#!/usr/bin/env python3
"""A: full-archive ac1 census of GE_ART payloads (MAP-independent).

Reads the 926,147-row index, keeps has_wave=1 & GE_ART, then for each file:
  - regex-extract every <mg name="GE_ART"> block
  - regex-extract <m name="Wave"> base64, decode to int16, compute lag-1 autocorr
  - file class = median packet ac1; A: < -0.5, N: |a|<=0.5, S: > 0.5
No numeric MAP is required, so this runs over the whole archive.

Usage:
  python3 ac1_census_20261030.py --limit 100 --workers 8
  python3 ac1_census_20261030.py --workers 24
"""
import argparse
import base64
import csv
import gzip
import os
import re
import statistics
import struct
from concurrent.futures import ProcessPoolExecutor

BASE = os.environ.get("MOVER_WAVE_ROOT",
                       "/mnt/share/ljfang/exp/raw/wave_pilot_20260909")
INDEX = os.environ.get(
    "MOVER_WAVE_INDEX",
    "/mnt/share/ljfang/exp/results/run-20260911-113000/"
    "11_wave_gates/wave_file_index_20260911.csv.gz")

MG_RE = re.compile(r'<mg name="GE_ART">(.*?)</mg>', re.S)
WAVE_RE = re.compile(r'<m name="Wave">([^<]+)</m>', re.S)


def ac1(codes):
    n = len(codes)
    if n < 3:
        return None
    m = sum(codes) / n
    x = [c - m for c in codes]
    d = sum(v * v for v in x)
    if d <= 0:
        return None
    return sum(x[i] * x[i + 1] for i in range(n - 1)) / d


def classify(a):
    if a is None:
        return "NA"
    if a < -0.5:
        return "A"
    if a > 0.5:
        return "S"
    return "N"


def process(relpath):
    path = os.path.join(BASE, relpath)
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
    except Exception:
        return (relpath, 0, None, "NA")
    acs = []
    for mg in MG_RE.finditer(text):
        w = WAVE_RE.search(mg.group(1))
        if not w:
            continue
        try:
            raw = base64.b64decode(w.group(1))
            codes = list(struct.unpack("<%dh" % (len(raw) // 2), raw))
        except Exception:
            continue
        a = ac1(codes)
        if a is not None:
            acs.append(a)
    if not acs:
        return (relpath, 0, None, "NA")
    med = statistics.median(acs)
    return (relpath, len(acs), med, classify(med))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="0 = all")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--out", default=os.environ.get(
        "MOVER_OUT", "/mnt/share/ljfang/exp/results/run-20261030-gate"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    rows = []
    seen = set()
    with gzip.open(INDEX, "rt") as fh:
        for row in csv.DictReader(fh):
            if row["has_wave"] == "1" and "GE_ART" in row["channels"]:
                key = (row["pid"], row["suffix"], row["fn_ts"])
                if key in seen:
                    continue
                seen.add(key)
                rows.append(row["relpath"])
    print("GE_ART unique files after de-dup:", len(rows), flush=True)
    if a.limit:
        rows = rows[: a.limit]
        print("limited to:", len(rows), flush=True)

    tag = "all" if not a.limit else "n%d" % len(rows)
    outcsv = os.path.join(a.out, "ac1_census_%s.csv" % tag)
    done = 0
    with open(outcsv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["relpath", "n_packets", "median_ac1", "class"])
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            for res in ex.map(process, rows, chunksize=200):
                if res:
                    w.writerow(res)
                done += 1
                if done % 5000 == 0:
                    print(done, flush=True)
    print("done:", outcsv, "rows:", done, flush=True)


if __name__ == "__main__":
    main()
