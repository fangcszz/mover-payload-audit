#!/usr/bin/env python3
"""B: byte-level spot-check of GE_ART payloads (no ac1 classifier).

Randomly samples N GE_ART files, decodes the first packet's base64 to int16,
and prints the first codes with even/odd subsequences, so a human can see
whether the period-2 alternation is present.
"""
import base64
import csv
import gzip
import os
import random
import re
import struct
import sys

BASE = os.environ.get("MOVER_WAVE_ROOT",
                       "/mnt/share/ljfang/exp/raw/wave_pilot_20260909")
INDEX = os.environ.get(
    "MOVER_WAVE_INDEX",
    "/mnt/share/ljfang/exp/results/run-20260911-113000/"
    "11_wave_gates/wave_file_index_20260911.csv.gz")
MG_RE = re.compile(r'<mg name="GE_ART">(.*?)</mg>', re.S)
WAVE_RE = re.compile(r'<m name="Wave">([^<]+)</m>', re.S)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    rows = []
    with gzip.open(INDEX, "rt") as fh:
        for row in csv.DictReader(fh):
            if row["has_wave"] == "1" and "GE_ART" in row["channels"]:
                rows.append(row["relpath"])
    random.seed(20261030)
    sample = random.sample(rows, min(n, len(rows)))
    print("sampled %d of %d GE_ART files" % (len(sample), len(rows)), flush=True)
    for rp in sample:
        path = os.path.join(BASE, rp)
        try:
            with open(path, encoding="utf-8", errors="replace") as fh:
                text = fh.read()
        except Exception as e:
            print(os.path.basename(rp), "READ_ERR", e)
            continue
        mg = MG_RE.search(text)
        if not mg:
            print(os.path.basename(rp), "NO_GE_ART")
            continue
        w = WAVE_RE.search(mg.group(1))
        if not w:
            print(os.path.basename(rp), "NO_WAVE")
            continue
        raw = base64.b64decode(w.group(1))
        codes = list(struct.unpack("<%dh" % (len(raw) // 2), raw))
        even = codes[0:12:2]
        odd = codes[1:13:2]
        print("%s len=%d first12=%s even=%s odd=%s" % (
            os.path.basename(rp), len(codes), codes[:12], even, odd), flush=True)


if __name__ == "__main__":
    main()
