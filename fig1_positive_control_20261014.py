#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""审计线论文 Fig 1 —— 中心图：阳性对照与三类载荷

投稿规格重制（原分析级图在 gpu1 `/mnt/share/ljfang/exp/results/run-20260913-abpcontrol/`）。

本脚本**在 gpu1 上运行**（原始波形只在 gpu1 的 /mnt/share 下，本机无挂载）。

输出：Fig1（A 阳性对照 / B 三类载荷代表波形 / C 增益分布 + 生理窗）
      —— 对应 `plans/配图方案_20261014.md` Figure 1。

口径（2026-10-14 已定，见 `plans/口径修正_阳性对照分母_20261014.md`）：
  - S 类 24 例（检出）
  - 其中 21 例有正的增益拟合 → final_stats.json 的 S_subset.n
  - 其中 23 例有配对 charted MAP → 正文 r 的 n
  - 90.5% = 19/21（生理窗内）
本图只画 21 例口径的增益与 SBP/DBP/PP，图注须写明。

用法（gpu1）：
  python3 fig1_positive_control_20261014.py --out /mnt/share/ljfang/exp/results/run-20261014-figs
"""
import argparse
import collections
import csv
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 分析级产出目录（gpu1）
SRC = "/mnt/share/ljfang/exp/results/run-20260913-abpcontrol"
# 生成器与配对辅助模块（gpu1 上 _probe/ 内的工作副本）
HELPERS = f"{SRC}/_probe"

# 分类代表例：从 S/A 类名单取首例（A 类名单在 abp_final_figs.py 中给定）
S_REP = "e17daa9fac8a5985"
A_REP = "d93535750d53d8fe"
NUM_SHIFT = 15.0 * 3600.0  # 数值通道与波形的时基偏移（沿用分析脚本）


def load_tables(src):
    tbl = {r["log_id"]: r for r in json.load(open(f"{src}/Sfinal_table.json"))}
    stats = json.load(open(f"{src}/final_stats.json"))
    return tbl, stats


def panel_a(ax, tbl, stats):
    """阳性对照：拟合增益 vs 声明 0.25（21 例）。"""
    gains = np.array([r["gain_file"] for r in tbl.values()
                      if r.get("gain_file") == r.get("gain_file") and r["gain_file"] > 0])
    n = len(gains)
    assert n == stats["S_subset"]["n"], f"n mismatch: {n} vs {stats['S_subset']['n']}"

    med = float(np.median(gains))
    ax.hist(gains, bins=12, color="#4C72B0", edgecolor="white", linewidth=0.6)
    ax.axvline(0.25, color="#C44E52", lw=1.8, ls="--", label="declared gain = 0.25")
    ax.axvline(med, color="black", lw=1.4, ls="-",
               label=f"recovered median = {med:.4f}")
    ax.set_xlabel("recovered gain (per case)", fontsize=9)
    ax.set_ylabel("cases", fontsize=9)
    ax.set_title(f"A  Gain concordance, readable class (n = {n})", fontsize=10, loc="left")
    ax.legend(fontsize=7.5, frameon=False)
    ax.grid(alpha=0.22, lw=0.4)
    ax.tick_params(labelsize=8)
    return n, med


def _a_class_pid(src):
    """A 类代表例的 pid（A 类不在 Sfinal_table 中，从 demod 表读）。"""
    with open(f"{src}/demod_variants_per_case.csv") as fh:
        for x in csv.DictReader(fh):
            if x["log_id"] == A_REP:
                return x["pid"]
    return None


def _load_helpers():
    """载入 gpu1 _probe/ 下的原始波形读取辅助模块。"""
    sys.path.insert(0, HELPERS)
    from p1_final import packets_for          # noqa: E402
    from abp_pair import read_numeric, files_by_pid   # noqa: E402
    return packets_for, read_numeric, files_by_pid


def _window_from_files(pid, fbp, dt_utc, span=600.0):
    """由文件索引自身的时标取一个窗口（不依赖数值通道）。"""
    fl = fbp.get(pid, [])
    if not fl:
        return None
    times = []
    for f in fl:
        try:
            times.append((dt_utc(f["first_dt_utc"]).timestamp(),
                          dt_utc(f["last_dt_utc"]).timestamp()))
        except Exception:                          # noqa: BLE001
            continue
    if not times:
        return None
    times.sort()
    lo = times[0][0]
    return lo, lo + span


def _best_segment(codes, want=1440):
    """选一段**连续**的高幅度码流。

    关键：包之间可能不连续（不同文件/时段），必须只在连续高 ptp 的包里取，
    否则拼出的是跳跃的假波形。先按 ptp 找连续 run，再取该 run 内最长的段。
    """
    if not codes:
        return None
    code_list = [np.asarray(c, dtype=float) for c in codes]
    ptp = np.array([float(np.ptp(c)) if len(c) else 0.0 for c in code_list])
    thr = np.percentile(ptp[ptp > 0], 60) if (ptp > 0).any() else 0.0

    # 找连续的高幅度 run
    best_run, run = None, []
    for i, p in enumerate(ptp):
        if p >= thr:
            if run and i == run[-1] + 1:
                run.append(i)
            else:
                if best_run is None or len(run) > len(best_run):
                    best_run = list(run)
                run = [i]
    if best_run is None or len(run) > len(best_run):
        best_run = list(run)
    if not best_run:
        return None

    # 在 run 内拼接至 want 长度
    out = []
    total = 0
    for i in best_run:
        out.append(code_list[i])
        total += len(code_list[i])
        if total >= want:
            break
    if not out:
        return None
    seg = np.concatenate(out)
    return seg[:want] if len(seg) >= want else seg


def panel_b(ax, src, rows):
    """三类载荷代表波形：S（可读脉动）/ A（周期 2 交替）各取真实码流。

    A 类在 8 s 尺度上是 Nyquist 交替的实心带（raw_ac1 ≈ −0.99），
    因此**右半幅改为 ~0.15 s 放大窗**，让"逐样本交替、包络仍有压力波"
    这一事实可见——这正是 A 类被判定为不可读的原因。
    """
    try:
        sys.path.insert(0, HELPERS)
        from p1_final import packets_for                 # noqa: E402
        from abp_pair import read_numeric, files_by_pid, dt_utc   # noqa: E402
    except Exception as exc:                              # noqa: BLE001
        print("panel B: helper import failed:", exc, file=sys.stderr)
        return False
    fbp = files_by_pid()
    return _panel_b_impl(ax, src, rows, packets_for, read_numeric, files_by_pid, dt_utc, fbp)


def _panel_b_impl(ax, src, rows, packets_for, read_numeric, files_by_pid, dt_utc, fbp):
    """在给定轴上绘 S（全窗）+ A（放大窗）双面板。"""
    parent = ax.get_subplotspec()
    fig = ax.figure
    ax.remove()
    inner = parent.subgridspec(1, 2, wspace=0.28)
    axL = fig.add_subplot(inner[0])
    axR = fig.add_subplot(inner[1])

    drewL = drewR = False

    # --- 左：S 类，8 s 全窗 ---
    if S_REP in rows:
        r = rows[S_REP]
        try:
            arr, _s = read_numeric(S_REP)
            pk = packets_for(r["pid"], arr[0, 0] + NUM_SHIFT, arr[-1, 0] + NUM_SHIFT, fbp)
            seg = _best_segment([c for _, c in pk], want=1440)
            if seg is not None:
                t = np.arange(len(seg)) / 180.0
                axL.plot(t, seg, lw=0.85, color="#2E7D32")
                drewL = True
        except Exception as exc:                          # noqa: BLE001
            print(f"panel B: S failed: {exc}", file=sys.stderr)

    # --- 右：A 类，放大窗 ~0.15 s（27 samples）---
    try:
        pid_a = _a_class_pid(src)
        if pid_a:
            fl = fbp.get(pid_a, [])
            if fl:
                t0 = min(dt_utc(f["first_dt_utc"]).timestamp() for f in fl)
                pk = packets_for(pid_a, t0, t0 + 600, fbp)
                if pk:
                    c = np.asarray(pk[0][1], dtype=float)
                    n = min(27, len(c))
                    t = np.arange(n) / 180.0 * 1000.0   # ms
                    axR.plot(t, c[:n], lw=1.0, marker="o", ms=3.0,
                             color="#C62828", mfc="white", mew=0.9, label="raw samples")
                    axR.plot(t, c[:n], lw=0.7, color="#C62828", alpha=0.5)
                    # 偶/奇子序列包络
                    axR.plot(t[0::2], c[0:n:2], lw=0.5, ls=":", color="#1565C0",
                             label="even samples")
                    axR.plot(t[1::2], c[1:n:2], lw=0.5, ls=":", color="#EF6C00",
                             label="odd samples")
                    drewR = True
    except Exception as exc:                              # noqa: BLE001
        print(f"panel B: A failed: {exc}", file=sys.stderr)

    if not drewL and not drewR:
        return False

    axL.set_xlabel("time within segment (s)", fontsize=8.5)
    axL.set_ylabel("raw code", fontsize=8.5)
    axL.set_title("B1  Class S: readable pulsatile (8 s)", fontsize=9, loc="left")
    axL.grid(alpha=0.22, lw=0.4)
    axL.tick_params(labelsize=7.5)

    axR.set_xlabel("time (ms)", fontsize=8.5)
    axR.set_title("B2  Class A: period-2 alternation (144 ms)", fontsize=9, loc="left")
    axR.grid(alpha=0.22, lw=0.4)
    axR.tick_params(labelsize=7.5)
    axR.set_ylim(-560, 620)
    axR.legend(fontsize=6.0, frameon=False, loc="upper center", ncol=3,
               bbox_to_anchor=(0.5, 1.0), handlelength=1.1,
               columnspacing=0.8, borderaxespad=0.15)
    return True


def panel_c(ax, src):
    """生理窗落点：SBP / DBP / PP 箱线 vs 预定窗（21 例）。"""
    rows = json.load(open(f"{src}/Sfinal_table.json"))
    sbp = np.array([r["g0.25_SBP"] for r in rows if r.get("g0.25_SBP")])
    dbp = np.array([r["g0.25_DBP"] for r in rows if r.get("g0.25_DBP")])
    pp = np.array([r["g0.25_PP"] for r in rows if r.get("g0.25_PP")])
    n = len(sbp)
    if n == 0:
        return None

    inside = ((sbp >= 90) & (sbp <= 180) & (dbp >= 40) & (dbp <= 100) &
              (pp >= 30) & (pp <= 100))
    frac = float(inside.mean())

    bp = ax.boxplot([sbp, dbp, pp], tick_labels=["SBP", "DBP", "PP"],
                    widths=0.5, patch_artist=True, showfliers=True,
                    medianprops=dict(color="black", lw=1.4))
    for patch, col in zip(bp["boxes"], ["#4C72B0", "#55A868", "#C44E52"]):
        patch.set_facecolor(col)
        patch.set_alpha(0.65)

    # 预定生理窗（阴影）
    for i, (lo, hi) in enumerate([(90, 180), (40, 100), (30, 100)], start=1):
        ax.add_patch(plt.Rectangle((i - 0.32, lo), 0.64, hi - lo,
                                   facecolor="#F0F0F0", edgecolor="none", zorder=0))

    ax.set_ylabel("mmHg", fontsize=9)
    ax.set_title(f"C  Decoded pressures vs window ({frac * 100:.1f}% inside, n = {n})",
                 fontsize=10, loc="left")
    ax.grid(alpha=0.22, lw=0.4, axis="y")
    ax.tick_params(labelsize=8)
    return frac, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="输出目录（gpu1 /mnt/share 下）")
    ap.add_argument("--src", default=SRC)
    ap.add_argument("--dpi", type=int, default=400)
    args = ap.parse_args()

    if not os.path.isdir(args.src):
        sys.exit(f"SRC not found: {args.src} (must run on gpu1)")
    os.makedirs(args.out, exist_ok=True)

    tbl, stats = load_tables(args.src)

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig = plt.figure(figsize=(7.2, 8.6))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 1.0], hspace=0.42)

    ax_a = fig.add_subplot(gs[0])
    n_a, med_a = panel_a(ax_a, tbl, stats)

    ax_b = fig.add_subplot(gs[1])
    has_b = panel_b(ax_b, args.src, tbl)
    if not has_b:
        ax_b.text(0.5, 0.5,
                  "Panel B: representative payloads\n"
                  "(source npz not found — generate from raw packets)",
                  ha="center", va="center", fontsize=8, color="#888888",
                  transform=ax_b.transAxes)
        ax_b.set_title("B  Representative payloads (pending data)", fontsize=10, loc="left")
        ax_b.set_xticks([]); ax_b.set_yticks([])

    ax_c = fig.add_subplot(gs[2])
    res_c = panel_c(ax_c, args.src)

    out_png = os.path.join(args.out, "Fig1_positive_control_20261014.png")
    out_pdf = os.path.join(args.out, "Fig1_positive_control_20261014.pdf")
    fig.savefig(out_png, dpi=args.dpi, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")

    summary = {
        "panel_a": {"n": n_a, "recovered_median": round(med_a, 6)},
        "panel_b": {"rendered": bool(has_b)},
        "panel_c": ({"frac_inside": round(res_c[0], 6), "n": res_c[1]}
                    if res_c else None),
        "denominators": {"S_class": 24, "gain_fit": 21, "paired_map": 23},
        "source": args.src,
    }
    with open(os.path.join(args.out, "fig1_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary, indent=1))
    print("WROTE", out_png)


if __name__ == "__main__":
    main()
