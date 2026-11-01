#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""审计线论文 Fig 3 —— 影响量化：绝对单位在任何估计量类下都不可恢复（Step 3）

对应 `plans/配图方案_20261014.md` **Figure 3**。
数据源：`mirror/results/run-20260913-audit/22_wp4_wcal/stats_3arms_FULL.json`
（六臂恢复统计，n=1,058 例 / 41,047 分箱）+ `per_case_arms_FULL.csv`。

**本脚本不需原始波形**，可在本机运行。

口径要点（如实呈现，不美化）：
  · 乐观臂（file_is / case_is）在同一批分箱上拟合与评分 ⇒ **不是可部署估计**
  · 天花板 file_is LoA 半宽 25.07 mmHg = 阈值 10 的 2.5×，且 ≈ charted MAP 自身边缘 SD 20.79
  · 可部署臂（file_oos / case_oos / case_oos60）LoA 半宽 691–876 mmHg，**比不标定还差一个量级**
  · raw 臂 bias −62.2、LoA 半宽 61.4、ICC 0.0125、93.9% 分箱 <65 mmHg ⇒ **归谬而非一致性研究**

用法：
  python3 fig3_recoverability_20261014.py --out <dir> [--stats <json>] [--percase <csv>]
"""
import argparse
import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

THRESHOLD = 10.0          # 预设验收阈值：LoA 半宽 ≤ 10 mmHg
MARGINAL_SD = 20.79       # charted MAP 自身边缘 SD

ARM_LABELS = {
    "raw":        "released decoder\n(no recalibration)",
    "file_is":    "per-file, in-sample\n(ceiling)",
    "case_is":    "per-case, in-sample",
    "file_oos":   "per-file, out-of-sample",
    "case_oos":   "per-case, out-of-sample",
    "case_oos60": "per-case, out-of-sample\n60 s anchor",
}
DEPLOYABLE = {"file_oos", "case_oos", "case_oos60"}


def fig3(out_dir, stats_path, percase_path, dpi=400):
    S = json.load(open(stats_path))
    arms = S["arms"]
    order = ["raw", "file_is", "case_is", "file_oos", "case_oos", "case_oos60"]

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    fig = plt.figure(figsize=(7.2, 9.2))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.15, 1.0, 1.0], hspace=0.62)

    # ---------- A：LoA 半宽（对数轴，阈值 10 与 MAP 边缘 SD 并列）----------
    axA = fig.add_subplot(gs[0])
    hw = [arms[a]["overall"]["loa_halfwidth"] for a in order]
    cols = ["#C62828" if a in DEPLOYABLE or a == "raw" else "#1565C0" for a in order]
    ypos = np.arange(len(order))[::-1]
    bars = axA.barh(ypos, hw, color=cols, height=0.6, edgecolor="white", lw=0.8)
    axA.set_xscale("log")
    axA.set_yticks(ypos)
    axA.set_yticklabels([ARM_LABELS[a] for a in order], fontsize=7.6)
    for y, a, v in zip(ypos, order, hw):
        axA.text(v * 1.12, y, f"{v:,.1f}" if v < 100 else f"{v:,.0f}",
                 va="center", fontsize=8, weight="bold")
    axA.axvline(THRESHOLD, color="#2E7D32", ls="--", lw=1.4)
    axA.axvline(MARGINAL_SD, color="#616161", ls=":", lw=1.3)
    axA.text(THRESHOLD * 1.12, -0.68, f"acceptance threshold {THRESHOLD:.0f}",
             fontsize=7, color="#2E7D32", ha="left")
    axA.text(MARGINAL_SD * 1.02, -1.22, f"charted MAP marginal SD {MARGINAL_SD}",
             fontsize=7, color="#616161", ha="left")
    axA.set_xlim(5, 3000)
    axA.set_ylim(-1.35, len(order) - 0.2)
    axA.set_xlabel("limits-of-agreement half-width (mmHg, log scale)", fontsize=8.5)
    axA.set_title("A  Recovery of absolute pressure, by estimator class",
                  fontsize=10, loc="left")
    axA.grid(alpha=0.22, lw=0.4, axis="x")

    # ---------- B：raw 臂 Bland–Altman ----------
    axB = fig.add_subplot(gs[1])
    raw = arms["raw"]["overall"]
    # 用 per-case 表绘制逐例分布（避免依赖巨大 npz）
    md, lo, hi = raw["bias"], raw["loa_lo"], raw["loa_hi"]
    if os.path.exists(percase_path):
        vals = []
        with open(percase_path) as fh:
            for r in csv.DictReader(fh):
                try:
                    vals.append(float(r.get("raw_bias", "nan")))
                except ValueError:
                    pass
        vals = np.array([v for v in vals if np.isfinite(v)])
        if len(vals) > 10:
            axB.scatter(vals, np.random.default_rng(0).normal(0, 0.35, len(vals)),
                        s=5, alpha=0.35, color="#C62828", edgecolors="none")
    axB.axvline(md, color="black", lw=1.4, label=f"bias = {md:.1f}")
    axB.axvline(lo, color="#C62828", ls="--", lw=1.1, label=f"LoA = {lo:.1f}")
    axB.axvline(hi, color="#C62828", ls="--", lw=1.1, label=f"LoA = {hi:.1f}")
    axB.set_xlabel("waveform − numeric, per 5-min bin (mmHg)", fontsize=8.5)
    axB.set_yticks([])
    axB.set_title(f"B  The released decoder against the numeric record "
                  f"(bias {md:.1f} mmHg, half-width {raw['loa_halfwidth']:.1f}, "
                  f"ICC {raw['icc21']:.3f})", fontsize=10, loc="left")
    axB.legend(fontsize=7.5, frameon=False, loc="upper left")
    axB.grid(alpha=0.22, lw=0.4, axis="x")

    # ---------- C：5-min 分箱压力分布 + 65 mmHg 阈值 ----------
    axC = fig.add_subplot(gs[2])
    wq, eq = raw["w_quantiles"], raw["e_quantiles"]
    xs = [0, 1]
    for i, (q, col, lab) in enumerate([(wq, "#C62828", "waveform arm"),
                                       (eq, "#1565C0", "numeric arm")]):
        axC.plot([i, i], [q["p25"], q["p75"]], color=col, lw=7, solid_capstyle="butt",
                 alpha=0.55)
        axC.plot([i, i], [q["p1"], q["p99"]], color=col, lw=1.4)
        axC.plot(i, q["p50"], "o", color=col, ms=7, zorder=5)
        axC.text(i + 0.14, q["p50"], f"median {q['p50']:.1f}", fontsize=8,
                 va="center", color=col)
    axC.axhline(65, color="#2E7D32", ls="--", lw=1.4)
    axC.text(0.02, 0.97, f"hypotension threshold 65 mmHg; "
                         f"{raw['frac_w_below_65'] * 100:.1f}% of waveform bins below",
             transform=axC.transAxes, fontsize=7.2, color="#2E7D32",
             ha="left", va="top")
    axC.set_xticks(xs)
    axC.set_xticklabels(["waveform arm\n(released decoder)", "numeric arm\n(reference)"],
                        fontsize=8)
    axC.set_xlim(-0.5, 1.5)
    axC.set_ylabel("5-min bin pressure (mmHg)", fontsize=8.5)
    axC.set_title("C  What a non-auditing user inherits (reductio, not agreement)",
                  fontsize=10, loc="left")
    axC.grid(alpha=0.22, lw=0.4, axis="y")

    os.makedirs(out_dir, exist_ok=True)
    png = os.path.join(out_dir, "Fig3_recoverability_20261014.png")
    pdf = os.path.join(out_dir, "Fig3_recoverability_20261014.pdf")
    fig.savefig(png, dpi=dpi, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")

    summary = {
        "threshold_mmHg": THRESHOLD,
        "charted_map_marginal_sd": MARGINAL_SD,
        "arms": {a: {"bias": round(arms[a]["overall"]["bias"], 3),
                     "loa_halfwidth": round(arms[a]["overall"]["loa_halfwidth"], 3),
                     "icc21": round(arms[a]["overall"]["icc21"], 5),
                     "n_bins": arms[a]["overall"]["n"],
                     "deployable": a in DEPLOYABLE} for a in order},
        "raw_frac_below_65": round(raw["frac_w_below_65"], 5),
    }
    with open(os.path.join(out_dir, "fig3_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary, indent=1))
    print("WROTE", png)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--stats", default="mirror/results/run-20260913-audit/22_wp4_wcal/stats_3arms_FULL.json")
    ap.add_argument("--percase", default="mirror/results/run-20260913-audit/22_wp4_wcal/per_case_arms_FULL.csv")
    ap.add_argument("--dpi", type=int, default=400)
    a = ap.parse_args()
    fig3(a.out, a.stats, a.percase, a.dpi)
