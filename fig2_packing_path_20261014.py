#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""审计线论文 Fig 2 —— 定位：缺陷归属单一交付路径（Step 2）

对应 `plans/配图方案_20261014.md` **Figure 2**。
数据源（均已实测，见 `scripts/chanpred_20261007/RESULTS.md`）：

  · 200 个 `IP-` 文件内，`GE_ART` 与 `GE_ECG` 的分类交叉表
      两者都受损 124 / 都干净 70 / 仅 ART 5 / 仅 ECG 1
      Po = 0.970, Pe = 0.536, κ = +0.935
  · 60 个 `CB-` 文件（另一条通道集，94.3% 为空 ⇒ 3 例有波形）
      ECG1  ac1 中位 +0.975, 受损 0
      PLETH ac1 中位 +0.988, 受损 0

**本脚本不需原始波形**（全部数字来自已冻结的分析产出），可在本机运行。

用法：
  python3 fig2_packing_path_20261014.py --out <dir>
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

# ---- 冻结的分析结果（scripts/chanpred_20261007/RESULTS.md）----
N_BOTH_AFFLICTED = 124
N_BOTH_CLEAN = 70
N_ART_ONLY = 5
N_ECG_ONLY = 1
N_FILES = N_BOTH_AFFLICTED + N_BOTH_CLEAN + N_ART_ONLY + N_ECG_ONLY   # 200
PO = 0.970
PE = 0.536
KAPPA = 0.935

# CB- 负对照通道集（n=3 有波形）
CB_CHANNELS = [("ECG1", 3, 0.975, 0), ("PLETH", 3, 0.988, 0)]


def fig2(out_dir, dpi=400):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    fig = plt.figure(figsize=(7.2, 8.0))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.05, 1.0, 0.95], hspace=0.55)

    # ---------- A：2×2 交叉表（mosaic 式，面积∝计数）----------
    # 坐标：x=0 表示 ECG 受损，x=1 表示 ECG 干净；y=0 表示 ART 受损，y=1 表示 ART 干净
    axA = fig.add_subplot(gs[0])
    total = N_FILES
    cells = [
        # (x, y, count, colour, label)
        (0, 0, N_BOTH_AFFLICTED, "#C62828", "both affected"),
        (1, 1, N_BOTH_CLEAN, "#2E7D32", "both clean"),
        (0, 1, N_ECG_ONLY, "#9E9E9E", "ECG only"),
        (1, 0, N_ART_ONLY, "#9E9E9E", "ART only"),
    ]
    for x, y, n, col, lab in cells:
        w = h = np.sqrt(n / total)          # 面积 ∝ 计数
        axA.add_patch(Rectangle((x - w / 2, y - h / 2), w, h,
                                facecolor=col, edgecolor="white", lw=1.4, alpha=0.85))
        txtcol = "white" if n > 20 else "black"
        axA.text(x, y, f"{lab}\n{n}", ha="center", va="center",
                 fontsize=8.5, color=txtcol, weight="bold")

    axA.set_xlim(-0.8, 1.8)
    axA.set_ylim(-0.8, 1.9)
    axA.set_xticks([0, 1]); axA.set_xticklabels(["ECG affected", "ECG clean"], fontsize=8)
    axA.set_yticks([0, 1]); axA.set_yticklabels(["ART\naffected", "ART\nclean"], fontsize=8)
    axA.set_title(f"A  Co-affliction within the same files (n = {N_FILES})",
                  fontsize=10, loc="left")
    axA.grid(alpha=0.15, lw=0.4)
    axA.tick_params(length=0)
    axA.text(0.5, 1.85, f"Po = {PO:.3f}   Pe = {PE:.3f}   κ = {KAPPA:+.3f}    "
                        f"discordant {N_ART_ONLY + N_ECG_ONLY} of {N_FILES}",
             ha="center", va="top", fontsize=8.5,
             bbox=dict(boxstyle="round,pad=0.35", fc="#F5F5F5", ec="#BDBDBD", lw=0.7))

    # ---------- B：逐文件类别一致率（比 Po/Pe 对比更直观）----------
    axB = fig.add_subplot(gs[1])
    cats = ["concordant", "discordant"]
    vals = [N_BOTH_AFFLICTED + N_BOTH_CLEAN, N_ART_ONLY + N_ECG_ONLY]
    cols = ["#2E7D32", "#9E9E9E"]
    bars = axB.barh(cats, vals, color=cols, height=0.45, edgecolor="white", lw=0.8)
    for i, (b, v) in enumerate(zip(bars, vals)):
        axB.text(v + 2.5, b.get_y() + b.get_height() / 2,
                 f"{v} of {N_FILES}  ({v / N_FILES * 100:.1f}%)",
                 va="center", fontsize=9, weight="bold")
    axB.set_xlim(0, 232)
    axB.set_xlabel("files (of 200)", fontsize=8.5)
    axB.set_title("B  The two channels agree in 97% of files", fontsize=10, loc="left")
    axB.grid(alpha=0.22, lw=0.4, axis="x")
    axB.tick_params(labelsize=8.5)

    # ---------- C：另一通道集（负对照）----------
    axC = fig.add_subplot(gs[2])
    names = [c[0] for c in CB_CHANNELS]
    acs = [c[2] for c in CB_CHANNELS]
    axC.bar(names, acs, color="#2E7D32", width=0.4, edgecolor="white", lw=0.8)
    axC.axhline(0.5, color="#C62828", ls="--", lw=1.2)
    axC.text(2.32, 0.5, "affliction threshold (+0.5)", fontsize=7.5,
             color="#C62828", ha="left", va="center")
    axC.set_xlim(-0.55, 3.6)
    for i, (nm, n, ac, aff) in enumerate(CB_CHANNELS):
        axC.text(i, ac + 0.03, f"median ac1 = {ac:+.3f}\nn = {n}, affected = {aff}",
                 ha="center", fontsize=8)
    axC.set_ylabel("median lag-1\nautocorrelation", fontsize=8.5)
    axC.set_ylim(0, 1.35)
    axC.set_title("C  The other channel set is clean (negative control)",
                  fontsize=10, loc="left")
    axC.grid(alpha=0.22, lw=0.4, axis="y")
    axC.tick_params(labelsize=8.5)

    os.makedirs(out_dir, exist_ok=True)
    png = os.path.join(out_dir, "Fig2_packing_path_20261014.png")
    pdf = os.path.join(out_dir, "Fig2_packing_path_20261014.pdf")
    fig.savefig(png, dpi=dpi, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")

    summary = {
        "panel_a": {"n_files": N_FILES, "both_afflicted": N_BOTH_AFFLICTED,
                    "both_clean": N_BOTH_CLEAN, "art_only": N_ART_ONLY,
                    "ecg_only": N_ECG_ONLY, "Po": PO, "Pe": PE, "kappa": KAPPA},
        "panel_b": {"concordant": N_BOTH_AFFLICTED + N_BOTH_CLEAN,
                    "discordant": N_ART_ONLY + N_ECG_ONLY,
                    "concordant_pct": round((N_BOTH_AFFLICTED + N_BOTH_CLEAN)
                                            / N_FILES * 100, 2)},
        "chance_terms": {"Po": PO, "Pe": PE, "kappa": KAPPA},
        "panel_c": [{"channel": n, "n": k, "median_ac1": a, "affected": f}
                    for n, k, a, f in CB_CHANNELS],
    }
    with open(os.path.join(out_dir, "fig2_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary, indent=1))
    print("WROTE", png)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--dpi", type=int, default=400)
    a = ap.parse_args()
    fig2(a.out, a.dpi)
