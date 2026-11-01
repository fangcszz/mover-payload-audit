#!/usr/bin/env python3
"""Fig 4 — Design-level bound and its failure on a second corpus.

Dual-curve log-log plot: MOVER relSE rises as anchor windows narrow;
VitalDB relSE falls.  The bound is conditional, not general.

Data sources (all local in mirror/):
  MOVER:  mirror/results/run-20260913-anchor/40_resolution/wp4b_identifiability_by_anchor.json
          (ART channel, two measured grids: 300 s and 60 s native)
  VitalDB: manuscript §3.5 (0.095 at 300 s, 0.026 at 30 s)

Power-law exponents (manuscript §3.5):
  sigma_res  ∝ Delta^(-0.331)
  S_xx * a^2 ∝ n^(0.264)
  ⇒ relSE = sigma_res / sqrt(S_xx * a^2) ∝ Delta^(-0.331) / n^(0.132)
  with n ∝ T/Delta ⇒ relSE ∝ Delta^(-0.199)   [MOVER, rising as Delta narrows]

Usage:
  python3 fig4_design_bound_20261014.py
Produces:
  scratch/figs_20261014/Fig4_design_bound_20261014.{png,pdf}
  scratch/figs_20261014/fig4_summary.json
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── paths ────────────────────────────────────────────────────────────────
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("MOVER_PROJECT_ROOT",
                      os.path.dirname(os.path.dirname(HERE)))
ANCHOR_JSON = os.path.join(
    ROOT,
    "mirror/results/run-20260913-anchor/40_resolution/wp4b_identifiability_by_anchor.json",
)
OUTDIR = os.environ.get("MOVER_FIG_OUT", os.path.join(ROOT, "scratch/figs_20261014"))
os.makedirs(OUTDIR, exist_ok=True)

# ── load MOVER measured relSE ────────────────────────────────────────────
with open(ANCHOR_JSON) as fh:
    idf = json.load(fh)

art = idf["ART"]

mover_300 = dict(
    delta_s=300.0,
    relse_p50=art["rel_se_a"]["p50"],
    relse_p25=art["rel_se_a"]["p25"],
    relse_p75=art["rel_se_a"]["p75"],
    r2_p50=art["r2"]["p50"],
    n_pts_p50=art["n_bins_fit"]["p50"],
    n=art["rel_se_a"]["n"],
)
mover_60 = dict(
    delta_s=60.0,
    relse_p50=art["rel_se_a60f"]["p50"],
    relse_p25=art["rel_se_a60f"]["p25"],
    relse_p75=art["rel_se_a60f"]["p75"],
    r2_p50=art["r2_60f"]["p50"],
    n_pts_p50=art["n1_fit"]["p50"],
    n=art["rel_se_a60f"]["n"],
)

# ── VitalDB (from manuscript §3.5) ───────────────────────────────────────
vitaldb_300 = dict(delta_s=300.0, relse=0.095)
vitaldb_30 = dict(delta_s=30.0, relse=0.026)

# ── power-law fits ───────────────────────────────────────────────────────
# MOVER:  exponent from sub-exponents (sigma_res ∝ Δ^-0.331, S_xx·a² ∝ n^0.264)
#   relSE ∝ Δ^(-0.331 + 0.132) = Δ^(-0.199)
#   Verify from two measured points: log(1.788/1.299)/log(60/300) = -0.199  ✓
mover_exp = -0.199
mover_A = mover_300["relse_p50"] / (mover_300["delta_s"] ** mover_exp)

# VitalDB: fit from two points
vitaldb_exp = np.log(vitaldb_300["relse"] / vitaldb_30["relse"]) / np.log(
    vitaldb_300["delta_s"] / vitaldb_30["delta_s"]
)
vitaldb_A = vitaldb_300["relse"] / (vitaldb_300["delta_s"] ** vitaldb_exp)

# curve x-range
delta_curve = np.logspace(np.log10(20), np.log10(400), 200)
mover_curve = mover_A * delta_curve ** mover_exp
vitaldb_curve = vitaldb_A * delta_curve ** vitaldb_exp

# ── plot ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7.0, 5.2))

# MOVER curve + measured points
ax.plot(delta_curve, mover_curve, color="#c0392b", lw=2.2, zorder=3,
        label="MOVER (this archive)")
ax.errorbar(
    [mover_300["delta_s"], mover_60["delta_s"]],
    [mover_300["relse_p50"], mover_60["relse_p50"]],
    yerr=[
        [mover_300["relse_p50"] - mover_300["relse_p25"],
         mover_60["relse_p50"] - mover_60["relse_p25"]],
        [mover_300["relse_p75"] - mover_300["relse_p50"],
         mover_60["relse_p75"] - mover_60["relse_p50"]],
    ],
    fmt="o", color="#c0392b", ms=7, capsize=4, zorder=4,
)

# VitalDB curve + measured points
ax.plot(delta_curve, vitaldb_curve, color="#2471a3", lw=2.2, ls="--", zorder=3,
        label="VitalDB (second corpus)")
ax.plot(
    [vitaldb_300["delta_s"], vitaldb_30["delta_s"]],
    [vitaldb_300["relse"], vitaldb_30["relse"]],
    "s", color="#2471a3", ms=7, zorder=4,
)

# direction arrows
ax.annotate(
    "rises as\nwindows narrow",
    xy=(45, mover_A * 45 ** mover_exp),
    xytext=(70, 4.5),
    fontsize=8.5, color="#c0392b", ha="center",
    arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.3),
)
ax.annotate(
    "falls as\nwindows narrow",
    xy=(40, vitaldb_A * 40 ** vitaldb_exp),
    xytext=(120, 0.012),
    fontsize=8.5, color="#2471a3", ha="center",
    arrowprops=dict(arrowstyle="->", color="#2471a3", lw=1.3),
)

# exponent annotations
ax.text(0.97, 0.97,
        r"$\mathbf{MOVER:}$" + "\n" +
        r"relSE $\propto \Delta^{%.3f}$" % mover_exp + "\n" +
        r"($\sigma_{res}\propto\Delta^{-.331}$, "
        r"$S_{xx}a^2\propto n^{.264}$)",
        transform=ax.transAxes, fontsize=8, color="#c0392b",
        va="top", ha="right",
        bbox=dict(fc="white", ec="#c0392b", alpha=0.9, lw=0.8))
ax.text(0.97, 0.20,
        r"$\mathbf{VitalDB:}$" + "\n" +
        r"relSE $\propto \Delta^{%.3f}$" % vitaldb_exp,
        transform=ax.transAxes, fontsize=8, color="#2471a3",
        va="top", ha="right",
        bbox=dict(fc="white", ec="#2471a3", alpha=0.9, lw=0.8))

# ~50x gap annotation
gap_300 = mover_300["relse_p50"] / vitaldb_300["relse"]
ax.annotate(
    "",
    xy=(300, mover_300["relse_p50"]),
    xytext=(300, vitaldb_300["relse"]),
    arrowprops=dict(arrowstyle="<->", color="grey", lw=1.0),
)
ax.text(285, np.sqrt(mover_300["relse_p50"] * vitaldb_300["relse"]),
        "%.0f× gap" % gap_300, fontsize=8, color="grey", va="center", ha="right")

# axes
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("anchor bin width  $\\Delta$  (s)", fontsize=11)
ax.set_ylabel("relative slope SE   relSE($\\Delta$)", fontsize=11)
ax.set_xlim(17, 460)
ax.set_xticks([20, 30, 60, 100, 300, 400])
ax.set_xticklabels(["20", "30", "60", "100", "300", "400"])
ax.set_ylim(0.008, 12)
ax.legend(loc="center left", fontsize=9, framealpha=0.9)
ax.tick_params(labelsize=9)

# grid
ax.grid(True, which="major", ls="-", lw=0.4, alpha=0.3)
ax.grid(True, which="minor", ls=":", lw=0.3, alpha=0.2)

fig.tight_layout()

# ── save ─────────────────────────────────────────────────────────────────
png = os.path.join(OUTDIR, "Fig4_design_bound_20261014.png")
pdf = os.path.join(OUTDIR, "Fig4_design_bound_20261014.pdf")
fig.savefig(png, dpi=400)
fig.savefig(pdf)
plt.close(fig)

# ── summary JSON ─────────────────────────────────────────────────────────
summary = {
    "mover": {
        "measured": [mover_300, mover_60],
        "power_law": {"exponent": mover_exp, "A": mover_A},
    },
    "vitaldb": {
        "measured": [vitaldb_300, vitaldb_30],
        "power_law": {"exponent": vitaldb_exp, "A": vitaldb_A},
    },
    "gap_at_300s": gap_300,
    "interpretation": "MOVER relSE rises as windows narrow (exponent < 0); "
                      "VitalDB relSE falls (exponent > 0). "
                      "The bound is conditional, not general.",
}
sj = os.path.join(OUTDIR, "fig4_summary.json")
with open(sj, "w") as fh:
    json.dump(summary, fh, indent=2)

print("OK")
print("  PNG:", png)
print("  PDF:", pdf)
print("  JSON:", sj)
print("  MOVER exponent: %.4f" % mover_exp)
print("  VitalDB exponent: %.4f" % vitaldb_exp)
print("  gap at 300s: %.1fx" % gap_300)
