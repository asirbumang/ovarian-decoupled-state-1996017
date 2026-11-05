# -*- coding: utf-8 -*-
"""Regenerate Fig6b / Fig8a / Fig8b at their final physical panel sizes (fonts 7 pt)."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

ROOT = r"D:\子宫内膜异位症2"
FIGS = os.path.join(ROOT, "analysis", "L1_scRNA", "results", "figs")
plt.rcParams.update({"font.family": "Arial", "font.size": 7,
                     "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                     "pdf.fonttype": 42, "ps.fonttype": 42})
C = {"blue": "#3C5488", "red": "#E64B35", "grey": "#8491B4"}

# ---------------- Fig 6b: two forests, final width ~70 mm (2.75 in)
tier = pd.read_csv(os.path.join(ROOT, "analysis", "L0_bulk", "results", "survival_meta",
                                "survival_meta_results.csv"))
cont = pd.read_csv(os.path.join(ROOT, "analysis", "L0_bulk", "results", "track2", "decreg_cox_meta.csv"))
order = ["GSE26193", "GSE14764", "GSE13876", "GSE63885", "TCGA-OV"]
tier = tier.set_index("cohort").loc[order].reset_index()
cont = cont.set_index("cohort").loc[order].reset_index()
fig, ax = plt.subplots(1, 2, figsize=(2.75, 1.05), sharey=True)
ys = np.arange(len(order))[::-1]
for a, df, hr, lo, hi, pcol in [(ax[0], tier, "hr_tier", "ci_lo", "ci_hi", "p_trend"),
                                (ax[1], cont, "hr_cont", "ci_lo", "ci_hi", "p_cont")]:
    for y, (_, r) in zip(ys, df.iterrows()):
        a.plot([r[lo], r[hi]], [y, y], color="black", lw=0.6)
        a.plot(r[hr], y, "o", color=C["blue"], ms=2)
        a.text(2.55, y, f"{r[pcol]:.3f}".lstrip("0"), va="center", fontsize=5.2)
    a.axvline(1, color="black", lw=0.4, ls="--")
    a.set_xscale("log"); a.set_xlim(0.45, 2.3)
    a.set_xticks([0.5, 1, 2]); a.set_xticklabels(["0.5", "1", "2"], fontsize=6)
    a.tick_params(axis="y", length=0)
    for s in ["top", "right"]:
        a.spines[s].set_visible(False)
ax[0].set_yticks(ys); ax[0].set_yticklabels(order, fontsize=5.6)
ax[0].set_title("raw (groups)", fontsize=6.2, loc="left")
ax[1].set_title("adjusted (per SD)", fontsize=6.2, loc="left")
fig.tight_layout(pad=0.2)
fig.savefig(os.path.join(FIGS, "F6b_forests.png"), dpi=600)
plt.close(fig)

# ---------------- Fig 8a: KM, final width ~45 mm (1.8 in)
d = pd.read_csv(os.path.join(ROOT, "analysis", "L0_bulk", "results", "track2", "TCGA_OV_residual.csv"))
d = d[(d.os_days > 0) & d.event.notna()].copy()
q = d.dec_res.quantile([0, .25, .75, 1]).values
d["grp"] = pd.cut(d.dec_res, q, labels=["Low", "Intermediate", "High"], include_lowest=True)
lr = multivariate_logrank_test(d.os_days, d.grp, d.event)
fig, ax = plt.subplots(figsize=(2.3, 1.95))
for g, col in [("Low", C["blue"]), ("Intermediate", C["grey"]), ("High", C["red"])]:
    s = d[d.grp == g]
    KaplanMeierFitter().fit(s.os_days / 365.25, s.event, label=f"{g} (n={len(s)})").plot_survival_function(
        ax=ax, ci_show=False, color=col, lw=0.9)
ax.set_xlabel("Years", fontsize=6.5); ax.set_ylabel("Overall survival", fontsize=6.2, labelpad=3.0)
ax.set_ylim(0, 1.02); ax.set_xlim(0, 12)
ax.tick_params(labelsize=6)
ax.text(0.97, 0.96, f"P = {lr.p_value:.3f}", transform=ax.transAxes, ha="right", va="top", fontsize=6)
ax.legend(frameon=False, fontsize=4.4, loc="upper center", bbox_to_anchor=(0.5, -0.40), ncol=3, handlelength=1.2, columnspacing=0.8)
ax.set_title("TCGA-OV (n = 423)", fontsize=6.2, loc="left")
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
fig.subplots_adjust(left=0.24, right=0.98, bottom=0.42, top=0.90)
fig.savefig(os.path.join(FIGS, "F8a_km.png"), dpi=600)
plt.close(fig)

# ---------------- Fig 8b: P-value bars, final width ~48 mm (1.9 in)
items = [("OS raw", 0.096), ("OS adj", 0.060), ("Resp", 0.412),
         ("Resp adj", 0.923), ("ER/PR", 0.699), ("ER/PR mal", 0.394)]
fig, ax = plt.subplots(figsize=(2.2, 1.2))
x = np.arange(len(items))
ax.bar(x, [-np.log10(p) for _, p in items], color=C["grey"], width=0.66)
ax.axhline(-np.log10(0.05), color=C["red"], lw=0.5, ls="--")
for i, (_, p) in enumerate(items):
    ax.text(i, -np.log10(p) + 0.02, f"{p:.3f}", ha="center", fontsize=5)
ax.set_xticks(x); ax.set_xticklabels([l for l, _ in items], fontsize=5.2, rotation=45, ha="right")
ax.set_ylabel("-log10(P)", fontsize=6.2)
ax.tick_params(axis="y", labelsize=5.5)
ax.set_ylim(0, 1.35)
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
fig.tight_layout(pad=0.2)
fig.savefig(os.path.join(FIGS, "F8b_layers.png"), dpi=600)
plt.close(fig)
print("panels regenerated at final size")
