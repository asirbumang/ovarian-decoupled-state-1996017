# -*- coding: utf-8 -*-
"""Redraw main-figure panels flagged as hard to read (Fig. 6b, 8a, 8b)."""
import json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

ROOT = r"D:\子宫内膜异位症2"
RES = os.path.join(ROOT, "R2_revision", "04_analysis", "results")
FIG = os.path.join(ROOT, "R2_revision", "06_figures")
plt.rcParams.update({"font.family": "Arial", "font.size": 9, "axes.linewidth": 0.7,
                     "xtick.major.width": 0.7, "ytick.major.width": 0.7,
                     "pdf.fonttype": 42, "ps.fonttype": 42})
C = {"blue": "#3C5488", "red": "#E64B35", "green": "#00A087", "grey": "#8491B4"}

# ---------------------------------------------------------------- Fig 6b forest
meta = pd.read_csv(os.path.join(ROOT, "analysis", "L0_bulk", "results", "track2", "decreg_cox_meta.csv"))
a5 = json.load(open(os.path.join(RES, "A5_clinical.json")))
fig, ax = plt.subplots(figsize=(4.4, 2.6))
ys = np.arange(len(meta))[::-1]
for y, (_, c) in zip(ys, meta.iterrows()):
    ax.plot([c.ci_lo, c.ci_hi], [y, y], color="black", lw=1.0)
    ax.plot(c.hr_cont, y, "s", color=C["blue"], ms=5)
    ax.text(2.55, y, f"{c.hr_cont:.2f} ({c.ci_lo:.2f}–{c.ci_hi:.2f})", va="center", fontsize=8)
m = a5["meta_external"]
y0 = -1.3
ax.plot([m["RE_CI"][0], m["RE_CI"][1]], [y0, y0], color=C["red"], lw=1.4)
ax.plot(m["RE_HR"], y0, "D", color=C["red"], ms=6)
ax.text(2.55, y0, f"{m['RE_HR']:.2f} ({m['RE_CI'][0]:.2f}–{m['RE_CI'][1]:.2f})", va="center",
        fontsize=8, color=C["red"])
ax.axvline(1, color="black", lw=0.7, ls="--")
ax.set_yticks(list(ys) + [y0])
ax.set_yticklabels(list(meta.cohort) + ["External, pooled"], fontsize=8)
ax.set_xscale("log"); ax.set_xlim(0.45, 2.2)
ax.set_xticks([0.5, 0.75, 1, 1.5, 2]); ax.set_xticklabels(["0.5", "0.75", "1.0", "1.5", "2.0"])
import matplotlib.ticker as mticker
ax.xaxis.set_minor_locator(mticker.NullLocator())
ax.set_xlabel("Hazard ratio per 1 SD (95% CI)", fontsize=9)
ax.set_title("Overall survival, composition-adjusted score", fontsize=9, loc="left")
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(FIG, f"Fig6b_redrawn.{ext}"), dpi=400, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Fig 8a KM
d = pd.read_csv(os.path.join(ROOT, "analysis", "L0_bulk", "results", "track2", "TCGA_OV_residual.csv"))
d = d[(d.os_days > 0) & d.event.notna()].copy()
q = d.dec_res.quantile([0, .25, .75, 1]).values
d["grp"] = pd.cut(d.dec_res, q, labels=["Low", "Intermediate", "High"], include_lowest=True)
lr = multivariate_logrank_test(d.os_days, d.grp, d.event)
fig, ax = plt.subplots(figsize=(3.6, 3.1))
cols = {"Low": C["blue"], "Intermediate": C["grey"], "High": C["red"]}
for g in ["Low", "Intermediate", "High"]:
    s = d[d.grp == g]
    km = KaplanMeierFitter().fit(s.os_days / 365.25, s.event, label=f"{g} (n={len(s)})")
    km.plot_survival_function(ax=ax, ci_show=False, color=cols[g], lw=1.3)
ax.set_xlabel("Years", fontsize=9); ax.set_ylabel("Overall survival", fontsize=9)
ax.set_ylim(0, 1.02); ax.set_xlim(0, 12)
ax.text(0.98, 0.96, f"Log-rank P = {lr.p_value:.3f}", transform=ax.transAxes,
        ha="right", va="top", fontsize=9)
ax.legend(frameon=False, fontsize=7, loc="upper right", bbox_to_anchor=(1.0, 0.88))
ax.set_title("TCGA-OV (n = 423)", fontsize=9, loc="left")
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(FIG, f"Fig8a_redrawn.{ext}"), dpi=400, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Fig 8b p bars
items = [("OS, crude", 0.096), ("OS, adjusted", 0.060), ("Response, crude", 0.412),
         ("Response, adjusted", 0.923), ("ER vs PR, all spots", 0.70),
         ("ER vs PR, malignant", 0.39)]
fig, ax = plt.subplots(figsize=(3.8, 2.4))
y = np.arange(len(items))[::-1]
ax.barh(y, [-np.log10(p) for _, p in items], color=C["grey"], height=0.6)
ax.axvline(-np.log10(0.05), color=C["red"], lw=0.8, ls="--")
for yy, (lab, p) in zip(y, items):
    ax.text(-np.log10(p) + 0.02, yy, f"P = {p:.3f}".replace("P = 0.060", "P = 0.060"),
            va="center", fontsize=7)
ax.set_yticks(y); ax.set_yticklabels([l for l, _ in items], fontsize=7.5)
ax.set_xlabel("−log10(P)", fontsize=9)
ax.set_xlim(0, 1.3)
ax.text(-np.log10(0.05), len(items) - 0.4, "P = 0.05", color=C["red"], fontsize=7, ha="center")
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(FIG, f"Fig8b_redrawn.{ext}"), dpi=400, bbox_inches="tight")
plt.close(fig)
print("main panels redrawn")
