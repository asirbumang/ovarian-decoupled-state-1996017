# -*- coding: utf-8 -*-
"""Generate new supplementary figures for the revision (S31-S33, S36-S37)."""
import json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = r"D:\子宫内膜异位症2"
RES = os.path.join(ROOT, "R2_revision", "04_analysis", "results")
FIG = os.path.join(ROOT, "R2_revision", "06_figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "pdf.fonttype": 42, "ps.fonttype": 42})
C = {"blue": "#3C5488", "red": "#E64B35", "green": "#00A087", "grey": "#8491B4",
     "orange": "#F39B7F", "teal": "#4DBBD5"}

# ---------------------------------------------------------------- S31 null models
a1 = json.load(open(os.path.join(RES, "A1_null_models.json")))
a1b = json.load(open(os.path.join(RES, "A1b_loglinear.json")))
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw={"width_ratios": [1.1, 1]})
labels = ["Observed", "Independent\n(0.25³)", "Permutation\n(marginals)",
          "Gaussian\ncopula", "Pairwise\nlog-linear"]
vals = [a1["state_pct"], a1["independence_expected_pct"],
        a1["permutation_null"]["null_pct_mean"], a1["copula_null"]["expected_pct"],
        a1b["models"][0]["expected_ABC_pairwise"] / a1["n_cells"] * 100]
cols = [C["red"], C["grey"], C["grey"], C["teal"], C["blue"]]
bars = ax[0].bar(range(5), vals, color=cols, width=0.65, edgecolor="black", linewidth=0.4)
ax[0].axhline(a1["state_pct"], color=C["red"], lw=0.8, ls="--")
ax[0].set_xticks(range(5)); ax[0].set_xticklabels(labels, fontsize=6.5)
ax[0].set_ylabel("State frequency (% of cells)")
for i, v in enumerate(vals):
    ax[0].text(i, v + 0.03, f"{v:.2f}", ha="center", va="bottom", fontsize=6.5)
ax[0].set_ylim(0, 1.85)
ax[0].set_title("Frequency vs null models (observed/independent = %.2f)" % a1["observed_over_independence"], fontsize=8, loc="left")


mods = a1b["models"]
names = ["Pairwise\nonly", "+ cell type", "+ patient", "+ both"]
beta = [m["ABC_coef"] for m in mods]
lo = [m["ABC_coef"] - 1.96 * m["ABC_se"] for m in mods]
hi = [m["ABC_coef"] + 1.96 * m["ABC_se"] for m in mods]
ax[1].axhline(0, color="black", lw=0.6)
ax[1].errorbar(range(len(mods)), beta, yerr=[np.array(beta) - lo, np.array(hi) - np.array(beta)],
               fmt="o", color=C["blue"], ecolor=C["blue"], capsize=2, ms=4, lw=1)
for i, m in enumerate(mods):
    ax[1].text(i, hi[i] + 0.015, f"O/E {m['OE']:.2f}", ha="center", fontsize=6, clip_on=False)
ax[1].set_xticks(range(len(mods))); ax[1].set_xticklabels(names, fontsize=6.5)
ax[1].set_ylabel("Three-way interaction (log scale)")
ax[1].set_ylim(-0.05, 0.66)
ax[1].set_title("Excess of the joint state over pairwise structure", fontsize=8, loc="left")
for s in ["top", "right"]:
    ax[0].spines[s].set_visible(False); ax[1].spines[s].set_visible(False)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(FIG, f"S31_null_models.{ext}"), dpi=400, bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------- S32 threshold sensitivity
th = pd.read_csv(os.path.join(RES, "A2_threshold_grid.csv"))
alt = pd.read_csv(os.path.join(RES, "A2_alternative_definitions.csv"))
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
piv = th.pivot(index="M1_pct", columns="M2M3_pct", values="pct_cells")
im = ax[0].imshow(piv.values, cmap="Reds", aspect="auto", vmin=0)
ax[0].set_xticks(range(len(piv.columns))); ax[0].set_xticklabels(piv.columns)
ax[0].set_yticks(range(len(piv.index))); ax[0].set_yticklabels(piv.index)
ax[0].set_xlabel("M2/M3 percentile cutoff"); ax[0].set_ylabel("M1 percentile cutoff")
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        ax[0].text(j, i, f"{piv.values[i, j]:.2f}", ha="center", va="center", fontsize=6,
                   color="white" if piv.values[i, j] > piv.values.max() * 0.6 else "black")
cb = fig.colorbar(im, ax=ax[0], fraction=0.046, pad=0.02)
cb.set_label("% state cells", fontsize=7); cb.ax.tick_params(labelsize=6)
ax[0].set_title("Cell frequency across cutoffs", fontsize=8, loc="left")

for key, col, lab in [("plasma_OR", C["red"], "Plasma"), ("fibroblast_OR", C["blue"], "Fibroblast"),
                      ("epithelial_OR", C["grey"], "Epithelial")]:
    ax[1].plot(range(len(th)), th[key], "-o", ms=2.5, lw=0.9, color=col, label=lab)
ax[1].axhline(1, color="black", lw=0.6, ls="--")
ax[1].set_xticks(range(len(th))); ax[1].set_xticklabels([f"{a}/{b}" for a, b in zip(th.M1_pct, th.M2M3_pct)],
                                                       rotation=90, fontsize=5)
ax[1].set_ylabel("State OR vs other lineages")
ax[1].legend(frameon=False, fontsize=6, loc="upper right")
ax[1].set_title("Carrier enrichment across cutoffs (25 settings)", fontsize=8, loc="left")
for s in ["top", "right"]:
    ax[1].spines[s].set_visible(False)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(FIG, f"S32_threshold_sensitivity.{ext}"), dpi=400, bbox_inches="tight")
plt.close(fig)

# ----------------------------------------------------------- S33 technical effects
a3 = json.load(open(os.path.join(RES, "A3_technical.json")))
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8))
q = a3["state_rate_by_nFeature_quintile"]
ax[0].bar(range(len(q)), [x["mean"] * 100 for x in q], color=C["blue"], edgecolor="black", lw=0.4)
ax[0].set_xticks(range(len(q))); ax[0].set_xticklabels([f"Q{i+1}" for i in range(len(q))])
ax[0].set_xlabel("Detected-gene quintile (low → high)"); ax[0].set_ylabel("State frequency (%)")
ax[0].set_title("State frequency falls with library complexity", fontsize=8, loc="left")
for i, x in enumerate(q):
    ax[0].text(i, x["mean"] * 100 + 0.05, f"{x['mean']*100:.2f}", ha="center", fontsize=6)

sp = a3["module_spearman"]
x = np.arange(3); w = 0.25
for i, (k, lab) in enumerate([("M1z", "M1"), ("M2z", "M2"), ("M3z", "M3")]):
    vals = [sp[v][k]["rho"] for v in ["nCount_RNA", "nFeature_RNA", "percent.mt"]]
    ax[1].bar(x + (i - 1) * w, vals, width=w, label=lab, color=[C["red"], C["green"], C["blue"]][i])
ax[1].axhline(0, color="black", lw=0.6)
ax[1].set_xticks(x); ax[1].set_xticklabels(["nCount", "nFeature", "% mito"])
ax[1].set_ylabel("Spearman ρ with module z")
ax[1].legend(frameon=False, fontsize=6, ncol=3, loc="upper right")
ax[1].set_title("Module scores are depth-dependent", fontsize=8, loc="left")
for s in ["top", "right"]:
    ax[0].spines[s].set_visible(False); ax[1].spines[s].set_visible(False)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(FIG, f"S33_technical_effects.{ext}"), dpi=400, bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------ S36 pseudotime roots
a7 = json.load(open(os.path.join(RES, "A7_pseudotime.json")))
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
for lin, col in [("fib", C["red"]), ("plasma", C["blue"])]:
    for root, ls in [("min", "-"), ("max", "--")]:
        if root in a7[lin]:
            bins = pd.DataFrame(a7[lin][root]["bins"])
            ax[0].plot(bins["bin"], bins["frac"] * 100, ls, color=col, lw=1.1,
                       marker="o", ms=2.5, label=f"{lin} ({root}-root)")
ax[0].set_xlabel("Pseudotime bin"); ax[0].set_ylabel("State frequency (%)")
ax[0].legend(frameon=False, fontsize=5.5, loc="upper left")
ax[0].set_title("Root dependence of state-along-pseudotime", fontsize=8, loc="left")

rows = []
for lin, col in [("fib", C["red"]), ("plasma", C["blue"])]:
    for root in ["min", "max"]:
        if root in a7[lin] and a7[lin][root]["per_patient"]:
            for r in a7[lin][root]["per_patient"]:
                rows.append({"lin": lin, "root": root, "rho": r["rho"]})
rr = pd.DataFrame(rows)
pos = {"fib": 0, "plasma": 1}
for i, (lin, col) in enumerate([("fib", C["red"]), ("plasma", C["blue"])]):
    sub = rr[(rr.lin == lin) & (rr.root == "min")]
    ax[1].scatter(np.full(len(sub), i) + np.random.default_rng(1).normal(0, 0.04, len(sub)),
                  sub.rho, s=8, color=col, alpha=0.7, linewidth=0)
ax[1].axhline(0, color="black", lw=0.6)
ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(["Fibroblast", "Plasma"])
ax[1].set_ylabel("Per-patient Spearman ρ (state vs pseudotime)")
ax[1].set_title("No consistent within-patient direction", fontsize=8, loc="left")
for s in ["top", "right"]:
    ax[0].spines[s].set_visible(False); ax[1].spines[s].set_visible(False)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(FIG, f"S36_pseudotime_roots.{ext}"), dpi=400, bbox_inches="tight")
plt.close(fig)

# ----------------------------------------------------------------- S37 forest plot
a5 = json.load(open(os.path.join(RES, "A5_clinical.json")))
co = a5["meta_external"]["cohorts"]
fig, ax = plt.subplots(figsize=(4.6, 2.8))
ys = np.arange(len(co))[::-1]
for y, c in zip(ys, co):
    ax.plot([c["ci_lo"], c["ci_hi"]], [y, y], color="black", lw=0.9)
    ax.plot(c["hr_cont"], y, "s", color=C["blue"], ms=4)
    ax.text(2.35, y, f"{c['hr_cont']:.2f} ({c['ci_lo']:.2f}-{c['ci_hi']:.2f})", va="center", fontsize=6)
y0 = -1.2
m = a5["meta_external"]
ax.plot([m["RE_CI"][0], m["RE_CI"][1]], [y0, y0], color=C["red"], lw=1.2)
ax.plot(m["RE_HR"], y0, "D", color=C["red"], ms=5)
ax.text(2.35, y0, f"{m['RE_HR']:.2f} ({m['RE_CI'][0]:.2f}-{m['RE_CI'][1]:.2f})", va="center", fontsize=6,
        color=C["red"])
ax.axvline(1, color="black", lw=0.6, ls="--")
ax.set_yticks(list(ys) + [y0]); ax.set_yticklabels([c["cohort"] for c in co] + ["Random effects"], fontsize=6.5)
ax.set_xscale("log"); ax.set_xlim(0.4, 2.4)
ax.set_xlabel("HR per 1 SD (95% CI)")
ax.set_title(f"External cohorts: I² = {m['I2_pct']:.0f}%, τ² = {m['tau2']:.3f}", fontsize=8, loc="left")
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(FIG, f"S37_forest_random_effects.{ext}"), dpi=400, bbox_inches="tight")
plt.close(fig)
print("figures written to", FIG)
