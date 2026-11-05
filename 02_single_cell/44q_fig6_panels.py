# -*- coding: utf-8 -*-
"""44q_fig6_panels.py: 图6 验证+阴性报告——a 效应方向矩阵 + b 双森林(raw vs adjusted) + c crude/resid 对照。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"

## a: 15 结果方向×显著矩阵
rows = [
    ("Metastasis (epi, 27a)", -1, "cell p~0; sample-level NS"),
    ("Persister C vs T (27b)", +1, "sample-level p=0.70 NS"),
    ("Stage (epi within, 26b)", -1, "within-epi lower post-chemo"),
    ("Sites (27a)", -1, "cell-level sig; sample NS"),
    ("Spatial Moran (28h)", +1, "14/14 p<0.05"),
    ("Neighbourhood OR (28d)", +1, "7/12 p<0.05"),
    ("COMMOT axis (28f2)", +1, "CCL2/TIM3/CCL5/HGF perm p<0.05"),
    ("NNLS plasma delta", +1, "1/3 CI excl 0"),
    ("GeoMx compartment", -1, "Tumor vs Immune p=3.5e-15"),
    ("TLS/CD20 (28e3)", 0, "NS"),
    ("ER vs PR (26e)", 0, "NS both"),
    ("Drug screen (CTRPv2)", -1, "277/544 q<0.25"),
    ("CNV-malignant bridge", -1, "AUC 0.377 negative"),
    ("OS adjusted TCGA", +1, "HR 1.17 p=0.025 (exploratory)"),
    ("OS external (dec_res)", 0, "all NS; pooled p=0.71"),
]
fig, ax = plt.subplots(figsize=(105/25.4, 78/25.4))
Y = np.arange(len(rows))
for i, (lab, d, note) in enumerate(rows):
    col = {"-1": "#DC0000", "0": "#B0B0B0", "1": "#3C5488"}[str(d)]
    ax.barh(i, 0.35, left=0, color=col, height=0.55)
    ax.text(0.45, i, f"{lab}: {note}", va="center", ha="left", fontsize=4.6)
ax.set_xlim(0, 100); ax.axis("off")
ax.invert_yaxis()
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in ("#DC0000", "#B0B0B0", "#3C5488")]
ax.legend(handles, ["opposite", "no change/NS", "positive/consistent"], fontsize=4.8, loc="lower left", ncol=3, frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F6a_direction.png"), dpi=600)
plt.close(fig)
print("F6a saved", flush=True)

## b: 双森林（raw tier HR vs dec_res cont HR）
raw = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L0_bulk\results\survival_meta\survival_meta_results.csv")
adj = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L0_bulk\results\track2\decreg_cox_meta.csv")
fig, axes = plt.subplots(1, 2, figsize=(110/25.4, 54/25.4), sharey=True)
for ax, d, title, hr, lo, hi, pv in [
        (axes[0], raw.iloc[::-1], "raw decouple (tertile)", "hr_tier", "ci_lo", "ci_hi", "p_trend"),
        (axes[1], adj.iloc[::-1], "adjusted dec_res (continuous)", "hr_cont", "ci_lo", "ci_hi", "p_cont")]:
    ys = np.arange(len(d))
    ax.errorbar(d[hr], ys, xerr=[d[hr] - d[lo], d[hi] - d[hr]], fmt="o", ms=3.5, color="#3C5488",
                ecolor="#3C5488", elinewidth=0.6, capsize=1.5)
    ax.axvline(1, color="grey", linestyle=":", linewidth=0.8)
    ax.set_yticks(ys); ax.set_yticklabels(d["cohort"], fontsize=5)
    ax.set_xscale("log"); ax.set_xlim(0.45, 2.6)
    ax.set_title(title, size=6); ax.tick_params(labelsize=5)
    for i, (_, r) in enumerate(d.iterrows()):
        star = "*" if r[pv] < 0.05 else ""
        ax.text(2.7, i, f"p={r[pv]:.3f}{star}", fontsize=4.5, va="center")
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F6b_forests.png"), dpi=600)
plt.close(fig)
print("F6b saved", flush=True)

## c: crude vs residual 对照（TCGA OS / GSE156699 反应）
crd = pd.DataFrame({"metric": ["TCGA OS (3-tier)", "GSE156699 Resp vs NonResp"],
                    "crude": [0.096, 0.412], "resid": [0.060, 0.923]})
fig, ax = plt.subplots(figsize=(62/25.4, 40/25.4))
x = np.arange(2)
w = 0.3
ax.bar(x - w / 2, -np.log10(crd["crude"]), w, color="#B9B9B9", label="crude")
ax.bar(x + w / 2, -np.log10(crd["resid"]), w, color="#3C5488", label="composition-adjusted")
ax.axhline(-np.log10(0.05), color="#DC0000", linestyle="--", linewidth=0.5)
ax.text(1.45, -np.log10(0.05) + 0.02, "p=0.05", fontsize=4.5, color="#DC0000")
for i, (_, r) in enumerate(crd.iterrows()):
    ax.text(i - w / 2, -np.log10(r["crude"]) + 0.03, f"{r['crude']:.3f}", ha="center", fontsize=4.5)
    ax.text(i + w / 2, -np.log10(r["resid"]) + 0.03, f"{r['resid']:.3f}", ha="center", fontsize=4.5)
ax.set_xticks(x); ax.set_xticklabels(crd["metric"], fontsize=4.8, rotation=12)
ax.set_ylabel("-log10(p)", size=6); ax.tick_params(labelsize=5)
ax.legend(fontsize=4.8, frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F6c_crude_resid.png"), dpi=600)
plt.close(fig)
print("F6c saved", flush=True)
print("DONE 44q")
