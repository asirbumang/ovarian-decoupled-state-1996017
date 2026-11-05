# -*- coding: utf-8 -*-
"""44e_fig3_panels.py: 图3 b(CNV 分数分布) + c(横向三组：P-M / C-T / chemo-IDS)。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
QRY = "D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs"

## --- b: per-cell CNV 分数分布（恶性 vs 非恶 上皮）---
cnv = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\percell_cnv_score_v2.csv")
cnv = cnv[cnv["group"] == "Epithelial"].copy()
bins = np.linspace(cnv["cnv_top"].min(), cnv["cnv_top"].max(), 40)
fig, ax = plt.subplots(figsize=(62/25.4, 54/25.4))
ax.hist(cnv.loc[~cnv["malignant_like"], "cnv_top"], bins=bins, color="#BBBBBB", alpha=0.8, label="non-malignant-like")
ax.hist(cnv.loc[cnv["malignant_like"], "cnv_top"], bins=bins, color="#DC0000", alpha=0.8, label="malignant-like (n=%d)" % int(cnv["malignant_like"].sum()))
ax.axvline(cnv.loc[~cnv["malignant_like"], "cnv_top"].quantile(.99), color="#3C5488", linestyle="--", linewidth=0.5,
           label="Ref q99")
ax.set_xlabel("CNV score (top-3000 genes |dev|)", size=6)
ax.set_ylabel("cells", size=6)
ax.tick_params(labelsize=5.5)
leg = ax.legend(fontsize=5, frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F3b_cnv_dist.png"), dpi=600)
plt.close(fig)
print("F3b saved; malignant:", int(cnv["malignant_like"].sum()), flush=True)

## --- c: 横向三组（decouple 箱线）---
def box_panel(ax, vals, labels, colors, ylab):
    bp = ax.boxplot(vals, tick_labels=labels, patch_artist=True, widths=0.55, showfliers=False,
                    medianprops=dict(color="black", linewidth=0.6),
                    boxprops=dict(linewidth=0.5), whiskerprops=dict(linewidth=0.5),
                    capprops=dict(linewidth=0.5))
    for p, c in zip(bp["boxes"], colors):
        p.set_facecolor(c); p.set_alpha(0.75)
    ax.set_ylabel(ylab, size=6)
    ax.tick_params(labelsize=5.5)

fig, axes = plt.subplots(1, 3, figsize=(150/25.4, 58/25.4), sharey=True)
## 1) GSE222556 Primary vs Metastatic（恶性上皮）
d1 = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\GSE222556_type_annotated.csv",
                 usecols=["decouple", "type"], low_memory=False)
print("222556 cols ok; n =", len(d1), flush=True)
d1 = d1[d1["type"].notna() & d1["type"].isin(["Primary", "Metastatic"])]
d1 = d1[d1["decouple"].notna()]
box_panel(axes[0], [d1.loc[d1["type"] == "Primary", "decouple"], d1.loc[d1["type"] == "Metastatic", "decouple"]],
          ["Primary", "Metastatic"], ["#4DBBD5", "#E64B35"], "decouple (z)")
axes[0].set_title("GSE222556 (metastasis)", size=6)
## 2) GSE293459 C vs T
d2 = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE293459.csv",
                 usecols=["sample", "decouple", "sample_type"], low_memory=False)
d2 = d2[d2["sample_type"].isin(["C", "T"]) & d2["decouple"].notna()]
box_panel(axes[1], [d2.loc[d2["sample_type"] == "C", "decouple"], d2.loc[d2["sample_type"] == "T", "decouple"]],
          ["C (n=%d)" % int((d2["sample_type"] == "C").sum()), "T"], ["#4DBBD5", "#E64B35"], "")
axes[1].set_title("GSE293459 (persister)", size=6)
## 3) GSE266577 chemo-naive vs IDS
d3 = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577.csv",
                 usecols=["decouple", "treatment_stage"], low_memory=False)
d3 = d3[d3["treatment_stage"].isin(["chemo-naive", "IDS"]) & d3["decouple"].notna()]
box_panel(axes[2], [d3.loc[d3["treatment_stage"] == "chemo-naive", "decouple"], d3.loc[d3["treatment_stage"] == "IDS", "decouple"]],
          ["chemo-naive", "IDS"], ["#4DBBD5", "#E64B35"], "")
axes[2].set_title("GSE266577 (stage)", size=6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F3c_horizontal.png"), dpi=600)
plt.close(fig)
print("F3c saved; n1:", len(d1), " n2:", len(d2), " n3:", len(d3), flush=True)
import pandas as pd
print("medians:", d1.groupby("type")["decouple"].median().round(3).to_dict(),
      d2.groupby("sample_type")["decouple"].median().round(3).to_dict(),
      d3.groupby("treatment_stage")["decouple"].median().round(3).to_dict(), flush=True)
print("DONE 44e")
