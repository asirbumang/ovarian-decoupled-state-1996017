# -*- coding: utf-8 -*-
"""44h_fig3_v4.py: 图3 v4 —— a(染色体×细胞 CNV 谱热图) + b(含 AUC) + c(含 n/p) + d(CNV×decouple 桥接)。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
QRY = "D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs"

## a: 染色体×细胞 CNV 谱（按 malignant 排序）
prof = pd.read_csv(os.path.join(QRY, "F3a_chr_profile.csv.gz"), index_col=0)
cnv = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\percell_cnv_score_v2.csv")
mal_map = dict(zip(cnv["cell"], cnv["malignant_like"]))
cells = prof.columns.tolist()
order = sorted(cells, key=lambda c: (not mal_map.get(c, False), c))
prof2 = prof[order]
flag = np.array([mal_map.get(c, False) for c in order])
fig, ax = plt.subplots(figsize=(95/25.4, 60/25.4))
## 每 chr 中位居中（相对全部细胞中位；显示"恶性 vs 非恶"的染色体级偏差）
P = prof2.values.copy()
P = P - np.median(P, axis=1, keepdims=True)
# 再按行稳健缩放（除以行 MAD）避免 chr 间量纲差
mad = np.median(np.abs(P - np.median(P, axis=1, keepdims=True)), axis=1, keepdims=True)
P = P / (mad + 1e-6)
im = ax.imshow(P, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3,
               extent=[0, len(cells), 22, 0], interpolation="nearest")
ax.set_yticks(np.arange(1, 23)); ax.set_yticklabels([f"chr{i}" for i in range(1, 23)], fontsize=4)
ax.set_xticks([])
ax.axvline(np.sum(~flag), color="white", linewidth=1.2)
ax.text(len(cells)*0.5, -0.8, f"non-mal-like ({np.sum(~flag)}) | malignant-like ({np.sum(flag)})",
        ha="center", fontsize=5.5)
cb = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label="chr-median-centered CNV (MAD)")
cb.ax.tick_params(labelsize=5)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F3a_chr_heatmap.png"), dpi=600)
plt.close(fig)
print("F3a saved", prof2.shape, flush=True)

## b: 更新（AUC）
sc_e = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577_epithelial.csv",
                   usecols=["cell", "decouple"], low_memory=False)
cnv_e = cnv[cnv["group"] == "Epithelial"]
mrg = pd.merge(cnv_e, sc_e, on="cell")
auc = roc_auc_score(mrg["malignant_like"], mrg["decouple"])
print("AUC(decouple | malignant-like):", round(auc, 3), flush=True)
fig, ax = plt.subplots(figsize=(62/25.4, 54/25.4))
bins = np.linspace(cnv_e["cnv_top"].min(), cnv_e["cnv_top"].max(), 40)
ax.hist(cnv_e.loc[~cnv_e["malignant_like"], "cnv_top"], bins=bins, color="#BBBBBB", alpha=0.8,
        label="non-malignant-like")
ax.hist(cnv_e.loc[cnv_e["malignant_like"], "cnv_top"], bins=bins, color="#DC0000", alpha=0.8,
        label="malignant-like (n=%d)" % int(cnv_e["malignant_like"].sum()))
q99 = cnv_e.loc[~cnv_e["malignant_like"], "cnv_top"].quantile(.99)
ax.axvline(q99, color="#3C5488", linestyle="--", linewidth=0.5, label="Ref q99=%.3f" % q99)
ax.set_xlabel("CNV score (top-3000 genes |dev|)", size=6)
ax.set_ylabel("cells", size=6)
ax.tick_params(labelsize=5.5)
ax.legend(fontsize=5, frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F3b_cnv_dist.png"), dpi=600)
plt.close(fig)
print("F3b saved; AUC:", round(auc, 3), flush=True)

## d: CNV×decouple 散点（桥接；负向=恶性样解耦更低）
mrg2 = mrg.dropna(subset=["decouple"])
from scipy.stats import spearmanr, mannwhitneyu
rho, pp = spearmanr(mrg2["cnv_top"], mrg2["decouple"])
print("spearman(cnv, decouple):", round(rho, 3), "p:", pp, flush=True)
fig, ax = plt.subplots(figsize=(62/25.4, 54/25.4))
ax.scatter(mrg2.loc[~mrg2["malignant_like"], "cnv_top"], mrg2.loc[~mrg2["malignant_like"], "decouple"],
           s=0.3, c="#BBBBBB", alpha=0.3, linewidths=0)
ax.scatter(mrg2.loc[mrg2["malignant_like"], "cnv_top"], mrg2.loc[mrg2["malignant_like"], "decouple"],
           s=0.5, c="#DC0000", alpha=0.6, linewidths=0)
ax.set_xlabel("CNV score", size=6); ax.set_ylabel("decouple", size=6)
ax.tick_params(labelsize=5.5)
ax.text(0.03, 0.97, f"rho={rho:.2f} (p={pp:.1e})", transform=ax.transAxes, fontsize=5.5, va="top")
mwu = mannwhitneyu(mrg2.loc[mrg2["malignant_like"], "decouple"], mrg2.loc[~mrg2["malignant_like"], "decouple"])
ax.text(0.03, 0.90, f"MWU p={mwu.pvalue:.1e}", transform=ax.transAxes, fontsize=5.5, va="top")
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F3d_cnv_decouple.png"), dpi=600)
plt.close(fig)
print("F3d saved", flush=True)

## c: n/p 补充
from scipy.stats import mannwhitneyu as mwu2
d1 = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\GSE222556_type_annotated.csv",
                 usecols=["decouple", "type"], low_memory=False)
d1 = d1[d1["type"].isin(["Primary", "Metastatic"]) & d1["decouple"].notna()]
p1 = mwu2(d1.loc[d1["type"] == "Metastatic", "decouple"], d1.loc[d1["type"] == "Primary", "decouple"]).pvalue
d2 = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE293459.csv",
                 usecols=["sample_type", "decouple"], low_memory=False)
d2 = d2[d2["sample_type"].isin(["C", "T"]) & d2["decouple"].notna()]
p2 = mwu2(d2.loc[d2["sample_type"] == "T", "decouple"], d2.loc[d2["sample_type"] == "C", "decouple"]).pvalue
d3 = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577.csv",
                 usecols=["treatment_stage", "decouple"], low_memory=False)
d3 = d3[d3["treatment_stage"].isin(["chemo-naive", "IDS"]) & d3["decouple"].notna()]
p3 = mwu2(d3.loc[d3["treatment_stage"] == "IDS", "decouple"], d3.loc[d3["treatment_stage"] == "chemo-naive", "decouple"]).pvalue
print("Wilcoxon p:", "P-M %.1e (n=%d/%d)" % (p1, (d1['type']=='Primary').sum(), (d1['type']=='Metastatic').sum()),
      "| C-T %.1e (n=%d/%d)" % (p2, (d2['sample_type']=='C').sum(), (d2['sample_type']=='T').sum()),
      "| naive-IDS %.1e (n=%d/%d)" % (p3, (d3['treatment_stage']=='chemo-naive').sum(), (d3['treatment_stage']=='IDS').sum()), flush=True)
print("DONE 44h")
