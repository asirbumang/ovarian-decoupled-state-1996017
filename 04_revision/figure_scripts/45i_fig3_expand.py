# -*- coding: utf-8 -*-
"""45i_fig3_expand.py: 图3 补 d(恶性 vs 非恶 三模块) + e(各谱系 CNV 分数)。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"

cnv = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\percell_cnv_score_v2.csv")
cnv = cnv[cnv["group"] == "Epithelial"].copy()
sc = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577_epithelial.csv",
                 usecols=["cell", "M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"], low_memory=False)
for m in ["M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"]:
    sc[m] = pd.to_numeric(sc[m], errors="coerce").astype(float)
mrg = cnv.merge(sc, on="cell")
fig, axes = plt.subplots(1, 3, figsize=(90/25.4, 46/25.4))
for ax, m, t in zip(axes, ["M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"],
                    ["M1 immune", "M2 metab/GSH", "M3 mito-cyto"]):
    a = mrg.loc[mrg["malignant_like"], m].values
    b = mrg.loc[~mrg["malignant_like"], m].values
    vp = ax.violinplot([b, a], positions=[1, 2], showmedians=True, widths=0.7)
    for pc in vp["bodies"]:
        pc.set_alpha(0.5)
    p = mannwhitneyu(a, b).pvalue
    pstr = "P < 0.001" if p < 0.001 else f"P = {p:.3f}"
    ax.set_xticks([1, 2]); ax.set_xticklabels(["non-mal", "malignant"], fontsize=5.4)
    ax.set_title(f"{t} ({pstr})", size=5.2); ax.tick_params(labelsize=5.4)
    ax.set_xlim(0.4, 2.6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F3d_modules_mal.png"), dpi=600)
plt.close(fig)
print("F3d saved; n mal/non:", int(mrg['malignant_like'].sum()), int((~mrg['malignant_like']).sum()), flush=True)

## e: 各谱系 CNV 分数箱线（42,719 全有 CNV：上皮 22,719 + 参考 20,000）
cnv_e = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\percell_cnv_score_v2.csv", low_memory=False)
scf = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577.csv",
                  usecols=["cell", "cell_type"], low_memory=False)
LIN = {"Epithelial cells":"Epi","Fibroblasts":"Fibro","Mesothelial cells":"Meso","Endothelial cells":"Endo",
 "Macrophages":"Myeloid","Classical monocytes":"Myeloid","Migratory DCs":"DC","DC1":"DC","DC2":"DC","pDC":"DC",
 "Mast cells":"Mast","Tem/Trm cytotoxic T cells":"T","Tcm/Naive helper T cells":"T","Regulatory T cells":"T",
 "Type 17 helper T cells":"T","NK cells":"NK","CD16- NK cells":"NK","CD16+ NK cells":"NK","ILC":"NK",
 "Memory B cells":"B","Naive B cells":"B","Plasma cells":"Plasma"}
scf["lin"] = scf["cell_type"].map(LIN)
m2 = cnv_e.merge(scf[["cell", "lin"]], on="cell")
order = ["Epi", "Fibro", "Meso", "Endo", "Myeloid", "DC", "Mast", "T", "NK", "B", "Plasma"]
order = [l for l in order if l in set(m2["lin"].dropna())]
fig, ax = plt.subplots(figsize=(80/25.4, 50/25.4))
data = [m2.loc[m2["lin"] == l, "cnv_top"] for l in order]
bp = ax.boxplot(data, tick_labels=order, widths=0.55, patch_artist=True, showfliers=False)
for p in bp["boxes"]:
    p.set_facecolor("#4DBBD5"); p.set_alpha(0.7)
ax.set_ylabel("CNV score", size=6); ax.tick_params(labelsize=5.4, labelrotation=30)
ax.set_title("CNV score by lineage (epithelial highest)", size=6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F3e_cnv_lineage.png"), dpi=600)
plt.close(fig)
print("F3e saved", flush=True)
print("DONE 45i")
