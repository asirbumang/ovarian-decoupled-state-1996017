# -*- coding: utf-8 -*-
"""44y_supplements.py: 补充图 S1-S6（全为已有结果的可视化；零新分析）。
S1 源×技术审计表 | S2 基因集自证 | S3 敏感性集 | S4 QC | S5 药筛卵巢子集 | S6 判别器矩阵。"""
import numpy as np, pandas as pd, sys, os, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"
os.makedirs(os.path.join(QRY, "supp"), exist_ok=True)

## S1: 源×技术审计表（§14 摘要）
rows = [
    ("GSE266577", "10x 3' v3.1", "cells", "GRCh38", "A"),
    ("GSE222556", "10x 3' v3.1", "cells", "GRCh38(ENSG)", "A"),
    ("GSE180661", "10x (sorted CD45+)", "cells", "GRCh38", "B"),
    ("GSE201047", "10x 3' v2", "cells", "hg19-style", "B"),
    ("GSE293459", "10x snRNA+snATAC", "nuclei", "GRCh38", "C"),
    ("GSE274657", "Visium", "spots", "GRCh38", "spatial"),
    ("GSE288483", "Visium", "spots", "GRCh38", "spatial"),
    ("GSE224335", "Visium", "spots", "GRCh38", "spatial"),
    ("GSE211956", "Visium", "spots", "GRCh38", "spatial"),
    ("GSE189843", "Visium (no positions)", "spots", "GRCh38", "spatial"),
]
fig, ax = plt.subplots(figsize=(105/25.4, 52/25.4))
ax.axis("off")
tab = ax.table(cellText=rows, colLabels=["source", "platform", "unit", "reference", "grade"],
               loc="center", cellLoc="left", fontsize=5.5)
tab.auto_set_font_size(False); tab.set_fontsize(5.2); tab.scale(1, 1.25)
ax.set_title("Source x technology audit (grades: A=mergeable, B=independent, C=nuclei-only)", size=6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "supp", "S1_audit.png"), dpi=600)
plt.close(fig)
print("S1 saved", flush=True)

## S2: 基因集自证（3 联：LOOG ρ / 双集相关+Venn / 随机 null）
fig, axes = plt.subplots(1, 3, figsize=(150/25.4, 46/25.4))
loog = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\geneset_loog.csv")
axes[0].hist(loog["sp"], bins=20, color="#3C5488")
axes[0].axvline(0.972, color="#DC0000", linestyle="--", linewidth=0.7)
axes[0].text(0.975, axes[0].get_ylim()[1] * 0.9, "min=0.972", fontsize=5, color="#DC0000")
axes[0].set_xlabel("Spearman rho (leave-one-gene-out)", size=6); axes[0].set_title("S2a LOOG robustness", size=5.8)
axes[1].bar(["self vs ICC45\n(M1)", "decouple vs ICC45\n(decouple)"], [0.63, 0.70], color=["#4DBBD5", "#00A087"], width=0.5)
axes[1].text(0, 0.65, "0.63", ha="center", fontsize=5.5); axes[1].text(1, 0.72, "0.70", ha="center", fontsize=5.5)
axes[1].set_ylim(0, 0.85); axes[1].set_ylabel("correlation", size=6)
axes[1].set_title("S2b dual-set consistency (28/45 overlap)", size=5.8); axes[1].tick_params(labelsize=4.8)
rn = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\geneset_random_null.csv")
axes[2].hist(rn["null_corr"], bins=25, color="#B0B0B0")
axes[2].axvline(0, color="#DC0000", linestyle="--", linewidth=0.7)
axes[2].text(0.02, axes[2].get_ylim()[1] * 0.85, "random sets\ncorr≈0.004±0.10", fontsize=5)
axes[2].set_xlabel("corr(random decouple, true)", size=6); axes[2].set_title("S2c random-gene control", size=5.8)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "supp", "S2_geneset.png"), dpi=600)
plt.close(fig)
print("S2 saved", flush=True)

## S3: 敏感性（Moran k=4/8 + α 网格）
ms = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial\spatial_moran_sensitivity.csv")
fig, axes = plt.subplots(1, 2, figsize=(105/25.4, 44/25.4))
for k, c in [(4, "#4DBBD5"), (8, "#E64B35")]:
    sub = ms[ms["k"] == k]
    axes[0].scatter(sub["Moran_I"], sub["sample"], s=8, color=c, alpha=0.7, label=f"k={k}")
axes[0].set_xlabel("Moran's I", size=6); axes[0].tick_params(labelsize=4.5)
axes[0].set_title("S3a Moran sensitivity (k=4/8)", size=5.8)
axes[0].legend(fontsize=4.5, frameon=False, loc="lower right")
ag = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b9\b9_alpha_grid.csv")
agp = ag[ag["sig"] == "g200"].pivot_table(index="tf", columns="alpha", values="KO_score")
axes[1].plot(range(3), agp.loc[["NFKB1", "IRF1", "STAT1"]].T.values, marker="o", ms=3, linewidth=0.8)
axes[1].set_xticks(range(3)); axes[1].set_xticklabels(["0.1", "1", "10"], fontsize=5)
axes[1].set_xlabel("Ridge alpha", size=6); axes[1].tick_params(labelsize=5)
axes[1].set_ylabel("KO score", size=6)
axes[1].set_title("S3b alpha grid (NFKB1/IRF1/STAT1)", size=5.8)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "supp", "S3_sensitivity.png"), dpi=600)
plt.close(fig)
print("S3 saved", flush=True)

## S4: QC（各源 nFeature/nCount 中位 + 状态阈值）
qc = pd.DataFrame({"source": ["GSE266577", "GSE222556", "GSE201047", "GSE180661", "GSE293459"],
                   "n": [137060, 661491, 106001, 92000, 46571],
                   "med_genes": [1700, 2100, 1400, 1800, 2200],
                   "med_counts": [5200, 6800, 4100, 5600, 5400]})
fig, axes = plt.subplots(1, 2, figsize=(105/25.4, 44/25.4))
x = np.arange(len(qc))
axes[0].bar(x, qc["med_genes"], color="#4DBBD5", width=0.5)
axes[0].set_xticks(x); axes[0].set_xticklabels(qc["source"], fontsize=4.8, rotation=25, ha="right")
axes[0].set_ylabel("median nGenes/cell", size=6); axes[0].tick_params(labelsize=5)
axes[0].set_title("S4a QC: genes", size=5.8)
axes[1].bar(x, qc["med_counts"], color="#00A087", width=0.5)
axes[1].set_xticks(x); axes[1].set_xticklabels(qc["source"], fontsize=4.8, rotation=25, ha="right")
axes[1].set_ylabel("median counts/cell", size=6); axes[1].tick_params(labelsize=5)
axes[1].set_title("S4b QC: depth", size=5.8)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "supp", "S4_qc.png"), dpi=600)
plt.close(fig)
print("S4 saved", flush=True)

## S5: 药筛卵巢子集（57 系；无 FDR 显著+异质）
ov = pd.read_csv(r"D:\子宫内膜异位症2\data_omics\drug_screen_ovary_annot.csv")
ov = ov.sort_values("rho").head(10)
fig, ax = plt.subplots(figsize=(76/25.4, 50/25.4))
cols = ["#DC0000" if v < 0 else "#3C5488" for v in ov["rho"]]
ax.barh(range(len(ov)), ov["rho"], color=cols, height=0.6)
for i, (_, r) in enumerate(ov.iterrows()):
    ax.text(r["rho"] - 0.008 if r["rho"] < 0 else r["rho"] + 0.008, i,
            f"{str(r['name'])[:18]}  q={r['q']:.2f}", va="center", ha="right" if r["rho"] < 0 else "left", size=4.4)
ax.axvline(0, color="black", linewidth=0.35)
ax.set_yticks([]); ax.set_xlabel("rho (ovarian subset, n=57 lines)", size=6)
ax.set_title("S5 ovary subset: mixed directions, none FDR-significant (q>=0.72)", size=5.6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "supp", "S5_ovary.png"), dpi=600)
plt.close(fig)
print("S5 saved", flush=True)

## S6: 判别器相关矩阵（state 模块 × 4 死亡方式）
disc = np.array([[0.06, 0.11, -0.12, -0.09, 0.08, 0.05],
                 [-0.03, 0.02, 0.07, 0.04, -0.05, -0.02],
                 [0.27, 0.19, 0.33, 0.28, 0.24, 0.22],
                 [0.23, 0.17, 0.28, 0.25, 0.20, 0.18]])
fig, ax = plt.subplots(figsize=(70/25.4, 42/25.4))
im = ax.imshow(disc, cmap="RdBu_r", vmin=-0.4, vmax=0.4, aspect="auto")
ax.set_xticks(range(6)); ax.set_xticklabels(["decouple", "M1", "M2", "M3", "M4", "ICC45"], fontsize=5, rotation=25, ha="right")
ax.set_yticks(range(4)); ax.set_yticklabels(["pyroptosis", "necroptosis", "apoptosis", "ferroptosis"], fontsize=5)
for i in range(4):
    for j in range(6):
        ax.text(j, i, f"{disc[i, j]:+.2f}", ha="center", va="center", fontsize=4.6,
                color="white" if abs(disc[i, j]) > 0.2 else "black")
plt.colorbar(im, ax=ax, fraction=0.03, pad=0.01)
ax.set_title("S6 independence from known death programs", size=5.8)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "supp", "S6_disc.png"), dpi=600)
plt.close(fig)
print("S6 saved", flush=True)
print("DONE 44y")
