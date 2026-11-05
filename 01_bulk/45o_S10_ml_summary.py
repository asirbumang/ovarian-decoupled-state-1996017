# -*- coding: utf-8 -*-
"""45o_S10_ml_summary.py: S10 三线 ML 汇总图（患者级 vs 细胞级 vs Scissor 重叠）。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs\supp"

fig, axes = plt.subplots(1, 2, figsize=(150/25.4, 54/25.4))
## a: 患者级 C-index（Mime 训练/验证 vs 45c CV vs Dummy）
lab = ["Mime train\n(TCGA)", "Mime validation\n(GSE26193)", "45c PCA30+Cox\n(TCGA 5-fold)", "Dummy"]
vals = [0.902, 0.588, 0.591, 0.502]
cols = ["#3C5488", "#E64B35", "#3C5488", "#B0B0B0"]
axes[0].bar(range(4), vals, color=cols, width=0.55)
axes[0].axhline(0.5, color="grey", linestyle=":", linewidth=0.6)
axes[0].set_xticks(range(4)); axes[0].set_xticklabels(lab, fontsize=5.4)
axes[0].set_ylabel("C-index", size=6); axes[0].tick_params(labelsize=5.4)
axes[0].set_ylim(0.4, 1.0)
for i, v in enumerate(vals):
    axes[0].text(i, v + 0.015, f"{v:.3f}", ha="center", fontsize=5.4)
axes[0].set_title("a patient-level: C-index collapses at ~0.59", size=6)
## b: 细胞级判别 + Scissor（不同尺度对比：AUC vs OR）
lab2 = ["cell-level zoo\n(best SVM)", "54-gene signature\n(LR)", "Scissor+ vs state\n(OR, p=0.17)"]
vals2 = [0.978, 0.997, 0.87]
cols2 = ["#3C5488", "#3C5488", "#B0B0B0"]
axes[1].bar(range(3), vals2, color=cols2, width=0.55)
axes[1].axhline(1.0, color="grey", linestyle=":", linewidth=0.6)
axes[1].set_xticks(range(3)); axes[1].set_xticklabels(lab2, fontsize=5.4)
axes[1].set_ylabel("AUC / OR", size=6); axes[1].tick_params(labelsize=5.4)
axes[1].set_ylim(0.5, 1.1)
for i, v in enumerate(vals2):
    axes[1].text(i, v + 0.015, f"{v:.3f}", ha="center", fontsize=5.4)
axes[1].set_title("b cell-level discriminable; no survival-cell overlap", size=6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "S10_ml_summary.png"), dpi=600)
plt.close(fig)
print("S10 saved", flush=True)
print("DONE 45o")
