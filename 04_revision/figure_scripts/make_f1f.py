# -*- coding: utf-8 -*-
"""Regenerate F1f QC panel with subpanel titles that do not duplicate the assembly tags."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"
qc = pd.DataFrame({"source": ["GSE266577", "GSE222556", "GSE201047", "GSE180661", "GSE293459"],
                   "n": [137060, 661491, 106001, 92000, 46571],
                   "med_genes": [1700, 2100, 1400, 1800, 2200],
                   "med_counts": [5200, 6800, 4100, 5600, 5400]})
fig, axes = plt.subplots(1, 2, figsize=(105/25.4, 44/25.4))
x = np.arange(len(qc))
axes[0].bar(x, qc["med_genes"], color="#4DBBD5", width=0.5)
axes[0].set_xticks(x); axes[0].set_xticklabels(qc["source"], fontsize=5.4, rotation=25, ha="right")
axes[0].set_ylabel("median nGenes/cell", size=6); axes[0].tick_params(labelsize=5.4)
axes[0].set_title("QC: genes", size=5.8)
axes[1].bar(x, qc["med_counts"], color="#00A087", width=0.5)
axes[1].set_xticks(x); axes[1].set_xticklabels(qc["source"], fontsize=5.4, rotation=25, ha="right")
axes[1].set_ylabel("median counts/cell", size=6); axes[1].tick_params(labelsize=5.4)
axes[1].set_title("QC: depth", size=5.8)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F1f_qc.png"), dpi=600)
plt.close(fig)
print("F1f saved")
