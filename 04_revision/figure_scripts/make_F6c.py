# -*- coding: utf-8 -*-
"""Clean replacement for Fig. 6c (crude vs composition-adjusted P values)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "Arial", "font.size": 9})
fig, ax = plt.subplots(figsize=(1464/400, 944/400), dpi=400)
labels = ["TCGA OS\n(3-tier)", "GSE156699\nresponse"]
crude = [-np.log10(0.096), -np.log10(0.412)]
adj = [-np.log10(0.060), -np.log10(0.923)]
x = np.arange(2); w = 0.34
b1 = ax.bar(x - w/2, crude, w, color="#B0B0B0", label="crude")
b2 = ax.bar(x + w/2, adj, w, color="#3C5488", label="composition-adjusted")
for bars, vals in [(b1, [0.096, 0.412]), (b2, [0.060, 0.923])]:
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.03,
                f"{v:.3f}", ha="center", va="bottom", fontsize=7.5)
ax.axhline(-np.log10(0.05), color="#E64B35", lw=0.8, ls="--")
ax.text(1.48, -np.log10(0.05) + 0.02, "P = 0.05", color="#E64B35", fontsize=7, ha="right")
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("-log10(P)", fontsize=9)
ax.set_ylim(0, 1.35)
ax.legend(frameon=False, fontsize=7, loc="upper right")
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs\F6c_crude_resid.png", dpi=400)
print("saved")
