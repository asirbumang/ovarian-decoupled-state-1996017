# -*- coding: utf-8 -*-
"""44k_fig4_cd.py: 图4 c(Moran I 点图) + d(邻域 OR)（3 层定量）。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
OUT = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial"
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"

m = pd.read_csv(os.path.join(OUT, "spatial_moran_summary.csv"))
n = pd.read_csv(os.path.join(OUT, "spatial_neighborhood_summary.csv"))
fig, axes = plt.subplots(1, 2, figsize=(120/25.4, 52/25.4))
m["lab"] = [f"{r.GSE[:7]} {r.sample}" for _, r in m.iterrows()]
order = m.sort_values("Moran_I")["lab"].tolist()
axes[0].scatter(m["Moran_I"], range(len(m)), s=8, c="#3C5488")
axes[0].set_yticks(range(len(m))); axes[0].set_yticklabels(order, fontsize=4.5)
axes[0].set_xlabel("Moran's I (decouple, k=6)", size=6); axes[0].tick_params(labelsize=5.5)
axes[0].axvline(m["Moran_I"].median(), color="#E64B35", linestyle="--", linewidth=0.5)
axes[0].text(m["Moran_I"].median() + 0.01, len(m) - 1.2, "median=%.3f" % m["Moran_I"].median(), fontsize=5, color="#E64B35")
axes[0].set_title("Spatial autocorrelation (14/14 p<0.05)", size=6)

n = n[n["fisher_p"].notna()]
n["logOR"] = np.log(n["neigh_OR"])
n["lab"] = [f"{r.GSE[:7]} {r.sample}" for _, r in n.iterrows()]
order2 = n.sort_values("neigh_OR")["lab"].tolist()
axes[1].scatter(n["neigh_OR"], range(len(n)), s=8, c=("#DC0000" if p2 < 0.05 else "#B0B0B0") if False else "#3C5488")
lab2 = [r["neigh_OR"] if r["fisher_p"] < 0.05 else None for _, r in n.iterrows()]
for i, (_, r) in enumerate(n.iterrows()):
    if r["fisher_p"] < 0.05:
        axes[1].scatter(r["neigh_OR"], i, s=10, c="#DC0000")
axes[1].axvline(1, color="grey", linestyle=":", linewidth=0.5)
axes[1].set_yticks(range(len(n))); axes[1].set_yticklabels(order2, fontsize=4.5)
axes[1].set_xlabel("neighbourhood OR (state-high)", size=6); axes[1].tick_params(labelsize=5.5)
axes[1].set_title("Neighbourhood enrichment (red p<0.05)", size=6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F4cd_spatial_stats.png"), dpi=600)
plt.close(fig)
print("F4cd saved", flush=True)
print("DONE 44k")
