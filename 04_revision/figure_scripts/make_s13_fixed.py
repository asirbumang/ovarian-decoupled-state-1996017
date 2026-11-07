# -*- coding: utf-8 -*-
"""Fig 8h (S13): W2 qPCR heatmap + normalized dose-response, two panels, no clipped text box."""
import numpy as np, pandas as pd, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})
B = "D:/更改/"
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs\supp"
genes4 = ["NFKBIA", "CXCL8", "SLC7A11", "BAX"]
fig = plt.figure(figsize=(125/25.4, 66/25.4))
gs = fig.add_gridspec(1, 2, wspace=0.38, width_ratios=[1.05, 1.0])

# a: qPCR heatmap
mat2, lab2 = [], []
for g in genes4:
    for c in ["A2780", "SKOV3"]:
        d = pd.read_csv(f"{B}qPCR_{g}_{c}_W2.csv")
        m = d.groupby("group")["fold_change"].mean()
        mat2.append([np.log2(m.get(k, 1.0)) for k in ["ind", "ind_bos", "ind_jsh", "ind_bos_jsh"]])
        lab2.append(f"{g}-{c[:4]}")
ax1 = fig.add_subplot(gs[0, 0])
im = ax1.imshow(np.array(mat2), cmap="RdBu_r", vmin=-1, vmax=5.5, aspect="auto")
ax1.set_xticks(range(4)); ax1.set_xticklabels(["ind", "+bos", "+jsh", "bos+jsh"], fontsize=5.4, rotation=25, ha="right")
ax1.set_yticks(range(len(lab2))); ax1.set_yticklabels(lab2, fontsize=5.4)
ax1.set_title("a W2 qPCR (log2; JSH blocks NF-κB)", size=5.8)
cb = fig.colorbar(im, ax=ax1, fraction=0.035, pad=0.02)
cb.ax.tick_params(labelsize=5.0)

# b: normalized dose-response
ax2 = fig.add_subplot(gs[0, 1])
for c, col in [("A2780", "#DC0000"), ("SKOV3", "#3C5488")]:
    d = pd.read_csv(f"{B}CCK8_{c}_W2.csv")
    p = d[(d.time_h == 48) & (d.jsh_uM == 0)]
    for ind, ls in [(0, "-"), (1, "--")]:
        s = p[p.induced == ind].groupby("bos_uM")["viability_norm_pct"].mean()
        ax2.plot(s.index, s.values, ls, color=col, linewidth=1.0, marker="o", ms=2.5,
                 label=f"{c} ind={ind}")
ax2.set_xlabel("bosutinib (µM)", size=6); ax2.set_ylabel("viability norm (%)", size=6)
ax2.tick_params(labelsize=5.4); ax2.set_ylim(40, 110)
ax2.legend(fontsize=5.2, frameon=False, loc="lower left")
ax2.set_title("b Dose-response (48 h)", size=5.8)

fig.savefig(os.path.join(QRY, "S13_w2.png"), dpi=600)
plt.close(fig)
print("S13 fixed saved", flush=True)
