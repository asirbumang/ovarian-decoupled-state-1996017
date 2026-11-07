# -*- coding: utf-8 -*-
"""Regenerate only F7e, F7f and F8f (avoids the slow F3g loop in 46c)."""
import numpy as np, pandas as pd, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"

## F7e: alpha grid
ag = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b9\b9_alpha_grid.csv")
ag = ag[ag["sig"] == "g200"]
fig, ax = plt.subplots(figsize=(70/25.4, 46/25.4))
for tf, c in [("NFKB1", "#E64B35"), ("IRF1", "#3C5488"), ("STAT1", "#00A087")]:
    sub = ag[ag["tf"] == tf].groupby("alpha")["KO_score"].mean()
    ax.plot(range(3), sub.values, marker="o", ms=3, color=c, label=tf, linewidth=0.9)
ax.set_xticks(range(3)); ax.set_xticklabels(["0.1", "1", "10"], fontsize=5.4)
ax.set_xlabel("Ridge alpha", size=6); ax.set_ylabel("KO score", size=6); ax.tick_params(labelsize=5.4)
ax.legend(fontsize=5.4, frameon=False)
ax.set_title("Ridge alpha grid (3 hosts merged)", size=5.8)
fig.tight_layout(); fig.savefig(os.path.join(QRY, "F7e_alpha.png"), dpi=600); plt.close(fig)
print("F7e saved", flush=True)

## F7f: 6-gene evidence matrix
genes_k = ["IRF1", "NFKBIA", "CXCL8", "LDHA", "RHOA", "ACTG1"]
sig_f = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b12\sig_freq.csv")
freq = dict(zip(sig_f["gene"], sig_f["freq"]))
tile = np.zeros((6, 4))
b9 = {"IRF1": 1, "NFKBIA": 1}
for i, g in enumerate(genes_k):
    tile[i, 0] = b9.get(g, 0)
    tile[i, 1] = freq.get(g, 0)
fig, ax = plt.subplots(figsize=(72/25.4, 46/25.4))
ax.imshow(tile, cmap="Reds", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(4)); ax.set_xticklabels(["B9 KO*", "ML freq", "drug target", "COMMOT axis"],
                                            fontsize=5.4, rotation=20, ha="right")
ax.set_yticks(range(6)); ax.set_yticklabels(genes_k, fontsize=5.4)
for i in range(6):
    for j in range(4):
        ax.text(j, i, "%.1f" % tile[i, j], ha="center", va="center", fontsize=5.4,
                color="white" if tile[i, j] > 0.6 else "black")
ax.set_title("6-gene evidence matrix (convergence)", size=5.8)
fig.tight_layout(); fig.savefig(os.path.join(QRY, "F7f_evidence.png"), dpi=600); plt.close(fig)
print("F7f saved", flush=True)

## F8f: residual disease stacked bar
d = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L0_bulk\results\track2\tcga_subtype_consensus.csv", low_memory=False)
d = d[d["subtype"].notna()]
tab = pd.crosstab(d["subtype"], d["residual_tumor"]).reindex(["high", "low"])
cats = ["No Macroscopic disease", "1-10 mm", "11-20 mm", ">20 mm"]
norm = tab[cats].div(tab[cats].sum(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(62/25.4, 42/25.4))
bot = np.zeros(2)
cols = ["#00A087", "#4DBBD5", "#3C5488", "#E64B35"]
for i, c in enumerate(cats):
    v = norm[c].values
    ax.bar([0, 1], v, bottom=bot, color=cols[i], width=0.55, label=c)
    bot += v
ax.set_xticks([0, 1])
ax.set_xticklabels(["high\n(n=%d)" % len(d[d.subtype == "high"]),
                    "low\n(n=%d)" % len(d[d.subtype == "low"])], fontsize=5.4)
ax.set_ylabel("fraction", size=6); ax.tick_params(labelsize=5.4)
ax.legend(fontsize=5.4, frameon=False, loc="upper left", bbox_to_anchor=(1.05, 1.0))
ax.set_title("residual tumour by subtype (χ² P = 0.025)", size=5.8)
fig.tight_layout(); fig.savefig(os.path.join(QRY, "F8f_residual.png"), dpi=600); plt.close(fig)
print("F8f saved", flush=True)
