# -*- coding: utf-8 -*-
"""Regenerate F7f (evidence matrix) and F8f (residual disease) with corrected labels."""
import numpy as np, pandas as pd, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"

## ---------------- F7f: 6-gene evidence matrix ----------------
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
ax.set_xticks(range(4))
ax.set_xticklabels(["TF\nknockout", "54-gene\nsignature", "Drug\ntarget", "COMMOT\naxis"],
                   fontsize=5.0, rotation=0, ha="center")
ax.set_yticks(range(6)); ax.set_yticklabels(genes_k, fontsize=5.4)
for i in range(6):
    for j in range(4):
        ax.text(j, i, "%.1f" % tile[i, j], ha="center", va="center", fontsize=5.4,
                color="white" if tile[i, j] > 0.6 else "black")
ax.set_title("6-gene evidence matrix (four lines)", size=5.8)
fig.tight_layout(); fig.savefig(os.path.join(QRY, "F7f_evidence.png"), dpi=600); plt.close(fig)
print("F7f saved", flush=True)

## ---------------- F8f: residual disease ----------------
d = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L0_bulk\results\track2\tcga_subtype_consensus.csv",
                low_memory=False)
d = d[d["subtype"].notna()]
tab = pd.crosstab(d["subtype"], d["residual_tumor"]).reindex(["high", "low"])
cats = ["No Macroscopic disease", "1-10 mm", "11-20 mm", ">20 mm"]
den = tab[cats].sum(axis=1)            # known residual disease only (R0 excluded)
norm = tab[cats].div(den, axis=0)
print("denominators:", den.to_dict(), flush=True)
print(">20 mm:", (tab[">20 mm"] / den).round(4).to_dict(), flush=True)
fig, ax = plt.subplots(figsize=(62/25.4, 42/25.4))
bot = np.zeros(2)
cols = ["#00A087", "#4DBBD5", "#3C5488", "#E64B35"]
for i, c in enumerate(cats):
    v = norm[c].values
    ax.bar([0, 1], v, bottom=bot, color=cols[i], width=0.55, label=c)
    bot += v
ax.set_xticks([0, 1])
ax.set_xticklabels(["high\n(n=%d)" % int(den["high"]),
                    "low\n(n=%d)" % int(den["low"])], fontsize=5.4)
ax.set_ylabel("fraction", size=6); ax.tick_params(labelsize=5.4)
ax.legend(fontsize=5.4, frameon=False, loc="upper left", bbox_to_anchor=(1.05, 1.0))
ax.set_title("residual tumour by subtype (P = 0.026)", size=5.8)
fig.tight_layout(); fig.savefig(os.path.join(QRY, "F8f_residual.png"), dpi=600); plt.close(fig)
print("F8f saved", flush=True)
