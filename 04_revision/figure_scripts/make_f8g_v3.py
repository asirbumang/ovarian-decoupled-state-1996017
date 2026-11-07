# -*- coding: utf-8 -*-
"""Fig 8g v3: 9-gene x (2 cell lines x 4 groups) heatmap + decouple index + LDH 24 h.
CCK-8 curves remain in Fig. S12, which keeps this panel readable at final width."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "svg.fonttype": "none"})

WET = r"D:\wetlab_results_"
AUD = r"D:\子宫内膜异位症2\_audit_r3"
GENES = [("NFKBIA", "M1"), ("CXCL8", "M1"), ("IL1B", "M1"), ("IRF1", "M1"),
         ("SLC7A11", "M2"), ("GCLC", "M2"), ("GCLM", "M2"),
         ("BAX", "M3"), ("RHOA", "M3")]
GROUPS = ["ctrl", "bso", "tnfa", "bso_tnfa"]
GLAB = ["ctrl", "BSO", "TNFα", "B+T"]

fc = {}
for gname, _ in GENES:
    for cell in ["A2780", "SKOV3"]:
        d = pd.read_csv(os.path.join(WET, f"qPCR_{gname}_{cell}.csv"))
        sub = d[d.time == "24h"]
        fc[(gname, cell)] = [sub[sub.group == g].fold_change.mean() for g in GROUPS]

mat = np.log2(np.array([[fc[(g, c)][i] for i in range(4)] for c in ["A2780", "SKOV3"] for g, _ in GENES]))
mat = mat.reshape(2, 9, 4).transpose(1, 0, 2).reshape(9, 8)

def geo(vals):
    return float(np.exp(np.mean(np.log(np.array(vals)))))

idx = {}
for cell in ["A2780", "SKOV3"]:
    vals = []
    for i in range(4):
        m1 = geo([fc[(g, cell)][i] for g, m in GENES if m == "M1"])
        m2 = geo([fc[(g, cell)][i] for g, m in GENES if m == "M2"])
        m3 = geo([fc[(g, cell)][i] for g, m in GENES if m == "M3"])
        vals.append(np.log2(m1) - (np.log2(m2) + np.log2(m3)) / 2)
    idx[cell] = vals

FW, FH = 1300 / 600.0, 1290 / 600.0
fig = plt.figure(figsize=(FW, FH))
gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1.0], hspace=0.80, wspace=0.45,
                      left=0.295, right=0.935, top=0.86, bottom=0.12)
fig.text(0.5, 0.965, "a W1 qPCR: joint induction (log2, 24 h)",
         ha="center", va="top", fontsize=6.2)

axa = fig.add_subplot(gs[0, :])
im = axa.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=5)
axa.set_yticks(range(9))
axa.set_yticklabels([g for g, _ in GENES], fontsize=5.4)
axa.set_xticks(range(8))
axa.set_xticklabels(GLAB * 2, fontsize=5.4, rotation=30, ha="right")
axa.axvline(3.5, color="white", lw=3)
axa.text(1.5, -1.25, "A2780", ha="center", va="bottom", fontsize=5.6, fontweight="bold")
axa.text(5.5, -1.25, "SKOV3", ha="center", va="bottom", fontsize=5.6, fontweight="bold")
for i in range(9):
    for j in range(8):
        axa.text(j, i, f"{mat[i, j]:.1f}", ha="center", va="center", fontsize=5.4,
                 color="white" if mat[i, j] > 3 else "#222222")
cb = fig.colorbar(im, ax=axa, fraction=0.018, pad=0.006)
cb.ax.tick_params(labelsize=5.4)
cb.set_label("log2 fold", fontsize=5.4, labelpad=2)
cb.outline.set_linewidth(0.4)

axb = fig.add_subplot(gs[1, 0])
x = np.arange(4); w = 0.36
axb.bar(x - w / 2, idx["A2780"], w, color="#DC0000", label="A2780")
axb.bar(x + w / 2, idx["SKOV3"], w, color="#3C5488", label="SKOV3")
axb.set_xticks(x); axb.set_xticklabels(["ctrl", "BSO", "TNFα", "B+T"], fontsize=5.4, rotation=30, ha="right")
axb.tick_params(labelsize=5.4)
axb.set_ylabel("decouple index (log2)", fontsize=5.4, labelpad=1)
axb.set_title("b decouple index", fontsize=5.4, pad=2)
axb.axhline(0, color="k", lw=0.4)
axb.set_ylim(-1.0, 4.15)
for s in ("top", "right"):
    axb.spines[s].set_visible(False)

axd = fig.add_subplot(gs[1, 1])
for i, cell in enumerate(["A2780", "SKOV3"]):
    d = pd.read_csv(os.path.join(WET, f"LDH_{cell}.csv"))
    d = d[d.time_h == 24]
    gm = d.groupby("group")["ldh_pct"].mean()
    axd.bar(np.arange(4) + i * 0.38, [gm.get(k, 0) for k in GROUPS], 0.36,
            color="#DC0000" if cell == "A2780" else "#3C5488", label=cell)
axd.set_xticks(np.arange(4) + 0.19)
axd.set_xticklabels(["ctrl", "BSO", "TNFα", "B+T"], fontsize=5.4, rotation=30, ha="right")
axd.tick_params(labelsize=5.4)
axd.set_ylabel("LDH release (%)", fontsize=5.4, labelpad=1)
axd.set_title("c LDH at 24 h", fontsize=5.4, pad=2)
for s in ("top", "right"):
    axd.spines[s].set_visible(False)

pg = os.path.join(AUD, "new_fig8g_w1_v3.png")
fig.savefig(pg, dpi=600); plt.close(fig)
print("saved", pg)
