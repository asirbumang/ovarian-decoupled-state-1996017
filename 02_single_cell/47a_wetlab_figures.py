# -*- coding: utf-8 -*-
"""47a_wetlab_figures.py: S12(W1: qPCR 热图/指数/CCK8/LDH) + S13(W2: qPCR 热图/IC50/CI)。"""
import numpy as np, pandas as pd, glob, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
B = "D:/wetlab_results_/"
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs\supp"

def load(gene, cell, w2=False):
    f = f"{B}qPCR_{gene}_{cell}_W2.csv" if w2 else f"{B}qPCR_{gene}_{cell}.csv"
    assert os.path.exists(f), f
    d = pd.read_csv(f)
    d["bio_fold"] = d.groupby(["group"] + (["time"] if "time" in d.columns else []))["fold_change"].transform(
        lambda x: x.mean() * 0 + x)  # noop
    return d

## 数据装载：W1 24h 每基因-细胞：组→bio 级（9 复 per group? 每 group 3 bio×3 tech → bio fold = tech mean）
def bio_fold(d, which="mean"):
    grp = ["group", "bio_rep", "time"] if "time" in d.columns else ["group", "bio_rep"]
    b = d.groupby(grp)["fold_change"].mean().reset_index()
    return b

genes9 = ["NFKBIA", "CXCL8", "IL1B", "IRF1", "SLC7A11", "GCLC", "GCLM", "BAX", "RHOA"]
## S12-a: 24h log2 fold heatmap（gene × [ctrl,bso,tnfa,bso_tnfa]（两系分块））
vals = {}
for g in genes9:
    for c in ["A2780", "SKOV3"]:
        d = pd.read_csv(f"{B}qPCR_{g}_{c}.csv")
        b = bio_fold(d)
        b24 = b[b.time == "24h"]
        m = b24.groupby("group")["fold_change"].mean()
        vals[(g, c)] = {k: np.log2(m.get(k, 1.0)) for k in ["ctrl", "bso", "tnfa", "bso_tnfa"]}
mat = np.array([[vals[(g, c)][k] for k in ["ctrl", "bso", "tnfa", "bso_tnfa"]] for c in ["A2780", "SKOV3"] for g in genes9])
fig = plt.figure(figsize=(105/25.4, 84/25.4))
gs = fig.add_gridspec(2, 3, height_ratios=[1.35, 1.0], hspace=0.55)
ax = fig.add_subplot(gs[0, :])
im = ax.imshow(mat, cmap="RdBu_r", vmin=-1, vmax=5, aspect="auto")
ax.set_xticks(range(4)); ax.set_xticklabels(["ctrl", "BSO", "TNFa", "BSO+TNFa"], fontsize=5.4)
ax.set_yticks(range(18)); ax.set_yticklabels([f"{g.upper()[:7]}·{c[:4]}" for c in ["A2780", "SKOV3"] for g in genes9], fontsize=5.4)
for i in range(18):
    for j in range(4):
        ax.text(j, i, "%.1f" % mat[i, j], ha="center", va="center", fontsize=5.4,
                color="white" if abs(mat[i, j]) > 2 else "black")
plt.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label="log2 fold (24h)")
ax.set_title("a W1 qPCR: joint induction (log2, 24h)", size=5.8)

## S12-b: decouple-like 指数（每组）
idx = {}
for c in ["A2780", "SKOV3"]:
    for gset, nm in [(["NFKBIA", "CXCL8", "IL1B", "IRF1"], "M1"), (["SLC7A11", "GCLC", "GCLM"], "M2"), (["BAX", "RHOA"], "M3")]:
        for k in ["ctrl", "bso", "tnfa", "bso_tnfa"]:
            idx[(c, k, nm)] = np.mean([vals[(g, c)][k] for g in gset])
axb = fig.add_subplot(gs[1, 0])
x = np.arange(4); w = 0.38
for i, (c, col) in enumerate([("A2780", "#DC0000"), ("SKOV3", "#3C5488")]):
    vals_idx = [idx[(c, k, "M1")] - (idx[(c, k, "M2")] + idx[(c, k, "M3")]) / 2 for k in ["ctrl", "bso", "tnfa", "bso_tnfa"]]
    axb.bar(x + (i - 0.5) * w, vals_idx, w, color=col, label=c)
axb.set_xticks(x); axb.set_xticklabels(["ctrl", "BSO", "TNFa", "B+T"], fontsize=5.4)
axb.set_ylabel("decouple-like index (log2)", size=6); axb.tick_params(labelsize=5.4)
axb.legend(fontsize=5.4, frameon=False)
axb.set_title("b decouple index (M1-(M2+M3)/2)", size=5.8)

## S12-c: CCK8（两组曲线合并放在 c；LDH 放 d）
axc = fig.add_subplot(gs[1, 1])
for c, col in [("A2780", "#DC0000"), ("SKOV3", "#3C5488")]:
    d = pd.read_csv(f"{B}CCK8_{c}_W1.csv")
    for gr, ls in [("ctrl", "-"), ("bso_tnfa", "-.")]:
        s = d[d.group == gr].groupby("time_h")["viability_pct"].mean()
        axc.plot(s.index, s.values, ls, color=col, linewidth=0.9, label=f"{c} {gr}")
axc.set_xlabel("hours", size=6); axc.set_ylabel("viability (%)", size=6); axc.tick_params(labelsize=5.4)
axc.set_ylim(0, 115)
axc.legend(fontsize=5.4, frameon=False, ncol=2, loc="lower left")
axc.set_title("c W1 CCK8 (ctrl vs B+T)", size=5.8)

## S12-d: LDH
axd = fig.add_subplot(gs[1, 2])
for i, c in enumerate(["A2780", "SKOV3"]):
    d = pd.read_csv(f"{B}LDH_{c}.csv")
    gm = d.groupby("group")["ldh_pct"].mean()
    axd.bar(np.arange(4) + i * 0.38, [gm.get(k, 0) for k in ["ctrl", "bso", "tnfa", "bso_tnfa"]], 0.36,
            color="#DC0000" if c == "A2780" else "#3C5488", label=c)
axd.set_xticks(np.arange(4) + 0.19); axd.set_xticklabels(["ctrl", "BSO", "TNFa", "B+T"], fontsize=5.4)
axd.set_ylabel("LDH release (%)", size=6); axd.tick_params(labelsize=5.4)
axd.legend(fontsize=5.4, frameon=False); axd.set_title("d W1 LDH (24h)", size=5.8)
fig.savefig(os.path.join(QRY, "S12_w1.png"), dpi=600)
plt.close(fig)
print("S12 saved", flush=True)

## ---- S13: W2 ----
genes4 = ["NFKBIA", "CXCL8", "SLC7A11", "BAX"]
fig = plt.figure(figsize=(110/25.4, 66/25.4))
gs2 = fig.add_gridspec(1, 3, wspace=0.5)
ax1 = fig.add_subplot(gs2[0, 0])
mat2 = []
lab2 = []
for g in genes4:
    for c in ["A2780", "SKOV3"]:
        d = pd.read_csv(f"{B}qPCR_{g}_{c}_W2.csv")
        b = d.groupby(["group", "bio_rep"])["fold_change"].median().reset_index()
        m = b.groupby("group")["fold_change"].mean()
        mat2.append([np.log2(m.get(k, 1.0)) for k in ["ctrl", "ind", "bos", "tbso" if False else "jsh", "ind_bos", "ind_jsh", "ind_bos_jsh"] if k in m.index])
        lab2.append(f"{g} ({c})")
ax1.imshow(np.array(mat2), cmap="RdBu_r", vmin=-1, vmax=4, aspect="auto")
ax1.set_xticks(range(len(["ind", "bos", "jsh", "ind_bos", "ind_jsh", "ind_bos_jsh"])))
ax1.set_xticklabels(["ind", "bos", "jsh", "ind_bos", "ind_jsh", "ind_bos_jsh"], fontsize=5.4, rotation=30, ha="right")
ax1.set_yticks(range(len(lab2))); ax1.set_yticklabels(lab2, fontsize=5.4)
ax1.set_title("a W2 qPCR (log2 vs ctrl)", size=5.8)
ax2 = fig.add_subplot(gs2[0, 1])
for c, col in [("A2780", "#DC0000"), ("SKOV3", "#3C5488")]:
    d = pd.read_csv(f"{B}CCK8_{c}_W2.csv")
    for ind, ls in [(0, "-"), (1, "--")]:
        s = d[(d.induced == ind) & (d.jsh_uM == 0)].groupby("bos_uM")["viability_pct"].mean()
        ax2.plot(s.index, s.values, ls, color=col, linewidth=0.9, label=f"{c} ind={ind}")
ax2.set_xlabel("bosutinib (uM)", size=6); ax2.set_ylabel("viability (%)", size=6); ax2.tick_params(labelsize=5.4)
ax2.legend(fontsize=5.4, frameon=False); ax2.set_title("b W2 bos IC50 (no JSH)", size=5.8)
ax3 = fig.add_subplot(gs2[0, 2])
ax3.axis("off")
ax3.text(0.02, 0.92, "W2 notes:\n- 4 genes qPCR (NFKBIA/CXCL8/SLC7A11/BAX)\n- IC50/CI: next step (per-condition)\n- induced groups: 24-25% viability", fontsize=5.4, va="top")
fig.savefig(os.path.join(QRY, "S13_w2.png"), dpi=600)
plt.close(fig)
print("S13 saved", flush=True)
print("DONE 47a")
