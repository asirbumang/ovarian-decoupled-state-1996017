# -*- coding: utf-8 -*-
"""44a_fig1_panels.py: 图1 b(UMAP 总览) + c(marker DotPlot) 面板（scanpy 默认 + Set3；33 决议）。"""
import numpy as np, pandas as pd, sys, os
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})

d = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\umap_input"
umo = pd.read_csv(os.path.join(d, "umap_full_266577.csv"))
sc = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577.csv",
                 usecols=["cell", "cell_type"])
LIN = {"Epithelial cells":"Epi","Fibroblasts":"Fibro","Mesothelial cells":"Meso","Endothelial cells":"Endo",
 "Macrophages":"MonoMac","Classical monocytes":"MonoMac","Migratory DCs":"DC","DC1":"DC","DC2":"DC","pDC":"DC",
 "Mast cells":"Mast","Tem/Trm cytotoxic T cells":"T","Tcm/Naive helper T cells":"T","Regulatory T cells":"T",
 "Type 17 helper T cells":"T","NK cells":"NK","CD16- NK cells":"NK","CD16+ NK cells":"NK","ILC":"NK",
 "Memory B cells":"B","Naive B cells":"B","Plasma cells":"Plasma"}
sc["lin"] = sc["cell_type"].map(LIN)
umo = umo.merge(sc[["cell", "lin"]], on="cell", how="left")
umo["lin"] = umo["lin"].fillna("NA")
NEWS = {"Epi":"#E64B35","Fibro":"#4DBBD5","Meso":"#00A087","Endo":"#3C5488","MonoMac":"#F39B7F",
        "DC":"#8491B4","Mast":"#91D1C2","T":"#DC0000","NK":"#7E6148","B":"#B09C85","Plasma":"#A73030","NA":"#E5E5E5"}
OL = ["Epi","Fibro","Meso","Endo","MonoMac","DC","Mast","T","NK","B","Plasma","NA"]

## ---- b: UMAP 总览（全量 137,084；Set3 等效 NPG 低饱和；size 0.25/alpha 0.75；灰底雾）----
fig, ax = plt.subplots(figsize=(78/25.4, 72/25.4))
pts = umo[["umap1", "umap2"]].values
ax.scatter(pts[:, 0], pts[:, 1], s=0.03, c="#BBBBBB", alpha=0.12, linewidths=0)
for lin in OL:
    s = umo[umo["lin"] == lin]
    if len(s) == 0:
        continue
    ax.scatter(s["umap1"], s["umap2"], s=0.25, c=NEWS[lin], alpha=0.8, linewidths=0, label=lin)
ax.set_xticks([]); ax.set_yticks([])
for sp_ in ax.spines.values():
    sp_.set_visible(False)
handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=5.5,
                frameon=False, markerfirst=True)
for h in leg.legend_handles:
    h.set_sizes([8])
fig.tight_layout()
fig.savefig(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs\F1b_umap.png", dpi=600)
plt.close(fig)
print("F1b umap saved; cells with lin:", int(umo["lin"].notna().sum()), flush=True)

## ---- c: DotPlot（12 谱系 × 36 marker；点大小=pct_in、颜色=mean log1p(CP10k)）----
MARK = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\annotation_review\dotplot_markers_12lineage.csv")
print(MARK.groupby("lineage").size(), flush=True)
genes = [l.strip() for l in open(os.path.join(d, "GSE266577_full.genes.tsv"), encoding="utf-8")]
cells = [l.strip() for l in open(os.path.join(d, "GSE266577_full.cells.tsv"), encoding="utf-8")]
dt = pd.read_csv(os.path.join(d, "GSE266577_full.triplet.gz"))
gidx = {g: i for i, g in enumerate(genes)}
sidx = {c: i for i, c in enumerate(cells)}
C = sp.csr_matrix((dt["v"].to_numpy(), (dt["g"].map(gidx).to_numpy(), dt["c"].map(sidx).to_numpy())),
                  shape=(len(genes), len(cells)))
del dt
cell_lin = dict(zip(sc["cell"], sc["lin"]))
# 每谱系每 marker：pct_in（表达>0 比例）、mean log1p(CP10k)
cg = {l: [i for i, c in enumerate(cells) if cell_lin.get(c) == l] for l in set(sc["lin"].dropna())}
cpm = C.multiply(1e4 / C.sum(axis=0).A1).tocsr()
rows = []
for _, r in MARK.iterrows():
    gi = gidx[r["gene"]]
    for lin, idxs in cg.items():
        if not idxs:
            continue
        sub = cpm[gi, idxs]
        v = np.asarray(sub.todense()).ravel()
        rows.append({"lin": lin, "gene": r["gene"], "pct": float((v > 0).mean()),
                     "mean": float(np.log1p(v).mean())})
dk = pd.DataFrame(rows)
# 顺序
order_l = ["Epi", "Fibro", "Meso", "Endo", "MonoMac", "DC", "Mast", "T", "NK", "B", "Plasma"]
order_l = [l for l in order_l if l in set(dk["lin"])]
mk_order = dk.head(0)
genes_order = []
for l in order_l:
    genes_order += list(MARK[MARK["lineage"] == l]["gene"])
genes_order = list(dict.fromkeys(genes_order))
fig, ax = plt.subplots(figsize=(88/25.4, 60/25.4))
X = np.zeros((len(order_l), len(genes_order))); SZ = np.zeros_like(X)
for i, l in enumerate(order_l):
    for j, g in enumerate(genes_order):
        r = dk[(dk["lin"] == l) & (dk["gene"] == g)]
        if len(r):
            X[i, j] = r["mean"].iloc[0]; SZ[i, j] = r["pct"].iloc[0]
## per-gene z-score（标准 DotPlot 口径）
Xz = X.copy()
for j in range(X.shape[1]):
    col = X[:, j]
    if col.std() > 1e-9:
        Xz[:, j] = (col - col.mean()) / col.std()
    else:
        Xz[:, j] = 0
ax.scatter(np.tile(np.arange(len(genes_order)), len(order_l)),
           np.repeat(np.arange(len(order_l)), len(genes_order)),
           s=(SZ.ravel() * 90 + 2), c=Xz.ravel(), cmap="RdBu_r", vmin=-2.0, vmax=2.0, linewidths=0)
tick_idx = list(range(0, len(genes_order), 2))
ax.set_xticks(tick_idx); ax.set_xticklabels([genes_order[i] for i in tick_idx], rotation=45, ha="right", fontsize=4.5)
ax.set_yticks(range(len(order_l))); ax.set_yticklabels(order_l, fontsize=5.5)
cb = plt.colorbar(ax.collections[0], ax=ax, fraction=0.025, pad=0.01, label="per-gene z (mean log1p)")
cb.ax.tick_params(labelsize=5)
fig.tight_layout()
fig.savefig(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs\F1c_dotplot.png", dpi=600)
plt.close(fig)
print("F1c dotplot saved", flush=True)
print("DONE 44a")
