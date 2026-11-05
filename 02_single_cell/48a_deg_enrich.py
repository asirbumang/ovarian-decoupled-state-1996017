# -*- coding: utf-8 -*-
"""48a_deg_enrich.py: 补齐 DEG+火山+GSEA 三件套（state vs rest，GSE266577 全细胞抽样 22k）。
① 读全量 triplet → 子集（state 全部 + 随机非 state 20,000）② Wilcoxon DEG（scanpy，BH）
③ 火山图（top 标注）④ gseapy prerank（Hallmark；如可用）气泡图 → S14。"""
import numpy as np, pandas as pd, sys, os, warnings
warnings.filterwarnings("ignore")
import scipy.sparse as sp
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding="utf-8")
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
D = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\umap_input"
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"
OUTD = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\deg"
os.makedirs(OUTD, exist_ok=True)

genes = [l.strip() for l in open(os.path.join(D, "GSE266577_full.genes.tsv"), encoding="utf-8")]
cells = [l.strip() for l in open(os.path.join(D, "GSE266577_full.cells.tsv"), encoding="utf-8")]
dt = pd.read_csv(os.path.join(D, "GSE266577_full.triplet.gz"))
gidx = {g: i for i, g in enumerate(genes)}
sidx = {c: i for i, c in enumerate(cells)}
sc_meta = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577.csv",
                      usecols=["cell", "M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"], low_memory=False)
for m in ["M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"]:
    sc_meta[m] = pd.to_numeric(sc_meta[m], errors="coerce").astype(float)
z = lambda x: (x - x.mean()) / (x.std() + 1e-9)
z1, z2, z3 = z(sc_meta["M1_immune_activation"]), z(sc_meta["M2_metabolic_GSH"]), z(sc_meta["M3_mito_cytoskeleton"])
sc_meta["state"] = (z1 >= np.quantile(z1, .75)) & (z2 <= np.quantile(z2, .25)) & (z3 <= np.quantile(z3, .25))
state_cells = set(sc_meta.loc[sc_meta["state"], "cell"])
keep_state = [c for c in cells if c in state_cells]
rng = np.random.default_rng(42)
others = [c for c in cells if c not in state_cells]
keep_other = list(rng.choice(others, size=20000, replace=False))
keep = keep_state + keep_other
print("state:", len(keep_state), "other:", len(keep_other), flush=True)
sel = set(keep)
dt = dt[dt["c"].isin(sel)]
print("triplet subset rows:", len(dt), flush=True)
c2i = {c: i for i, c in enumerate(keep)}
C = sp.csr_matrix((dt["v"].to_numpy(), (dt["g"].map(gidx).to_numpy(), dt["c"].map(c2i).to_numpy())),
                  shape=(len(genes), len(keep)))
del dt
cpm = C.multiply(1e4 / C.sum(axis=0).A1).tocsr()
X = cpm.T.tocsr()
ad = sc.AnnData(X=X.astype(np.float32))
ad.var_names = genes
ad.obs_names = keep
ad.obs["state"] = [c in state_cells for c in keep]
sc.pp.log1p(ad)
sc.pp.filter_genes(ad, min_cells=10)
print("adata:", ad.shape, flush=True)
sc.tl.rank_genes_groups(ad, "state", method="wilcoxon", pts=True, n_genes=ad.shape[1])
res = sc.get.rank_genes_groups_df(ad, group=True)
res = res.rename(columns={"names": "gene", "logfoldchanges": "logFC", "pvals_adj": "padj", "pct_nz_group": "pct_in", "pct_nz_reference": "pct_out"})
res = res.sort_values("logFC", ascending=False)
res.to_csv(os.path.join(OUTD, "deg_state_vs_rest.csv"), index=False)
print("DEG rows:", len(res), "| up(padj<0.05,logFC>0.5):", int(((res.padj < 0.05) & (res.logFC > 0.5)).sum()),
      "| down:", int(((res.padj < 0.05) & (res.logFC < -0.5)).sum()), flush=True)

## 火山图
r = res.dropna(subset=["padj"]).copy()
r["neglog"] = -np.log10(r["padj"].clip(lower=1e-300))
fig, ax = plt.subplots(figsize=(80/25.4, 64/25.4))
sig = (r.padj < 0.05) & (r.logFC.abs() > 0.5)
ax.scatter(r.loc[~sig, "logFC"], r.loc[~sig, "neglog"], s=1.5, c="#BBBBBB", alpha=0.5, linewidths=0)
ax.scatter(r.loc[sig & (r.logFC > 0), "logFC"], r.loc[sig & (r.logFC > 0), "neglog"], s=2.5, c="#DC0000", alpha=0.8, linewidths=0, label="up in state")
ax.scatter(r.loc[sig & (r.logFC < 0), "logFC"], r.loc[sig & (r.logFC < 0), "neglog"], s=2.5, c="#3C5488", alpha=0.8, linewidths=0, label="down in state")
top = pd.concat([r[sig].nlargest(10, "neglog"), r[sig & (r.logFC < 0)].nsmallest(6, "logFC")]).drop_duplicates("gene")
for _, t in top.iterrows():
    ax.annotate(t["gene"], (t["logFC"], t["neglog"]), fontsize=5.4, ha="center", va="bottom")
ax.axvline(0.5, ls=":", lw=0.4, c="grey"); ax.axvline(-0.5, ls=":", lw=0.4, c="grey")
ax.axhline(-np.log10(0.05), ls=":", lw=0.4, c="grey")
ax.set_xlabel("log2 fold change (state vs rest)", size=6); ax.set_ylabel("-log10 adjusted p", size=6)
ax.tick_params(labelsize=5.4); ax.legend(fontsize=5.4, frameon=False)
ax.set_title("State vs rest DEG (GSE266577; %d state + 20,000 sampled)" % len(keep_state), size=6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "S14a_volcano.png"), dpi=600)
plt.close(fig)
print("volcano saved", flush=True)

## GSEA（Hallmark via gseapy prerank；ranking=logFC*sign）
try:
    import gseapy as gp
    rk = res.dropna(subset=["logFC"])[["gene", "logFC"]].drop_duplicates("gene").sort_values("logFC", ascending=False)
    pre = gp.prerank(rnk=rk, gene_sets="MSigDB_Hallmark_2020", min_size=10, max_size=500, permutation_num=200,
                     seed=42, outdir=None, verbose=False)
    g = pre.res2d.copy()
    g["FDR q-val"] = pd.to_numeric(g["FDR q-val"], errors="coerce")
    g["NES"] = pd.to_numeric(g["NES"], errors="coerce")
    g = g.dropna(subset=["FDR q-val"]).sort_values("NES")
    g.to_csv(os.path.join(OUTD, "gsea_hallmark_state.csv"), index=False)
    gg = pd.concat([g.head(10), g.tail(10)])
    fig, ax = plt.subplots(figsize=(90/25.4, 70/25.4))
    cols = ["#3C5488" if v < 0 else "#DC0000" for v in gg["NES"]]
    ax.barh(range(len(gg)), gg["NES"], color=cols, height=0.65)
    ax.set_yticks(range(len(gg))); ax.set_yticklabels([t[:46] for t in gg["Term"]], fontsize=5.0)
    ax.axvline(0, c="k", lw=0.4); ax.set_xlabel("NES (Hallmark)", size=6); ax.tick_params(labelsize=5.2)
    ax.set_title("GSEA Hallmark (state vs rest; prerank)", size=6)
    fig.tight_layout(); fig.savefig(os.path.join(QRY, "S14b_gsea.png"), dpi=600); plt.close(fig)
    print("GSEA saved; top:", g.tail(3)[["Term","NES"]].to_dict("records"), flush=True)
except Exception as e:
    print("GSEA skipped/failed:", repr(e)[:200], flush=True)
print("DONE 48a")
