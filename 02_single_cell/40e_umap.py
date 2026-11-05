# -*- coding: utf-8 -*-
"""40e: 全量 UMAP（GSE266577 137k；scanpy 默认参数 + n=30/min_dist 0.3/seed42；33 决议）。"""
import numpy as np, pandas as pd, sys, os
import scipy.sparse as sp
import scanpy as sc
sys.stdout.reconfigure(encoding='utf-8')
d = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\umap_input"
genes = [l.strip() for l in open(os.path.join(d, "GSE266577_full.genes.tsv"), encoding="utf-8")]
cells = [l.strip() for l in open(os.path.join(d, "GSE266577_full.cells.tsv"), encoding="utf-8")]
dt = pd.read_csv(os.path.join(d, "GSE266577_full.triplet.gz"))
gidx = {g: i for i, g in enumerate(genes)}
sidx = {c: i for i, c in enumerate(cells)}
C = sp.csr_matrix((dt["v"].to_numpy(), (dt["g"].map(gidx).to_numpy(), dt["c"].map(sidx).to_numpy())),
                  shape=(len(genes), len(cells)))
print("counts:", C.shape, flush=True)
## 仅取蛋白编码/表达基因过滤（scanpy 常规）
ad = sc.AnnData(X=C.T.astype(np.float32), obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))
del C, dt
sc.pp.filter_genes(ad, min_cells=3)
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)
sc.pp.highly_variable_genes(ad, n_top_genes=3000)
ad = ad[:, ad.var["highly_variable"]].copy()
print("hvg:", ad.shape, flush=True)
sc.tl.pca(ad, n_comps=30, random_state=42)
sc.pp.neighbors(ad, n_neighbors=30, random_state=42)
sc.tl.umap(ad, min_dist=0.3, random_state=42)
out = pd.DataFrame({"cell": ad.obs_names, "umap1": ad.obsm["X_umap"][:, 0], "umap2": ad.obsm["X_umap"][:, 1]})
out.to_csv(os.path.join(d, "umap_full_266577.csv"), index=False)
print("saved umap:", out.shape, flush=True)
print("DONE 40e")
