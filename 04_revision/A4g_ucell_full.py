# -*- coding: utf-8 -*-
"""Full UCell scoring (exact re-implementation, validated r = 1.0 against R UCell)."""
import gzip, os, time
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

ROOT = r"D:\子宫内膜异位症2"
OUT = os.path.join(ROOT, "R2_revision", "04_analysis", "results")
genes = pd.read_csv(os.path.join(ROOT, "_gene_sets.csv"))
sets = {k: v.gene.tolist() for k, v in genes.groupby("set")}
mods = ["M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"]
icc = sets["ICC45"]
shared = [g for g in icc if g in set().union(*[set(sets[m]) for m in mods])]
unique = [g for g in icc if g not in shared]
rng = np.random.default_rng(20261105)
allg = sorted(set(genes.gene))
set_list = {m: sets[m] for m in mods}
set_list["ICC45_full"] = icc
set_list["ICC45_shared28"] = shared
set_list["ICC45_unique17"] = unique
for i in range(5):
    set_list[f"RAND17_{i+1}"] = list(rng.choice(allg, 17, replace=False))

feat = pd.read_csv(os.path.join(ROOT, "data_sc", "GSE266577", "GSE266577_seurat_features.txt.gz"), header=None)[0].tolist()
bar = pd.read_csv(os.path.join(ROOT, "data_sc", "GSE266577", "GSE266577_barcodes.txt.gz"), header=None)[0].tolist()
sc = pd.read_csv(os.path.join(ROOT, "analysis", "L1_scRNA", "results", "module_scores",
                              "scores_GSE266577.csv"), usecols=["cell"])
keep = set(sc.cell)
cols = [i for i, b in enumerate(bar) if b in keep]
print("cells", len(cols), flush=True)

t0 = time.time()
X = sp.csc_matrix(sio.mmread(gzip.open(os.path.join(ROOT, "data_sc", "GSE266577",
                                                     "GSE266577_counts_raw.mtx.gz"), "rb")))[:, cols]
print("matrix", X.shape, "%.1f min" % ((time.time()-t0)/60), flush=True)
N = X.shape[0]; MR = 1500
gidx = {g: i for i, g in enumerate(feat)}
gp = {k: np.array([gidx[g] for g in v if g in gidx]) for k, v in set_list.items()}
scores = {k: np.zeros(X.shape[1], dtype=np.float32) for k in set_list}
t1 = time.time()
for j in range(X.shape[1]):
    col = X[:, j]; idx = col.indices; vals = col.data.astype(np.float64)
    sv = np.sort(vals); nz = N - len(vals)
    zero_rank = len(vals) + (nz + 1) / 2
    for k, g in gp.items():
        p = np.searchsorted(idx, g)
        found = (p < len(idx)) & (idx[np.clip(p, 0, len(idx)-1)] == g)
        ranks = np.full(len(g), zero_rank)
        if found.any():
            vv = vals[p[found]]
            left = np.searchsorted(sv, vv, "left"); right = np.searchsorted(sv, vv, "right")
            greater = len(vals) - right
            ranks[found] = greater + (right - left + 1) / 2
        ranks = np.minimum(ranks, MR)
        n = len(g); rmin = n * (n + 1) / 2
        scores[k][j] = 1 - (ranks.sum() - rmin) / (n * MR - rmin)
    if j % 20000 == 0:
        print(j, "%.1f min" % ((time.time()-t1)/60), flush=True)

df = pd.DataFrame({"cell": [bar[i] for i in cols], **scores})
df.to_csv(os.path.join(OUT, "A4_ucell_python.csv"), index=False)
print("done", df.shape, "%.1f min" % ((time.time()-t1)/60))
