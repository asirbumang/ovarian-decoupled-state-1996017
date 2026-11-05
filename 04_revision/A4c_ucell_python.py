# -*- coding: utf-8 -*-
"""
A4c: UCell score implemented directly (Mann-Whitney U on within-cell ranks),
validated against the R UCell output on 3,000 cells.

UCell: rank genes within a cell (ascending, ties averaged, zeros tied at the
bottom), U = sum(rank of set genes) - n(n+1)/2, Umax = n(N-n),
score = 1 - U/Umax.
"""
import gzip, json, os, sys, time
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
rand17 = [list(rng.choice(allg, 17, replace=False)) for _ in range(5)]
set_list = {m: sets[m] for m in mods}
set_list["ICC45_full"] = icc
set_list["ICC45_shared28"] = shared
set_list["ICC45_unique17"] = unique
for i, r in enumerate(rand17):
    set_list[f"RAND17_{i+1}"] = list(r)

t0 = time.time()
feat = pd.read_csv(os.path.join(ROOT, "data_sc", "GSE266577", "GSE266577_seurat_features.txt.gz"),
                   header=None)[0].tolist()
bar = pd.read_csv(os.path.join(ROOT, "data_sc", "GSE266577", "GSE266577_barcodes.txt.gz"),
                  header=None)[0].tolist()
sc = pd.read_csv(os.path.join(ROOT, "analysis", "L1_scRNA", "results", "module_scores",
                              "scores_GSE266577.csv"), usecols=["cell"])
keep = set(sc.cell)
cols = [i for i, b in enumerate(bar) if b in keep]
print("matrix cols kept", len(cols), flush=True)

X = sio.mmread(gzip.open(os.path.join(ROOT, "data_sc", "GSE266577",
                                      "GSE266577_counts_raw.mtx.gz"), "rb"))
X = sp.csc_matrix(X)
print("matrix", X.shape, "%.1f min" % ((time.time()-t0)/60), flush=True)
X = X[:, cols]
bar = [bar[i] for i in cols]

N = X.shape[0]
gidx = {g: i for i, g in enumerate(feat)}
gene_pos = {k: np.array([gidx[g] for g in v if g in gidx]) for k, v in set_list.items()}

scores = {k: np.zeros(X.shape[1], dtype=np.float32) for k in set_list}
t1 = time.time()
for j in range(X.shape[1]):
    col = X[:, j]
    idx = col.indices
    vals = col.data.astype(np.float64)
    if len(vals) == 0:
        continue
    order = np.argsort(vals, kind="stable")
    sv = vals[order]
    nz = N - len(vals)
    zero_rank = (nz + 1) / 2.0
    for k, gp in gene_pos.items():
        pos = np.searchsorted(idx, gp)
        found = (pos < len(idx)) & (idx[np.clip(pos, 0, len(idx)-1)] == gp)
        ranks = np.full(len(gp), zero_rank)
        if found.any():
            vv = vals[pos[found]]
            left = np.searchsorted(sv, vv, side="left")
            right = np.searchsorted(sv, vv, side="right")
            ranks[found] = nz + left + (right - left + 1) / 2.0
        n = len(gp)
        U = ranks.sum() - n * (n + 1) / 2.0
        scores[k][j] = 1.0 - U / (n * (N - n))
    if j % 20000 == 0:
        print(j, "%.1f min" % ((time.time()-t1)/60), flush=True)

df = pd.DataFrame({"cell": bar, **scores})
df.to_csv(os.path.join(OUT, "A4_ucell_python.csv"), index=False)
print("done", df.shape, "%.1f min" % ((time.time()-t1)/60))
