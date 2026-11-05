# -*- coding: utf-8 -*-
"""Exact re-implementation of UCell::u_stat (descending ranks, cap at maxRank)."""
import gzip, os
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

ROOT = r"D:\子宫内膜异位症2"
test = pd.read_csv(os.path.join(ROOT, "R2_revision", "04_analysis", "results", "A4_ucell_TEST.csv"))
cells = test.cell.tolist()
genes = pd.read_csv(os.path.join(ROOT, "_gene_sets.csv"))
sets = {k: v.gene.tolist() for k, v in genes.groupby("set")}
mods = ["M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"]
feat = pd.read_csv(os.path.join(ROOT, "data_sc", "GSE266577", "GSE266577_seurat_features.txt.gz"), header=None)[0].tolist()
bar = pd.read_csv(os.path.join(ROOT, "data_sc", "GSE266577", "GSE266577_barcodes.txt.gz"), header=None)[0].tolist()
pos = {b: i for i, b in enumerate(bar)}
X = sp.csc_matrix(sio.mmread(gzip.open(os.path.join(ROOT, "data_sc", "GSE266577",
                                                     "GSE266577_counts_raw.mtx.gz"), "rb")))[:, [pos[c] for c in cells]]
N = X.shape[0]; MR = 1500
gidx = {g: i for i, g in enumerate(feat)}
gp = {m: np.array([gidx[g] for g in sets[m] if g in gidx]) for m in mods}
out = {m: np.zeros(len(cells)) for m in mods}
for j in range(len(cells)):
    col = X[:, j]; idx = col.indices; vals = col.data.astype(float)
    sv = np.sort(vals)
    nz = N - len(vals)
    zero_rank = len(vals) + (nz + 1) / 2
    for m in mods:
        g = gp[m]; p = np.searchsorted(idx, g)
        found = (p < len(idx)) & (idx[np.clip(p, 0, len(idx)-1)] == g)
        ranks = np.full(len(g), zero_rank)
        if found.any():
            vv = vals[p[found]]
            left = np.searchsorted(sv, vv, "left"); right = np.searchsorted(sv, vv, "right")
            greater = len(vals) - right
            ranks[found] = greater + (right - left + 1) / 2
        ranks = np.minimum(ranks, MR)
        n = len(g)
        rmin = n * (n + 1) / 2
        out[m][j] = 1 - (ranks.sum() - rmin) / (n * MR - rmin)
for m in mods:
    r = np.corrcoef(out[m], test[m + "_UCell"].values)[0, 1]
    print(m, "pearson", round(r, 6), "maxdiff", round(float(np.abs(out[m]-test[m+"_UCell"].values).max()), 6))
