# -*- coding: utf-8 -*-
"""Validate the Python UCell implementation against the R UCell test output."""
import gzip, os, time
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
cols = [pos[c] for c in cells]
X = sp.csc_matrix(sio.mmread(gzip.open(os.path.join(ROOT, "data_sc", "GSE266577",
                                                     "GSE266577_counts_raw.mtx.gz"), "rb")))[:, cols]
N = X.shape[0]
gidx = {g: i for i, g in enumerate(feat)}
gp = {m: np.array([gidx[g] for g in sets[m] if g in gidx]) for m in mods}
out = {m: np.zeros(len(cells)) for m in mods}
for j in range(len(cells)):
    col = X[:, j]; idx = col.indices; vals = col.data.astype(float)
    sv = np.sort(vals); nz = N - len(vals); zr = len(vals) + (nz + 1) / 2
    for m in mods:
        g = gp[m]; p = np.searchsorted(idx, g)
        found = (p < len(idx)) & (idx[np.clip(p, 0, len(idx)-1)] == g)
        ranks = np.full(len(g), zr)
        if found.any():
            vv = vals[p[found]]
            left = np.searchsorted(sv, vv, "left"); right = np.searchsorted(sv, vv, "right")
            greater = len(vals) - right
            ranks[found] = greater + (right - left + 1) / 2
        n = len(g); U = ranks.sum() - n*(n+1)/2
        out[m][j] = 1 - U/(n*(N-n))
for m in mods:
    r = np.corrcoef(out[m], test[m+"_UCell"].values)[0, 1]
    d = np.abs(out[m] - test[m+"_UCell"].values).max()
    print(m, "pearson", round(r, 5), "maxdiff", round(float(d), 5),
          "spearman", round(pd.Series(out[m]).corr(pd.Series(test[m+"_UCell"].values), method="spearman"), 5))
