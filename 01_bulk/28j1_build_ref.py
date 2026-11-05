# -*- coding: utf-8 -*-
"""28j1: cell2location 参考构建+训练（GSE266577 12 谱系 × ≤300 细胞）。"""
import numpy as np, pandas as pd, sys, gzip, os
import scipy.sparse as sp
import anndata as ad
sys.stdout.reconfigure(encoding='utf-8')

LIN = {"Epithelial cells": "Epithelial", "Fibroblasts": "Fibroblast", "Mesothelial cells": "Mesothelial",
  "Endothelial cells": "Endothelial", "Macrophages": "MonoMac", "Classical monocytes": "MonoMac",
  "Migratory DCs": "DC", "DC1": "DC", "DC2": "DC", "pDC": "DC", "Mast cells": "Mast",
  "Tem/Trm cytotoxic T cells": "T", "Tcm/Naive helper T cells": "T", "Regulatory T cells": "T",
  "Type 17 helper T cells": "T", "NK cells": "NK", "CD16- NK cells": "NK", "CD16+ NK cells": "NK",
  "ILC": "NK", "Memory B cells": "B", "Naive B cells": "B", "Plasma cells": "Plasma"}

## 读取 GSE266577 原始计数（triplet gz；CRLF 容错）
p = r"D:\子宫内膜异位症2\data_sc\GSE266577\GSE266577_counts_raw.mtx.gz"
print("reading counts...", flush=True)
with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
    data = fh.read()
lines = data.split("\n")
# 跳过 % 行（banner/dims）
i = 0
while lines[i].startswith("%"):
    if lines[i].startswith("%%MatrixMarket"):
        i += 1
        break
    i += 1
dims = lines[i].split()
nrow, ncol, nnz = int(dims[0]), int(dims[1]), int(dims[2])
body = "\n".join(lines[i+1:])
arr = np.fromstring(body, dtype=np.int64, sep=" ").reshape(nnz, 3)
C = sp.csr_matrix((arr[:, 2], (arr[:, 0] - 1, arr[:, 1] - 1)), shape=(nrow, ncol))
del arr, body, data; 
feat = [l.strip() for l in gzip.open(r"D:\子宫内膜异位症2\data_sc\GSE266577\GSE266577_seurat_features.txt.gz", "rt", encoding="utf-8")]
bar = [l.strip() for l in gzip.open(r"D:\子宫内膜异位症2\data_sc\GSE266577\GSE266577_barcodes.txt.gz", "rt", encoding="utf-8")]
print("counts:", C.shape, "cells:", len(bar), "genes:", len(feat), flush=True)

## 细胞→谱系
sc = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577.csv", usecols=["cell", "cell_type"])
sc["lineage"] = sc["cell_type"].map(LIN)
print("lineage counts:", sc.groupby("lineage").size().to_dict(), flush=True)

## 参考子集 ≤300/谱系
rng = np.random.default_rng(2026)
keep = []
for lin, g in sc.groupby("lineage"):
    cells = g["cell"].sample(min(300, len(g)), random_state=2026).tolist()
    keep.extend(cells)
print("ref cells:", len(keep), flush=True)
bar_idx = {c: i for i, c in enumerate(bar)}
idx = [bar_idx[c] for c in keep]
Csub = C[:, idx].tocsr()
annt = sc.set_index("cell").loc[keep, "lineage"].to_numpy()
ref = ad.AnnData(X=Csub.T.tocsr().astype(np.float32),
                 obs=pd.DataFrame({"cell_type": annt}, index=keep),
                 var=pd.DataFrame(index=feat))
print("ref adata:", ref.shape, flush=True)
ref.write_h5ad(r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial\deconv_input\reference_lineage12.h5ad")
print("DONE 28j1 (ref script follows separately)")
