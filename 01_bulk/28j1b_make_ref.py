# -*- coding: utf-8 -*-
"""28j1b: 参考 AnnData 构建（读 28j1R 紧凑 triplet）。"""
import numpy as np, pandas as pd, sys, os
import scipy.sparse as sp
import anndata as ad
sys.stdout.reconfigure(encoding='utf-8')
DE = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial\deconv_input"
genes = [l.strip() for l in open(os.path.join(DE, "reference_lineage12.genes.tsv"), encoding="utf-8")]
spots = pd.read_csv(os.path.join(DE, "reference_lineage12.spots.tsv"))
dt = pd.read_csv(os.path.join(DE, "reference_lineage12.triplet.gz"))
gidx = {g: i for i, g in enumerate(genes)}
sidx = {s: i for i, s in enumerate(spots["spot"])}
C = sp.csr_matrix((dt["v"].to_numpy(), (dt["c"].map(sidx).to_numpy(), dt["g"].map(gidx).to_numpy())),
                  shape=(len(spots), len(genes)))
ref = ad.AnnData(X=C.astype(np.float32), obs=spots.set_index("spot"),
                 var=pd.DataFrame(index=genes))
print("ref adata:", ref.shape, flush=True)
ref.write_h5ad(os.path.join(DE, "reference_lineage12.h5ad"))
print("DONE 28j1b")
