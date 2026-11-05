# -*- coding: utf-8 -*-
"""28j2 v2: cell2location 0.1.5 新 API（RegressionModel → cell_state_df → Cell2location）。"""
import numpy as np, pandas as pd, sys, os, gzip, glob
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
sys.stdout.reconfigure(encoding='utf-8')
from cell2location.models import RegressionModel, Cell2location

OUT = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial"
DE = os.path.join(OUT, "deconv_input")
GSE = sys.argv[1]

ref = ad.read_h5ad(os.path.join(DE, "reference_lineage12.h5ad"))
ref = ref[:, ~ref.var_names.duplicated()].copy()
print("ref:", ref.shape, flush=True)
state_csv = os.path.join(DE, "cell_state_df.csv")
if os.path.exists(state_csv):
    state_df = pd.read_csv(state_csv, index_col=0)
    print("loaded cached cell_state_df:", state_df.shape, flush=True)
    genes = state_df.index.to_numpy()
else:
    ## 细胞类型 × 基因 平均原始计数（计数尺度，LocationModel 所需）
    genes = ref.var_names.to_numpy()
    ct = ref.obs["cell_type"].astype(str).to_numpy()
    ctype_list = np.unique(ct)
    Xr = ref.X.toarray() if hasattr(ref.X, "toarray") else np.asarray(ref.X)
    state_df = pd.DataFrame(0.0, index=genes, columns=ctype_list)
    for ct0 in ctype_list:
        sel = np.where(ct == ct0)[0]
        state_df[ct0] = Xr[sel].mean(axis=0)
    state_df = state_df.loc[:, state_df.sum(axis=0) > 0]
    print("cell_state_df (raw-count means):", state_df.shape, flush=True)
    state_df.to_csv(state_csv)

## 逐样本
for fp in sorted(glob.glob(os.path.join(DE, f"{GSE}_*.triplet.gz"))):
    code = os.path.basename(fp).replace(f"{GSE}_", "").replace(".triplet.gz", "")
    genes_s = [l.strip() for l in open(os.path.join(DE, f"{GSE}_{code}.genes.tsv"), encoding="utf-8")]
    spots_s = [l.strip() for l in open(os.path.join(DE, f"{GSE}_{code}.spots.tsv"), encoding="utf-8")]
    dt = pd.read_csv(fp)
    gidx = {g: i for i, g in enumerate(genes_s)}
    sidx = {s: i for i, s in enumerate(spots_s)}
    C = sp.csr_matrix((dt["v"].to_numpy(), (dt["c"].map(sidx).to_numpy(), dt["g"].map(gidx).to_numpy())),
                      shape=(len(spots_s), len(genes_s)))
    adata = ad.AnnData(X=C.astype(np.float32), var=pd.DataFrame(index=genes_s),
                       obs=pd.DataFrame(index=spots_s))
    adata = adata[:, ~adata.var_names.duplicated()].copy()
    adata.obs["n_counts"] = np.asarray(C.sum(axis=1)).ravel()
    adata.layers["counts"] = adata.X.copy()
    adata = adata[adata.obs["n_counts"] >= 200].copy()
    if adata.n_obs < 50:
        print("skip", code, adata.n_obs, flush=True)
        continue
    common = [g for g in genes if g in adata.var_names]
    adata2 = adata[:, common].copy()
    csd = state_df.loc[common]
    Cell2location.setup_anndata(adata2, layer="counts")
    lm = Cell2location(adata2, cell_state_df=csd)
    lm.train(max_epochs=100)
    adata2 = lm.export_posterior(adata2, sample_kwargs={"num_samples": 50})
    frac = adata2.obsm["q05_cell_abundance_w_sf"]
    tot = np.asarray(frac).sum(axis=1, keepdims=True)
    frac = np.asarray(frac) / tot
    names = [str(c) for c in csd.columns]
    out = pd.DataFrame(frac, index=adata2.obs_names, columns=names)
    out["sample"] = code
    out = out.reset_index().rename(columns={"index": "spot"})
    out.to_csv(os.path.join(OUT, f"{GSE}_{code}_deconv.csv"), index=False)
    print("saved", code, out.shape, flush=True)

## state vs 非 state
sc_scores = pd.read_csv(os.path.join(OUT, f"{GSE}_domains.csv"))
allf = []
for fp in sorted(glob.glob(os.path.join(OUT, f"{GSE}_*_deconv.csv"))):
    allf.append(pd.read_csv(fp))
allf = pd.concat(allf, ignore_index=True)
merged = allf.merge(sc_scores[["sample", "spot", "state_high", "decouple", "domain"]], on=["sample", "spot"], how="inner")
comp_cols = [c for c in merged.columns if c not in ("sample", "spot", "state_high", "decouple", "domain")]
mc = merged.groupby("state_high")[comp_cols].mean()
print("\nstate vs non-state 平均组成:\n", mc.round(4).T.to_string(), flush=True)
merged.to_csv(os.path.join(OUT, f"{GSE}_deconv_state.csv"), index=False)
print("DONE 28j2", GSE)
