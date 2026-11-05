# -*- coding: utf-8 -*-
"""28k: NNLS 反卷积（参考=GSE266577 21 型计数均值；spot 用 log1p(CP10k)）。
输入：deconv_input 的 triplet + cell_state_df 原始计数均值。
输出：每 spot 21 型组成 → state vs 非 state 组成对比。
"""
import numpy as np, pandas as pd, sys, os, gzip, glob
import scipy.sparse as sp
from scipy.optimize import nnls
sys.stdout.reconfigure(encoding='utf-8')
OUT = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial"
DE = os.path.join(OUT, "deconv_input")
GSE = sys.argv[1]

state_df = pd.read_csv(os.path.join(DE, "cell_state_df.csv"), index_col=0)
print("cell_state_df:", state_df.shape, flush=True)
# 参考 profile：log1p(CPK) per gene per type（用计数均值归一）
ref = state_df.values  # genes × types
ref_cpk = ref / ref.sum(axis=0, keepdims=True) * 1e4
ref_log = np.log1p(ref_cpk)

rows = []
for fp in sorted(glob.glob(os.path.join(DE, f"{GSE}_*.triplet.gz"))):
    code = os.path.basename(fp).replace(f"{GSE}_", "").replace(".triplet.gz", "")
    genes = [l.strip() for l in open(os.path.join(DE, f"{GSE}_{code}.genes.tsv"), encoding="utf-8")]
    spots = [l.strip() for l in open(os.path.join(DE, f"{GSE}_{code}.spots.tsv"), encoding="utf-8")]
    dt = pd.read_csv(fp)
    gidx = {g: i for i, g in enumerate(genes)}
    sidx = {s: i for i, s in enumerate(spots)}
    C = sp.csr_matrix((dt["v"].to_numpy(), (dt["c"].map(sidx).to_numpy(), dt["g"].map(gidx).to_numpy())),
                      shape=(len(spots), len(genes)))
    # 共同基因，去重复
    keep = [g for g in state_df.index if g in genes]
    keep = [g for g in dict.fromkeys(keep)]
    Ck = C[:, [genes.index(g) for g in keep]]
    R = ref_log[[list(state_df.index).index(g) for g in keep], :]
    # spot log1p(CP10k)
    cpm = Ck.multiply(1e4 / Ck.sum(axis=1).A1[:, None]).toarray()
    X = np.log1p(cpm)
    # NNLS per spot（21 型）
    n = X.shape[0]
    frac = np.zeros((n, R.shape[1]))
    for i in range(n):
        sol, _ = nnls(R, X[i])
        tot = sol.sum()
        frac[i] = sol / tot if tot > 0 else np.ones(R.shape[1]) / R.shape[1]
    out = pd.DataFrame(frac, index=spots, columns=state_df.columns)
    out["sample"] = code
    out = out.reset_index().rename(columns={"index": "spot"})
    out.to_csv(os.path.join(OUT, f"{GSE}_{code}_nnls.csv"), index=False)
    print("saved", code, out.shape, flush=True)

# state vs 非 state
sc_scores = pd.read_csv(os.path.join(OUT, f"{GSE}_domains.csv"))
allf = []
for fp in sorted(glob.glob(os.path.join(OUT, f"{GSE}_*_nnls.csv"))):
    allf.append(pd.read_csv(fp))
allf = pd.concat(allf, ignore_index=True)
merged = allf.merge(sc_scores[["sample", "spot", "state_high", "decouple", "domain"]], on=["sample", "spot"], how="inner")
comp_cols = [c for c in merged.columns if c not in ("sample", "spot", "state_high", "decouple", "domain")]
mc = merged.groupby("state_high")[comp_cols].mean()
mc["delta"] = mc.loc[True] - mc.loc[False]
print("\nstate vs non-state 平均组成:\n", mc.round(4).T.to_string(), flush=True)
merged.to_csv(os.path.join(OUT, f"{GSE}_nnls_state.csv"), index=False)
print("DONE 28k", GSE)
