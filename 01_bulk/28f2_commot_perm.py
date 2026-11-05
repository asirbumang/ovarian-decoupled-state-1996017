# -*- coding: utf-8 -*-
"""28f2: COMMOT 置换显著性——保存每 LR spot 信号矩阵 → 1000× 标签置换 → 每 LR p。
（对 GSE274657 4 样本；ratio 尺度上的条件置换）"""
import numpy as np, pandas as pd, sys, os, glob
import anndata as ad
import scipy.sparse as sp
sys.stdout.reconfigure(encoding='utf-8')
np.Inf = np.inf
import commot, commot.tl as ctl

OUT = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial"
GSE = sys.argv[1]
db = pd.read_csv(os.path.join(OUT, "cellchatdb_filtered.csv"))
rng = np.random.default_rng(2026)

out_rows = []
for fp in sorted(glob.glob(os.path.join(OUT, "lr_input", f"{GSE}_*.csv.gz"))):
    code = os.path.basename(fp).replace(f"{GSE}_", "").replace(".csv.gz", "")
    X = pd.read_csv(fp, index_col=0)
    sc = pd.read_csv(os.path.join(OUT, f"{GSE}_domains.csv"))
    sc = sc[(sc["sample"] == code) & sc["array_row"].notna()].set_index("spot")
    common = [b for b in sc.index if b in X.columns]
    if len(common) < 50: continue
    Xs = X[common].T
    sub = sc.loc[common]
    adata = ad.AnnData(X=Xs.values.astype(np.float32), obs=sub[["state_high"]])
    adata.var_names = Xs.columns.tolist()
    adata.obsm["spatial"] = sub[["array_row", "array_col"]].values.astype(float)
    import sklearn.neighbors as nb
    nn = nb.NearestNeighbors(n_neighbors=6).fit(adata.obsm["spatial"])
    _, idx = nn.kneighbors(adata.obsm["spatial"])
    conn = sp.csr_matrix((np.ones(idx.size), (np.repeat(np.arange(len(idx)), 6), idx.ravel())), shape=(len(idx), len(idx)))
    adata.obsp["spatial_connectivities"] = ((conn + conn.T) > 0).astype(float)
    ctl.spatial_communication(adata, database_name="df_ligrec", df_ligrec=db, dis_thr=250)
    st = sub["state_high"].values.astype(bool)
    n_st = int(st.sum())
    if n_st < 5:
        print("skip", code, "state n:", n_st, flush=True); continue
    # 每 LR 信号矩阵（obsp：commot-df_ligrec-<lr>）→ 每 spot 的发送总量（行和）
    lr_mat = {}
    for k in adata.obsp.keys():
        if k.startswith("commot-df_ligrec-") and k not in ("commot-df_ligrec-sum-sender", "commot-df_ligrec-sum-receiver"):
            mt = adata.obsp[k].toarray()
            lr_mat[k.replace("commot-df_ligrec-", "")] = mt
    for lr, mt in lr_mat.items():
        send = mt.sum(axis=1)
        s_obs = send[st].mean() / (send[~st].mean() + 1e-9)
        # 置换：打乱标签 1000×
        cnt = 0
        for _ in range(1000):
            perm = rng.permutation(len(send))
            ps = send[perm][st].mean() / (send[perm][~st].mean() + 1e-9)
            if np.abs(ps) >= np.abs(s_obs):  # 方向性：观察为正时看 ≥
                cnt += 1
        p = (cnt + 1) / 1001
        out_rows.append({"GSE": GSE, "sample": code, "lr": lr, "ratio_obs": round(s_obs, 3), "perm_p": round(p, 4)})
    print(code, "done LRs:", len(lr_mat), flush=True)
res = pd.DataFrame(out_rows)
res.to_csv(os.path.join(OUT, f"{GSE}_commot_perm.csv"), index=False)
print(res.sort_values("perm_p").head(15).to_string(index=False), flush=True)
print("DONE 28f2")
