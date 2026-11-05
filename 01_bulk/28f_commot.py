# -*- coding: utf-8 -*-
"""28f: COMMOT 空间通讯（250µm, 精选 LR 集）→ state spot 的发送/接收信号对比。"""
import numpy as np
np.Inf = np.inf  # 必须先于 commot 导入（其默认参数在 import 时求值）
import pandas as pd, sys, os, glob
import anndata as ad
import commot
import commot.tl as ctl
sys.stdout.reconfigure(encoding='utf-8')

LRD = {
    "FASLG-FAS": ("FASLG", "FAS"), "TNF-TNFRSF1A": ("TNF", "TNFRSF1A"),
    "TRAIL-TNFRSF10B": ("TNFSF10", "TNFRSF10B"), "IL1B-IL1R1": ("IL1B", "IL1R1"),
    "IL18-IL18R1": ("IL18", "IL18R1"), "IFNG-IFNGR1": ("IFNG", "IFNGR1"),
    "TGFB1-TGFBR1": ("TGFB1", "TGFBR1"), "CSF1-CSF1R": ("CSF1", "CSF1R"),
    "SPP1-CD44": ("SPP1", "CD44"), "MDK-SDC2": ("MDK", "SDC2"), "MDK-SDC1": ("MDK", "SDC1"),
    "CXCL12-CXCR4": ("CXCL12", "CXCR4"), "CXCL16-CXCR6": ("CXCL16", "CXCR6"),
    "CCL2-CCR2": ("CCL2", "CCR2"), "CCL5-CCR5": ("CCL5", "CCR5"),
    "PDL1-PD1": ("CD274", "PDCD1"), "NECTIN2-TIGIT": ("NECTIN2", "TIGIT"),
    "HGF-MET": ("HGF", "MET"), "EGF-EGFR": ("EGF", "EGFR"),
    "VEGFA-KDR": ("VEGFA", "KDR"), "IL6-IL6R": ("IL6", "IL6R"), "IL10-IL10RA": ("IL10", "IL10RA"),
    "LGALS9-HAVCR2": ("LGALS9", "HAVCR2"), "IDO1": ("IDO1", "IDO1"),
}
OUT = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial"

## 官方 CellChatDB（commot 内置）过滤到我们导出的 LR 基因
import commot.pp as cpp, os as _os
DB = os.path.join(_os.path.dirname(commot.__file__), "preprocessing", "_data", "LRdatabase", "CellChat", "CellChatDB.ligrec.human.csv")
full_db = pd.read_csv(DB, header=None, names=["ligand", "receptor", "pathway_name", "annotation"])
print("full CellChatDB pairs:", len(full_db), flush=True)
LR_GENES = set(pd.read_csv(os.path.join(OUT, "lr_input", sorted(glob.glob(os.path.join(OUT, "lr_input", f"{sys.argv[1]}_*.csv.gz")))[0].split("\\")[-1].split("_")[0] + "_" if False else sorted(glob.glob(os.path.join(OUT, "lr_input", f"{sys.argv[1]}_*.csv.gz")))[0]), index_col=0).index)
print("exported LR genes:", len(LR_GENES), flush=True)
db = full_db[full_db["ligand"].isin(LR_GENES) & full_db["receptor"].isin(LR_GENES)]
print("filtered pairs:", len(db), flush=True)
db = db.rename(columns={"ligand": 0, "receptor": 1, "pathway_name": 2, "annotation": 3})
db.to_csv(os.path.join(OUT, "cellchatdb_filtered.csv"), index=False)

GSE = sys.argv[1]
scores = pd.read_csv(os.path.join(OUT, f"{GSE}_domains.csv"))
rows = []
for fp in sorted(glob.glob(os.path.join(OUT, "lr_input", f"{GSE}_*.csv.gz"))):
    code = os.path.basename(fp).replace(f"{GSE}_", "").replace(".csv.gz", "")
    X = pd.read_csv(fp, index_col=0)  # genes × spots
    sc = scores[(scores["sample"] == code) & scores["array_row"].notna() & scores["array_col"].notna()]
    sc = sc.set_index("spot")
    common = [b for b in sc.index if b in X.columns]
    if len(common) < 50:
        print("skip", code, len(common), flush=True)
        continue
    Xs = X[common].T  # spots × genes
    sub = sc.loc[common]
    adata = ad.AnnData(X=Xs.values.astype(np.float32), obs=sub[["state_high", "decouple"]])
    adata.var_names = Xs.columns.tolist()
    adata.obsm["spatial"] = sub[["array_row", "array_col"]].values.astype(float)
    # 邻接（kNN=6, array 坐标）
    import sklearn.neighbors as nb
    nn = nb.NearestNeighbors(n_neighbors=6).fit(adata.obsm["spatial"])
    dist, idx = nn.kneighbors(adata.obsm["spatial"])
    import scipy.sparse as sp
    conn = sp.csr_matrix((np.ones(idx.size), (np.repeat(np.arange(len(idx)), 6), idx.ravel())), shape=(len(idx), len(idx)))
    conn = (conn + conn.T) > 0
    adata.obsp["spatial_connectivities"] = conn.astype(float)
    # COMMOT（官方 CellChatDB 过滤版）
    ctl.spatial_communication(adata, database_name="df_ligrec", df_ligrec=db, dis_thr=250)
    s_sum = adata.obsm.get("commot-df_ligrec-sum-sender")
    r_sum = adata.obsm.get("commot-df_ligrec-sum-receiver")
    if s_sum is None or r_sum is None:
        print("no summary keys!", list(adata.obsm.keys()), flush=True)
        continue
    s_sum = pd.DataFrame(s_sum) if not isinstance(s_sum, pd.DataFrame) else s_sum
    r_sum = pd.DataFrame(r_sum) if not isinstance(r_sum, pd.DataFrame) else r_sum
    st = sub["state_high"].values.astype(bool)
    if st.sum() < 5:
        print("skip", code, "state n:", st.sum(), flush=True)
        continue
    per_lr = {}
    for k in adata.obsp.keys():
        if k.startswith("commot-df_ligrec-") and k not in ("commot-df_ligrec-sum-sender", "commot-df_ligrec-sum-receiver"):
            mt = adata.obsp[k].toarray()
            st_v = mt[st].sum() / max(st.sum(), 1)
            ns_v = mt[~st].sum() / max((~st).sum(), 1)
            per_lr[k] = (float(st_v) + 1e-9) / (float(ns_v) + 1e-9)
    row = {"GSE": GSE, "sample": code, "n_spots": len(sub), "n_state": int(st.sum()),
           "state_send": float(s_sum[st].values.mean().mean()), "state_recv": float(r_sum[st].values.mean().mean()),
           "nonstate_send": float(s_sum[~st].values.mean().mean()), "nonstate_recv": float(r_sum[~st].values.mean().mean())}
    row["top_lr"] = sorted(per_lr.items(), key=lambda kv: -kv[1])[:6]
    rows.append(row)
    print(code, "state", st.sum(), "send ratio(state/non):", round(row["state_send"] / (row["nonstate_send"] + 1e-9), 2),
          "recv ratio:", round(row["state_recv"] / (row["nonstate_recv"] + 1e-9), 2), flush=True)
    print("   top LR ratio:", [(k.replace("commot-df_ligrec-", ""), round(v, 2)) for k, v in row["top_lr"]], flush=True)
res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, f"{GSE}_commot_summary.csv"), index=False)
print("DONE 28f", GSE)
