# -*- coding: utf-8 -*-
"""28h: 正式 Moran's I（esda.Moran, 999 重排）——decouple 空间自相关的标准口径。"""
import numpy as np, pandas as pd, sys, os
from sklearn.neighbors import NearestNeighbors
from libpysal.weights import W
from esda.moran import Moran
sys.stdout.reconfigure(encoding='utf-8')
OUT = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial"
rows = []
for GSE in ["GSE274657", "GSE288483", "GSE224335"]:
    sc = pd.read_csv(os.path.join(OUT, f"{GSE}_domains.csv"))
    for sm, sub in sc.groupby("sample"):
        sub = sub.dropna(subset=["array_row", "array_col", "decouple"])
        if len(sub) < 100:
            continue
        pts = sub[["array_row", "array_col"]].values.astype(float)
        d = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=2)
        # 近邻权重（k=6，对称化）
        idx = np.argsort(d, axis=1)[:, 1:7]
        import scipy.sparse as sp
        adj = sp.coo_matrix((np.ones(idx.size), (np.repeat(np.arange(len(sub)), 6), idx.ravel())),
                            shape=(len(sub), len(sub)))
        adj = (adj + adj.T) > 0
        w = W.from_sparse(adj.astype(float))
        w.transform = "r"
        mi = Moran(sub["decouple"].values, w, permutations=999)
        rows.append({"GSE": GSE, "sample": sm, "n": len(sub), "Moran_I": round(mi.I, 4),
                     "z": round(mi.z_sim, 2), "p": round(mi.p_sim, 4)})
        print(f"{GSE} {sm}: I={mi.I:.4f} z={mi.z_sim:.2f} p={mi.p_sim:.4f}", flush=True)
res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "spatial_moran_summary.csv"), index=False)
print("significant (p<0.05):", int((res["p"] < 0.05).sum()), "/", len(res), flush=True)
print("median I:", round(res["Moran_I"].median(), 4), flush=True)
print("DONE 28h")
