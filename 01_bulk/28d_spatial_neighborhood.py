# -*- coding: utf-8 -*-
"""28d: 空间邻域分析（kNN=6, array coords）：
- state_high 邻域富集（OR + Fisher per 样本，跨样本汇总）
- 连续 decouple 的空间自相关（邻域-自身差）
- 代表样本：状态地图 + 域地图（PNG 规范件初稿）
"""
import os, sys
import numpy as np, pandas as pd
from sklearn.neighbors import NearestNeighbors
from scipy.stats import fisher_exact, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.size": 6.5, "figure.dpi": 600, "font.family": "Arial",
                     "axes.linewidth": 0.35})

OUT = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial"
COILS = ["array_row", "array_col", "pxl_row", "pxl_col"]
rows = []
fig_samples = {"GSE274657": "Patient_16", "GSE288483": "Pt2-4", "GSE224335": "OCCR2v"}

for GSE in sys.argv[1:]:
    sc = pd.read_csv(os.path.join(OUT, f"{GSE}_domains.csv"))
    have_pos = all(c in sc.columns for c in COILS)
    print("==", GSE, "positions:", have_pos, flush=True)
    if not have_pos:
        continue
    for sm, sub in sc.groupby("sample"):
        sub = sub.dropna(subset=["array_row", "array_col"])
        if len(sub) < 100:
            continue
        pts = sub[["array_row", "array_col"]].values.astype(float)
        nn = NearestNeighbors(n_neighbors=6).fit(pts)
        dist, idx = nn.kneighbors(pts)
        neigh = idx[:, 1:]
        # 1) state_high 邻域富集
        st = sub["state_high"].fillna(False).values
        n_st = int(st.sum())
        if n_st < 5 or n_st > len(st) * 0.9:
            print("  skip", sm, "state n too small/large:", n_st, flush=True)
            continue
        # 邻居对：state spot 的邻居数 vs 期望
        obs_st_neigh = sum(st[neigh[np.where(st)[0]]])
        exp_st_neigh = st.mean() * neigh[np.where(st)[0]].size
        # OR：以"spot 为 state 的邻居比"构造 2×2：state spot 的邻居 state 比例 vs 非 state spot 的邻居 state 比例
        st_idx = np.where(st)[0]; nst_idx = np.where(~st)[0]
        a = int(st[neigh[st_idx]].sum()); b = int(neigh[st_idx].size - a)
        c = int(st[neigh[nst_idx]].sum()); d = int(neigh[nst_idx].size - c)
        orr, p = fisher_exact([[a, b], [c, d]])
        rows.append({"GSE": GSE, "sample": sm, "n_spots": len(sub), "n_state": n_st,
                     "state_rate": round(st.mean(), 4), "neigh_OR": round(orr, 2), "fisher_p": p,
                     "obs_frac": round(a / (a + b), 4), "exp_frac": round(c / (c + d), 4)})
        # 2) 连续自相关：邻居平均 decouple - 自身
        dco = sub["decouple"].values
        nb_mean = dco[neigh].mean(axis=1)
        rho, pp = spearmanr(dco, nb_mean)
        rows[-1]["moran_like_rho"] = round(rho, 3)
        rows[-1]["rho_p"] = pp
        # figure: state map + domain map
        if fig_samples.get(GSE) == sm:
            fig, ax = plt.subplots(1, 3, figsize=(180/25.4, 60/25.4))
            cols = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000"]
            ax0 = ax[0].scatter(pts[:, 1], pts[:, 0], c=sub["decouple"], cmap="RdBu_r", s=0.6, linewidths=0)
            plt.colorbar(ax0, ax=ax[0], shrink=0.8); ax[0].set_title(f"{sm} decouple")
            im1 = ax[1].scatter(pts[:, 1], pts[:, 0], c=sub["domain"], cmap=matplotlib.colors.ListedColormap(cols), s=0.6, linewidths=0)
            ax[1].set_title("NMF domain")
            ax[2].scatter(pts[:, 1], pts[:, 0], c="#DDDDDD", s=0.6, linewidths=0)
            ax[2].scatter(pts[st_idx, 1], pts[st_idx, 0], c="#E64B35", s=1.2, linewidths=0)
            ax[2].set_title("state_high")
            for a in ax:
                a.set_aspect("equal"); a.axis("off")
            fig.tight_layout()
            fig.savefig(os.path.join(OUT, f"{GSE}_{sm}_maps.png"), dpi=600)
            plt.close(fig)
            print("  fig saved", f"{GSE}_{sm}_maps.png", flush=True)
res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "spatial_neighborhood_summary.csv"), index=False)
print(res.to_string(index=False), flush=True)
print("DONE 28d")
