# -*- coding: utf-8 -*-
"""40b3: B9 P1——alpha 网格 {0.1,1,10} + 签名替代 {100/300 基因}（响应复审 P1 项）。"""
import numpy as np, pandas as pd, sys, os
import scipy.sparse as sp
import decoupler as dc
from sklearn.linear_model import Ridge
sys.stdout.reconfigure(encoding='utf-8')
OUT = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b9"
net = dc.op.collectri()[["source", "target"]].drop_duplicates()
print("collectri edges:", len(net), flush=True)

def load_tri(nm):
    d = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b9_input"
    genes = [l.strip() for l in open(os.path.join(d, f"{nm}.genes.tsv"), encoding="utf-8")]
    dt = pd.read_csv(os.path.join(d, f"{nm}.triplet.gz"))
    gidx = {g: i for i, g in enumerate(genes)}
    cells = dt["c"].unique(); sidx = {v: i for i, v in enumerate(cells)}
    C = sp.csr_matrix((dt["v"].to_numpy(), (dt["g"].map(gidx).to_numpy(), dt["c"].map(sidx).to_numpy())),
                      shape=(len(genes), len(cells)))
    return genes, C

HOSTS = {"GSE293459_all": ("GSE293459_all", None), "GSE266577_epithelial": ("GSE266577_epi", None),
         "GSE266577_malignant": ("GSE266577_mal", None)}
sig_up = set(open(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\state_sig_up.txt", encoding="utf-8").read().split())
sig_dn = set(open(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\state_sig_dn.txt", encoding="utf-8").read().split())
TF = ["RELA", "NFKB1", "IRF1", "STAT1"]

rows = []
for alpha in [0.1, 1.0, 10.0]:
    for host, (nm, _) in HOSTS.items():
        genes, C = load_tri(nm)
        gidx = {g: i for i, g in enumerate(genes)}
        net_h = net[net["source"].isin(gidx) & net["target"].isin(gidx)]
        # 限本宿主数据拟合（Ridge）
        cpm = C.multiply(1e4 / C.sum(axis=0).A1).tocsr()
        X = np.log1p(cpm.toarray().T)
        del cpm, C
        edges = []
        for tgt, sub in net_h.groupby("target"):
            parents = sub["source"].unique()
            cols = [gidx[p] for p in parents if p in gidx]
            if not cols: continue
            m = Ridge(alpha=alpha)
            try: m.fit(X[:, cols], X[:, gidx[tgt]])
            except Exception: continue
            for p, coef in zip(parents, m.coef_):
                edges.append({"tf": p, "target": tgt, "beta": coef})
        hb = pd.DataFrame(edges)
        del X
        # 分数（每签名版；α 网格只对 200 基线）
        for nsig, (U, D) in {"g200": (sig_up, sig_dn), "g100": (set(list(sig_up)[:100]), set(list(sig_dn)[:100])),
                             "g300": (set(list(sig_up)[:300]), set(list(sig_dn)[:300]))}.items():
            svec = np.zeros(len(genes))
            for g in U:
                if g in gidx: svec[gidx[g]] = 1.0
            for g in D:
                if g in gidx: svec[gidx[g]] = -1.0
            sn = np.linalg.norm(svec) or 1.0
            for tf in TF:
                cfs = hb[hb["tf"] == tf]
                if len(cfs) == 0: continue
                delta = np.zeros(len(genes))
                for _, r in cfs.iterrows():
                    if r["target"] in gidx:
                        delta[gidx[r["target"]]] += r["beta"] * (-1.0)
                sc = float(np.dot(delta, svec) / (np.linalg.norm(delta) * sn + 1e-12))
                rows.append({"alpha": alpha, "host": host, "sig": nsig, "tf": tf, "KO_score": sc})
        print(host, "alpha", alpha, "edges:", len(hb), flush=True)
res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "b9_alpha_grid.csv"), index=False)
piv = res[res["sig"] == "g200"].pivot_table(index=["alpha"], columns=["host", "tf"], values="KO_score")
print(piv.round(3).to_string(), flush=True)
print("\n== sig 敏感性（alpha=1）==")
piv2 = res[res["alpha"] == 1.0].pivot_table(index=["host"], columns=["sig", "tf"], values="KO_score")
print(piv2.round(3).to_string(), flush=True)
print("DONE 40b3")
