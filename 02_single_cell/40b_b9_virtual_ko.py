# -*- coding: utf-8 -*-
"""40b-b9 同构版（卵巢癌）：CollecTRI base GRN + Ridge host GRN + 一步传播 → KO/OE 分数。
hosts: GSE293459 (persister, 全部细胞) / GSE266577 (上皮+恶性肿瘤细胞, 有注释)
签名 sig = 解耦态 vs 其余 的 logFC 向量（GSE266577 上皮内）+ 反向（逆转=返回代谢偶联）。
"""
import sys, os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import scipy.sparse as sp

OUT = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b9"
os.makedirs(OUT, exist_ok=True)
HOSTS = ["GSE293459_all", "GSE266577_epithelial", "GSE266577_malignant"]
TFS = ["RELA", "NFKB1", "IRF1", "IRF3", "STAT1", "STAT3", "IRF7", "HIF1A", "SP1", "JUN"]
ALPHA = 1.0

print("== CollecTRI ==", flush=True)
import decoupler as dc
net = dc.op.collectri()
net = net[["source", "target"]].drop_duplicates()
print("collectri edges:", len(net), flush=True)

## ---- 每个 host 构建 X（cells × genes, log1p CP1e4）----
def load_293459():
    import gzip
    d = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b9_input"
    genes = [l.strip() for l in open(d + r"\GSE293459_all.genes.tsv", encoding="utf-8")]
    dt = pd.read_csv(d + r"\GSE293459_all.triplet.gz")
    gidx = {g: i for i, g in enumerate(genes)}
    cells = dt["c"].unique()
    sidx = {v: i for i, v in enumerate(cells)}
    C = sp.csr_matrix((dt["v"].to_numpy(), (dt["g"].map(gidx).to_numpy(), dt["c"].map(sidx).to_numpy())),
                      shape=(len(genes), len(cells)))
    return genes, C

def load_266577(subset_lineage=None, subset_malign=None):
    ## 快速路径：R 已导出紧凑 triplet（40b_prep）
    def load_tri(nm):
        d = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b9_input"
        genes = [l.strip() for l in open(d + rf"\GSE266577_{nm}.genes.tsv", encoding="utf-8")]
        dt = pd.read_csv(d + rf"\GSE266577_{nm}.triplet.gz")
        gidx = {g: i for i, g in enumerate(genes)}
        cells = dt["c"].unique()
        sidx = {v: i for i, v in enumerate(cells)}
        C = sp.csr_matrix((dt["v"].to_numpy(), (dt["g"].map(gidx).to_numpy(), dt["c"].map(sidx).to_numpy())),
                          shape=(len(genes), len(cells)))
        return genes, C
    g2, Cepi = load_tri("epi")
    _, Cmal = load_tri("mal")
    bars2 = None
    return g2, Cepi, Cmal

print("== loading GSE293459 ==", flush=True)
genes1, C1 = load_293459()
cpm1 = C1.multiply(1e4 / C1.sum(axis=0).A1).tocsr()
print("GSE293459 host cells:", cpm1.shape[1], flush=True)

print("== loading GSE266577 ==", flush=True)
genes2, Cepi, Cmal = load_266577()
print("GSE266577 epi cells:", Cepi.shape[1], " mal cells:", Cmal.shape[1], flush=True)
cpm2 = Cepi.multiply(1e4 / Cepi.sum(axis=0).A1).tocsr()
cpm3 = Cmal.multiply(1e4 / Cmal.sum(axis=0).A1).tocsr()

HOSTS_DATA = [("GSE293459_all", genes1, cpm1), ("GSE266577_epithelial", genes2, cpm2),
              ("GSE266577_malignant", genes2, cpm3)]
## 签名向量：解耦态 logFC（来自 GSE266577 上皮 top/bottom quartile DE，已算过→ 用 state_sig 向量近似：
## sig = +1 for up genes, -1 for dn genes（二值向量；cos 对方向敏感）
sig_up = set(open(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\state_sig_up.txt", encoding="utf-8").read().split())
sig_dn = set(open(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\state_sig_dn.txt", encoding="utf-8").read().split())

rows = []
for host, genes, cpm in HOSTS_DATA:
    gidx = {g: i for i, g in enumerate(genes)}
    net_h = net[net["source"].isin(gidx) & net["target"].isin(gidx)]
    print(host, "genes:", len(genes), "cells:", cpm.shape[1], "edges:", len(net_h), flush=True)
    X = np.log1p(cpm.toarray().T)  # cells × genes
    del cpm; 
    for tgt, sub in net_h.groupby("target"):
        parents = sub["source"].unique()
        cols = [gidx[p] for p in parents if p in gidx]
        if not cols:
            continue
        Xp = X[:, cols]
        y = X[:, gidx[tgt]]
        from sklearn.linear_model import Ridge
        m = Ridge(alpha=ALPHA)
        try:
            m.fit(Xp, y)
        except Exception:
            continue
        for p, coef in zip(parents, m.coef_):
            rows.append({"host": host, "tf": p, "target": tgt, "beta": coef})
    del X; 
hb = pd.DataFrame(rows)
hb.to_csv(os.path.join(OUT, "grn_ridge.csv"), index=False)
print("GRN edges:", len(hb), flush=True)

## ---- 扰动分数 ----
sig_vec = np.zeros(len(genes1))
tbl = []
for host, genes, _ in HOSTS_DATA:
    gidx = {g: i for i, g in enumerate(genes)}
    svec = np.zeros(len(genes))
    for g in sig_up:
        if g in gidx: svec[gidx[g]] = 1.0
    for g in sig_dn:
        if g in gidx: svec[gidx[g]] = -1.0
    sn = np.linalg.norm(svec) or 1.0
    hb_h = hb[hb["host"] == host]
    for tf in TFS:
        coefs = hb_h[hb_h["tf"] == tf]
        if len(coefs) == 0:
            continue
        for mode, sign in [("KO", -1.0), ("OE", 1.0)]:
            delta = np.zeros(len(genes))
            for _, r in coefs.iterrows():
                if r["target"] in gidx:
                    delta[gidx[r["target"]]] += r["beta"] * sign * 1.0
            score = float(np.dot(delta, svec) / (np.linalg.norm(delta) * sn + 1e-12))
            tbl.append({"host": host, "tf": tf, "mode": mode, "pert_score": score, "n_edges": len(coefs)})
            print(host, tf, mode, "score=%.4f" % score, flush=True)
summ = pd.DataFrame(tbl)
summ.to_csv(os.path.join(OUT, "b9_summary.csv"), index=False)

## ---- 热图 ----
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42})
piv = summ.pivot_table(index=["tf", "mode"], columns="host", values="pert_score")
fig, ax = plt.subplots(figsize=(150 / 25.4, 100 / 25.4))
vmax = np.abs(piv.values).max()
im = im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns, fontsize=6.5, rotation=30, ha="right")
ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels([f"{i}|{m}" for i, m in piv.index], fontsize=6.5)
ax.set_title("Virtual perturbation score (cos ΔX vs mitoxy signature)\n>0 mimic decoupled | <0 revert", fontsize=7)
plt.colorbar(im, ax=ax, fraction=0.03)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "F_B9a_score_heatmap.png"), dpi=600)
plt.savefig(os.path.join(OUT, "F_B9a_score_heatmap.svg"))
plt.close()
print("DONE b9")
