# -*- coding: utf-8 -*-
"""40b2 v2: 单 TF 随机对照（1000 随机 TF 的 |score| 零分布）+ 集合统计（抽10取max 1000×）。"""
import numpy as np, pandas as pd, sys, os
sys.stdout.reconfigure(encoding='utf-8')
OUT = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b9"
hb = pd.read_csv(os.path.join(OUT, "grn_ridge.csv"))
TFS = ["RELA", "NFKB1", "IRF1", "IRF3", "STAT1", "STAT3", "IRF7", "HIF1A", "SP1", "JUN"]
sig_up = set(open(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\state_sig_up.txt", encoding="utf-8").read().split())
sig_dn = set(open(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\state_sig_dn.txt", encoding="utf-8").read().split())

def load_genes(nm):
    d = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b9_input"
    return [l.strip() for l in open(os.path.join(d, f"{nm}.genes.tsv"), encoding="utf-8")]

HOSTS = {"GSE293459_all": load_genes("GSE293459_all"),
         "GSE266577_epithelial": load_genes("GSE266577_epi"),
         "GSE266577_malignant": load_genes("GSE266577_mal")}

rng = np.random.default_rng(2026)
rows = []
for host, genes in HOSTS.items():
    h = hb[hb["host"] == host]
    gidx = {g: i for i, g in enumerate(genes)}
    svec = np.zeros(len(genes))
    for g in sig_up:
        if g in gidx: svec[gidx[g]] = 1.0
    for g in sig_dn:
        if g in gidx: svec[gidx[g]] = -1.0
    sn = np.linalg.norm(svec) or 1.0
    coef = {}
    for (tf, tgt), g in h.groupby(["tf", "target"]):
        coef.setdefault(tf, {})[tgt] = float(g["beta"].sum())
    all_tf = sorted(coef.keys())
    def score_one(tf, sign):
        delta = np.zeros(len(genes))
        for tgt, b in coef.get(tf, {}).items():
            if tgt in gidx:
                delta[gidx[tgt]] += b * sign
        return float(np.dot(delta, svec) / (np.linalg.norm(delta) * sn + 1e-12))
    real_tf = [t for t in TFS if t in coef]
    real_s = {tf: score_one(tf, -1.0) for tf in real_tf}
    ## 零分布：单 TF 随机（1000 个）
    pool = [t for t in all_tf if t not in real_tf]
    null1 = np.array([abs(score_one(rng.choice(pool), -1.0)) for _ in range(2000)])
    ## 集合统计：随机抽与真实 TF 同样数量(10)的 TF，取 max|score|（1000×）
    nullm = np.array([max(abs(score_one(t, -1.0)) for t in rng.choice(pool, size=min(10, len(pool)), replace=False))
                      for _ in range(1000)])
    obs_max = max(abs(real_s[tf]) for tf in real_tf)
    p_set = float((nullm >= obs_max).mean()) + 1e-4
    print(f"== {host}: obs_max={obs_max:.4f} nullm95={np.quantile(nullm, .95):.4f} p_set={p_set:.4f}", flush=True)
    for tf in real_tf:
        p = float((null1 >= abs(real_s[tf])).mean()) + 1e-4
        print(f"   {tf}: KO {real_s[tf]:+.4f}  p(perm)={p:.4f}", flush=True)
        rows.append({"host": host, "tf": tf, "KO": real_s[tf], "p_perm": p, "p_set": p_set})
    # 组合
    combos = [("NFKB1", "IRF1"), ("NFKB1", "STAT1"), ("IRF1", "STAT1"), ("NFKB1", "IRF1", "STAT1")]
    for cb in combos:
        ok = all(t in coef for t in cb)
        if ok:
            delta = np.zeros(len(genes))
            for tf in cb:
                for tgt, b in coef[tf].items():
                    if tgt in gidx: delta[gidx[tgt]] += b * (-1.0)
            sc = float(np.dot(delta, svec) / (np.linalg.norm(delta) * sn + 1e-12))
            print(f"   combo {'+'.join(cb)} KO: {sc:+.4f}", flush=True)
            rows.append({"host": host, "tf": "+".join(cb) + "_KO", "KO": sc, "p_perm": None, "p_set": p_set})
res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "b9_summary_perm_combo.csv"), index=False)
print("DONE 40b2 v2")
