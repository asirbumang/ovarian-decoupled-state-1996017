# -*- coding: utf-8 -*-
"""40d: 基因集自证 B——随机 60 基因集零分布（1000×）：
① cor(随机M1,随机M2)/cor(随机M1,随机M3) 与真实（负）对照；② cor(随机decouple, 真decouple)。"""
import numpy as np, pandas as pd, sys, os
import scipy.sparse as sp
from scipy.stats import spearmanr
sys.stdout.reconfigure(encoding='utf-8')
d = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b9_input"
genes = [l.strip() for l in open(os.path.join(d, "GSE266577_epi.genes.tsv"), encoding="utf-8")]
dt = pd.read_csv(os.path.join(d, "GSE266577_epi.triplet.gz"))
gidx = {g: i for i, g in enumerate(genes)}
cells = dt["c"].unique(); sidx = {v: i for i, v in enumerate(cells)}
C = sp.csr_matrix((dt["v"].to_numpy(), (dt["g"].map(gidx).to_numpy(), dt["c"].map(sidx).to_numpy())),
                  shape=(len(genes), len(cells)))
print("matrix:", C.shape, flush=True)
cpm = C.multiply(1e4 / C.sum(axis=0).A1).tocsr()
X = np.log1p(cpm.toarray().T)  # cells × genes
col_of = {g: i for i, g in enumerate(genes)}
mods = {}
for line in open(r"D:\子宫内膜异位症2\pipeline\mitoxy_four_modules.gmt", encoding="utf-8"):
    p = line.strip().split("\t")
    if p[0].startswith("MITOXY"):
        mods[p[0]] = [g for g in p[2:] if g in col_of]
m1 = np.mean(X[:, [col_of[g] for g in mods["MITOXY_M1_immune_activation"]]], axis=1)
m2 = np.mean(X[:, [col_of[g] for g in mods["MITOXY_M2_metabolic_GSH"]]], axis=1)
m3 = np.mean(X[:, [col_of[g] for g in mods["MITOXY_M3_mito_cytoskeleton"]]], axis=1)
z = lambda a: (a - a.mean()) / (a.std() + 1e-9)
true_dec = z(m1) - (z(m2) + z(m3)) / 2
true_rho12, true_rho13 = spearmanr(m1, m2)[0], spearmanr(m1, m3)[0]
print(f"TRUE: rho(M1,M2)={true_rho12:.3f}  rho(M1,M3)={true_rho13:.3f}", flush=True)
pool = [g for g in genes if g in col_of]
pool = [g for g in pool if g not in set().union(*[set(v) for v in mods.values()])]
rng = np.random.default_rng(2026)
n = 1000
null12, null13, null_corr = [], [], []
for _ in range(n):
    s = rng.choice(pool, size=60, replace=False)
    a, b, cc = s[:20], s[20:40], s[40:]
    v1 = X[:, [col_of[g] for g in a]].mean(axis=1)
    v2 = X[:, [col_of[g] for g in b]].mean(axis=1)
    v3 = X[:, [col_of[g] for g in cc]].mean(axis=1)
    null12.append(spearmanr(v1, v2)[0])
    null13.append(spearmanr(v1, v3)[0])
    rd = z(v1) - (z(v2) + z(v3)) / 2
    null_corr.append(spearmanr(rd, true_dec)[0])
null12, null13, null_corr = map(np.array, (null12, null13, null_corr))
p12 = float((np.abs(null12) >= abs(true_rho12)).mean())
p13 = float((np.abs(null13) >= abs(true_rho13)).mean())
print(f"NULL: rho12 mean={null12.mean():+.4f} sd={null12.std():.4f}  p={p12:.6f}", flush=True)
print(f"NULL: rho13 mean={null13.mean():+.4f} sd={null13.std():.4f}  p={p13:.6f}", flush=True)
print(f"NULL: cor(rand_dec, true_dec) mean={null_corr.mean():+.4f} sd={null_corr.std():.4f}", flush=True)
out = pd.DataFrame({"null12": null12, "null13": null13, "null_corr": null_corr,
                    "true12": true_rho12, "true13": true_rho13})
out.to_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\geneset_random_null.csv", index=False)
print("DONE 40d")
