# -*- coding: utf-8 -*-
"""45m v2: Scissor 等价 → 相关矩阵 + Cox-Lasso + 阈值 + 重叠 + S9。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from sksurv.linear_model import CoxnetSurvivalAnalysis
sys.stdout.reconfigure(encoding='utf-8')
BASE = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\track2\mime"
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs\supp"

## scRNA 60k×2900（基因×细胞）→ 细胞×基因
sc = pd.read_csv(os.path.join(BASE, "scissor_expr60k.csv.gz"), index_col=0)
cells = pd.read_csv(os.path.join(BASE, "scissor_cells60k.csv"))["cell"].tolist()
Xc = sc.T.values.astype(np.float32)  # cells x genes
genes = sc.index.tolist()
print("sc:", Xc.shape, flush=True)
## bulk 基因×425
bt = pd.read_csv(os.path.join(BASE, "tcga_exp.csv.gz"), index_col=0)
cl = pd.read_csv(os.path.join(BASE, "tcga_clin.csv"))
cl = cl[cl["sample"].isin(bt.columns)]
bt = bt[cl["sample"].tolist()]
common = [g for g in genes if g in bt.index]
print("common genes:", len(common), flush=True)
Xc = Xc[:, [genes.index(g) for g in common]]
Xb = bt.loc[common].values.astype(np.float32)  # genes x patients
## 相关矩阵 R = cells × patients（跨基因 pearson）
def zrow(a):
    m = a.mean(axis=1, keepdims=True); s = a.std(axis=1, keepdims=True) + 1e-9
    return (a - m) / s
Zc = zrow(Xc)                       # cells x genes
Zb = zrow(Xb.T)                     # patients x genes
R = (Zc @ Zb.T) / (len(common) - 1)  # cells x patients
R = np.nan_to_num(R.astype(np.float32))
print("R:", R.shape, flush=True)
## Cox-Lasso：特征=细胞（60k），样本=患者（425）
T = cl["futime"].values.astype(float); E = cl["fustat"].values.astype(int)
y = np.array([(e, t) for e, t in zip(E, T)], dtype=[("event", bool), ("time", float)])
X = R.T  # patients x cells
mdl = CoxnetSurvivalAnalysis(alpha_min_ratio=0.05, n_alphas=25, max_iter=1000, tol=1e-6)
alphas = np.geomspace(0.01, 1.0, 25)
mdl.fit(X, y)
coefs = np.array([c for c in mdl.coef_])  # alphas x cells
print("coefs:", coefs.shape, flush=True)
## 稳定非零：coef_ = (cells, alphas)；跨 alpha 非零率≥50%（轴 1）
nz = (coefs != 0).sum(axis=1)
stable = nz >= coefs.shape[1] * 0.5
mcoef = coefs.mean(axis=1)
mcoef[~stable] = 0.0
idx_pos = np.argsort(-mcoef)[: int(0.2 * len(mcoef))]
idx_neg = np.argsort(mcoef)[: int(0.2 * len(mcoef))]
pos = np.zeros(len(mcoef), bool); neg = np.zeros(len(mcoef), bool)
pos[idx_pos] = True; neg[idx_neg] = True
print("Scissor+:", int(pos.sum()), "Scissor-:", int(neg.sum()), flush=True)
print("lens:", len(cells), len(mcoef), len(pos), len(neg), flush=True)
res = pd.DataFrame({"cell": np.array(cells), "coef": mcoef, "scissor_pos": pos, "scissor_neg": neg})
res.to_csv(os.path.join(BASE, "scissor_equiv_coefs.csv"), index=False)

## 与 state 定义重叠（GSE266577 全细胞四分位；本次采样）
scores = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577.csv",
                     usecols=["cell", "M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"], low_memory=False)
for m_ in ["M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"]:
    scores[m_] = pd.to_numeric(scores[m_], errors="coerce").astype(float)
z = lambda x: (x - x.mean()) / (x.std() + 1e-9)
dz = scores.set_index("cell")
z1, z2, z3 = z(dz["M1_immune_activation"]), z(dz["M2_metabolic_GSH"]), z(dz["M3_mito_cytoskeleton"])
dz["state"] = (z1 >= np.quantile(z1, .75)) & (z2 <= np.quantile(z2, .25)) & (z3 <= np.quantile(z3, .25))
dz = dz.reset_index()
m = res.merge(dz[["cell", "state"]], on="cell")
a = int(((m["scissor_pos"]) & (m["state"])).sum()); b = int((m["scissor_pos"] & ~m["state"]).sum())
c_ = int((~m["scissor_pos"] & m["state"]).sum()); d_ = int((~m["scissor_pos"] & ~m["state"]).sum())
orr, p = fisher_exact([[a, b], [c_, d_]])
print("Scissor+ vs state: OR=%.2f p=%.2e (a=%d b=%d c=%d d=%d)" % (orr, p, a, b, c_, d_), flush=True)
## S9 图：UMAP 叠加（Scissor+ 橙 / Scissor- 蓝）
umo = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\umap_input\umap_full_266577.csv")
mu = umo.merge(m[["cell", "scissor_pos", "scissor_neg"]], on="cell", how="left").fillna({"scissor_pos": False, "scissor_neg": False})
fig, ax = plt.subplots(figsize=(70/25.4, 64/25.4))
ax.scatter(mu["umap1"], mu["umap2"], s=0.04, c="#BBBBBB", alpha=0.12, linewidths=0)
sp = mu[mu["scissor_pos"]]; sn = mu[mu["scissor_neg"]]
ax.scatter(sn["umap1"], sn["umap2"], s=0.3, c="#3C5488", alpha=0.8, linewidths=0, label=f"Scissor- (n={len(sn)})")
ax.scatter(sp["umap1"], sp["umap2"], s=0.3, c="#DC0000", alpha=0.8, linewidths=0, label=f"Scissor+ (n={len(sp)})")
ax.set_xticks([]); ax.set_yticks([])
for s_ in ax.spines.values():
    s_.set_visible(False)
leg = ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=5.5, frameon=False)
for h in leg.legend_handles:
    h.set_sizes([8])
fig.tight_layout()
fig.savefig(os.path.join(QRY, "S9_scissor.png"), dpi=600)
plt.close(fig)
print("S9 saved", flush=True)
print("DONE 45m")
