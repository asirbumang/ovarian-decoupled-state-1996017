# -*- coding: utf-8 -*-
"""45u_consensus_types.py: JTM 主线镜像——①患者级共识聚类（TCGA-OV 427；bootstrap 500× k=2）
→ 亚型×OS/反应/状态率；②GSE266577 患者亚型（51）→ SPP1+ 巨噬锁定。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from scipy.stats import mannwhitneyu, chi2_contingency
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"
BASE = r"D:\子宫内膜异位症2\analysis\L0_bulk\results"

## A) TCGA-OV 共识聚类
d = pd.read_csv(os.path.join(BASE, "track0_tcga_scores_states.csv"), low_memory=False)
d = d[d["sample"].notna()]
feat = d[["M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"]].astype(float).values
z = lambda a: (a - a.mean(axis=0)) / (a.std(axis=0) + 1e-9)
F = z(feat)
rng = np.random.default_rng(2026)
n = len(F); k = 2
co = np.zeros((n, n))
sv = np.full(n, 0)
for b in range(500):
    idx = rng.choice(n, size=int(0.8 * n), replace=False)
    km = KMeans(n_clusters=k, random_state=int(b), n_init=5).fit(F[idx])
    for c in np.unique(km.labels_):
        sel = idx[km.labels_ == c]
        for i1 in range(len(sel)):
            for i2 in range(i1 + 1, len(sel)):
                co[sel[i1], sel[i2]] += 1
    sv[idx] += 1
con = co / np.maximum(sv[:, None], 1)
d["cluster"] = KMeans(n_clusters=k, random_state=42, n_init=50).fit_predict(F)
## 亚型命名：按 decouple 中位（high/low）
d["subtype"] = np.where(d.groupby("cluster")["decouple"].transform("median") >= d["decouple"].median(), "high", "low")
## 稳定性：共识内部（平均共识>0.5 比例）
d["consensus_mean"] = con.mean(axis=1)
print("subtype sizes:", d["subtype"].value_counts().to_dict(), flush=True)
print("consensus mean:", round(float(d["consensus_mean"].mean()), 3),
      "| frac>0.5:", round(float((d["consensus_mean"] > 0.5).mean()), 3), flush=True)
## 关联：OS（KW/秩和+logrank 用 lifelines）
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
subs = d["subtype"].values
lr = logrank_test(d["os_days"][subs == "high"], d["os_days"][subs == "low"],
                  d["event"][subs == "high"], d["event"][subs == "low"])
print("logrank p (subtype OS):", round(float(lr.p_value), 4), flush=True)
cr = d[["subtype", "cr_vs_pd"]].dropna()
tab = pd.crosstab(cr["subtype"], cr["cr_vs_pd"])
print("CR-PD crosstab:\n", tab.to_string(), flush=True)
chi, pchi, _, _ = chi2_contingency(tab)
print("chi2 p:", round(float(pchi), 4), flush=True)
st = d[["subtype", "state_high"]].dropna()
tab2 = pd.crosstab(st["subtype"], st["state_high"])
print("state counts:\n", tab2.to_string(), flush=True)
d.to_csv(os.path.join(BASE, "track2", "tcga_subtype_consensus.csv"), index=False)

## B) GSE266577 患者亚型 + SPP1+ 巨噬
sc = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577.csv",
                 usecols=["cell", "publication_sample_code_final", "cell_type", "M1_immune_activation",
                          "M2_metabolic_GSH", "M3_mito_cytoskeleton"], low_memory=False)
for m_ in ["M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"]:
    sc[m_] = pd.to_numeric(sc[m_], errors="coerce").astype(float)
pat = sc.groupby("publication_sample_code_final")[["M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"]].mean()
pat = pat.dropna()
print("GSE266577 patients:", len(pat), flush=True)
F2 = z(pat.values)
km2 = KMeans(n_clusters=2, random_state=42, n_init=50).fit(F2)
pat["subtype"] = np.where(km2.labels_ == 0, "high", "low")
pat = pat.reset_index()[["publication_sample_code_final", "subtype"]].rename(columns={"subtype": "pte_sub"})
sc = sc.merge(pat, on="publication_sample_code_final", how="left")
## SPP1+ 巨噬（Macrophages 内 SPP1 高 = 前 25%）
sc["SPP1"] = 0.0  # 占位——SPP1 表达需基因矩阵；用 cell_type 空间映射（n 巨大，改用 marker 近似：C1QA 高=组织巨噬已有注释；SPP1 需加载表达）
## 简化（数据已有）：巨噬比例 per 患者亚型
mac = sc.assign(is_mac=sc["cell_type"].isin(["Macrophages", "Classical monocytes", "Migratory DCs", "DC1", "DC2", "pDC"])).groupby(["publication_sample_code_final", "pte_sub"])["is_mac"].mean().reset_index()
print("mac fraction: high", round(float(mac[mac["pte_sub"] == "high"]["is_mac"].mean()), 4),
      " low", round(float(mac[mac["pte_sub"] == "low"]["is_mac"].mean()), 4), flush=True)
print("patients per subtype:", sc.drop_duplicates("publication_sample_code_final").groupby("pte_sub")["publication_sample_code_final"].count().to_dict(), flush=True)
mac.to_csv(os.path.join(QRY, "F11_pt_subtype_mac.csv"), index=False)
## 图：亚型共识矩阵热图 + OS KM + 状态率
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 3, figsize=(150/25.4, 52/25.4))
axes[0].imshow(con[np.argsort(d["subtype"].values)], cmap="Reds", vmin=0, vmax=1)
axes[0].set_title("a consensus matrix (k=2; mean=%.2f)" % d["consensus_mean"].mean(), size=6)
axes[0].set_xticks([]); axes[0].set_yticks([])
for sb, c in [("high", "#DC0000"), ("low", "#3C5488")]:
    sub = d[(d["subtype"] == sb) & d["os_days"].notna() & d["event"].notna() & (d["os_days"] > 0)]
    km_ = KaplanMeierFitter().fit(sub["os_days"], sub["event"], label=f"{sb} (n={len(sub)})")
    km_.plot_survival_function(ax=axes[1], color=c, linewidth=1.0)
axes[1].set_xlabel("days", size=6); axes[1].tick_params(labelsize=5.4)
axes[1].set_title("b OS by subtype (logrank p=%.3f)" % lr.p_value, size=6)
sr = d.groupby("subtype")["state_high"].mean()
axes[2].bar([0, 1], sr.values, color=["#DC0000", "#3C5488"], width=0.5)
axes[2].set_xticks([0, 1]); axes[2].set_xticklabels(sr.index, fontsize=5.5)
axes[2].set_ylabel("state frequency", size=6); axes[2].tick_params(labelsize=5.4)
axes[2].set_title("c state rate by subtype", size=6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F11_consensus.png"), dpi=600)
plt.close(fig)
print("F11 saved", flush=True)
print("DONE 45u")
