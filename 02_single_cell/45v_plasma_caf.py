# -*- coding: utf-8 -*-
"""45v_plasma_caf.py: I5 创新执行——浆细胞 & CAF 亚群刻画（方法借 JTM：亚群划分→marker→与状态/亚型关联；
对象=我们的富集群）。①全量矩阵重载 ②两谱系各 kmeans k=2 ③marker+state%+患者亚型关联 ④F12 面板。"""
import numpy as np, pandas as pd, sys, os
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from scipy.stats import mannwhitneyu
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"
d = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\umap_input"
genes = [l.strip() for l in open(os.path.join(d, "GSE266577_full.genes.tsv"), encoding="utf-8")]
cells = [l.strip() for l in open(os.path.join(d, "GSE266577_full.cells.tsv"), encoding="utf-8")]
dt = pd.read_csv(os.path.join(d, "GSE266577_full.triplet.gz"))
gidx = {g: i for i, g in enumerate(genes)}
sidx = {c: i for i, c in enumerate(cells)}
C = sp.csr_matrix((dt["v"].to_numpy(), (dt["g"].map(gidx).to_numpy(), dt["c"].map(sidx).to_numpy())),
                  shape=(len(genes), len(cells)))
del dt
## 谱系标签 + state + 患者亚型
sc = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577.csv",
                 usecols=["cell", "cell_type", "publication_sample_code_final", "M1_immune_activation",
                          "M2_metabolic_GSH", "M3_mito_cytoskeleton"], low_memory=False)
for m_ in ["M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"]:
    sc[m_] = pd.to_numeric(sc[m_], errors="coerce").astype(float)
z = lambda x: (x - x.mean()) / (x.std() + 1e-9)
sc["state"] = (z(sc["M1_immune_activation"]) >= np.quantile(z(sc["M1_immune_activation"]), .75)) & \
              (z(sc["M2_metabolic_GSH"]) <= np.quantile(z(sc["M2_metabolic_GSH"]), .25)) & \
              (z(sc["M3_mito_cytoskeleton"]) <= np.quantile(z(sc["M3_mito_cytoskeleton"]), .25))
sub = pd.read_csv(os.path.join(QRY, "F11_pt_subtype_mac.csv"))
## 患者亚型表（从 45u 落盘：cl = sc groupby -> 重新由 sc 聚合即可；直接再算）
if os.path.exists(os.path.join(QRY, "F11_pt_subtype_mac.csv")):
    psm = pd.read_csv(os.path.join(QRY, "F11_pt_subtype_mac.csv"))
else:
    psm = None
meta = sc.set_index("cell")
ct_map = dict(zip(sc["cell"], sc["cell_type"]))
st_map = dict(zip(sc["cell"], sc["state"]))
res_all = []
for LIN, key in [("Plasma", "Plasma cells"), ("CAF", "Fibroblasts")]:
    sel = [c for c in cells if ct_map.get(c) == key]
    if len(sel) < 300:
        print("skip", key, len(sel), flush=True); continue
    Xc = np.log1p(C[:, [sidx[c] for c in sel]].multiply(1e4 / C[:, [sidx[c] for c in sel]].sum(axis=0).A1).toarray().T)
    v = Xc.var(axis=0); hvg = np.argsort(-v)[:1500]
    P = PCA(n_components=20, random_state=42).fit_transform(Xc[:, hvg])
    lab = KMeans(n_clusters=2, random_state=42, n_init=20).fit_predict(P)
    print(key, "n=", len(sel), "sizes:", pd.Series(lab).value_counts().to_dict(), flush=True)
    ps = pd.DataFrame({"cell": sel, "sub": lab})
    ps["state"] = [st_map.get(c, False) for c in sel]
    print("  state% by sub:", ps.groupby("sub")["state"].agg(["mean", "count"]).round(3).to_dict(), flush=True)
    ## marker top6
    outm = []
    for s_ in np.unique(lab):
        a = Xc[lab == s_][:, hvg]; b = Xc[lab != s_][:, hvg]
        lf = a.mean(axis=0) - b.mean(axis=0)
        top = np.argsort(-lf)[:6]
        outm.append({"sub": s_, "markers": ", ".join(genes[hvg[j]] for j in top)})
    print(pd.DataFrame(outm).to_string(index=False), flush=True)
    res_all.append(pd.DataFrame({"cell": sel, "lineage": key, "sub": lab,
                                 "state": [st_map.get(c, False) for c in sel]}))
res = pd.concat(res_all)
res.to_csv(os.path.join(QRY, "F12_plasma_caf_sub.csv"), index=False)
## F12 图（每谱系：子群占比×state% 双条）
fig, axes = plt.subplots(1, 2, figsize=(110/25.4, 46/25.4))
for ax, key in zip(axes, ["Plasma cells", "Fibroblasts"]):
    r = res[res["lineage"] == key]
    gs = r.groupby("sub")["state"].mean()
    ns = r.groupby("sub")["state"].count()
    ax.bar(range(len(gs)), gs.values * 100, color=["#3C5488", "#E64B35"], width=0.5)
    for i, (g_, n_) in enumerate(zip(gs.values, ns.values)):
        ax.text(i, g_ * 100 + 0.3, f"n={n_}", ha="center", fontsize=5.4)
    ax.set_xticks(range(len(gs))); ax.set_xticklabels([f"{key} sub{s}" for s in gs.index], fontsize=5.4)
    ax.set_ylabel("state fraction (%)", size=6); ax.tick_params(labelsize=5.4)
    ax.set_title(f"{key} subclusters (state enrichment)", size=6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F12_plasma_caf.png"), dpi=600)
plt.close(fig)
print("F12 saved", flush=True)
print("DONE 45v")
