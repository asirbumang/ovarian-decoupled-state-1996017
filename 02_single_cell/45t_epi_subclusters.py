# -*- coding: utf-8 -*-
"""45t_epi_subclusters.py: 上皮内无监督亚群（JTM 同构）——HVG2000→PCA30→kmeans k=4；
UMAP(n=30/min_dist 0.3)；每亚群：CNV/state/模块/marker；输出 F3f。"""
import numpy as np, pandas as pd, sys, os
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"
d = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b9_input"
genes = [l.strip() for l in open(os.path.join(d, "GSE266577_epi.genes.tsv"), encoding="utf-8")]
dt = pd.read_csv(os.path.join(d, "GSE266577_epi.triplet.gz"))
gidx = {g: i for i, g in enumerate(genes)}
cells = dt["c"].unique(); sidx = {v: i for i, v in enumerate(cells)}
C = sp.csr_matrix((dt["v"].to_numpy(), (dt["g"].map(gidx).to_numpy(), dt["c"].map(sidx).to_numpy())),
                  shape=(len(genes), len(cells)))
cpm = C.multiply(1e4 / C.sum(axis=0).A1).tocsr()
X = np.log1p(cpm.toarray().T)  # cells x genes
print("X:", X.shape, flush=True)
## HVG2000（方差）
var = X.var(axis=0)
hvg = np.argsort(-var)[:2000]
Xh = X[:, hvg]
## PCA30
pca = PCA(n_components=30, random_state=42)
P = pca.fit_transform(Xh)
k = KMeans(n_clusters=4, random_state=42, n_init=10).fit(P)
lab = k.labels_
print("cluster sizes:", pd.Series(lab).value_counts().to_dict(), flush=True)
## UMAP
try:
    import umap as umap_
    emb = umap_.UMAP(n_neighbors=30, min_dist=0.3, random_state=42).fit_transform(P)
except Exception:
    from sklearn.manifold import TSNE
    emb = TSNE(n_components=2, random_state=42, perplexity=30).fit_transform(P[:20000])
np.save(os.path.join(QRY, "F3f_epi_umap.npy"), emb)
np.save(os.path.join(QRY, "F3f_epi_lab.npy"), lab)
## 每亚群统计
cnv = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\percell_cnv_score_v2.csv",
                  usecols=["cell", "cnv_top", "malignant_like", "group"])
cnv = cnv[cnv["group"].astype(str) == "Epithelial"]
sc = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577_epithelial.csv",
                 usecols=["cell", "decouple", "M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"], low_memory=False)
mm = pd.DataFrame({"cell": cells, "lab": lab})
mm = mm.merge(cnv, on="cell", how="left").merge(sc, on="cell", how="left")
gs = mm.groupby("lab").agg(n=("cell", "count"), cnv=("cnv_top", "median"),
                            mal=("malignant_like", "mean"), dec=("decouple", "median"))
print(gs.round(3).to_string(), flush=True)
mm.to_csv(os.path.join(QRY, "F3f_epi_subcluster.csv"), index=False)
## 每亚群 marker（top2 平均 logFC）
out = []
for c0 in sorted(mm["lab"].unique()):
    a = X[mm["lab"] == c0]; b = X[mm["lab"] != c0]
    lf = a.mean(axis=0) - b.mean(axis=0)
    top = np.argsort(-lf)[:6]
    out.append({"lab": c0, "genes": ", ".join(genes[hvg[i]] for i in top)})
print(pd.DataFrame(out).to_string(index=False), flush=True)
print("DONE 45t")
