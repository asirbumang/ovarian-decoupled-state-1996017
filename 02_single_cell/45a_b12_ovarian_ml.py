# -*- coding: utf-8 -*-
"""45a_b12_ovarian_ml.py: 卵巢癌 ML（镜像母版 B12 三段式）。
A) 细胞级动物园（10FS×6 分类器 × 5 折 GroupKFold(样本)；细胞=top2000 状态 DE 基因；AUC）
B) bootstrap L1 最小签名 → 频率
C) 简化说明：bulk 外验单独 R 脚本（45b）
"""
import numpy as np, pandas as pd, sys, os, warnings
warnings.filterwarnings("ignore")
import scipy.sparse as sp
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.feature_selection import SelectFromModel, SelectKBest, mutual_info_classif
from sklearn.decomposition import PCA
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
from sklearn.dummy import DummyClassifier
sys.stdout.reconfigure(encoding='utf-8')

OUT = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b12"
os.makedirs(OUT, exist_ok=True)
d = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b9_input"
genes = [l.strip() for l in open(os.path.join(d, "GSE266577_epi.genes.tsv"), encoding="utf-8")]
dt = pd.read_csv(os.path.join(d, "GSE266577_epi.triplet.gz"))
gidx = {g: i for i, g in enumerate(genes)}
cells = dt["c"].unique(); sidx = {v: i for i, v in enumerate(cells)}
C = sp.csr_matrix((dt["v"].to_numpy(), (dt["g"].map(gidx).to_numpy(), dt["c"].map(sidx).to_numpy())),
                  shape=(len(genes), len(cells)))
sc = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577.csv",
                 usecols=["cell", "M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton",
                          "publication_sample_code_final"], low_memory=False)
for m in ["M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"]:
    sc[m] = pd.to_numeric(sc[m], errors="coerce").astype(float)
z = lambda x: (x - x.mean()) / (x.std() + 1e-9)
z1, z2, z3 = z(sc["M1_immune_activation"]), z(sc["M2_metabolic_GSH"]), z(sc["M3_mito_cytoskeleton"])
sc["state"] = (z1 >= np.quantile(z1, .75)) & (z2 <= np.quantile(z2, .25)) & (z3 <= np.quantile(z3, .25))
## 匹配细胞
cell_info = sc[["cell", "state", "publication_sample_code_final"]].set_index("cell")
keep = [c for c in cells if c in cell_info.index]
idx_ = [sidx[c] for c in keep]
Ck = C[:, idx_].tocsr()
cpm = Ck.multiply(1e4 / Ck.sum(axis=0).A1).tocsr()
X = np.log1p(cpm.toarray().T)  # cells × genes
y = cell_info.loc[keep, "state"].values.astype(int)
grp = cell_info.loc[keep, "publication_sample_code_final"].values
del C, Ck, cpm
## 特征：top2000 状态 DE 简单口径（state vs rest 的 mean 差绝对 top2000）
ms = X[y == 1].mean(axis=0) - X[y == 0].mean(axis=0)
top = np.argsort(-np.abs(ms))[:2000]
Xf = X[:, top]
print("X:", Xf.shape, "state:", y.sum(), "groups:", len(np.unique(grp)), flush=True)
## 抽 6 万（保持 state/non 比例、按组抽取）
rng = np.random.default_rng(2026)
per_g = 60000 // len(np.unique(grp))
sel = []
for g in np.unique(grp):
    ii = np.where(grp == g)[0]
    sel.extend(rng.choice(ii, size=min(per_g, len(ii)), replace=False))
sel = np.array(sel)
Xf, y, grp = Xf[sel], y[sel], grp[sel]
print("sub:", Xf.shape, "state:", y.sum(), "groups:", len(np.unique(grp)), flush=True)

CLFS = {
    "LR_l1": lambda: LogisticRegression(penalty="l1", solver="liblinear", C=0.5, class_weight="balanced", max_iter=3000),
    "LR_l2": lambda: LogisticRegression(penalty="l2", C=0.5, class_weight="balanced", max_iter=3000),
    "SVM": lambda: SVC(kernel="rbf", C=0.5, probability=True, class_weight="balanced"),
    "RF": lambda: RandomForestClassifier(n_estimators=200, max_depth=8, random_state=42, n_jobs=8),
    "ET": lambda: ExtraTreesClassifier(n_estimators=200, max_depth=8, random_state=42, n_jobs=8),
    "GBM": lambda: GradientBoostingClassifier(n_estimators=150, max_depth=3, random_state=42),
}
SELS = {
    "None": lambda: None,
    "L1": lambda: SelectFromModel(LogisticRegression(penalty="l1", solver="liblinear", C=0.3, class_weight="balanced", max_iter=3000), threshold="median"),
    "RFimp": lambda: SelectFromModel(RandomForestClassifier(n_estimators=200, max_depth=8, random_state=42, n_jobs=8), threshold="median"),
    "MI200": lambda: SelectKBest(mutual_info_classif, k=200),
    "PCA30": lambda: PCA(n_components=30),
}
gkf = GroupKFold(n_splits=5)
res = []
for sn, sf in SELS.items():
    for cn, cf in CLFS.items():
        preds = np.zeros(len(y))
        for tr, te in gkf.split(Xf, y, grp):
            Xtr, Xte = Xf[tr], Xf[te]
            selm = sf() if sf is not None else None
            if selm is not None:
                selm.fit(Xtr, y[tr]); Xtr = selm.transform(Xtr); Xte = selm.transform(Xte)
            m = cf()
            m.fit(Xtr, y[tr])
            pr = m.predict_proba(Xte)[:, 1] if hasattr(m, "predict_proba") else m.decision_function(Xte)
            preds[te] = pr
        auc = roc_auc_score(y, preds)
        res.append({"fs": sn, "clf": cn, "auc": round(auc, 4)})
        print(sn, cn, round(auc, 4), flush=True)
## Dummy 基线
preds = np.zeros(len(y))
for tr, te in gkf.split(Xf, y, grp):
    preds[te] = y[tr].mean()
res.append({"fs": "None", "clf": "Dummy", "auc": round(roc_auc_score(y, preds), 4)})
rr = pd.DataFrame(res)
rr.to_csv(os.path.join(OUT, "modelzoo_cell_auc.csv"), index=False)
print(rr.sort_values("auc", ascending=False).head(10).to_string(index=False), flush=True)

## B) bootstrap L1 最小签名（用 LR_l1 系数频率；30 次 bootstrap）
from sklearn.linear_model import LogisticRegression as LR
coef_sum = np.zeros(Xf.shape[1]); freq = np.zeros(Xf.shape[1])
for b in range(30):
    boot = rng.choice(len(y), size=len(y), replace=True)
    Xb, yb = Xf[boot], y[boot]
    m = LR(penalty="l1", solver="liblinear", C=0.3, class_weight="balanced", max_iter=3000)
    m.fit(Xb, yb)
    coef_sum += (m.coef_[0] != 0).astype(int)
freq = coef_sum / 30
sig_genes = [genes[top[i]] for i in np.where(freq > 0.5)[0]]
print("bootstrap signature genes (freq>0.5):", len(sig_genes), sig_genes[:20], flush=True)
pd.DataFrame({"gene": [genes[top[i]] for i in range(len(freq))], "freq": freq}).to_csv(os.path.join(OUT, "sig_freq.csv"), index=False)
## 该签名上的 LR 5 折 AUC
if len(sig_genes) > 0:
    sel_sig = [i for i, g in enumerate(genes) if g in set(sig_genes)]
    Xs = X[:, [top[i] for i, g in enumerate(genes) if False]]
    # top 数组索引→基因
    pos = [i for i in np.where(freq > 0.5)[0]]
    Xs = Xf[:, pos]
    preds = np.zeros(len(y))
    for tr, te in gkf.split(Xs, y, grp):
        m = LR(penalty="l2", C=0.5, class_weight="balanced", max_iter=3000)
        m.fit(Xs[tr], y[tr]); preds[te] = m.predict_proba(Xs[te])[:, 1]
    print("signature LR-boot AUC:", round(roc_auc_score(y, preds), 4), flush=True)
print("DONE 45a")
