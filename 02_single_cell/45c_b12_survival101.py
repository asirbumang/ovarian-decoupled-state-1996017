# -*- coding: utf-8 -*-
"""45c v2: TCGA-OV 生存模型动物园 101-范式（C-index；lifelines 生存专属算法保留）。"""
import numpy as np, pandas as pd, sys, os, warnings
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectFromModel, SelectKBest, mutual_info_classif
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.model_selection import StratifiedKFold
from lifelines import CoxPHFitter, WeibullAFTFitter, LogNormalAFTFitter, LogLogisticAFTFitter
from lifelines.utils import concordance_index
sys.stdout.reconfigure(encoding='utf-8')
OUT = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b12"
os.makedirs(OUT, exist_ok=True)

mx = pd.read_csv(r"D:\子宫内膜异位症2\data_omics\tcga_ov_counts_matrix.csv.gz", index_col=0)
sm = pd.read_csv(r"D:\子宫内膜异位症2\data_omics\tcga_ov_sample_map.csv")
clin = pd.read_csv(r"D:\子宫内膜异位症2\data_omics\tcga_ov_clinical.csv")
clin["sample"] = clin["barcode"].str.replace(r"-\d{2}[A-Z]$", "", regex=True)
print("mapped:", clin["sample"].notna().sum(), flush=True)
clin = clin.dropna(subset=["sample"])
clin["os_days"] = np.where(clin["vital_status"] == "Dead", clin["days_to_death"], clin["days_to_last_followup"])
clin["event"] = (clin["vital_status"] == "Dead").astype(int)
clin = clin[clin["os_days"].notna() & (clin["os_days"] > 0)]
samples = [c for c in clin["sample"] if c in mx.columns]
X = mx[samples].T.values.astype(np.float32)  # patients × genes
genes = mx.index.tolist()
print("patients:", len(samples), " genes:", len(genes), flush=True)
cpm = X / X.sum(axis=1, keepdims=True) * 1e6
Xl = np.log2(cpm + 1)
y_t = clin.set_index("sample").loc[samples, "os_days"].values
y_e = clin.set_index("sample").loc[samples, "event"].values
# 特征：top2000 基因（全部/按方差；含状态签名）
sig = set(open(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\state_sig_up.txt").read().split()) | set(open(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\state_sig_dn.txt").read().split())
gi = {g: i for i, g in enumerate(genes)}
sig_pos = sorted([gi[g] for g in sig if g in gi])
var = Xl.var(axis=0)
top_var = np.argsort(-var)[:2000]
feat = np.unique(np.concatenate([top_var, sig_pos]))[:2500]
print("features:", len(feat), " deaths:", int(y_e.sum()), flush=True)
Xf = Xl[:, feat]

FS = {
    "None": lambda: None,
    "L1med": lambda: SelectFromModel(LogisticRegression(penalty="l1", solver="liblinear", C=0.4, class_weight="balanced", max_iter=3000), threshold="median"),
    "RFimp": lambda: SelectFromModel(ExtraTreesClassifier(n_estimators=200, max_depth=8, random_state=42, n_jobs=8), threshold="median"),
    "MI200": lambda: SelectKBest(mutual_info_classif, k=200),
    "PCA30": lambda: PCA(n_components=30),
}
MODELS = {
    "Cox": lambda: CoxPHFitter(penalizer=0.0),
    "Cox_l2": lambda: CoxPHFitter(penalizer=1.0),
    "Weibull": lambda: WeibullAFTFitter(penalizer=0.0),
    "LogNormal": lambda: LogNormalAFTFitter(penalizer=0.0),
    "LogLogistic": lambda: LogLogisticAFTFitter(penalizer=0.0),
}
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
res = []
for fn, ff in FS.items():
    for mn, mf in MODELS.items():
        cidx = []
        for tr, te in skf.split(Xf, y_e):
            selm = ff() if ff is not None else None
            Xtr, Xte = Xf[tr], Xf[te]
            if selm is not None:
                selm.fit(Xtr, y_e[tr])
                Xtr = selm.transform(Xtr).astype(np.float32); Xte = selm.transform(Xte).astype(np.float32)
            d = pd.DataFrame(Xtr, columns=[str(i) for i in range(Xtr.shape[1])])
            d["T"] = y_t[tr]; d["E"] = y_e[tr]
            try:
                m = mf(); m.fit(d, duration_col="T", event_col="E")
                dte = pd.DataFrame(Xte, columns=d.columns[:-2]).copy()
                dte["T"] = y_t[te]; dte["E"] = y_e[te]
                if hasattr(m, "predict_partial_hazard"):
                    risq = m.predict_partial_hazard(dte).values
                    c = concordance_index(y_t[te], -risq, y_e[te])
                else:
                    risq = m.predict_median(dte).values
                    c = concordance_index(y_t[te], -risq, y_e[te])
            except Exception as e2:
                c = np.nan
            cidx.append(c)
        c0 = float(np.nanmean(cidx))
        res.append({"fs": fn, "model": mn, "cindex": round(c0, 4)})
        print(fn, mn, round(c0, 4), flush=True)
rr = pd.DataFrame(res)
rr.to_csv(os.path.join(OUT, "survival101_cindex.csv"), index=False)
# Dummy
c = []
for tr, te in skf.split(Xf, y_e):
    c.append(concordance_index(y_t[te], np.random.default_rng(1).random(len(te)), y_e[te]))
rr = pd.concat([rr, pd.DataFrame([{"fs": "None", "model": "Dummy", "cindex": round(float(np.mean(c)), 4)}])])
rr.to_csv(os.path.join(OUT, "survival101_cindex.csv"), index=False)
print(rr.sort_values("cindex", ascending=False).head(12).to_string(index=False), flush=True)
print("DONE 45c")
