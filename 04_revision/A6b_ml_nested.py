# -*- coding: utf-8 -*-
"""
A6b: nested 5-fold CV for the patient-level survival model (top-500 model).

Outer 5-fold cross-validation; the top-500 variance filter and feature
standardisation are fitted on the training fold only; the ridge penalty is
fixed at alpha = 1.0. This is the model reported in Methods Section 8 and
Results Section 7 (mean C-index 0.46, SD 0.02).

A second nested run with the top-2,000 variance filter and an unstandardised
ridge penalty of 1.0 is reported by A6_ml.py (mean C-index 0.584, SD 0.038).
"""
import json, os, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.metrics import concordance_index_censored

ROOT = r"D:\子宫内膜异位症2"
OUT = os.path.join(ROOT, "R2_revision", "04_analysis", "results")

mx = pd.read_csv(os.path.join(ROOT, "data_omics", "tcga_ov_counts_matrix.csv.gz"), index_col=0)
clin = pd.read_csv(os.path.join(ROOT, "data_omics", "tcga_ov_clinical.csv"))
clin["sample"] = clin["barcode"].str.replace(r"-\d{2}[A-Z]$", "", regex=True)
clin = clin.dropna(subset=["sample"])
clin["os_days"] = np.where(clin.vital_status == "Dead", clin.days_to_death, clin.days_to_last_followup)
clin["event"] = (clin.vital_status == "Dead").astype(int)
clin = clin[clin.os_days.notna() & (clin.os_days > 0)]
samples = [c for c in clin["sample"] if c in mx.columns]
X = mx[samples].T.values.astype(np.float64)
X = np.log2(X / X.sum(1, keepdims=True) * 1e6 + 1)
t = clin.set_index("sample").loc[samples, "os_days"].values.astype(float)
e = clin.set_index("sample").loc[samples, "event"].values.astype(int)
y = np.array([(bool(a), float(b)) for a, b in zip(e, t)],
             dtype=[("event", bool), ("time", float)])

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
cs = []
for tr, te in skf.split(X, e):
    v = X[tr].var(0)
    f = np.argsort(-v)[:500]
    mu = X[tr][:, f].mean(0)
    sd = X[tr][:, f].std(0) + 1e-8
    m = CoxPHSurvivalAnalysis(alpha=1.0).fit((X[tr][:, f] - mu) / sd, y[tr])
    p = m.predict((X[te][:, f] - mu) / sd)
    cs.append(concordance_index_censored(e[te].astype(bool), t[te], -p)[0])

res = {"n": int(len(t)), "events": int(e.sum()),
       "nested_cindex_mean": float(np.mean(cs)),
       "sd": float(np.std(cs, ddof=1)),
       "folds": [float(x) for x in cs],
       "model": "ridge Cox alpha=1, top-500 variance genes fitted inside each training fold"}
json.dump(res, open(os.path.join(OUT, "A6b_ml_nested.json"), "w"), indent=2)
print(json.dumps(res, indent=2))
