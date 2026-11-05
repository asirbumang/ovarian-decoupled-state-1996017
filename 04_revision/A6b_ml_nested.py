# -*- coding: utf-8 -*-
"""
A6b: fast nested CV for the patient-level survival model using scikit-survival.
Outer 5-fold CV; variance filter and scaling fitted on the training fold only;
elastic-net Cox (Coxnet) with the alpha path selected by inner 3-fold CV.
Also reports the full-cohort variance-filter variant for comparison.
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
cpm = X / X.sum(axis=1, keepdims=True) * 1e6
Xl = np.log2(cpm + 1)
t = clin.set_index("sample").loc[samples, "os_days"].values.astype(float)
e = clin.set_index("sample").loc[samples, "event"].values.astype(bool)
y = np.array([(bool(a), float(b)) for a, b in zip(e, t)],
             dtype=[("event", bool), ("time", float)])

outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)


def fit_predict(Xtr, ytr, Xte):
    # standardise on training fold
    mu = Xtr.mean(0); sd = Xtr.std(0) + 1e-8
    Xtr = (Xtr - mu) / sd; Xte = (Xte - mu) / sd
    # inner CV over a small ridge-alpha grid
    inner = StratifiedKFold(n_splits=3, shuffle=True, random_state=7)
    grid = [0.1, 1.0, 10.0]
    scores = np.full((len(grid), 3), np.nan)
    for k, (itr, ite) in enumerate(inner.split(Xtr, ytr["event"])):
        for ai, a in enumerate(grid):
            try:
                m = CoxPHSurvivalAnalysis(alpha=a).fit(Xtr[itr], ytr[itr])
                pred = m.predict(Xtr[ite])
                scores[ai, k] = concordance_index_censored(ytr[ite]["event"], ytr[ite]["time"], -pred)[0]
            except Exception:
                pass
    best = grid[int(np.nanargmax(np.nanmean(scores, 1)))]
    m = CoxPHSurvivalAnalysis(alpha=best).fit(Xtr, ytr)
    return m.predict(Xte)


def run(nested):
    cs = []
    for tr, te in outer.split(Xl, e):
        v = Xl[tr].var(0) if nested else Xl.var(0)
        feat = np.argsort(-v)[:1000]
        pred = fit_predict(Xl[tr][:, feat], y[tr], Xl[te][:, feat])
        cs.append(concordance_index_censored(e[te], t[te], -pred)[0])
    return float(np.mean(cs)), float(np.std(cs, ddof=1)), [float(c) for c in cs]


nested = run(True)
res = {"n": int(len(t)), "events": int(e.sum()),
       "model": "ridge Cox; top-1000 variance genes; alpha by inner 3-fold CV",
       "nested_cindex_mean": nested[0], "nested_cindex_sd": nested[1], "nested_folds": nested[2]}
with open(os.path.join(OUT, "A6b_ml_nested.json"), "w", encoding="utf-8") as f:
    json.dump(res, f, indent=2)
print(json.dumps(res, indent=2))
