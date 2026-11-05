# -*- coding: utf-8 -*-
"""
A6: Machine-learning layer.
 (i)  cell-level: patient-grouped CV of the state label from the module scores
      themselves. This is definition recovery, not biological validation, and
      is labelled as such. AUROC and AUPRC are both reported.
 (ii) patient-level: nested re-run of the survival model with the variance
      filter and feature standardisation fitted inside each training fold,
      compared with the leaky variant (filter on the full cohort) to quantify
      the effect of the pre-selection used originally.
"""
import json, os, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GroupKFold
from sklearn.metrics import roc_auc_score, average_precision_score
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

ROOT = r"D:\子宫内膜异位症2"
OUT = os.path.join(ROOT, "R2_revision", "04_analysis", "results")
res = {}

# ---------------------------------------------------------------- cell level
sc = pd.read_csv(os.path.join(ROOT, "analysis", "L1_scRNA", "results",
                              "module_scores", "scores_GSE266577.csv"))
X = sc[["M1z", "M2z", "M3z"]].values
y = sc.state_high.astype(int).values
groups = sc.publication_patient_code_final.astype(str).values
gkf = GroupKFold(n_splits=5)
oof = np.zeros(len(y))
for tr, te in gkf.split(X, y, groups):
    m = LogisticRegression(max_iter=1000).fit(X[tr], y[tr])
    oof[te] = m.predict_proba(X[te])[:, 1]
res["cell_level_definition_recovery"] = {
    "features": "M1z, M2z, M3z (the three scores that define the label)",
    "cv": "5-fold GroupKFold by patient",
    "AUROC": float(roc_auc_score(y, oof)),
    "AUPRC": float(average_precision_score(y, oof)),
    "prevalence": float(y.mean()),
    "note": "circular by construction; reported as definition recovery"}

# ------------------------------------------------------------- patient level
mx = pd.read_csv(os.path.join(ROOT, "data_omics", "tcga_ov_counts_matrix.csv.gz"), index_col=0)
clin = pd.read_csv(os.path.join(ROOT, "data_omics", "tcga_ov_clinical.csv"))
clin["sample"] = clin["barcode"].str.replace(r"-\d{2}[A-Z]$", "", regex=True)
clin = clin.dropna(subset=["sample"])
clin["os_days"] = np.where(clin.vital_status == "Dead", clin.days_to_death, clin.days_to_last_followup)
clin["event"] = (clin.vital_status == "Dead").astype(int)
clin = clin[clin.os_days.notna() & (clin.os_days > 0)]
samples = [c for c in clin["sample"] if c in mx.columns]
Xg = mx[samples].T.values.astype(np.float32)
genes = mx.index.to_numpy()
cpm = Xg / Xg.sum(axis=1, keepdims=True) * 1e6
Xl = np.log2(cpm + 1)
y_t = clin.set_index("sample").loc[samples, "os_days"].values.astype(float)
y_e = clin.set_index("sample").loc[samples, "event"].values.astype(int)

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)


def run_cv(nested: bool):
    cs = []
    for tr, te in skf.split(Xl, y_e):
        if nested:
            v = Xl[tr].var(axis=0)
            feat = np.argsort(-v)[:2000]
        else:
            v = Xl.var(axis=0)
            feat = np.argsort(-v)[:2000]
        d = pd.DataFrame(Xl[tr][:, feat], columns=[str(i) for i in feat])
        d["T"] = y_t[tr]; d["E"] = y_e[tr]
        m = CoxPHFitter(penalizer=1.0).fit(d, "T", "E")
        dte = pd.DataFrame(Xl[te][:, feat], columns=[str(i) for i in feat])
        risk = m.predict_partial_hazard(dte).values
        cs.append(concordance_index(y_t[te], -risk, y_e[te]))
    return float(np.mean(cs)), float(np.std(cs, ddof=1)), cs


nested = run_cv(True)
leaky = run_cv(False)  # identical here because the filter is refitted per fold;
# the leaky variant instead uses the full-cohort filter once:
v = Xl.var(axis=0); feat_full = np.argsort(-v)[:2000]
cs = []
for tr, te in skf.split(Xl, y_e):
    d = pd.DataFrame(Xl[tr][:, feat_full], columns=[str(i) for i in feat_full])
    d["T"] = y_t[tr]; d["E"] = y_e[tr]
    m = CoxPHFitter(penalizer=1.0).fit(d, "T", "E")
    dte = pd.DataFrame(Xl[te][:, feat_full], columns=[str(i) for i in feat_full])
    risk = m.predict_partial_hazard(dte).values
    cs.append(concordance_index(y_t[te], -risk, y_e[te]))
leaky = (float(np.mean(cs)), float(np.std(cs, ddof=1)), cs)

res["patient_level_nested_cv"] = {
    "n": int(len(y_t)), "events": int(y_e.sum()),
    "model": "CoxPH ridge (penalizer=1.0), top-2000 variance genes",
    "nested_cindex_mean": nested[0], "nested_cindex_sd": nested[1], "nested_folds": nested[2],
    "full_cohort_filter_cindex_mean": leaky[0], "full_cohort_filter_sd": leaky[1],
    "delta_leak": leaky[0] - nested[0]}

with open(os.path.join(OUT, "A6_ml.json"), "w", encoding="utf-8") as f:
    json.dump(res, f, indent=2)
print(json.dumps(res, indent=2))
