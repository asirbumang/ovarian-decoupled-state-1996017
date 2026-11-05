# -*- coding: utf-8 -*-
"""
A2: sensitivity of the state definition to percentile cutoffs and to
    alternative operational definitions (carrier enrichment must be stable).
A3: technical covariates (depth, detected genes, mitochondrial fraction) and
    composition effects on module scores and on state membership.
"""
import json, os, itertools
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

ROOT = r"D:\子宫内膜异位症2"
OUT = os.path.join(ROOT, "R2_revision", "04_analysis", "results")
os.makedirs(OUT, exist_ok=True)
df = pd.read_csv(os.path.join(ROOT, "analysis", "L1_scRNA", "results",
                              "module_scores", "scores_GSE266577.csv"))
df["patient"] = df["publication_patient_code_final"].astype(str)
n = len(df)


def state_from(m1, m2, m3, p_hi=75, p_lo=25):
    return ((m1 >= np.quantile(m1, p_hi / 100)) &
            (m2 <= np.quantile(m2, p_lo / 100)) &
            (m3 <= np.quantile(m3, p_lo / 100)))


def carrier_stats(mask):
    out = {}
    for name, sel in [("plasma", df.cell_type == "Plasma cells"),
                      ("fibroblast", df.cell_type == "Fibroblasts"),
                      ("epithelial", df.cell_type == "Epithelial cells")]:
        a = int((sel & mask).sum()); b = int((sel & ~mask).sum())
        c = int((~sel & mask).sum()); d = int((~sel & ~mask).sum())
        if min(a, b, c, d) == 0:
            a2, b2, c2, d2 = a + .5, b + .5, c + .5, d + .5
        else:
            a2, b2, c2, d2 = a, b, c, d
        out[name] = {"OR": (a2 * d2) / (b2 * c2), "state": a}
    return out


base = df.state_high.values
sens = []
for p_hi, p_lo in itertools.product([70, 72.5, 75, 77.5, 80], [20, 22.5, 25, 27.5, 30]):
    m = state_from(df.M1z.values, df.M2z.values, df.M3z.values, p_hi, p_lo)
    cs = carrier_stats(m)
    sens.append({"M1_pct": p_hi, "M2M3_pct": p_lo, "pct_cells": 100 * m.mean(),
                 "n": int(m.sum()), "jaccard_vs_base": float((m & base).sum() / (m | base).sum()),
                 "plasma_OR": cs["plasma"]["OR"], "fibroblast_OR": cs["fibroblast"]["OR"],
                 "epithelial_OR": cs["epithelial"]["OR"]})
sens = pd.DataFrame(sens)
sens.to_csv(os.path.join(OUT, "A2_threshold_grid.csv"), index=False)

# alternative definitions
alts = {}
d = df.decouple.values
alts["decouple_top10pct"] = d >= np.quantile(d, .90)
alts["decouple_top5pct"] = d >= np.quantile(d, .95)
alts["median_split"] = (df.M1z.values > 0) & (df.M2z.values < 0) & (df.M3z.values < 0)
alts["decouple_gt_1"] = d > 1
# residual definition: M1 unexplained by M2 and M3 (orthogonalised immune arm)
X = sm.add_constant(np.column_stack([df.M2z.values, df.M3z.values]))
resid = sm.OLS(df.M1z.values, X).fit().resid
alts["M1_residual_top25"] = resid >= np.quantile(resid, .75)
alt_rows = []
for k, m in alts.items():
    cs = carrier_stats(np.asarray(m))
    alt_rows.append({"definition": k, "pct_cells": 100 * np.mean(m), "n": int(np.sum(m)),
                     "jaccard_vs_base": float((np.asarray(m) & base).sum() / (np.asarray(m) | base).sum()),
                     "plasma_OR": cs["plasma"]["OR"], "fibroblast_OR": cs["fibroblast"]["OR"],
                     "epithelial_OR": cs["epithelial"]["OR"]})
pd.DataFrame(alt_rows).to_csv(os.path.join(OUT, "A2_alternative_definitions.csv"), index=False)

# ------------------------------------------------------------------ A3 technical
tech = {}
tech["state_median"] = {k: {"state": float(df.loc[df.state_high, k].median()),
                            "rest": float(df.loc[~df.state_high, k].median())}
                        for k in ["nCount_RNA", "nFeature_RNA", "percent.mt"]}
for k in ["nCount_RNA", "nFeature_RNA", "percent.mt"]:
    u = stats.mannwhitneyu(df.loc[df.state_high, k], df.loc[~df.state_high, k])
    auc = u.statistic / (df.state_high.sum() * (~df.state_high).sum())
    tech.setdefault("state_mannwhitney", {})[k] = {"p": float(u.pvalue), "AUC": float(auc)}

# module z vs technical covariates
tech["module_spearman"] = {}
for k in ["nCount_RNA", "nFeature_RNA", "percent.mt"]:
    tech["module_spearman"][k] = {m: {"rho": float(stats.spearmanr(df[m], df[k]).statistic),
                                      "p": float(stats.spearmanr(df[m], df[k]).pvalue)}
                                  for m in ["M1z", "M2z", "M3z", "decouple"]}

# state ~ technical covariates, adjusted for cell type (logistic)
keep_ct = [c for c in df.cell_type.unique()
           if (df.cell_type == c).sum() >= 100 and df.loc[df.cell_type == c, "state_high"].sum() > 0]
sub = df[df.cell_type.isin(keep_ct)].copy()
ct = pd.get_dummies(sub.cell_type, drop_first=True).astype(float)
z_feat = (sub["nFeature_RNA"] - sub["nFeature_RNA"].mean()) / sub["nFeature_RNA"].std()
z_mt = (sub["percent.mt"] - sub["percent.mt"].mean()) / sub["percent.mt"].std()
X = pd.concat([pd.DataFrame({"z_nFeature": z_feat, "z_percent_mt": z_mt}), ct], axis=1)
X = sm.add_constant(X)
X.columns = [str(c) for c in X.columns]
logit = sm.Logit(sub.state_high.astype(int), X).fit(disp=0, method="bfgs", maxiter=500)
tech["logistic_adjusted"] = {
    "n_cells": int(len(sub)), "n_celltypes": len(keep_ct),
    "nFeature_per_SD_OR": float(np.exp(logit.params["z_nFeature"])),
    "nFeature_p": float(logit.pvalues["z_nFeature"]),
    "percent_mt_per_SD_OR": float(np.exp(logit.params["z_percent_mt"])),
    "percent_mt_p": float(logit.pvalues["z_percent_mt"])}

# state rate across nFeature quintiles
df["nFeature_q"] = pd.qcut(df.nFeature_RNA, 5, labels=False)
tech["state_rate_by_nFeature_quintile"] = df.groupby("nFeature_q").state_high.agg(["mean", "size"]).reset_index().to_dict("records")

# residualised modules: remove technical + composition effects, recompute state
ct_full = pd.get_dummies(df.cell_type, drop_first=True).astype(float)
Z = pd.concat([df[["nCount_RNA", "nFeature_RNA", "percent.mt"]].astype(float), ct_full], axis=1)
Z = sm.add_constant(Z); Z.columns = [str(c) for c in Z.columns]
res_mod = {}
for m in ["M1z", "M2z", "M3z"]:
    res_mod[m] = sm.OLS(df[m].values, Z).fit().resid
m_res = state_from(res_mod["M1z"], res_mod["M2z"], res_mod["M3z"])
cs = carrier_stats(m_res)
tech["residualised_definition"] = {
    "pct_cells": float(100 * m_res.mean()), "n": int(m_res.sum()),
    "jaccard_vs_base": float((m_res & base).sum() / (m_res | base).sum()),
    "plasma_OR": cs["plasma"]["OR"], "fibroblast_OR": cs["fibroblast"]["OR"],
    "epithelial_OR": cs["epithelial"]["OR"],
    "corr_resid_vs_orig": {m: float(np.corrcoef(res_mod[m], df[m])[0, 1]) for m in res_mod}}

with open(os.path.join(OUT, "A3_technical.json"), "w", encoding="utf-8") as f:
    json.dump(tech, f, indent=2)

print(sens.head(10).to_string())
print(pd.DataFrame(alt_rows).to_string())
print(json.dumps(tech, indent=2)[:3000])
