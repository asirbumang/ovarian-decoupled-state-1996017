# -*- coding: utf-8 -*-
"""
A1b: Is the three-way (M1-high, M2-low, M3-low) pattern present beyond the
three pairwise associations, and does it survive adjustment for cell-type
composition and for patient identity?

Method: Poisson log-linear models on the 2x2x2 table (A = M1 >= P75,
B = M2 <= P25, C = M3 <= P25). The ABC coefficient is the three-way
interaction; expected count in the ABC cell comes from the model without ABC.
Models:
  M0  A + B + C + AB + AC + BC                     (pairwise-only null)
  M1  M0 + ABC
  M2  M0 + cell_type + cell_type:A/B/C + ABC       (composition-adjusted)
  M3  M0 + patient  + patient:A/B/C  + ABC         (patient-adjusted)
  M4  M2 + patient + patient:A/B/C                 (both)
"""
import json, os
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

ROOT = r"D:\子宫内膜异位症2"
OUT = os.path.join(ROOT, "R2_revision", "04_analysis", "results")
df = pd.read_csv(os.path.join(ROOT, "analysis", "L1_scRNA", "results",
                              "module_scores", "scores_GSE266577.csv"))
df["patient"] = df["publication_patient_code_final"].astype(str)
A = (df.M1z >= np.quantile(df.M1z, .75)).astype(int)
B = (df.M2z <= np.quantile(df.M2z, .25)).astype(int)
C = (df.M3z <= np.quantile(df.M3z, .25)).astype(int)
df["A"], df["B"], df["C"] = A, B, C
df["ABC"] = A * B * C


def fit_long(agg, extra_terms=(), weighted_cols=()):
    """agg: columns A,B,C + stratifier columns + y."""
    d = agg.copy()
    X = pd.DataFrame({"const": 1.0, "A": d.A.astype(float), "B": d.B.astype(float),
                      "C": d.C.astype(float)})
    X["AB"] = X.A * X.B; X["AC"] = X.A * X.C; X["BC"] = X.B * X.C
    terms = {}
    for col in extra_terms:
        levels = sorted(d[col].unique())
        for lev in levels[1:]:
            ind = (d[col] == lev).astype(float)
            terms[f"{col}={lev}"] = ind
            terms[f"{col}={lev}:A"] = ind * X.A
            terms[f"{col}={lev}:B"] = ind * X.B
            terms[f"{col}={lev}:C"] = ind * X.C
    for k, v in terms.items():
        X[k] = v
    X["ABC"] = X.A * X.B * X.C
    m = sm.GLM(d.y.astype(float), X, family=sm.families.Poisson()).fit()
    return m, X


def summarize(m, d, label):
    b = m.params.get("ABC", np.nan); se = m.bse.get("ABC", np.nan)
    # expected ABC count from the same model without the ABC term
    X = m.model.exog
    keep = [i for i, n in enumerate(m.model.exog_names) if n != "ABC"]
    m0 = sm.GLM(m.model.endog, X[:, keep], family=sm.families.Poisson()).fit()
    mask = (d.A == 1) & (d.B == 1) & (d.C == 1)
    obs = d.loc[mask, "y"].sum()
    exp = np.asarray(m0.fittedvalues)[mask.values].sum()
    return {"model": label, "ABC_coef": float(b), "ABC_se": float(se),
            "z": float(b / se), "p": float(2 * stats.norm.sf(abs(b / se))),
            "observed_ABC": int(obs), "expected_ABC_pairwise": float(exp),
            "OE": float(obs / exp), "n_cells": int(d.y.sum()),
            "deviance": float(m.deviance), "df_resid": int(m.df_resid)}


def aggregate(by):
    g = df.groupby(by + ["A", "B", "C"]).size().reset_index(name="y")
    full = pd.MultiIndex.from_product([sorted(df[c].unique()) for c in by] + [[0, 1], [0, 1], [0, 1]])
    g = g.set_index(by + ["A", "B", "C"]).reindex(full, fill_value=0).reset_index()
    g.columns = by + ["A", "B", "C", "y"]
    return g


res = {"models": []}
# overall
g0 = df.groupby(["A", "B", "C"]).size().reset_index(name="y")
full0 = pd.MultiIndex.from_product([[0, 1], [0, 1], [0, 1]])
g0 = g0.set_index(["A", "B", "C"]).reindex(full0, fill_value=0).reset_index()
g0.columns = ["A", "B", "C", "y"]
m0, _ = fit_long(g0); m1, _ = fit_long(g0)
res["models"].append(summarize(m1, g0, "pairwise-only null + ABC (all cells)"))

# composition-adjusted (cell types with >= 300 cells to avoid separation)
keep_ct = df.cell_type.value_counts()
keep_ct = keep_ct[keep_ct >= 300].index
d_ct = df[df.cell_type.isin(keep_ct)].copy()
g_ct = aggregate(["cell_type"]) if False else None
tmp = d_ct.groupby(["cell_type", "A", "B", "C"]).size().reset_index(name="y")
full = pd.MultiIndex.from_product([sorted(keep_ct), [0, 1], [0, 1], [0, 1]])
tmp = tmp.set_index(["cell_type", "A", "B", "C"]).reindex(full, fill_value=0).reset_index()
tmp.columns = ["cell_type", "A", "B", "C", "y"]
m2, _ = fit_long(tmp, extra_terms=["cell_type"])
res["models"].append(summarize(m2, tmp, "composition-adjusted (cell_type)"))

# patient-adjusted
tmp2 = df.groupby(["patient", "A", "B", "C"]).size().reset_index(name="y")
full2 = pd.MultiIndex.from_product([sorted(df.patient.unique()), [0, 1], [0, 1], [0, 1]])
tmp2 = tmp2.set_index(["patient", "A", "B", "C"]).reindex(full2, fill_value=0).reset_index()
tmp2.columns = ["patient", "A", "B", "C", "y"]
m3, _ = fit_long(tmp2, extra_terms=["patient"])
res["models"].append(summarize(m3, tmp2, "patient-adjusted"))

# both
tmp3 = df.groupby(["cell_type", "patient", "A", "B", "C"]).size().reset_index(name="y")
res["note_both"] = "combined model fitted on the cell-type x patient x ABC table"
m4, _ = fit_long(tmp3, extra_terms=["cell_type", "patient"]) if len(tmp3) > 100 else (None, None)
if m4 is not None:
    res["models"].append(summarize(m4, tmp3.assign(**{"cell_type": "all"}), "cell-type + patient adjusted"))

# within-type three-way tables (descriptive, for cell types with >=1000 cells)
within = []
for ct, g in df.groupby("cell_type"):
    if len(g) < 1000:
        continue
    n111 = int(((g.A == 1) & (g.B == 1) & (g.C == 1)).sum())
    n_expected = len(g) * ((g.A == 1).mean() * (g.B == 1).mean() * (g.C == 1).mean())
    within.append({"cell_type": ct, "n": len(g), "state": n111,
                   "pct": 100 * n111 / len(g), "expected_independent": n_expected,
                   "OE": n111 / n_expected if n_expected else np.nan})
res["within_celltype_descriptive"] = within

with open(os.path.join(OUT, "A1b_loglinear.json"), "w", encoding="utf-8") as f:
    json.dump(res, f, indent=2)
print(json.dumps(res, indent=2))
