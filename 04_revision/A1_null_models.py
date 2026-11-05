# -*- coding: utf-8 -*-
"""
A1: Null models, composition adjustment and patient-level accounting for the
immune-metabolic decoupled state (GSE266577 primary reference).

Inputs : analysis/L1_scRNA/results/module_scores/scores_GSE266577.csv
Outputs: R2_revision/04_analysis/results/A1_*.csv/json/md

All tests are two-sided unless stated. Random seed = 20261105.
"""
import json, os, math
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.contingency_tables import StratifiedTable

SEED = 20261105
rng = np.random.default_rng(SEED)
ROOT = r"D:\子宫内膜异位症2"
OUT = os.path.join(ROOT, "R2_revision", "04_analysis", "results")
os.makedirs(OUT, exist_ok=True)

df = pd.read_csv(os.path.join(ROOT, "analysis", "L1_scRNA", "results",
                              "module_scores", "scores_GSE266577.csv"))
df["patient"] = df["publication_patient_code_final"].astype(str)
df["sample"] = df["publication_sample_code_final"].astype(str)

# state_high in the source table is the P75/P25/P25 rule on raw module scores
state = df["state_high"].astype(bool).values
n = len(df)
obs = int(state.sum())

# ---------------------------------------------------------------- marginals
q1 = np.quantile(df.M1z, 0.75)
q2 = np.quantile(df.M2z, 0.25)
q3 = np.quantile(df.M3z, 0.25)
p1 = float((df.M1z >= q1).mean())
p2 = float((df.M2z <= q2).mean())
p3 = float((df.M3z <= q3).mean())
p_ind = p1 * p2 * p3

res = {
    "n_cells": n, "n_state": obs, "state_pct": 100 * obs / n,
    "n_patients": int(df["patient"].nunique()), "n_samples": int(df["sample"].nunique()),
    "marginals": {"M1_ge_P75": p1, "M2_le_P25": p2, "M3_le_P25": p3},
    "independence_expected_pct": 100 * p_ind,
    "observed_over_independence": (obs / n) / p_ind,
    "module_correlations": df[["M1z", "M2z", "M3z"]].corr().round(4).to_dict(),
}

# exact binomial test of observed count vs independence expectation
bt = stats.binomtest(obs, n, p_ind)
res["binomial_vs_independence"] = {"p_two_sided": bt.pvalue,
                                   "ci95": list(bt.proportion_ci(0.95))}

# ---------------------------------------------------------- permutation null
# Shuffle M2 and M3 values independently across cells: preserves each module's
# marginal distribution and cell-type composition of M1, destroys cell-level
# co-occurrence.
B = 10000
null_counts = np.empty(B, dtype=np.int64)
m1 = df.M1z.values; m2 = df.M2z.values; m3 = df.M3z.values
hi1 = m1 >= q1; lo2 = m2 <= q2; lo3 = m3 <= q3
for b in range(B):
    null_counts[b] = int((hi1 & rng.permutation(lo2) & rng.permutation(lo3)).sum())
p_perm_low = (1 + int((null_counts <= obs).sum())) / (B + 1)
p_perm_high = (1 + int((null_counts >= obs).sum())) / (B + 1)
res["permutation_null"] = {
    "B": B, "mean_count": float(null_counts.mean()),
    "sd_count": float(null_counts.std(ddof=1)),
    "p_depletion": p_perm_low, "p_enrichment": p_perm_high,
    "null_pct_mean": 100 * float(null_counts.mean()) / n}

# ------------------------------------------- copula null (marginals + corr)
# Gaussian copula with the observed correlation matrix and the empirical
# marginal quantile functions: preserves both marginal distributions and
# pairwise correlations, i.e. tests whether the region frequency follows from
# a continuous correlated distribution rather than from a distinct state.
R = df[["M1z", "M2z", "M3z"]].corr().values
L = np.linalg.cholesky(R)
M = 400000
Z = rng.standard_normal((M, 3)) @ L.T
U = stats.norm.cdf(Z)
sim = np.column_stack([
    np.quantile(m1, U[:, 0]), np.quantile(m2, U[:, 1]), np.quantile(m3, U[:, 2])])
p_cop = float(((sim[:, 0] >= q1) & (sim[:, 1] <= q2) & (sim[:, 2] <= q3)).mean())
res["copula_null"] = {"M": M, "expected_pct": 100 * p_cop,
                      "observed_over_expected": (obs / n) / p_cop}

# ---------------------------------------- composition: observed vs expected
tab = []
for ct, g in df.groupby("cell_type"):
    n_t = len(g); k_t = int(g.state_high.sum())
    e_glob = n_t * obs / n                       # global-rate expectation
    # within-type independence expectation (type-specific marginals)
    pt = ((g.M1z >= np.quantile(g.M1z, .75)).mean() *
          (g.M2z <= np.quantile(g.M2z, .25)).mean() *
          (g.M3z <= np.quantile(g.M3z, .25)).mean())
    e_type = n_t * pt
    p_binom = stats.binomtest(k_t, n_t, obs / n).pvalue
    tab.append({"cell_type": ct, "n": n_t, "n_state": k_t, "pct": 100 * k_t / n_t,
                "expected_global": e_glob, "OE_global": k_t / e_glob if e_glob else np.nan,
                "expected_within_type": e_type,
                "OE_within_type": k_t / e_type if e_type else np.nan,
                "p_vs_global_binomial": p_binom})
comp = pd.DataFrame(tab).sort_values("pct", ascending=False)
comp.to_csv(os.path.join(OUT, "A1_composition_by_celltype.csv"), index=False)

# carrier enrichment (Fisher) + CMH across patients
def fisher_or(mask):
    a = int((mask & state).sum()); b = int((mask & ~state).sum())
    c = int((~mask & state).sum()); d = int((~mask & ~state).sum())
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    lo, hi = _or_ci(a, b, c, d)
    return {"n_carrier": int(mask.sum()), "state_in": a, "pct_in": 100 * a / mask.sum(),
            "state_out": c, "pct_out": 100 * c / (~mask).sum(),
            "OR": orr, "CI_low": lo, "CI_high": hi, "p_fisher": p}

def _or_ci(a, b, c, d):
    if min(a, b, c, d) == 0:
        a, b, c, d = a + .5, b + .5, c + .5, d + .5
    lor = math.log(a * d / (b * c)); se = math.sqrt(1/a + 1/b + 1/c + 1/d)
    return math.exp(lor - 1.96 * se), math.exp(lor + 1.96 * se)

res["carrier_enrichment"] = {
    "plasma": fisher_or((df.cell_type == "Plasma cells").values),
    "fibroblast": fisher_or((df.cell_type == "Fibroblasts").values),
    "epithelial": fisher_or((df.cell_type == "Epithelial cells").values)}

# CMH across patients (stratified by patient)
cmh = {}
for name, mask in [("plasma", df.cell_type == "Plasma cells"),
                   ("fibroblast", df.cell_type == "Fibroblasts")]:
    tabs = []
    for pt, g in df.groupby("patient"):
        m = mask.loc[g.index].values; s = g.state_high.values
        tabs.append([[int((m & s).sum()), int((m & ~s).sum())],
                     [int((~m & s).sum()), int((~m & ~s).sum())]])
    st = StratifiedTable(tabs)
    cmh[name] = {"OR_MH": float(st.oddsratio_pooled),
                 "CI_low": float(st.oddsratio_pooled_confint()[0]),
                 "CI_high": float(st.oddsratio_pooled_confint()[1]),
                 "p": float(st.test_null_odds().pvalue),
                 "n_patients": len(tabs)}
res["carrier_CMH_patient_stratified"] = cmh

# ------------------------------------------------------- patient-level rates
pat = df.groupby("patient").agg(n=("state_high", "size"), k=("state_high", "sum"),
                                plasma_frac=("cell_type", lambda s: (s == "Plasma cells").mean()),
                                fibro_frac=("cell_type", lambda s: (s == "Fibroblasts").mean()))
pat["rate"] = pat.k / pat.n
pat = pat.sort_values("rate", ascending=False)
pat.to_csv(os.path.join(OUT, "A1_patient_state_rates.csv"))
pbar = obs / n
chi2 = float((((pat.k - pat.n * pbar) ** 2) / (pat.n * pbar * (1 - pbar))).sum())
dfree = len(pat) - 1
res["patient_overdispersion"] = {
    "chi2": chi2, "df": dfree, "p": float(stats.chi2.sf(chi2, dfree)),
    "patients_with_state": int((pat.k > 0).sum()), "n_patients": len(pat),
    "rate_min": float(pat.rate.min()), "rate_max": float(pat.rate.max()),
    "rate_median": float(pat.rate.median()),
    "spearman_rate_vs_plasma": stats.spearmanr(pat.rate, pat.plasma_frac).statistic,
    "spearman_rate_p_plasma": stats.spearmanr(pat.rate, pat.plasma_frac).pvalue,
    "spearman_rate_vs_fibro": stats.spearmanr(pat.rate, pat.fibro_frac).statistic,
    "spearman_rate_p_fibro": stats.spearmanr(pat.rate, pat.fibro_frac).pvalue}

# variance components (one-way random effects on the logit-ish scale via ANOVA
# of patient rates; ICC approximated as between-patient variance / total)
rates = pat.rate.values
grand = rates.mean()
var_between = max(0.0, rates.var(ddof=1) - grand * (1 - grand) / pat.n.mean())
var_total = var_between + grand * (1 - grand) / pat.n.mean()
res["patient_ICC_approx"] = float(var_between / var_total)

with open(os.path.join(OUT, "A1_null_models.json"), "w", encoding="utf-8") as f:
    json.dump(res, f, indent=2, default=float)

print(json.dumps(res, indent=2, default=float))
