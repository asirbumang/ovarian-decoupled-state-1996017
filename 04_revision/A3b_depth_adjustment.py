# -*- coding: utf-8 -*-
"""
A3b: is carrier enrichment explained by library complexity?
 (i)  median nFeature by lineage;
 (ii) state OR for plasma/fibroblast within nFeature deciles (Mantel-Haenszel);
 (iii) state re-derived after residualising modules on technical covariates ONLY
       (no cell type), so that between-lineage signal is preserved.
"""
import json, os
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.contingency_tables import StratifiedTable

ROOT = r"D:\子宫内膜异位症2"
OUT = os.path.join(ROOT, "R2_revision", "04_analysis", "results")
df = pd.read_csv(os.path.join(ROOT, "analysis", "L1_scRNA", "results",
                              "module_scores", "scores_GSE266577.csv"))
res = {}

res["median_nFeature_by_lineage"] = (df.groupby("cell_type").nFeature_RNA.median()
                                     .sort_values(ascending=False).to_dict())
df["nf_decile"] = pd.qcut(df.nFeature_RNA, 10, labels=False)

def or_ci(a, b, c, d):
    a2, b2, c2, d2 = (a + .5, b + .5, c + .5, d + .5) if min(a, b, c, d) == 0 else (a, b, c, d)
    lor = np.log(a2 * d2 / (b2 * c2)); se = np.sqrt(1/a2 + 1/b2 + 1/c2 + 1/d2)
    return float(np.exp(lor)), float(np.exp(lor - 1.96*se)), float(np.exp(lor + 1.96*se))

# (ii) CMH within nFeature deciles
cmh = {}
for name, mask in [("plasma", df.cell_type == "Plasma cells"),
                   ("fibroblast", df.cell_type == "Fibroblasts")]:
    tabs = []
    for dec, g in df.groupby("nf_decile"):
        m = mask.loc[g.index].values; s = g.state_high.values
        tabs.append([[int((m & s).sum()), int((m & ~s).sum())],
                     [int((~m & s).sum()), int((~m & ~s).sum())]])
    st = StratifiedTable(tabs)
    cmh[name] = {"OR_MH": float(st.oddsratio_pooled),
                 "CI": [float(st.oddsratio_pooled_confint()[0]), float(st.oddsratio_pooled_confint()[1])],
                 "p": float(st.test_null_odds().pvalue)}
res["carrier_CMH_within_nFeature_deciles"] = cmh

# stratified state rates by nFeature decile for carriers vs rest
rows = []
for dec, g in df.groupby("nf_decile"):
    rows.append({"decile": int(dec), "n": len(g),
                 "plasma_state": int(((g.cell_type == "Plasma cells") & g.state_high).sum()),
                 "plasma_n": int((g.cell_type == "Plasma cells").sum()),
                 "fibro_state": int(((g.cell_type == "Fibroblasts") & g.state_high).sum()),
                 "fibro_n": int((g.cell_type == "Fibroblasts").sum()),
                 "all_state": int(g.state_high.sum())})
res["nFeature_decile_table"] = rows

# (iii) technical-only residualisation
Z = df[["nCount_RNA", "nFeature_RNA", "percent.mt"]].astype(float)
Z = sm.add_constant((Z - Z.mean()) / Z.std())
resid = {m: sm.OLS(df[m].values, Z).fit().resid for m in ["M1z", "M2z", "M3z"]}
mask = ((resid["M1z"] >= np.quantile(resid["M1z"], .75)) &
        (resid["M2z"] <= np.quantile(resid["M2z"], .25)) &
        (resid["M3z"] <= np.quantile(resid["M3z"], .25)))
base = df.state_high.values
cs = {}
for name, sel in [("plasma", df.cell_type == "Plasma cells"),
                  ("fibroblast", df.cell_type == "Fibroblasts"),
                  ("epithelial", df.cell_type == "Epithelial cells")]:
    a = int((sel & mask).sum()); b = int((sel & ~mask).sum())
    c = int((~sel & mask).sum()); d = int((~sel & ~mask).sum())
    o, lo, hi = or_ci(a, b, c, d)
    cs[name] = {"OR": o, "CI": [lo, hi], "state_in": a, "n": int(sel.sum())}
res["technical_residualised_definition"] = {
    "pct": float(100 * mask.mean()), "n": int(mask.sum()),
    "jaccard_vs_base": float((mask & base).sum() / (mask | base).sum()),
    "carrier": cs,
    "corr_resid_vs_orig": {m: float(np.corrcoef(resid[m], df[m])[0, 1]) for m in resid}}
res["state_rate_in_plasma_by_nFeature_tertile"] = (
    df[df.cell_type == "Plasma cells"].assign(t=pd.qcut(df[df.cell_type == "Plasma cells"].nFeature_RNA, 3, labels=False))
      .groupby("t").state_high.agg(["mean", "size"]).to_dict("index"))
res["state_rate_in_fibro_by_nFeature_tertile"] = (
    df[df.cell_type == "Fibroblasts"].assign(t=pd.qcut(df[df.cell_type == "Fibroblasts"].nFeature_RNA, 3, labels=False))
      .groupby("t").state_high.agg(["mean", "size"]).to_dict("index"))

with open(os.path.join(OUT, "A3b_depth_adjustment.json"), "w", encoding="utf-8") as f:
    json.dump(res, f, indent=2, default=float)
print(json.dumps(res, indent=2, default=float))
