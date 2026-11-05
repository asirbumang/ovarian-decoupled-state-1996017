# -*- coding: utf-8 -*-
"""A8: UCell-based re-scoring and ICC45 deconvolution (full/shared/unique)."""
import json, os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = r"D:\子宫内膜异位症2"
OUT = os.path.join(ROOT, "R2_revision", "04_analysis", "results")
sc = pd.read_csv(os.path.join(ROOT, "analysis", "L1_scRNA", "results",
                              "module_scores", "scores_GSE266577.csv"))
uc = pd.read_csv(os.path.join(OUT, "A4_ucell_python.csv"))
uc = uc.rename(columns={c: c + "_UCell" for c in uc.columns if c != "cell"})
df = sc.merge(uc, on="cell", how="inner")
res = {"n_cells": int(len(df))}

# concordance with AddModuleScore
res["spearman_ucell_vs_addmodulescore"] = {
    m: {"rho": float(stats.spearmanr(df[m], df[m + "_UCell"]).statistic),
        "pearson": float(np.corrcoef(df[m], df[m + "_UCell"])[0, 1])}
    for m in ["M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"]}

# state from UCell
q1 = np.quantile(df["M1_immune_activation_UCell"], .75)
q2 = np.quantile(df["M2_metabolic_GSH_UCell"], .25)
q3 = np.quantile(df["M3_mito_cytoskeleton_UCell"], .25)
state_u = ((df["M1_immune_activation_UCell"] >= q1) &
           (df["M2_metabolic_GSH_UCell"] <= q2) &
           (df["M3_mito_cytoskeleton_UCell"] <= q3)).values
base = df.state_high.values


def orci(mask):
    out = {}
    for name, sel in [("plasma", df.cell_type == "Plasma cells"),
                      ("fibroblast", df.cell_type == "Fibroblasts"),
                      ("epithelial", df.cell_type == "Epithelial cells")]:
        a = int((sel & mask).sum()); b = int((sel & ~mask).sum())
        c = int((~sel & mask).sum()); d = int((~sel & ~mask).sum())
        orr, p = stats.fisher_exact([[a, b], [c, d]])
        out[name] = {"OR": float(orr), "p": float(p), "state": a, "n": int(sel.sum())}
    return out


res["ucell_state"] = {
    "pct": float(100 * state_u.mean()), "n": int(state_u.sum()),
    "jaccard_vs_primary": float((state_u & base).sum() / (state_u | base).sum()),
    "carrier": orci(state_u)}
res["primary_state_carrier_for_comparison"] = orci(base)
df["nf_q"] = pd.qcut(df.nFeature_RNA, 5, labels=False)
res["ucell_state_by_nFeature_quintile"] = (df.groupby("nf_q")
                                           .apply(lambda g: float(state_u[g.index].mean()), include_groups=False)
                                           .to_dict())
res["addmodulescore_state_by_nFeature_quintile"] = df.groupby("nf_q").state_high.mean().to_dict()

# ---------------------------------------------------------------- ICC45
icc_sets = ["ICC45_full", "ICC45_shared28", "ICC45_unique17"] + [f"RAND17_{i}" for i in range(1, 6)]
icc_res = {}
for s in icc_sets:
    r = stats.spearmanr(df.decouple, df[s + "_UCell"])
    # AUC for the primary state label
    u = stats.mannwhitneyu(df.loc[base, s + "_UCell"], df.loc[~base, s + "_UCell"])
    auc = u.statistic / (base.sum() * (~base).sum())
    icc_res[s] = {"spearman_with_decouple": float(r.statistic), "p": float(r.pvalue),
                  "AUC_for_state": float(auc)}
res["icc45"] = icc_res
res["icc45_random_summary"] = {
    "spearman_mean": float(np.mean([icc_res[f"RAND17_{i}"]["spearman_with_decouple"] for i in range(1, 6)])),
    "spearman_sd": float(np.std([icc_res[f"RAND17_{i}"]["spearman_with_decouple"] for i in range(1, 6)], ddof=1))}

with open(os.path.join(OUT, "A8_ucell_icc45.json"), "w", encoding="utf-8") as f:
    json.dump(res, f, indent=2, default=float)
print(json.dumps(res, indent=2, default=float))
