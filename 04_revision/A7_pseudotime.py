# -*- coding: utf-8 -*-
"""
A7: Trajectory robustness.
 (i)  compare the state-along-pseudotime profile under the two roots
      (min-decouple vs max-decouple cluster) for fibroblasts and plasma cells;
 (ii) test the trend formally (Cochran-Armitage / logistic on bin);
 (iii) check consistency across patients (per-patient Spearman of state vs dpt,
      and per-patient peak bin).
"""
import json, os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = r"D:\子宫内膜异位症2"
PT = os.path.join(ROOT, "analysis", "L1_scRNA", "results", "pseudotime")
OUT = os.path.join(ROOT, "R2_revision", "04_analysis", "results")
scores = pd.read_csv(os.path.join(ROOT, "analysis", "L1_scRNA", "results",
                                  "module_scores", "scores_GSE266577.csv"),
                     usecols=["cell", "publication_patient_code_final"])
scores["patient"] = scores["publication_patient_code_final"].astype(str)

res = {}
for lin in ["fib", "plasma", "epi"]:
    curves = {}
    for root in ["min", "max"]:
        f = os.path.join(PT, f"PT_{lin}_{root}_decouple_cells.csv")
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f)
        d = d.merge(scores[["cell", "patient"]], on="cell", how="left")
        d["bin"] = pd.qcut(d.dpt, 7, labels=False, duplicates="drop")
        g = d.groupby("bin").agg(n=("state", "size"), k=("state", "sum"),
                                 dpt=("dpt", "median"))
        g["frac"] = g.k / g.n
        # Wilson CI
        z = 1.96
        lo, hi = [], []
        for k, n in zip(g.k, g.n):
            p = k / n
            den = 1 + z**2 / n
            c = (p + z**2 / (2*n)) / den
            h = z * np.sqrt(p*(1-p)/n + z**2/(4*n*n)) / den
            lo.append(c-h); hi.append(c+h)
        g["lo"], g["hi"] = lo, hi
        # trend: logistic regression state ~ bin
        import statsmodels.api as sm
        X = sm.add_constant(g.index.values.astype(float))
        y = g.k.values / g.n.values
        try:
            w = sm.WLS(y, X, weights=g.n.values).fit()
            trend_p = float(w.pvalues[1]); trend_beta = float(w.params[1])
        except Exception:
            trend_p, trend_beta = np.nan, np.nan
        # per-patient trend
        pt = []
        for pid, sub in d.groupby("patient"):
            if sub.state.sum() < 3 or len(sub) < 100:
                continue
            rho = stats.spearmanr(sub.dpt, sub.state.astype(int))
            pt.append({"patient": pid, "n": len(sub), "k": int(sub.state.sum()),
                       "rho": float(rho.statistic), "p": float(rho.pvalue)})
        pt = pd.DataFrame(pt)
        curves[root] = {"bins": g.reset_index().to_dict("records"),
                        "trend_beta_per_bin": trend_beta, "trend_p": trend_p,
                        "peak_bin": int(g.frac.idxmax()), "peak_frac": float(g.frac.max()),
                        "global_frac": float(d.state.mean()),
                        "per_patient": pt.to_dict("records") if len(pt) else [],
                        "n_patients_testable": int(len(pt)),
                        "n_patients_positive": int((pt.rho > 0).sum()) if len(pt) else 0,
                        "n_patients_negative": int((pt.rho < 0).sum()) if len(pt) else 0,
                        "median_patient_rho": float(pt.rho.median()) if len(pt) else np.nan}
    if len(curves) == 2:
        a = pd.DataFrame(curves["min"]["bins"]).set_index("bin").frac
        b = pd.DataFrame(curves["max"]["bins"]).set_index("bin").frac
        common = a.index.intersection(b.index)
        rho = stats.spearmanr(a.loc[common], b.loc[common])
        curves["root_concordance"] = {"spearman": float(rho.statistic),
                                      "p": float(rho.pvalue), "n_bins": int(len(common))}
    res[lin] = curves

with open(os.path.join(OUT, "A7_pseudotime.json"), "w", encoding="utf-8") as f:
    json.dump(res, f, indent=2, default=float)

for lin, c in res.items():
    print("==", lin)
    for root in ["min", "max"]:
        if root in c:
            x = c[root]
            print(f"  {root}: peak bin {x['peak_bin']} ({x['peak_frac']*100:.1f}%), "
                  f"trend p={x['trend_p']:.3g}, per-patient +/{x['n_patients_positive']} "
                  f"-/{x['n_patients_negative']} median rho={x['median_patient_rho']:.2f}")
    if "root_concordance" in c:
        print("  root concordance:", c["root_concordance"])
