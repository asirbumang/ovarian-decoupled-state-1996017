# -*- coding: utf-8 -*-
"""
A5: Clinical layer clarification.
 (i)  continuous Cox (per 1 SD of composition-residualised decouple score) with
      and without covariate adjustment in TCGA-OV;
 (ii) tertile log-rank kept separate from the continuous model;
 (iii) fixed-effect vs random-effects (DL and REML) pooling of the four
      external cohorts, with I2, tau2 and a 95% prediction interval, given
      I2 = 99.5%.
"""
import json, os, math
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test, logrank_test
from scipy import stats

ROOT = r"D:\子宫内膜异位症2"
OUT = os.path.join(ROOT, "R2_revision", "04_analysis", "results")
meta = pd.read_csv(os.path.join(ROOT, "analysis", "L0_bulk", "results", "track2",
                                "decreg_cox_meta.csv"))
tcga = pd.read_csv(os.path.join(ROOT, "analysis", "L0_bulk", "results", "track2",
                                "TCGA_OV_residual.csv"))

res = {}

# ---------------------------------------------------------------- TCGA models
d = tcga[(tcga.os_days > 0) & tcga.event.notna()].copy()
d["dec_res_z"] = (d.dec_res - d.dec_res.mean()) / d.dec_res.std()
d["age10"] = d.age_at_dx / 10.0
d["resid_gt1cm"] = d.residual_tumor.isin([">20 mm", "1-10 mm", "11-20 mm"]).astype(int)
d["resid_macro"] = d.residual_tumor.isin([">20 mm"]).astype(int)
d["age10"] = d["age10"].fillna(d["age10"].median())
d["resid_gt1cm"] = d["resid_gt1cm"].fillna(0)


def cox(df, cols, label):
    c = CoxPHFitter()
    c.fit(df[["os_days", "event"] + cols], duration_col="os_days", event_col="event")
    s = c.summary.loc[cols[0]]
    return {"model": label, "n": int(len(df)), "events": int(df.event.sum()),
            "HR_per_SD": float(s["exp(coef)"]), "CI_low": float(s["exp(coef) lower 95%"]),
            "CI_high": float(s["exp(coef) upper 95%"]), "p": float(s["p"]),
            "covariates": cols[1:]}


res["tcga_continuous"] = [
    cox(d, ["dec_res_z"], "unadjusted"),
    cox(d, ["dec_res_z", "age10"], "age-adjusted"),
    cox(d, ["dec_res_z", "age10", "resid_gt1cm"], "age + residual disease adjusted"),
    cox(d, ["dec_res_z", "age10", "resid_macro"], "age + gross residual (>20 mm) adjusted")]

# tertile analysis, kept separate from the continuous model
d["tier3"] = pd.qcut(d.dec_res, 3, labels=["Low", "Mid", "High"])
lr = multivariate_logrank_test(d.os_days, d.tier3, d.event)
res["tcga_tertile"] = {"logrank_p": float(lr.p_value), "n": int(len(d)),
                       "n_by_tier": d.tier3.value_counts().sort_index().to_dict()}
ct = CoxPHFitter().fit(d.assign(high=(d.tier3 == "High").astype(int))[
    ["os_days", "event", "high"]], "os_days", "event")
s = ct.summary.loc["high"]
res["tcga_tertile_cox_high_vs_rest"] = {"HR": float(s["exp(coef)"]), "CI_low": float(s["exp(coef) lower 95%"]),
                                        "CI_high": float(s["exp(coef) upper 95%"]), "p": float(s["p"])}

# proportional hazards check for the adjusted model
from lifelines.statistics import proportional_hazard_test
cph = CoxPHFitter().fit(d[["os_days", "event", "dec_res_z", "age10", "resid_gt1cm"]], "os_days", "event")
ph = proportional_hazard_test(cph, d[["os_days", "event", "dec_res_z", "age10", "resid_gt1cm"]], time_transform="rank")
res["tcga_ph_test_p"] = {str(k): float(v) for k, v in ph.summary["p"].items()}

# ------------------------------------------------------------ meta-analysis
ext = meta[meta.cohort != "TCGA-OV"].copy()
ext["logHR"] = np.log(ext.hr_cont)
ext["se"] = (np.log(ext.ci_hi) - np.log(ext.ci_lo)) / 3.92
ext["w_fixed"] = 1 / ext.se ** 2


def dl_meta(loghr, se):
    w = 1 / se ** 2
    fe = (w * loghr).sum() / w.sum()
    Q = (w * (loghr - fe) ** 2).sum()
    k = len(loghr); C = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (Q - (k - 1)) / C)
    I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    wr = 1 / (se ** 2 + tau2)
    re = (wr * loghr).sum() / wr.sum(); sere = math.sqrt(1 / wr.sum())
    pred_se = math.sqrt(tau2 + sere ** 2)
    tcrit = stats.t.ppf(0.975, k - 2)
    return {"k": k, "FE_logHR": fe, "FE_HR": math.exp(fe),
            "FE_CI": (math.exp(fe - 1.96 * math.sqrt(1 / w.sum())), math.exp(fe + 1.96 * math.sqrt(1 / w.sum()))),
            "Q": Q, "Q_p": stats.chi2.sf(Q, k - 1), "I2_pct": I2, "tau2": tau2,
            "RE_logHR": re, "RE_HR": math.exp(re),
            "RE_CI": (math.exp(re - 1.96 * sere), math.exp(re + 1.96 * sere)),
            "RE_p": 2 * stats.norm.sf(abs(re / sere)),
            "prediction_interval": (math.exp(re - tcrit * pred_se), math.exp(re + tcrit * pred_se))}


res["meta_external"] = dl_meta(ext.logHR.values, ext.se.values)
res["meta_all5"] = dl_meta(np.log(meta.hr_cont.values), ((np.log(meta.ci_hi) - np.log(meta.ci_lo)) / 3.92).values)
res["meta_external"]["cohorts"] = ext[["cohort", "n", "events", "hr_cont", "ci_lo", "ci_hi"]].to_dict("records")

with open(os.path.join(OUT, "A5_clinical.json"), "w", encoding="utf-8") as f:
    json.dump(res, f, indent=2, default=float)
print(json.dumps(res, indent=2, default=float))
