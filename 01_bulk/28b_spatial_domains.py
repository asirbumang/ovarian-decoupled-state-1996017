# -*- coding: utf-8 -*-
"""28b v4: NMF max_iter=1000 + 样本级域-状态一致性汇总。"""
import os, sys, gzip, glob
import numpy as np, pandas as pd
from sklearn.decomposition import NMF
from scipy.stats import chi2_contingency, fisher_exact
sys.stdout.reconfigure(encoding='utf-8')

OUT = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial"
NMF_IN = os.path.join(OUT, "nmf_input")
K = 8
rows_all = []
for GSE in sys.argv[1:]:
    scores = pd.read_csv(os.path.join(OUT, f"{GSE}_spot_scores.csv"))
    print("==", GSE, flush=True)
    domrows = []
    for fp in sorted(glob.glob(os.path.join(NMF_IN, f"{GSE}_*.csv.gz"))):
        code = os.path.basename(fp).replace(f"{GSE}_", "").replace(".csv.gz", "")
        logc = pd.read_csv(fp, index_col=0)
        nmf = NMF(n_components=K, init="nndsvda", random_state=2026, max_iter=1000)
        W = nmf.fit_transform(logc.values.T)
        dom = W.argmax(axis=1)
        for b, dd in zip(logc.columns, dom):
            domrows.append({"sample": code, "spot": b, "domain": int(dd)})
    d = pd.DataFrame(domrows)
    sc = scores.merge(d, on=["sample", "spot"], how="left")
    sc.to_csv(os.path.join(OUT, f"{GSE}_domains.csv"), index=False)
    # 样本级：每样本域大小 + 状态率；检查"最大域=最低状态率?"和"小域富集"
    for sm, sub in sc.groupby("sample"):
        dd = sub.groupby("domain")["state_high"].agg(["size", "mean"]).reset_index()
        dd["frac"] = dd["size"] / dd["size"].sum()
        dd = dd.sort_values("size", ascending=False)
        big = dd.iloc[0]
        small_enriched = dd[(dd["frac"] < 0.15) & (dd["mean"] > 0.03)]
        rows_all.append({"GSE": GSE, "sample": sm, "n_spots": len(sub),
                         "big_domain_frac": round(big["frac"], 3), "big_domain_rate": round(big["mean"], 4),
                         "small_enrich_n": len(small_enriched),
                         "max_small_rate": round(float(dd["mean"].max()), 4)})
        print(f"  {sm}: n={len(sub)} big_domain {big['frac']:.2f}rate{big['mean']:.3f} | small-enriched domains {len(small_enriched)}", flush=True)
res = pd.DataFrame(rows_all)
res.to_csv(os.path.join(OUT, "spatial_domain_state_summary_all.csv"), index=False)
print(res.groupby("GSE")["small_enrich_n"].agg(["mean", "sum"]), flush=True)
print("ALL DONE")
