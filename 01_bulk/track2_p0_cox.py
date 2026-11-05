# >=- coding: utf-8 -*-
"""P0 收尾：TCGA-OV deconvolution-adjusted（dec_res）× OS。
比较：raw decouple vs dec_res（肿瘤组分校正后）。"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import pandas as pd, numpy as np
from lifelines import CoxPHFitter
d = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L0_bulk\results\track2\TCGA_OV_residual.csv")
print("samples:", len(d), flush=True)
d = d[d["event"].notna() & d["os_days"].notna()].copy()
print("with outcome:", len(d), flush=True)
for col, label in [("decouple", "RAW decouple"), ("dec_res", "ADJUSTED dec_res")]:
    x = d[[col, "os_days", "event"]].copy()
    x = x[x[col].notna()]
    ## z 标准化后 Cox（每 SD 的 HR）
    x[col] = (x[col] - x[col].mean()) / x[col].std()
    cph = CoxPHFitter()
    try:
        cph.fit(x, duration_col="os_days", event_col="event", formula=f"Q('{col}')")
        print(f"{label}: HR={cph.hazard_ratios_.iloc[0]:.3f} p={cph.summary['p'].iloc[0]:.4f} n={len(x)}", flush=True)
    except Exception as e:
        print(label, "ERR", repr(e)[:80], flush=True)
## 高低分组（median）简化逻辑
for col, label in [("decouple", "RAW"), ("dec_res", "ADJ")]:
    x = d[["sample", col, "os_days", "event"]].dropna()
    hi = x[col] > x[col].median()
    x["grp"] = np.where(hi, "High", "Low")
    cph = CoxPHFitter()
    try:
        cph.fit(x[["grp", "os_days", "event"]], duration_col="os_days", event_col="event", formula="C(grp, Treatment('Low'))")
        print(f"{label} med-split: HR(High/Low)={cph.hazard_ratios_.iloc[0]:.3f} p={cph.summary['p'].iloc[0]:.4f}", flush=True)
    except Exception as e:
        print(label, "ERR2", repr(e)[:80], flush=True)
