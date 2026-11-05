# -*- coding: utf-8 -*-
"""GeoMx 段占比相关 + 汇总（28e2 补）。"""
import pandas as pd, sys
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
for p in [r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial\GSE276935_aoi_scores_annot.csv",
          r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial\GSE276935_aoi_scores.csv"]:
    try:
        d = pd.read_csv(p)
        print("loaded", p, d.shape)
        break
    except Exception as e:
        print("fail", e)
for frac in ["CD20", "Immune", "Stroma", "Tumor", "other"]:
    c = "meta_" + frac
    if c in d.columns:
        rho, p = stats.spearmanr(d["decouple"], d[c])
        print(f"decouple vs {frac}: rho={rho:.3f} p={p:.2e}", flush=True)
## 患者级总结已存；补：state_high 定义阈值敏感（宽松：decouple>P75）
d["state_hi_loose"] = d["decouple"] > d["decouple"].quantile(.75)
print("\nloose-state by seg:", dict(d.groupby("seg")["state_hi_loose"].sum()), flush=True)
ktab = pd.crosstab(d["seg"], d["state_hi_loose"])
print(ktab, flush=True)
