# -*- coding: utf-8 -*-
"""28e2: GeoMx AOI 级：段类型 × 解耦 + CD20(TLS代理)/免疫占比 相关 + 患者级。"""
import pandas as pd, numpy as np, sys
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
d = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial\GSE276935_aoi_scores.csv")
print("AOIs:", len(d), flush=True)
d["seg"] = d["meta_SegmentLabel"]
d["patient"] = d["aoi"].str.extract(r"^(\d+)")[0]
d["site"] = d["aoi"].str.extract(r"(\w+?)0\d+(CD20|Immune|Stroma|Tumor|other|)$")[0]
print(d.groupby("seg").size(), flush=True)

## 1. decouple by segment (Kruskal + Tumor vs Immune/Stroma)
print("\n== decouple by segment ==")
for s, sub in d.groupby("seg"):
    print(f"  {s:8s} n={len(sub):3d} med={sub['decouple'].median():.3f}", flush=True)
k = stats.kruskal(*[sub["decouple"].values for _, sub in d.groupby("seg")])
print("Kruskal p =", k.pvalue, flush=True)
w = stats.mannwhitneyu(d.loc[d.seg == "Tumor", "decouple"], d.loc[d.seg == "Immune", "decouple"])
print("Tumor vs Immune MWU p =", w.pvalue, flush=True)

## 2. CD20 fraction (TLS proxy) and Immune fraction vs decouple (hepatocytes: 全段/肿瘤段)
for frac in ["CD20", "Immune", "Stroma", "Tumor"]:
    if frac in d.columns:
        rho, p = stats.spearmanr(d["decouple"], d[frac])
        print(f"decouple vs {frac}-fraction: rho={rho:.3f} p={p:.2e}", flush=True)

## 3. state_high 分布
print("\nstate_high:", int(d["state_high"].sum()), " by seg:", dict(d.groupby("seg")["state_high"].sum()), flush=True)

## 4. 患者级（多 AOI/患者 → 患者中位；Tumor 段只算肿瘤患者）
pat = d.groupby("patient")["decouple"].agg(["size", "median"]).reset_index()
print("\npatients:", len(pat), " median n AOI:", pat["size"].median(), flush=True)
print(pat.sort_values("median").head(10).to_string(index=False), flush=True)
d.to_csv(r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial\GSE276935_aoi_scores_annot.csv", index=False)
print("DONE 28e2")
