# -*- coding: utf-8 -*-
"""GeoMx：AOI 级四模块打分 + 解耦状态 × ROI 类型/TLS（Cancer Cell 2024 范式对接）。"""
import pandas as pd, numpy as np, sys
sys.stdout.reconfigure(encoding='utf-8')
path = r"D:\子宫内膜异位症2\data_spatial\GSE276935\GSE276935_Processed.xlsx"
xl = pd.ExcelFile(path, engine="openpyxl")

## 1. 元数据（SegmentProperties 全列）
meta = xl.parse("SegmentProperties")
print("all cols:", list(meta.columns), flush=True)
print("n AOI:", len(meta), flush=True)

## 2. 靶点 → HUGO
tp = xl.parse("TargetProperties")
print("n targets:", len(tp), flush=True)

## 3. 计数矩阵
cm = xl.parse("TargetCountMatrix")
print("count matrix:", cm.shape, flush=True)
aoi_names = cm.columns[1:].tolist()
print("AOIs:", len(aoi_names), flush=True)

## 4. 模块基因集
sets = {}
with open(r"D:\子宫内膜异位症2\pipeline\mitoxy_four_modules.gmt", encoding="utf-8") as f:
    for line in f:
        p = line.strip().split("\t")
        sets[p[0]] = set(p[2:])
print("modules:", {k: len(v) for k, v in sets.items()}, flush=True)

## 靶点映射
sym2t = {}
for _, r in tp.iterrows():
    s = str(r.get("HUGOSymbol", "")).strip()
    if s and s != "nan":
        sym2t.setdefault(s, []).append(r["TargetName"])

def score_set(syms):
    tns = set()
    for g in syms:
        tns.update(sym2t.get(g, []))
    return sorted([i for i, t in enumerate(cm["TargetName"]) if t in tns])

sel_all = []
mod_idx = {}
for m, genes in sets.items():
    idx = score_set(genes)
    mod_idx[m] = idx
    sel_all.extend(idx)
sel_all = sorted(set(sel_all))
print("sel targets:", len(sel_all), flush=True)

X = cm.iloc[sel_all, 1:].values.astype(float)   # targets × AOI
names = [cm["TargetName"].iloc[i] for i in sel_all]
# AOI→gene 映射（score_set 顺序与 mod_idx 一致；简化：每模块 avg of its idx）
## log2 CPM
cpm = X / X.sum(axis=0, keepdims=True) * 1e6
logc = np.log2(cpm + 1)
res = pd.DataFrame({"aoi": aoi_names})
for m, idx in mod_idx.items():
    pos = [i for i, n in enumerate(names) if n in set(cm["TargetName"].iloc[idx])]
    # 直接按原始索引
    pos = [j for j, g in enumerate(sel_all) if g in idx]
    scores = logc[pos, :].mean(axis=0)
    res[m] = scores
    print(m, "targets:", len(pos), flush=True)
# decouple（AOI 间 z）
for m in list(sets.keys()):
    z = (res[m] - res[m].mean()) / res[m].std()
    res[f"z_{m}"] = z
res["decouple"] = res["z_MITOXY_M1_immune_activation"] - (res["z_MITOXY_M2_metabolic_GSH"] + res["z_MITOXY_M3_mito_cytoskeleton"]) / 2
q75 = res["MITOXY_M1_immune_activation"].quantile(.75); q25 = res["MITOXY_M2_metabolic_GSH"].quantile(.25); q25b = res["MITOXY_M3_mito_cytoskeleton"].quantile(.25)
res["state_high"] = (res["MITOXY_M1_immune_activation"] >= q75) & (res["MITOXY_M2_metabolic_GSH"] <= q25) & (res["MITOXY_M3_mito_cytoskeleton"] <= q25b)

## 5. 合并元数据（ROILabel 等）
meta_map = meta.set_index(meta.columns[0])
mm = meta_map.loc[res["aoi"]] if res["aoi"].isin(meta_map.index).all() else None
if mm is not None:
    for c in meta_map.columns:
        res["meta_" + str(c)] = meta_map[c].values
res.to_csv(r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial\GSE276935_aoi_scores.csv", index=False)
print("saved aoi scores:", res.shape, flush=True)
print("state_high:", int(res["state_high"].sum()), "/", len(res), flush=True)
print(res[["aoi", "decouple", "state_high"] + [c for c in res.columns if c.startswith("meta_ROI") or c.startswith("meta_RL") or c.startswith("meta_Se")]].head(8).to_string(), flush=True)
print("DONE geomx")
