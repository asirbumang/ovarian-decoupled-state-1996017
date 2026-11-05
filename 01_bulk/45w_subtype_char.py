# -*- coding: utf-8 -*-
"""45w_subtype_char.py: 亚型特征刻画（JTM 式步骤①）——亚型×6 特征 heatmap + 临床表 + 状态/轴关联。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, chi2_contingency
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"
d = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L0_bulk\results\track2\tcga_subtype_consensus.csv", low_memory=False)
d = d[d["subtype"].notna()]
feats = ["decouple", "M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton", "M4_lipid_perox_rupture", "ICC45"]
rows = []
for f in feats:
    a = d.loc[d["subtype"] == "high", f].dropna()
    b = d.loc[d["subtype"] == "low", f].dropna()
    p = mannwhitneyu(a, b).pvalue
    rows.append({"feature": f, "high_med": a.median(), "low_med": b.median(), "p": p, "direction": "+" if a.median() > b.median() else "-"})
fd = pd.DataFrame(rows)
print(fd.round(4).to_string(index=False), flush=True)
fd.to_csv(os.path.join(QRY, "S11_feat_diff.csv"), index=False)
## 临床特征表（x²）
clin_rows = []
for col, label in [("histological_type", "histology"), ("residual_tumor", "residual")]:
    tab = pd.crosstab(d["subtype"], d[col])
    chi, p, _, _ = chi2_contingency(tab)
    clin_rows.append({"clinical": label, "p": p})
    print(label, "chi2 p=", round(p, 4), flush=True)
age = mannwhitneyu(d.loc[d["subtype"] == "high", "age_at_dx"].dropna(),
                   d.loc[d["subtype"] == "low", "age_at_dx"].dropna()).pvalue
print("age p=", round(age, 4), flush=True)
## 图：亚型×特征差异热图（标准化 z 后 high-low）
fig, axes = plt.subplots(1, 3, figsize=(150/25.4, 50/25.4))
mz = d[feats].astype(float)
zs = (mz - mz.mean()) / (mz.std() + 1e-9)
order = np.argsort(d["subtype"].values)
im = axes[0].imshow(zs.values[order].T, cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
axes[0].set_yticks(range(len(feats))); axes[0].set_yticklabels(feats, fontsize=5.4)
axes[0].set_xticks([])
axes[0].axvline(np.sum(d["subtype"].values[order] == "high") - 0.5, color="white", linewidth=1.0)
axes[0].set_title("a features by subtype (left=high)", size=6)
plt.colorbar(im, ax=axes[0], fraction=0.03, pad=0.01)
fb = fd.sort_values("direction", key=lambda s: -s.map({"high": 1, "low": 1}).astype(float) * 0 + fd["p"]).iloc[::-1]
axes[1].barh(range(len(fb)), -np.log10(fb["p"].clip(lower=1e-300)), color=["#DC0000" if p < 0.05 else "#B0B0B0" for p in fb["p"]], height=0.6)
axes[1].set_yticks(range(len(fb))); axes[1].set_yticklabels(fb["feature"], fontsize=5.4)
axes[1].set_xlabel("-log10(p)", size=6); axes[1].tick_params(labelsize=5.4)
axes[1].set_title("b feature tests (high vs low)", size=6)
clin = pd.DataFrame(clin_rows + [{"clinical": "age", "p": age}])
axes[2].bar(range(len(clin)), -np.log10(clin["p"].clip(lower=1e-300)), color=["#3C5488"] * len(clin), width=0.5)
axes[2].set_xticks(range(len(clin))); axes[2].set_xticklabels(clin["clinical"], fontsize=5.4)
axes[2].set_ylabel("-log10(p)", size=6); axes[2].tick_params(labelsize=5.4)
axes[2].axhline(-np.log10(0.05), color="#DC0000", linestyle="--", linewidth=0.5)
axes[2].set_title("c clinical features (x2 / MWU)", size=6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "S11_subtype_char.png"), dpi=600)
plt.close(fig)
print("S11 saved", flush=True)
print("DONE 45w")
