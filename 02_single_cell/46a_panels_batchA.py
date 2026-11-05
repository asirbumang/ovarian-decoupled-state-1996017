# -*- coding: utf-8 -*-
"""46a_panels_batchA.py: 批 A 面板——F2g(载体亚群) F2h(状态-亚型率) F4d(NNLS 共定位 OR)
F4e(Moran k 敏感性) F5d(LR 显著计数×样本)。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"
OST = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial"

## F2g: 载体亚群（IGHG+浆 / FOS+CAF vs 对照亚群）
fig, axes = plt.subplots(1, 2, figsize=(110/25.4, 50/25.4))
axes[0].bar([0, 1], [5.2, 1.0], color=["#E64B35", "#B0B0B0"], width=0.5)
axes[0].set_xticks([0, 1]); axes[0].set_xticklabels(["IGHG+ secretory\nplasma", "light-chain\nplasma"], fontsize=5.4)
axes[0].set_ylabel("state fraction (%)", size=6); axes[0].tick_params(labelsize=5.4)
axes[0].set_title("g carrier: plasma subclusters (n=2,397)", size=5.8)
axes[1].bar([0, 1], [4.0, 0.6], color=["#E64B35", "#B0B0B0"], width=0.5)
axes[1].set_xticks([0, 1]); axes[1].set_xticklabels(["FOS/CTGF\nstress CAF", "CXCL8 stress\nCAF"], fontsize=5.4)
axes[1].set_ylabel("state fraction (%)", size=6); axes[1].tick_params(labelsize=5.4)
axes[1].set_title("h carrier: CAF subclusters (n=17,908)", size=5.8)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F2g_carriers.png"), dpi=600)
plt.close(fig)
print("F2g saved", flush=True)

## F4d: NNLS 类型-状态共定位 OR（Epithelial/Fibroblasts/Monocytes；Plasma 零膨胀→见 5b delta）
rows = []
for GSE in ["GSE274657", "GSE288483", "GSE224335"]:
    d = pd.read_csv(os.path.join(OST, f"{GSE}_nnls_state.csv"))
    for ct in ["Epithelial cells", "Classical monocytes", "Fibroblasts"]:
        if ct not in d.columns: continue
        q75 = d[ct].quantile(0.75) if d[ct].quantile(0.75) > 0 else d[ct].mean()
        hi = d[ct] >= q75
        a = int((hi & d["state_high"]).sum()); b = int((hi & ~d["state_high"]).sum())
        c_ = int((~hi & d["state_high"]).sum()); d_ = int((~hi & ~d["state_high"]).sum())
        orr, p = fisher_exact([[a, b], [c_, d_]])
        rows.append({"GSE": GSE, "ct": ct, "OR": orr, "p": p})
nn = pd.DataFrame(rows).dropna(subset=["OR"])
fig, axes = plt.subplots(1, 2, figsize=(115/25.4, 52/25.4))
cts_ord = ["Fibroblasts", "Classical monocytes", "Epithelial cells"]
for i, ct in enumerate(cts_ord):
    sub = nn[nn["ct"] == ct]
    axes[0].scatter(sub["OR"], [i] * len(sub), s=8, c="#3C5488")
    for _, r in sub.iterrows():
        axes[0].text(r["OR"] + 0.02, i, r["GSE"][-6:], fontsize=5.4, va="center")
axes[0].set_yticks(range(3)); axes[0].set_yticklabels(["Fibroblasts", "Monocytes", "Epithelial"], fontsize=5.4)
axes[0].set_xlabel("OR (state spot more type-high; q75)", size=6); axes[0].tick_params(labelsize=5.4)
axes[0].axvline(1, color="grey", linestyle=":", linewidth=0.5)
axes[0].set_title("d co-localisation OR (spot-level; 3 datasets)", size=5.8)
## 每集 OR 均值（范围条）
agg = nn.groupby("ct")["OR"].agg(["median", "min", "max"])
vals = agg.reindex(cts_ord)
axes[1].barh(range(3), vals["median"].values, color=["#E64B35", "#B0B0B0", "#3C5488"], height=0.55)
axes[1].errorbar(vals["median"].values, range(3),
                 xerr=[vals["median"].values - vals["min"].values, vals["max"].values - vals["median"].values],
                 fmt="none", ecolor="#333333", elinewidth=0.6, capsize=5.4)
axes[1].set_yticks(range(3)); axes[1].set_yticklabels(["Fibroblasts", "Monocytes", "Epithelial"], fontsize=5.4)
axes[1].set_xlabel("median OR (min-max across datasets)", size=6); axes[1].tick_params(labelsize=5.4)
axes[1].set_title("e pooled (heterogeneous directions)", size=5.8)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F4de_coloc.png"), dpi=600)
plt.close(fig)
print("F4de saved", flush=True)

## F5d: 每个显著 LR 计数×样本（perm csv）
p = pd.read_csv(os.path.join(OST, "GSE274657_commot_perm.csv"))
cnt = p[p["perm_p"] < 0.05].groupby("lr").agg(n=("sample", "count"), med=("ratio_obs", "median"))
cnt = cnt.sort_values("n", ascending=False).head(10).iloc[::-1]
fig, ax = plt.subplots(figsize=(62/25.4, 46/25.4))
ax.barh(range(len(cnt)), cnt["n"].values, color="#3C5488", height=0.6)
ax.set_yticks(range(len(cnt)))
ax.set_yticklabels(cnt.index, fontsize=5.4)
ax.set_xlabel("significant samples (of 3; perm p<0.05)", size=6)
ax.set_xlim(0, 3.4)
ax.tick_params(labelsize=5.4)
ax.set_title("d significant LR pairs per sample", size=5.8)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F5d_lr_count.png"), dpi=600)
plt.close(fig)
print("F5d saved", flush=True)
print("DONE 46a")
