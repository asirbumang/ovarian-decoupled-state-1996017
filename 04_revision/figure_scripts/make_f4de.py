# -*- coding: utf-8 -*-
"""Fig 4d: regenerate the two-subplot panel without the stray 'e ' title prefix."""
import numpy as np, pandas as pd, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"
OST = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial"

rows = []
for GSE in ["GSE274657", "GSE288483", "GSE224335"]:
    d = pd.read_csv(os.path.join(OST, f"{GSE}_nnls_state.csv"))
    for ct in ["Epithelial cells", "Classical monocytes", "Fibroblasts"]:
        if ct not in d.columns:
            continue
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
axes[0].set_title("co-localisation OR (spot-level; 3 datasets)", size=5.8)
agg = nn.groupby("ct")["OR"].agg(["median", "min", "max"])
vals = agg.reindex(cts_ord)
axes[1].barh(range(3), vals["median"].values, color=["#E64B35", "#B0B0B0", "#3C5488"], height=0.55)
axes[1].errorbar(vals["median"].values, range(3),
                 xerr=[vals["median"].values - vals["min"].values, vals["max"].values - vals["median"].values],
                 fmt="none", ecolor="#333333", elinewidth=0.6, capsize=5.4)
axes[1].set_yticks(range(3)); axes[1].set_yticklabels(["Fibroblasts", "Monocytes", "Epithelial"], fontsize=5.4)
axes[1].set_xlabel("median OR (min-max across datasets)", size=6); axes[1].tick_params(labelsize=5.4)
axes[1].set_title("pooled (heterogeneous directions)", size=5.8)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F4de_coloc.png"), dpi=600)
plt.close(fig)
print("F4de saved", flush=True)
