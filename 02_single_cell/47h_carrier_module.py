# -*- coding: utf-8 -*-
"""47h_carrier_module.py: 载体亚群内模块结构（堵"定义偏倚"）——IGHG+浆 vs 轻链浆；FOS+CAF vs CXCL8+CAF。"""
import pandas as pd, numpy as np, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
sys.stdout.reconfigure(encoding="utf-8")
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
sub = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs\F12_plasma_caf_sub.csv")
sc = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577.csv",
                 usecols=["cell", "M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"], low_memory=False)
for m in ["M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"]:
    sc[m] = pd.to_numeric(sc[m], errors="coerce").astype(float)
m = sub.merge(sc, on="cell")
print(m.columns.tolist(), m["sub"].dtype, sorted(m["sub"].unique()), flush=True)
fig, axes = plt.subplots(1, 2, figsize=(115/25.4, 50/25.4))
sets = [("Plasma cells", 0, "IGHG+ sec. plasma", "light-chain plasma"),
        ("Fibroblasts", 1, "FOS/CTGF stress CAF", "CXCL8 stress CAF")]
for ax, (lin, s_carrier, lab1, lab2) in zip(axes, sets):
    d = m[m["lineage"] == lin]
    a = d[np.asarray(d["sub"]) == s_carrier]
    b = d[np.asarray(d["sub"]) != s_carrier]
    vals = []
    for mm in ["M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"]:
        p = mannwhitneyu(a[mm], b[mm]).pvalue
        vals.append((a[mm].median(), b[mm].median(), "<1e-300" if p == 0 else "%.1e" % p))
    x = np.arange(3); w = 0.36
    ax.bar(x - w / 2, [v[0] for v in vals], w, color="#DC0000", label=lab1)
    ax.bar(x + w / 2, [v[1] for v in vals], w, color="#B0B0B0", label=lab2)
    print(lin, "| carrier-med/other-med:", ["%.3f/%.3f (%s)" % (v[0], v[1], v[2]) for v in vals], flush=True)
    for i, v in enumerate(vals):
        ax.text(i, max(v[0], v[1]) + 0.02, v[2], ha="center", fontsize=5.4)
    ax.set_xticks(x); ax.set_xticklabels(["M1", "M2", "M3"], fontsize=6)
    ax.set_ylabel("module score (median)", size=6); ax.tick_params(labelsize=5.4)
    ax.legend(fontsize=5.4, frameon=False, loc="upper right")
    ax.set_title(f"{lab1} vs {lab2} (n={len(a)}/{len(b)})", size=6)
fig.tight_layout()
fig.savefig(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs\F2g2_module_structure.png", dpi=600)
plt.close(fig)
print("F2g2 saved", flush=True)
print("DONE 47h")
