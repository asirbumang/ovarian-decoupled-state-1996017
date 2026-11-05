# -*- coding: utf-8 -*-
"""45p_fig6_main_panels.py: 图6 正文化面板 d(Mime) e(45c+Dummy) f(Scissor OR)。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"

## d: Mime 训练 vs 验证（7 个代表模型条形对）
p = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L0_bulk\results\track2\mime\mime101_summary.csv").set_index("Model")
sel = ["RSF", "RSF + CoxBoost", "Lasso + survival-SVM", "CoxBoost + GBM", "StepCox[backward] + survival-SVM", "GBM"]
sel = [s for s in sel if s in p.index]
ss = p.loc[sel]
fig, axes = plt.subplots(1, 3, figsize=(150/25.4, 56/25.4))
x = np.arange(len(sel)); w = 0.36
axes[0].bar(x - w/2, ss["train"], w, color="#3C5488", label="train (TCGA 425)")
axes[0].bar(x + w/2, ss["val"], w, color="#E64B35", label="validation (GSE26193 107)")
axes[0].axhline(0.5, color="grey", linestyle=":", linewidth=0.5)
axes[0].set_xticks(x); axes[0].set_xticklabels([s[:24] for s in sel], fontsize=5.4, rotation=30, ha="right")
axes[0].set_ylabel("C-index", size=6); axes[0].tick_params(labelsize=5.4)
axes[0].set_ylim(0.45, 0.95)
axes[0].legend(fontsize=5.4, frameon=False, loc="upper left")
axes[0].set_title("d Mime 101: train vs validation (representative)", size=5.8)
## e: 45c 交叉验证（5 模型 + Dummy）
c = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b12\survival101_cindex.csv")
sel2 = ["PCA30", "RFimp", "MI200", "None"]
s2 = c[c["fs"].isin(sel2) & (c["model"] == "Cox_l2")].set_index("fs").reindex(sel2)
db = c[(c["model"] == "Dummy")]
v2 = [s2.loc[s2.index[i], "cindex"] for i in range(len(sel2))] + [db["cindex"].iloc[0]]
lab2 = [f"{s}+Cox_l2" for s in sel2] + ["Dummy"]
axes[1].bar(range(5), v2, color=["#3C5488"]*4 + ["#B0B0B0"], width=0.55)
axes[1].axhline(0.5, color="grey", linestyle=":", linewidth=0.5)
axes[1].set_xticks(range(5)); axes[1].set_xticklabels(lab2, fontsize=5.4, rotation=30, ha="right")
axes[1].set_ylabel("C-index (TCGA 5-fold)", size=6); axes[1].tick_params(labelsize=5.4)
axes[1].set_ylim(0.45, 0.65)
axes[1].set_title("e independent implementation cross-check", size=5.8)
## f: Scissor 重叠（OR + CI）
import math
a, b, cc, d_ = 134, 11860, 614, 47376
orr = (a * d_) / (b * cc)
se = math.sqrt(1/a + 1/b + 1/cc + 1/d_)
lo, hi = math.exp(math.log(orr) - 1.96*se), math.exp(math.log(orr) + 1.96*se)
axes[2].bar(0, orr, color="#B0B0B0", width=0.3)
axes[2].errorbar(0, orr, yerr=[[orr-lo], [hi-orr]], fmt="none", ecolor="#333333", elinewidth=0.7, capsize=5.4)
axes[2].axhline(1, color="#DC0000", linestyle="--", linewidth=0.6)
axes[2].set_xticks([0]); axes[2].set_xticklabels(["Scissor+ vs state"], fontsize=5.4)
axes[2].set_ylabel("OR (95% CI)", size=6); axes[2].tick_params(labelsize=5.4)
axes[2].set_ylim(0.4, 1.6)
axes[2].text(0.08, orr + 0.06, f"OR={orr:.2f} [%.2f-%.2f] p=0.17" % (lo, hi), fontsize=5.4)
axes[2].set_title("f Scissor (cell-OS association): no overlap", size=5.8)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F6def_main.png"), dpi=600)
plt.close(fig)
print("F6def saved", flush=True)
print("DONE 45p")
