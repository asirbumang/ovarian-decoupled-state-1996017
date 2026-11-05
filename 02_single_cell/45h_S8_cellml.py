# -*- coding: utf-8 -*-
"""45h_S8_cellml.py: S8 补充图——细胞级判别 ML：a 30 组合 AUC、b 签名频率 top15、c 54 基因签名交叉注。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
OUT = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b12"
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs\supp"

z = pd.read_csv(os.path.join(OUT, "modelzoo_cell_auc.csv"))
z = z.sort_values("auc").iloc[::-1]
fig, axes = plt.subplots(1, 2, figsize=(150/25.4, 56/25.4))
cols = ["#DC0000" if m in ("SVM", "LR_l2", "LR_l1") else "#B0B0B0" for m in z["clf"]]
axes[0].barh(range(len(z)), z["auc"], color="#3C5488", height=0.7)
axes[0].set_yticks(range(len(z)))
axes[0].set_yticklabels([f"{r.fs}|{r.clf}" for _, r in z.iterrows()], fontsize=5.4)
axes[0].set_xlabel("AUC (GroupKFold by patient, 5-fold)", size=6); axes[0].tick_params(labelsize=5.4)
axes[0].axvline(z["auc"].min(), color="grey", linestyle=":", linewidth=0.5)
axes[0].set_title("a cell-level zoo (30 combos; Dummy 0.5)", size=6)
sm = pd.read_csv(os.path.join(OUT, "sig_freq.csv"))
sm = sm.sort_values("freq", ascending=False).head(15).iloc[::-1]
axes[1].barh(range(len(sm)), sm["freq"], color="#E64B35", height=0.6)
axes[1].set_yticks(range(len(sm)))
axes[1].set_yticklabels(sm["gene"], fontsize=5.4)
axes[1].set_xlabel("bootstrap frequency (30x)", size=6); axes[1].tick_params(labelsize=5.4)
axes[1].set_title("b bootstrap L1 signature (54 genes; LR AUC=0.997)", size=6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "S8_cellml.png"), dpi=600)
plt.close(fig)
print("S8 saved", flush=True)
print("DONE 45h")
