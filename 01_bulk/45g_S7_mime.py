# -*- coding: utf-8 -*-
"""45g_S7_mime.py: S7 补充图——Mime 101（117 模型）：a 训练/验证 C-index 排序条形；b gap 直方；c 验证 top 条。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
OUT = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\track2\mime"
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs\supp"
p = pd.read_csv(os.path.join(OUT, "mime101_summary.csv")).set_index("Model")
p = p.sort_values("train")

fig, axes = plt.subplots(1, 3, figsize=(160/25.4, 58/25.4))
## a: 训练 vs 验证（按训练排序）
y = np.arange(len(p))
axes[0].scatter(p["train"], y, s=2.5, c="#3C5488", label="train (TCGA, n=425)")
axes[0].scatter(p["val"], y, s=2.5, c="#E64B35", label="validation (GSE26193, n=107)")
axes[0].set_yticks([]); axes[0].set_xlabel("C-index", size=6); axes[0].tick_params(labelsize=5.4)
axes[0].set_title("a 117 models: train vs validation", size=6)
axes[0].legend(fontsize=5.4, frameon=False, loc="lower right")
axes[0].axvline(0.5, color="grey", linestyle=":", linewidth=0.5)
## b: gap 分布
axes[1].hist(p["gap"], bins=25, color="#B0B0B0")
axes[1].axvline(0.14, color="#E64B35", linestyle="--", linewidth=0.7)
axes[1].text(0.15, axes[1].get_ylim()[1] * 0.9, "mean gap 0.140", fontsize=5.4, color="#E64B35")
axes[1].set_xlabel("train - validation C-index", size=6); axes[1].tick_params(labelsize=5.4)
axes[1].set_title("b generalisation gap", size=6)
## c: 验证 top10
tp = p.sort_values("val", ascending=False).head(10).iloc[::-1]
axes[2].barh(range(len(tp)), tp["val"], color="#E64B35", height=0.6)
for i, (_, r) in enumerate(tp.iterrows()):
    axes[2].text(r["val"] + 0.005, i, f"{_:30s} train={r['train']:.2f}", va="center", ha="left", fontsize=5.4)
axes[2].set_yticks([]); axes[2].set_xlim(0.45, 0.95)
axes[2].set_xlabel("validation C-index", size=6); axes[2].tick_params(labelsize=5.4)
axes[2].set_title("c best validation (max 0.588)", size=6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "S7_mime101.png"), dpi=600)
plt.close(fig)
print("S7 saved", flush=True)
print("DONE 45g")
