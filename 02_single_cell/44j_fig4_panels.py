# -*- coding: utf-8 -*-
"""44j_fig4_panels.py: 图4 空间三层——a 代表样本 triptych(黑底+比例尺+红框zoom) + b 样本×域 state 率热图。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
OUT = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial"
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"

## a: 代表样本（P16 / Pt2-4）三联动：域 + state + zoom
def draw_panel(ax, GSE, sm, mode, um=100.0):
    d = pd.read_csv(os.path.join(OUT, f"{GSE}_domains.csv"))
    sub = d[(d["sample"] == sm)].dropna(subset=["array_row", "array_col"])
    ax.scatter(sub["array_col"], -sub["array_row"], s=0.4, c="#333333", linewidths=0)
    if mode == "domain":
        cmap = plt.cm.tab10
        ax.scatter(sub["array_col"], -sub["array_row"], s=1.2, c=sub["domain"], cmap=cmap, vmin=0, vmax=9, linewidths=0)
    else:
        st = sub[sub["state_high"]]
        ax.scatter(st["array_col"], -st["array_row"], s=2.0, c="#DC0000", linewidths=0, alpha=0.9)
    ax.set_aspect("equal"); ax.axis("off")
    # 1mm 比例尺（Visium spot pitch = 100µm；1mm = 10 spots）
    x1 = sub["array_col"].min() + 0.04 * (sub["array_col"].max() - sub["array_col"].min())
    x2 = x1 + 10
    y1 = sub["array_row"].min() + 0.92 * (sub["array_row"].max() - sub["array_row"].min())
    ax.plot([x1, x2], [-y1, -y1], color="white", linewidth=0.8)
    ax.text((x1 + x2) / 2, -y1 - 1.5, "1 mm", color="white", fontsize=4.5, ha="center")

fig, axes = plt.subplots(2, 4, figsize=(170/25.4, 62/25.4))
for j, (GSE, sm) in enumerate([("GSE274657", "Patient_16"), ("GSE288483", "Pt2-4")]):
    d = pd.read_csv(os.path.join(OUT, f"{GSE}_domains.csv"))
    sub = d[(d["sample"] == sm)].dropna(subset=["array_row", "array_col"])
    # 奇数列=域 / state；末列=zoom(取 state 高密度区)
    draw_panel(axes[j, 0], GSE, sm, "domain")
    draw_panel(axes[j, 1], GSE, sm, "state")
    axes[j, 0].set_title(f"{sm}: NMF domains", size=5.5)
    axes[j, 1].set_title(f"{sm}: state (orange-red)", size=5.5)
    # zoom：以 state 细胞包围盒（扩大 8 倍）绘制红框+子图
    st = sub[sub["state_high"]]
    if len(st) > 10:
        x0, x1 = st["array_col"].quantile(0.1), st["array_col"].quantile(0.9)
        y0, y1 = st["array_row"].quantile(0.1), st["array_row"].quantile(0.9)
        axp = axes[j, 2]
        axp.scatter(sub["array_col"], -sub["array_row"], s=0.3, c="#555555", linewidths=0)
        axp.scatter(st["array_col"], -st["array_row"], s=1.4, c="#DC0000", linewidths=0)
        axp.set_xlim(x0 - (x1 - x0) * 0.3, x1 + (x1 - x0) * 0.3)
        axp.set_ylim(-(y1 + (y1 - y0) * 0.3), -(y0 - (y1 - y0) * 0.3))
        axp.set_aspect("equal"); axp.axis("off")
        axes[j, 2].set_title("zoom (state)", size=5.5)
        axes[j, 0].add_patch(Rectangle((x0, -y1), x1 - x0, y1 - y0, fill=False, edgecolor="white", linewidth=0.7))
    axes[j, 3].axis("off")
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F4a_triptych.png"), dpi=600)
plt.close(fig)
print("F4a saved", flush=True)

## b: 14 样本×8 域 state 率热图
rows = []
for GSE in ["GSE274657", "GSE288483", "GSE224335"]:
    d = pd.read_csv(os.path.join(OUT, f"{GSE}_domains_state.csv"))
    for _, r in d.iterrows():
        rows.append({"GSE": GSE, "sample": r["sample"], "domain": int(r["domain"]),
                     "rate": r["mean"] * 100, "n": int(r["size"])})
df = pd.DataFrame(rows)
piv = df.pivot_table(index=["GSE", "sample"], columns="domain", values="rate")
piv = piv.fillna(0)
fig, ax = plt.subplots(figsize=(80/25.4, 74/25.4))
im = ax.imshow(piv.values, cmap="YlOrRd", vmin=0, vmax=max(piv.values.max(), 6), aspect="auto")
ax.set_xticks(range(8)); ax.set_xticklabels([f"D{i}" for i in range(8)], fontsize=5)
ax.set_yticks(range(len(piv))); ax.set_yticklabels([f"{i} ({j})" for i, j in piv.index], fontsize=4.5)
cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.01, label="state rate (%)")
cb.ax.tick_params(labelsize=5)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F4b_domain_state.png"), dpi=600)
plt.close(fig)
print("F4b saved", piv.shape, flush=True)
print("DONE 44j")
