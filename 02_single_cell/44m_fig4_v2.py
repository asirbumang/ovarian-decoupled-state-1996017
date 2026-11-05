# -*- coding: utf-8 -*-
"""44m_fig4_v2.py: 图4 v2 —— a 黑底 triptych(2×3) + c 横排短标签统计（重制）。"""
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

def dark(ax):
    ax.set_facecolor("black")
    for s_ in ax.spines.values():
        s_.set_visible(False)

def draw(ax, sub, mode):
    ax.scatter(sub["array_col"], -sub["array_row"], s=0.35, c="#444444", linewidths=0)
    if mode == "domain":
        ax.scatter(sub["array_col"], -sub["array_row"], s=1.2, c=sub["domain"], cmap=plt.cm.tab10, vmin=0, vmax=9, linewidths=0)
    else:
        st = sub[sub["state_high"]]
        ax.scatter(st["array_col"], -st["array_row"], s=2.0, c="#DC0000", linewidths=0, alpha=0.95)
    ax.set_aspect("equal"); ax.axis("off")
    dark(ax)
    x1 = sub["array_col"].min() + 0.04 * (sub["array_col"].max() - sub["array_col"].min())
    x2 = x1 + 10
    y1 = sub["array_row"].min() + 0.90 * (sub["array_row"].max() - sub["array_row"].min())
    ax.plot([x1, x2], [-y1, -y1], color="white", linewidth=0.8)
    ax.text((x1 + x2) / 2, -y1 - 1.5, "1 mm", color="white", fontsize=4.2, ha="center")

fig, axes = plt.subplots(2, 3, figsize=(120/25.4, 74/25.4))
from matplotlib.patches import Patch
leg_handles = [Patch(facecolor=plt.cm.tab10(i / 10), edgecolor="none", label=f"D{i}") for i in range(8)]
for j, (GSE, sm) in enumerate([("GSE274657", "Patient_16"), ("GSE288483", "Pt2-4")]):
    d = pd.read_csv(os.path.join(OUT, f"{GSE}_domains.csv"))
    sub = d[(d["sample"] == sm)].dropna(subset=["array_row", "array_col"])
    draw(axes[j, 0], sub, "domain")
    draw(axes[j, 1], sub, "state")
    for ax, t in [(axes[j, 0], f"{sm} | domains"), (axes[j, 1], f"{sm} | state")]:
        ax.text(0.5, 0.97, t, transform=ax.transAxes, ha="center", va="top", color="white", fontsize=5.2)
    st = sub[sub["state_high"]]
    if len(st) > 10:
        x0, x1 = st["array_col"].quantile(0.1), st["array_col"].quantile(0.9)
        y0, y1 = st["array_row"].quantile(0.1), st["array_row"].quantile(0.9)
        axz = axes[j, 2]
        axz.scatter(sub["array_col"], -sub["array_row"], s=0.25, c="#666666", linewidths=0)
        axz.scatter(st["array_col"], -st["array_row"], s=1.4, c="#DC0000", linewidths=0)
        axz.set_xlim(x0 - (x1 - x0) * 0.5, x1 + (x1 - x0) * 0.5)
        axz.set_ylim(-(y1 + (y1 - y0) * 0.5), -(y0 - (y1 - y0) * 0.5))
        axz.set_aspect("equal"); axz.axis("off"); dark(axz)
        axz.text(0.5, 0.97, "zoom (state)", transform=axz.transAxes, ha="center", va="top", color="white", fontsize=5.2)
        axes[j, 0].add_patch(Rectangle((x0, -y1), x1 - x0, y1 - y0, fill=False, edgecolor="white", linewidth=0.7))
fig.legend(handles=leg_handles, loc="upper right", ncol=4, fontsize=4.5, frameon=False, title="NMF domains",
           title_fontsize=5)
fig.savefig(os.path.join(QRY, "F4a_triptych_v2.png"), dpi=600)
plt.close(fig)
print("F4a v2 saved", flush=True)

## c: 横排点图（Moran + OR；样本缩写）
m = pd.read_csv(os.path.join(OUT, "spatial_moran_summary.csv"))
n = pd.read_csv(os.path.join(OUT, "spatial_neighborhood_summary.csv"))
n = n[n["fisher_p"].notna()]
short = {"GSE274657": "V", "GSE288483": "G", "GSE224335": "X"}
def mk(s):
    g = str(s["GSE"])[:6]
    sm = str(s["sample"]).replace("Patient_", "P").replace("Pt", "P")
    return short.get(str(s["GSE"]), g[:6]) + "-" + sm
m["ab"] = m.apply(mk, axis=1)
n["ab"] = n.apply(mk, axis=1)
fig, axes = plt.subplots(1, 2, figsize=(115/25.4, 58/25.4), sharey=True)
ypos = np.arange(len(m))
axes[0].scatter(m["Moran_I"], ypos, s=9, c="#3C5488")
axes[0].set_yticks(ypos); axes[0].set_yticklabels(m["ab"], fontsize=4.5)
axes[0].set_xlabel("Moran's I (decouple)", size=6); axes[0].tick_params(labelsize=5.5)
axes[0].axvline(m["Moran_I"].median(), color="#E64B35", linestyle="--", linewidth=0.5)
axes[0].text(m["Moran_I"].median() + 0.015, len(m) * 0.97, "median %.3f" % m["Moran_I"].median(), fontsize=5, color="#E64B35")
axes[0].set_title("14/14 samples p<0.05", size=6)
ny = [i for i in range(len(n))]
axes[1].scatter(n["neigh_OR"], ny, s=9, c="#3C5488")
for i, (_, r) in enumerate(n.iterrows()):
    if r["fisher_p"] < 0.05:
        axes[1].scatter(r["neigh_OR"], i, s=12, c="#DC0000")
axes[1].axvline(1, color="grey", linestyle=":", linewidth=0.5)
axes[1].set_xlabel("neighbourhood OR", size=6); axes[1].tick_params(labelsize=5.5)
axes[1].set_title("red: p<0.05 (7/12); all samples n_spots>=898", size=6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F4cd_stats_v2.png"), dpi=600)
plt.close(fig)
print("F4cd v2 saved", flush=True)
print("DONE 44m")
