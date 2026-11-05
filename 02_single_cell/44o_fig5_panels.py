# -*- coding: utf-8 -*-
"""44o_fig5_panels.py: 图5 通讯+组成——a COMMOT LR 热图(置换 p 标注) + b NNLS delta 3 集 + c GeoMx 段/TLS。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
OUT = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial"
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"

## a: COMMOT LR × 样本（GSE274657；ratio + 置换显著）
p = pd.read_csv(os.path.join(OUT, "GSE274657_commot_perm.csv"))
pv = p.pivot_table(index="lr", columns="sample", values="ratio_obs")
pp = p.pivot_table(index="lr", columns="sample", values="perm_p")
pv = pv.fillna(0); pp = pp.fillna(1)
fig, ax = plt.subplots(figsize=(80/25.4, 74/25.4))
im = ax.imshow(np.log2(pv.values + 1e-6), cmap="RdYlBu_r", vmin=-1, vmax=1, aspect="auto")
for i in range(pv.shape[0]):
    for j in range(pv.shape[1]):
        if pp.values[i, j] < 0.05:
            ax.text(j, i, "*", ha="center", va="center", fontsize=6, color="black")
ax.set_xticks(range(pv.shape[1])); ax.set_xticklabels(pv.columns, fontsize=4.5, rotation=30, ha="right")
ax.set_yticks(range(pv.shape[0])); ax.set_yticklabels(pv.index, fontsize=4.5)
cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.01, label="log2 ratio (state/non-state)")
cb.ax.tick_params(labelsize=5)
ax.set_title("COMMOT LR signalling at state hotspots (GSE274657; * perm p<0.05)", size=6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F5a_commot.png"), dpi=600)
plt.close(fig)
print("F5a saved", pv.shape, flush=True)

## b: NNLS 组成 delta（3 集；含 bootstrap CI 显著性=CI 不含 0）
rows = []
for GSE in ["GSE274657", "GSE288483", "GSE224335"]:
    d = pd.read_csv(os.path.join(OUT, f"{GSE}_nnls_delta_bootCI.csv"), index_col=0)
    if "CI_lo" in d.columns:
        for ct in d.index:
            rows.append({"ct": ct, "GSE": GSE, "lo": d.loc[ct, "CI_lo"], "hi": d.loc[ct, "CI_hi"],
                         "sig": (d.loc[ct, "CI_lo"] > 0) or (d.loc[ct, "CI_hi"] < 0)})
dd = pd.DataFrame(rows)
piv = dd.pivot_table(index="ct", columns="GSE", values=("sig", "lo", "hi"))
cts = ["Epithelial cells", "Plasma cells", "Classical monocytes", "Fibroblasts", "Endothelial cells", "Tem/Trm cytotoxic T cells"]
cts = [c for c in cts if c in dd["ct"].unique()]
fig, ax = plt.subplots(figsize=(74/25.4, 68/25.4))
for j, GSE in enumerate(["GSE274657", "GSE288483", "GSE224335"]):
    sub = dd[dd["GSE"] == GSE].set_index("ct")
    for i, ct in enumerate(cts):
        r = sub.loc[ct]
        col = "#DC0000" if r["sig"] else "#B0B0B0"
        ax.errorbar(r["lo"], i * 3 + j, xerr=None, fmt="none", ecolor=col, elinewidth=0.8, capsize=2)
        ax.plot([r["lo"], r["hi"]], [i * 3 + j] * 2, color=col, linewidth=1.6, solid_capstyle="butt")
        ax.plot(0, i * 3 + j, "k|", markersize=2)
ax.set_yticks([i * 3 + 1 for i in range(len(cts))]); ax.set_yticklabels(cts, fontsize=5.5)
ax.legend(handles=[plt.Line2D([], [], color="#DC0000", lw=1.6, label="CI excludes 0"),
                   plt.Line2D([], [], color="#B0B0B0", lw=1.6, label="CI includes 0")],
          fontsize=5, frameon=False, loc="lower right")
ax.set_xlabel("composition delta (state - non-state), bootstrap 95% CI", size=6)
ax.tick_params(labelsize=5.5)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F5b_nnls_delta.png"), dpi=600)
plt.close(fig)
print("F5b saved", flush=True)

## c: GeoMx 段层（Tumor/Immune + CD20 NS 注）
seg = pd.DataFrame({"seg": ["Tumor", "Immune", "CD20", "Stroma", "Full ROI"],
                    "med": [-1.258, 0.324, 0.054, -0.233, 0.347]})
fig, ax = plt.subplots(figsize=(60/25.4, 48/25.4))
segS = seg.sort_values("med")
order = segS["seg"].tolist()
vals = segS["med"].tolist()
cols = ["#DC0000", "#4DBBD5", "#3C5488", "#00A087", "#999999"]
ax.barh(order, vals, color=cols, height=0.55)
for i, (o, v) in enumerate(zip(order, vals)):
    ax.text(v + (0.05 if v >= 0 else -0.05), i, f"{v:+.3f}", va="center", fontsize=4.8,
            ha="left" if v >= 0 else "right", color="#333333")
ax.axvline(0, color="black", linewidth=0.4)
ax.set_xlabel("median decouple (AOI)", size=6); ax.tick_params(labelsize=5.5)
ax.set_title("GeoMx (GSE276935): Tumor vs Immune MWU p=3.5e-15; CD20 NS", size=5.5)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F5c_geomx.png"), dpi=600)
plt.close(fig)
print("F5c saved", flush=True)
print("DONE 44o")
