# -*- coding: utf-8 -*-
"""44c_fig2_panels.py: 图2 a(UMAP 状态高亮) + c(四象限) + d(稀有度跨源)。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
z = lambda x: (x - x.mean()) / (x.std() + 1e-9)
QRY = "D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs"

## --- a: UMAP + state 高亮（GSE266577 数据集内全细胞定义一致性；§17 相同口径）---
umo = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\umap_input\umap_full_266577.csv")
scf = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577.csv",
                  usecols=["cell", "M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"],
                  low_memory=False)
for m in ["M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"]:
    scf[m] = pd.to_numeric(scf[m], errors="coerce").astype(float)
z1 = z(scf["M1_immune_activation"]); z2 = z(scf["M2_metabolic_GSH"]); z3 = z(scf["M3_mito_cytoskeleton"])
scf["state_high"] = (z1 >= np.quantile(z1, .75)) & (z2 <= np.quantile(z2, .25)) & (z3 <= np.quantile(z3, .25))
umo = umo.merge(scf[["cell", "state_high"]], on="cell", how="left")
umo["state_high"] = umo["state_high"].fillna(False)
fig, ax = plt.subplots(figsize=(72/25.4, 66/25.4))
ax.scatter(umo["umap1"], umo["umap2"], s=0.06, c="#AAAAAA", alpha=0.22, linewidths=0)
## 12 型着色（图1b 同色；证明跨类型）
sc_l = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577.csv",
                   usecols=["cell", "cell_type"], low_memory=False)
LIN = {"Epithelial cells":"Epi","Fibroblasts":"Fibro","Mesothelial cells":"Meso","Endothelial cells":"Endo",
 "Macrophages":"MonoMac","Classical monocytes":"MonoMac","Migratory DCs":"DC","DC1":"DC","DC2":"DC","pDC":"DC",
 "Mast cells":"Mast","Tem/Trm cytotoxic T cells":"T","Tcm/Naive helper T cells":"T","Regulatory T cells":"T",
 "Type 17 helper T cells":"T","NK cells":"NK","CD16- NK cells":"NK","CD16+ NK cells":"NK","ILC":"NK",
 "Memory B cells":"B","Naive B cells":"B","Plasma cells":"Plasma"}
sc_l["lin"] = sc_l["cell_type"].map(LIN)
umo = umo.merge(sc_l[["cell", "lin"]], on="cell", how="left")
NEWS = {"Epi":"#E64B35","Fibro":"#4DBBD5","Meso":"#00A087","Endo":"#3C5488","MonoMac":"#F39B7F",
        "DC":"#8491B4","Mast":"#91D1C2","T":"#DC0000","NK":"#7E6148","B":"#B09C85","Plasma":"#A73030"}
ums = umo[umo["lin"].notna()]
for lin, col in NEWS.items():
    s = ums[ums["lin"] == lin]
    if len(s):
        ax.scatter(s["umap1"], s["umap2"], s=0.18, c=col, alpha=0.5, linewidths=0)
st = umo[umo["state_high"]]
ax.scatter(st["umap1"], st["umap2"], s=0.35, c="#000000", alpha=0.9, linewidths=0,
           label="decoupled state (n=%d; 1.28%% of cells)" % len(st))
ax.set_xticks([]); ax.set_yticks([])
for s_ in ax.spines.values():
    s_.set_visible(False)
handles, labels = [], []
for lin, col in NEWS.items():
    handles.append(plt.Line2D([], [], marker="o", color=col, linestyle="None", markersize=2.2, alpha=0.7))
    labels.append(lin)
handles.append(plt.Line2D([], [], marker="o", color="black", linestyle="None", markersize=3.0))
labels.append("decoupled state")
leg = ax.legend(handles, labels, loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=5.0, frameon=False, ncol=1)
for h in leg.legend_handles:
    h.set_markersize(2.5)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F2a_state_umap.png"), dpi=600)
plt.close(fig)
print("F2a saved; state cells (all):", int(st["state_high"].sum()), flush=True)

## --- c: 解耦四象限（全细胞；M1z vs mean(M2z,M3z)）---
x = z(scf["M2_metabolic_GSH"] + scf["M3_mito_cytoskeleton"])
y = z1
stt = scf["state_high"].values.astype(bool)
fig, ax = plt.subplots(figsize=(62/25.4, 62/25.4))
ax.scatter(x[~stt], y[~stt], s=0.4, c="#BBBBBB", alpha=0.35, linewidths=0)
ax.scatter(x[stt], y[stt], s=1.4, c="#DC0000", alpha=0.8, linewidths=0)
ax.axvline(0, linestyle=":", linewidth=0.35, color="grey")
ax.axhline(0, linestyle=":", linewidth=0.35, color="grey")
ax.set_xlabel("M2+M3 (gsh/metabolic · mito-cyto), z")
ax.set_ylabel("M1 (immune activation), z")
ax.tick_params(labelsize=5.5)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F2c_quadrant.png"), dpi=600)
plt.close(fig)
print("F2c saved; state frac:", round(stt.mean() * 100, 2), "% of all cells", flush=True)

## --- d: 稀有度跨源（各 scRNA 源内 state%）---
rows = []
for f, tag in [("scores_GSE266577.csv", "GSE266577"), ("scores_GSE222556.csv", "GSE222556"),
               ("scores_GSE201047.csv", "GSE201047"), ("scores_GSE180661_subset.csv", "GSE180661"),
               ("scores_GSE293459.csv", "GSE293459")]:
    p = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores" + "\\" + f
    if not os.path.exists(p):
        print("skip", tag, p); continue
    d = pd.read_csv(p, low_memory=False)
    for m in ["M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"]:
        d[m] = pd.to_numeric(d[m], errors="coerce").astype(float)
    z1 = z(d["M1_immune_activation"]); z2 = z(d["M2_metabolic_GSH"]); z3 = z(d["M3_mito_cytoskeleton"])
    stt = (z1 >= np.quantile(z1, .75)) & (z2 <= np.quantile(z2, .25)) & (z3 <= np.quantile(z3, .25))
    rows.append({"source": tag, "n": len(d), "state_n": int(stt.sum()), "pct": round(100 * stt.mean(), 3)})
rr = pd.DataFrame(rows)
rr.to_csv(os.path.join(QRY, "F2d_rarity.csv"), index=False)
print(rr.to_string(index=False), flush=True)
fig, ax = plt.subplots(figsize=(72/25.4, 52/25.4))
P = rr["pct"] / 100.0
N = rr["n"].astype(float)
se = np.sqrt(P * (1 - P) / N) * 100
ax.scatter(rr["pct"], range(len(rr)), s=14, c="#3C5488")
ax.errorbar(rr["pct"], range(len(rr)), xerr=1.96 * se, fmt="none", ecolor="#3C5488", elinewidth=0.5, capsize=1.5)
for i, r in rr.iterrows():
    ax.text(r["pct"] + 0.38, i, f"{r['source']} (n={r['n']:,}) | {r['pct']:.2f}% [±{1.96*se[i]:.2f}]", va="center", fontsize=5)
ax.set_yticks([]); ax.set_xlabel("state frequency (% of cells)", size=6)
ax.set_xlim(0, max(rr["pct"]) * 1.9)
ax.tick_params(labelsize=5.5)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F2d_rarity.png"), dpi=600)
plt.close(fig)
print("F2d saved")
print("DONE 44c")
