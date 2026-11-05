# -*- coding: utf-8 -*-
"""44w_fig2_expand.py: 图2 扩充面板 d(宿主/谱系 OR) e(三模块小提琴) f(组成堆叠)。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact, mannwhitneyu
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"

## d: 谱系 OR（GSE266577 全细胞 state × 12 谱系；Fisher + CI 近似）
sc = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577.csv",
                 usecols=["cell_type", "M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"], low_memory=False)
for m in ["M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"]:
    sc[m] = pd.to_numeric(sc[m], errors="coerce").astype(float)
z = lambda x: (x - x.mean()) / (x.std() + 1e-9)
z1, z2, z3 = z(sc["M1_immune_activation"]), z(sc["M2_metabolic_GSH"]), z(sc["M3_mito_cytoskeleton"])
st = (z1 >= np.quantile(z1, .75)) & (z2 <= np.quantile(z2, .25)) & (z3 <= np.quantile(z3, .25))
LIN = {"Epithelial cells":"Epithelial","Fibroblasts":"Fibroblast","Mesothelial cells":"Mesothelial","Endothelial cells":"Endothelial",
 "Macrophages":"Myeloid","Classical monocytes":"Myeloid","Migratory DCs":"DC","DC1":"DC","DC2":"DC","pDC":"DC",
 "Mast cells":"Mast","Tem/Trm cytotoxic T cells":"T","Tcm/Naive helper T cells":"T","Regulatory T cells":"T",
 "Type 17 helper T cells":"T","NK cells":"NK","CD16- NK cells":"NK","CD16+ NK cells":"NK","ILC":"NK",
 "Memory B cells":"B","Naive B cells":"B","Plasma cells":"Plasma"}
sc["lin"] = sc["cell_type"].map(LIN)
rows = []
for lin, g in sc.groupby("lin"):
    a = int((st[g.index]).sum()); b = int(len(g) - a)
    c_ = int(st.sum() - a); d_ = int((~st).sum() - b)
    orr, p = fisher_exact([[a, b], [c_, d_]])
    se = np.sqrt(1 / max(a, 1) + 1 / max(b, 1) + 1 / max(c_, 1) + 1 / max(d_, 1))
    rows.append({"lin": lin, "OR": orr, "lo": np.exp(np.log(orr) - 1.96 * se), "hi": np.exp(np.log(orr) + 1.96 * se), "p": p, "n": len(g)})
dd = pd.DataFrame(rows).sort_values("OR")
fig, ax = plt.subplots(figsize=(60/25.4, 58/25.4))
y = np.arange(len(dd))
ax.errorbar(dd["OR"], y, xerr=[dd["OR"] - dd["lo"], dd["hi"] - dd["OR"]], fmt="o", ms=3, color="#3C5488",
            ecolor="#3C5488", elinewidth=0.6, capsize=1.5)
ax.axvline(1, color="grey", linestyle=":", linewidth=0.5)
ax.set_yticks(y); ax.set_yticklabels(dd["lin"], fontsize=4.8)
ax.set_xscale("log"); ax.set_xlim(0.05, 40)
ax.set_xlabel("OR (state vs rest per lineage)", size=6); ax.tick_params(labelsize=5)
ax.set_title("State enrichment by lineage (GSE266577)", size=5.8)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F2d_or.png"), dpi=600)
plt.close(fig)
print("F2d OR saved", flush=True)
print(dd[["lin", "OR", "p"]].round(3).to_string(index=False), flush=True)

## e: 三模块小提琴（state vs 其余）
fig, axes = plt.subplots(1, 3, figsize=(90/25.4, 46/25.4))
for ax, m, t in zip(axes, ["M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton"],
                    ["M1 immune", "M2 metab/GSH", "M3 mito-cyto"]):
    a = sc.loc[st, m].values; b = sc.loc[~st, m].values
    vp = ax.violinplot([b, a], positions=[1, 2], showmedians=True, widths=0.7)
    for pc in vp["bodies"]:
        pc.set_alpha(0.5)
    p = mannwhitneyu(a, b).pvalue
    pstr = "<1e-300" if p == 0 else f"{p:.1e}"
    ax.set_xticks([1, 2]); ax.set_xticklabels(["rest", "state"], fontsize=5)
    ax.set_title(f"{t} (MWU p={pstr})", size=5.5); ax.tick_params(labelsize=5)
    ax.set_xlim(0.4, 2.6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F2e_modules.png"), dpi=600)
plt.close(fig)
print("F2e saved", flush=True)

## f: 组成堆叠（state vs 全体）
comp_order = ["T", "Myeloid", "NK", "B", "Plasma", "DC", "Mast", "Endothelial", "Fibroblast", "Mesothelial", "Epithelial"]
cnt_all = sc["lin"].value_counts(normalize=True).reindex(comp_order).fillna(0)
cnt_st = sc.loc[st, "lin"].value_counts(normalize=True).reindex(comp_order).fillna(0)
fig, ax = plt.subplots(figsize=(64/25.4, 46/25.4))
cols = plt.cm.tab20(np.linspace(0, 1, len(comp_order)))
bottom = np.zeros(2)
for i, c in enumerate(comp_order):
    vals = [cnt_all[c], cnt_st[c]]
    ax.bar([0, 1], vals, bottom=bottom, color=cols[i], width=0.6, label=c)
    bottom += vals
ax.set_xticks([0, 1]); ax.set_xticklabels(["all cells", "state"], fontsize=5.5)
ax.set_ylabel("fraction", size=6); ax.tick_params(labelsize=5)
ax.legend(fontsize=4.5, frameon=False, loc="upper right", bbox_to_anchor=(1.55, 1.0))
ax.set_title("Composition: state vs all cells", size=5.8)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F2f_comp.png"), dpi=600)
plt.close(fig)
print("F2f saved", flush=True)
print("DONE 44w")
