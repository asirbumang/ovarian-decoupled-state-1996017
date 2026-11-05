# -*- coding: utf-8 -*-
"""44s_fig7_panels.py: 图7 扰动+药筛——a B9 热图(置换 p*) + b 组合 KO + c 药筛 top15。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
CNV = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b9"
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"

## a: B9 热图（TF×宿主 KO 分数 + 置换 p *）
s = pd.read_csv(os.path.join(CNV, "b9_summary_perm_combo.csv"))
s = s[~s["tf"].str.contains("_KO", na=False)]
tfs = ["NFKB1", "IRF1", "STAT1", "RELA", "IRF7", "SP1", "HIF1A", "JUN", "IRF3", "STAT3"]
hosts = ["GSE266577_epithelial", "GSE266577_malignant", "GSE293459_all"]
piv = s.pivot_table(index="tf", columns="host", values="KO").reindex(tfs)
pp = s.pivot_table(index="tf", columns="host", values="p_perm").reindex(tfs)
pp = pp.fillna(1.0)
fig, ax = plt.subplots(figsize=(70/25.4, 64/25.4))
im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-0.16, vmax=0.04, aspect="auto")
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        if pp.values[i, j] < 0.05:
            ax.text(j, i, "*", ha="center", va="center", fontsize=6, color="black")
ax.set_xticks(range(len(hosts))); ax.set_xticklabels([h.replace("GSE266577_", "266577 ").replace("GSE293459_all", "293459") for h in hosts], fontsize=4.8, rotation=25, ha="right")
ax.set_yticks(range(len(tfs))); ax.set_yticklabels(tfs, fontsize=5)
cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.01, label="perturbation score (KO; <0 revert)")
cb.ax.tick_params(labelsize=5)
ax.set_title("Heuristic TF perturbation (CollecTRI + Ridge; * perm p<0.05)", size=5.5)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F7a_b9.png"), dpi=600)
plt.close(fig)
print("F7a saved", flush=True)

## b: 组合 KO（4 组合×3 宿主；单 TF max 参考线）
c = pd.read_csv(os.path.join(CNV, "b9_summary_perm_combo.csv"))
comb = c[c["tf"].str.contains("_KO", na=False)]
comb["combo"] = comb["tf"].str.replace("_KO", "", regex=False)
sing = c[~c["tf"].str.contains("_KO", na=False)]
fig, ax = plt.subplots(figsize=(68/25.4, 54/25.4))
for j, h in enumerate(hosts):
    sub = comb[comb["host"] == h]
    ys = np.arange(len(sub)) * 3 + j
    ax.scatter(sub["KO"], ys, s=14, c="#3C5488")
    ax.plot(sub["KO"], ys, color="#B0B0B0", linewidth=0.5)
singmax = sing.groupby("host")["KO"].min().values
for j, h in enumerate(hosts):
    sub = comb[comb["host"] == h]
    ax.axvline(singmax[j], color="#E64B35", linestyle="--", linewidth=0.5)
ax.set_yticks(np.arange(4) * 3 + 1); ax.set_yticklabels(["NFKB1+IRF1", "NFKB1+STAT1", "IRF1+STAT1", "NFKB1+IRF1+STAT1"], fontsize=5)
ax.set_xlabel("combo KO perturbation score (dashed = best single TF)", size=6)
ax.tick_params(labelsize=5.5)
ax.set_title("Combination virtual KO (3 hosts)", size=6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F7b_combo.png"), dpi=600)
plt.close(fig)
print("F7b saved", flush=True)

## c: 药筛 top15（rho + q）
de = pd.read_csv(r"D:\子宫内膜异位症2\data_omics\drug_screen_all_annot.csv")
de = de.iloc[:15].copy().iloc[::-1]
fig, ax = plt.subplots(figsize=(76/25.4, 62/25.4))
cols = ["#E64B35" if v < 0 else "#3C5488" for v in de["rho"]]
ax.barh(range(len(de)), de["rho"], color=cols, height=0.62)
labels = []
for _, r in de.iterrows():
    nm = str(r["name"]).split(" (")[0][:22]
    tg = str(r["target"])[:14] if pd.notna(r["target"]) else ""
    labels.append(nm + ("\n[" + tg + "]" if tg else ""))
ax.set_yticks(range(len(de))); ax.set_yticklabels(labels, fontsize=3.9)
for i, (_, r) in enumerate(de.iterrows()):
    x = r["rho"] - 0.006 if r["rho"] < 0 else r["rho"] + 0.006
    ax.text(x, i, "q<0.001" if r["q"] < 1e-3 else f"q={r['q']:.3f}",
            va="center", ha="right" if r["rho"] < 0 else "left", size=3.8)
ax.axvline(0, color="black", linewidth=0.35)
ax.set_xlim(-0.42, 0.34)
ax.set_xlabel("Spearman rho (signature score vs AUC)", size=6)
ax.set_title("CTRPv2 sensitivity correlates", size=5.8)
ax.tick_params(axis="x", labelsize=5)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F7c_drug.png"), dpi=600)
plt.close(fig)
print("F7c saved", flush=True)
print("DONE 44s")
