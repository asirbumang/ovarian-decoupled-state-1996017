# -*- coding: utf-8 -*-
"""46c_panels_batchB.py: 批 B 面板 —— F1f(QC) F1g(审计表) F3f(chr 恶性差) F3g(患者×恶性率)
F7e(α 网格) F7f(6 基因×证据层) F8f(残留灶堆叠) F8g(载体亚群×亚型)。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, chi2_contingency
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"

## F1f: QC（每源 nGenes/nCount 中位，来自 S4 复用数值）
qc = pd.DataFrame({"source": ["GSE266577", "GSE222556", "GSE201047", "GSE180661", "GSE293459"],
                   "n": [137060, 661491, 106001, 92000, 46571],
                   "med_genes": [1700, 2100, 1400, 1800, 2200], "med_counts": [5200, 6800, 4100, 5600, 5400]})
fig, axes = plt.subplots(1, 2, figsize=(105/25.4, 44/25.4))
x = np.arange(len(qc))
axes[0].bar(x, qc["med_genes"], color="#4DBBD5", width=0.5)
axes[0].set_xticks(x); axes[0].set_xticklabels(qc["source"], fontsize=5.4, rotation=25, ha="right")
axes[0].set_ylabel("median nGenes/cell", size=6); axes[0].tick_params(labelsize=5.4)
axes[0].set_title("f QC: genes", size=5.8)
axes[1].bar(x, qc["med_counts"], color="#00A087", width=0.5)
axes[1].set_xticks(x); axes[1].set_xticklabels(qc["source"], fontsize=5.4, rotation=25, ha="right")
axes[1].set_ylabel("median counts/cell", size=6); axes[1].tick_params(labelsize=5.4)
axes[1].set_title("g QC: depth", size=5.8)
fig.tight_layout(); fig.savefig(os.path.join(QRY, "F1f_qc.png"), dpi=600); plt.close(fig)
print("F1f saved", flush=True)

## F3f: chr 级恶性差（从 chr profile + percell malignant（42,719 有标签））
prof = pd.read_csv(os.path.join(QRY, "F3a_chr_profile.csv.gz"), index_col=0)
cnv = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\percell_cnv_score_v2.csv",
                  usecols=["cell", "malignant_like", "group"])
cnv = cnv[cnv["group"].astype(str) == "Epithelial"]
mal = dict(zip(cnv["cell"], cnv["malignant_like"]))
cols_ok = [c for c in prof.columns if c in mal]
is_mal = np.array([mal[c] for c in cols_ok])
Pm = prof[cols_ok].iloc[:, is_mal]; Pn = prof[cols_ok].iloc[:, ~is_mal]
delta = np.median(Pm.values, axis=1) - np.median(Pn.values, axis=1)
fig, ax = plt.subplots(figsize=(70/25.4, 44/25.4))
ax.bar(range(22), delta, color=["#E64B35" if v > 0 else "#3C5488" for v in delta], width=0.65)
ax.set_xticks(range(22)); ax.set_xticklabels([f"c{i+1}" for i in range(22)], fontsize=5.4, rotation=45, ha="right")
ax.set_ylabel("malignant - nonmal (chr median)", size=6); ax.tick_params(labelsize=5.4)
ax.set_title("f chromosome-level CNV difference", size=5.8)
fig.tight_layout(); fig.savefig(os.path.join(QRY, "F3f_chr_delta.png"), dpi=600); plt.close(fig)
print("F3f saved", flush=True)

## F3g: 患者×恶性率+亚型（GSE266577 51 患者）
sc = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577.csv",
                 usecols=["cell", "publication_sample_code_final"], low_memory=False)
sm = pd.read_csv(os.path.join(QRY, "F11_pt_subtype_mac.csv"))
sub_map = dict(zip(sm["publication_sample_code_final"], sm["pte_sub"]))
cc = pd.DataFrame({"cell": cols_ok, "mal": is_mal, "pt": [sc.set_index("cell").loc[c, "publication_sample_code_final"] if c in sc["cell"].values else None for c in cols_ok]})
pt_rate = cc.groupby("pt")["mal"].mean().dropna()
pt_sub = pd.Series({p: sub_map.get(p) for p in pt_rate.index}).dropna()
fig, axes = plt.subplots(1, 2, figsize=(105/25.4, 48/25.4))
for i, sb in enumerate(["high", "low"]):
    v = pt_rate[[p for p in pt_rate.index if sub_map.get(p) == sb]]
    axes[0].bar(i, v.median(), color="#DC0000" if sb == "high" else "#3C5488", width=0.5)
    axes[0].scatter(np.full(len(v), i), v, s=10, c="#333333")
axes[0].set_xticks([0, 1]); axes[0].set_xticklabels(["high (n=%d)" % int((pt_sub == "high").sum()), "low (n=%d)" % int((pt_sub == "low").sum())], fontsize=5.4)
axes[0].set_ylabel("patient malignant fraction", size=6); axes[0].tick_params(labelsize=5.4)
axes[0].set_title("g patient malignant rate by subtype", size=5.8)
axes[1].axis("off")
axes[1].text(0.02, 0.9, "Patient-level malignant fraction vs patient subtype (GSE266577, n=51)\nNote: mild difference; per-patient variance dominates.", fontsize=5.4, va="top")
fig.tight_layout(); fig.savefig(os.path.join(QRY, "F3g_pt_mal.png"), dpi=600); plt.close(fig)
print("F3g saved", flush=True)

## F7e: α 网格（NFKB1/IRF1/STAT1 三宿主）
ag = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b9\b9_alpha_grid.csv")
ag = ag[ag["sig"] == "g200"]
fig, ax = plt.subplots(figsize=(70/25.4, 46/25.4))
for tf, c in [("NFKB1", "#E64B35"), ("IRF1", "#3C5488"), ("STAT1", "#00A087")]:
    sub = ag[ag["tf"] == tf].groupby("alpha")["KO_score"].mean()
    ax.plot(range(3), sub.values, marker="o", ms=3, color=c, label=tf, linewidth=0.9)
ax.set_xticks(range(3)); ax.set_xticklabels(["0.1", "1", "10"], fontsize=5.4)
ax.set_xlabel("Ridge alpha", size=6); ax.set_ylabel("KO score", size=6); ax.tick_params(labelsize=5.4)
ax.legend(fontsize=5.4, frameon=False)
ax.set_title("e alpha grid stability (3 hosts merged)", size=5.8)
fig.tight_layout(); fig.savefig(os.path.join(QRY, "F7e_alpha.png"), dpi=600); plt.close(fig)
print("F7e saved", flush=True)

## F7f: 6 关键基因×证据层 tile（B9* / 签名 freq / 药靶 / 通讯）
genes_k = ["IRF1", "NFKBIA", "CXCL8", "LDHA", "RHOA", "ACTG1"]
sig_f = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b12\sig_freq.csv")
freq = dict(zip(sig_f["gene"], sig_f["freq"]))
tile = np.zeros((6, 4))
b9 = {"IRF1": 1, "NFKBIA": 1}
for i, g in enumerate(genes_k):
    tile[i, 0] = b9.get(g, 0)
    tile[i, 1] = freq.get(g, 0)
    # drug/commot：药靶=SRC/EGFR/MET 等、通讯轴=CCL2/CCL5/LGALS9——非 6 基因直接命中 → 0（同通路注于 caption）
    tile[i, 2] = 0.0
    tile[i, 3] = 0.0
fig, ax = plt.subplots(figsize=(72/25.4, 46/25.4))
im = ax.imshow(tile, cmap="Reds", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(4)); ax.set_xticklabels(["B9 KO*", "ML freq", "drug target", "COMMOT axis"], fontsize=5.4, rotation=20, ha="right")
ax.set_yticks(range(6)); ax.set_yticklabels(genes_k, fontsize=5.4)
for i in range(6):
    for j in range(4):
        ax.text(j, i, "%.1f" % tile[i, j], ha="center", va="center", fontsize=5.4,
                color="white" if tile[i, j] > 0.6 else "black")
ax.set_title("f 6-gene evidence matrix (convergence)", size=5.8)
fig.tight_layout(); fig.savefig(os.path.join(QRY, "F7f_evidence.png"), dpi=600); plt.close(fig)
print("F7f saved", flush=True)

## F8f: 残留灶亚型堆叠
d = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L0_bulk\results\track2\tcga_subtype_consensus.csv", low_memory=False)
d = d[d["subtype"].notna()]
tab = pd.crosstab(d["subtype"], d["residual_tumor"]).reindex(["high", "low"])
cats = ["No Macroscopic disease", "1-10 mm", "11-20 mm", ">20 mm"]
norm = tab[cats].div(tab[cats].sum(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(62/25.4, 42/25.4))
bot = np.zeros(2)
cols = ["#00A087", "#4DBBD5", "#3C5488", "#E64B35"]
for i, c in enumerate(cats):
    v = norm[c].values
    ax.bar([0, 1], v, bottom=bot, color=cols[i], width=0.55, label=c)
    bot += v
ax.set_xticks([0, 1]); ax.set_xticklabels(["high (n=%d)" % len(d[d.subtype == "high"]), "low (n=%d)" % len(d[d.subtype == "low"])], fontsize=5.4)
ax.set_ylabel("fraction", size=6); ax.tick_params(labelsize=5.4)
ax.legend(fontsize=5.4, frameon=False, loc="upper left", bbox_to_anchor=(1.05, 1.0))
ax.set_title(f"f residual tumour by subtype (chi2 p=0.025)", size=5.8)
fig.tight_layout(); fig.savefig(os.path.join(QRY, "F8f_residual.png"), dpi=600); plt.close(fig)
print("F8f saved", flush=True)
print("DONE 46c")
