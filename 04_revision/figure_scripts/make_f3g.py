# -*- coding: utf-8 -*-
"""Fig 3g: relabel from 'patient-level' to 'sample-level' and use the correct n.
Reproduces the original 46c computation (cell-level malignant flags grouped by sample).
"""
import numpy as np, pandas as pd, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"

prof = pd.read_csv(os.path.join(QRY, "F3a_chr_profile.csv.gz"), index_col=0)
cnv = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\percell_cnv_score_v2.csv",
                  usecols=["cell", "malignant_like", "group"])
cnv = cnv[cnv["group"].astype(str) == "Epithelial"]
mal = dict(zip(cnv["cell"], cnv["malignant_like"]))
cols_ok = [c for c in prof.columns if c in mal]
is_mal = np.array([mal[c] for c in cols_ok])

sc = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\module_scores\scores_GSE266577.csv",
                 usecols=["cell", "publication_sample_code_final"], low_memory=False)
smap = dict(zip(sc["cell"], sc["publication_sample_code_final"]))
sm = pd.read_csv(os.path.join(QRY, "F11_pt_subtype_mac.csv"))
sub_map = dict(zip(sm["publication_sample_code_final"], sm["pte_sub"]))

cc = pd.DataFrame({"cell": cols_ok, "mal": is_mal,
                   "pt": [smap.get(c) for c in cols_ok]})
pt_rate = cc.dropna(subset=["pt"]).groupby("pt")["mal"].mean()
pt_sub = pd.Series({p: sub_map.get(p) for p in pt_rate.index}).dropna()
n_high = int((pt_sub == "high").sum()); n_low = int((pt_sub == "low").sum())
print("samples plotted:", len(pt_rate), "patients:", len(set(p.split('_')[0] for p in pt_rate.index)),
      "high:", n_high, "low:", n_low, flush=True)
print("median high/low:", pt_rate[[p for p in pt_rate.index if sub_map.get(p) == 'high']].median(),
      pt_rate[[p for p in pt_rate.index if sub_map.get(p) == 'low']].median(), flush=True)

fig, ax = plt.subplots(figsize=(70/25.4, 52/25.4))
for i, sb in enumerate(["high", "low"]):
    v = pt_rate[[p for p in pt_rate.index if sub_map.get(p) == sb]]
    ax.bar(i, v.median(), color="#DC0000" if sb == "high" else "#3C5488", width=0.5)
    ax.scatter(np.full(len(v), i), v, s=10, c="#333333")
ax.set_xticks([0, 1])
ax.set_xticklabels(["high (n=%d)" % n_high, "low (n=%d)" % n_low], fontsize=5.4)
ax.set_ylabel("sample malignant fraction", size=6)
ax.tick_params(labelsize=5.4)
ax.set_title("sample malignant fraction by sample subtype", size=5.8)
ax.text(0.02, 0.02, "%d samples from 29 patients" % len(pt_rate),
        transform=ax.transAxes, fontsize=4.8, ha="left", va="bottom")
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F3g_pt_mal.png"), dpi=600)
plt.close(fig)
print("F3g saved", flush=True)
