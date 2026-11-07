# -*- coding: utf-8 -*-
"""Fig 3f: chromosome-level CNV difference, without the stray 'f ' title prefix."""
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
Pm = prof[cols_ok].iloc[:, is_mal]; Pn = prof[cols_ok].iloc[:, ~is_mal]
delta = np.median(Pm.values, axis=1) - np.median(Pn.values, axis=1)

fig, ax = plt.subplots(figsize=(70/25.4, 44/25.4))
ax.bar(range(22), delta, color=["#E64B35" if v > 0 else "#3C5488" for v in delta], width=0.65)
ax.set_xticks(range(22)); ax.set_xticklabels([f"c{i+1}" for i in range(22)], fontsize=5.4, rotation=45, ha="right")
ax.set_ylabel("malignant - nonmal (chr median)", size=6); ax.tick_params(labelsize=5.4)
ax.set_title("chromosome-level CNV difference", size=5.8)
fig.tight_layout(); fig.savefig(os.path.join(QRY, "F3f_chr_delta.png"), dpi=600); plt.close(fig)
print("F3f saved", flush=True)
