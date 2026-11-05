# -*- coding: utf-8 -*-
"""45r_fig7_panel_d.py: 图7 d 面板（54 基因签名：AUC + bootstrap 频率 top12）。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
OUT = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\cnv\b12"
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"
sm = pd.read_csv(os.path.join(OUT, "sig_freq.csv")).sort_values("freq", ascending=False).head(12).iloc[::-1]
fig, axes = plt.subplots(1, 2, figsize=(95/25.4, 50/25.4))
axes[0].barh(range(len(sm)), sm["freq"], color="#E64B35", height=0.6)
axes[0].set_yticks(range(len(sm)))
axes[0].set_yticklabels(sm["gene"], fontsize=5.4)
axes[0].set_xlabel("bootstrap frequency (30x L1)", size=6)
axes[0].tick_params(labelsize=5.4)
axes[0].set_title("d 54-gene signature (freq>0.5); LR AUC=0.997", size=5.8)
mech = {"IRF1": True, "NFKBIA": True, "PKM": True, "LDHA": True, "CXCL8": True, "CXCL3": True, "RHOA": True, "ACTG1": True}
axes[1].axis("off")
axes[1].text(0.02, 0.95, "Mechanistic convergence of the signature:\n"
                        "• IRF1, NFKBIA — B9 consensus TFs (KO revert)\n"
                        "• CXCL8, CXCL3 — chemokine axis (COMMOT CCL2-CCR2; drug screen)\n"
                        "• PKM, LDHA — M2 metabolic/GSH module\n"
                        "• RHOA, ACTG1, CD9 — M3 mito-cytoskeleton module\n"
                        "• LDHA, IRF1, ACTG1 — overlap with master 7-gene signature",
             transform=axes[1].transAxes, fontsize=5.4, va="top")
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F7d_signature.png"), dpi=600)
plt.close(fig)
print("F7d saved", flush=True)
print("DONE 45r")
