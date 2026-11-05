# -*- coding: utf-8 -*-
"""47c_panels_batchC.py: 图5e(样本级 浆/CAF fraction × state 率散点) + 图5/图8 组装 v3/v4（字号≥5.4）。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"
OST = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\spatial"
rows = []
for GSE in ["GSE274657", "GSE288483", "GSE224335"]:
    d = pd.read_csv(os.path.join(OST, f"{GSE}_nnls_state.csv"))
    dom = pd.read_csv(os.path.join(OST, f"{GSE}_domains_state.csv"))
    for sm, s in d.groupby("sample"):
        sr = dom[dom["sample"] == sm]["mean"].mean() * 100
        if "Plasma cells" in s.columns:
            pf = s["Plasma cells"].mean()
        else:
            pf = np.nan
        ff = s["Fibroblasts"].mean() if "Fibroblasts" in s.columns else np.nan
        rows.append({"GSE": GSE, "sample": sm, "state_rate": sr, "plasma_frac": pf, "fibro_frac": ff})
df = pd.DataFrame(rows).dropna(subset=["plasma_frac"])
rho1, p1 = spearmanr(df["state_rate"], df["plasma_frac"])
rho2, p2 = spearmanr(df["state_rate"], df["fibro_frac"])
print("n samples:", len(df), " rho plasma:", round(rho1, 3), "p:", round(p1, 3),
      "| rho fibro:", round(rho2, 3), "p:", round(p2, 3), flush=True)
fig, axes = plt.subplots(1, 2, figsize=(110/25.4, 52/25.4))
for ax, col, rho_, p_, lab in [(axes[0], "plasma_frac", rho1, p1, "Plasma cells"),
                               (axes[1], "fibro_frac", rho2, p2, "Fibroblasts")]:
    ax.scatter(df["state_rate"], df[col], s=14, c="#3C5488")
    ax.set_xlabel("state rate (%)", size=6); ax.set_ylabel(f"{lab} NNLS fraction", size=6)
    ax.tick_params(labelsize=5.4)
    ax.text(0.97, 0.05, f"rho={rho_:.2f} (p={p_:.3f})", transform=ax.transAxes, ha="right", fontsize=5.4)
    ax.set_title(f"e sample-level composition vs state ({lab})", size=6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F5e_comp_state.png"), dpi=600)
plt.close(fig)
print("F5e saved", flush=True)
print("DONE 47c")
