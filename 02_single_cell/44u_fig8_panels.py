# -*- coding: utf-8 -*-
"""44u_fig8_panels.py: 图8 —— a KM(dec_res 三分位) + b 临床各层 p 汇总 + c 诚实声明卡。"""
import numpy as np, pandas as pd, sys, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf-8')
plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "figure.dpi": 600, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
QRY = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\figs"

## a: KM（TCGA-OV dec_res 三分位；logrank p=0.060（26d 已算））
d = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L0_bulk\results\track2\TCGA_OV_residual.csv",
                usecols=["os_days", "event", "dec_res"], low_memory=False)
d = d[d["event"].notna() & d["os_days"].notna() & d["os_days"] > 0 & d["dec_res"].notna()]
q1, q2 = d["dec_res"].quantile([1 / 3, 2 / 3])
d["grp"] = np.where(d["dec_res"] > q2, "High", np.where(d["dec_res"] > q1, "Mid", "Low"))

def km(t, e):
    idx = np.argsort(t)
    t, e = t[idx], e[idx]
    times = np.unique(t)
    S = [1.0]; at = [0.0]
    for tt in times:
        n = np.sum(t >= tt)
        d_ = np.sum((e == 1) & (t == tt))
        if n > 0:
            S.append(S[-1] * (1 - d_ / n)); at.append(tt)
    return np.array(at), np.array(S)

fig, ax = plt.subplots(figsize=(66/25.4, 54/25.4))
for g, c in [("High", "#DC0000"), ("Mid", "#3C5488"), ("Low", "#B0B0B0")]:
    sub = d[d["grp"] == g]
    at, S = km(sub["os_days"].values, sub["event"].values)
    ax.step(at, S, where="post", color=c, linewidth=1.0, label=f"{g} (n={len(sub)})")
ax.set_xlabel("days", size=6); ax.set_ylabel("OS probability", size=6)
ax.tick_params(labelsize=5)
ax.legend(fontsize=4.8, frameon=False, loc="lower left")
ax.text(0.98, 0.97, "logrank p=0.060 (adjusted dec_res; 3-tier)", transform=ax.transAxes, ha="right",
        va="top", fontsize=5)
ax.set_title("TCGA-OV OS by composition-adjusted decouple", size=6)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F8a_km.png"), dpi=600)
plt.close(fig)
print("F8a saved; n:", len(d), flush=True)

## b: 临床各层 p（−log10；NS 界线）
lay = pd.DataFrame({
    "layer": ["TCGA OS (raw 3-tier, n=423)", "TCGA OS (adjusted, n=423)", "GSE156699 response (raw, n=88)",
              "GSE156699 response (adjusted, n=88)", "ER vs PR (all spots, n=12)", "ER vs PR (malignant spots, n=12)"],
    "p": [0.096, 0.060, 0.412, 0.923, 0.350, 0.394]})
fig, ax = plt.subplots(figsize=(66/25.4, 46/25.4))
y = np.arange(len(lay))[::-1]
bars = ax.barh(y, -np.log10(lay["p"]), color=["#3C5488" if p < 0.05 else "#B0B0B0" for p in lay["p"]], height=0.6)
ax.axvline(-np.log10(0.05), color="#DC0000", linestyle="--", linewidth=0.6)
ax.text(-np.log10(0.05) + 0.03, len(lay) - 0.4, "p=0.05", fontsize=4.6, color="#DC0000")
for y_, (_, r) in zip(y, lay.iterrows()):
    ax.text(-np.log10(r["p"]) + 0.02, y_, f"p={r['p']:.3f}", va="center", fontsize=4.6)
ax.set_yticks(y); ax.set_yticklabels(lay["layer"], fontsize=4.8)
ax.set_xlabel("-log10(p)", size=6); ax.tick_params(labelsize=5)
ax.set_title("Clinical layers: all p>0.05", size=5.8)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F8b_layers.png"), dpi=600)
plt.close(fig)
print("F8b saved", flush=True)

## c: 诚实声明卡（文字）
text = ("Honest conclusions\n"
        "• OS association: TCGA-only exploratory (continuous HR 1.17, p=0.025 adjusted;\n"
        "  3-tier logrank p=0.060); 4 external cohorts NS (pooled HR 0.91, p=0.14; I2=0%) - null\n"
        "• Treatment response: NS across 6 layers (chemo-naive/IDS, platinum, PARPi, ER/PR)\n"
        "• Composition adjustment: direction stable in OS (0.096->0.060) but reversed/weakened in response layers\n"
        "• What is robust: rare (1.2-2.9%), spatially-clustered, TF-controlled (NF-kB/IRF1/STAT1), drug-sensitive axis")
fig, ax = plt.subplots(figsize=(150/25.4, 40/25.4))
ax.axis("off")
ax.text(0.01, 0.95, text, fontsize=5.2, va="top", family="Arial", transform=ax.transAxes)
fig.tight_layout()
fig.savefig(os.path.join(QRY, "F8c_honest.png"), dpi=600)
plt.close(fig)
print("F8c saved", flush=True)
print("DONE 44u")
