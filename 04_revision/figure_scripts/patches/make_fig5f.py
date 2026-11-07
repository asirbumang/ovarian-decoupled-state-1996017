# -*- coding: utf-8 -*-
"""Redraw Fig 5 panel f (GeoMx) without scientific notation and without the clipped label."""
import os, sys, io
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
plt.rcParams.update({'font.family': 'Arial', 'pdf.fonttype': 42})

LAB = ['Tumor', 'Stroma', 'CD20', 'Immune', 'Full ROI']
VAL = [-1.26, -0.233, 0.054, 0.324, 0.347]
COL = ['#C0392B', '#56B4E9', '#3B4F8C', '#009E73', '#999999']

W, H = 1380 / 600.0, 1160 / 600.0
fig, ax = plt.subplots(figsize=(W, H))
y = np.arange(len(LAB))
ax.barh(y, VAL, height=0.62, color=COL)
ax.axvline(0, color='k', lw=0.6)
for yi, v in zip(y, VAL):
    ax.text(v + (0.05 if v >= 0 else -0.05), yi, '%+.2f' % v,
            va='center', ha='left' if v >= 0 else 'right', fontsize=5.4)
ax.set_yticks(y); ax.set_yticklabels(LAB, fontsize=5.6)
ax.set_xlim(-1.55, 0.62); ax.tick_params(labelsize=5.4)
ax.set_xlabel('median decouple score per AOI', fontsize=5.8)
ax.set_title('GeoMx compartments (GSE276935, 566 AOIs)', fontsize=6.0, loc='left', pad=3)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
fig.tight_layout(pad=0.3)
fig.savefig(r'_audit_r3\new_fig5f.png', dpi=600)
plt.close(fig)
print('panel f saved')
