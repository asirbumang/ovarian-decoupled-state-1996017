# -*- coding: utf-8 -*-
"""Rebuild Fig 1 panel a with complete accession labels and a de-duplicated legend.
Cell colours are sampled from the existing strip so the palette is preserved."""
import os, sys, io
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
plt.rcParams.update({'font.family': 'Arial', 'pdf.fonttype': 42})

SRC = r'最终\02_Figures_Main\Fig1.png'
AUD = r'_audit_r3'
a = np.asarray(Image.open(SRC).convert('RGB')).astype(int)

X0, X1 = 521, 3796
NC = 18
rows_y = {'tech': (412, 593), 'role': (601, 782), 'layer': (789, 970)}
cw = (X1 - X0) / NC


def sample(y0, y1, i):
    x = int(X0 + (i + 0.5) * cw)
    y = int((y0 + y1) / 2)
    return tuple(a[y, x] / 255.0)


TECH = [sample(*rows_y['tech'], i) for i in range(NC)]
ROLE = [sample(*rows_y['role'], i) for i in range(NC)]
LAYER = [sample(*rows_y['layer'], i) for i in range(NC)]
for i in range(NC):
    print(i, 'tech', tuple(round(v * 255) for v in TECH[i]),
          'role', tuple(round(v * 255) for v in ROLE[i]),
          'layer', tuple(round(v * 255) for v in LAYER[i]))

ACC = ['GSE266577', 'GSE222556', 'GSE201047', 'GSE180661', 'GSE293459',
       'GSE156699', 'GSE26193', 'GSE13876', 'GSE14764', 'GSE63885',
       'GSE26712', 'TCGA-OV', 'GSE274657', 'GSE288483', 'GSE224335',
       'GSE211956', 'GSE189843', 'GSE276935']

# legend keys in the same order as the columns they first appear in
LAYER_KEYS = [('scRNA', 0), ('bulk', 5), ('Visium', 12), ('GeoMx', 17)]
ROLE_KEYS = [('State/cancer', 0), ('Metastasis', 1), ('Chemo-naive', 2), ('Site', 3),
             ('Persister', 4), ('Response', 5),
             ('OS', 6), ('OS+Platinum', 9), ('Tumour/Normal', 10), ('OS+Response', 11),
             ('Spatial', 12), ('TLS', 17)]
TECH_KEYS = [("10x 3' v3.1", 0), ('10x v2', 2), ('10x sorted', 3), ('snRNA', 4),
             ('RNA-seq', 5), ('GPL570', 6), ('GPL96', 7), ('10x Visium', 12),
             ('GeoMx DSP', 17)]

W_MM, H_MM = 183.0, 67.3
fig = plt.figure(figsize=(W_MM / 25.4, H_MM / 25.4))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W_MM); ax.set_ylim(0, H_MM)
ax.axis('off')

# geometry in mm
x_left = 22.0
strip_w = 148.0
label_base = 55.0          # strip top
row_h = 9.0
col_w = strip_w / NC

for r, (name, ys) in enumerate([('tech', TECH), ('role', ROLE), ('layer', LAYER)]):
    y = label_base - (r + 1) * row_h
    for i in range(NC):
        ax.add_patch(Rectangle((x_left + i * col_w, y), col_w * 0.93, row_h * 0.93,
                               facecolor=ys[i], edgecolor='none'))
    ax.text(x_left - 1.5, y + row_h / 2, name, ha='right', va='center', fontsize=5.6)

for i, acc in enumerate(ACC):
    ax.text(x_left + (i + 0.5) * col_w, label_base + 0.8, acc, rotation=45,
            rotation_mode='anchor', ha='left', va='bottom', fontsize=5.4)

# legends: layer | role (6) | role (6) | tech
def legend(items, x0, y0, title):
    ax.text(x0, y0, title, fontsize=5.8, fontweight='bold', va='top')
    yy = y0 - 2.2
    for lab, col in items:
        ax.add_patch(Rectangle((x0, yy - 1.5), 2.4, 1.5, facecolor=col, edgecolor='none'))
        ax.text(x0 + 3.0, yy - 0.8, lab, fontsize=5.4, va='center')
        yy -= 2.55

yb = label_base - 3 * row_h - 1.8
legend([(k, LAYER[c]) for k, c in LAYER_KEYS], 1.0, yb, 'layer')
legend([(k, ROLE[c]) for k, c in ROLE_KEYS[:6]], 32.0, yb, '')
legend([(k, ROLE[c]) for k, c in ROLE_KEYS[6:]], 62.0, yb, '')
ax.text(47.0, yb, 'role', fontsize=5.8, fontweight='bold', va='top', ha='center')
legend([(k, TECH[c]) for k, c in TECH_KEYS], 106.0, yb, 'tech')

ax.text(1.5, H_MM - 1.2, 'a', fontsize=8, fontweight='bold', ha='left', va='top')

out = os.path.join(AUD, 'new_fig1a.png')
fig.savefig(out, dpi=600)
plt.close(fig)
print('saved', out, Image.open(out).size)
