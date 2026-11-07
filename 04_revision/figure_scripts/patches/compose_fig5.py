# -*- coding: utf-8 -*-
"""Fig 5: replace panels b, f, g and the c table with verified-isolated redraws."""
import os, sys, io
from PIL import Image
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

FIGDIR = r'最终\02_Figures_Main'
SRC = os.path.join(FIGDIR, 'Fig5.png')
AUD = r'_audit_r3'
im = Image.open(SRC).convert('RGB')

# c table
im.paste((255, 255, 255), (3355, 100, 4200, 806))
im.paste(Image.open(os.path.join(AUD, 'new_fig5c_table.png')).convert('RGB'), (3364, 117))

# b (panel d is at y < 2300, so the whole band below 1245 is free to the right of 1460)
im.paste((255, 255, 255), (8, 1245, 1450, 2312))
im.paste((255, 255, 255), (8, 2313, 1460, 2438))
pb = Image.open(os.path.join(AUD, 'new_fig5b.png')).convert('RGB').resize((1406, 1044), Image.LANCZOS)
im.paste(pb, (19, 1266))

# f (bottom-left block; panel d ends above y=2300)
im.paste((255, 255, 255), (8, 2439, 1720, 3615))
pf = Image.open(os.path.join(AUD, 'new_fig5f.png')).convert('RGB').resize((1400, 1155), Image.LANCZOS)
im.paste(pf, (18, 2448))

# g
im.paste((255, 255, 255), (2095, 2415, 4315, 3720))
pg = Image.open(os.path.join(AUD, 'new_fig5g.png')).convert('RGB').resize((2139, 1286), Image.LANCZOS)
im.paste(pg, (2161, 2420))

im.save(SRC, dpi=(600, 600))
im.save(os.path.join(FIGDIR, 'Fig5.tiff'), dpi=(600, 600), compression='tiff_lzw')
im.save(os.path.join(FIGDIR, 'Fig5.pdf'), 'PDF', resolution=600.0)
print('saved')
