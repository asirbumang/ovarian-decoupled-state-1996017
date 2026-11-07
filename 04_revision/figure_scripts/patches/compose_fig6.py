# -*- coding: utf-8 -*-
"""Composite the new Fig 6b forest into Fig 6."""
import os, sys, io
from PIL import Image
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

FIGDIR = r'最终\02_Figures_Main'
SRC = os.path.join(FIGDIR, 'Fig6.png')
AUD = r'_audit_r3'
im = Image.open(SRC).convert('RGB')
print('Fig6', im.size)
im.paste((255, 255, 255), (1512, 320, 3266, 1012))
nb = Image.open(os.path.join(AUD, 'new_fig6b.png')).convert('RGB').resize((1460, 660), Image.LANCZOS)
im.paste(nb, (1522, 342))
im.save(SRC, dpi=(600, 600))
im.save(os.path.join(FIGDIR, 'Fig6.tiff'), dpi=(600, 600), compression='tiff_lzw')
im.save(os.path.join(FIGDIR, 'Fig6.pdf'), 'PDF', resolution=600.0)
print('saved')
