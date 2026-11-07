# -*- coding: utf-8 -*-
"""Composite the rebuilt Fig 1 panel a."""
import os, sys, io
from PIL import Image
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

FIGDIR = r'最终\02_Figures_Main'
SRC = os.path.join(FIGDIR, 'Fig1.png')
AUD = r'_audit_r3'
im = Image.open(SRC).convert('RGB')
print('Fig1', im.size)
im.paste((255, 255, 255), (0, 0, 4322, 1590))
pa = Image.open(os.path.join(AUD, 'new_fig1a.png')).convert('RGB')
print('panel a', pa.size)
im.paste(pa, (0, 0))
im.save(SRC, dpi=(600, 600))
im.save(os.path.join(FIGDIR, 'Fig1.tiff'), dpi=(600, 600), compression='tiff_lzw')
im.save(os.path.join(FIGDIR, 'Fig1.pdf'), 'PDF', resolution=600.0)
print('saved')
