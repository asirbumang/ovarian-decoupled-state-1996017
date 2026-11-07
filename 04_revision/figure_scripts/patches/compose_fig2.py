# -*- coding: utf-8 -*-
"""Paste regenerated panels h and j into Fig 2, preserving panel i's title tail."""
import os, sys, io
from PIL import Image
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

FIGDIR = r'最终\02_Figures_Main'
SRC = os.path.join(FIGDIR, 'Fig2.png')
AUD = r'_audit_r3'

im = Image.open(SRC).convert('RGB')
print('Fig2', im.size)

# --- panel i title tail that must survive the erase of the j row
tail = im.crop((1240, 3830, 1470, 3928))
tail.save(os.path.join(AUD, 'fig2_i_title_tail.png'))

# --- erase + paste panel j
im.paste((255, 255, 255), (1380, 3820, 3630, 5075))
im.paste(tail, (1240, 3830))
pj = Image.open(os.path.join(AUD, 'new_fig2j.png')).convert('RGB').resize((2160, 1255), Image.LANCZOS)
im.paste(pj, (1470, 3820))
print('panel j pasted at x=1470..3630')

# --- erase + paste panel h
im.paste((255, 255, 255), (2035, 2900, 3810, 3525))
ph = Image.open(os.path.join(AUD, 'new_fig2h.png')).convert('RGB').resize((1775, 625), Image.LANCZOS)
im.paste(ph, (2035, 2900))
print('panel h pasted at x=2035..3810')

im.save(SRC, dpi=(600, 600))
im.save(os.path.join(FIGDIR, 'Fig2.tiff'), dpi=(600, 600), compression='tiff_lzw')
im.save(os.path.join(FIGDIR, 'Fig2.pdf'), 'PDF', resolution=600.0)
print('saved')
