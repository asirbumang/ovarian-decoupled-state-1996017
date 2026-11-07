# -*- coding: utf-8 -*-
"""Rebuild Fig 8 panel c card at the exact block bbox (1727 x 917 px) and zero out
the leftover 'see text / Table 7' placeholder."""
import os, sys, io
from PIL import Image, ImageDraw, ImageFont
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
AUD = r'_audit_r3'

CW, CH = 1727, 917
img = Image.new('RGB', (CW, CH), 'white')
d = ImageDraw.Draw(img)
FB = r'C:\Windows\Fonts\arialbd.ttf'
FR = r'C:\Windows\Fonts\arial.ttf'
f_head = ImageFont.truetype(FB, 54)
f_item = ImageFont.truetype(FR, 45)

GREEN, GREEN_BG = (0x0F, 0x8F, 0x7F), (0xE6, 0xF3, 0xEF)
RED, RED_BG = (0xC0, 0x39, 0x2B), (0xFB, 0xE9, 0xE5)

blocks = [
    ('ROBUST / reproducible', GREEN, GREEN_BG,
     ['Carrier-enriched (adjusted)', 'Spatially clustered (14/14)', 'NF-\u03baB/IRF1 axis (computational)']),
    ('NOT REPRODUCED', RED, RED_BG,
     ['External-cohort survival null', 'Response: not significant',
      'ML: no generalisation', 'Ovarian drug signal: none']),
]
pad, row_h, bar_w, bar_h = 20, 97, 58, 26
y = 0
for title, col, bg, items in blocks:
    h = row_h * (len(items) + 1) + 8
    d.rectangle([0, y, CW - 1, y + h], fill=bg)
    d.text((pad, y + 6), title, font=f_head, fill=col)
    yy = y + row_h + 2
    for it in items:
        d.rounded_rectangle([pad, yy + (row_h - bar_h) // 2 + 6,
                             pad + bar_w, yy + (row_h + bar_h) // 2 + 6],
                            radius=11, fill=col)
        d.text((pad + bar_w + 18, yy + 10), it, font=f_item, fill=col)
        yy += row_h
    y += h + 6
img.save(os.path.join(AUD, 'new_fig8c_card.png'))
print('card', img.size)
