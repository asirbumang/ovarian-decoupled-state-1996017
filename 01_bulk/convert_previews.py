# -*- coding: utf-8 -*-
from PIL import Image
import os, glob

d = r"D:\子宫内膜异位症2\analysis\L0_bulk\results"
for t in sorted(glob.glob(os.path.join(d, "*.tiff"))):
    im = Image.open(t)
    im.thumbnail((1400, 1400), Image.LANCZOS)
    out = os.path.join(d, os.path.basename(t).replace(".tiff", "_preview.png"))
    im.convert("RGB").save(out, "PNG")
    print(os.path.basename(out), im.size)
