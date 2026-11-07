# Figure assembly (final revision)

Main figures were assembled at 183 mm width and 600 dpi. The final TIFF/PDF/PNG files are in
`05_Figures/` of the submission package.

## Main assembly scripts

- `makeFig1.R` ... `makeFig7.R`
- `makeFig8_fixed.R` (Fig 8 with the corrected panels c, f, g and h)

These write `analysis/figs_final/_v9/v9_Fig*`.

## Panel patch scripts (`patches/`)

- `make_fig1a.py` + `compose_fig1.py` — Fig 1a (complete accession labels, single role legend)
- `compose_fig2.py` — Fig 2h/2j rebuilt panels
- `make_fig5f.py` + `compose_fig5.py` — Fig 5b/5f/5g and the 5c table
- `compose_fig6.py` — Fig 6b forest plot
- `make_fig8_card.py` — Fig 8c evidence card

## New panels

- `make_f8g_v3.py` — Fig 8g (9-gene heatmap, decouple index, LDH at 24 h)
- `make_s13_fixed.py` — Fig 8h (two-panel W2 qPCR and dose-response)
- `make_f7f8f_panels.py` — Fig 7e/7f and Fig 8f

## Order of execution

1. Run the panel scripts needed for the figure.
2. Run the corresponding `makeFig*.R` assembly script.
3. Run the corresponding compose script (if any) to paste rebuilt panels into the final figure.
4. Export PNG/TIFF/PDF at 600 dpi.

The scripts use absolute paths from the analysis machine; the input tables are described in the
repository README and in Supplementary Tables S1/S6/S7.
