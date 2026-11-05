
suppressMessages({library(ggplot2); library(patchwork); library(png); library(ragg)})
save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 90, dpi = 600) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi); print(plot); dev.off()
  ggsave(paste0(filename, ".png"), plot, width = w, height = h, dpi = 600)
}
FIG <- "D:/子宫内膜异位症2/analysis/figs_final/fig8"
rd <- function(f) readPNG(f)
mk <- function(f, cap) ggplot() + annotation_raster(rd(f), 0, 1, 0, 1) + theme_void() +
  labs(caption = cap) + theme(plot.caption = element_text(size = 4.6, hjust = 0))
pa <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F8a_km.png", "a: TCGA-OV OS (tertiles)")
pb <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F8b_layers.png", "b: clinical layers (all p>0.05)")
pc <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F8c_honest.png", "c: honest conclusions")
pd_ <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F11_consensus.png", "d: consensus subtyping (k=2; residual p=0.025: high more >20mm)")
pe_ <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/S11_subtype_char.png", "e: subtype features (M1/M3/M4/ICC45 high; M2 NS) + carrier subpopulations (IGHG+ plasma 5.2%, FOS+ CAF 4.0%)")
layout <- (pa + pb) / (pc + pd_) / (pe_) +
  plot_layout(heights = c(1.0, 0.9, 0.85)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(layout, file.path(FIG, "F8_combined_v2"), width_mm = 183, height_mm = 138)
cat("saved F8_combined_v2\nDONE 45x\n")
