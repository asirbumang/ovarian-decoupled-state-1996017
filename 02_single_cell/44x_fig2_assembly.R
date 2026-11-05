
suppressMessages({library(ggplot2); library(patchwork); library(png); library(ragg)})
theme_set(theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35)))
save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 90, dpi = 600) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi); print(plot); dev.off()
  ggsave(paste0(filename, ".png"), plot, width = w, height = h, dpi = 600)
}
FIG <- "D:/子宫内膜异位症2/analysis/figs_final/fig2"
rd <- function(f) readPNG(f)
mk <- function(f, cap) ggplot() + annotation_raster(rd(f), 0, 1, 0, 1) + theme_void() +
  labs(caption = cap) + theme(plot.caption = element_text(size = 4.6, hjust = 0))
pa <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F2a_state_umap.png", "a: UMAP with state (black)")
pa <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F2a_state_umap.png", "UMAP with state highlighted (black)")
pb <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F2c_quadrant.png", "Definition: M1z >= P75 & M2z/M3z <= P25 (within source)")
pc <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F2d_rarity.png", "Rarity across 5 sources (95% CI)")
pd_ <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F2d_or.png", "Lineage OR (Fisher, 95% CI)")
pe_ <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F2e_modules.png", "Module distributions (state vs rest)")
pf <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F2f_comp.png", "Composition: state vs all cells")
top <- (pa + pb + pc)
bottom <- (pd_ + pe_ + pf)
layout <- (top / bottom) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(layout, file.path(FIG, "F2_combined_v2"), width_mm = 183, height_mm = 96)
cat("saved F2_combined_v2\nDONE 44x\n")
