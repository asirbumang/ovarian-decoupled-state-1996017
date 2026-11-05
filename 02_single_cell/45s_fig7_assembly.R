
suppressMessages({library(ggplot2); library(patchwork); library(png); library(ragg)})
save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 90, dpi = 600) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi); print(plot); dev.off()
  ggsave(paste0(filename, ".png"), plot, width = w, height = h, dpi = 600)
}
FIG <- "D:/子宫内膜异位症2/analysis/figs_final/fig7"
rd <- function(f) readPNG(f)
mk <- function(f, cap) ggplot() + annotation_raster(rd(f), 0, 1, 0, 1) + theme_void() +
  labs(caption = cap) + theme(plot.caption = element_text(size = 4.6, hjust = 0))
pa <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F7a_b9.png", "a: heuristic TF perturbation (* permutation p<0.05)")
pb <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F7b_combo.png", "b: combination KO (dashed=best single)")
pc <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F7c_drug.png", "c: CTRPv2 sensitivity (rho; q<0.25 n=277)")
pd_ <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F7d_signature.png", "d: 54-gene ML signature (mechanistic convergence)")
layout <- (pa + pb) / (pc + pd_) +
  plot_layout(heights = c(1.15, 1.0)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(layout, file.path(FIG, "F7_combined_v3"), width_mm = 183, height_mm = 112)
cat("saved F7_combined_v3\nDONE 45s\n")
