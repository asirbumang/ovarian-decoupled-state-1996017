suppressMessages({library(ggplot2); library(patchwork); library(png); library(ragg)})
save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 128, dpi = 600) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi); print(plot); dev.off()
  ggsave(paste0(filename, ".png"), plot, width = w, height = h, dpi = 600)
}
rd <- function(f) readPNG(f)
mk <- function(f, cap) ggplot() + annotation_raster(rd(f), 0, 1, 0, 1) + theme_void() +
  labs(caption = cap) + theme(plot.caption = element_text(size = 5.4, hjust = 0))
Q <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/"
p <- list(
  mk(paste0(Q, "F2a_state_umap.png"), "a UMAP state"),
  mk(paste0(Q, "F2c_quadrant.png"), "b definition"),
  mk(paste0(Q, "F2d_rarity.png"), "c rarity (5 sources)"),
  mk(paste0(Q, "F2d_or.png"), "d lineage OR"),
  mk(paste0(Q, "F2e_modules.png"), "e modules (state vs rest)"),
  mk(paste0(Q, "F2f_comp.png"), "f composition"),
  mk(paste0(Q, "F2g_carriers.png"), "g carrier subpopulations (state fraction)"),
  mk(paste0(Q, "F2g2_module_structure.png"), "h module structure within carriers (IGHG+ plasma: M1-driven; FOS/CTGF+ CAF: M2/M3-driven)"))
layout <- (p[[1]] + p[[2]] + p[[3]]) / (p[[4]] + p[[5]] + p[[6]]) / (p[[7]] + p[[8]]) +
  plot_layout(heights = c(1, 1, 0.42)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(layout, "D:/子宫内膜异位症2/analysis/figs_final/fig2/F2_combined_v6", width_mm = 183, height_mm = 128)
cat("DONE 47i fix\n")
