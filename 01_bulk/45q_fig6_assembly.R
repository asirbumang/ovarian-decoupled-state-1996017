
suppressMessages({library(ggplot2); library(patchwork); library(png); library(ragg)})
save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 90, dpi = 600) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi); print(plot); dev.off()
  ggsave(paste0(filename, ".png"), plot, width = w, height = h, dpi = 600)
}
FIG <- "D:/子宫内膜异位症2/analysis/figs_final/fig6"
rd <- function(f) readPNG(f)
mk <- function(f, cap) ggplot() + annotation_raster(rd(f), 0, 1, 0, 1) + theme_void() +
  labs(caption = cap) + theme(plot.caption = element_text(size = 4.6, hjust = 0))
pa <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F6a_direction.png", "a: direction summary (15 results)")
pb <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F6b_forests.png", "b: raw vs adjusted OS forests (TCGA only *)")
pc <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F6c_crude_resid.png", "c: crude vs adjusted p (-log10)")
pd_ <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F6def_main.png", "d-f: ML layer (Mime 101 / cross-check / Scissor)")
layout <- (pa / (pb + pc) / (pd_)) +
  plot_layout(heights = c(1.1, 0.75, 0.9)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(layout, file.path(FIG, "F6_combined_v3"), width_mm = 183, height_mm = 138)
cat("saved F6_combined_v3\nDONE 45q\n")
