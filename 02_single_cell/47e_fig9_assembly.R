
suppressMessages({library(ggplot2); library(patchwork); library(png); library(ragg)})
save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 90, dpi = 600) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi); print(plot); dev.off()
  ggsave(paste0(filename, ".png"), plot, width = w, height = h, dpi = 600)
}
rd <- function(f) readPNG(f)
# S12 v3 四联裁剪为 2x2 组件（直接整图作为 a-d 组图槽位不可；改为整图分槽）
s12 <- rd("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/supp/S12_w1.png")
s13 <- rd("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/supp/S13_w2.png")
pa <- ggplot() + annotation_raster(s12, 0, 1, 0, 1) + theme_void() +
  labs(caption = "a-d: W1 joint induction (BSO+TNFa: M1 >> M2/M3; decouple index +2.6-2.8; LDH 24%)") +
  theme(plot.caption = element_text(size = 5.4, hjust = 0))
pe <- ggplot() + annotation_raster(s13, 0, 1, 0, 1) + theme_void() +
  labs(caption = "e-g: W2 (JSH blocks NF-kB induction 90-93%; bos sensitivity not enhanced)") +
  theme(plot.caption = element_text(size = 5.4, hjust = 0))
layout <- (pa / pe) + plot_layout(heights = c(1.15, 0.85)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
FIG <- "D:/子宫内膜异位症2/analysis/figs_final/fig9"
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
save_pub_r(layout, file.path(FIG, "F9_combined_v1"), width_mm = 183, height_mm = 128)
cat("saved F9_combined_v1\nDONE 47e\n")
