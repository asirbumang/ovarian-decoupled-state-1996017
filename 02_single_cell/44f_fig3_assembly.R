
suppressMessages({library(ggplot2); library(patchwork); library(png); library(ragg)})
theme_set(theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        legend.title = element_text(size = 6.2), legend.text = element_text(size = 5.8)))
save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 90, dpi = 600) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi); print(plot); dev.off()
  ggsave(paste0(filename, ".png"), plot, width = w, height = h, dpi = 600)
}
FIG <- "D:/子宫内膜异位症2/analysis/figs_final/fig3"
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
img_a <- readPNG("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F3a_cropped.png")
img_b <- readPNG("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F3b_cnv_dist.png")
img_c <- readPNG("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F3c_horizontal.png")
pa <- ggplot() + annotation_raster(img_a, xmin = 0, xmax = 1, ymin = 0, ymax = 1) + theme_void() +
  labs(caption = "inferCNV subclusters (epi 22,719 + ref 20,000; no-HMM)") +
  theme(plot.caption = element_text(size = 4.8, hjust = 0))
pb <- ggplot() + annotation_raster(img_b, xmin = 0, xmax = 1, ymin = 0, ymax = 1) + theme_void() +
  labs(caption = "all-cell decouple; within-epi dirs in text") +
  theme(plot.caption = element_text(size = 4.8, hjust = 0))
pc <- ggplot() + annotation_raster(img_c, xmin = 0, xmax = 1, ymin = 0, ymax = 1) + theme_void()
layout <- (pa + plot_spacer()) / (pb + pc) +
  plot_layout(heights = c(1.55, 1.0), widths = c(1.05, 1.0)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(layout, file.path(FIG, "F3_combined_v1"), width_mm = 183, height_mm = 110)
cat("saved F3_combined_v1\nDONE 44f\n")
