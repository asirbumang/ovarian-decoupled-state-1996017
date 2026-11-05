
suppressMessages({library(ggplot2); library(patchwork); library(png); library(ragg)})
theme_set(theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35)))
save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 95, dpi = 600) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi); print(plot); dev.off()
  ggsave(paste0(filename, ".png"), plot, width = w, height = h, dpi = 600)
}
FIG <- "D:/子宫内膜异位症2/analysis/figs_final/fig6"
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
rd <- function(f) readPNG(f)
pa <- ggplot() + annotation_raster(rd("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F6a_direction.png"), 0, 1, 0, 1) + theme_void() +
  labs(caption = "a: evidence direction matrix (15 results; red=opposite/negative, grey=NS, blue=consistent)") + theme(plot.caption = element_text(size = 4.8, hjust = 0))
pb <- ggplot() + annotation_raster(rd("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F6b_forests.png"), 0, 1, 0, 1) + theme_void() +
  labs(caption = "b: raw tertile vs composition-adjusted dec_res (5 cohorts; TCGA only significant: HR 1.17, p=0.025, exploratory)") + theme(plot.caption = element_text(size = 4.8, hjust = 0))
pc <- ggplot() + annotation_raster(rd("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F6c_crude_resid.png"), 0, 1, 0, 1) + theme_void() +
  labs(caption = "c: crude vs composition-adjusted association p-values") + theme(plot.caption = element_text(size = 4.8, hjust = 0))
layout <- (pa / plot_spacer()) / (pb + pc) +
  plot_layout(heights = c(1.3, 0.16, 0.9), widths = c(1.35, 0.85)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(layout, file.path(FIG, "F6_combined_v1"), width_mm = 183, height_mm = 122)
cat("saved F6_combined_v1\nDONE 44r\n")
