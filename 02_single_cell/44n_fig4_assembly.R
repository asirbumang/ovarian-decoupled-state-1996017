
suppressMessages({library(ggplot2); library(patchwork); library(png); library(ragg)})
theme_set(theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        legend.title = element_text(size = 6.2), legend.text = element_text(size = 5.8)))
save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 100, dpi = 600) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi); print(plot); dev.off()
  ggsave(paste0(filename, ".png"), plot, width = w, height = h, dpi = 600)
}
FIG <- "D:/子宫内膜异位症2/analysis/figs_final/fig4"
rd <- function(f) readPNG(f)
pa <- ggplot() + annotation_raster(rd("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F4a_triptych_v2.png"), 0, 1, 0, 1) + theme_void() +
  labs(caption = "black background; 1 mm = 10 spots (100 um pitch); red box -> zoom") +
  theme(plot.caption = element_text(size = 4.8, hjust = 0))
pb <- ggplot() + annotation_raster(rd("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F4b_domain_state.png"), 0, 1, 0, 1) + theme_void() +
  labs(caption = "state rate per NMF domain (sample-specific domains)") + theme(plot.caption = element_text(size = 4.8, hjust = 0))
pc <- ggplot() + annotation_raster(rd("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F4cd_stats_v2.png"), 0, 1, 0, 1) + theme_void() +
  labs(caption = "Moran's I (k=6, row-std) 14/14 p<0.05; neighbourhood OR 7/12 p<0.05") + theme(plot.caption = element_text(size = 4.8, hjust = 0))
layout <- (pa + pb) / (pc + plot_spacer()) +
  plot_layout(heights = c(1.2, 0.9), widths = c(1.5, 0.85)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(layout, file.path(FIG, "F4_combined_v2"), width_mm = 183, height_mm = 105)
cat("saved F4_combined_v2\nDONE 44n\n")
