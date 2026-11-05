
suppressMessages({library(ggplot2); library(patchwork); library(png); library(ragg)})
theme_set(theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        legend.title = element_text(size = 6.2), legend.text = element_text(size = 5.8)))
save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 95, dpi = 600) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi); print(plot); dev.off()
  ggsave(paste0(filename, ".png"), plot, width = w, height = h, dpi = 600)
}
FIG <- "D:/子宫内膜异位症2/analysis/figs_final/fig5"
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
rd <- function(f) readPNG(f)
pa <- ggplot() + annotation_raster(rd("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F5a_commot.png"), 0, 1, 0, 1) + theme_void() +
  labs(caption = "a: COMMOT (GSE274657; 3 samples; CellChatDB 19 pairs; 1000x label permutation; * p<0.05; Patient_2 skipped: no state)") + theme(plot.caption = element_text(size = 4.8, hjust = 0))
pb <- ggplot() + annotation_raster(rd("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F5b_nnls_delta.png"), 0, 1, 0, 1) + theme_void() +
  labs(caption = "b: NNLS delta (three Visium datasets incl. GSE274657; bootstrap 95% CI)") + theme(plot.caption = element_text(size = 4.8, hjust = 0))
pc <- ggplot() + annotation_raster(rd("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F5c_geomx.png"), 0, 1, 0, 1) + theme_void() +
  labs(caption = "c: GeoMx GSE276935 (566 AOIs; CD20 n=83, B-enriched, no histologic TLS confirmation)") + theme(plot.caption = element_text(size = 4.8, hjust = 0))
layout <- (pa + pb) / (pc + plot_spacer()) +
  plot_layout(heights = c(1.1, 0.9), widths = c(1.1, 0.95)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(layout, file.path(FIG, "F5_combined_v1"), width_mm = 183, height_mm = 105)
cat("saved F5_combined_v1\nDONE 44p\n")
