
suppressMessages({library(ggplot2); library(patchwork); library(png); library(ragg)})
save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 90, dpi = 600) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi); print(plot); dev.off()
  ggsave(paste0(filename, ".png"), plot, width = w, height = h, dpi = 600)
}
FIG <- "D:/子宫内膜异位症2/analysis/figs_final/fig3"
rd <- function(f) readPNG(f)
mk <- function(f, cap) ggplot() + annotation_raster(rd(f), 0, 1, 0, 1) + theme_void() +
  labs(caption = cap) + theme(plot.caption = element_text(size = 4.6, hjust = 0))
pa <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F3a_chr_heatmap.png", "chr-level CNV profile (median per cell; sorted)")
pb <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F3b_cnv_dist.png", "CNV score distribution; AUC=0.377")
pc <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F3c_horizontal.png", "three-context decouple (all-cell)")
pd_ <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F3d_modules_mal.png", "modules by malignancy")
pe_ <- mk("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F3e_cnv_lineage.png", "CNV by lineage")
layout <- (pa / (pb + pc)) / (pd_ + pe_) +
  plot_layout(heights = c(1.25, 0.75, 0.85)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(layout, file.path(FIG, "F3_combined_v6"), width_mm = 183, height_mm = 128)
cat("saved F3_combined_v6\nDONE 45j\n")
