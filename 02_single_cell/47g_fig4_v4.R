
suppressMessages({library(ggplot2); library(patchwork); library(png); library(ragg)})
save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 90, dpi = 600) {
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
pa <- mk(paste0(Q, "F4a_triptych_v2.png"), "a triptych (1mm scale; red box zoom)")
pb <- mk(paste0(Q, "F4b_domain_state.png"), "b domain x state")
pc <- mk(paste0(Q, "F4cd_stats_v2.png"), "c Moran + neighbourhood OR")
pd_ <- mk(paste0(Q, "F4de_coloc.png"), "d-e co-localisation (heterogeneous)")
pf <- mk(paste0(Q, "F5c_geomx.png"), "e GeoMx segments (566 AOIs; p=3.5e-15)")
layout <- (pa + pb) / (pc + pd_) / (pf + plot_spacer()) +
  plot_layout(heights = c(1.05, 0.95, 0.75)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
FIG <- "D:/子宫内膜异位症2/analysis/figs_final/fig4"
save_pub_r(layout, file.path(FIG, "F4_combined_v4"), width_mm = 183, height_mm = 128)
cat("saved F4_combined_v4\nDONE 47g v2\n")
