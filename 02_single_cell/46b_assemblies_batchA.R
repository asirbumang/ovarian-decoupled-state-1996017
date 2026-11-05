
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
  labs(caption = cap) + theme(plot.caption = element_text(size = 4.6, hjust = 0))
Q <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/"

## 图2 v3（a-h 8 面板）
p2 <- list(
  mk(paste0(Q, "F2a_state_umap.png"), "a UMAP state"),
  mk(paste0(Q, "F2c_quadrant.png"), "b definition"),
  mk(paste0(Q, "F2d_rarity.png"), "c rarity"),
  mk(paste0(Q, "F2d_or.png"), "d lineage OR"),
  mk(paste0(Q, "F2e_modules.png"), "e modules"),
  mk(paste0(Q, "F2f_comp.png"), "f composition"),
  mk(paste0(Q, "F2g_carriers.png"), "g-h carrier subpopulations"))
f2 <- (p2[[1]] + p2[[2]] + p2[[3]]) / (p2[[4]] + p2[[5]] + p2[[6]]) / (p2[[7]] + plot_spacer() + plot_spacer()) +
  plot_layout(heights = c(1, 1, 0.85)) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(f2, "D:/子宫内膜异位症2/analysis/figs_final/fig2/F2_combined_v3", width_mm = 183, height_mm = 128)

## 图4 v3（a-f 6 面板）
p4 <- list(
  mk(paste0(Q, "F4a_triptych_v2.png"), "a triptych"),
  mk(paste0(Q, "F4b_domain_state.png"), "b domain x state"),
  mk(paste0(Q, "F4cd_stats_v2.png"), "c Moran + OR"),
  mk(paste0(Q, "F4de_coloc.png"), "d-e co-localisation OR"))
f4 <- (p4[[1]] + p4[[2]]) / (p4[[3]] + p4[[4]]) +
  plot_layout(heights = c(1.1, 1.0)) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(f4, "D:/子宫内膜异位症2/analysis/figs_final/fig4/F4_combined_v3", width_mm = 183, height_mm = 118)

## 图5 v2（a-d 4 面板：a + b / c + d）
p5 <- list(
  mk(paste0(Q, "F5a_commot.png"), "a COMMOT heatmap"),
  mk(paste0(Q, "F5b_nnls_delta.png"), "b NNLS delta"),
  mk(paste0(Q, "F5c_geomx.png"), "c GeoMx segments"),
  mk(paste0(Q, "F5d_lr_count.png"), "d significant LR counts"))
f5 <- (p5[[1]] + p5[[2]]) / (p5[[3]] + p5[[4]]) +
  plot_layout(heights = c(1.05, 0.85)) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(f5, "D:/子宫内膜异位症2/analysis/figs_final/fig5/F5_combined_v2", width_mm = 183, height_mm = 112)
cat("DONE 46b\n")
