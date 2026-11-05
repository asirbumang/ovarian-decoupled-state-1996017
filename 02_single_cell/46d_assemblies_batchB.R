
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
FIG <- "D:/子宫内膜异位症2/analysis/figs_final"
p3 <- list(
  mk(paste0(Q, "F3a_chr_heatmap.png"), "a chr profile"),
  mk(paste0(Q, "F3b_cnv_dist.png"), "b distribution (AUC 0.377)"),
  mk(paste0(Q, "F3c_horizontal.png"), "c contexts"),
  mk(paste0(Q, "F3d_modules_mal.png"), "d modules"),
  mk(paste0(Q, "F3e_cnv_lineage.png"), "e lineage CNV"),
  mk(paste0(Q, "F3f_chr_delta.png"), "f chr delta (mal vs non)"),
  mk(paste0(Q, "F3g_pt_mal.png"), "g patient rate"))
f3 <- (p3[[1]] + p3[[2]]) / (p3[[3]] + p3[[4]]) / (p3[[5]] + p3[[6]] + p3[[7]]) +
  plot_layout(heights = c(1.2, 0.9, 0.9)) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(f3, file.path(FIG, "fig3", "F3_combined_v3"), width_mm = 183, height_mm = 138)
p7 <- list(
  mk(paste0(Q, "F7a_b9.png"), "a B9 matrix"),
  mk(paste0(Q, "F7b_combo.png"), "b combos"),
  mk(paste0(Q, "F7c_drug.png"), "c drugs"),
  mk(paste0(Q, "F7d_signature.png"), "d signature"),
  mk(paste0(Q, "F7e_alpha.png"), "e alpha grid"),
  mk(paste0(Q, "F7f_evidence.png"), "f evidence matrix"))
f7 <- (p7[[1]] + p7[[2]]) / (p7[[3]] + p7[[4]]) / (p7[[5]] + p7[[6]]) +
  plot_layout(heights = c(1.05, 1.0, 0.85)) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(f7, file.path(FIG, "fig7", "F7_combined_v4"), width_mm = 183, height_mm = 140)
p8 <- list(
  mk(paste0(Q, "F8a_km.png"), "a KM"),
  mk(paste0(Q, "F8b_layers.png"), "b layers"),
  mk(paste0(Q, "F8c_honest.png"), "c honest"),
  mk(paste0(Q, "F11_consensus.png"), "d consensus"),
  mk(paste0(Q, "S11_subtype_char.png"), "e features"),
  mk(paste0(Q, "F8f_residual.png"), "f residual"))
f8 <- (p8[[1]] + p8[[2]]) / (p8[[3]] + p8[[4]]) / (p8[[5]] + p8[[6]]) +
  plot_layout(heights = c(1.0, 0.9, 0.85)) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(f8, file.path(FIG, "fig8", "F8_combined_v3"), width_mm = 183, height_mm = 138)
cat("DONE 46d v2\n")
