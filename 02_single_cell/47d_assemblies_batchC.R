
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
FIG <- "D:/子宫内膜异位症2/analysis/figs_final"
p5 <- list(
  mk(paste0(Q, "F5a_commot.png"), "a COMMOT heatmap (log2 ratio; * perm p<0.05)"),
  mk(paste0(Q, "F5b_nnls_delta.png"), "b NNLS composition delta (95% CI)"),
  mk(paste0(Q, "F5d_lr_count.png"), "c significant LR pairs per sample"),
  mk(paste0(Q, "F5c_geomx.png"), "d GeoMx segment decouple"),
  mk(paste0(Q, "F5e_comp_state.png"), "e sample-level composition vs state (rho=0.34/0.19 NS, n=14)"))
f5 <- (p5[[1]] + p5[[2]]) / (p5[[3]] + p5[[4]] + p5[[5]]) +
  plot_layout(heights = c(1.1, 0.8), widths = c(1, 0.9, 0.9)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(f5, file.path(FIG, "fig5", "F5_combined_v3"), width_mm = 183, height_mm = 118)
p8 <- list(
  mk(paste0(Q, "F8a_km.png"), "a OS KM"),
  mk(paste0(Q, "F8b_layers.png"), "b clinical layers"),
  mk(paste0(Q, "F8c_honest.png"), "c honest conclusions"),
  mk(paste0(Q, "F11_consensus.png"), "d consensus subtyping"),
  mk(paste0(Q, "S11_subtype_char.png"), "e subtype features"),
  mk(paste0(Q, "F8f_residual.png"), "f residual tumour"),
  mk(paste0(Q, "supp/S12_w1.png"), "g W1: joint induction (qPCR/CCK8/LDH)"),
  mk(paste0(Q, "supp/S13_w2.png"), "h W2: JSH blocks NF-kB; bos sensitivity not enhanced"))
f8 <- (p8[[1]] + p8[[2]] + p8[[3]]) / (p8[[4]] + p8[[5]] + p8[[6]]) / (p8[[7]] + p8[[8]] + plot_spacer()) +
  plot_layout(heights = c(0.95, 0.95, 1.25)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(f8, file.path(FIG, "fig8", "F8_combined_v4"), width_mm = 183, height_mm = 150)
cat("DONE 47d\n")
