
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
pa <- mk(paste0(Q, "F2a_state_umap.png"), "a UMAP state")
pb <- mk(paste0(Q, "F2c_quadrant.png"), "b definition")
pc <- mk(paste0(Q, "F2d_rarity.png"), "c rarity (5 sources)")
pd_ <- mk(paste0(Q, "F2d_or.png"), "d lineage OR")
pe_ <- mk(paste0(Q, "F2e_modules.png"), "e modules (state vs rest)")
pf <- mk(paste0(Q, "F2f_comp.png"), "f composition")
## g 双联：占比（F2g_carriers）+ 模块结构（F2g2）
g1 <- rd(paste0(Q, "F2g_carriers.png"))
g2 <- rd(paste0(Q, "F2g2_module_structure.png"))
pg <- (ggplot() + annotation_raster(g1, 0, 1, 0, 1) + theme_void() +
       labs(caption = "g carrier subpopulations (state fraction)") + theme(plot.caption = element_text(size = 5.4, hjust = 0)) +
       ggplot() + annotation_raster(g2, 0, 1, 0, 1) + theme_void() +
       labs(caption = "module structure within carriers: IGHG+ plasma M1-driven; FOS/CTGF+ CAF M2/M3-driven (IAMD two arms)") +
       theme(plot.caption = element_text(size = 5.4, hjust = 0))) +
  plot_layout(ncol = 2) + plot_annotation(tag_levels = "a")
layout <- (pa + pb + pc) / (pd_ + pe_ + pf) / (pg + plot_spacer()) +
  plot_layout(heights = c(1, 1, 0.9)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
FIG <- "D:/子宫内膜异位症2/analysis/figs_final/fig2"
save_pub_r(layout, file.path(FIG, "F2_combined_v6"), width_mm = 183, height_mm = 128)
cat("saved F2_combined_v6\nDONE 47i\n")
