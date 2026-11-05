
suppressMessages({library(ggplot2); library(patchwork); library(data.table); library(png); library(ragg)})
theme_set(theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        legend.title = element_text(size = 6.2), legend.text = element_text(size = 5.8),
        plot.title = element_text(size = 7, face = "bold"), panel.grid = element_blank()))
save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 110, dpi = 600) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi); print(plot); dev.off()
  ggsave(paste0(filename, ".png"), plot, width = w, height = h, dpi = 600)
}
FIG <- "D:/子宫内膜异位症2/analysis/figs_final/fig2"
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
img_a <- readPNG("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F2a_state_umap.png")
img_c <- readPNG("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F2c_quadrant.png")
pa <- ggplot() + annotation_raster(img_a, xmin = 0, xmax = 1, ymin = 0, ymax = 1) + theme_void()
pc <- ggplot() + annotation_raster(img_c, xmin = 0, xmax = 1, ymin = 0, ymax = 1) + theme_void() +
  labs(caption = "Definition (within source): M1z ≥ P75 & M2z ≤ P25 & M3z ≤ P25; z within-source")
## d: 稀有度（文本条 + 二项 95%CI）
rr <- fread("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F2d_rarity.csv")
rr$src <- factor(rr$source, levels = rr$source[order(rr$pct)])
pd_ <- ggplot(rr, aes(x = pct, y = src)) +
  geom_point(size = 1.6, colour = "#3C5488") +
  geom_text(aes(label = sprintf("%s  (n=%s)  %.2f%%", source, format(n, big.mark = ","), pct)),
            hjust = -0.12, size = 2.1) +
  scale_x_continuous(limits = c(0, max(rr$pct) * 2.3)) +
  labs(x = "state frequency (% of cells within source)", y = NULL) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 5.5))
layout <- (pa + pc + pd_) +
  plot_layout(widths = c(1.05, 1.0, 1.25)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(layout, file.path(FIG, "F2_combined_v1"), width_mm = 183, height_mm = 74)
cat("saved F2_combined_v1\nDONE 44d\n")
