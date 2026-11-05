
suppressMessages({library(ggplot2); library(patchwork); library(data.table); library(png); library(ragg)})
theme_set(theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        legend.title = element_text(size = 6.2), legend.text = element_text(size = 5.8),
        plot.title = element_text(size = 7, face = "bold"), panel.grid = element_blank()))
save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 90, dpi = 600) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi); print(plot); dev.off()
  ggsave(paste0(filename, ".png"), plot, width = w, height = h, dpi = 600)
}
meta <- data.table(
  set = c("GSE266577","GSE222556","GSE201047","GSE180661*","GSE293459#","GSE156699","GSE26193","GSE13876","GSE14764","GSE63885","GSE26712","TCGA-OV","GSE274657","GSE288483","GSE224335","GSE211956","GSE189843","GSE276935"),
  layer = c(rep("scRNA",5), rep("bulk",7), rep("Visium",5), "GeoMx"),
  tech = c("[10x v3.1]","[10x v3.1]","[10x v2]","[10x sorted]","[snRNA]","[RNA-seq]","[GPL570]","[GPL96]","[GPL96]","[GPL570]","[GPL570]","[RNA-seq]","[Visium]","[Visium]","[Visium]","[Visium]","[Visium]","[GeoMx]"))
meta$set <- factor(meta$set, levels = meta$set)
long <- rbind(data.table(set = meta$set, ann = "layer", val = meta$layer),
              data.table(set = meta$set, ann = "tech", val = meta$tech))
long$val <- factor(long$val, levels = unique(long$val))
pa <- ggplot(long, aes(x = set, y = ann, fill = val)) + geom_tile(width = 1, height = 0.92, colour = "white", linewidth = 0.2) +
  scale_x_discrete(position = "top") +
  scale_fill_manual(values = c(scRNA="#E64B35", bulk="#4DBBD5", Visium="#00A087", GeoMx="#3C5488",
    "[10x v3.1]"="#FFF5F0","[10x v2]"="#FEE0D2","[10x sorted]"="#C6DBEF","[snRNA]"="#9ECAE1","[RNA-seq]"="#E5F5E0","[GPL570]"="#A1D99B","[GPL96]"="#74C476","[Visium]"="#BCBDDC","[GeoMx]"="#807DBA")) +
  theme_minimal(base_size = 6.5) + theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 5.4),
        axis.text.y = element_text(size = 6), panel.grid = element_blank(), legend.position = "right",
        legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 5.4), legend.title = element_blank()) +
  labs(x = NULL, y = NULL)
gd <- fread("D:/子宫内膜异位症2/analysis/L0_bulk/results/survival_meta/module_effects_hedgesg.csv")
gd <- gd[cohort == "GSE26712" & module != "M4_lipid_perox_rupture"]
lab_map <- c(M1_immune_activation = "M1 immune", M2_metabolic_GSH = "M2 metab/GSH", M3_mito_cytoskeleton = "M3 mito-cyto")
gd$m2 <- lab_map[gd$module]
pd_ <- ggplot(gd, aes(x = g, y = factor(m2, levels = c("M3 mito-cyto", "M2 metab/GSH", "M1 immune")))) +
  geom_point(size = 1.4, colour = "#E64B35") +
  geom_segment(aes(x = ci_lo, xend = ci_hi, y = y, yend = y),
               data = transform(gd, y = as.numeric(factor(m2, levels = c("M3 mito-cyto", "M2 metab/GSH", "M1 immune")))), linewidth = 0.4, colour = "#E64B35") +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35) +
  annotate("text", x = -0.9, y = 3.75, label = "Tumor vs normal OSE (n=185 vs 10)", size = 5.0, fontface = "bold", hjust = 0, vjust = 1) +
  coord_cartesian(xlim = c(-1.0, 2.8), ylim = c(0.6, 3.9)) + labs(x = "Hedges' g [95% CI] (Welch)", y = NULL) +
  theme(axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 6))
disc <- data.frame(x = factor(c("Pyroptosis","Necroptosis","Apoptosis","Ferroptosis"), levels = c("Pyroptosis","Necroptosis","Apoptosis","Ferroptosis")), r = c(-0.03,-0.11,0.27,0.23))
pe_ <- ggplot(disc, aes(x = x, y = r)) + geom_col(width = 0.6, fill = "#3C5488") + geom_hline(yintercept = 0, linewidth = 0.35) +
  labs(x = NULL, y = "partial r", title = "Partial correlation with death-program scores") +
  theme(axis.text.x = element_text(size = 5.4, angle = 30, hjust = 1), axis.text.y = element_text(size = 5.4))
img_b <- readPNG("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F1b_umap.png")
img_c <- readPNG("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F1c_dotplot.png")
img_f <- readPNG("D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs/F1f_qc.png")
pb <- ggplot() + annotation_raster(img_b, 0, 1, 0, 1) + theme_void()
pc <- ggplot() + annotation_raster(img_c, 0, 1, 0, 1) + theme_void()
pf <- ggplot() + annotation_raster(img_f, 0, 1, 0, 1) + theme_void()
layout <- (pa / (pb + pc) / (pd_ + pe_ + pf)) + plot_layout(heights = c(0.9, 1.5, 1.0)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 8, face = "bold"))
FIG <- "D:/子宫内膜异位症2/analysis/figs_final/fig1"
save_pub_r(layout, file.path(FIG, "F1_combined_v2"), width_mm = 183, height_mm = 135)
cat("saved F1_combined_v2\nDONE 47f\n")
