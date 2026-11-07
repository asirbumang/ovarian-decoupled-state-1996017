## makeFig4.R -- Fig4 v9: (1) EVERY panel white background (the F4a triptych keeps its
## internal dark imaging content but is no longer framed by a black rectangle); (2) panel c
## is REDRAWN from CSV: left = 14-sample Moran's I dot plot, right = 12-sample neighbourhood
## OR bar (red = Fisher p<0.05) -- the old raster F4cd_spatial_stats.png is deleted;
## (3) global compactness: row heights = content + ~0mm, spacing <=1mm, canvas 183 x 208 mm.
options(warn = 1)
source("analysis/figs_final/_v9/R_scripts/common.R")
theme_set(base_theme())

FIGS <- "analysis/L1_scRNA/results/figs"
NEW  <- "analysis/figs_new"
SPAT <- "data_spatial/GSE211956/raw"
SPC  <- "analysis/L0_bulk/results/spatial"

## ---- panel a: dark imaging plate -> WHITE outer background (no black rectangle) ----
p_a <- panel_img(file.path(FIGS, "F4a_triptych.png"), bg = "white") + labs(tag = "a")

## ---- white maps (cropped + enlarged to fill their cell) ----
p_e <- panel_img(file.path(NEW, "SPT1_celltype_maps.png")) + labs(tag = "e")
p_f <- panel_img(file.path(NEW, "SPT2_state_maps.png"))   + labs(tag = "f")
p_b <- panel_img(file.path(FIGS, "F4b_domain_state.png")) + labs(tag = "b")
p_d <- panel_img(file.path(FIGS, "F4de_coloc.png"))       + labs(tag = "d")
p_g <- panel_img(file.path(NEW, "SPT3_domain_composition.png")) + labs(tag = "g")
p_h <- panel_img(file.path(NEW, "SPT4_radial.png"))       + labs(tag = "h")

## ---- panel c: REDRAW from CSV (raster F4cd_spatial_stats.png removed) ----
## short sample names
sn <- c("Patient_16"="P16","Patient_17"="P17","Patient_2"="P2","Patient_8"="P8",
        "Pt1-1"="Pt1-1","Pt1-2"="Pt1-2","Pt1-3"="Pt1-3","Pt1-4"="Pt1-4",
        "Pt2-1"="Pt2-1","Pt2-2"="Pt2-2","Pt2-3"="Pt2-3","Pt2-4"="Pt2-4",
        "OCCR2v"="OCCR2v","OCCS3v"="OCCS3v")
## source of Moran's I: glob *moran* once; prefer the 14-sample summary.
moran_f <- list.files(SPC, pattern = "moran", full.names = TRUE)
moran_f <- moran_f[grepl("summary", moran_f, ignore.case = TRUE)]
if (length(moran_f) == 0) moran_f <- list.files(SPC, pattern = "moran", full.names = TRUE)[1]
m <- read.csv(moran_f[1])
m$short <- sn[m$sample]
m$short[is.na(m$short)] <- m$sample[is.na(m$short)]
m <- m[order(m$Moran_I), ]
m$short <- factor(m$short, levels = m$short)
md <- median(m$Moran_I)
cat(sprintf("[Fig4c] moran n=%d  computed median=%.4f\n", nrow(m), md))
md_line <- 0.160   # per user decision (data median = 0.1596)
p_cm <- ggplot(m, aes(x = Moran_I, y = short)) +
  geom_vline(xintercept = md_line, linetype = 2, linewidth = 0.4, colour = "grey45") +
  geom_point(size = 2.0, colour = "#2C3E50") +
  scale_x_continuous(limits = c(0, 0.60), breaks = seq(0, 0.6, 0.2), expand = c(0.02, 0)) +
  scale_y_discrete(expand = c(0.03, 0.03)) +
  labs(x = "Moran's I", y = NULL, title = "Spatial autocorrelation (KNN)") +
  base_theme() +
  theme(plot.title = element_text(size = 6.5, face = "bold", hjust = 0),
        axis.text.y = element_text(size = 6.0),
        axis.text.x = element_text(size = 6.0),
        axis.title.x = element_text(size = 6.5),
        plot.margin = margin(1, 1, 1, 1)) +
  labs(tag = "c")
## right / OR bar (12 samples; red = Fisher p<0.05)
nb <- read.csv(file.path(SPC, "spatial_neighborhood_summary.csv"))
nb$short <- sn[nb$sample]
nb$sig <- nb$fisher_p < 0.05
nb <- nb[order(nb$neigh_OR), ]
nb$short <- factor(nb$short, levels = nb$short)
p_cr <- ggplot(nb, aes(x = neigh_OR, y = short, fill = sig)) +
  geom_col(width = 0.64) +
  geom_text(aes(x = neigh_OR, label = sprintf("%.2f", neigh_OR)),
            hjust = -0.12, size = 2.12, family = "Arial") +
  scale_fill_manual(values = c("TRUE" = SIG, "FALSE" = GRA), name = NULL,
                    breaks = c(TRUE, FALSE), labels = c("p<0.05", "n.s.")) +
  scale_x_continuous(limits = c(0, 8.2), expand = c(0.02, 0)) +
  scale_y_discrete(expand = c(0.03, 0.03)) +
  labs(x = "neighbourhood OR", y = NULL, title = "State-spot neighbourhood enrichment") +
  base_theme() +
  theme(plot.title = element_text(size = 6.5, face = "bold", hjust = 0),
        axis.text.y = element_text(size = 6.0),
        axis.text.x = element_text(size = 6.0),
        axis.title.x = element_text(size = 6.5),
        legend.position = "right",
        legend.text = element_text(size = 6.0),
        legend.key.size = unit(0.24, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(1, 1, 1, 1))
## composite c = Moran|OR (flattened into the row to avoid patchwork nesting issues)

## ---- panel i : domain PCA (redrawn from domain*.csv) ----
domf <- list.files(SPC, pattern = "_domains\\.csv$", full.names = TRUE)
if (length(domf) > 0) {
  cols <- c("M1_immune_activation","M2_metabolic_GSH",
            "M3_mito_cytoskeleton","M4_lipid_perox_rupture","domain")
  d <- do.call(rbind, lapply(domf, function(f) read.csv(f)[, cols, drop = FALSE]))
  ml <- d[, cols[1:4]]
  pc <- prcomp(ml, scale. = TRUE)
  xy <- data.frame(pc$x[, 1:2], domain = factor(d$domain))
  p_i <- ggplot(xy, aes(x = PC1, y = PC2, colour = domain)) +
    geom_point(size = 0.30, alpha = 0.55) +
    scale_colour_viridis_d(option = "magma", name = "domain") +
    labs(title = "Domain PCA (M1-M4 module scores)",
         caption = sprintf("n=%d spatial spots", nrow(xy))) +
    theme(axis.text = element_text(size = 6.0), axis.title = element_text(size = 6.5),
          plot.title = element_text(size = 6.5, face = "bold", hjust = 0),
          plot.caption = element_text(size = 6.0, hjust = 0),
          legend.position = "right", legend.key.size = unit(0.22, "cm"),
          legend.text = element_text(size = 6.0), legend.title = element_text(size = 6.0),
          plot.margin = margin(1, 1, 1, 1)) +
    labs(tag = "i")
} else {
  p_i <- ggplot() + annotate("label", x = 0.5, y = 0.5, label = "N/A",
                             size = 4) + theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = NA)) +
    labs(tag = "i")
}

## ---- panel j / k : H&E histology + spot overlay (GSE211956) ----
he_panel <- function(gsm, sp, scale, tag = NULL) {
  img <- png::readPNG(file.path(SPAT, gsm, "spatial/tissue_hires_image.png"))
  W <- ncol(img); H <- nrow(img)
  pos <- read.csv(file.path(SPAT, gsm, "spatial/tissue_positions_list.csv"), header = FALSE)
  names(pos) <- c("barcode","in_tissue","array_row","array_col","pxl_row","pxl_col")
  pos <- pos[pos$in_tissue == 1, ]
  pos$x <- pos$pxl_col * scale
  pos$y <- H - pos$pxl_row * scale
  sc <- read.csv(file.path(SPC, "GSE211956_domains.csv"))
  sc <- sc[sc$sample == sp, c("spot","domain")]
  m <- merge(pos, sc, by.x = "barcode", by.y = "spot")
  m <- m[!is.na(m$domain), ]
  m$domain <- factor(m$domain)
  p <- ggplot(m, aes(x = x, y = y, colour = domain)) +
    annotation_custom(grid::rasterGrob(img, interpolate = TRUE),
                      xmin = 0, xmax = W, ymin = 0, ymax = H) +
    geom_point(alpha = 0.8, size = 0.55) +
    scale_colour_viridis_d(option = "magma", name = "domain") +
    coord_fixed(ratio = H / W) +
    scale_x_continuous(limits = c(0, W), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, H), expand = c(0, 0)) +
    theme_void() +
    theme(legend.position = "right",
          legend.text = element_text(size = 6.0),
          legend.title = element_text(size = 6.0),
          legend.key.size = unit(0.2, "cm"),
          plot.title = element_text(size = 6.0, face = "bold", family = "Arial", hjust = 0),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(1, 1, 1, 1)) +
    annotate("label", x = 0.02 * W, y = 0.98 * H, hjust = 0, vjust = 1,
             label = sprintf("GSE211956 %s\nspots=%d", gsm, nrow(m)),
             size = 2.12, family = "Arial", fill = "white", colour = "#222222") +
    labs(title = sprintf("%s (SP%d)", gsm, as.integer(sub("SP", "", sp))))
  if (!is.null(tag)) p <- p + labs(tag = tag)
  p
}
p_j <- he_panel("GSM6506110", "SP1", 0.1906396,  tag = "j")
p_k <- he_panel("GSM6506112", "SP3", 0.17675652, tag = "k")

## ---- assemble: image rows set width ~ aspect so they fill; data row c|i ----
row1 <- (p_a | p_b) + plot_layout(widths = c(0.66, 0.34))
row2 <- (p_e | p_f) + plot_layout(widths = c(0.44, 0.56))
row3 <- (p_d | p_g | p_h) + plot_layout(widths = c(0.39, 0.34, 0.27))
row4 <- (p_cm | p_cr | p_i) + plot_layout(widths = c(0.29, 0.29, 0.42))
row5 <- (p_j | p_k) + plot_layout(widths = c(0.48, 0.52))

fig <- row1 / row2 / row3 / row4 / row5 +
  plot_layout(heights = c(0.50, 0.28, 0.28, 0.38, 0.48)) +
  tag_theme() +
  theme(plot.margin = margin(0, 0, 0, 0),
        plot.background = element_rect(fill = "white", colour = NA))

save_pub_r(fig, "analysis/figs_final/_v9/v9_Fig4", width_mm = 183, height_mm = 192, dpi = 600)
cat("Fig4 OK\n")
cat("median 0.160 per user decision（data 0.1596）\n")
