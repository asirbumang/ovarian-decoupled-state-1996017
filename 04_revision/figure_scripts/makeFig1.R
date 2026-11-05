## makeFig1.R -- Fig1 v8: 7 panels (a-g). No text-dominated panel; a = annotation
## strip (full-width, top = reading start), then b|c, d|e, f|g.  BLANK spacer rows
## REMOVED (no empty band); panel gaps ~1mm (<=1.5mm); text >=6pt; tags bold 7.5pt;
## content-proportional layout; canvas 183mm x 210mm.
options(warn = 1)
source("analysis/figs_final/_v8/R_scripts/common.R")

FIGS <- "analysis/L1_scRNA/results/figs"

## ---------------- panel a: annotation strip ----------------
coh <- data.frame(
  ds    = c("GSE266577","GSE222556","GSE201047","GSE180661*","GSE293459#","GSE156699",
            "GSE26193","GSE13876","GSE14764","GSE63885","GSE26712","TCGA-OV",
            "GSE274657","GSE288483","GSE224335","GSE211956","GSE189843","GSE276935"),
  layer = c("scRNA","scRNA","scRNA","scRNA","scRNA","bulk","bulk","bulk","bulk","bulk",
            "bulk","bulk","Visium","Visium","Visium","Visium","Visium","GeoMx"),
  role  = c("State/cancer","Metastasis","Chemo-naive","Site","Persister","Response","OS","OS",
            "OS","OS+Platinum","Tumour/Normal","OS+Response","Spatial","Spatial","Spatial",
            "Spatial","Spatial","TLS"),
  tech  = c("10x 3' v3.1","10x 3' v3.1","10x v2","10x sorted","snRNA","RNA-seq","GPL570",
            "GPL96","GPL96","GPL570","GPL570","RNA-seq","10x Visium","10x Visium",
            "10x Visium","10x Visium","10x Visium","GeoMx DSP"),
  stringsAsFactors = FALSE)

layer_col <- c(scRNA = "#E64B35", bulk = "#4DBBD5", Visium = "#00A087", GeoMx = "#3C5488")
role_col  <- c("State/cancer"="#E64B35","Metastasis"="#4DBBD5","Chemo-naive"="#00A087",
               "Site"="#3C5488","Persister"="#EE4C97","Response"="#F39B7F","OS"="#91D1C2",
               "OS+Platinum"="#8491B4","Tumour/Normal"="#DC0000","OS+Response"="#FFDC91",
               "Spatial"="#B8B8B8","TLS"="#7E6148")
tech_col  <- c("10x 3' v3.1"="#FFF5F0","10x v2"="#FEE0D2","10x sorted"="#C6DBEF",
               "snRNA"="#9ECAE1","RNA-seq"="#E5F5E0","GPL570"="#A1D99B","GPL96"="#74C476",
               "10x Visium"="#FFE4C4","GeoMx DSP"="#807DBA")

long <- data.frame(
  ds  = rep(coh$ds, 3),
  attr = rep(c("layer","role","tech"), each = nrow(coh)),
  val = c(as.character(coh$layer), as.character(coh$role), as.character(coh$tech)),
  stringsAsFactors = FALSE)
long$fill <- mapply(function(a, v) {
  if (a == "layer") layer_col[[v]] else if (a == "role") role_col[[v]] else tech_col[[v]]
}, long$attr, long$val)
long$attr <- factor(long$attr, levels = c("layer","role","tech"))
long$ds   <- factor(long$ds,   levels = coh$ds)

a_strip <- ggplot(long, aes(x = ds, y = attr)) +
  geom_tile(aes(fill = fill), colour = "white", linewidth = 0.4) +
  scale_fill_identity() +
  scale_x_discrete(position = "top") +
  scale_y_discrete(expand = c(0, 0)) +
  theme_void() +
  theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 1,
                                       size = 6.0, family = "Arial"),
        axis.text.y = element_text(size = 6.0, hjust = 1, family = "Arial"),
        plot.margin = margin(2, 2, 0, 0)) +
  labs(tag = "a")

## ---- legend: full-width, four columns (role split), 5.5pt ----
leg_groups <- list(
  layer = list(labs = names(layer_col), cols = unname(layer_col)),
  role1 = list(labs = names(role_col)[1:6], cols = unname(role_col)[1:6]),
  role2 = list(labs = names(role_col)[7:12], cols = unname(role_col)[7:12]),
  tech  = list(labs = names(tech_col), cols = unname(tech_col)))
x0 <- c(layer = 0.0, role1 = 2.7, role2 = 5.4, tech = 8.4)
leg_df <- data.frame()
for (s in names(leg_groups)) {
  labs <- leg_groups[[s]]$labs; cols <- leg_groups[[s]]$cols
  leg_df <- rbind(leg_df, data.frame(x = x0[[s]], y = 9, lab = ifelse(s %in% c("role1","role2"), "role", s),
                                     col = NA, header = TRUE, stringsAsFactors = FALSE))
  ys <- seq(8, by = -1, length.out = length(labs))
  leg_df <- rbind(leg_df, data.frame(x = x0[[s]], y = ys, lab = labs, col = cols,
                                     header = FALSE, stringsAsFactors = FALSE))
}
a_leg <- ggplot(leg_df) +
  geom_text(data = leg_df[leg_df$header, ], aes(x = x, y = y, label = lab),
            hjust = 0, size = 1.95, fontface = "bold", family = "Arial") +
  geom_tile(data = leg_df[!leg_df$header, ], aes(x = x + 0.14, y = y, fill = col),
            width = 0.26, height = 0.82, colour = NA) +
  geom_text(data = leg_df[!leg_df$header, ], aes(x = x + 0.34, y = y, label = lab),
            hjust = 0, size = 1.85, family = "Arial") +
  scale_fill_identity() +
  scale_x_continuous(limits = c(-0.1, 11.6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.3, 9.7), expand = c(0, 0)) +
  theme_void() + theme(plot.margin = margin(0, 1, 0, 1, "mm"))

## ---------------- panel b / c / f : images ----------------
p_b <- panel_img(file.path(FIGS, "F1b_umap.png")) + labs(tag = "b")
p_c <- panel_img(file.path(FIGS, "F1c_dotplot.png")) + labs(tag = "c")
p_f <- panel_img(file.path(FIGS, "F1f_qc.png")) + labs(tag = "f")

## ---------------- panel d: Hedges' g forest ----------------
hedge <- data.frame(lab = c("M3 mito-cyto", "M2 metab/GSH", "M1 immune"),
                    g = c(0.66, 1.79, -0.14), lo = c(0.02, 1.13, -0.78),
                    hi = c(1.30, 2.45, 0.50), stringsAsFactors = FALSE)
hedge$y <- rev(seq_len(nrow(hedge)))
p_d <- ggplot(hedge, aes(y = y)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = "grey55") +
  geom_segment(aes(x = lo, xend = hi, y = y, yend = y), linewidth = 0.5, colour = SIG) +
  geom_point(aes(x = g), size = 1.3, colour = SIG) +
  geom_text(aes(x = hi + 0.12, label = sprintf("%.2f [%.2f-%.2f]", g, lo, hi)),
            size = 2.12, hjust = 0, family = "Arial") +
  scale_y_continuous(breaks = hedge$y, labels = hedge$lab, limits = c(0.4, 3.6)) +
  scale_x_continuous(limits = c(-1.1, 3.6),
                     breaks = c(-1, 0, 1, 2), labels = c("-1", "0", "1", "2")) +
  labs(x = "Hedges' g [95% CI]", y = NULL,
       title = "Module contrast (transformation axis)") +
  annotate("text", x = -1.0, y = 3.45, label = "GSE26712 tumour vs normal OSE",
           size = 2.12, hjust = 0, fontface = "bold", family = "Arial") +
  annotate("text", x = -1.0, y = 3.25, label = "Welch, 95% CI; n=185 vs 10",
           size = 2.12, hjust = 0, family = "Arial", colour = "#555555") +
  base_theme() + labs(tag = "d")

## ---------------- panel e: discriminator ----------------
disc <- data.frame(lab = c("Pyroptosis","Necroptosis","Apoptosis","Ferroptosis"),
                   v = c(-0.03, -0.11, 0.27, 0.23), stringsAsFactors = FALSE)
disc$y <- rev(seq_len(nrow(disc)))
p_e <- ggplot(disc, aes(y = y)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey55") +
  geom_segment(aes(x = 0, xend = v, y = y, yend = y), linewidth = 0.9, colour = BLU) +
  geom_point(aes(x = v, y = y), size = 1.2, colour = BLU) +
  geom_text(aes(x = v, label = sprintf("%+.2f", v), hjust = ifelse(v >= 0, -0.15, 1.15)),
            size = 2.12, family = "Arial") +
  scale_y_continuous(breaks = disc$y, labels = disc$lab, limits = c(0.4, 4.8)) +
  scale_x_continuous(limits = c(-0.20, 0.36), breaks = c(-0.1, 0, 0.1, 0.2)) +
  labs(x = "partial r", y = NULL, title = "Decouple vs death programmes") +
  annotate("text", x = -0.19, y = 4.65, label = "composition-adjusted",
           size = 2.12, hjust = 0, family = "Arial", colour = "#444444") +
  base_theme() + labs(tag = "e")

## ---------------- panel g: sample-size log bar ----------------
samp <- data.frame(
  ds  = c("GSE266577","GSE222556","GSE201047","GSE180661*","GSE293459#","GSE26193",
          "GSE14764","GSE13876","GSE63885","GSE26712","TCGA-OV","GSE274657","GSE288483",
          "GSE224335","GSE189843"),
  n   = c(137084,661491,106001,92000,62286,107,80,157,75,195,427,13629,15322,2410,4801),
  lay = c("scRNA","scRNA","scRNA","scRNA","scRNA","bulk","bulk","bulk","bulk","bulk",
          "bulk","Vis","Vis","Vis","Vis"), stringsAsFactors = FALSE)
samp$y <- rev(seq_len(nrow(samp)))
samp$lx <- samp$n * 1.2
lay_all <- c(scRNA = "#E64B35", bulk = "#4DBBD5", Vis = "#00A087")
p_g <- ggplot(samp, aes(y = y)) +
  geom_segment(aes(x = 40, xend = n, colour = lay), linewidth = 0.9) +
  geom_point(aes(x = n, colour = lay), size = 1.1) +
  scale_colour_manual(values = lay_all, name = NULL) +
  scale_x_log10(limits = c(40, 3e6),
                breaks = c(1e2,1e3,1e4,1e5,1e6),
                labels = expression(10^2, 10^3, 10^4, 10^5, 10^6)) +
  scale_y_continuous(breaks = samp$y, labels = samp$ds, limits = c(0.4, 15.6)) +
  geom_text(aes(x = lx, label = format(n, big.mark = ",")),
            size = 2.12, hjust = 0, family = "Arial") +
  labs(x = "cells / samples (log10)", y = NULL, title = "Cohort size (log scale)") +
  theme(legend.position = "bottom", legend.text = element_text(size = 6.0),
        legend.key.size = unit(0.25, "cm")) +
  base_theme() + labs(tag = "g")

## ---------------- assemble (7 panels, content-proportional, NO blank spacer) ----------------
row_bc  <- (p_b | p_c) + plot_layout(widths = c(0.43, 0.57))
row_de  <- (p_d | p_e) + plot_layout(widths = c(0.50, 0.50))
row_fg  <- (p_f | p_g) + plot_layout(widths = c(0.58, 0.42))

fig <- (a_strip / a_leg) / row_bc / row_de / row_fg +
  plot_layout(heights = c(1.35, 1.45, 2.6, 2.0, 2.4)) +
  tag_theme() +
  theme(plot.margin = margin(0, 0, 0, 0),
        plot.background = element_rect(fill = "white", colour = NA))

save_pub_r(fig, "analysis/figs_final/_v9/v9_Fig1", width_mm = 183, height_mm = 225, dpi = 600)
cat("Fig1 OK\n")
