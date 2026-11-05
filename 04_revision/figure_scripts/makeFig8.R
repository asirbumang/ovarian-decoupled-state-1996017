## makeFig8.R -- Fig8 v8: 8 panels (a-h). c honest-conclusions graphical bar (kept,
## no long text); d/e wide images; f residual hero; g W1 (S12) & h W2 (S13) -- these
## carry the content of the DELETED Fig9 (main-text refs point to Fig8g/h).  Blank
## spacer rows REMOVED; 183mm x 210mm.
options(warn = 1)
source("analysis/figs_final/_v8/R_scripts/common.R")

FIGS <- "analysis/L1_scRNA/results/figs"
SUP  <- file.path(FIGS, "supp")

p_a <- panel_img(file.path(FIGS, "F8a_km.png")) + labs(tag = "a")
p_b <- panel_img(file.path(FIGS, "F8b_layers.png")) + labs(tag = "b")
p_d <- panel_img(file.path(FIGS, "F11_consensus.png")) + labs(tag = "d")   ## wide
p_e <- panel_img(file.path(FIGS, "S11_subtype_char.png")) + labs(tag = "e") ## wide
p_f <- panel_img(file.path(FIGS, "F8f_residual.png")) + labs(tag = "f")    ## hero
p_g <- panel_img(file.path(SUP, "S12_w1.png")) + labs(tag = "g")
p_h <- panel_img(file.path(SUP, "S13_w2.png")) + labs(tag = "h")

## ---------------- panel c: honest conclusions -> graphical summary (kept) ----------------
rob <- data.frame(
  lab = c("Carrier-enriched", "Spatially clustered", "TF-controlled", "Drug-axis (CTRPv2)"),
  val = c(1.0, 1.0, 1.0, 1.0), stringsAsFactors = FALSE)
not <- data.frame(
  lab = c("External I\u00b2 = 0%", "Response NS", "ML no generalisation"),
  val = c(0.50, 0.50, 0.50), stringsAsFactors = FALSE)
rob$y <- c(6, 7, 8, 9)
not$y <- c(1, 2, 3)
hc  <- rbind(transform(rob, col = GRN), transform(not, col = SIG))
band <- data.frame(ymin = c(5.45, 0.45), ymax = c(9.6, 4.55),
                   fill = c("#E6F3EF", "#FBE9E4"), stringsAsFactors = FALSE)
p_c <- ggplot(hc, aes(y = y)) +
  geom_rect(data = band, inherit.aes = FALSE,
            aes(xmin = -0.18, xmax = 2.75, ymin = ymin, ymax = ymax),
            fill = band$fill, colour = NA) +
  geom_segment(aes(x = -0.08, xend = val, y = y, yend = y, colour = col),
               linewidth = 3.4, lineend = "butt") +
  geom_point(aes(x = val, y = y, colour = col), size = 1.2) +
  geom_text(aes(x = val + 0.06, y = y, label = lab, colour = col),
            hjust = 0, vjust = 0.5, size = 2.12, family = "Arial") +
  scale_colour_identity() +
  scale_x_continuous(limits = c(-0.20, 3.2), breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.2, 10.0), breaks = NULL, expand = c(0, 0)) +
  annotate("text", x = -0.20, y = 9.85, label = "ROBUST / reproducible",
           hjust = 0, size = 2.28, fontface = "bold", colour = GRN, family = "Arial") +
  annotate("text", x = -0.20, y = 0.60, label = "NOT REPRODUCIBLE",
           hjust = 0, size = 2.28, fontface = "bold", colour = SIG, family = "Arial") +
  annotate("text", x = 3.15, y = 9.85, label = "see text / Table 7",
           hjust = 1, size = 2.12, family = "Arial", colour = "#666666") +
  theme_void() +
  theme(plot.margin = margin(3, 3, 3, 3),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA)) +
  labs(tag = "c")

row1 <- (p_a | p_b | p_c) + plot_layout(widths = c(0.24, 0.26, 0.50))
row2 <- (p_d | p_e) + plot_layout(widths = c(0.50, 0.50))
row3 <- (p_f | p_g | p_h) + plot_layout(widths = c(0.36, 0.30, 0.34))
fig <- row1 / row2 / row3 +
  plot_layout(heights = c(0.34, 0.27, 0.39)) +   # content-proportional (35/30/43mm)
  tag_theme() +
  theme(plot.margin = margin(0, 0, 0, 0),
        plot.background = element_rect(fill = "white", colour = NA))

save_pub_r(fig, "analysis/figs_final/_v9/v9_Fig8", width_mm = 183, height_mm = 140, dpi = 600)
cat("Fig8 OK\n")
