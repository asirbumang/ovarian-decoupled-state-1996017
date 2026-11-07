## makeFig8_fixed.R -- rebuild Fig 8 with the corrected panels:
## c = corrected evidence card, f = two-line residual labels, g = complete row labels
## and LDH at 24 h, h = two-panel W2 (no clipped conclusion box).
options(warn = 1)
source("analysis/figs_final/_v8/R_scripts/common.R")

FIGS <- "analysis/L1_scRNA/results/figs"
SUP  <- file.path(FIGS, "supp")
AUD  <- "_audit_r3"

p_a <- panel_img(file.path(FIGS, "F8a_km.png")) + labs(tag = "a")
p_b <- panel_img(file.path(FIGS, "F8b_layers.png")) + labs(tag = "b")
p_c <- panel_img(file.path(AUD, "new_fig8c_card.png")) + labs(tag = "c")
p_d <- panel_img(file.path(FIGS, "F11_consensus.png")) + labs(tag = "d")
p_e <- panel_img(file.path(FIGS, "S11_subtype_char.png")) + labs(tag = "e")
p_f <- panel_img(file.path(FIGS, "F8f_residual.png")) + labs(tag = "f")
p_g <- panel_img(file.path(AUD, "new_fig8g_w1_v3.png")) + labs(tag = "g")
p_h <- panel_img(file.path(SUP, "S13_w2.png")) + labs(tag = "h")

row1 <- (p_a | p_b | p_c) + plot_layout(widths = c(0.24, 0.26, 0.50))
row2 <- (p_d | p_e) + plot_layout(widths = c(0.50, 0.50))
row3 <- (p_f | p_g | p_h) + plot_layout(widths = c(0.36, 0.30, 0.34))
fig <- row1 / row2 / row3 +
  plot_layout(heights = c(0.34, 0.27, 0.39)) +
  tag_theme() +
  theme(plot.margin = margin(0, 0, 0, 0),
        plot.background = element_rect(fill = "white", colour = NA))

save_pub_r(fig, "analysis/figs_final/_v9/v9_Fig8_fixed", width_mm = 183, height_mm = 140, dpi = 600)
cat("Fig8 fixed OK\n")
