## makeFig5.R -- Fig5 v8: 7 panels (a-g).  READING-ORDER FIX: previously c (network
## hero) was LEFT of a (commot) -- violates "no c-left/a-right".  Now a is top-left,
## c is the right/centre hero.  REDUNDANCY FIX: b and h were the SAME image
## (F5d_lr_count) -- h removed (7 panels).  c = pure-graphical comm network (no param
## table); no blank rows; 183mm x 220mm.
options(warn = 1)
source("analysis/figs_final/_v8/R_scripts/common.R")
theme_set(base_theme())

FIGS <- "analysis/L1_scRNA/results/figs"
NEW  <- "analysis/figs_new"

p_a <- panel_img(file.path(FIGS, "F5a_commot.png")) + labs(tag = "a")
p_b <- panel_img(file.path(FIGS, "F5d_lr_count.png")) + labs(tag = "b")
p_c <- panel_img(file.path(NEW, "SPT5_comm_network.png")) + labs(tag = "c")  ## hero, graphical
p_d <- panel_img(file.path(FIGS, "F5e_comp_state.png")) + labs(tag = "d")
p_e <- panel_img(file.path(FIGS, "F5b_nnls_delta.png")) + labs(tag = "e")
p_f <- panel_img(file.path(FIGS, "F5c_geomx.png")) + labs(tag = "f")

## ---- panel g: LR significance (7 pairs), proper horizontal bars (>=6pt) ----
lr_g <- data.frame(
  pair = paste("LR", 1:7),
  p = c(0.002, 0.001, 0.004, 0.001, 0.023, 0.018, 0.027))
lr_g$pair <- factor(lr_g$pair, levels = rev(paste("LR", 1:7)))
p_g <- ggplot(lr_g, aes(y = pair, x = -log10(p))) +
  geom_col(fill = SIG, width = 0.62) +
  geom_text(aes(label = sprintf("p=%.3f", p), x = -log10(p) + 0.15),
            hjust = 0, size = 2.12, family = "Arial") +
  scale_x_continuous(limits = c(0, 4.1), expand = c(0.02, 0.02)) +
  labs(x = "-log10(p) [COMMOT permutation]", y = NULL,
       title = "LR pair significance (7 pairs)") +
  theme(axis.text = element_text(size = 6.0), axis.title = element_text(size = 6.5),
        plot.title = element_text(size = 7, face = "bold", hjust = 0)) +
  labs(tag = "g")

## ---- assemble: a top-left, c (network) hero to its RIGHT; no blank rows ----
row1 <- (p_a | p_c) + plot_layout(widths = c(0.40, 0.60))   # a then c (hero) - reading ok
row2 <- (p_b | p_d | p_e) + plot_layout(widths = c(0.34, 0.33, 0.33))
row3 <- (p_f | p_g) + plot_layout(widths = c(0.50, 0.50))
fig <- row1 / row2 / row3 +
  plot_layout(heights = c(0.42, 0.30, 0.28)) +   # content-proportional (67/45/35mm)
  tag_theme() +
  theme(plot.margin = margin(0, 0, 0, 0),
        plot.background = element_rect(fill = "white", colour = NA))

save_pub_r(fig, "analysis/figs_final/_v9/v9_Fig5", width_mm = 183, height_mm = 158, dpi = 600)
cat("Fig5 OK\n")
