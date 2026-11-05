## makeFig7.R -- Fig7 v8: 8 panels (a-h); c drug-sens hero (row1 right); g/h R-drawn
## (>=6pt).  Blank spacer rows REMOVED; no blank rows; 183mm x 220mm.
options(warn = 1)
source("analysis/figs_final/_v8/R_scripts/common.R")
theme_set(base_theme())

FIGS <- "analysis/L1_scRNA/results/figs"

p_a <- panel_img(file.path(FIGS, "F7a_b9.png")) + labs(tag = "a")
p_b <- panel_img(file.path(FIGS, "F7b_combo.png")) + labs(tag = "b")
p_c <- panel_img(file.path(FIGS, "F7c_drug.png")) + labs(tag = "c")   ## hero
p_d <- panel_img(file.path(FIGS, "F7d_signature.png")) + labs(tag = "d")
p_e <- panel_img(file.path(FIGS, "F7e_alpha.png")) + labs(tag = "e")
p_f <- panel_img(file.path(FIGS, "F7f_evidence.png")) + labs(tag = "f")

## ---- panel g: permutation p (5 drugs), proper horizontal bars (>=6pt) ----
perm <- data.frame(drug = c("bos","ner","afa","lap","erlo"),
                   p = c(0.0021,0.0011,0.0111,0.0141,0.0221))
perm$drug <- factor(perm$drug, levels = rev(c("bos","ner","afa","lap","erlo")))
p_g <- ggplot(perm, aes(y = drug, x = -log10(p))) +
  geom_col(fill = BLU, width = 0.6) +
  geom_text(aes(label = sprintf("p=%.4f", p), x = -log10(p) + 0.15),
            hjust = 0, size = 2.12, family = "Arial") +
  scale_x_continuous(limits = c(0, 4.4), expand = c(0.02, 0.02)) +
  labs(x = "-log10(p) permutation", y = NULL, title = "Drug permutation") +
  theme(axis.text = element_text(size = 6.0), axis.title = element_text(size = 6.5),
        plot.title = element_text(size = 7, face = "bold", hjust = 0)) +
  labs(tag = "g")

## ---- panel h: stratified effect (>=6pt) ----
strat <- data.frame(grp = c("breast","other","LUNG"),
                    mid = c(-0.44, -0.295, NA),
                    lo  = c(-0.47, -0.34, NA),
                    hi  = c(-0.41, -0.25, NA))
p_h <- ggplot(strat, aes(x = grp, y = mid)) +
  geom_segment(aes(x = grp, xend = grp, y = lo, yend = hi), linewidth = 0.7, colour = SIG) +
  geom_point(size = 1.4, colour = SIG) +
  annotate("text", x = 3, y = 0.02, label = "LUNG NS", hjust = 1,
           size = 2.12, family = "Arial") +
  annotate("text", x = 1.2, y = -0.30, label = "-0.41 to -0.47", hjust = 0,
           size = 2.12, family = "Arial") +
  annotate("text", x = 2.2, y = -0.30, label = "-0.25 to -0.34", hjust = 0,
           size = 2.12, family = "Arial") +
  scale_y_continuous(limits = c(-0.62, 0.12), breaks = c(-0.6, -0.4, -0.2, 0.0),
                     labels = c("-0.6", "-0.4", "-0.2", "0.0")) +
  labs(x = NULL, y = "effect size", title = "Stratified effect") +
  theme(axis.text.x = element_text(size = 6.0)) +
  labs(tag = "h")

row1 <- (p_a | p_b | p_c) + plot_layout(widths = c(0.31, 0.31, 0.38))
row2 <- (p_d | p_e | p_f) + plot_layout(widths = c(0.34, 0.33, 0.33))
row3 <- (p_g | p_h) + plot_layout(widths = c(0.5, 0.5))
fig <- row1 / row2 / row3 +
  plot_layout(heights = c(0.42, 0.28, 0.30)) +   # content-proportional (57/39/35mm)
  tag_theme() +
  theme(plot.margin = margin(0, 0, 0, 0),
        plot.background = element_rect(fill = "white", colour = NA))

save_pub_r(fig, "analysis/figs_final/_v9/v9_Fig7", width_mm = 183, height_mm = 148, dpi = 600)
cat("Fig7 OK\n")
