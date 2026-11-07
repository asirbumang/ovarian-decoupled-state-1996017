## makeFig6.R -- Fig6 v8: 8 panels (a-h). panel a = 15-result direction summary
## graphic (no text wall); d MIME heatmap hero; b/c images; e/f/g/h R-drawn (>=6pt).
## Blank spacer rows REMOVED; no blank rows; 183mm x 230mm.
options(warn = 1)
source("analysis/figs_final/_v8/R_scripts/common.R")

FIGS <- "analysis/L1_scRNA/results/figs"
NEW  <- "analysis/figs_new"

## ---------------- panel a: 15-result direction summary graphic ----------------
summ <- data.frame(
  grp  = c("consistent","opposite","NS"),
  n    = c(7, 5, 3),
  col  = c(BLU, "#DC0000", GRA),
  y    = c(3, 2, 1),
  stringsAsFactors = FALSE)
p_a <- ggplot(summ, aes(y = y)) +
  geom_segment(aes(x = 0, xend = n, y = y, yend = y, colour = grp),
               linewidth = 4.6, lineend = "butt") +
  geom_text(aes(x = n + 0.20, y = y, label = paste0(n, " / 15"), colour = grp),
            hjust = 0, size = 2.28, fontface = "bold", family = "Arial") +
  scale_colour_manual(values = c(consistent = BLU, opposite = "#DC0000", NS = GRA),
                      name = NULL, breaks = c("consistent","opposite","NS"),
                      labels = c("expected","opposite","no difference")) +
  scale_x_continuous(limits = c(0, 9.2), breaks = c(0, 2, 4, 6, 8), expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = c(0.4, 3.8), breaks = NULL) +
  labs(x = NULL, y = NULL, title = "15 results: direction summary") +
  annotate("text", x = 0.6, y = 3.55, hjust = 0, size = 2.28, colour = "#333333",
           label = "OS/Resp/Space/Comm/ML/CNV all canonical", family = "Arial") +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.title = element_text(size = 7, face = "bold", hjust = 0),
        legend.position = "bottom", legend.text = element_text(size = 6.0),
        legend.key.size = unit(0.26, "cm"), legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(2, 2, 2, 2)) +
  labs(tag = "a")

p_b <- panel_img(file.path(FIGS, "F6b_forests.png")) + labs(tag = "b")
p_c <- panel_img(file.path(FIGS, "F6c_crude_resid.png")) + labs(tag = "c")
p_d <- panel_img(file.path(NEW, "F6d2_mime101_heatmap.png")) + labs(tag = "d")  ## hero
p_e <- panel_img(file.path(NEW, "F6d2_mime101_trainval_scatter.png")) + labs(tag = "e")

## ---- panel f: zoo bar (4 values) ----
zoo <- data.frame(lab = c("Z1","Z2","Z3","Z4"), v = c(0.978, 0.896, 0.857, 0.335))
zoo$lab <- factor(zoo$lab, levels = zoo$lab[order(-zoo$v)])
p_f <- ggplot(zoo, aes(x = lab, y = v)) +
  geom_col(fill = BLU, width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", v)), vjust = -0.4, size = 2.12, family = "Arial") +
  scale_y_continuous(limits = c(0, 1.12), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = NULL, y = "score", title = "Zoo (4 models)") +
  theme(axis.text.x = element_text(size = 6.0)) +
  base_theme() + labs(tag = "f")

## ---- panel g: Scissor HR forest ----
p_g <- ggplot(data.frame()) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.35, colour = "grey55") +
  geom_segment(aes(x = 0.72, xend = 1.05, y = 1, yend = 1), linewidth = 0.8, colour = SIG) +
  geom_point(aes(x = 0.87, y = 1), size = 1.4, colour = SIG) +
  annotate("text", x = 1.07, y = 1, hjust = 0, size = 2.12, family = "Arial",
           label = "0.87 [0.72-1.05]\np = 0.17") +
  scale_x_continuous(limits = c(0.4, 1.70), breaks = c(0.5, 1.0, 1.5)) +
  scale_y_continuous(limits = c(0.55, 1.45)) +
  labs(x = "Odds ratio", y = NULL, title = "Scissor") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  base_theme() + labs(tag = "g")

## ---- panel h: 45c bar + range ----
d45 <- data.frame(lab = c("A", "B"), v = c(0.591, 0.502))
p_h <- ggplot(d45, aes(x = lab, y = v)) +
  geom_col(fill = GRA, width = 0.55) +
  geom_text(aes(label = sprintf("%.3f", v)), vjust = -0.4, size = 2.12, family = "Arial") +
  annotate("text", x = 2.55, y = 0.60, hjust = 1, label = "range 0.40-0.47",
           size = 2.12, family = "Arial") +
  scale_y_continuous(limits = c(0, 0.72), breaks = c(0, 0.2, 0.4, 0.6)) +
  labs(x = NULL, y = "score", title = "45c") +
  base_theme() + labs(tag = "h")

top <- (p_a | p_b | p_c) + plot_layout(widths = c(0.36, 0.38, 0.26))
mid <- (p_d | p_e) + plot_layout(widths = c(0.68, 0.32))
bot <- (p_f | p_g | p_h) + plot_layout(widths = c(0.34, 0.33, 0.33))
fig <- top / mid / bot + plot_layout(heights = c(0.32, 0.42, 0.26)) + tag_theme() +
  theme(plot.margin = margin(0, 0, 0, 0),
        plot.background = element_rect(fill = "white", colour = NA))

save_pub_r(fig, "analysis/figs_final/_v9/v9_Fig6", width_mm = 183, height_mm = 162, dpi = 600)
cat("Fig6 OK\n")
