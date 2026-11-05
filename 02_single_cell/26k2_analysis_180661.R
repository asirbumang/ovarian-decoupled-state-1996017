suppressMessages({library(data.table); library(ggplot2); library(ggpubr); library(patchwork); library(ggsci); library(ragg)})
df <- fread("D:/子宫内膜异位症2/analysis/L1_scRNA/results/module_scores/scores_GSE180661_subset.csv", data.table = FALSE)
cat("cells:", nrow(df), "\n")
cat("cell types:\n"); print(sort(table(df$cell_type)))
cat("\nsupersite:\n"); print(sort(table(df$tumor_supersite)))
cat("\nsite (exact):\n"); print(sort(table(df$tumor_site))[1:30])

canc <- df[grepl("cancer\\.cell", df$cell_type, ignore.case = TRUE), ]
cat("\ncancer cells:", nrow(canc), "\n")
## 位点级（tumor_supersite）
s <- aggregate(decouple ~ tumor_supersite, data = canc, FUN = function(x) c(n = length(x), med = median(x)))
print(s)
## 卵巢 vs 网膜（supersite = Ovary/Omentum?）
for (grp in unique(canc$tumor_supersite)) {
  v <- canc$decouple[canc$tumor_supersite == grp]
  cat(sprintf("%-14s n=%5d med=%.3f\n", grp, length(v), median(v)))
}
w <- tryCatch(wilcox.test(decouple ~ tumor_supersite,
                          data = canc[canc$tumor_supersite %in% c("Ovary", "Omentum"), ]), error = function(e) NULL)
if (!is.null(w)) cat("\nOvary vs Omentum (cancer) p =", format(w$p.value, digits = 4), "\n")
## 病人级（患者中位）——避免伪重复
pw <- aggregate(decouple ~ patient_id, data = canc, FUN = median)
cat("patients with cancer cells:", nrow(pw), "\n")
pw$sup <- aggregator <- NULL
pw_exact <- aggregate(decouple ~ patient_id + tumor_supersite, data = canc, FUN = median)
cat("\npatient x supersite medians (top):\n")
print(head(pw_exact, 20))
## 图：癌细胞的部位分布（快速）
theme_set(theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35), panel.grid = element_blank()))
p1 <- ggplot(canc, aes(tumor_supersite, decouple, fill = tumor_supersite)) + geom_boxplot(outlier.size = .3, width = .5) +
  stat_compare_means(comparisons = list(c("Ovary", "Omentum")), size = 2.2) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cancer-cell decoupling by site (pre-Rx)", x = NULL, y = "decoupling")
out <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/module_scores"
svglite::svglite(file.path(out, "GSE180661_A_site.svg"), width = 90/25.4, height = 70/25.4); print(p1); dev.off()
ragg::agg_tiff(file.path(out, "GSE180661_A_site.tiff"), width = 90/25.4, height = 70/25.4, units = "in", res = 600); print(p1); dev.off()
cat("DONE follow-up\n")
