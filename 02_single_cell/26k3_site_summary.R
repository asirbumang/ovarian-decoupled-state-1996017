suppressMessages({library(data.table); library(ggplot2); library(ggpubr); library(ggsci); library(ragg)})
df <- fread("D:/子宫内膜异位症2/analysis/L1_scRNA/results/module_scores/scores_GSE180661_subset.csv", data.table = FALSE)
cat("cell types:\n"); print(sort(table(df$cell_type)))
canc <- df[grepl("cancer\\.cell", df$cell_type, ignore.case = TRUE), ]
cat("\ncancer:", nrow(canc), "\n")
cat("supersites (cancer):\n"); print(table(canc$tumor_supersite))
cat("\nmedians by supersite (cancer):\n")
for (grp in unique(canc$tumor_supersite)) {
  v <- canc$decouple[canc$tumor_supersite == grp]
  cat(sprintf("%-16s n=%5d med=%.3f\n", grp, length(v), median(v)))
}
## 患者级中位（supersite × patient）
pw <- aggregate(decouple ~ patient_id + tumor_supersite, data = canc, FUN = median)
cat("\npatient-level by supersite:\n")
for (grp in unique(pw$tumor_supersite)) {
  v <- pw$decouple[pw$tumor_supersite == grp]
  cat(sprintf("%-16s n_pat=%3d med=%.3f\n", grp, length(v), median(v)))
}
w <- NULL
## Adnexa vs Omentum（患者级）
pw2 <- pw[pw$tumor_supersite %in% c("Adnexa", "Omentum"), ]
wo <- wilcox.test(decouple ~ tumor_supersite, data = pw2)
cat("\nPATIENT-LEVEL Adnexa vs Omentum p =", format(wo$p.value, digits = 4), "\n")

theme_set(theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35), panel.grid = element_blank()))
p1 <- ggplot(canc, aes(tumor_supersite, decouple, fill = tumor_supersite)) +
  geom_boxplot(outlier.size = .3, width = .5) +
  stat_compare_means(comparisons = list(c("Adnexa", "Omentum")), size = 2.2) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cancer-cell decoupling by site (pre-Rx)", x = NULL, y = "decoupling")
out <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/module_scores"
svglite::svglite(file.path(out, "GSE180661_A_site.svg"), width = 90/25.4, height = 70/25.4); print(p1); dev.off()
ragg::agg_tiff(file.path(out, "GSE180661_A_site.tiff"), width = 90/25.4, height = 70/25.4, units = "in", res = 600); print(p1); dev.off()
cat("DONE\n")
