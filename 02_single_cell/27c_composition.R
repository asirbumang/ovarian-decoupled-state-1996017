
suppressMessages({library(data.table); library(ggplot2); library(ggsci); library(ragg)})
sdir <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/module_scores"
sc <- fread(file.path(sdir, "scores_GSE266577.csv"), data.table = TRUE)
## 12 谱系映射（03 决议；ILC 并入 NK；Th17/DC2/pDC 并入相应大群）
LIN <- c(
  "Epithelial cells" = "Epithelial", "Fibroblasts" = "Fibroblast", "Mesothelial cells" = "Mesothelial",
  "Endothelial cells" = "Endothelial", "Macrophages" = "MonoMac", "Classical monocytes" = "MonoMac",
  "Migratory DCs" = "DC", "DC1" = "DC", "DC2" = "DC", "pDC" = "DC", "Mast cells" = "Mast",
  "Tem/Trm cytotoxic T cells" = "T", "Tcm/Naive helper T cells" = "T", "Regulatory T cells" = "T",
  "Type 17 helper T cells" = "T", "NK cells" = "NK", "CD16- NK cells" = "NK", "CD16+ NK cells" = "NK",
  "ILC" = "NK", "Memory B cells" = "B", "Naive B cells" = "B", "Plasma cells" = "Plasma")
sc$lineage <- LIN[sc$cell_type]
cat("lineage table:\n"); print(table(sc$lineage))
## 样本级组成 × 阶段
comp <- sc[, .N, by = .(sample = publication_sample_code_final, stage = treatment_stage, lineage)]
comp <- comp[, frac := N / sum(N), by = sample]
cat("\nsamples:", comp[, uniqueN(sample)], "\n")
## 组间组成差异（各谱系：chemo-naive vs IDS 样本级 Wilcoxon）
setorder(comp, stage, sample)
res <- list()
for (ln in unique(comp$lineage)) {
  sub <- comp[lineage == ln]
  w <- wilcox.test(frac ~ stage, data = sub)
  res[[ln]] <- data.frame(lineage = ln, med_naive = median(sub$frac[sub$stage == "chemo-naive"]),
                          med_ids = median(sub$frac[sub$stage == "IDS"]), p = w$p.value)
}
res <- rbindlist(res)
print(res[order(p)])
fwrite(res, file.path(sdir, "GSE266577_composition_by_stage.csv"))
cat("DONE 27c\n")
