
suppressMessages({library(data.table); library(ggplot2); library(ggpubr); library(ggsci); library(ragg)})
sdir <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/module_scores"
sc <- fread(file.path(sdir, "scores_GSE222556.csv"), data.table = TRUE)
cat("cells:", nrow(sc), "\n")
## 样本类型：sample 形如 P01_P1-1 / P02_M1-1 → 字母 P/M
sc$samp_letter <- sub("^[A-Z0-9]+_([A-Z])[0-9]+.*", "\\1", sc$sample)
sc$type <- ifelse(sc$samp_letter == "P", "Primary", ifelse(sc$samp_letter == "M", "Metastatic", NA))
sc$decouple <- as.numeric(scale(sc$M1_immune_activation)) - (as.numeric(scale(sc$M2_metabolic_GSH)) + as.numeric(scale(sc$M3_mito_cytoskeleton))) / 2
cat("type counts(cells):\n"); print(table(sc$type, useNA = "ifany"))
## 样本级中位
samp <- sc[, .(n_cells = .N, med = median(decouple)), by = .(sample, type)]
w1 <- wilcox.test(med ~ type, data = samp)
cat(sprintf("SAMPLE-level Primary vs Metastatic (n=%d vs %d): p=%.4f ; med %.3f vs %.3f\n",
            sum(samp$type == "Primary"), sum(samp$type == "Metastatic"), w1$p.value,
            median(samp$med[samp$type == "Primary"]), median(samp$med[samp$type == "Metastatic"])))
## 细胞级（全部 + 恶性/上皮亚群）
for (ct in c("Epithelial cell", "Myeloid cell", "T & NK cell")) {
  sub <- sc[sc$`Major cell type` == ct & !is.na(sc$type), ]
  if (nrow(sub) < 30) next
  w <- tryCatch(wilcox.test(decouple ~ type, data = sub), error = function(e) NULL)
  if (!is.null(w)) cat(sprintf("CELL-level %-16s: n=%d  p=%.4g  med P %.3f vs M %.3f\n", ct, nrow(sub),
      w$p.value, median(sub$decouple[sub$type == "Primary"]), median(sub$decouple[sub$type == "Metastatic"])))
}
fwrite(sc, file.path(sdir, "GSE222556_type_annotated.csv"))
fwrite(samp, file.path(sdir, "GSE222556_type_sample_summary.csv"))
cat("DONE 27a\n")
