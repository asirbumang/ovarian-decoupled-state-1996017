
suppressMessages({library(data.table)})
sdir <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/module_scores"
df <- fread(file.path(sdir, "scores_GSE293459.csv"), data.table = TRUE)
## 清掉旧错误列（08_merge 解析 bug）；gsm 从细胞名前缀提取
keep <- setdiff(colnames(df), c("sample", "sample_type", "gsm"))
df <- df[, ..keep]
df$gsm <- sub("^(GSM[0-9]+)_.*", "\\1", df$cell)
## 用修好的 meta（gsm → sample/sample_type）
sampmap <- fread("D:/子宫内膜异位症2/data_sc/GSE293459/meta_samples.csv", data.table = TRUE)
df <- merge(df, sampmap, by = "gsm", all.x = TRUE)
fwrite(df, file.path(sdir, "scores_GSE293459.csv"))
cat("cells:", nrow(df), " C-cell:", sum(df$sample_type == "C", na.rm = TRUE),
    " T-cell:", sum(df$sample_type == "T", na.rm = TRUE), "\n")
## 样本级中位
samp <- df[, .(n_cells = .N, med_dec = median(decouple)), by = .(sample, sample_type)]
print(as.data.frame(samp), nrow = 20)
## 配对：3 对患者（PH235/PH27/PH626）
pair <- samp[sample_type %in% c("C", "T")]
pair$patient <- sub("_[CT]$", "", pair$sample)
w <- wilcox.test(med_dec ~ sample_type, data = as.data.frame(pair))
cat("\nC vs T (6 samples, 3 pairs): p=", format(w$p.value, digits = 4), "\n", sep = "")
print(as.data.frame(pair), nrow = 10)
## 细胞级（全）
wc <- wilcox.test(decouple ~ sample_type, data = as.data.frame(df[sample_type %in% c("C", "T")]))
cat("CELL-level C vs T: p=", format(wc$p.value, digits = 4), "\n", sep = "")
cat("med C:", round(median(df$decouple[df$sample_type == "C"], na.rm = TRUE), 3),
    " vs T:", round(median(df$decouple[df$sample_type == "T"], na.rm = TRUE), 3), "\n")
cat("DONE 27b-supp\n")
