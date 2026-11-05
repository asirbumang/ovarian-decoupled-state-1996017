
suppressMessages({library(infercnv); library(data.table)})
rd <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/cnv/run"
obj <- readRDS(file.path(rd, "15_tumor_subclusters.leiden.infercnv_obj"))
expr <- as.matrix(obj@expr.data)   # 基因×细胞（ref 居中）
gord <- obj@gene_order
cat("genes:", nrow(expr), " cells:", ncol(expr), "\n")
gn <- rownames(expr)
chr_of <- gord$chr[match(gn, rownames(gord))]
chr_of <- as.character(chr_of)
## 每细胞每 chr 均值（均值对观测/参考敏感度 OK——用中位数更稳）
mat <- matrix(NA, 24, ncol(expr), dimnames = list(c(as.character(1:22), "X", "Y"), colnames(expr)))
for (ch in c(as.character(1:22), "X", "Y")) {
  sel <- which(chr_of == ch)
  if (length(sel) >= 3) mat[ch, ] <- colMeans(expr[sel, , drop = FALSE])
}
keep <- c(as.character(1:22))
m2 <- mat[keep, ]
## 保存
out <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/figs"
fwrite(data.table(chr = rownames(m2), m2), file.path(out, "F3a_chr_profile.csv.gz"))
cat("saved chr profile:", dim(m2), "\nDONE 44g\n")
