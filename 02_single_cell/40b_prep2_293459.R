
suppressMessages({library(Matrix); library(data.table)})
cat("reading matrix_lf...\n")
mtx <- readMM("D:/子宫内膜异位症2/data_sc/GSE293459/merged10x/matrix_lf.mtx")
mtx <- as(mtx, "dgCMatrix")
fea <- fread("D:/子宫内膜异位症2/data_sc/GSE293459/merged10x/features.tsv", header = FALSE, data.table = FALSE)
bar <- fread("D:/子宫内膜异位症2/data_sc/GSE293459/merged10x/barcodes.tsv", header = FALSE, data.table = FALSE)
rownames(mtx) <- fea$V1; colnames(mtx) <- bar$V1
rm(fea, bar); gc()
s <- summary(mtx)
dt <- data.table(g = rownames(mtx)[s$i], c = colnames(mtx)[s$j], v = s$x)
outdir <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/cnv/b9_input"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
fwrite(dt, file.path(outdir, "GSE293459_all.triplet.gz"))
fwrite(data.table(gene = rownames(mtx)), file.path(outdir, "GSE293459_all.genes.tsv"))
cat("saved GSE293459_all nnz:", nrow(dt), "\nDONE prep2\n")
