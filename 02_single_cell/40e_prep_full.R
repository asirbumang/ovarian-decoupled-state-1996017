
suppressMessages({library(Matrix); library(data.table)})
cat("reading...\n")
mtx <- readMM(gzfile("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_counts_raw.mtx.gz"))
mtx <- as(mtx, "dgCMatrix")
feat <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_seurat_features.txt.gz", header = FALSE, data.table = FALSE)
bar <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_barcodes.txt.gz", header = FALSE, data.table = FALSE)
rownames(mtx) <- feat$V1; colnames(mtx) <- bar$V1
rm(feat, bar); gc()
cat("full:", dim(mtx), "\n")
s <- summary(mtx)
dt <- data.table(g = rownames(mtx)[s$i], c = colnames(mtx)[s$j], v = s$x)
outdir <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/cnv/umap_input"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
fwrite(dt, file.path(outdir, "GSE266577_full.triplet.gz"))
fwrite(data.table(gene = rownames(mtx)), file.path(outdir, "GSE266577_full.genes.tsv"))
fwrite(data.table(cell = colnames(mtx)), file.path(outdir, "GSE266577_full.cells.tsv"))
cat("saved nnz:", nrow(dt), "\nDONE 40e-prep\n")
