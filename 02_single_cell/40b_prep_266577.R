
suppressMessages({library(Matrix); library(data.table)})
cat("reading...\n")
mtx <- readMM(gzfile("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_counts_raw.mtx.gz"))
mtx <- as(mtx, "dgCMatrix")
feat <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_seurat_features.txt.gz", header = FALSE, data.table = FALSE)
bar <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_barcodes.txt.gz", header = FALSE, data.table = FALSE)
rownames(mtx) <- feat$V1; colnames(mtx) <- bar$V1
rm(feat, bar); gc()
sc <- fread("D:/子宫内膜异位症2/analysis/L1_scRNA/results/module_scores/scores_GSE266577.csv",
            select = c("cell", "cell_type"), data.table = TRUE)
cnv <- fread("D:/子宫内膜异位症2/analysis/L1_scRNA/results/cnv/percell_cnv_score_v2.csv", data.table = TRUE)
mal <- cnv[malignant_like == TRUE, cell]
epi <- sc[cell_type == "Epithelial cells", cell]
mal <- intersect(mal, colnames(mtx)); epi <- intersect(epi, colnames(mtx))
cat("epi:", length(epi), " mal:", length(mal), "\n")
outdir <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/cnv/b9_input"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
for (nm in c("epi", "mal")) {
  cells <- if (nm == "epi") epi else mal
  sub <- mtx[, cells]
  s <- summary(sub)
  dt <- data.table(g = rownames(sub)[s$i], c = colnames(sub)[s$j], v = s$x)
  fwrite(dt, file.path(outdir, paste0("GSE266577_", nm, ".triplet.gz")))
  fwrite(data.table(gene = rownames(sub)), file.path(outdir, paste0("GSE266577_", nm, ".genes.tsv")))
  cat("saved", nm, "nnz:", nrow(dt), "\n")
}
cat("DONE prep\n")
