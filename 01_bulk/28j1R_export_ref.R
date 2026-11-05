
suppressMessages({library(Matrix); library(data.table)})
set.seed(2026)
LIN <- c("Epithelial cells" = "Epithelial", "Fibroblasts" = "Fibroblast", "Mesothelial cells" = "Mesothelial",
  "Endothelial cells" = "Endothelial", "Macrophages" = "MonoMac", "Classical monocytes" = "MonoMac",
  "Migratory DCs" = "DC", "DC1" = "DC", "DC2" = "DC", "pDC" = "DC", "Mast cells" = "Mast",
  "Tem/Trm cytotoxic T cells" = "T", "Tcm/Naive helper T cells" = "T", "Regulatory T cells" = "T",
  "Type 17 helper T cells" = "T", "NK cells" = "NK", "CD16- NK cells" = "NK", "CD16+ NK cells" = "NK",
  "ILC" = "NK", "Memory B cells" = "B", "Naive B cells" = "B", "Plasma cells" = "Plasma")
cat("reading counts...\n")
mtx <- readMM(gzfile("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_counts_raw.mtx.gz"))
mtx <- as(mtx, "dgCMatrix")
feat <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_seurat_features.txt.gz", header = FALSE, data.table = FALSE)
bar <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_barcodes.txt.gz", header = FALSE, data.table = FALSE)
rownames(mtx) <- feat$V1; colnames(mtx) <- bar$V1
rm(feat, bar); gc()
cat("counts:", dim(mtx), "\n")
sc <- fread("D:/子宫内膜异位症2/analysis/L1_scRNA/results/module_scores/scores_GSE266577.csv",
            select = c("cell", "cell_type"), data.table = TRUE)
sc$lineage <- unname(LIN[sc$cell_type])
set.seed(2026)
keep <- sc[, .(cell = sample(cell, min(300, .N))), by = lineage]$cell
cat("ref cells:", length(keep), "\n")
sub <- mtx[, keep]
cat("sub:", dim(sub), "\n")
s <- summary(sub)
dt <- data.table(g = rownames(sub)[s$i], c = colnames(sub)[s$j], v = s$x)
outdir <- "D:/子宫内膜异位症2/analysis/L0_bulk/results/spatial/deconv_input"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
fwrite(dt, file.path(outdir, "reference_lineage12.triplet.gz"))
fwrite(data.table(gene = rownames(sub)), file.path(outdir, "reference_lineage12.genes.tsv"))
fwrite(data.table(spot = colnames(sub), lineage = unname(LIN[sc$lineage[match(colnames(sub), sc$cell)]]),
                  cell_type = sc$cell_type[match(colnames(sub), sc$cell)]), file.path(outdir, "reference_lineage12.spots.tsv"))
cat("saved ref triplet nnz:", nrow(dt), "\nDONE 28j1R\n")
