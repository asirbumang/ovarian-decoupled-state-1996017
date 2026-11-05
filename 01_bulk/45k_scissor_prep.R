
suppressMessages({library(data.table)})
d <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/cnv/umap_input"
genes <- readLines(file.path(d, "GSE266577_full.genes.tsv"))
cells <- readLines(file.path(d, "GSE266577_full.cells.tsv"))
dt <- fread(file.path(d, "GSE266577_full.triplet.gz"))
gidx <- setNames(seq_along(genes), genes)
sidx <- setNames(seq_along(cells), cells)
dt[, gi := gidx[g]]; dt[, ci := sidx[c]]
cat("triplet:", nrow(dt), "\n")
sig <- unique(c(readLines("D:/子宫内膜异位症2/analysis/L1_scRNA/results/cnv/state_sig_up.txt"),
                readLines("D:/子宫内膜异位症2/analysis/L1_scRNA/results/cnv/state_sig_dn.txt")))
sig <- intersect(sig, genes)
cat("sig genes:", length(sig), "\n")
mat <- Matrix::sparseMatrix(i = gidx[dt$g], j = dt$ci, x = dt$v,
                            dims = c(length(genes), length(cells)))
rm(dt); gc()
xs <- mat; xs@x <- xs@x^2
varr <- Matrix::rowMeans(xs) - Matrix::rowMeans(mat)^2
rm(xs); gc()
ordv <- order(-varr)
topv <- setdiff(genes[ordv], sig)[1:2500]
sel <- intersect(unique(c(sig, topv)), genes)
cat("feature genes (with topvar):", length(sel), "\n")
m2 <- Matrix::Matrix(mat[gidx[sel], , drop = FALSE], sparse = TRUE)
rm(mat); gc()
m2 <- m2 * (1e4 / Matrix::colSums(m2)[col(m2)])
m2 <- log1p(m2)
out <- "D:/子宫内膜异位症2/analysis/L0_bulk/results/track2/mime"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(expr = m2, genes = sel, cells = cells), file.path(out, "scissor_scRNA.rds"))
cat("saved scRNA:", dim(m2), "\nDONE 45k\n")
