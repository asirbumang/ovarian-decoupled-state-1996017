
suppressMessages({library(data.table); library(Matrix)})
basedir <- "D:/子宫内膜异位症2/analysis/L0_bulk/results/track2/mime"
sc <- readRDS(file.path(basedir, "scissor_scRNA.rds"))
set.seed(42)
sel <- sample(sc$cells, 60000)
pos <- match(sel, sc$cells)
m <- as.matrix(sc$expr[, pos])
rownames(m) <- sc$genes; colnames(m) <- sel
out <- data.table(gene = rownames(m), m)
fwrite(out, file.path(basedir, "scissor_expr60k.csv.gz"))
fwrite(data.table(cell = sel), file.path(basedir, "scissor_cells60k.csv"))
cat("saved", dim(m), "\nDONE 45n v2\n")
