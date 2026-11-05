# A4b: UCell scores for the essential signatures only (fast mode)
suppressMessages({library(Matrix); library(data.table); library(UCell)})
root <- "D:/子宫内膜异位症2"
outdir <- file.path(root, "R2_revision", "04_analysis", "results")
gs <- fread(file.path(root, "_gene_sets.csv"), data.table = FALSE)
sets <- split(gs$gene, gs$set)
mods <- sets[c("M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton")]
icc45 <- sets[["ICC45"]]
shared <- intersect(icc45, unique(unlist(mods)))
icc_unique <- setdiff(icc45, shared)
set.seed(20261105)
rand17 <- replicate(10, sample(unique(gs$gene), 17), simplify = FALSE)
names(rand17) <- paste0("RAND17_", 1:10)
features <- c(mods, list(ICC45_full = icc45, ICC45_shared28 = shared,
                         ICC45_unique17 = icc_unique), rand17)

cat("reading matrix ...\n"); t0 <- Sys.time()
counts <- readMM(gzfile(file.path(root, "data_sc", "GSE266577",
                                 "GSE266577_counts_raw.mtx.gz")))
counts <- as(counts, "dgCMatrix")
feat <- fread(file.path(root, "data_sc", "GSE266577", "GSE266577_seurat_features.txt.gz"), header = FALSE)[[1]]
bar <- fread(file.path(root, "data_sc", "GSE266577", "GSE266577_barcodes.txt.gz"), header = FALSE)[[1]]
rownames(counts) <- feat; colnames(counts) <- bar
cat(sprintf("matrix %d x %d (%.1f min)\n", nrow(counts), ncol(counts),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
sc <- fread(file.path(root, "analysis", "L1_scRNA", "results", "module_scores",
                      "scores_GSE266577.csv"), select = c("cell"))
counts <- counts[, sc$cell]
cat("cells:", ncol(counts), "\n")
features <- lapply(features, function(g) intersect(g, rownames(counts)))
print(sapply(features, length))
t1 <- Sys.time()
scores <- ScoreSignatures_UCell(counts, features = features, maxRank = 1500,
                                ncores = 1, chunk.size = 500)
cat(sprintf("UCell done in %.1f min\n", as.numeric(difftime(Sys.time(), t1, units = "mins"))))
scores <- data.frame(cell = rownames(scores), scores, check.names = FALSE)
fwrite(scores, file.path(outdir, "A4_ucell.csv"))
cat("saved\n")
