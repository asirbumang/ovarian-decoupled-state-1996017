# A4d: UCell on a stratified subsample (rank-based scoring is the bottleneck)
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

sc <- fread(file.path(root, "analysis", "L1_scRNA", "results", "module_scores",
                      "scores_GSE266577.csv"), select = c("cell", "cell_type", "state_high"))
set.seed(20261105)
idx <- unlist(lapply(split(seq_len(nrow(sc)), sc$cell_type), function(ix)
  sample(ix, max(20, round(length(ix) * 30000 / nrow(sc))))))
idx <- sort(idx)
cat("subsample cells:", length(idx), "state:", sum(sc$state_high[idx]), "\n")
fwrite(data.frame(cell = sc$cell[idx]), file.path(outdir, "A4_ucell_subsample_cells.csv"))

cat("reading matrix ...\n"); t0 <- Sys.time()
counts <- readMM(gzfile(file.path(root, "data_sc", "GSE266577", "GSE266577_counts_raw.mtx.gz")))
counts <- as(counts, "dgCMatrix")
feat <- fread(file.path(root, "data_sc", "GSE266577", "GSE266577_seurat_features.txt.gz"), header = FALSE)[[1]]
bar <- fread(file.path(root, "data_sc", "GSE266577", "GSE266577_barcodes.txt.gz"), header = FALSE)[[1]]
rownames(counts) <- feat; colnames(counts) <- bar
counts <- counts[, sc$cell[idx]]
cat(sprintf("matrix %d x %d (%.1f min)\n", nrow(counts), ncol(counts),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
features <- lapply(features, function(g) intersect(g, rownames(counts)))
t1 <- Sys.time()
scores <- ScoreSignatures_UCell(counts, features = features, maxRank = 1500,
                                ncores = 1, chunk.size = 500)
cat(sprintf("UCell done in %.1f min\n", as.numeric(difftime(Sys.time(), t1, units = "mins"))))
scores <- data.frame(cell = rownames(scores), scores, check.names = FALSE)
fwrite(scores, file.path(outdir, "A4_ucell_subsample.csv"))
cat("saved\n")
