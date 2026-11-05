# A4: UCell scores (rank-based, depth-robust) for the module and ICC45 gene sets
# Usage: Rscript A4_ucell.R [n_cells_for_test]
suppressMessages({
  library(Matrix); library(data.table); library(UCell)
})
args <- commandArgs(trailingOnly = TRUE)
ntest <- if (length(args) >= 1) as.integer(args[1]) else 0L

root <- "D:/子宫内膜异位症2"
outdir <- file.path(root, "R2_revision", "04_analysis", "results")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## gene sets
gs <- fread(file.path(root, "_gene_sets.csv"), data.table = FALSE)
sets <- split(gs$gene, gs$set)
mods <- sets[c("M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton",
               "M4_lipid_perox_rupture")]
icc45 <- sets[["ICC45"]]
shared <- intersect(icc45, unique(unlist(mods)))
icc_unique <- setdiff(icc45, shared)
cat("ICC45", length(icc45), "shared", length(shared), "unique", length(icc_unique), "\n")

set.seed(20261105)
rand17 <- replicate(20, sample(unique(gs$gene), 17), simplify = FALSE)
rand45 <- replicate(20, sample(unique(gs$gene), 45), simplify = FALSE)
names(rand17) <- paste0("RAND17_", 1:20)
names(rand45) <- paste0("RAND45_", 1:20)

features <- c(mods, list(ICC45_full = icc45, ICC45_shared28 = shared,
                         ICC45_unique17 = icc_unique), rand17, rand45)

## matrix
cat("reading matrix ...\n"); t0 <- Sys.time()
counts <- readMM(gzfile(file.path(root, "data_sc", "GSE266577",
                                 "GSE266577_counts_raw.mtx.gz")))
counts <- as(counts, "dgCMatrix")
feat <- fread(file.path(root, "data_sc", "GSE266577", "GSE266577_seurat_features.txt.gz"),
              header = FALSE)[[1]]
bar <- fread(file.path(root, "data_sc", "GSE266577", "GSE266577_barcodes.txt.gz"),
             header = FALSE)[[1]]
rownames(counts) <- feat; colnames(counts) <- bar
cat(sprintf("matrix %d x %d (%.1f min)\n", nrow(counts), ncol(counts),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

## keep QC-passed cells from the scoring table
sc <- fread(file.path(root, "analysis", "L1_scRNA", "results", "module_scores",
                      "scores_GSE266577.csv"), select = c("cell", "cell_type"))
counts <- counts[, sc$cell]
if (ntest > 0) {
  idx <- sort(sample(ncol(counts), min(ntest, ncol(counts))))
  counts <- counts[, idx]
}
cat("cells:", ncol(counts), "\n")

features <- lapply(features, function(g) intersect(g, rownames(counts)))
cat("set sizes after intersect:\n"); print(sapply(features, length))

t1 <- Sys.time()
scores <- ScoreSignatures_UCell(counts, features = features, maxRank = 1500,
                                ncores = 8, chunk.size = 200)
cat(sprintf("UCell done in %.1f min\n", as.numeric(difftime(Sys.time(), t1, units = "mins"))))
scores <- data.frame(cell = rownames(scores), scores, check.names = FALSE)
fwrite(scores, file.path(outdir, if (ntest > 0) "A4_ucell_TEST.csv" else "A4_ucell.csv"))
cat("saved\n")
