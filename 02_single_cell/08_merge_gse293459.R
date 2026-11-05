
suppressMessages({library(Matrix); library(hdf5r); library(data.table)})
base <- "D:/子宫内膜异位症2/data_sc/GSE293459/rna"
files <- list.files(base, pattern = "filtered_feature_bc_matrix.h5", recursive = TRUE, full.names = TRUE)
cat("files:", length(files), "\n")
common <- NULL
mat_list <- list(); meta_list <- list()
for (fp in files) {
  gsm <- sub(".*/(GSM[0-9]+)/.*", "\\1", fp)
  base_name <- sub(".*_", "", sub("_filtered.*", "", fp))
  h5 <- hdf5r::H5File$new(fp, mode = "r")
  genes <- h5[["matrix/features/name"]][]
  barv <- h5[["matrix/barcodes"]][]
  dat <- h5[["matrix/data"]][]; idx <- h5[["matrix/indices"]][]
  ptr <- h5[["matrix/indptr"]][]; shp <- h5[["matrix/shape"]][]
  h5$close_all()
  gens <- sapply(genes, function(x) if (is.raw(x)) rawToChar(x) else as.character(x))
  bars <- sapply(barv, function(x) if (is.raw(x)) rawToChar(x) else as.character(x))
  m <- Matrix::sparseMatrix(i = idx + 1, j = rep(seq_len(shp[2]), diff(ptr)), x = dat,
                            dims = c(shp[1], shp[2]))
  rownames(m) <- makesymbols <- make.unique(gens)
  colnames(m) <- paste0(gsm, "_", bars)
  mat_list[[fp]] <- m
  meta_list[[fp]] <- data.frame(cell = colnames(m), sample = base_name, gsm = gsm, stringsAsFactors = FALSE)
  common <- if (is.null(common)) rownames(m) else intersect(common, rownames(m))
  cat("  ", gsm, base_name, "cells:", ncol(m), " genes:", nrow(m), " common:", length(common), "\n")
  rm(m); gc()
}
cat("common genes:", length(common), "\n")
## cbind 对齐
mats2 <- lapply(mat_list, function(m) m[common, ])
mm <- do.call(cbind, mats2)
cat("merged:", dim(mm), "\n")
meta_all <- rbindlist(meta_list)
meta_all$sample_type <- ifelse(grepl("_C$", meta_all$sample), "C",
                        ifelse(grepl("_T$", meta_all$sample), "T", "R"))
fwrite(meta_all, "D:/子宫内膜异位症2/data_sc/GSE293459/meta_cells.csv")
dir.create("D:/子宫内膜异位症2/data_sc/GSE293459/merged10x", showWarnings = FALSE)
mc <- as(mm, "CsparseMatrix")
writeMM(mc, file = "D:/子宫内膜异位症2/data_sc/GSE293459/merged10x/matrix.mtx")
writeLines(rownames(mc), "D:/子宫内膜异位症2/data_sc/GSE293459/merged10x/features.tsv")
writeLines(colnames(mc), "D:/子宫内膜异位症2/data_sc/GSE293459/merged10x/barcodes.tsv")
cat("DONE merge\n")
