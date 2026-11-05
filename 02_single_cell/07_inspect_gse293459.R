
suppressMessages({library(Seurat)})
cat("loading rds...\n")
t0 <- Sys.time()
obj <- readRDS("D:/子宫内膜异位症2/data_sc/GSE293459/GSE293459_patient_seurat.rds")
cat("load:", round(as.numeric(difftime(Sys.time(), t0, units="mins")), 1), "min\n")
cat("class:", class(obj)[1], "\n")
if (inherits(obj, "Seurat")) {
  cat("assays:", paste(names(obj@assays), collapse=", "), "\n")
  cat("cells:", ncol(obj), "\n")
  for (a in names(obj@assays)) {
    cat("  assay", a, ":", class(obj@assays[[a]])[1], " genes:", nrow(obj@assays[[a]]), "\n")
  }
  cat("meta cols:", paste(head(colnames(obj@meta.data), 25), collapse=", "), "\n")
} else {
  cat("names:", names(obj), "\n")
  for (n in names(obj)) cat("  ", n, ":", class(obj[[n]])[1], "\n")
}
cat("DONE\n")
