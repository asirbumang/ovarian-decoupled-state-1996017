# ============================================================
# 26j: GSE201047（Cyprotex? —— "Distinct transcriptional programs of ascitic and solid
# cancer cells induce a different response to chemotherapy in HGSOC"）
# 22 个 10x h5：命名 N_/T_ × 部位(OT/A/PER/OM/BL) × PT#（N=化疗前 naive, T=treated）
# 逐样本 AddModuleScore → decouple → N vs T × 部位
# ============================================================
suppressMessages({
  library(Matrix); library(Seurat); library(data.table); library(ggplot2)
  library(ggpubr); library(patchwork); library(ggsci); library(ragg); library(hdf5r)
})
set.seed(42)
rawdir <- "D:/子宫内膜异位症2/data_sc/GSE201047/raw"
outdir <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/module_scores"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sets <- readRDS("D:/子宫内膜异位症2/pipeline/mitoxy_sets.RDS")
modules <- sets$modules; icc45 <- sets$icc45
mitoxy_genes <- unique(unlist(modules))
death_sets_clean <- list(
  pyroptosis  = setdiff(c("GSDMD","GSDMC","NLRP3","CASP1","CASP4","CASP5"), mitoxy_genes),
  necroptosis = setdiff(c("MLKL","RIPK3","RIPK1"), mitoxy_genes),
  apoptosis   = setdiff(c("CASP3","CASP7","CASP8","CASP9","BCL2","BAX"), mitoxy_genes),
  ferroptosis = setdiff(c("TFRC","NCOA4","FTH1","FTL","HMOX1","CHAC1","PTGS2","SLC3A2"), mitoxy_genes))
score_sets <- c(modules, list(ICC45 = icc45), list(pyroptosis = death_sets_clean$pyroptosis,
                 necroptosis = death_sets_clean$necroptosis,
                 apoptosis = death_sets_clean$apoptosis,
                 ferroptosis = death_sets_clean$ferroptosis))

files <- list.files(rawdir, pattern = "h5$")
cat("samples:", length(files), "\n")
res_list <- list()
for (fp in files) {
  sample_name <- sub("^GSM[0-9]+_[0-9]+_", "", sub("_filtered_gene_bc_matrices_h5\\.h5$", "", fp))  # N_OT_PT1
  cat("  ", sample_name, "\n", sep = "")
  ## 旧版 10x h5（hg19 组）手动读
  h5 <- hdf5r::H5File$new(file.path(rawdir, fp), mode = "r")
  g <- if ("hg19" %in% h5$ls()$name) "hg19" else h5$ls()$name[1]
  genes <- h5[[paste0(g, "/gene_names")]][]
  barv <- h5[[paste0(g, "/barcodes")]][]
  dat <- h5[[paste0(g, "/data")]][]
  idx <- h5[[paste0(g, "/indices")]][]
  ptr <- h5[[paste0(g, "/indptr")]][]
  shp <- h5[[paste0(g, "/shape")]][]
  h5$close_all()
  mtx <- Matrix::sparseMatrix(i = idx + 1,
                              j = rep(seq_len(shp[2]), diff(ptr)),
                              x = dat, dims = c(shp[1], shp[2]))
  rownames(mtx) <- genes; colnames(mtx) <- barv
  if (anyDuplicated(rownames(mtx))) mtx <- mtx[!duplicated(rownames(mtx)), ]
  nf <- Matrix::colSums(mtx > 0); nc <- Matrix::colSums(mtx)
  keep <- nf >= 100 & nc >= 200
  mtx <- mtx[, keep]
  seu <- CreateSeuratObject(counts = mtx, min.cells = 0, min.features = 0)
  seu <- NormalizeData(seu, verbose = FALSE)
  si <- lapply(score_sets, function(g) intersect(g, rownames(seu)))
  seu <- AddModuleScore(seu, features = si, name = "sc", nbin = 24, ctrl = 100)
  cols <- grep("^sc[0-9]+$", colnames(seu@meta.data), value = TRUE)
  stopifnot(length(cols) == length(score_sets))
  df <- data.frame(sample = sample_name, cell = colnames(seu), stringsAsFactors = FALSE)
  for (k in seq_along(score_sets)) df[[names(score_sets)[k]]] <- seu@meta.data[[cols[k]]]
  res_list[[fp]] <- df
  rm(seu, mtx); gc()
}
all <- rbindlist(res_list, fill = TRUE)
cat("total cells:", nrow(all), "\n")
mkz <- function(x) as.numeric(scale(x))
all[, M1z := mkz(M1_immune_activation)]
all[, M2z := mkz(M2_metabolic_GSH)]
all[, M3z := mkz(M3_mito_cytoskeleton)]
all[, decouple := M1z - (M2z + M3z) / 2]
all[, patient := sub("^.*_", "", sample)]
all[, treatment := sub("_.*", "", sample)]
all[, site := sub("^(N|T)_", "", sample)]
fwrite(all, file.path(outdir, "scores_GSE201047.csv"))

## 样本级
samp <- all[, .(n_cells = .N, med_dec = median(decouple)), by = .(sample, treatment, site, patient)]
fwrite(samp, file.path(outdir, "GSE201047_sample_summary.csv"))
cat("\n== samples ==\n"); print(as.data.frame(samp))
## 细胞级 N vs T
w <- wilcox.test(decouple ~ treatment, data = as.data.frame(all))
cat("\ncell-level N vs T (overall) p =", format(w$p.value, digits = 4), "\n")
## 样本级 N vs T（配对患者）
w1 <- tryCatch(wilcox.test(med_dec ~ treatment, data = as.data.frame(samp)), error = function(e) NULL)
cat("sample-level N vs T p =", if (!is.null(w1)) format(w1$p.value, digits = 4) else NA, "\n")
for (st in c("N", "T")) {
  med <- median(samp$med_dec[samp$treatment == st])
  cat(sprintf("  %s: n_samples=%d med=%.3f\n", st, sum(samp$treatment == st), med))
}
cat("DONE 26j\n")
