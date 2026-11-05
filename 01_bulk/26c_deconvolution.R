# ============================================================
# 26c: NNLS 解卷积（母版 28b 口径）—— TCGA-OV(427) + GSE156699(88)
# 参考：26a 卵巢 6 谱系伪bulk（GSE266577）；marker=每谱系 top90 特异性
# 输出: analysis/L0_bulk/results/track2/deconv_<cohort>.csv
# ============================================================
suppressMessages({library(data.table); library(Matrix); library(nnls)})
set.seed(20260831)
dir.create("D:/子宫内膜异位症2/analysis/L0_bulk/results/track2", recursive = TRUE, showWarnings = FALSE)

ref <- readRDS("D:/子宫内膜异位症2/analysis/L1_scRNA/results/track2/pb_ovarian_celltype_counts.RDS")
ref_cpm <- ref$cpm; ct_lev <- ref$lineage
cat("reference:", dim(ref_cpm), " lineages:", ct_lev, "\n")

## 特异性 marker（logCPM sdiff top90）
rl <- log2(ref_cpm + 1)
spec <- rl - apply(rl, 1, function(x) sort(x, decreasing = TRUE)[2])
markers <- lapply(colnames(rl), function(ct) names(head(sort(spec[, ct], decreasing = TRUE), 90)))
names(markers) <- colnames(rl)
cat("markers:", sapply(markers, length), "\n")
mark_all <- sort(unique(unlist(markers)))
rl_m <- rl[mark_all, , drop = FALSE]
## 剔除任何谱系零方差（z 标准化会产出 NA）
keep_g <- apply(rl_m, 1, function(x) stats::var(x) > 1e-12)
mark_all <- mark_all[keep_g]
ref_z <- t(scale(t(rl[mark_all, ])))
stopifnot(!any(is.na(ref_z)))
cat("ref_z:", dim(ref_z), "\n")

deconv <- function(bulk_log2, prefix) {
  ## bulk_log2: 基因 × 样本（log2 尺度）
  mk <- intersect(mark_all, rownames(bulk_log2))
  cat(prefix, "marker overlap:", length(mk), "/", length(mark_all), "\n")
  b <- bulk_log2[mk, , drop = FALSE]
  ## 剔除 bulk 中零方差基因（z 标准化 NA）
  keep_b <- apply(b, 1, function(x) stats::var(x) > 1e-12)
  mk <- mk[keep_b]
  b <- b[mk, , drop = FALSE]
  bz <- t(scale(t(b)))   # 每基因×样本 z
  stopifnot(!any(is.na(bz)))
  A <- ref_z[mk, , drop = FALSE]
  stopifnot(!any(is.na(A)), nrow(A) >= 50)
  est <- matrix(NA, ncol(bz), ncol(A), dimnames = list(colnames(bz), colnames(A)))
  fitv <- numeric(ncol(bz))
  for (j in seq_len(ncol(bz))) {
    sol <- nnls(A, bz[, j])
    est[j, ] <- sol$x
    fit <- 1 - sum((bz[, j] - A %*% sol$x)^2) / sum((bz[, j] - mean(bz[, j]))^2)
    fitv[j] <- fit
  }
  fc <- est / rowSums(est)
  out <- data.frame(sample = rownames(fc), fc, fit = fitv, stringsAsFactors = FALSE)
  out$immune_frac <- (out$Myeloid + out$T_NK + out$B_Plasma) / rowSums(fc)
  fwrite(out, paste0("D:/子宫内膜异位症2/analysis/L0_bulk/results/track2/deconv_", prefix, ".csv"))
  cat(prefix, "deconv saved, immune_frac median:", round(median(out$immune_frac), 3), "\n")
  invisible(out)
}

## ---- TCGA-OV ----
cat("loading TCGA counts...\n")
mx <- fread("D:/子宫内膜异位症2/data_omics/tcga_ov_counts_matrix.csv.gz", data.table = FALSE)
rownames(mx) <- mx$gene; mx$gene <- NULL
expr <- as.matrix(mx); storage.mode(expr) <- "double"; rm(mx); gc()
## 病例重复去重（与 track0 相同策略：保留文库最大）
lib <- colSums(expr)
dup <- unique(colnames(expr)[duplicated(colnames(expr))])
for (d in dup) { cols <- which(colnames(expr) == d); expr <- expr[, -setdiff(cols, cols[which.max(lib[cols])])] }
expr <- expr[!duplicated(rownames(expr)), ]
cpm <- sweep(expr, 2, colSums(expr) / 1e6, "/")
logt <- log2(cpm + 1); rm(expr, cpm); gc()
cat("TCGA log2:", dim(logt), "\n")
t1 <- deconv(logt, "TCGA_OV")

## ---- GSE156699 ----
cat("loading GSE156699...\n")
mx2 <- fread("D:/子宫内膜异位症2/data_bulk/GSE156699_UI_mRNA.txt.gz", data.table = FALSE)
colnames(mx2)[1] <- "gene"; mx2 <- mx2[!duplicated(mx2$gene), ]
rownames(mx2) <- mx2$gene; mx2$gene <- NULL
e2 <- as.matrix(mx2); storage.mode(e2) <- "double"; rm(mx2); gc()
t2 <- deconv(e2, "GSE156699")
cat("DONE 26c\n")
