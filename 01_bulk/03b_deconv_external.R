# ============================================================
# 03b: P1-1 外队 dec_res × OS 复算（26c 同口径 NNLS + 26d 残差 + Cox 连续/分层）
# 覆盖 GSE26193/14764/13876/63885 + TCGA-OV（用现有 deconv/dec_res 一致化）
# ============================================================
suppressMessages({library(data.table); library(Matrix); library(nnls); library(survival)
  library(AnnotationDbi); library(hgu133plus2.db); library(hgu133a2.db)})
set.seed(20260831)
tr <- "D:/子宫内膜异位症2/analysis/L0_bulk/results/track2"

parse_sm <- function(gz) {
  con <- gzfile(gz, "rt"); lines <- readLines(con, warn = FALSE); close(con)
  beg <- which(startsWith(lines, "!series_matrix_table_begin"))
  end <- which(startsWith(lines, "!series_matrix_table_end"))
  hdr <- strsplit(sub('"', "", lines[beg + 1]), "\t")[[1]]
  dat <- read.table(text = lines[(beg + 2):(end - 1)], sep = "\t", quote = "\"",
                    header = FALSE, row.names = 1, check.names = FALSE)
  colnames(dat) <- gsub('"', "", hdr[-1])
  char <- list()
  for (i in seq_len(beg)) {
    if (startsWith(lines[i], "!Sample_characteristics_ch1")) {
      parts <- gsub('"', "", strsplit(lines[i], "\t")[[1]][-1])
      if (length(parts) != ncol(dat)) next
      k <- sub(":.*$", "", parts[1])
      v <- sub("^[^:]+: ?", "", parts)
      char[[k]] <- v
    }
  }
  list(expr = as.matrix(dat), char = as.data.frame(char, stringsAsFactors = FALSE, check.names = FALSE))
}
p2g_db <- function(expr, db) {
  syms <- mapIds(db, keys = rownames(expr), keytype = "PROBEID", column = "SYMBOL", multiVals = "first")
  ok <- !is.na(syms)
  cat("  probe->symbol:", sum(ok), "/", nrow(expr), "\n")
  expr <- expr[ok, ]; rownames(expr) <- syms[ok]
  expr
}
p2g_map <- function(expr, mapfile) {
  mp <- fread(mapfile, data.table = FALSE)
  idx <- match(rownames(expr), as.character(mp$ID))
  syms <- mp$Gene_Symbol[idx]
  ok <- !is.na(syms) & syms != ""
  cat("  probe->symbol(custom):", sum(ok), "/", nrow(expr), "\n")
  expr <- expr[ok, ]; rownames(expr) <- syms[ok]
  expr
}

## 参考（26c 口径：6 谱系伪bulk marker z）
ref <- readRDS("D:/子宫内膜异位症2/analysis/L1_scRNA/results/track2/pb_ovarian_celltype_counts.RDS")
ref_cpm <- ref$cpm
rl <- log2(ref_cpm + 1)
spec <- rl - apply(rl, 1, function(x) sort(x, decreasing = TRUE)[2])
markers <- lapply(colnames(rl), function(ct) names(head(sort(spec[, ct], decreasing = TRUE), 90)))
mark_all <- sort(unique(unlist(markers)))
keep_g <- apply(rl[mark_all, ], 1, function(x) stats::var(x) > 1e-12)
mark_all <- mark_all[keep_g]
ref_z <- t(scale(t(rl[mark_all, ])))

deconv_fit <- function(bulk_log2, label) {
  mk <- intersect(mark_all, rownames(bulk_log2))
  b <- bulk_log2[mk, , drop = FALSE]
  keep_b <- apply(b, 1, function(x) stats::var(x) > 1e-12)
  mk <- mk[keep_b]; b <- b[mk, , drop = FALSE]
  bz <- t(scale(t(b)))
  A <- ref_z[mk, , drop = FALSE]
  est <- matrix(NA, ncol(bz), ncol(A), dimnames = list(colnames(bz), colnames(A)))
  for (j in seq_len(ncol(bz))) est[j, ] <- nnls(A, bz[, j])$x
  fc <- est / rowSums(est)
  out <- data.frame(sample = rownames(fc), fc, stringsAsFactors = FALSE)
  out$immune_frac <- (out$Myeloid + out$T_NK + out$B_Plasma) / rowSums(fc)
  cat("  ", label, "deconv n=", nrow(out), " median immune_frac:",
      round(median(out$immune_frac), 3), "\n", sep = "")
  out
}

cox_res <- function(df, lab) {
  df$z <- as.numeric(scale(df$dec_res))
  fx <- coxph(Surv(os_days, event) ~ z, data = df)
  co <- summary(fx)$coefficients
  df$tier <- cut(df$dec_res, quantile(df$dec_res, c(0, .25, .75, 1)), labels = c("Low", "Mid", "High"), include.lowest = TRUE)
  lr <- survdiff(Surv(os_days, event) ~ tier, data = df)
  out <- data.frame(cohort = lab, n = nrow(df), events = sum(df$event),
                    hr_cont = exp(co[1]), ci_lo = exp(co[1] - 1.96 * co[3]), ci_hi = exp(co[1] + 1.96 * co[3]),
                    p_cont = co[5],
                    logrank_p = pchisq(lr$chisq, 2, lower.tail = FALSE))
  cat(sprintf("  [%s] dec_res continuous HR=%.3f p=%.4f | logrank p=%.4f (n=%d ev=%d)\n",
              lab, out$hr_cont, out$p_cont, out$logrank_p, nrow(df), sum(df$event)))
  out
}

res_all <- list()
## ---------- GSE26193 ----------
d <- parse_sm("D:/子宫内膜异位症2/data_bulk/GSE26193_series_matrix.txt.gz")
e <- p2g_db(d$expr, hgu133plus2.db)
e <- e[!duplicated(rownames(e)), ]
if (max(e) > 50) e <- log2(e + 1)
dv <- deconv_fit(e, "26193")
sc <- fread("D:/子宫内膜异位症2/analysis/L0_bulk/results/survival_meta/scores_GSE26193.csv", data.table = FALSE)
df <- merge(sc, dv, by = "sample")
df$os_days <- df$os_days; df$event <- df$os_event
df <- df[!is.na(df$os_days) & df$os_days > 0 & !is.na(df$event), ]
df$dec_res <- resid(lm(decouple ~ immune_frac, data = df))
res_all[[1]] <- cox_res(df, "GSE26193")

## ---------- GSE14764 ----------
d <- parse_sm("D:/子宫内膜异位症2/data_bulk/GSE14764_series_matrix.txt.gz")
e <- p2g_db(d$expr, hgu133a2.db)
e <- e[!duplicated(rownames(e)), ]
if (max(e) > 50) e <- log2(e + 1)
dv <- deconv_fit(e, "14764")
sc <- fread("D:/子宫内膜异位症2/analysis/L0_bulk/results/survival_meta/scores_GSE14764.csv", data.table = FALSE)
df <- merge(sc, dv, by = "sample")
df$os_days <- df$os_days; df$event <- df$os_event
df <- df[!is.na(df$os_days) & df$os_days > 0 & !is.na(df$event), ]
df$dec_res <- resid(lm(decouple ~ immune_frac, data = df))
res_all[[2]] <- cox_res(df, "GSE14764")

## ---------- GSE13876（患者级折叠）----------
d <- parse_sm("D:/子宫内膜异位症2/data_bulk/GSE13876_series_matrix.txt.gz")
e <- p2g_map(d$expr, "D:/子宫内膜异位症2/data_omics/GSE13876_platform_table.csv.gz")
e <- e[!duplicated(rownames(e)), ]
if (max(e) > 50) e <- log2(e + 1)
dv <- deconv_fit(e, "13876")
sc <- fread("D:/子宫内膜异位症2/analysis/L0_bulk/results/survival_meta/scores_GSE13876_patient.csv", data.table = FALSE)
pat_id <- d$char[["assigned unique patient id"]]
dv$patient <- pat_id[match(dv$sample, colnames(d$expr))]
dvp <- as.data.table(dv)[, lapply(.SD, mean), by = patient, .SDcols = c("Epithelial", "Stromal", "Myeloid", "T_NK", "B_Plasma", "immune_frac")]
dvp <- as.data.frame(dvp)
df <- merge(sc, dvp, by = "patient")
df$os_days <- df$fumnd * 30.44
df$event <- df$status
df <- df[!is.na(df$os_days) & df$os_days > 0 & !is.na(df$event), ]
df$dec_res <- resid(lm(decouple ~ immune_frac, data = df))
res_all[[3]] <- cox_res(df, "GSE13876")

## ---------- GSE63885 ----------
d <- parse_sm("D:/子宫内膜异位症2/data_bulk/GSE63885_series_matrix.txt.gz")
e <- p2g_db(d$expr, hgu133plus2.db)
e <- e[!duplicated(rownames(e)), ]
if (max(e) > 50) e <- log2(e + 1)
dv <- deconv_fit(e, "63885")
sc <- fread("D:/子宫内膜异位症2/analysis/L0_bulk/results/survival_meta/scores_GSE63885.csv", data.table = FALSE)
cn <- function(pat) { hit <- grep(pat, names(d$char), value = TRUE)[1]; d$char[[hit]] }
df <- merge(sc, dv, by = "sample")
df$os_days <- as.numeric(cn("os - overall survival"))
df$status <- cn("clinical status at last follow-up")
df$event <- ifelse(df$status == "DOD", 1, 0)
df <- df[!is.na(df$os_days) & df$os_days > 0 & !is.na(df$event), ]
df$dec_res <- resid(lm(decouple ~ immune_frac, data = df))
res_all[[4]] <- cox_res(df, "GSE63885")

## ---------- TCGA-OV（复用既有 dec_res）----------
d <- fread(file.path(tr, "TCGA_OV_residual.csv"), data.table = FALSE)
df <- d[!is.na(d$event) & d$os_days > 0 & !is.na(d$dec_res), ]
df <- df[!is.na(df$dec_res), ]
res_all[[5]] <- cox_res(df, "TCGA-OV")

RR <- rbindlist(res_all, fill = TRUE)
fwrite(RR, file.path(tr, "decreg_cox_meta.csv"))
cat("\n==== dec_res×OS 汇总 ====\n"); print(as.data.frame(RR))
## 逆方差固定效应合并（连续 HR）
lg <- log(RR$hr_cont); se <- (log(RR$ci_hi) - log(RR$ci_lo)) / 3.92
pool_w <- sum(lg / se^2) / sum(1 / se^2); pool_se <- sqrt(1 / sum(1 / se^2))
cat(sprintf("\nPOOLED continuous HR=%.3f (%.3f–%.3f) p=%.4f  I2=%.1f%%\n",
            exp(pool_w), exp(pool_w - 1.96 * pool_se), exp(pool_w + 1.96 * pool_se),
            2 * pnorm(-abs(pool_w / pool_se)),
            100 * max(0, (sum(1 / se^2) - sum(1 / se^2 * (lg - pool_w)^2) / (length(lg) - 1)) / sum(1 / se^2))))
cat("DONE 03b\n")
