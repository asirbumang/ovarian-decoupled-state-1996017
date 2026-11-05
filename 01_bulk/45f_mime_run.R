
suppressMessages({library(Mime1); library(data.table); library(dplyr)})
basedir <- "D:/子宫内膜异位症2/analysis/L0_bulk/results/track2/mime"
## ---- 构建两个队列 ----
mk <- function(exp_path, clin_path, name) {
  e <- as.data.frame(fread(exp_path, data.table = FALSE))
  rownames(e) <- e[[1]]; e <- e[, -1]
  c0 <- as.data.frame(fread(clin_path))
  rownames(c0) <- c0$sample
  samples <- intersect(rownames(c0), colnames(e))
  x <- t(e[, samples])
  d <- data.frame(ID = rownames(x), OS.time = c0[samples, "futime"], OS = c0[samples, "fustat"], x)
  rownames(d) <- NULL
  d
}
load_train <- tryCatch(mk(file.path(basedir, "tcga_exp.csv.gz"), file.path(basedir, "tcga_clin.csv"), "TCGA"),
                       error = function(e) NULL)
load_test <- tryCatch(mk(file.path(basedir, "GSE26193_exp.csv.gz"), file.path(basedir, "GSE26193_clin.csv"), "GSE26193"),
                      error = function(e) NULL)
if (is.null(load_train) || is.null(load_test)) { cat("INPUT_MISSING\n"); quit(status = 1) }
cat("train:", nrow(load_train), " test:", nrow(load_test), "\n")
genes_common <- intersect(colnames(load_train), colnames(load_test))
genes_common <- setdiff(genes_common, c("ID", "OS.time", "OS"))
cat("common genes:", length(genes_common), "\n")
## 候选基因 = 状态签名 ∩ common（400 签名在两个平台都有——若 <30 用 top-var common）
sig <- unique(c(readLines("D:/子宫内膜异位症2/analysis/L1_scRNA/results/cnv/state_sig_up.txt"),
                readLines("D:/子宫内膜异位症2/analysis/L1_scRNA/results/cnv/state_sig_dn.txt")))
cand <- intersect(sig, genes_common)
cat("candidate genes:", length(cand), "\n")
if (length(cand) < 30) stop("too few candidates")
sel <- c("ID", "OS.time", "OS", cand)
l <- list(Dataset1 = load_train[, sel], Dataset2 = load_test[, sel])
## OS.time 单位统一为月（253 天标准）→ Mime 内部用月/年换算；用天（文档默认天）
## rfsrc/兼容垫片：v3.7 已移除 var.select（原=按 VIMP 选基因）；按同语义实现
var.select <- function(object, conservative = "high") {
  imp <- object$importance
  if (is.null(imp)) return(rownames(object$variable.importance) %||% names(object$variable.importance))
  if (is.vector(imp) || is.null(dim(imp))) {
    nm <- names(imp); v <- as.numeric(imp)
  } else {
    nm <- rownames(imp); v <- as.numeric(imp[, 1])
  }
  if (is.null(nm)) return(character(0))
  sel <- nm[v > 0]
  sel <- sel[order(-v[match(sel, nm)])]
  if (length(sel) > 30) sel <- sel[1:30]
  if (length(sel) < 3) {
    ord <- nm[order(-v)]
    sel <- ord[1:min(5, length(ord))]
  }
  list(topvars = sel, vars = sel)
}
res <- ML.Dev.Prog.Sig(train_data = l[["Dataset1"]], list_train_vali_Data = l,
                       unicox.filter.for.candi = TRUE, unicox_p_cutoff = 0.05,
                       candidate_genes = cand, mode = "all", nodesize = 5, seed = 5201314)
res$Cindex.res %>% fwrite(file.path(basedir, "mime101_cindex.csv"))
saveRDS(res, file.path(basedir, "mime101_res.rds"))
print(head(res$Cindex.res, 10))
cat("DONE 45f\n")
