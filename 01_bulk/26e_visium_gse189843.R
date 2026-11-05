# ============================================================
# 26e: GSE189843 Visium spot 级解耦 × NACT 反应（ER 6 vs PR 6）
# 每 spot: AddModuleScore（母版 27a 参数）；样本级=spot 中位 + 高解耦 spot 比例
# ============================================================
suppressMessages({
  library(Matrix); library(Seurat); library(data.table); library(ggplot2)
  library(ggpubr); library(patchwork); library(ggsci); library(ragg)
})
set.seed(42)
sd <- "D:/子宫内膜异位症2/data_sc/GSE189843/raw"
outdir <- "D:/子宫内膜异位症2/analysis/L0_bulk/results/spatial"
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

## 样本分组表（从 GEO 元数据：ER 1-6 / PR 1-6）
grp_map <- data.frame(gsm = sprintf("GSM57084%02d", 85:96),
                      sample = c(paste0("ER", 1:6), paste0("PR", 1:6)),
                      stringsAsFactors = FALSE)
stopifnot(nrow(grp_map) == 12)

files <- list.files(sd, pattern = "matrix_.*\\.mtx\\.gz$")
res_list <- list()
for (fi in files) {
  gsm <- sub("_matrix_.*", "", fi)
  ii <- sub(".*_matrix_(II\\d+)\\.mtx\\.gz", "\\1", fi)
  gsm_num <- as.integer(sub("GSM", "", gsm))
  grp <- grp_map$sample[match(gsm, grp_map$gsm)]
  cat("==", gsm, grp, "==\n")
  mtx <- readMM(gzfile(file.path(sd, fi))); mtx <- as(mtx, "dgCMatrix")
  fea <- fread(file.path(sd, sub("\\.mtx\\.gz", "", fi) %>% sub("matrix", "features", .) %>% paste0(".tsv.gz")), header = FALSE, data.table = FALSE)
  bar <- fread(file.path(sd, sub("\\.mtx\\.gz", "", fi) %>% sub("matrix", "barcodes", .) %>% paste0(".tsv.gz")), header = FALSE, data.table = FALSE)
  rownames(mtx) <- fea$V2   # 10x v3 features: gene symbol 第2列
  colnames(mtx) <- bar$V1
  if (anyDuplicated(rownames(mtx))) {
    mtx <- mtx[!duplicated(rownames(mtx)), ]
    cat("   dedup genes:", nrow(mtx), "\n")
  }
  ## QC：>200 基因 & >500 counts 的 spot
  nf <- Matrix::colSums(mtx > 0); nc <- Matrix::colSums(mtx)
  keep <- nf >= 100 & nc >= 200
  mtx <- mtx[, keep]
  cat("   spots:", ncol(mtx), " genes:", nrow(mtx), "\n")
  seu <- CreateSeuratObject(counts = mtx, min.cells = 0, min.features = 0)
  seu <- NormalizeData(seu, verbose = FALSE)
  si <- lapply(score_sets, function(g) intersect(g, rownames(seu)))
  cov <- sapply(si, length)
  seu <- AddModuleScore(seu, features = si, name = "sc", nbin = 24, ctrl = 100)
  cols <- grep("^sc[0-9]+$", colnames(seu@meta.data), value = TRUE)
  stopifnot(length(cols) == length(score_sets))
  df <- data.frame(gsm = gsm, sample = grp, spot = colnames(seu), stringsAsFactors = FALSE)
  for (k in seq_along(score_sets)) df[[names(score_sets)[k]]] <- seu@meta.data[[cols[k]]]
  res_list[[fi]] <- df
  rm(seu, mtx); gc()
}
all <- rbindlist(res_list, fill = TRUE)
cat("total spots:", nrow(all), "\n")

## 解耦 + 严格态（全部 spot 合并标准化 —— 与 bulk 定义一致）
mkz <- function(x) as.numeric(scale(x))
all[, M1z := mkz(M1_immune_activation)]
all[, M2z := mkz(M2_metabolic_GSH)]
all[, M3z := mkz(M3_mito_cytoskeleton)]
all[, decouple := M1z - (M2z + M3z) / 2]
all[, state_high := M1_immune_activation >= quantile(M1_immune_activation, .75) &
     M2_metabolic_GSH <= quantile(M2_metabolic_GSH, .25) &
     M3_mito_cytoskeleton <= quantile(M3_mito_cytoskeleton, .25)]
fwrite(all, file.path(outdir, "GSE189843_spot_scores.csv"))

## 样本级汇总 + 组间检验
samp <- all[, .(n_spots = .N, med_dec = median(decouple), pct_state = 100 * mean(state_high)),
            by = .(sample)]
samp[, grp := ifelse(grepl("^ER", sample), "ER", "PR")]
cat("\n== 样本级 ==\n"); print(samp)
w1 <- tryCatch(wilcox.test(med_dec ~ grp, data = samp), error = function(e) NULL)
w2 <- tryCatch(wilcox.test(pct_state ~ grp, data = samp), error = function(e) NULL)
cat("样本级 Wilcoxon med_dec: ER vs PR p =", if (!is.null(w1)) format(w1$p.value, digits = 4) else NA, "\n")
cat("样本级 Wilcoxon pct_state: ER vs PR p =", if (!is.null(w2)) format(w2$p.value, digits = 4) else NA, "\n")
cat("ER med:", round(median(samp$med_dec[samp$grp == "ER"]), 3),
    " PR med:", round(median(samp$med_dec[samp$grp == "PR"]), 3), "\n")

## 图
theme_set(theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        legend.title = element_text(size = 6.2), legend.text = element_text(size = 5.8),
        plot.title = element_text(size = 7, face = "bold"), panel.grid = element_blank()))
save_pub_r <- function(p, fn, width_mm = 90, height_mm = 70, dpi = 600) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(fn, ".svg"), width = w, height = h); print(p); dev.off()
  grDevices::cairo_pdf(paste0(fn, ".pdf"), width = w, height = h, family = "Arial"); print(p); dev.off()
  ragg::agg_tiff(paste0(fn, ".tiff"), width = w, height = h, units = "in", res = dpi); print(p); dev.off()
}
p1 <- ggplot(samp, aes(grp, med_dec, fill = grp)) + geom_boxplot(width = .5, outlier.size = .5) +
  geom_jitter(width = .08, size = .6) + stat_compare_means(size = 2.2) +
  scale_fill_manual(values = c(ER = "#00A087", PR = "#DC0000"), guide = "none") +
  labs(title = "Spot-level decoupling: ER vs PR", x = "NACT response", y = "median spot decoupling")
p2 <- ggplot(all, aes(x = decouple, fill = sample)) + geom_density(alpha = .3, linewidth = .2) +
  scale_fill_manual(values = rep(NA, 12)) +
  labs(title = "Spot decoupling density", x = "decoupling", y = "density") +
  theme(legend.position = "none")
cp <- (p1 + p2) + plot_layout(ncol = 2) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))
save_pub_r(p1, file.path(outdir, "GSE189843_A_box"))
save_pub_r(p2, file.path(outdir, "GSE189843_B_density"))
save_pub_r(cp, file.path(outdir, "GSE189843_composite"), 180, 60)
fwrite(samp, file.path(outdir, "GSE189843_sample_summary.csv"))
cat("DONE 26e\n")
