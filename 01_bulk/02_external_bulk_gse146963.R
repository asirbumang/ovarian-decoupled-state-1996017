# ============================================================
# 外部 bulk P0 · GSE146963（NAC 前后 28 对）/ GSE146964（初治多部位 38）
# 与 00_track0 同口径；GSE146963 主分析 = 配对 NACT 前后解耦变化
# ============================================================
suppressMessages({
  library(data.table); library(GSVA); library(ggplot2); library(ggpubr)
  library(patchwork); library(ggsci); library(ragg)
})
res_dir <- "D:/子宫内膜异位症2/analysis/L0_bulk/results/external"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

sets <- readRDS("D:/子宫内膜异位症2/pipeline/mitoxy_sets.RDS")
modules <- sets$modules; icc45 <- sets$icc45
mitoxy_genes <- unique(unlist(modules))
death_sets_clean <- list(
  pyroptosis  = setdiff(c("GSDMD","GSDMC","NLRP3","CASP1","CASP4","CASP5"), mitoxy_genes),
  necroptosis = setdiff(c("MLKL","RIPK3","RIPK1"), mitoxy_genes),
  apoptosis   = setdiff(c("CASP3","CASP7","CASP8","CASP9","BCL2","BAX"), mitoxy_genes),
  ferroptosis = setdiff(c("TFRC","NCOA4","FTH1","FTL","HMOX1","CHAC1","PTGS2","SLC3A2"), mitoxy_genes))
score_sets <- c(modules, list(ICC45 = icc45),
                list(pyroptosis = death_sets_clean$pyroptosis,
                     necroptosis = death_sets_clean$necroptosis,
                     apoptosis = death_sets_clean$apoptosis,
                     ferroptosis = death_sets_clean$ferroptosis))

score_cohort <- function(gse, title_tag) {
  expdf <- fread(paste0("D:/子宫内膜异位症2/data_omics/", gse, "_expr_symbol.csv.gz"), data.table = FALSE)
  rownames(expdf) <- expdf$gene; expdf$gene <- NULL
  expr <- as.matrix(expdf); rm(expdf); gc()
  storage.mode(expr) <- "double"
  meta <- fread(file.path("D:/子宫内膜异位症2/data_omics", paste0(gse, "_meta.csv")), data.table = FALSE)
  sc <- lapply(score_sets, function(g) intersect(g, rownames(expr)))
  cat(sprintf("\n[%s] gene coverage: %s\n", gse,
              paste(sprintf("%s=%d", names(sc), sapply(sc, length)), collapse = " | ")))
  g <- GSVA::gsva(GSVA::ssgseaParam(expr, sc))
  s <- t(as.matrix(g)); rm(g); gc()
  dec <- as.numeric(scale(s[, "M1_immune_activation"])) -
         (as.numeric(scale(s[, "M2_metabolic_GSH"])) + as.numeric(scale(s[, "M3_mito_cytoskeleton"]))) / 2
  out <- data.frame(sample = colnames(expr), decouple = dec, stringsAsFactors = FALSE)
  for (cn in colnames(s)) out[[cn]] <- s[, cn]
  out <- merge(out, meta, by.x = "sample", by.y = "gsm", all.x = TRUE)
  out
}

d1 <- score_cohort("GSE146963")
d2 <- score_cohort("GSE146964")

## ---------- GSE146963：NAC 前后配对 ----------
cat("\n== GSE146963 paired NACT ==\n")
d1$treat <- d1$`treatment status`
d1$pid <- sub("([A-Z]+[0-9]+)[A-Z]?$", "\\1", d1$title)  # NAC1A/NAC1B -> NAC1
d1$timepoint <- ifelse(grepl("A$", d1$title), "pre", ifelse(grepl("B$", d1$title), "post", NA))
cat("pairs:", table(d1$pid, useNA = "ifany"), "\n")
cat("pre n:", sum(d1$timepoint == "pre"), " post n:", sum(d1$timepoint == "post"), "\n")
d1w <- reshape(d1[, c("pid", "timepoint", "decouple")], idvar = "pid", timevar = "timepoint", direction = "wide")
d1w <- d1w[!is.na(d1w$decouple.pre) & !is.na(d1w$decouple.post), ]
cat("complete pairs:", nrow(d1w), "\n")
wp <- wilcox.test(d1w$decouple.post, d1w$decouple.pre, paired = TRUE)
cat("paired Wilcoxon decouple pre vs post: p =", format(wp$p.value, digits = 4),
    " med pre =", round(median(d1w$decouple.pre), 3), " med post =", round(median(d1w$decouple.post), 3), "\n")
for (m in names(modules)) {
  w1 <- d1[d1$timepoint == "pre", m]; w2 <- d1[d1$timepoint == "post", m]
  wm <- tryCatch(wilcox.test(w2, w1, paired = FALSE), error = function(e) NULL)
  if (!is.null(wm)) cat(sprintf("  %-24s p=%.4g  (pre %.3f vs post %.3f)\n", m, wm$p.value, median(w1), median(w2)))
}
fwrite(d1, file.path(res_dir, "ext_GSE146963_scores.csv"))
fwrite(d1w, file.path(res_dir, "ext_GSE146963_pairs.csv"))

## ---------- GSE146964：初治多部位 ----------
cat("\n== GSE146964 naive multi-site ==\n")
cat("tissues:\n"); print(table(d2$tissue))
kr <- kruskal.test(decouple ~ tissue, data = d2)
cat("Kruskal decouple x tissue: p =", format(kr$p.value, digits = 4), "\n")
for (tt in sort(unique(d2$tissue))) {
  s <- d2$decouple[d2$tissue == tt]
  cat(sprintf("  %-22s n=%2d median=%.3f\n", tt, length(s), median(s)))
}
fwrite(d2, file.path(res_dir, "ext_GSE146964_scores.csv"))

## ---------- 图（规范件）----------
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
npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488")
## 配对图：pre/post 连线 + 箱线
d1w$delta <- d1w$decouple.post - d1w$decouple.pre
p_pairs <- ggplot(d1w, aes(x = factor(1), y = delta)) + geom_boxplot(width = .3, outlier.size = 0.5, fill = npg[1], alpha = .7) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = .3) +
  annotate("text", x = 1, y = max(d1w$delta) * 1.02, label = paste0("paired Wilcoxon p = ", format(wp$p.value, digits = 3)), size = 2.4) +
  labs(title = "NACT-induced decoupling change (post-pre)", x = NULL, y = "Δ decoupling") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
npg_d <- c(pre = "#4DBBD5", post = "#E64B35")
p_pair2 <- ggplot(d1w, aes(x = decouple.pre, y = decouple.post)) + geom_point(size = .8, color = npg[4]) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = .3) +
  geom_smooth(method = "lm", se = TRUE, linewidth = .5, color = npg[1]) +
  labs(title = "Pre vs post decoupling (paired)", x = "pre-treatment decoupling", y = "post-treatment decoupling") +
  stat_cor(size = 2.4, method = "spearman")
p_box <- ggplot(d1, aes(x = timepoint, y = decouple, fill = timepoint)) + geom_boxplot(width = .5, outlier.size = .5) +
  scale_fill_manual(values = npg_d, guide = "none") +
  labs(title = "Decoupling by NACT timepoint", x = NULL, y = "decoupling")
cp <- (p_pairs + p_pair2 + p_box) + plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 10, face = "bold"))
save_pub_r(p_pairs, file.path(res_dir, "ext_GSE146963_A_delta"))
save_pub_r(p_pair2, file.path(res_dir, "ext_GSE146963_B_scatter"))
save_pub_r(p_box, file.path(res_dir, "ext_GSE146963_C_box"))
save_pub_r(cp, file.path(res_dir, "ext_GSE146963_composite"), width_mm = 180, height_mm = 60)
cat("\nDONE GSE146963/4\n")
