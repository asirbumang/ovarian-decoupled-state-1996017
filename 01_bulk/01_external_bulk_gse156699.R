# ============================================================
# 外部 bulk P0 · GSE156699（88 例 HGSC，反应者=PFS≥6mo）
# 与 00_track0_tcga.R 同口径：同一基因集/同一 ssGSEA 参数/同一状态与解耦分数定义
# 输入：data_bulk/GSE156699_UI_mRNA.txt.gz（log2 表达，Symbol × 111 患者）
#       data_omics/gse156699_meta.csv（88 例：Stage + response）
# 输出：analysis/L0_bulk/results/ext_GSE156699_*
# ============================================================
suppressMessages({
  library(data.table); library(GSVA); library(ggplot2); library(ggpubr)
  library(patchwork); library(survival); library(survminer); library(ggsci); library(ragg)
})
set.seed(20260831)
res_dir <- "D:/子宫内膜异位症2/analysis/L0_bulk/results"
out_pref <- "ext_GSE156699"
direct <- file.path(res_dir, "external")
dir.create(direct, showWarnings = FALSE, recursive = TRUE)

sets <- readRDS("D:/子宫内膜异位症2/pipeline/mitoxy_sets.RDS")
modules <- sets$modules; icc45 <- sets$icc45
mitoxy_genes <- unique(unlist(modules))
death_sets_clean <- list(
  pyroptosis  = setdiff(c("GSDMD","GSDMC","NLRP3","CASP1","CASP4","CASP5"), mitoxy_genes),
  necroptosis = setdiff(c("MLKL","RIPK3","RIPK1"), mitoxy_genes),
  apoptosis   = setdiff(c("CASP3","CASP7","CASP8","CASP9","BCL2","BAX"), mitoxy_genes),
  ferroptosis = setdiff(c("TFRC","NCOA4","FTH1","FTL","HMOX1","CHAC1","PTGS2","SLC3A2"), mitoxy_genes))

## 1. 读表达（Symbol × 患者；已是 log2 尺度）
cat("loading GSE156699 UI_mRNA...\n")
mx <- fread("D:/子宫内膜异位症2/data_bulk/GSE156699_UI_mRNA.txt.gz", data.table = FALSE)
colnames(mx)[1] <- "gene"
mx <- mx[!duplicated(mx$gene), ]
genes <- mx$gene
expr <- as.matrix(mx[, -1]); rownames(expr) <- mx$gene; rm(mx); gc()
cat("dim:", dim(expr), "\n")
## 质检：全部有限值、无全 0
cat("NA values:", sum(is.na(expr)), " zero rows:", sum(rowSums(expr, na.rm = TRUE) == 0), "\n")

## 2. 元数据（88 例）对齐
meta <- fread("D:/子宫内膜异位症2/data_omics/gse156699_meta.csv", data.table = FALSE)
meta$patient <- meta$`subect id`
meta$resp <- meta$`response to chemotherapy`
meta$stage <- meta$Stage
cat("meta n:", nrow(meta), "\n")
com <- intersect(colnames(expr), meta$patient)
cat("expr x meta intersect:", length(com), "\n")
expr <- expr[, com]
meta <- meta[match(com, meta$patient), ]
meta$device <- colnames(expr)

## 3. 打分（与 TCGA 同参数）
score_sets <- c(modules, list(ICC45 = icc45),
                list(pyroptosis = death_sets_clean$pyroptosis,
                     necroptosis = death_sets_clean$necroptosis,
                     apoptosis = death_sets_clean$apoptosis,
                     ferroptosis = death_sets_clean$ferroptosis))
score_sets <- lapply(score_sets, function(g) intersect(g, rownames(expr)))
cat("gene coverage:\n"); print(sapply(score_sets, length))
gsva_res <- GSVA::gsva(GSVA::ssgseaParam(expr, score_sets))
scores <- t(as.matrix(gsva_res)); rm(gsva_res); gc()

## 4. 状态与解耦分数（同定义）
M1 <- scores[, "M1_immune_activation"]; M2 <- scores[, "M2_metabolic_GSH"]; M3 <- scores[, "M3_mito_cytoskeleton"]
state_strict <- M1 >= quantile(M1, 0.75) & M2 <= quantile(M2, 0.25) & M3 <= quantile(M3, 0.25)
state_relax  <- M1 >= quantile(M1, 0.60) & M2 <= quantile(M2, 0.40) & M3 <= quantile(M3, 0.40)
decouple <- as.numeric(scale(M1)) - (as.numeric(scale(M2)) + as.numeric(scale(M3))) / 2
dq <- cut(decouple, quantile(decouple, c(0, 0.25, 0.75, 1)), labels = c("Low", "Mid", "High"), include.lowest = TRUE)
res <- data.frame(sample = colnames(expr), state_high = state_strict, state_relax = state_relax,
                  decouple = decouple, decouple_tier = dq)
for (cn in colnames(scores)) res[[cn]] <- scores[, cn]
res <- merge(res, meta[, c("patient", "resp", "stage")], by.x = "sample", by.y = "patient", all.x = TRUE)

## 5. 关联：Resp vs NonResp（Wilcoxon, decouple + 各模块）、Fisher(严格态)、Stage Kruskal
cat("\n== KEY RESULTS ==\n")
cat("decoupled strict:", sum(res$state_high), "/", nrow(res), "\n")
w <- wilcox.test(decouple ~ resp, data = res)
cat("Wilcoxon decouple  Resp vs NonResp  p =", format(w$p.value, digits = 4), "\n")
for (m in names(modules)) {
  wm <- wilcox.test(res[[m]] ~ res$resp)
  cat(sprintf("  %-24s p=%.4g\n", m, wm$p.value))
}
tab <- table(res$state_high, res$resp)
fis <- tryCatch(fisher.test(tab), error = function(e) NA)
cat("Fisher strict x Resp  p =", ifelse(is.list(fis), format(fis$p.value, digits = 4), NA), "\n")
ks <- tryCatch(kruskal.test(decouple ~ stage, data = res[!is.na(res$stage), ]), error = function(e) NA)
cat("Kruskal decouple x Stage  p =", ifelse(is.list(ks), format(ks$p.value, digits = 4), NA), "\n")
cat("Resp median decouple: ", round(median(res$decouple[res$resp == "Responders"]), 3),
    " vs NonResp ", round(median(res$decouple[res$resp == "Non-responders"]), 3), "\n")

## 6. 判别器（同口径 ppcor）
if (requireNamespace("ppcor", quietly = TRUE)) {
  scd <- data.frame(mitoxy_total = rowMeans(scores[, names(modules)]), scores[, c("M1_immune_activation", "M4_lipid_perox_rupture", names(death_sets_clean))])
  pc <- ppcor::pcor(as.matrix(scd), method = "spearman")
  partials <- data.frame(death = names(death_sets_clean),
                         partial_r = pc$estimate["mitoxy_total", names(death_sets_clean)],
                         partial_p = pc$p.value["mitoxy_total", names(death_sets_clean)],
                         simple_r = sapply(names(death_sets_clean), function(d) cor(scd$mitoxy_total, scd[[d]], method = "spearman")))
  print(partials)
  fwrite(partials, file.path(direct, paste0(out_pref, "_partial.csv")))
}

## 7. 落盘
fwrite(res, file.path(direct, paste0(out_pref, "_scores.csv")))
saveRDS(list(wilcox = w$p.value, fisher = if (is.list(fis)) fis$p.value else NA, kruskal = ks$p.value,
             n = nrow(res)), file.path(direct, paste0(out_pref, "_stats.RDS")))

## 8. 图（NPG + Arial + 6.5pt 规范；单面板 + 组图；组图带面板字母）
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
npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")
res$resp2 <- factor(res$resp, levels = c("Responders", "Non-responders"))
p1 <- ggplot(res, aes(x = resp2, y = decouple, fill = resp2)) + geom_boxplot(width = .55, outlier.size = 0.5) +
  stat_compare_means(size = 2.2) +
  scale_fill_manual(values = c(Responders = "#00A087", `Non-responders` = "#DC0000"), guide = "none") +
  labs(title = "Decoupling score by platinum response", x = NULL, y = "Decoupling score")
p2 <- ggplot(res, aes(x = resp2, fill = resp2)) + geom_bar(stat = "count", position = "dodge") +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 2.4) +
  scale_fill_manual(values = c(Responders = "#00A087", `Non-responders` = "#DC0000")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Decoupled state frequency", x = NULL, y = "n") + theme(legend.position = "none")
p3 <- ggplot(res, aes(x = resp2, y = M1_immune_activation, fill = resp2)) + geom_boxplot(width = .55, outlier.size = 0.5) +
  stat_compare_means(size = 2.2) +
  scale_fill_manual(values = c(Responders = "#00A087", `Non-responders` = "#DC0000"), guide = "none") +
  labs(title = "M1 immune activation", x = NULL, y = "ssGSEA")
cp <- (p1 + p2 + p3) + plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 10, face = "bold"))
save_pub_r(p1, file.path(direct, paste0(out_pref, "_A_box")))
save_pub_r(p2, file.path(direct, paste0(out_pref, "_B_freq")))
save_pub_r(p3, file.path(direct, paste0(out_pref, "_C_m1")))
save_pub_r(cp, file.path(direct, paste0(out_pref, "_composite")), width_mm = 180, height_mm = 60)
cat("\nDONE GSE156699 external bulk\n")
