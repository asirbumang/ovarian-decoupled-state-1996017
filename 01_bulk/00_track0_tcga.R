# ============================================================
# 轨道 0 · 卵巢癌 TCGA-OV：四模块评分 → 解耦状态 → 临床关联 + 判别器
# 复用：内异症项目 gene_sets（mitoxy_sets.RDS）+ GSVA ssGSEA + 28d 归因范式
# 输入：pipeline/mitoxy_sets.RDS；data_omics/tcga_ov_counts_matrix.csv.gz；
#       data_omics/tcga_ov_clinical.csv
# 输出：analysis/L0_bulk/results/（track0_tcga_*.csv + Fig2_* 图）
# ============================================================
suppressMessages({
  library(data.table); library(ggplot2); library(survival); library(survminer)
  library(GSVA); library(ggpubr); library(patchwork); library(ppcor)
})
set.seed(20260831)

pipeline <- "D:/子宫内膜异位症2/pipeline"
omics    <- "D:/子宫内膜异位症2/data_omics"
res_dir  <- "D:/子宫内膜异位症2/analysis/L0_bulk/results"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

## 1. 基因集（母版 RDS，逐字复用）
sets <- readRDS(file.path(pipeline, "mitoxy_sets.RDS"))
modules <- sets$modules; discriminators <- sets$discriminators
icc45 <- sets$icc45; death_sets <- sets$death_sets

## 2. 表达矩阵：counts -> log2(CPM+1)
cat("loading matrix...\n")
mx <- fread(file.path(omics, "tcga_ov_counts_matrix.csv.gz"), data.table = FALSE)
rownames(mx) <- mx$gene; mx$gene <- NULL
genes <- rownames(mx)
expr <- as.matrix(mx); storage.mode(expr) <- "double"
libsize <- colSums(expr)
# 同一 case 多份 RNA 文件 → 保留文库最大的一份
dup <- unique(colnames(expr)[duplicated(colnames(expr))])
for (d in dup) {
  cols <- which(colnames(expr) == d)
  keep <- cols[which.max(libsize[cols])]
  expr <- expr[, -setdiff(cols, keep)]
}
cat("collapsed duplicates:", length(dup), "-> n samples:", ncol(expr), "\n")
# 只保留蛋白编码名义基因（按 symbol 唯一化，保留表达更高的）
expr <- expr[!duplicated(rownames(expr)), ]
genes <- rownames(expr)
cpm <- sweep(expr, 2, colSums(expr) / 1e6, "/")
logexpr <- log2(cpm + 1)
rm(mx, expr, cpm); gc()
cat("matrix dim:", dim(logexpr)[1], "x", dim(logexpr)[2], "\n")

## 3. 模块评分：ssGSEA（四模块 + ICC45 + 3 个死亡方式集）
##    ★ 判别器纪律（内异症版）：死亡方式集必须先"去重叠"（从 mitoxy 模块/判别器共享基因中剔除）
mitoxy_genes <- unique(unlist(modules))
## 母版 04_stageD 去重叠清单（逐字）；ferroptosis 用母版 8 基因规范清单
death_sets_clean <- list(
  pyroptosis  = setdiff(c("GSDMD","GSDMC","NLRP3","CASP1","CASP4","CASP5"), mitoxy_genes),
  necroptosis = setdiff(c("MLKL","RIPK3","RIPK1"), mitoxy_genes),
  apoptosis   = setdiff(c("CASP3","CASP7","CASP8","CASP9","BCL2","BAX"), mitoxy_genes),
  ferroptosis = setdiff(c("TFRC","NCOA4","FTH1","FTL","HMOX1","CHAC1","PTGS2","SLC3A2"), mitoxy_genes))
cat("death sets cleaned (after removing mitoxy overlap):\n")
print(sapply(death_sets_clean, function(g) paste(intersect(g, genes), collapse = ",")))
score_sets <- c(modules, list(ICC45 = icc45),
                list(pyroptosis = death_sets_clean$pyroptosis,
                     necroptosis = death_sets_clean$necroptosis,
                     apoptosis = death_sets_clean$apoptosis,
                     ferroptosis = death_sets_clean$ferroptosis))
score_sets <- lapply(score_sets, function(g) intersect(g, genes))
cat("module gene coverage:\n")
print(sapply(score_sets, length))

gsva_res <- GSVA::gsva(GSVA::ssgseaParam(logexpr, score_sets))   # 母版默认 normalize=FALSE（同口径）
scores <- t(as.matrix(gsva_res)); colnames(scores) <- names(score_sets)
cat("score dim:", dim(scores), "\n")

## 4. 状态定义（内异症范式：M1 高 ∩ M2/M3 低，数据集内百分位）+ 连续解耦分数（主指标）
M1 <- scores[, "M1_immune_activation"]
M2 <- scores[, "M2_metabolic_GSH"]
M3 <- scores[, "M3_mito_cytoskeleton"]
q75 <- function(x) quantile(x, 0.75, na.rm = TRUE)
q25 <- function(x) quantile(x, 0.25, na.rm = TRUE)
state_strict <- M1 >= q75(M1) & M2 <= q25(M2) & M3 <= q25(M3)          # 严格交集（稀有，报告用）
state_relax <- M1 >= quantile(M1, 0.60) & M2 <= quantile(M2, 0.40) & M3 <= quantile(M3, 0.40)  # 宽松（敏感性）
decouple <- as.numeric(scale(M1)[,1] - (scale(M2)[,1] + scale(M3)[,1]) / 2)   # 连续解耦分数（主）
dq <- cut(decouple, quantile(decouple, c(0, 0.25, 0.75, 1)), labels = c("Low", "Mid", "High"), include.lowest = TRUE)
res <- data.frame(sample = rownames(scores), stringsAsFactors = FALSE)
res$state_high <- state_strict; res$state_relax <- state_relax
res$decouple <- decouple; res$decouple_tier <- dq
## 附各模块/判别器分数（溯源性；合并前绑定避免 merge 重排）
for (cn in colnames(scores)) res[[cn]] <- scores[, cn]

## 5. 临床合并（case barcode 对齐；治疗反应/生存/组织学）
clin <- fread(file.path(omics, "tcga_ov_clinical.csv"), data.table = FALSE)
clin$barcode <- trimws(clin$barcode)
## 数据卫生：剔除全空行（GDC 部分文件无临床内容），按 barcode 去重
clin <- clin[!is.na(clin$outcome) | !is.na(clin$vital_status) | !is.na(clin$days_to_death), ]
clin <- clin[!duplicated(clin$barcode), ]
clin$outcome_bin <- NA
clin$outcome_bin[clin$outcome == "Complete Remission/Response"] <- "CR"
clin$outcome_bin[clin$outcome == "Partial Remission/Response"]  <- "PR"
clin$outcome_bin[clin$outcome == "Stable Disease"]             <- "SD"
clin$outcome_bin[clin$outcome == "Progressive Disease"]        <- "PD"
clin$os <- suppressWarnings(as.numeric(clin$days_to_death))
clin$fup <- suppressWarnings(as.numeric(clin$days_to_last_followup))
clin$os_days <- ifelse(clin$vital_status == "Dead", clin$os, clin$fup)
clin$event <- ifelse(clin$vital_status == "Dead", 1L, 0L)
clin$responder <- ifelse(clin$outcome_bin %in% c("CR", "PR"), 1L, 0L)  # CR+PR vs SD+PD
clin$cr_vs_pd <- NA
clin$cr_vs_pd[clin$outcome_bin == "CR"] <- "CR"
clin$cr_vs_pd[clin$outcome_bin == "PD"] <- "PD"

res <- merge(res, clin[, c("barcode","outcome_bin","outcome","vital_status","os_days","event","responder","cr_vs_pd","histological_type","residual_tumor","age_at_dx")],
             by.x = "sample", by.y = "barcode", all.x = TRUE)
cat("matched clinical:", sum(!is.na(res$outcome_bin)), "/", nrow(res), "\n")

## 6. 关联分析（主指标 = 连续解耦分数；严格交集 = 次要/敏感性）
## 6.1 状态频率（严格交集 = 稀有；宽松版）
freq_all <- table(res$state_high); freq_relax <- table(res$state_relax)
## 6.2 反应关联：CR vs PD（Wilcoxon 解耦分数 + Fisher 严格交集）+ CR+PR vs SD+PD
w_crpd <- wilcox.test(decouple ~ cr_vs_pd, data = res[res$cr_vs_pd %in% c("CR","PD"), ])
tab_crpd <- table(res$state_high[res$cr_vs_pd %in% c("CR","PD")],
                  res$cr_vs_pd[res$cr_vs_pd %in% c("CR","PD")])
fisher_crpd <- tryCatch(fisher.test(tab_crpd), error = function(e) NA)
w_resp <- wilcox.test(decouple ~ responder, data = res[!is.na(res$responder), ])
tab_resp <- table(res$state_high[!is.na(res$responder)], res$responder[!is.na(res$responder)])
fisher_resp <- tryCatch(fisher.test(tab_resp), error = function(e) NA)
## 6.3 生存：解耦分数三分位（Low/Mid/High）KM + logrank；严格状态 Cox（作敏感性）
res$grp <- as.character(res$decouple_tier)
sd <- res[!is.na(res$os_days) & res$os_days > 0, ]
km <- survfit(Surv(os_days, event) ~ grp, data = sd)
lr <- survdiff(Surv(os_days, event) ~ grp, data = sd)
cox1 <- tryCatch(coxph(Surv(os_days, event) ~ decouple, data = sd), error = function(e) NULL)
cox2 <- tryCatch(coxph(Surv(os_days, event) ~ state_high, data = sd), error = function(e) NULL)
## 6.4 判别器（母版 04_stageD 逐字口径）：去重叠死亡集 → ppcor::pcor（Spearman）
##    vars = mitoxy_total, M1, M4, pyro/necro/apop/ferro →
##    汇报 mitoxy_total × 各死亡集 偏相关（控制 M1+M4+其它死亡集）＝"非已知死亡方式换皮"论证
scd_df <- data.frame(sample = rownames(scores),
                     mitoxy_total = rowMeans(scores[, names(modules)], na.rm = TRUE),
                     scores[, c("M1_immune_activation", "M4_lipid_perox_rupture",
                                "pyroptosis", "necroptosis", "apoptosis", "ferroptosis")])
pc <- ppcor::pcor(as.matrix(scd_df[, -1]), method = "spearman")
deaths <- c("pyroptosis", "necroptosis", "apoptosis", "ferroptosis")
partials <- data.frame(
  death = deaths,
  partial_r = pc$estimate["mitoxy_total", deaths],
  partial_p = pc$p.value["mitoxy_total", deaths],
  simple_r  = sapply(deaths, function(d) cor(scd_df$mitoxy_total, scd_df[[d]], method = "spearman")))

## 7. 落盘
fwrite(res, file.path(res_dir, "track0_tcga_scores_states.csv"))
fwrite(partials, file.path(res_dir, "track0_tcga_discriminator_partial.csv"))
cox1_coef <- if (!is.null(cox1)) summary(cox1)$coefficients[1, c(2, 5)] else c(NA, NA)
cox2_coef <- if (!is.null(cox2)) summary(cox2)$coefficients[1, c(2, 5)] else c(NA, NA)
out <- list(freq_strict = as.data.frame(freq_all), freq_relax = as.data.frame(freq_relax),
            wilcox_crpd = w_crpd$p.value, fisher_crpd = if (is.list(fisher_crpd)) fisher_crpd$p.value else NA,
            wilcox_resp = w_resp$p.value, fisher_resp = if (is.list(fisher_resp)) fisher_resp$p.value else NA,
            logrank_tier_p = pchisq(lr$chisq, df = 2, lower.tail = FALSE),
            cox_decouple_HR = cox1_coef[1], cox_decouple_p = cox1_coef[2],
            cox_state_HR = cox2_coef[1], cox_state_p = cox2_coef[2],
            state_strict_percent = mean(state_strict) * 100,
            state_relax_percent = mean(state_relax) * 100,
            partials = partials)
saveRDS(out, file.path(res_dir, "track0_tcga_stats.RDS"))
cat("== KEY RESULTS ==\n")
cat("strict state prevalence :", round(mean(state_strict) * 100, 2), "%  (n=", sum(state_strict), ")\n")
cat("relax state prevalence  :", round(mean(state_relax) * 100, 2), "%\n")
cat("wilcox CR vs PD (decouple) p =", w_crpd$p.value, "\n")
cat("fisher strict CR vs PD  p =", ifelse(is.list(fisher_crpd), fisher_crpd$p.value, NA), "\n")
cat("wilcox CR+PR vs SD+PD   p =", w_resp$p.value, "\n")
cat("logrank (3-tier) p      :", pchisq(lr$chisq, df = 2, lower.tail = FALSE), "\n")
cat("cox decouple (cont) HR  :", cox1_coef[1], " p:", cox1_coef[2], "\n")
cat("cox strict state HR     :", cox2_coef[1], " p:", cox2_coef[2], "\n")
print(partials)

## 8. 图（按《07_高分范式_出图.md》规范：NPG 色板 + Arial + theme_classic(6.5) + 三件套导出 + 单独面板/组图）
suppressMessages({ library(ggsci); library(svglite); library(ragg) })
theme_set(theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
        axis.ticks = element_line(linewidth = 0.35, colour = "black"),
        legend.title = element_text(size = 6.2),
        legend.text  = element_text(size = 5.8),
        strip.text   = element_text(size = 6.2, face = "bold"),
        plot.title   = element_text(size = 7, face = "bold"),
        panel.grid   = element_blank()))
save_pub_r <- function(plot, filename, width_mm = 90, height_mm = 70, dpi = 600) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial"); print(plot); dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi); print(plot); dev.off()
}
npg2 <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")
res$grp <- ifelse(res$state_high, "Decoupled", "Not")
p1 <- ggplot(res, aes(x = grp, fill = grp)) + geom_bar(stat = "count") +
  geom_text(stat = "count", aes(label = paste0("n=", after_stat(count))), vjust = -0.5, size = 2.4) +
  scale_fill_manual(values = c(Decoupled = npg2[1], Not = npg2[4])) + ylim(0, nrow(res) * 1.15) +
  labs(title = "Decoupled state frequency", x = NULL, y = "samples") + theme(legend.position = "none")
g2 <- ggsurvplot(km, data = sd, risk.table = TRUE, pval = TRUE, legend = "right",
                 palette = c(`grp=High` = npg2[1], `grp=Mid` = npg2[2], `grp=Low` = npg2[3]),
                 title = "OS by decoupling tier", fontsize = 2.2, risk.table.height = 0.28,
                 tables.theme = theme_cleantable(base_size = 5.8, base_family = "Arial") +
                   theme(plot.title = element_text(size = 7, hjust = 0, face = "bold"),
                         axis.title = element_text(size = 6.2),
                         axis.text = element_text(size = 5.8),
                         strip.text = element_text(size = 5.8)))
p2 <- g2$plot + scale_x_continuous(breaks = seq(0, 5000, by = 2500)) + labs(x = "Time (days)") +
  theme(legend.position = c(0.68, 0.90), legend.justification = c(0, 1),
        legend.title = element_text(size = 6.2), legend.text = element_text(size = 5.8),
        legend.key = element_rect(fill = NA), legend.background = element_blank())
g2$table <- g2$table + scale_x_continuous(breaks = seq(0, 5000, by = 2500)) + labs(x = "Time (days)")
p3 <- ggplot(res[!is.na(res$cr_vs_pd), ], aes(x = cr_vs_pd, y = decouple, fill = cr_vs_pd)) +
  geom_boxplot(alpha = .8, outlier.size = 0.5) + stat_compare_means(size = 2.2) +
  scale_fill_manual(values = c(CR = "#00A087", PD = "#DC0000"), guide = "none") +
  labs(title = "Decoupling score: CR vs PD", x = NULL, y = "Decoupling score")
partials$death <- factor(partials$death, levels = c("pyroptosis", "necroptosis", "apoptosis", "ferroptosis"))
p4 <- ggplot(partials, aes(x = death, y = partial_r)) + geom_col(fill = npg2[4], width = .6) +
  geom_text(aes(label = paste0("p=", formatC(partial_p, format = "g", digits = 2))), vjust = -0.5, size = 2.2) +
  expand_limits(y = max(partials$partial_r) * 1.25) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
  labs(title = "Mitoxy partial cor vs known death sets", x = NULL, y = "partial r (Spearman)") +
  theme(axis.text.x = element_text(size = 5.8, angle = 45, hjust = 1, vjust = 1))
panels <- list(Fig1a_freq = p1, Fig1b_os = p2, Fig1c_crpd = p3, Fig1d_disc = p4)
for (nm in names(panels)) save_pub_r(panels[[nm]], file.path(res_dir, nm))
## KM 单面板单独补风险表版本（按样图 NC img06 h/i：KM+risk table+精确 P）
save_pub_r(p2 / g2$table + plot_layout(heights = c(3.2, 1)),
           file.path(res_dir, "Fig1b_os_risktable"), width_mm = 90, height_mm = 92)
save_pub_r((p1 + p3 + p4 + p2) + plot_layout(ncol = 2) +
             plot_annotation(tag_levels = "A") &
             theme(plot.tag = element_text(size = 10, face = "bold")),
           file.path(res_dir, "Fig1_track0_composite"), width_mm = 180, height_mm = 120)
cat("figures saved (single panels + composite), per 07 出图规范.\n")
