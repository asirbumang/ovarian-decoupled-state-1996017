# ============================================================
# 26d: 免疫组成残差化解耦 —— 方向矛盾仲裁
# 问题: bulk decouple ~ 免疫构成? 残差解耦再关联临床(TCGA OS/反应, GSE156699 反应)
# 输入: track2/deconv_*.csv + 已有 scores csv
# ============================================================
suppressMessages({library(data.table); library(survival); library(ggplot2)
  library(ggpubr); library(survminer); library(patchwork); library(ggsci); library(ragg)})
set.seed(20260831)
tr <- "D:/子宫内膜异位症2/analysis/L0_bulk/results/track2"
out_dir <- file.path(tr, "figs"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
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

res_sum <- data.frame()
## ============ TCGA-OV ============
dv <- fread(file.path(tr, "deconv_TCGA_OV.csv"), data.table = FALSE)
sc <- fread("D:/子宫内膜异位症2/analysis/L0_bulk/results/track0_tcga_scores_states.csv", data.table = FALSE)
d <- merge(sc, dv[, c("sample", "immune_frac", "fit", "Epithelial", "Stromal", "Myeloid", "T_NK", "B_Plasma")], by = "sample")
cat("== TCGA: n =", nrow(d), " ==\n")
r2_imm <- summary(lm(decouple ~ immune_frac, data = d))$r.squared
rho <- suppressWarnings(cor(d$decouple, d$immune_frac, method = "spearman"))
cat(sprintf("decouple ~ immune_frac: R2=%.3f  rho=%.3f\n", r2_imm, rho))
d$dec_res <- resid(lm(decouple ~ immune_frac, data = d))
## OS：residual 三分位
sd <- d[!is.na(d$os_days) & d$os_days > 0, ]
sd$tier <- cut(sd$dec_res, quantile(sd$dec_res, c(0, .25, .75, 1)), labels = c("Low", "Mid", "High"), include.lowest = TRUE)
lr <- survdiff(Surv(os_days, event) ~ tier, data = sd)
lr_p <- pchisq(lr$chisq, 2, lower.tail = FALSE)
ct <- coxph(Surv(os_days, event) ~ as.numeric(tier), data = sd)
cat(sprintf("RESIDUAL OS: logrank p=%.4f ; trend HR=%.3f p=%.4f\n", lr_p, exp(coef(ct)), summary(ct)$coef[5]))
## 反应
w_res <- wilcox.test(dec_res ~ cr_vs_pd, data = d[d$cr_vs_pd %in% c("CR", "PD"), ])
cat(sprintf("RESIDUAL CR vs PD: p=%.4f (crude %.4f)\n", w_res$p.value, wilcox.test(decouple ~ cr_vs_pd, data = d[d$cr_vs_pd %in% c("CR", "PD"), ])$p.value))
w_resp <- wilcox.test(dec_res ~ responder, data = d[!is.na(d$responder), ])
cat(sprintf("RESIDUAL CR+PR vs SD+PD: p=%.4f\n", w_resp$p.value))
## 各模块 × immune_frac
cat("modules x immune_frac (Spearman):\n")
mod_r <- sapply(c("M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton", "M4_lipid_perox_rupture", "decouple"),
                function(m) suppressWarnings(cor(d[[m]], d$immune_frac, method = "spearman")))
print(round(mod_r, 3))
fwrite(d, file.path(tr, "TCGA_OV_residual.csv"))
res_sum <- rbind(res_sum, data.frame(cohort = "TCGA_OV", r2_imm = r2_imm, rho_imm = rho, metric = "OS-3tier logrank", crude_p = 0.09614245, residual_p = lr_p,
                                     crude_stat = NA, residual_stat = NA))

## 图（TCGA）：a scatter decouple~imm + b 残差KM + c CR/PD 残差箱线
p_a <- ggplot(d, aes(immune_frac, decouple)) + geom_point(size = .5, alpha = .5, color = npg[4]) +
  geom_smooth(method = "lm", se = TRUE, linewidth = .5, color = npg[1]) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = sprintf("R2=%.3f  rho=%.3f", r2_imm, rho), size = 2.4) +
  labs(title = "Decoupling vs immune fraction (TCGA)", x = "immune fraction (NNLS)", y = "decoupling")
g2 <- ggsurvplot(survfit(Surv(os_days, event) ~ tier, data = sd), data = sd, pval = TRUE, legend = "right",
                 palette = npg[c(1, 3, 2)], title = "Residual decoupling: OS", risk.table = FALSE,
                 legend.title = "", legend.labs = c("Low", "Mid", "High"))
p_b <- g2$plot + scale_x_continuous(breaks = seq(0, 5000, 2500)) +
  theme(legend.position = c(0.53, 0.92), legend.justification = c(0, 1),
        legend.key = element_rect(fill = NA), legend.key.size = unit(2.2, "mm"),
        legend.text = element_text(size = 5.6), legend.title = element_blank(),
        legend.background = element_blank(), plot.title = element_text(size = 6.8, face = "bold"))
p_c <- ggplot(d[d$cr_vs_pd %in% c("CR", "PD"), ], aes(cr_vs_pd, dec_res, fill = cr_vs_pd)) +
  geom_boxplot(width = .55, outlier.size = .5) + stat_compare_means(size = 2.2) +
  scale_fill_manual(values = c(CR = "#00A087", PD = "#DC0000"), guide = "none") +
  labs(title = "Residualized decoupling: CR vs PD", x = NULL, y = "decoupling (resid.)")
cp <- (p_a + p_c + p_b) + plot_layout(ncol = 3) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))
save_pub_r(p_a, file.path(out_dir, "TCGA_A_scatter"), 90, 70)
save_pub_r(p_c, file.path(out_dir, "TCGA_C_crpd"))
save_pub_r(p_b, file.path(out_dir, "TCGA_B_km"))
save_pub_r(cp, file.path(out_dir, "TCGA_composite"), 180, 60)

## ============ GSE156699 ============
dv2 <- fread(file.path(tr, "deconv_GSE156699.csv"), data.table = FALSE)
sc2 <- fread("D:/子宫内膜异位症2/analysis/L0_bulk/results/external/ext_GSE156699_scores.csv", data.table = FALSE)
d2 <- merge(sc2, dv2[, c("sample", "immune_frac", "fit")], by = "sample", all.x = TRUE)
cat("\n== GSE156699: n =", nrow(d2), " ==\n")
r2_imm2 <- summary(lm(decouple ~ immune_frac, data = d2))$r.squared
rho2 <- suppressWarnings(cor(d2$decouple, d2$immune_frac, method = "spearman"))
cat(sprintf("decouple ~ immune_frac: R2=%.3f rho=%.3f\n", r2_imm2, rho2))
d2$dec_res <- resid(lm(decouple ~ immune_frac, data = d2))
w2 <- wilcox.test(dec_res ~ resp, data = d2)
w2c <- wilcox.test(decouple ~ resp, data = d2)
cat(sprintf("RESIDUAL Resp vs NonResp: p=%.4f (crude %.4f)\n", w2$p.value, w2c$p.value))
cat("resp median res:", round(median(d2$dec_res[d2$resp == "Responders"]), 3),
    " vs ", round(median(d2$dec_res[d2$resp == "Non-responders"]), 3), "\n")
fwrite(d2, file.path(tr, "GSE156699_residual.csv"))
res_sum <- rbind(res_sum, data.frame(cohort = "GSE156699", r2_imm = r2_imm2, rho_imm = rho2, metric = "Resp vs NonResp", crude_p = 0.4115, residual_p = w2$p.value, crude_stat = NA, residual_stat = NA))
p_d <- ggplot(d2, aes(immune_frac, decouple)) + geom_point(size = .6, alpha = .6, color = npg[4]) +
  geom_smooth(method = "lm", se = TRUE, linewidth = .5, color = npg[1]) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, label = sprintf("R2=%.3f rho=%.3f", r2_imm2, rho2), size = 2.4) +
  labs(title = "Decoupling vs immune fraction (GSE156699)", x = "immune fraction", y = "decoupling")
p_e <- ggplot(d2, aes(resp, dec_res, fill = resp)) + geom_boxplot(width = .55, outlier.size = .5) +
  stat_compare_means(size = 2.2) +
  scale_fill_manual(values = c(Responders = "#00A087", `Non-responders` = "#DC0000"), guide = "none") +
  labs(title = "Residualized decoupling by response", x = NULL, y = "decoupling (resid.)")
cp2 <- (p_d + p_e) + plot_layout(ncol = 2) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))
save_pub_r(p_d, file.path(out_dir, "GSE156699_A_scatter"))
save_pub_r(p_e, file.path(out_dir, "GSE156699_B_resp"))
save_pub_r(cp2, file.path(out_dir, "GSE156699_composite"), 180, 60)

fwrite(res_sum, file.path(tr, "residual_assoc_summary.csv"))
cat("\n== SUMMARY ==\n"); print(res_sum)
cat("DONE 26d\n")
