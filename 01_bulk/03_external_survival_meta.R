# ============================================================
# 03: P1 外部队列生存 meta（GSE13876/26193/14764/63885 + TCGA）
# 同口径评分 + 每队列 Cox(每 tier) → 合并森林图（逆方差固定效应）
# 附：GSE63885 铂敏感三级 × 解耦
# ============================================================
suppressMessages({
  library(data.table); library(GSVA); library(survival); library(ggplot2)
  library(AnnotationDbi); library(hgu133plus2.db); library(hgu133a2.db)
  library(patchwork); library(ggsci); library(ragg)
})
set.seed(20260831)
res_dir <- "D:/子宫内膜异位症2/analysis/L0_bulk/results/survival_meta"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
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

## ---------- 通用：读 series matrix 表达表 + 元数据（char 多行展开为 data.frame） ----------
parse_sm <- function(gz) {
  con <- gzfile(gz, "rt"); lines <- readLines(con, warn = FALSE); close(con)
  beg <- which(startsWith(lines, "!series_matrix_table_begin"))
  end <- which(startsWith(lines, "!series_matrix_table_end"))
  hdr <- strsplit(sub('"', "", lines[beg + 1]), "\t")[[1]]
  dat <- read.table(text = lines[(beg + 2):(end - 1)], sep = "\t", quote = "\"",
                    header = FALSE, row.names = 1, check.names = FALSE)
  colnames(dat) <- gsub('"', "", hdr[-1])
  ## 元数据：characteristics 多行 → 每行一个字段列（列名 = 值中第一个 ":" 前）
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
  list(expr = as.matrix(dat), char = as.data.frame(char, stringsAsFactors = FALSE,
                                                   check.names = FALSE))
}
## probe → symbol（db 或 自建 map）
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
score_cohort <- function(expr, meta) {
  if (max(expr) > 50) { expr <- log2(expr + 1); cat("  -> log2\n") }
  expr <- expr[!duplicated(rownames(expr)), ]
  cat("  genes:", nrow(expr), " samples:", ncol(expr), "\n")
  sc <- lapply(score_sets, function(g) intersect(g, rownames(expr)))
  cat("  coverage:", paste(sprintf("%s=%d", names(sc), sapply(sc, length)), collapse = " | "), "\n")
  g <- GSVA::gsva(GSVA::ssgseaParam(expr, sc))
  s <- t(as.matrix(g)); rm(g); gc()
  dec <- as.numeric(scale(s[, "M1_immune_activation"])) -
         (as.numeric(scale(s[, "M2_metabolic_GSH"])) + as.numeric(scale(s[, "M3_mito_cytoskeleton"]))) / 2
  out <- data.frame(sample = colnames(expr), decouple = dec, stringsAsFactors = FALSE)
  for (cn in colnames(s)) out[[cn]] <- s[, cn]
  ## 样本序号顺序与 meta 同（series matrix 列序 = GSM 顺序）
  out
}
cohort_tab <- data.frame()
tcga_res <- NULL

## ========== GSE26193（OS 年 + OS event, GPL570） ==========
d <- parse_sm("D:/子宫内膜异位症2/data_bulk/GSE26193_series_matrix.txt.gz")
e <- p2g_db(d$expr, hgu133plus2.db)
sc <- score_cohort(e, d$meta)
m <- data.frame(sample = colnames(e), stringsAsFactors = FALSE)
m$os_years <- as.numeric(d$char[["os time (years)"]])
m$os_event <- as.numeric(d$char[["os event"]])
dd <- merge(sc, m, by = "sample")
cat("\n== GSE26193 n=", nrow(dd), " events:", sum(dd$os_event == 1), " ==\n")
dd$os_days <- dd$os_years * 365.25
dd <- dd[!is.na(dd$os_days) & dd$os_days > 0, ]
dd$tier <- cut(dd$decouple, quantile(dd$decouple, c(0, .25, .75, 1)), labels = c("Low", "Mid", "High"), include.lowest = TRUE)
fx <- coxph(Surv(os_days, os_event) ~ as.numeric(tier), data = dd)
sm <- summary(fx); co <- sm$coefficients
lr <- survdiff(Surv(os_days, os_event) ~ tier, data = dd)
cohort_tab <- rbindlist(list(cohort_tab, data.table(cohort = "GSE26193", n = nrow(dd), events = sum(dd$os_event == 1),
  hr_tier = exp(co[1]), ci_lo = exp(co[1] - 1.96 * co[3]), ci_hi = exp(co[1] + 1.96 * co[3]),
  p_trend = co[5], logrank_p = pchisq(lr$chisq, length(unique(dd$tier)) - 1, lower.tail = FALSE), endpoint = "OS")), fill = TRUE)
fwrite(dd, file.path(res_dir, "scores_GSE26193.csv"))

## ========== GSE14764（OS 月 + event, GPL96 → hgu133a2 兜底） ==========
d <- parse_sm("D:/子宫内膜异位症2/data_bulk/GSE14764_series_matrix.txt.gz")
e <- p2g_db(d$expr, hgu133a2.db)
sc <- score_cohort(e, d$meta)
dd <- sc
dd$os_months <- as.numeric(d$char[["overall survival time"]])
dd$os_event <- as.numeric(d$char[["overall survival event"]])
cat("\n== GSE14764 n=", nrow(dd), " events:", sum(dd$os_event == 1), " ==\n")
dd$os_days <- dd$os_months * 30.44
dd <- dd[!is.na(dd$os_days) & dd$os_days > 0, ]
dd$tier <- cut(dd$decouple, quantile(dd$decouple, c(0, .25, .75, 1)), labels = c("Low", "Mid", "High"), include.lowest = TRUE)
fx <- coxph(Surv(os_days, os_event) ~ as.numeric(tier), data = dd)
co <- summary(fx)$coefficients
lr <- survdiff(Surv(os_days, os_event) ~ tier, data = dd)
cohort_tab <- rbindlist(list(cohort_tab, data.table(cohort = "GSE14764", n = nrow(dd), events = sum(dd$os_event == 1),
  hr_tier = exp(co[1]), ci_lo = exp(co[1] - 1.96 * co[3]), ci_hi = exp(co[1] + 1.96 * co[3]),
  p_trend = co[5], logrank_p = pchisq(lr$chisq, length(unique(dd$tier)) - 1, lower.tail = FALSE), endpoint = "OS")), fill = TRUE)
fwrite(dd, file.path(res_dir, "scores_GSE14764.csv"))

## ========== GSE13876（患者级：多样本/患者 → 均值折叠; OS 月） ==========
d <- parse_sm("D:/子宫内膜异位症2/data_bulk/GSE13876_series_matrix.txt.gz")
e <- p2g_map(d$expr, "D:/子宫内膜异位症2/data_omics/GSE13876_platform_table.csv.gz")
sc <- score_cohort(e, d$meta)
i_pt <- grep("patient id", d$meta$characteristics)
i_st <- grep("^status", d$meta$characteristics)
i_fu <- grep("fumnd", d$meta$characteristics)
dd <- sc
dd$patient <- d$char[["assigned unique patient id"]]
dd$status <- as.numeric(d$char[["status"]])
dd$fumnd <- as.numeric(d$char[["fumnd"]])
cat("\n== GSE13876 samples:", nrow(dd), " patients:", length(unique(dd$patient)), " ==\n")
## 患者级折叠（decouple 取均值；生存字段一致）
dd <- as.data.table(dd)
pat <- dd[, .(decouple = mean(decouple), status = status[1], fumnd = fumnd[1]), by = patient]
pat <- as.data.frame(pat)
pat <- pat[!is.na(pat$fumnd) & pat$status %in% c(0, 1), ]
cat("patients with survival:", nrow(pat), " events:", sum(pat$status == 1), "\n")
pat$os_days <- pat$fumnd * 30.44
pat$tier <- cut(pat$decouple, quantile(pat$decouple, c(0, .25, .75, 1)), labels = c("Low", "Mid", "High"), include.lowest = TRUE)
fx <- coxph(Surv(os_days, status) ~ as.numeric(tier), data = pat)
co <- summary(fx)$coefficients
lr <- survdiff(Surv(os_days, status) ~ tier, data = pat)
cohort_tab <- rbindlist(list(cohort_tab, data.table(cohort = "GSE13876", n = nrow(pat), events = sum(pat$status == 1),
  hr_tier = exp(co[1]), ci_lo = exp(co[1] - 1.96 * co[3]), ci_hi = exp(co[1] + 1.96 * co[3]),
  p_trend = co[5], logrank_p = pchisq(lr$chisq, length(unique(dd$tier)) - 1, lower.tail = FALSE), endpoint = "OS")), fill = TRUE)
fwrite(pat, file.path(res_dir, "scores_GSE13876_patient.csv"))

## ========== GSE63885（OS 天 + DOD; DFS; 铂敏感 3 级; BRCA1） ==========
d <- parse_sm("D:/子宫内膜异位症2/data_bulk/GSE63885_series_matrix.txt.gz")
e <- p2g_db(d$expr, hgu133plus2.db)
sc <- score_cohort(e, d$meta)
cn <- function(pat) { hit <- grep(pat, names(d$char), value = TRUE)[1]; if (is.na(hit)) stop("no col: ", pat); d$char[[hit]] }
dd <- sc
dd$os_days <- as.numeric(cn("os - overall survival"))
dd$dfs_days <- as.numeric(cn("dfs - disease-free survival"))
dd$status <- cn("clinical status at last follow-up")
dd$chemo <- cn("adjuwant chemotherapy")
dd$brca1 <- cn("brca1 mutation")
dd$tp53 <- cn("somatic tp53")
dd$plat <- cn("platinium sensitivity")
cat("\n== GSE63885 n=", nrow(dd), " ==\n")
dd$os_event <- ifelse(dd$status == "DOD", 1, 0)
dd$os_event[is.na(dd$status)] <- NA
dd2 <- dd[!is.na(dd$os_days) & dd$os_days > 0 & !is.na(dd$os_event), ]
cat("OS cohort:", nrow(dd2), " events:", sum(dd2$os_event == 1), "\n")
dd2$tier <- cut(dd2$decouple, quantile(dd2$decouple, c(0, .25, .75, 1)), labels = c("Low", "Mid", "High"), include.lowest = TRUE)
fx <- coxph(Surv(os_days, os_event) ~ as.numeric(tier), data = dd2)
co <- summary(fx)$coefficients
lr <- survdiff(Surv(os_days, os_event) ~ tier, data = dd2)
cohort_tab <- rbindlist(list(cohort_tab, data.table(cohort = "GSE63885", n = nrow(dd2), events = sum(dd2$os_event == 1),
  hr_tier = exp(co[1]), ci_lo = exp(co[1] - 1.96 * co[3]), ci_hi = exp(co[1] + 1.96 * co[3]),
  p_trend = co[5], logrank_p = pchisq(lr$chisq, length(unique(dd$tier)) - 1, lower.tail = FALSE), endpoint = "OS")), fill = TRUE)
## 铂敏感：resistant vs highly sensitive（DFS 定义无歧义；任选 DFS>0 者）
## 铂敏感：resistant vs highly sensitive（标签在 "): " 之后）
plat_lab <- sub(".*\\): ", "", dd$plat)
cat("plat labels:", paste(names(table(plat_lab)), collapse = " | "), "\n")
dd$plat2 <- ifelse(plat_lab == "resistant", "Resistant",
            ifelse(plat_lab == "highly sensitive", "HighSens",
            ifelse(plat_lab == "moderately sensitive", "ModSens", NA_character_)))
wplat <- tryCatch(wilcox.test(decouple ~ plat2, data = dd[dd$plat2 %in% c("Resistant", "HighSens"), ]), error = function(e) NULL)
cat("decouple Resistant vs HighlySensitive p =", if (!is.null(wplat)) format(wplat$p.value, digits = 3) else NA, "\n")
if (!is.null(wplat)) cat("  med R:", round(median(dd$decouple[dd$plat2 == "Resistant"]), 3),
                         " med HS:", round(median(dd$decouple[dd$plat2 == "HighSens"]), 3), "\n")
## 临床反应 CR vs SD/P
dd$outcome <- cn("clinical status post 1st line")
dd$outcome2 <- ifelse(dd$outcome == "CR", "CR", ifelse(dd$outcome %in% c("PR", "SD", "P"), "NonCR", NA))
wo <- tryCatch(wilcox.test(decouple ~ outcome2, data = dd[!is.na(dd$outcome2), ]), error = function(e) NULL)
cat("decouple CR vs NonCR p =", if (!is.null(wo)) format(wo$p.value, digits = 3) else NA, "\n")
fwrite(dd, file.path(res_dir, "scores_GSE63885.csv"))

## ========== TCGA（同口径实时计算，供森林图） ==========
tc <- fread("D:/子宫内膜异位症2/analysis/L0_bulk/results/track0_tcga_scores_states.csv", data.table = FALSE)
tc <- tc[!is.na(tc$os_days) & tc$os_days > 0, ]
tc$tier <- cut(tc$decouple, quantile(tc$decouple, c(0, .25, .75, 1)), labels = c("Low", "Mid", "High"), include.lowest = TRUE)
fxt <- coxph(Surv(os_days, event) ~ as.numeric(tier), data = tc)
cot <- summary(fxt)$coefficients
lrt <- survdiff(Surv(os_days, event) ~ tier, data = tc)
cohort_tab <- rbindlist(list(cohort_tab, data.table(cohort = "TCGA-OV", n = nrow(tc), events = sum(tc$event),
  hr_tier = exp(cot[1]), ci_lo = exp(cot[1] - 1.96 * cot[3]), ci_hi = exp(cot[1] + 1.96 * cot[3]),
  p_trend = cot[5], logrank_p = pchisq(lrt$chisq, length(unique(tc$tier)) - 1, lower.tail = FALSE),
  endpoint = "OS")), fill = TRUE)
cohort_tab <- as.data.frame(cohort_tab)

fwrite(cohort_tab, file.path(res_dir, "survival_meta_results.csv"))
print(cohort_tab)

## ---------- 森林图（逆方差固定效应合并；geom_segment 版避免 errorbarh 弃用问题） ----------
ct <- cohort_tab
ct$loghr <- log(ct$hr_tier); ct$se <- (log(ct$ci_hi) - log(ct$ci_lo)) / (2 * 1.96)
ct$w <- 1 / ct$se^2
fe <- sum(ct$loghr * ct$w) / sum(ct$w)
fe_se <- sqrt(1 / sum(ct$w))
fe_lo <- exp(fe - 1.96 * fe_se); fe_hi <- exp(fe + 1.96 * fe_se); fe_hr <- exp(fe)
z <- fe / fe_se; fe_p <- 2 * pnorm(-abs(z))
cat(sprintf("\nFIXED-EFFECT pooled HR = %.3f [%.3f-%.3f] p=%.4f (I2=%.1f%%)\n",
            fe_hr, fe_lo, fe_hi, fe_p,
            100 * max(0, (sum(ct$w * ct$loghr^2) - sum(ct$loghr)^2 / sum(ct$w)) / (length(ct$w) - 1) / (sum(ct$w) - sum(ct$w^2) / sum(ct$w)))))
ct$y <- seq_len(nrow(ct))
ct$lab <- sprintf("HR=%.2f [%.2f-%.2f]  p=%.3f", ct$hr_tier, ct$ci_lo, ct$ci_hi, ct$p_trend)
pool <- data.frame(y = nrow(ct) + 1, hr = fe_hr, lo = fe_lo, hi = fe_hi,
                   lab = sprintf("POOLED  HR=%.2f [%.2f-%.2f]  p=%.3f", fe_hr, fe_lo, fe_hi, fe_p))
ct$cohort <- factor(ct$cohort, levels = rev(unique(as.character(ct$cohort))))
ct$lab <- sprintf("HR=%.2f [%.2f-%.2f]  p=%.3f", ct$hr_tier, ct$ci_lo, ct$ci_hi, ct$p_trend)
ct$is_pool <- FALSE
pool <- data.frame(cohort = factor("Pooled (fixed)", levels = c(levels(ct$cohort), "Pooled (fixed)")),
                   hr_tier = fe_hr, ci_lo = fe_lo, ci_hi = fe_hi,
                   lab = sprintf("POOLED  HR=%.2f [%.2f-%.2f]  p=%.3f", fe_hr, fe_lo, fe_hi, fe_p),
                   is_pool = TRUE, stringsAsFactors = FALSE)
p_f2 <- ggplot() +
  geom_vline(xintercept = 1, linetype = 2, linewidth = .3) +
  geom_segment(data = ct, aes(x = ci_lo, xend = ci_hi, y = cohort, yend = cohort), linewidth = .35) +
  geom_point(data = ct, aes(x = hr_tier, y = cohort), size = 1.3, color = "#3C5488") +
  geom_text(data = ct, aes(x = ci_hi * 1.04, y = cohort, label = lab), hjust = 0, size = 2.3) +
  geom_segment(data = pool, aes(x = ci_lo, xend = ci_hi, y = cohort, yend = cohort), linewidth = .6, color = "#E64B35") +
  geom_point(data = pool, aes(x = hr_tier, y = cohort), shape = 18, size = 2.6, color = "#E64B35") +
  geom_text(data = pool, aes(x = ci_hi * 1.04, y = cohort, label = lab), hjust = 0, size = 2.5, color = "#E64B35", fontface = "bold") +
  scale_x_log10(limits = c(0.3, max(ct$ci_hi) * 3.6), breaks = c(0.4, 0.6, 1, 1.5, 2.5, 4)) +
  labs(title = "Decoupling tier (Low \u2192 High) and overall survival",
       x = "HR per tier (95% CI)", y = NULL) +
  coord_cartesian(clip = "off") +
  theme(panel.grid = element_blank(), plot.margin = margin(4, 30, 4, 4),
        axis.text.y = element_text(size = 6.5, face = "bold"))
## ---------- 主题（森林图后置定义，避免覆盖主 theme） ----------
theme_set(theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        plot.title = element_text(size = 7, face = "bold"), panel.grid = element_blank()))
save_pub_r <- function(p, fn, width_mm = 180, height_mm = 70, dpi = 600) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(fn, ".svg"), width = w, height = h); print(p); dev.off()
  grDevices::cairo_pdf(paste0(fn, ".pdf"), width = w, height = h, family = "Arial"); print(p); dev.off()
  ragg::agg_tiff(paste0(fn, ".tiff"), width = w, height = h, units = "in", res = dpi); print(p); dev.off()
}
save_pub_r(p_f2, file.path(res_dir, "Fig_forest_decoupling_OS"), 180, 80)
cat("\nDONE 03 survival meta\n")
