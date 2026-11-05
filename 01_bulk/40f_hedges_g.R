
suppressMessages({library(data.table)})
hedges_g_ci <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  m1 <- mean(x); m2 <- mean(y)
  s1 <- sd(x); s2 <- sd(y)
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  g <- (m1 - m2) / sp
  J <- 1 - 3 / (4 * (n1 + n2) - 9)
  g <- g * J
  se <- sqrt((n1 + n2) / (n1 * n2) + g^2 / (2 * (n1 + n2)))
  c(g = as.numeric(g), lo = as.numeric(g - 1.96 * se), hi = as.numeric(g + 1.96 * se))
}
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
      k <- sub(":.*$", "", parts[1]); v <- sub("^[^:]+: ?", "", parts)
      char[[k]] <- v
    }
  }
  list(expr = as.matrix(dat), char = as.data.frame(char, stringsAsFactors = FALSE, check.names = FALSE))
}
p2g_db <- function(expr, db) {
  syms <- mapIds(db, keys = rownames(expr), keytype = "PROBEID", column = "SYMBOL", multiVals = "first")
  ok <- !is.na(syms); expr <- expr[ok, ]; rownames(expr) <- syms[ok]; expr
}
library(AnnotationDbi); library(hgu133plus2.db); library(hgu133a2.db)
sets <- readRDS("D:/子宫内膜异位症2/pipeline/mitoxy_sets.RDS")
mods <- sets$modules
score_c <- function(e, gset) {
  g <- intersect(gset, rownames(e))
  if (length(g) < 5) return(rep(NA, ncol(e)))
  GSVA::gsva(GSVA::ssgseaParam(e[, , drop = FALSE], list(g = g))) |> as.numeric()
}
res <- list()
run_cohort <- function(name, ez, char, grp_cols, tumor_vals) {
  e <- ez; 
  if (max(e) > 50) e <- log2(e + 1)
  e <- e[!duplicated(rownames(e)), ]
  grp <- rep(NA, ncol(e))
  for (gc_ in grp_cols) {
    if (gc_ %in% names(char)) {
      v <- tolower(as.character(char[[gc_]]))
      grp[grepl("tumor|cancer|carcinoma|malignant|tumour", v)] <- "Tumor"
      grp[grepl("normal|benign|control|surface", v)] <- "Normal"
    }
  }
  if (sum(grp == "Tumor", na.rm = TRUE) < 5 | sum(grp == "Normal", na.rm = TRUE) < 5) {
    cat(name, "groups insufficient:", table(grp), "\n"); return(invisible())
  }
  cat(name, "n tumor/normal:", sum(grp == "Tumor"), sum(grp == "Normal"), "\n")
  for (m in names(mods)) {
    s <- score_c(e, mods[[m]])
    d <- data.frame(grp = grp, s = s); d <- d[!is.na(d$s), ]
    if (length(unique(d$grp)) < 2) next
    gv <- hedges_g_ci(d$s[d$grp == "Tumor"], d$s[d$grp == "Normal"])
    res[[paste(name, m)]] <<- data.frame(cohort = name, module = m, g = gv["g"], 
                                         ci_lo = gv["lo"], ci_hi = gv["hi"],
                                         n_t = sum(d$grp == "Tumor"), n_n = sum(d$grp == "Normal"))
  }
}
libs <- list(
  GSE26712 = list("D:/子宫内膜异位症2/data_bulk/GSE26712_series_matrix.txt.gz", hgu133plus2.db, c("tissue")),
  GSE32062 = list("D:/子宫内膜异位症2/data_bulk/GSE32062-GPL570_series_matrix.txt.gz", hgu133plus2.db, c("tissue", "diagnosis")))
for (nm in names(libs)) {
  pdf <- libs[[nm]]; cat("==", nm, "==\n")
  d <- parse_sm(pdf[[1]]); e <- p2g_db(d$expr, pdf[[2]])
  print(names(d$char))
  run_cohort(nm, e, d$char, pdf[[3]], NULL)
}
RR <- rbindlist(res, fill = TRUE)
fwrite(RR, "D:/子宫内膜异位症2/analysis/L0_bulk/results/survival_meta/module_effects_hedgesg.csv")
print(RR)
cat("DONE 40f\n")
