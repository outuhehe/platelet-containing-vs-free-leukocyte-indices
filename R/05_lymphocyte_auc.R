# =============================================================================
# 05_lymphocyte_auc.R
# Lymphocyte-alone positive-control AUC (direction >), with neutrophil and NLR
# reference AUCs; cross-checks the ANC AUC against Table 2.
# Reads interim/cbc_export_v2.csv and tables/table2_single_auc.csv; writes
# tables/lymphocyte_auc_positive_control.csv. Logic unchanged from the original.
# =============================================================================
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
Sys.setenv(RA_ROOT = ra_root)
source(file.path(ra_root, "R", "00_common.R"))
load_analysis_pkgs("pROC")

d <- fread(file.path(P_INTERIM, "cbc_export_v2.csv"), check.names = FALSE)
stopifnot(nrow(d) == 879L, uniqueN(d$hadm_id) == 879L,
          sum(as.integer(d$hospital_expire_flag)) == 87L)

# The export aliases alc_f and raw alc to the same header; the first is alc_f.
alc_cols <- which(names(d) == "alc")
if (length(alc_cols) != 2L) stop("Expected resolved and raw alc columns in export")
alc_f <- d[[alc_cols[[1L]]]]
alc_zero_n <- sum(alc_f == 0, na.rm = TRUE)
alc_na_n <- sum(is.na(alc_f))
use_log <- alc_zero_n == 0L
lymph_score <- if (use_log) log(alc_f) else alc_f

as_roc <- function(y, score, direction) {
  pROC::roc(response = factor(y, levels = c(0,1)), predictor = score,
            levels = c("0","1"), direction = direction, quiet = TRUE)
}
y <- as.integer(d$hospital_expire_flag)
r_lym <- as_roc(y, lymph_score, ">")
ci_lym <- as.numeric(pROC::ci.auc(r_lym, method = "delong"))
r_anc <- as_roc(y, log(d$anc_single), "<")
ci_anc <- as.numeric(pROC::ci.auc(r_anc, method = "delong"))
r_nlr <- as_roc(y, log(d$nlr), "<")
nlr_auc <- as.numeric(pROC::auc(r_nlr))

t2 <- fread(file.path(P_TABLES, "table2_single_auc.csv"))
anc_table_auc <- t2[metric == "neutrophil (positive control)", AUC]
if (length(anc_table_auc) != 1L ||
    abs(as.numeric(pROC::auc(r_anc)) - anc_table_auc) > 1e-12) {
  stop("ANC AUC does not match the existing Table 2; stop for data audit")
}

out <- data.table(
  metric = c("lymphocyte_absolute", "neutrophil_absolute", "NLR"),
  AUC = c(as.numeric(pROC::auc(r_lym)), as.numeric(pROC::auc(r_anc)), nlr_auc),
  CI95_lower = c(ci_lym[[1L]], ci_anc[[1L]], NA_real_),
  CI95_upper = c(ci_lym[[3L]], ci_anc[[3L]], NA_real_),
  direction = c(">", "<", "<"),
  n = c(sum(complete.cases(y, lymph_score)), nrow(d), nrow(d)),
  events = c(sum(y[complete.cases(y, lymph_score)]), sum(y), sum(y))
)
fwrite(out, file.path(P_TABLES, "lymphocyte_auc_positive_control.csv"))

cat("Rows/events =", nrow(d), "/", sum(y), "\n")
cat("alc_f zero count =", alc_zero_n, "; NA count =", alc_na_n,
    "; use log =", use_log, "\n")
cat("Lymphocyte-alone AUC =", round(out$AUC[[1L]],3),
    "95% CI", sprintf("%.3f-%.3f", ci_lym[[1L]],ci_lym[[3L]]),
    "; direction >\n")
cat("NLR AUC =", round(nlr_auc,3), "\n")
print(out)
