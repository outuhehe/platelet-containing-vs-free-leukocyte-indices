# =============================================================================
# 06_lymphocyte_vs_nlr.R
# Paired DeLong comparison of NLR vs the lymphocyte-alone positive control;
# cross-checks both AUCs against the file written by 05.
# Reads interim/cbc_export_v2.csv and tables/lymphocyte_auc_positive_control.csv;
# writes tables/lymphocyte_vs_nlr_paired_delong.csv. Uses data.table instead of
# readr; statistical logic is unchanged.
# =============================================================================
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
Sys.setenv(RA_ROOT = ra_root)
source(file.path(ra_root, "R", "00_common.R"))
load_analysis_pkgs("pROC")

d <- fread(file.path(P_INTERIM, "cbc_export_v2.csv"),
           check.names = FALSE)
stopifnot(nrow(d) == 879L,
          sum(as.integer(d$hospital_expire_flag)) == 87L)
alc_cols <- which(names(d) == "alc")
if (length(alc_cols) != 2L) stop("Expected duplicated alc columns in export")
alc_f <- d[[alc_cols[[1L]]]]
if (anyNA(alc_f) || any(alc_f <= 0)) stop("alc_f contains missing or nonpositive values")

y <- as.integer(d$hospital_expire_flag)
r_lym <- pROC::roc(y, -log(alc_f), quiet = TRUE, direction = "<")
r_nlr <- pROC::roc(y, log(d$nlr), quiet = TRUE, direction = "<")
auc_lym <- as.numeric(pROC::auc(r_lym))
auc_nlr <- as.numeric(pROC::auc(r_nlr))
delta <- auc_nlr - auc_lym
tt <- pROC::roc.test(r_nlr, r_lym, paired = TRUE, method = "delong")

previous <- fread(file.path(P_TABLES, "lymphocyte_auc_positive_control.csv"))
old_lym <- previous$AUC[previous$metric == "lymphocyte_absolute"]
old_nlr <- previous$AUC[previous$metric == "NLR"]
if (length(old_lym) != 1L || length(old_nlr) != 1L ||
    abs(auc_lym-old_lym) > 1e-12 || abs(auc_nlr-old_nlr) > 1e-12) {
  stop("AUCs do not match the previously reported values; stop and audit")
}

out <- data.frame(n = nrow(d), events = sum(y), AUC_lymphocyte = auc_lym,
                  AUC_NLR = auc_nlr, delta_AUC_NLR_minus_lymphocyte = delta,
                  DeLong_p = unname(tt$p.value), paired = TRUE)
fwrite(out, file.path(P_TABLES, "lymphocyte_vs_nlr_paired_delong.csv"))

cat("n =", nrow(d), "events =", sum(y), "\n")
cat("Lymphocyte AUC =", round(auc_lym,3), "; NLR AUC =", round(auc_nlr,3), "\n")
cat("Delta AUC (NLR - lymphocyte) =", round(delta,3), "\n")
cat("Paired DeLong p =", format.pval(tt$p.value, digits = 3), "\n")
cat(if (tt$p.value < .05) "Statistically significant difference.\n"
    else "No statistically significant difference.\n")
