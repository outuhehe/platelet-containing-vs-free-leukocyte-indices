# =============================================================================
# 11_complete_submission_tables.R
# Assemble the paper-format Table 2 (single AUCs) and Table 3 (paired AUCs with
# DeLong/BH and bootstrap CIs), with hard cross-checks against the frozen deltas,
# plus the all-absolute subgroup table.
# Reads interim/cbc_export_v2.csv and tables produced by 03/04; writes to tables/.
# =============================================================================
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
Sys.setenv(RA_ROOT = ra_root)
source(file.path(ra_root, "R", "00_common.R"))
load_analysis_pkgs("pROC")

d <- fread(file.path(P_INTERIM, "cbc_export_v2.csv"))
t2_source <- fread(file.path(P_TABLES, "table2_single_auc.csv"))
main_prev <- fread(file.path(P_TABLES, "table3_paired_delta_v2.csv"))
main_boot <- fread(file.path(P_TABLES, "table3b_bootstrap_ci_v2.csv"))
stopifnot(nrow(d) == 879L, uniqueN(d$hadm_id) == 879L,
          sum(as.integer(d$hospital_expire_flag)) == 87L,
          "all_abs" %in% names(d))

metric_order <- c("NLR", "SII", "PIV", "SIRI", "PLR", "neutrophil (positive control)")
stopifnot(nrow(t2_source) == 6L, all(metric_order %in% t2_source$metric))
table2 <- t2_source[match(metric_order, metric),
  .(metric = metric, AUC = AUC,
    CI95 = sprintf("%.3f - %.3f", CI95_lower, CI95_upper), n, events)]
fwrite(table2, file.path(P_TABLES, "Table2_single_auc_paper.csv"))

y <- as.integer(d$hospital_expire_flag)
as_roc <- function(outcome, value) {
  pROC::roc(response = factor(outcome, levels = c(0,1)),
            predictor = log(value), levels = c("0","1"),
            direction = "<", quiet = TRUE)
}
pair_specs <- list(
  list(name = "SII vs NLR", with_platelet = "sii", control = "nlr"),
  list(name = "PIV vs SIRI", with_platelet = "piv", control = "siri")
)

main_rows <- lapply(pair_specs, function(p) {
  auc_a <- as.numeric(pROC::auc(as_roc(y, d[[p$with_platelet]])))
  auc_b <- as.numeric(pROC::auc(as_roc(y, d[[p$control]])))
  data.table(comparison = p$name, AUC_with_platelet = auc_a, AUC_control = auc_b,
             dAUC = auc_a - auc_b, n = nrow(d), events = sum(y))
})
main_table <- rbindlist(main_rows)
expected <- c(-0.108890630442355, -0.0874404969242428)
if (any(abs(main_table$dAUC - expected) > 1e-10) ||
    any(abs(main_table$dAUC - main_prev$delta_AUC) > 1e-10)) {
  stop(sprintf("Main analysis mismatch: SII-NLR=%.12f, PIV-SIRI=%.12f",
               main_table$dAUC[[1L]], main_table$dAUC[[2L]]))
}
main_table[, `:=`(direction = ifelse(dAUC < 0, "negative", "positive"),
                  p_DeLong = main_prev$DeLong_p,
                  p_BH = main_prev$p_BH)]
main_table[, `:=`(
  CI95_lower = main_boot$CI95_lower[match(comparison, main_boot$comparison)],
  CI95_upper = main_boot$CI95_upper[match(comparison, main_boot$comparison)])]
fwrite(main_table, file.path(P_TABLES, "Table3_paired_auc_paper.csv"))

sub <- d[all_abs %in% c(TRUE, 1L)]
ys <- as.integer(sub$hospital_expire_flag)
sub_rows <- lapply(pair_specs, function(p) {
  r_a <- as_roc(ys, sub[[p$with_platelet]])
  r_b <- as_roc(ys, sub[[p$control]])
  test <- pROC::roc.test(r_b, r_a, paired = TRUE, method = "delong")
  auc_a <- as.numeric(pROC::auc(r_a))
  auc_b <- as.numeric(pROC::auc(r_b))
  data.table(comparison = p$name, n = nrow(sub), events = sum(ys),
             AUC_with_platelet = auc_a, AUC_control = auc_b,
             dAUC = auc_a-auc_b,
             direction = if (auc_a-auc_b < 0) "negative" else "positive",
             p_DeLong = unname(test$p.value))
})
sub_table <- rbindlist(sub_rows)
fwrite(sub_table, file.path(P_TABLES, "tableS_all_abs_subgroup.csv"))

cat("Main cohort: n =", nrow(d), "events =", sum(y), "\n")
cat("Main pair cross-check passed: deltas =",
    paste(sprintf("%.6f", main_table$dAUC), collapse = ", "), "\n")
cat("All-absolute subgroup: n =", nrow(sub), "events =", sum(ys), "\n")
print(table2)
print(main_table)
