# =============================================================================
# 04_valid_paired_comparisons.R
# Final paired DeLong tests (SII vs NLR, PIV vs SIRI) with 2000-replicate
# stratified bootstrap CIs and the PLR single-marker AUC, on the 879-row panel.
# Reads interim/cbc_export_v2.csv; writes the v2 paper tables to tables/.
# Statistical logic is unchanged from outputs/run_valid_paired_comparisons_v2.R.
# =============================================================================
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
Sys.setenv(RA_ROOT = ra_root)
source(file.path(ra_root, "R", "00_common.R"))
load_analysis_pkgs("pROC")

d <- fread(file.path(P_INTERIM, "cbc_export_v2.csv"))
stopifnot(nrow(d) == 879L, uniqueN(d$hadm_id) == 879L,
          sum(as.integer(d$hospital_expire_flag)) == 87L)
y <- as.integer(d$hospital_expire_flag)

as_roc <- function(score) {
  pROC::roc(response = factor(y, levels = c(0,1)), predictor = log(score),
            levels = c("0","1"), direction = "<", quiet = TRUE)
}
pair_specs <- list(
  list(name = "SII vs NLR", with_platelet = "sii", control = "nlr"),
  list(name = "PIV vs SIRI", with_platelet = "piv", control = "siri")
)

pair_rows <- lapply(pair_specs, function(p) {
  r_a <- as_roc(d[[p$with_platelet]])
  r_b <- as_roc(d[[p$control]])
  test <- pROC::roc.test(r_b, r_a, paired = TRUE, method = "delong")
  auc_a <- as.numeric(pROC::auc(r_a))
  auc_b <- as.numeric(pROC::auc(r_b))
  data.table(comparison = p$name, AUC_with_platelet = auc_a,
             AUC_control = auc_b, delta_AUC = auc_a - auc_b,
             direction = if (auc_a - auc_b < 0) "negative" else "nonnegative",
             DeLong_p = unname(test$p.value), n = nrow(d), events = sum(y))
})
table3 <- rbindlist(pair_rows)
table3[, p_BH := p.adjust(DeLong_p, method = "BH")]

auc_rank <- function(outcome, score) {
  n1 <- sum(outcome == 1L)
  n0 <- sum(outcome == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[outcome == 1L]) - n1*(n1+1)/2) / (n1*n0)
}

set.seed(20261115)
B <- 2000L
workers <- min(4L, max(1L, parallel::detectCores(logical = FALSE)-1L))
idx0 <- which(y == 0L)
idx1 <- which(y == 1L)
boot_worker <- function(ids, outcome, scores, idx0, idx1) {
  out <- matrix(NA_real_, nrow = length(ids), ncol = length(scores))
  for (i in seq_along(ids)) {
    s <- c(sample(idx0, length(idx0), replace = TRUE),
           sample(idx1, length(idx1), replace = TRUE))
    yy <- outcome[s]
    for (j in seq_along(scores)) {
      out[i,j] <- auc_rank(yy, scores[[j]]$a[s]) -
                  auc_rank(yy, scores[[j]]$b[s])
    }
  }
  out
}
scores <- lapply(pair_specs, function(p) {
  list(a = log(d[[p$with_platelet]]), b = log(d[[p$control]]))
})
chunks <- split(seq_len(B), cut(seq_len(B), breaks = workers, labels = FALSE))
cl <- parallel::makeCluster(workers)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterSetRNGStream(cl, 20261115)
parallel::clusterExport(cl, c("y", "scores", "idx0", "idx1", "auc_rank",
                               "boot_worker"), envir = environment())
boot_parts <- parallel::parLapply(cl, chunks, boot_worker,
                                  outcome = y, scores = scores,
                                  idx0 = idx0, idx1 = idx1)
boot <- do.call(rbind, boot_parts)
boot_ci <- t(apply(boot, 2L, quantile, probs = c(.025,.975), na.rm = TRUE,
                   names = FALSE, type = 6))
table_boot <- data.table(comparison = vapply(pair_specs, `[[`, "", "name"),
                         delta_AUC = table3$delta_AUC,
                         CI95_lower = boot_ci[,1L], CI95_upper = boot_ci[,2L],
                         bootstrap_replicates = B)

r_plr <- as_roc(d$plr)
plr_ci <- as.numeric(pROC::ci.auc(r_plr, method = "delong"))
plr_row <- data.table(metric = "PLR", AUC = as.numeric(pROC::auc(r_plr)),
                      CI95_lower = plr_ci[[1L]], CI95_upper = plr_ci[[3L]],
                      n = nrow(d), events = sum(y))

fwrite(table3, file.path(P_TABLES, "table3_paired_delta_v2.csv"))
fwrite(table_boot, file.path(P_TABLES, "table3b_bootstrap_ci_v2.csv"))
fwrite(plr_row, file.path(P_TABLES, "plr_single_auc_v2.csv"))

cat("n =", nrow(d), "events =", sum(y), "workers =", workers,
    "bootstrap replicates =", B, "\n")
print(table3)
print(table_boot)
cat("PLR single-marker AUC =", round(plr_row$AUC,3),
    "(95% CI", round(plr_ci[[1L]],3), "to", round(plr_ci[[3L]],3), ")\n")
cat("Negative pair count =", sum(table3$delta_AUC < 0),
    "; BH-FDR < 0.05 count =", sum(table3$p_BH < .05), "\n")
