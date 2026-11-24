# Load shared paths/helpers; pROC is added for the ROC calculations.
.this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
.ra_root <- normalizePath(file.path(dirname(.this_file), "..", ".."), winslash = "/")
Sys.setenv(RA_ROOT = .ra_root)
suppressMessages(source(file.path(.ra_root, "R", "00_common.R")))
# Exploration can generate row-level records and observed ROC thresholds.
# Keep every exploratory artifact inside the restricted, git-ignored tree.
P_OUT <- file.path(P_INTERIM, "exploratory")
dir.create(P_OUT, recursive = TRUE, showWarnings = FALSE)
load_analysis_pkgs("pROC")

cat("Analysis source: outputs/analyze_eicu_sicu_external.R\n")
cat("Window: [-1440, 0] minutes; earliest valid panel per stay; patient-level main-stay roster; B=2000; seed=20261115\n")
cat("Filters: SICU first adult cohort, five CBC values from one offset, patient-level discharge death outcome; no imputation\n")

panels <- fread(file.path(P_OUT, "eicu_sicu_aligned_panels.csv"))
main_patient <- fread(file.path(P_OUT, "eicu_sicu_patient_index_roster.csv"))
value_cols <- c("polys_pct", "lymphs_pct", "monos_pct", "wbc_k", "platelet_k")
panels[, values_present := complete.cases(.SD), .SDcols = value_cols]
panels[, values_valid := values_present &
         polys_pct >= 0 & polys_pct <= 100 &
         lymphs_pct >= 0 & lymphs_pct <= 100 &
         monos_pct >= 0 & monos_pct <= 100 &
         wbc_k > 0 & platelet_k > 0]
panels[, `:=`(
  anc = wbc_k * polys_pct / 100,
  alc = wbc_k * lymphs_pct / 100,
  amc = wbc_k * monos_pct / 100
)]
panels[, values_valid := values_valid & is.finite(anc) & is.finite(alc) &
         is.finite(amc) & alc >= 0]

pre_candidates <- panels[labresultoffset >= -1440 & labresultoffset <= 0 & values_valid]
setorder(pre_candidates, patientunitstayid, labresultoffset)
pre_stay <- pre_candidates[, .SD[1L], by = patientunitstayid]
pre_stay[, `:=`(
  NLR = anc / alc,
  PLR = platelet_k / alc,
  SII = platelet_k * anc / alc,
  SIRI = anc * amc / alc,
  PIV = platelet_k * anc * amc / alc,
  neutrophil = anc,
  lymphocyte = alc
)]

pre_index <- pre_stay[main_patient[, .(uniquepid, patientunitstayid)],
                      on = .(uniquepid, patientunitstayid), nomatch = 0L]
pre_index[, `:=`(
  death = as.integer(tolower(trimws(hospitaldischargelocation)) == "death"),
  outcome_available = hospitaldischargeoffset > 0
)]
pre_analysis <- pre_index[alc > 0 & outcome_available == TRUE]
nonfinite_n <- sum(!is.finite(as.matrix(pre_analysis[, .(NLR, PLR, SII, SIRI, PIV,
                                                          neutrophil, lymphocyte)])))
if (nonfinite_n > 0L) stop("Non-finite index values remain in prewindow analysis")

n_pre <- nrow(pre_analysis)
events_pre <- sum(pre_analysis$death)
cat("\nSample-size check (prewindow): n =", n_pre, "; events =", events_pre,
    "; expected 1388 / 131\n")
cat("Prewindow valid panels before patient-stay matching =", nrow(pre_stay),
    "; matched patient stays =", nrow(pre_index),
    "; non-finite index rows excluded =", nonfinite_n, "\n")
if (n_pre != 1388L || events_pre != 131L) {
  cat("DIFFERENCE: n delta =", n_pre - 1388L,
      "; events delta =", events_pre - 131L, "\n")
  stop("Prewindow sample size does not match the required 1388/131; stopped")
}

roc_for <- function(outcome, score) {
  pROC::roc(response = factor(outcome, levels = c(0, 1)), predictor = score,
            levels = c("0", "1"), direction = ">", quiet = TRUE)
}
pair_specs <- list(
  list(name = "SII vs NLR", with_platelet = "SII", control = "NLR"),
  list(name = "PIV vs SIRI", with_platelet = "PIV", control = "SIRI")
)
pair_result <- function(dat, spec) {
  score_a <- -log1p(dat[[spec$with_platelet]])
  score_b <- -log1p(dat[[spec$control]])
  r_a <- roc_for(dat$death, score_a)
  r_b <- roc_for(dat$death, score_b)
  test <- pROC::roc.test(r_b, r_a, paired = TRUE, method = "delong")
  auc_a <- as.numeric(pROC::auc(r_a))
  auc_b <- as.numeric(pROC::auc(r_b))
  data.table(comparison = spec$name, AUC_a = auc_a, AUC_control = auc_b,
             delta_AUC = auc_a - auc_b,
             direction = if (auc_a - auc_b < 0) "negative" else "nonnegative",
             DeLong_p = unname(test$p.value), n = nrow(dat), events = sum(dat$death))
}

paired <- rbindlist(lapply(pair_specs, function(s) pair_result(pre_analysis, s)))

set.seed(20261115)
B <- 2000L
workers <- min(4L, max(1L, parallel::detectCores(logical = FALSE) - 1L))
idx0 <- which(pre_analysis$death == 0L)
idx1 <- which(pre_analysis$death == 1L)
score_a <- lapply(pair_specs, function(s) log1p(pre_analysis[[s$with_platelet]]))
score_b <- lapply(pair_specs, function(s) log1p(pre_analysis[[s$control]]))
auc_rank <- function(y, score) {
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
boot_worker <- function(ids, y, score_a, score_b, idx0, idx1) {
  out <- matrix(NA_real_, nrow = length(ids), ncol = length(score_a))
  for (i in seq_along(ids)) {
    ix <- c(sample(idx0, length(idx0), TRUE), sample(idx1, length(idx1), TRUE))
    for (j in seq_along(score_a)) {
      out[i, j] <- auc_rank(y[ix], score_a[[j]][ix]) -
        auc_rank(y[ix], score_b[[j]][ix])
    }
  }
  out
}
chunks <- split(seq_len(B), cut(seq_len(B), breaks = workers, labels = FALSE))
cl <- parallel::makeCluster(workers)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterSetRNGStream(cl, 20261115)
parallel::clusterExport(cl, c("pre_analysis", "score_a", "score_b", "idx0", "idx1",
                              "auc_rank", "boot_worker"), envir = environment())
boot_parts <- parallel::parLapply(cl, chunks, boot_worker,
                                  y = pre_analysis$death, score_a = score_a,
                                  score_b = score_b, idx0 = idx0, idx1 = idx1)
boot <- do.call(rbind, boot_parts)
boot_ci <- t(apply(boot, 2L, quantile, probs = c(.025, .975), na.rm = TRUE,
                   names = FALSE, type = 6))
paired[, `:=`(CI95_lower = boot_ci[, 1L], CI95_upper = boot_ci[, 2L],
              bootstrap_replicates = B)]

setcolorder(paired, c("comparison", "AUC_a", "AUC_control", "delta_AUC",
                      "direction", "DeLong_p", "n", "events", "CI95_lower",
                      "CI95_upper", "bootstrap_replicates"))
fwrite(paired, file.path(P_OUT, "eicu_prewindow_paired_delong_bootstrap.csv"))

cat("\nPrewindow paired DeLong + bootstrap:\n")
print(paired)
cat("\nSmall check table (single AUCs and delta AUC):\n")
print(paired[, .(comparison, AUC_a, AUC_control, delta_AUC, DeLong_p,
                 CI95_lower, CI95_upper)])
cat("\nwritten: outputs/eicu_prewindow_paired_delong_bootstrap.csv\n")
