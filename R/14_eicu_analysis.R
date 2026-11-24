# =============================================================================
# 14_eicu_analysis.R
# Direction-consistency analysis on the eICU SICU panels: numeric QC, index
# derivation, one index stay per patient, single AUCs, paired DeLong with 2000
# patient-stratified bootstrap CIs, pre-admission-window sensitivity, and
# comparison against the MIMIC point estimates. No claim of external validation.
# Reads interim/eicu_sicu_aligned_panels.csv plus MIMIC tables; writes tables/.
# =============================================================================
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
Sys.setenv(RA_ROOT = ra_root)
source(file.path(ra_root, "R", "00_common.R"))
load_analysis_pkgs("pROC")

panels <- fread(file.path(P_INTERIM, "eicu_sicu_aligned_panels.csv"))
stopifnot(uniqueN(panels$patientunitstayid) < nrow(panels))
setorder(panels, patientunitstayid, labresultoffset)

value_cols <- c("polys_pct","lymphs_pct","monos_pct","wbc_k","platelet_k")
qc <- rbindlist(lapply(value_cols, function(v) {
  x <- panels[[v]]
  data.table(variable = v, n_panel_rows = nrow(panels), n_missing = sum(is.na(x)),
             n_zero = sum(x == 0, na.rm = TRUE), n_negative = sum(x < 0, na.rm = TRUE),
             min = min(x, na.rm = TRUE), median = median(x, na.rm = TRUE),
             max = max(x, na.rm = TRUE))
}))
fwrite(qc, file.path(P_TABLES, "eicu_sicu_lab_value_qc.csv"))

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

main_stay <- panels[labresultoffset >= 0 & labresultoffset <= 1440 & values_valid]
setorder(main_stay, patientunitstayid, labresultoffset)
main_stay <- main_stay[, .SD[1L], by = patientunitstayid]
pre_stay <- panels[labresultoffset >= -1440 & labresultoffset <= 0 & values_valid]
setorder(pre_stay, patientunitstayid, labresultoffset)
pre_stay <- pre_stay[, .SD[1L], by = patientunitstayid]

derive_indices <- function(x) {
  x[, `:=`(
    NLR = anc / alc,
    PLR = platelet_k / alc,
    SII = platelet_k * anc / alc,
    SIRI = anc * amc / alc,
    PIV = platelet_k * anc * amc / alc,
    neutrophil = anc,
    lymphocyte = alc
  )]
  x
}
main_stay <- derive_indices(main_stay)
pre_stay <- derive_indices(pre_stay)

setorder(main_stay, uniquepid, hospitaldischargeyear, patientunitstayid)
main_patient <- main_stay[, .SD[1L], by = uniquepid]
main_patient[, death := as.integer(tolower(trimws(hospitaldischargelocation)) == "death")]
main_patient[, outcome_available := hospitaldischargeoffset > 0]

pre_index <- pre_stay[main_patient[, .(uniquepid, patientunitstayid)],
                      on = .(uniquepid, patientunitstayid), nomatch = 0L]
pre_index[, `:=`(
  death = as.integer(tolower(trimws(hospitaldischargelocation)) == "death"),
  outcome_available = hospitaldischargeoffset > 0)]
setorder(pre_index, uniquepid)

excluded_zero_lymph <- sum(main_patient$alc <= 0)
analysis <- main_patient[alc > 0 & outcome_available == TRUE]
if (any(!is.finite(as.matrix(analysis[, .(NLR,PLR,SII,SIRI,PIV,neutrophil,lymphocyte)])))) {
  stop("Non-finite index value remains after validity filtering")
}
if (uniqueN(analysis$uniquepid) != nrow(analysis)) stop("Patient-level index duplicated")

metric_names <- c("NLR","PLR","SII","SIRI","PIV","neutrophil","lymphocyte")
auc_direction <- function(v) {
  if (v == "lymphocyte") log1p(analysis[[v]]) else -log1p(analysis[[v]])
}
roc_for <- function(outcome, score) {
  pROC::roc(response = factor(outcome, levels = c(0,1)), predictor = score,
            levels = c("0","1"), direction = ">", quiet = TRUE)
}
single <- rbindlist(lapply(metric_names, function(v) {
  r <- roc_for(analysis$death, auc_direction(v))
  ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
  data.table(metric = v, AUC = as.numeric(pROC::auc(r)),
             CI95_lower = ci[[1L]], CI95_upper = ci[[3L]],
             direction = ">", n = nrow(analysis), events = sum(analysis$death))
}))

pair_specs <- list(
  list(name = "SII vs NLR", with_platelet = "SII", control = "NLR"),
  list(name = "PIV vs SIRI", with_platelet = "PIV", control = "SIRI"),
  list(name = "Lymphocyte vs NLR", with_platelet = "lymphocyte", control = "NLR")
)
pair_result <- function(dat, spec) {
  outcome <- dat$death
  score_a <- if (spec$with_platelet == "lymphocyte") log1p(dat[[spec$with_platelet]]) else -log1p(dat[[spec$with_platelet]])
  score_b <- if (spec$control == "lymphocyte") log1p(dat[[spec$control]]) else -log1p(dat[[spec$control]])
  r_a <- roc_for(outcome, score_a)
  r_b <- roc_for(outcome, score_b)
  test <- pROC::roc.test(r_b, r_a, paired = TRUE, method = "delong")
  data.table(comparison = spec$name, AUC_a = as.numeric(pROC::auc(r_a)),
             AUC_control = as.numeric(pROC::auc(r_b)),
             delta_AUC = as.numeric(pROC::auc(r_a)) - as.numeric(pROC::auc(r_b)),
             direction = if (as.numeric(pROC::auc(r_a))-as.numeric(pROC::auc(r_b)) < 0) "negative" else "nonnegative",
             DeLong_p = unname(test$p.value), n = nrow(dat), events = sum(outcome))
}
primary_specs <- pair_specs[1:2]
paired <- rbindlist(lapply(primary_specs, function(s) pair_result(analysis, s)))
lymph_nlr <- pair_result(analysis, pair_specs[[3L]])

set.seed(20261115)
B <- 2000L
workers <- min(4L, max(1L, parallel::detectCores(logical = FALSE)-1L))
idx0 <- which(analysis$death == 0L)
idx1 <- which(analysis$death == 1L)
score_a <- lapply(primary_specs, function(s) log1p(analysis[[s$with_platelet]]))
score_b <- lapply(primary_specs, function(s) log1p(analysis[[s$control]]))
auc_rank <- function(y, score) {
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n1*(n1+1)/2)/(n1*n0)
}
boot_worker <- function(ids, y, score_a, score_b, idx0, idx1) {
  m <- matrix(NA_real_, nrow = length(ids), ncol = length(score_a))
  for (i in seq_along(ids)) {
    ix <- c(sample(idx0,length(idx0),TRUE),sample(idx1,length(idx1),TRUE))
    for (j in seq_along(score_a)) m[i,j] <- auc_rank(y[ix],score_a[[j]][ix])-auc_rank(y[ix],score_b[[j]][ix])
  }
  m
}
chunks <- split(seq_len(B), cut(seq_len(B), breaks = workers, labels = FALSE))
cl <- parallel::makeCluster(workers)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterSetRNGStream(cl, 20261115)
parallel::clusterExport(cl, c("analysis","score_a","score_b","idx0","idx1","auc_rank","boot_worker"),
                        envir = environment())
parts <- parallel::parLapply(cl, chunks, boot_worker,
                             y = analysis$death, score_a = score_a, score_b = score_b,
                             idx0 = idx0, idx1 = idx1)
boot <- do.call(rbind, parts)
boot_ci <- t(apply(boot, 2L, quantile, c(.025,.975), na.rm = TRUE, names = FALSE, type = 6))
paired[, `:=`(CI95_lower = boot_ci[,1L], CI95_upper = boot_ci[,2L],
              bootstrap_replicates = B)]

pre_analysis <- pre_index[alc > 0 & outcome_available == TRUE]
pre_single <- rbindlist(lapply(metric_names, function(v) {
  score <- if (v == "lymphocyte") log1p(pre_analysis[[v]]) else -log1p(pre_analysis[[v]])
  r <- roc_for(pre_analysis$death, score)
  data.table(metric = v, AUC = as.numeric(pROC::auc(r)), n = nrow(pre_analysis),
             events = sum(pre_analysis$death), window = "[-1440,0]")
}))
pre_pairs <- rbindlist(lapply(primary_specs, function(s) {
  x <- pair_result(pre_analysis, s)
  x[, window := "[-1440,0]"]
  x
}))

# Compare point estimates to completed MIMIC tables.
mimic_auc <- fread(file.path(P_TABLES, "table2_single_auc.csv"))
mimic_lymph <- fread(file.path(P_TABLES, "lymphocyte_auc_positive_control.csv"))
mimic_auc <- rbind(mimic_auc[, .(metric = fifelse(metric == "neutrophil (positive control)","neutrophil",metric), AUC)],
                   mimic_lymph[metric == "lymphocyte_absolute", .(metric="lymphocyte",AUC)])
mimic_pairs <- fread(file.path(P_TABLES, "table3_paired_delta_v2.csv"))
mimic_pairs <- mimic_pairs[, .(comparison, mimic_delta_AUC = delta_AUC)]
compare_auc <- merge(mimic_auc, single[, .(metric,AUC_eICU=AUC)], by="metric", all=TRUE)
compare_auc[, delta_eICU_minus_MIMIC := AUC_eICU-AUC]
compare_pairs <- merge(mimic_pairs, paired[, .(comparison,eICU_delta_AUC=delta_AUC)], by="comparison", all=TRUE)
compare_pairs[, delta_difference := eICU_delta_AUC-mimic_delta_AUC]

main_flow <- data.table(
  stage = c("SICU first adult stays", "strict five-item same-offset stays",
            "stays with valid numeric main-window panel", "unique people with a valid main stay",
            "patient-level index stays with outcome offset", "patient-level analysis after positive ALC"),
  n = c(9824L, 3024L, nrow(main_stay), uniqueN(main_stay$uniquepid),
        sum(main_patient$outcome_available), nrow(analysis)),
  events = c(NA_integer_, 287L,
             sum(tolower(trimws(main_stay$hospitaldischargelocation)) == "death"),
             sum(main_patient$death), sum(main_patient$death[main_patient$outcome_available]),
             sum(analysis$death))
)

fwrite(main_flow, file.path(P_TABLES, "eicu_sicu_patient_level_flow.csv"))
fwrite(main_patient, file.path(P_TABLES, "eicu_sicu_patient_index_roster.csv"), na = "")
fwrite(analysis, file.path(P_TABLES, "eicu_sicu_external_validation_dataset.csv"), na = "")
fwrite(single, file.path(P_TABLES, "eicu_sicu_single_auc.csv"))
fwrite(paired, file.path(P_TABLES, "eicu_sicu_paired_delong_bootstrap.csv"))
fwrite(lymph_nlr, file.path(P_TABLES, "eicu_sicu_lymphocyte_vs_nlr.csv"))
fwrite(pre_single, file.path(P_TABLES, "eicu_sicu_prewindow_single_auc.csv"))
fwrite(pre_pairs, file.path(P_TABLES, "eicu_sicu_prewindow_pairs.csv"))
fwrite(compare_auc, file.path(P_TABLES, "eicu_vs_mimic_single_auc.csv"))
fwrite(compare_pairs, file.path(P_TABLES, "eicu_vs_mimic_paired_delta.csv"))

cat("Patient-level analytic n/events =", nrow(analysis), "/", sum(analysis$death),
    "; zero-ALC excluded =", excluded_zero_lymph, "\n")
print(paired)
print(compare_pairs)
