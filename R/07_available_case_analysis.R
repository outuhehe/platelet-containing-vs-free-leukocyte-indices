# =============================================================================
# 07_available_case_analysis.R
# Per-comparison available-case primary analyses, common-sample sensitivity,
# subject-clustered and first-stay bootstrap, admission-type stratification,
# platelet tertile/rcs, missingness mechanism and estimation accounting.
# Reads interim panels + roster (from 01); writes all results to tables/.
# No imputation; all bootstrap CIs are percentile intervals. Logic unchanged.
# =============================================================================
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
Sys.setenv(RA_ROOT = ra_root)
source(file.path(ra_root, "R", "00_common.R"))
load_analysis_pkgs("pROC")

set.seed(20261117)
B <- 2000L
dir_out <- P_TABLES

input_panels <- file.path(P_INTERIM, "cbc_available_case_panels.csv")
input_roster <- file.path(P_INTERIM, "cohort_roster_2376.csv")
if (!file.exists(input_panels) || !file.exists(input_roster)) {
  stop("Missing SQL exports. Required: ", input_panels, " and ", input_roster)
}
panels <- read.csv(input_panels,
                   stringsAsFactors = FALSE, check.names = TRUE)
roster <- read.csv(input_roster,
                   stringsAsFactors = FALSE, check.names = TRUE)

# 0. Reconstruct the five indices in R (Query A exports cell counts only).
mk_indices <- function(d) {
  d$nlr  <- d$anc_f / d$alc_f
  d$sii  <- d$plt * d$anc_f / d$alc_f
  d$siri <- d$anc_f * d$amc_f / d$alc_f
  d$piv  <- d$plt * d$anc_f * d$amc_f / d$alc_f
  d$plr  <- d$plt / d$alc_f
  d$inv_lymph <- -log(d$alc_f)
  d
}
panels <- mk_indices(panels)

# Query A exports all candidate panels; select one independently for each
# comparison so a missing monocyte cannot discard an otherwise usable panel.
best_panel <- function(d, require_monocyte = FALSE) {
  ok <- is.finite(d$anc_f) & d$anc_f > 0 &
         is.finite(d$alc_f) & d$alc_f > 0 &
         is.finite(d$plt) & d$plt > 0
  if (require_monocyte) ok <- ok & is.finite(d$amc_f) & d$amc_f > 0
  d <- d[ok, , drop = FALSE]
  if (!nrow(d)) return(d)
  d$all_abs <- as.logical(d$all_abs)
  panel_time_rank <- rank(as.character(d$panel_time), ties.method = "min")
  d <- d[order(d$hadm_id, -as.integer(d$all_abs), -panel_time_rank,
               is.na(d$specimen_id), d$specimen_id), , drop = FALSE]
  d[!duplicated(d$hadm_id), , drop = FALSE]
}

# Paired engine: AUCs, delta (containing - free), two-sided DeLong p, bootstrap CI.
paired_one <- function(d, free, containing, label) {
  y <- as.integer(d$hospital_expire_flag)
  rf <- roc(y, d[[free]], quiet = TRUE, direction = "<")
  rc <- roc(y, d[[containing]], quiet = TRUE, direction = "<")
  tt <- roc.test(rf, rc, paired = TRUE, method = "delong")
  i0 <- which(y == 0L); i1 <- which(y == 1L)
  deltas <- replicate(B, {
    ix <- c(sample(i0, length(i0), replace = TRUE),
            sample(i1, length(i1), replace = TRUE))
    auc(roc(y[ix], d[[containing]][ix], quiet = TRUE, direction = "<")) -
      auc(roc(y[ix], d[[free]][ix], quiet = TRUE, direction = "<"))
  })
  qs <- quantile(deltas, c(0.025, 0.975), na.rm = TRUE)
  data.frame(
    comparison = label,
    AUC_free = as.numeric(auc(rf)),
    AUC_containing = as.numeric(auc(rc)),
    delta_AUC = as.numeric(auc(rc) - auc(rf)),
    CI95_lower = as.numeric(qs[[1]]),
    CI95_upper = as.numeric(qs[[2]]),
    DeLong_p = tt$p.value,
    n = nrow(d),
    events = sum(y),
    stringsAsFactors = FALSE
  )
}

# 1. Available-case primary comparisons (per-comparison maximum samples).
s_nlr_sii <- best_panel(panels, require_monocyte = FALSE)
av_nlr_sii <- paired_one(s_nlr_sii, "nlr", "sii", "NLR vs SII (available-case)")
s_siri_piv <- best_panel(panels, require_monocyte = TRUE)
av_siri_piv <- paired_one(s_siri_piv, "siri", "piv", "SIRI vs PIV (available-case)")
av_primary <- rbind(av_nlr_sii, av_siri_piv)
av_primary$BH_FDR <- p.adjust(av_primary$DeLong_p, method = "BH")
write.csv(av_primary, file.path(dir_out, "available_case_primary.csv"),
          row.names = FALSE)
print(av_primary, digits = 4, row.names = FALSE)

av_plr <- paired_one(s_nlr_sii, "inv_lymph", "plr",
                     "Inverse lymphocyte vs PLR (available-case, descriptive)")
write.csv(av_plr, file.path(dir_out, "available_case_plr.csv"), row.names = FALSE)

# 2. Common-sample sensitivity (reproduces the frozen 879-stay manuscript analysis).
common <- best_panel(panels, require_monocyte = TRUE)
cat("Common five-index sample: n =", nrow(common),
    "events =", sum(common$hospital_expire_flag), "\n")
stopifnot(nrow(common) == 879L, sum(common$hospital_expire_flag) == 87L)
cm_nlr_sii <- paired_one(common, "nlr", "sii", "NLR vs SII (common sample)")
cm_siri_piv <- paired_one(common, "siri", "piv", "SIRI vs PIV (common sample)")
common_tbl <- rbind(cm_nlr_sii, cm_siri_piv)
common_tbl$BH_FDR <- p.adjust(common_tbl$DeLong_p, method = "BH")
write.csv(common_tbl, file.path(dir_out, "common_sample_sensitivity.csv"),
          row.names = FALSE)

# 3. Subject-clustered bootstrap: resample whole subjects, keeping all stays.
cluster_boot <- function(d, free, containing) {
  y <- as.integer(d$hospital_expire_flag)
  subjects <- unique(d$subject_id)
  replicate(B, {
    ss <- sample(subjects, length(subjects), replace = TRUE)
    ix <- unlist(lapply(ss, function(s) which(d$subject_id == s)),
                 use.names = FALSE)
    yy <- y[ix]
    if (length(unique(yy)) < 2L) {
      NA_real_
    } else {
      auc(roc(yy, d[[containing]][ix], quiet = TRUE, direction = "<")) -
        auc(roc(yy, d[[free]][ix], quiet = TRUE, direction = "<"))
    }
  })
}
clu_nlr <- cluster_boot(common, "nlr", "sii")
clu_siri <- cluster_boot(common, "siri", "piv")
cluster_tbl <- data.frame(
  comparison = c("NLR vs SII", "SIRI vs PIV"),
  delta_definition = rep("platelet-containing minus platelet-free", 2),
  cluster_CI95_lower = c(quantile(clu_nlr, 0.025, na.rm = TRUE),
                         quantile(clu_siri, 0.025, na.rm = TRUE)),
  cluster_CI95_upper = c(quantile(clu_nlr, 0.975, na.rm = TRUE),
                         quantile(clu_siri, 0.975, na.rm = TRUE)),
  n_subjects = length(unique(common$subject_id))
)
write.csv(cluster_tbl, file.path(dir_out, "cluster_bootstrap_ci.csv"),
          row.names = FALSE)
print(cluster_tbl, digits = 4, row.names = FALSE)

# 4. First-stay-per-subject sensitivity.
common <- common[order(common$admittime), ]
first_idx <- !duplicated(common$subject_id)
first <- common[first_idx, ]
cat("First-stay sample: n =", nrow(first),
    "events =", sum(first$hospital_expire_flag), "\n")
fs_nlr <- paired_one(first, "nlr", "sii", "NLR vs SII (first stay)")
fs_siri <- paired_one(first, "siri", "piv", "SIRI vs PIV (first stay)")
first_tbl <- rbind(fs_nlr, fs_siri)
write.csv(first_tbl, file.path(dir_out, "first_stay_sensitivity.csv"),
          row.names = FALSE)

# 5. Emergency vs elective stratified primary comparisons.
is_emerg <- grepl("EMERGENCY|URGENT|DIRECT EMER|EW EMER", common$admission_type,
                  ignore.case = TRUE)
is_elective <- toupper(trimws(common$admission_type)) == "ELECTIVE"
stratum_counts <- data.frame(
  admission_type_group = c("emergency_or_urgent", "elective", "other_or_unclassified"),
  n = c(sum(is_emerg), sum(is_elective), sum(!is_emerg & !is_elective)),
  events = c(sum(common$hospital_expire_flag[is_emerg]),
             sum(common$hospital_expire_flag[is_elective]),
             sum(common$hospital_expire_flag[!is_emerg & !is_elective]))
)
write.csv(stratum_counts, file.path(dir_out, "admission_type_stratum_counts.csv"),
          row.names = FALSE)
print(stratum_counts, row.names = FALSE)

strata <- list(emergency = common[is_emerg, , drop = FALSE],
               elective = common[is_elective, , drop = FALSE])
strat_rows <- do.call(rbind, lapply(names(strata), function(label) {
    dd <- strata[[label]]
    if (nrow(dd) < 20 || sum(dd$hospital_expire_flag) < 5) return(NULL)
    a <- paired_one(dd, "nlr", "sii", paste("NLR vs SII", label))
    b <- paired_one(dd, "siri", "piv", paste("SIRI vs PIV", label))
    rbind(a, b)
  }))
write.csv(strat_rows, file.path(dir_out, "emergency_stratified.csv"),
          row.names = FALSE)
print(strat_rows, digits = 4, row.names = FALSE)

# 6. Platelet non-linearity: mortality by tertile plus an rcs logistic fit.
tert <- cut(common$plt,
            breaks = quantile(common$plt, c(0, 1/3, 2/3, 1), na.rm = TRUE),
            include.lowest = TRUE, labels = c("T1(low)", "T2", "T3(high)"))
plt_tertile <- do.call(rbind, lapply(levels(tert), function(lv) {
  dd <- common[tert == lv, ]
  data.frame(platelet_tertile = lv, n = nrow(dd),
             events = sum(dd$hospital_expire_flag),
             mortality_pct = 100 * mean(dd$hospital_expire_flag),
             stringsAsFactors = FALSE)
}))
write.csv(plt_tertile, file.path(dir_out, "platelet_tertile_mortality.csv"),
          row.names = FALSE)
print(plt_tertile, digits = 4, row.names = FALSE)

if (requireNamespace("rms", quietly = TRUE)) {
  dd <- common
  dd$y <- as.integer(dd$hospital_expire_flag)
  fit <- rms::lrm(y ~ rms::rcs(log(plt), 3), data = dd, x = TRUE, y = TRUE)
  sink(file.path(dir_out, "platelet_rcs_logistic.txt"))
  print(fit)
  cat("\nWald test of nonlinear spline term:\n")
  print(anova(fit))
  sink()
}

# 7. Complete-case vs non-complete baseline + missingness mechanism (n=2,376).
as_flag <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "1", "t", "yes")
}
roster$complete <- as_flag(roster$complete_five_index)
roster$has_any_panel <- as_flag(roster$has_any_panel)
stopifnot(nrow(roster) == 2376L,
          sum(roster$hospital_expire_flag) == 181L,
          sum(roster$complete) == 879L,
          sum(roster$hospital_expire_flag[roster$complete]) == 87L)
mech <- data.frame(
  category = c("full eligible cohort",
               "no usable panel in window",
               "panel present but monocyte unavailable",
               "five-index complete case"),
  n = c(nrow(roster),
        sum(!roster$has_any_panel),
        sum(roster$has_any_panel & !roster$complete),
        sum(roster$complete)),
  events = c(sum(roster$hospital_expire_flag),
             sum(roster$hospital_expire_flag[!roster$has_any_panel]),
             sum(roster$hospital_expire_flag[roster$has_any_panel & !roster$complete]),
             sum(roster$hospital_expire_flag[roster$complete]))
)
mech$mortality_pct <- round(100 * mech$events / mech$n, 2)
write.csv(mech, file.path(dir_out, "missingness_mechanism.csv"),
          row.names = FALSE)
print(mech, row.names = FALSE)

fmt_med <- function(x) {
  x <- x[is.finite(x)]
  paste0(round(median(x), 2), " [", round(quantile(x, .25), 2), ", ",
         round(quantile(x, .75), 2), "]")
}
grp <- roster$complete
base_rows <- list()
for (v in c("age", "charlson")) {
  base_rows[[length(base_rows) + 1L]] <- data.frame(
    variable = v,
    complete = fmt_med(roster[[v]][grp]),
    non_complete = fmt_med(roster[[v]][!grp]),
    p = tryCatch(format(wilcox.test(roster[[v]] ~ grp)$p.value,
                        digits = 3), error = function(e) "NA"),
    test = "Wilcoxon")
}
for (v in c("sex", "admission_type")) {
  tab <- table(roster[[v]], grp)
  p <- tryCatch(format(fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value,
                       digits = 3), error = function(e) "NA")
  base_rows[[length(base_rows) + 1L]] <- data.frame(
    variable = v,
    complete = paste0(rownames(tab), ": ", tab[, "TRUE"], " (",
                      round(100 * tab[, "TRUE"] / sum(tab[, "TRUE"]), 1), "%)",
                      collapse = "; "),
    non_complete = paste0(rownames(tab), ": ", tab[, "FALSE"], " (",
                          round(100 * tab[, "FALSE"] / sum(tab[, "FALSE"]), 1), "%)",
                          collapse = "; "),
    p = p, test = "Fisher (Monte Carlo)")
}
base_rows[[length(base_rows) + 1L]] <- data.frame(
  variable = "in_hospital_mortality_pct",
  complete = round(100 * mean(roster$hospital_expire_flag[grp]), 2),
  non_complete = round(100 * mean(roster$hospital_expire_flag[!grp]), 2),
  p = format(prop.test(c(sum(roster$hospital_expire_flag[grp]),
                        sum(roster$hospital_expire_flag[!grp])),
                      c(sum(grp), sum(!grp)))$p.value, digits = 3),
  test = "two-sample proportion")
baseline_tbl <- do.call(rbind, base_rows)
write.csv(baseline_tbl, file.path(dir_out, "complete_vs_noncomplete_baseline.csv"),
          row.names = FALSE)
print(baseline_tbl, row.names = FALSE)

# 8. Percentage-estimation accounting.
est_tbl <- data.frame(
  all_absolute = sum(as_flag(common$all_abs)),
  at_least_one_estimated = sum(!as_flag(common$all_abs)),
  total = nrow(common)
)
write.csv(est_tbl, file.path(dir_out, "estimation_accounting.csv"),
          row.names = FALSE)
print(est_tbl, row.names = FALSE)

cat("\nAll available-case analyses complete. Outputs written to", dir_out, "\n")
