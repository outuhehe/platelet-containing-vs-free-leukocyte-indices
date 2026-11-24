# =============================================================================
# 03_paired_analysis.R
# Single-marker AUCs, paired DeLong differences vs NLR, 2000-replicate bootstrap
# CIs, Spearman correlation matrix and SII Q1/Q4 RMST, on the 879-row index panel.
# Reads interim/cbc_export_v2.csv (from 02); writes result tables to tables/.
# Statistical logic is unchanged from the validated outputs/run_cbc_paired_analysis.R.
# =============================================================================
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
Sys.setenv(RA_ROOT = ra_root)
source(file.path(ra_root, "R", "00_common.R"))
load_analysis_pkgs("pROC", "survRM2")

cat("Reading the validated specimen/charttime-aligned export...\n")
d <- fread(file.path(P_INTERIM, "cbc_export_v2.csv"))
n_raw <- nrow(d)
events_raw <- sum(as.integer(d$hospital_expire_flag))
if (uniqueN(d$hadm_id) != n_raw || abs(n_raw - 888L) >= 30L || events_raw != 87L) {
  stop(sprintf("Hard stop: export is %d rows, %d unique admissions, %d events; expected 888 +/- <30 and 87",
               n_raw, uniqueN(d$hadm_id), events_raw))
}
if (any(!is.finite(as.matrix(d[, .(nlr,sii,piv,siri,plr,anc_single)])))) {
  stop("Non-finite value detected in the exported index set")
}

selected_time <- d$dischtime
selected_time[d$hospital_expire_flag == 1L] <- d$deathtime[d$hospital_expire_flag == 1L]
d[, time_days := as.numeric(difftime(selected_time, admittime, units = "days"))]
d[, status := as.integer(hospital_expire_flag)]
bad_selected_time <- is.na(d$time_days) | !is.finite(d$time_days) | d$time_days < 0
cat("Export n =", n_raw, "events =", events_raw,
    "; unusable selected survival time =", sum(bad_selected_time), "\n")
# Audit copy with derived survival columns; never overwrites the raw interim extract.
fwrite(d, file.path(P_INTERIM, "cbc_export_with_survival.csv"), na = "")

as_roc <- function(y, score) {
  pROC::roc(response = factor(y, levels = c(0,1)), predictor = score,
            levels = c("0","1"), direction = "<", quiet = TRUE)
}
metric_map <- c(NLR = "nlr", SII = "sii", PIV = "piv", SIRI = "siri",
                PLR = "plr", "neutrophil (positive control)" = "anc_single")
auc_rows <- lapply(names(metric_map), function(nm) {
  r <- as_roc(d$hospital_expire_flag, log(d[[metric_map[[nm]]]]))
  ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
  data.table(metric = nm, AUC = as.numeric(pROC::auc(r)),
             CI95_lower = ci[[1L]], CI95_upper = ci[[3L]],
             n = n_raw, events = events_raw)
})
table2 <- rbindlist(auc_rows)
anc_roc <- as_roc(d$hospital_expire_flag, log(d$anc_single))
anc_auc <- as.numeric(pROC::auc(anc_roc))
anc_se <- sqrt(as.numeric(pROC::var(anc_roc, method = "delong")))
anc_p <- 2 * pnorm(-abs((anc_auc - 0.5) / anc_se))

markers <- c(SII = "sii", PIV = "piv", PLR = "plr", SIRI = "siri")
pair_rows <- lapply(names(markers), function(nm) {
  r_marker <- as_roc(d$hospital_expire_flag, log(d[[markers[[nm]]]]))
  r_nlr <- as_roc(d$hospital_expire_flag, log(d$nlr))
  test <- pROC::roc.test(r_marker, r_nlr, paired = TRUE, method = "delong")
  auc_marker <- as.numeric(pROC::auc(r_marker))
  auc_nlr <- as.numeric(pROC::auc(r_nlr))
  data.table(comparison = paste(nm, "vs NLR"), AUC_marker = auc_marker,
             AUC_NLR = auc_nlr, delta_AUC = auc_marker - auc_nlr,
             direction = if (auc_marker - auc_nlr < 0) "negative" else "nonnegative",
             DeLong_p = unname(test$p.value), n = n_raw, events = events_raw)
})
table3 <- rbindlist(pair_rows)
table3[, p_BH := p.adjust(DeLong_p, method = "BH")]

set.seed(20261115)
B <- 2000L
workers <- min(4L, max(1L, parallel::detectCores(logical = FALSE) - 1L))
idx_by_outcome <- split(seq_len(nrow(d)), d$hospital_expire_flag)
boot_worker <- function(ids, dt, markers, idx_by_outcome) {
  out <- matrix(NA_real_, nrow = length(ids), ncol = length(markers))
  for (i in seq_along(ids)) {
    ix <- c(sample(idx_by_outcome[["0"]], length(idx_by_outcome[["0"]]), TRUE),
            sample(idx_by_outcome[["1"]], length(idx_by_outcome[["1"]]), TRUE))
    yy <- dt$hospital_expire_flag[ix]
    nlr_auc <- as.numeric(pROC::auc(as_roc(yy, log(dt$nlr[ix]))))
    for (j in seq_along(markers)) {
      marker_auc <- as.numeric(pROC::auc(as_roc(yy, log(dt[[markers[[j]]]][ix]))))
      out[i,j] <- marker_auc - nlr_auc
    }
  }
  out
}
chunks <- split(seq_len(B), cut(seq_len(B), breaks = workers, labels = FALSE))
cl <- parallel::makeCluster(workers)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterSetRNGStream(cl, 20261115)
parallel::clusterExport(cl, c("d", "markers", "idx_by_outcome", "as_roc",
                               "boot_worker"), envir = environment())
boot_parts <- parallel::parLapply(cl, chunks, boot_worker,
                                  dt = d, markers = markers,
                                  idx_by_outcome = idx_by_outcome)
boot <- do.call(rbind, boot_parts)
boot_ci <- t(apply(boot, 2L, quantile, probs = c(.025,.975),
                   na.rm = TRUE, names = FALSE, type = 6))
table_boot <- data.table(comparison = table3$comparison,
                         delta_AUC = table3$delta_AUC,
                         CI95_lower = boot_ci[,1L], CI95_upper = boot_ci[,2L],
                         bootstrap_replicates = B)
table3 <- merge(table3, table_boot[, .(comparison, CI95_lower, CI95_upper)],
                by = "comparison", sort = FALSE)

cor_vars <- d[, .(sii, nlr, piv, siri, plr)]
cor_mat <- cor(cor_vars, method = "spearman", use = "complete.obs")
cor_table <- data.table(variable = rownames(cor_mat))
for (nm in colnames(cor_mat)) cor_table[, (nm) := cor_mat[,nm]]

d[, quartile := pmin(4L, ceiling(4 * frank(sii, ties.method = "first") / .N))]
surv_dat <- d[!bad_selected_time]
surv_n <- nrow(surv_dat)
if (surv_n == 0L || max(surv_dat$time_days) < 30) {
  stop("RMST cannot be estimated at tau=30 with the available follow-up")
}
rmst_dat <- surv_dat[quartile %in% c(1L,4L)]
rmst_dat[, arm := as.integer(quartile == 4L)]
if (uniqueN(rmst_dat$arm) != 2L) stop("SII Q1/Q4 RMST groups are not both present")
tau_median <- median(surv_dat$time_days)
rmst_at <- function(tau, label) {
  fit <- survRM2::rmst2(time = rmst_dat$time_days, status = rmst_dat$status,
                        arm = rmst_dat$arm, tau = tau)
  res <- fit$unadjusted.result
  diff_row <- grep("RMST.*arm=1.*arm=0", rownames(res), ignore.case = TRUE)
  if (!length(diff_row)) diff_row <- 3L
  vals <- res[diff_row[[1L]], ]
  data.table(analysis = label, tau_days = tau,
             rmst_Q1 = unname(fit$RMST.arm0$rmst[["Est."]]),
             rmst_Q4 = unname(fit$RMST.arm1$rmst[["Est."]]),
             difference_Q4_minus_Q1 = unname(vals["Est."]),
             CI95_lower = unname(vals["lower .95"]),
             CI95_upper = unname(vals["upper .95"]),
             p_value = unname(vals["p"]),
             n_Q1 = sum(rmst_dat$arm == 0L), n_Q4 = sum(rmst_dat$arm == 1L),
             events_Q1 = sum(rmst_dat$status[rmst_dat$arm == 0L]),
             events_Q4 = sum(rmst_dat$status[rmst_dat$arm == 1L]))
}
table_rmst <- rbind(rmst_at(30, "tau=30 days"),
                    rmst_at(tau_median, "tau=median observed stay"))

sub <- d[all_abs %in% c(TRUE, 1L)]
sub_delta <- as.numeric(pROC::auc(as_roc(sub$hospital_expire_flag, log(sub$sii)))) -
  as.numeric(pROC::auc(as_roc(sub$hospital_expire_flag, log(sub$nlr))))

fwrite(table2, file.path(P_TABLES, "table2_single_auc.csv"))
fwrite(table3, file.path(P_TABLES, "table3_paired_delta.csv"))
fwrite(table_boot, file.path(P_TABLES, "table3b_bootstrap_ci.csv"))
fwrite(cor_table, file.path(P_TABLES, "tableS_correlation.csv"))
fwrite(table_rmst, file.path(P_TABLES, "table4_rmst.csv"))

cat("n =", n_raw, "events =", events_raw,
    "EPV(7 parameters) =", round(events_raw/7,1),
    "; sample-size difference from 888 =", n_raw-888L, "\n")
cat("Spearman SII-NLR r =", round(cor_mat["sii","nlr"],3),
    "; PIV-SIRI r =", round(cor_mat["piv","siri"],3), "\n")
cat("ANC positive-control AUC =", round(anc_auc,3),
    "p(vs 0.5) =", format.pval(anc_p, digits = 3), "\n")
cat("Survival: n =", surv_n, "events =", sum(surv_dat$status),
    "; >30d =", round(mean(surv_dat$time_days > 30)*100,1),
    "% ; median =", round(tau_median,1), "days\n")
cat("RMST Q1/Q4 group counts and estimates:\n"); print(table_rmst)
cat("All-absolute subgroup: n =", nrow(sub),
    "events =", sum(sub$hospital_expire_flag),
    "; directional delta AUC (SII-NLR) =", round(sub_delta,3), "\n")
cat("Pairs with negative delta =", sum(table3$delta_AUC < 0),
    "; BH-FDR < 0.05 =", sum(table3$p_BH < 0.05), "\n")
cat("Tables 2, 3, 3b, S and 4 written to tables/.\n")
