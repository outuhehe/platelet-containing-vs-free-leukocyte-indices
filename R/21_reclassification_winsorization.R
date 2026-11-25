# Post-hoc exploratory NRI/IDI and winsorization sensitivity analyses.
# Reads restricted, locally regenerated intermediates; writes aggregate results only.
# Run after primary stages 01-15. No database connection or row-level export.

auc_rank <- function(y, score) {
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(score, ties.method = "average")[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

reclassification_metrics <- function(y, old, new) {
  stopifnot(length(y) == length(old), length(y) == length(new),
            all(is.finite(old)), all(is.finite(new)), all(old >= 0 & old <= 1),
            all(new >= 0 & new <= 1), all(y %in% 0:1))
  diff <- new - old
  movement <- ifelse(abs(diff) <= 1e-12, 0, sign(diff))
  event <- mean(movement[y == 1L]); nonevent <- -mean(movement[y == 0L])
  c(NRI_event = event, NRI_nonevent = nonevent, NRI_continuous = event + nonevent,
    IDI = mean(diff[y == 1L]) - mean(diff[y == 0L]))
}

make_folds <- function(y, subject, K = 10L) {
  ids <- unique(as.character(subject))
  event <- as.numeric(tapply(y, factor(as.character(subject), levels = ids), max))
  by_id <- setNames(integer(length(ids)), ids)
  for (s in 0:1) {
    eligible <- ids[event == s]
    eligible <- eligible[sample.int(length(eligible))]
    by_id[eligible] <- rep(seq_len(K), length.out = length(eligible))
  }
  unname(by_id[as.character(subject)])
}

cluster_resample <- function(d) {
  ids <- unique(as.character(d$subject_id))
  event <- as.numeric(tapply(d$death, factor(as.character(d$subject_id), levels = ids), max))
  sampled <- unlist(lapply(0:1, function(s) {
    z <- ids[event == s]; z[sample.int(length(z), length(z), replace = TRUE)]
  }), use.names = FALSE)
  rows <- split(seq_len(nrow(d)), as.character(d$subject_id))
  unlist(rows[sampled], use.names = FALSE)
}

cv_reclassification <- function(d, markers, repeats = 20L, K = 10L) {
  # Preserve script 16's prediction-model coding: EMER-labelled, ELECTIVE, other.
  # URGENT is included in other in this historical incremental-model specification.
  # This differs from M3 association models, which combine emergency and urgent.
  X0 <- model.matrix(~ age + sex + adm_prediction + charlson, data = d)
  X <- list(base = X0)
  for (m in markers) X[[paste0("add_", m)]] <- cbind(X0, log_index = log(d[[m]]))
  combos <- list(c("base", paste0("add_", markers[1])),
                c("base", paste0("add_", markers[2])),
                c(paste0("add_", markers[2]), paste0("add_", markers[1])))
  values <- array(NA_real_, c(length(combos), 4L, repeats))
  unconverged <- 0L
  for (r in seq_len(repeats)) {
    folds <- make_folds(d$death, d$subject_id, K)
    pred <- matrix(NA_real_, nrow(d), length(X), dimnames = list(NULL, names(X)))
    for (k in seq_len(K)) {
      train <- folds != k; test <- !train
      for (nm in names(X)) {
        fit <- suppressWarnings(glm.fit(X[[nm]][train, , drop = FALSE], d$death[train],
                                       family = binomial(), control = glm.control(maxit = 50L)))
        if (!isTRUE(fit$converged)) unconverged <- unconverged + 1L
        beta <- fit$coefficients; beta[is.na(beta)] <- 0
        pred[test, nm] <- plogis(drop(X[[nm]][test, , drop = FALSE] %*% beta))
      }
    }
    stopifnot(all(is.finite(pred)))
    for (j in seq_along(combos)) values[j, , r] <- reclassification_metrics(d$death, pred[, combos[[j]][1]], pred[, combos[[j]][2]])
  }
  point <- apply(values, c(1, 2), mean)
  split_sd <- apply(values, c(1, 2), sd)
  dimnames(point) <- dimnames(split_sd) <- list(vapply(combos, paste, collapse = "_to_", FUN.VALUE = ""),
    c("NRI_event", "NRI_nonevent", "NRI_continuous", "IDI"))
  list(point = point, split_sd = split_sd, nonconverged_fits = unconverged)
}

winsorize <- function(x, tail) {
  if (tail == 0) return(x)
  q <- quantile(x, c(tail, 1 - tail), type = 7, names = FALSE)
  pmin(pmax(x, q[1]), q[2])
}

apply_winsor <- function(d, scheme, tail) {
  z <- d
  if (scheme == "components") {
    for (m in c("anc", "alc", "amc", "plt")) z[[m]] <- winsorize(z[[m]], tail)
    z$nlr <- z$anc / z$alc; z$sii <- z$plt * z$nlr
    z$siri <- z$anc * z$amc / z$alc; z$piv <- z$plt * z$siri
    z$plr <- z$plt / z$alc
  } else if (scheme == "indices") {
    for (m in c("nlr", "sii", "siri", "piv", "plr")) z[[m]] <- winsorize(z[[m]], tail)
  } else if (scheme != "raw") stop("Unknown winsorization scheme")
  z
}

main <- function() {
  RA <- normalizePath(Sys.getenv("RA_ROOT", unset = "reproducible_analysis"), winslash = "/", mustWork = TRUE)
  OUT <- Sys.getenv("STUDY_NRI_OUT", unset = file.path(RA, "outputs", "reclassification_winsorization"))
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  B <- as.integer(Sys.getenv("STUDY_NRI_B", "500"))
  BW <- as.integer(Sys.getenv("STUDY_WINSOR_B", "2000"))
  REPS <- as.integer(Sys.getenv("STUDY_NRI_REPEATS", "20"))
  WORKERS <- max(1L, min(as.integer(Sys.getenv("STUDY_NRI_WORKERS", "4")), parallel::detectCores()))
  stopifnot(B >= 2L, BW >= 2L, REPS >= 2L)
  set.seed(20261121)
  roster <- read.csv(file.path(RA, "outputs", "interim", "cohort_roster_2376.csv"))
  aci <- read.csv(file.path(RA, "outputs", "interim", "cbc_available_case_index.csv"))
  aci <- merge(aci, roster[, c("hadm_id", "subject_id", "admission_type", "age", "sex", "charlson")], by = "hadm_id", sort = TRUE)
  aci$death <- as.integer(aci$hospital_expire_flag)
  aci$sex <- factor(aci$sex)
  aci$adm_prediction <- factor(ifelse(grepl("EMER", aci$admission_type, ignore.case = TRUE), "emergency",
    ifelse(grepl("ELECTIVE", aci$admission_type, ignore.case = TRUE), "elective", "other")))
  samples <- list(NLR_SII = aci[is.finite(aci$nlr) & is.finite(aci$sii), ],
                  SIRI_PIV = aci[is.finite(aci$siri) & is.finite(aci$piv), ])
  pairs <- list(NLR_SII = c("nlr", "sii"), SIRI_PIV = c("siri", "piv"))
  stopifnot(nrow(samples[[1]]) == 886L, nrow(samples[[2]]) == 879L,
            sum(samples[[1]]$death) == 87L, sum(samples[[2]]$death) == 87L)
  for (d in samples) stopifnot(!anyNA(d[, c("subject_id", "age", "sex", "adm_prediction", "charlson")]))
  point <- lapply(names(samples), function(nm) {
    set.seed(20261121L)
    cv_reclassification(samples[[nm]], pairs[[nm]], REPS)
  }); names(point) <- names(samples)
  message("Point estimates complete; full refit bootstrap B=", B, ", repeated CV=", REPS)

  boot_one <- function(b) {
    out <- list()
    for (nm in names(samples)) {
      set.seed(20261121L + 1000L + b)
      d <- samples[[nm]]; ix <- cluster_resample(d)
      out[[nm]] <- cv_reclassification(d[ix, , drop = FALSE], pairs[[nm]], REPS)
    }
    out
  }
  cl <- parallel::makePSOCKcluster(WORKERS)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, c("make_folds", "cluster_resample", "cv_reclassification", "reclassification_metrics"), envir = environment(main))
  parallel::clusterExport(cl, c("samples", "pairs", "REPS", "boot_one"), envir = environment())
  boot <- list()
  chunk <- max(4L, WORKERS * 2L)
  for (start in seq(1L, B, by = chunk)) {
    ids <- start:min(B, start + chunk - 1L)
    boot[ids] <- parallel::parLapply(cl, ids, boot_one)
    message("Bootstrap completed ", max(ids), "/", B)
  }
  rows <- list()
  for (nm in names(samples)) {
    for (j in seq_len(3L)) {
      row <- data.frame(sample = nm, comparison_old_to_new = rownames(point[[nm]]$point)[j],
        n = nrow(samples[[nm]]), events = sum(samples[[nm]]$death), unique_subjects = length(unique(samples[[nm]]$subject_id)))
      for (m in colnames(point[[nm]]$point)) {
        v <- vapply(boot, function(z) z[[nm]]$point[j, m], numeric(1))
        ci <- quantile(v, c(.025, .975), names = FALSE)
        row[[m]] <- point[[nm]]$point[j, m]; row[[paste0(m, "_CI_lower")]] <- ci[1]
        row[[paste0(m, "_CI_upper")]] <- ci[2]; row[[paste0(m, "_partition_SD")]] <- point[[nm]]$split_sd[j, m]
      }
      row$bootstrap_replicates <- B; row$cv_repeats <- REPS; row$cv_folds <- 10L
      rows[[length(rows) + 1L]] <- row
    }
  }
  write.csv(do.call(rbind, rows), file.path(OUT, "reclassification_NRI_IDI.csv"), row.names = FALSE)
  message("NRI/IDI aggregate table written")

  d <- read.csv(file.path(RA, "outputs", "interim", "cbc_export_v2.csv"), check.names = TRUE)
  d <- merge(d, roster[, c("hadm_id", "admission_type")], by = "hadm_id", sort = TRUE)
  d$death <- as.integer(d$hospital_expire_flag); d$sex <- factor(d$sex)
  d$adm_assoc <- factor(ifelse(grepl("EMER|URGENT", d$admission_type, ignore.case = TRUE), "emergency_urgent",
    ifelse(grepl("ELECTIVE", d$admission_type, ignore.case = TRUE), "elective", "other")))
  stopifnot(nrow(d) == 879L, sum(d$death) == 87L,
            all(abs(log(d$nlr / (d$anc / d$alc))) < 1e-8),
            all(abs(log(d$siri / (d$anc * d$amc / d$alc))) < 1e-8),
            all(abs(log(d$sii / (d$plt * d$nlr))) < 1e-8),
            all(abs(log(d$piv / (d$plt * d$siri))) < 1e-8))
  schemes <- data.frame(scheme = c("raw", "components", "components", "indices", "indices"),
                        tail = c(0, .01, .025, .01, .025))
  or_rows <- auc_rows <- list()
  for (s in seq_len(nrow(schemes))) {
    scheme <- schemes$scheme[s]; tail <- schemes$tail[s]
    z <- apply_winsor(d, scheme, tail)
    for (m in c("nlr", "sii", "siri", "piv", "plr")) {
      z$z <- as.numeric(scale(log(z[[m]])))
      fit <- glm(death ~ z + age + sex + adm_assoc + charlson, family = binomial(), data = z)
      est <- coef(summary(fit))["z", ]
      or_rows[[length(or_rows) + 1L]] <- data.frame(scheme = scheme, lower_quantile = tail, upper_quantile = 1 - tail,
        marker = toupper(m), n = nrow(z), events = sum(z$death), OR_per_SD_log_index = exp(est[1]),
        CI_lower = exp(est[1] - 1.96 * est[2]), CI_upper = exp(est[1] + 1.96 * est[2]), p = est[4])
    }
    set.seed(20261121L + 50000L)
    boot_auc <- replicate(BW, {
      ix <- cluster_resample(d); zz <- apply_winsor(d[ix, , drop = FALSE], scheme, tail)
      c(auc_rank(zz$death, zz$nlr) - auc_rank(zz$death, zz$sii),
        auc_rank(zz$death, zz$siri) - auc_rank(zz$death, zz$piv))
    })
    for (j in seq_along(pairs)) {
      m <- pairs[[j]]; ci <- quantile(boot_auc[j, ], c(.025, .975), names = FALSE)
      af <- auc_rank(z$death, z[[m[1]]]); ac <- auc_rank(z$death, z[[m[2]]])
      auc_rows[[length(auc_rows) + 1L]] <- data.frame(scheme = scheme, lower_quantile = tail, upper_quantile = 1 - tail,
        comparison = paste(toupper(m), collapse = "_vs_"), n = nrow(z), events = sum(z$death),
        AUC_platelet_free = af, AUC_platelet_containing = ac, delta_AUC_free_minus_containing = af - ac,
        CI_lower = ci[1], CI_upper = ci[2], bootstrap_replicates = BW)
    }
    message("Winsorization complete: ", scheme, " tail=", tail)
  }
  write.csv(do.call(rbind, auc_rows), file.path(OUT, "winsorized_paired_AUC.csv"), row.names = FALSE)
  write.csv(do.call(rbind, or_rows), file.path(OUT, "winsorized_adjusted_OR.csv"), row.names = FALSE)
  nonconv <- sum(vapply(point, function(z) z$nonconverged_fits, integer(1))) +
    sum(vapply(boot, function(z) sum(vapply(z, function(q) q$nonconverged_fits, integer(1))), integer(1)))
  writeLines(c("Exploratory post-hoc NRI/IDI and winsorization analysis.",
    paste("Full refit, subject-clustered, outcome-stratified bootstrap replicates:", B),
    paste("Cross-validation:", REPS, "repetitions of subject-grouped stratified 10-fold CV per bootstrap sample."),
    "Metrics are means across partitions; split SD is not a confidence interval.",
    "NRI is continuous (category-free); no clinical risk categories were specified.",
    "NRI/IDI compare old-to-new model risks; positive values favour the new model.",
    "IDI is on the absolute probability scale. CIs are bootstrap percentile intervals.",
    "No confirmatory p-values or clinical-utility claims are made for these post-hoc metrics.",
    "Prediction-model admission coding matches script 16 (EMER-labelled, ELECTIVE, other; URGENT in other).",
    "M3 association-model coding matches script 19 (EMER or URGENT combined, ELECTIVE, other).",
    "Winsorization uses the common five-index sample; no rows are removed.",
    "Components are capped before formulas are rebuilt; index-wise caps are a separate sensitivity.",
    "Winsorization quantiles use R type 7 and are refitted inside each paired-AUC bootstrap resample.",
    paste("Paired-AUC subject-clustered bootstrap replicates:", BW),
    paste("Nonconverged logistic fits:", nonconv),
    "Synthetic unit checks are provided separately; no restricted data are distributed."), file.path(OUT, "NRI_IDI_winsorization_methods.txt"))
  stopifnot(nonconv == 0L)
  message("Analysis complete; all logistic fits converged")
}

if (Sys.getenv("STUDY_SKIP_MAIN", "0") != "1") main()
