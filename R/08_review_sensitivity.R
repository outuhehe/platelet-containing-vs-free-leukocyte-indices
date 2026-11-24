# =============================================================================
# 08_review_sensitivity.R
# Sensitivity of the paired results to panel-selection rules (absolute-first vs
# latest eligible panel), separately for NLR-vs-SII (no monocyte) and
# SIRI-vs-PIV (monocyte required), with 2000-replicate bootstrap CIs.
# Reads interim/cbc_available_case_panels.csv; writes tables/. Logic unchanged.
# =============================================================================
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
Sys.setenv(RA_ROOT = ra_root)
source(file.path(ra_root, "R", "00_common.R"))
load_analysis_pkgs("pROC")

out_dir <- P_TABLES
panels <- read.csv(file.path(P_INTERIM, "cbc_available_case_panels.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE)

as_flag <- function(x) tolower(as.character(x)) %in% c("true", "1", "t", "yes")
panels$all_abs <- as_flag(panels$all_abs)
panels$panel_time <- as.character(panels$panel_time)
panels$nlr <- panels$anc_f / panels$alc_f
panels$sii <- panels$plt * panels$anc_f / panels$alc_f
panels$siri <- panels$anc_f * panels$amc_f / panels$alc_f
panels$piv <- panels$plt * panels$anc_f * panels$amc_f / panels$alc_f

pick_panel <- function(d, require_monocyte, prefer_absolute) {
  ok <- is.finite(d$anc_f) & d$anc_f > 0 & is.finite(d$alc_f) & d$alc_f > 0 &
    is.finite(d$plt) & d$plt > 0
  if (require_monocyte) ok <- ok & is.finite(d$amc_f) & d$amc_f > 0
  d <- d[ok, , drop = FALSE]
  if (!nrow(d)) return(d)
  if (prefer_absolute) {
    d <- d[order(d$hadm_id, -as.integer(d$all_abs), -xtfrm(d$panel_time),
                 is.na(d$specimen_id), d$specimen_id), , drop = FALSE]
  } else {
    d <- d[order(d$hadm_id, -xtfrm(d$panel_time), is.na(d$specimen_id),
                 d$specimen_id), , drop = FALSE]
  }
  d[!duplicated(d$hadm_id), , drop = FALSE]
}

roc_auc <- function(y, x) {
  as.numeric(auc(roc(y, x, quiet = TRUE, direction = "<")))
}
paired_result <- function(d, free, containing, cmp, selection) {
  y <- as.integer(d$hospital_expire_flag)
  a <- roc(y, d[[free]], quiet = TRUE, direction = "<")
  b <- roc(y, d[[containing]], quiet = TRUE, direction = "<")
  tt <- roc.test(a, b, paired = TRUE, method = "delong")
  set.seed(20261118)
  i0 <- which(y == 0L); i1 <- which(y == 1L)
  boot <- replicate(2000L, {
    ix <- c(sample(i0, length(i0), TRUE), sample(i1, length(i1), TRUE))
    roc_auc(y[ix], d[[free]][ix]) - roc_auc(y[ix], d[[containing]][ix])
  })
  data.frame(selection = selection, comparison = cmp, n = nrow(d),
             events = sum(y), AUC_free = as.numeric(auc(a)),
             AUC_containing = as.numeric(auc(b)), delta_AUC = as.numeric(auc(a)-auc(b)),
             DeLong_p = tt$p.value,
             CI95_lower = unname(quantile(boot, .025)),
             CI95_upper = unname(quantile(boot, .975)), bootstrap_replicates = length(boot))
}

res <- list()
for (require_mono in c(FALSE, TRUE)) {
  comp <- if (require_mono) "SIRI vs PIV" else "NLR vs SII"
  free <- if (require_mono) "siri" else "nlr"
  containing <- if (require_mono) "piv" else "sii"
  for (prefer_abs in c(TRUE, FALSE)) {
    selection <- if (prefer_abs) "absolute-first then latest" else "latest eligible panel"
    d <- pick_panel(panels, require_mono, prefer_abs)
    res[[length(res)+1L]] <- paired_result(d, free, containing, comp, selection)
  }
}
selection <- do.call(rbind, res)
write.csv(selection, file.path(out_dir, "cbc_panel_selection_sensitivity.csv"), row.names = FALSE)

chosen_a <- pick_panel(panels, FALSE, TRUE)
latest_a <- pick_panel(panels, FALSE, FALSE)
chosen_b <- pick_panel(panels, TRUE, TRUE)
latest_b <- pick_panel(panels, TRUE, FALSE)
switch_count <- function(a, b) {
  m <- merge(a[, c("hadm_id", "panel_time", "all_abs")],
             b[, c("hadm_id", "panel_time", "all_abs")], by = "hadm_id",
             suffixes = c("_absolute_first", "_latest"))
  sum(m$panel_time_absolute_first != m$panel_time_latest |
        m$all_abs_absolute_first != m$all_abs_latest)
}
switch_tbl <- data.frame(comparison = c("NLR vs SII", "SIRI vs PIV"),
                         n_absolute_first = c(nrow(chosen_a), nrow(chosen_b)),
                         n_latest_panel = c(nrow(latest_a), nrow(latest_b)),
                         n_hadm_selecting_different_panel =
                           c(switch_count(chosen_a, latest_a),
                             switch_count(chosen_b, latest_b)))
write.csv(switch_tbl, file.path(out_dir, "cbc_panel_selection_counts.csv"), row.names = FALSE)
cat("CBC selection sensitivity:\n"); print(selection, digits = 5, row.names = FALSE)
cat("Panel selection differences:\n"); print(switch_tbl, row.names = FALSE)
