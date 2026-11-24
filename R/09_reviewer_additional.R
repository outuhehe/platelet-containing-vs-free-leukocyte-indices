# =============================================================================
# 09_reviewer_additional.R
# PLR vs inverse-lymphocyte paired DeLong with bootstrap CI, platelet AUC, and
# the all-absolute subgroup (n=409) paired CIs. Uses data.table instead of readr.
# Reads interim/cbc_export_v2.csv; writes tables/. Statistical logic unchanged.
# =============================================================================
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
Sys.setenv(RA_ROOT = ra_root)
source(file.path(ra_root, "R", "00_common.R"))
load_analysis_pkgs("pROC")

d <- fread(file.path(P_INTERIM, "cbc_export_v2.csv"), check.names = FALSE)
stopifnot(nrow(d) == 879L,
          sum(as.integer(d$hospital_expire_flag)) == 87L)

y <- as.integer(d$hospital_expire_flag)
alc_cols <- which(names(d) == "alc")
alc <- d[[alc_cols[[1L]]]]
plr <- d$plr
plt <- d$plt
stopifnot(all(is.finite(plr)), all(is.finite(alc)), all(is.finite(plt)))
stopifnot(all(plr > 0), all(alc > 0), all(plt > 0))
stopifnot(isTRUE(all.equal(plr, plt / alc, tolerance = 1e-8)))

r_plr <- roc(y, log(plr), quiet = TRUE, direction = "<")
r_inv_alc <- roc(y, -log(alc), quiet = TRUE, direction = "<")
r_plt <- roc(y, -log(plt), quiet = TRUE, direction = "<")
test <- roc.test(r_plr, r_inv_alc, paired = TRUE, method = "delong")

set.seed(20261116)
B <- 2000L
i0 <- which(y == 0L)
i1 <- which(y == 1L)
boot_delta <- replicate(B, {
  ix <- c(sample(i0, length(i0), replace = TRUE), sample(i1, length(i1), replace = TRUE))
  auc(roc(y[ix], log(plr[ix]), quiet = TRUE, direction = "<")) -
      auc(roc(y[ix], -log(alc[ix]), quiet = TRUE, direction = "<"))
})
ci_delta <- unname(quantile(boot_delta, c(0.025, 0.975)))
ci_plt <- ci.auc(r_plt, method = "delong")

out <- data.frame(
  comparison = "PLR vs inverse lymphocyte",
  AUC_PLR = as.numeric(auc(r_plr)),
  AUC_inverse_lymphocyte = as.numeric(auc(r_inv_alc)),
  delta_AUC = as.numeric(auc(r_plr) - auc(r_inv_alc)),
  DeLong_p = test$p.value,
  CI95_lower = ci_delta[[1]],
  CI95_upper = ci_delta[[2]],
  bootstrap_replicates = B,
  platelet_AUC = as.numeric(auc(r_plt)),
  platelet_CI95_lower = as.numeric(ci_plt[1]),
  platelet_CI95_upper = as.numeric(ci_plt[3]),
  n = nrow(d),
  events = sum(y)
)
print(out, digits = 5, row.names = FALSE)
fwrite(out, file.path(P_TABLES, "reviewer_additional_results.csv"))

sub <- d[d$all_abs == 1, ]
ys <- as.integer(sub$hospital_expire_flag)
stopifnot(nrow(sub) == 409L, sum(ys) == 61L)
sub_pairs <- list(
  list(name = "SII vs NLR", a = "sii", b = "nlr"),
  list(name = "PIV vs SIRI", a = "piv", b = "siri")
)
set.seed(20261116)
sub_out <- do.call(rbind, lapply(sub_pairs, function(pair) {
  ra <- roc(ys, log(sub[[pair$a]]), quiet = TRUE, direction = "<")
  rb <- roc(ys, log(sub[[pair$b]]), quiet = TRUE, direction = "<")
  tt <- roc.test(rb, ra, paired = TRUE, method = "delong")
  s0 <- which(ys == 0L)
  s1 <- which(ys == 1L)
  deltas <- replicate(B, {
    ix <- c(sample(s0, length(s0), replace = TRUE), sample(s1, length(s1), replace = TRUE))
    auc(roc(ys[ix], log(sub[[pair$a]][ix]), quiet = TRUE, direction = "<")) -
      auc(roc(ys[ix], log(sub[[pair$b]][ix]), quiet = TRUE, direction = "<"))
  })
  qs <- unname(quantile(deltas, c(0.025, 0.975)))
  data.frame(
    comparison = pair$name,
    AUC_with_platelet = as.numeric(auc(ra)),
    AUC_control = as.numeric(auc(rb)),
    delta_AUC = as.numeric(auc(ra) - auc(rb)),
    DeLong_p = tt$p.value,
    CI95_lower = qs[[1]],
    CI95_upper = qs[[2]],
    n = nrow(sub),
    events = sum(ys),
    bootstrap_replicates = B
  )
}))
print(sub_out, digits = 5, row.names = FALSE)
fwrite(sub_out, file.path(P_TABLES, "reviewer_all_absolute_subgroup_ci.csv"))
