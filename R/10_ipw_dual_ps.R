# =============================================================================
# 10_ipw_dual_ps.R
# Dual-propensity-score IPW panel-availability analysis. Baseline models are
# primary; outcome-conditioned models are a sensitivity specification. Reports
# availability-model C statistics, PS overlap/extremes, SMD balance, weighted
# Mann-Whitney AUC and full-cohort-refit bootstrap CIs. Base R only.
# Reads roster/index/codebooks from interim; writes all results to tables/.
# Set IPW_RUN_BOOTSTRAP=false to refresh estimates without the 2000-rep bootstrap.
# =============================================================================
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
Sys.setenv(RA_ROOT = ra_root)
source(file.path(ra_root, "R", "00_common.R"))

## ----------------------------- paths / setup -------------------------------
out_dir <- P_TABLES
roster_f <- file.path(P_INTERIM, "cohort_roster_2376.csv")
index_candidates <- c(
  file.path(P_INTERIM, "cbc_available_case_index.csv"),
  file.path(P_INTERIM, "cbc_export_v2.csv")
)
icd9_cb  <- file.path(P_INTERIM, "D2_ICD9_45_46_codebook.csv")
icd10_cb <- file.path(P_INTERIM, "ICD10_actual_target_codes.csv")

B <- as.integer(Sys.getenv("IPW_BOOTSTRAP_REPLICATES", "2000"))
RUN_BOOTSTRAP <- tolower(Sys.getenv("IPW_RUN_BOOTSTRAP", "true")) %in% c("1", "true", "yes")
set.seed(20261117)

## ------------------------------ small helpers ------------------------------
pick_col <- function(df, candidates) {
  hits <- candidates[candidates %in% names(df)]
  if (length(hits) == 0) return(NULL)
  hits[1]
}
to01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  xl <- toupper(trimws(as.character(x)))
  if (all(xl %in% c("TRUE", "FALSE", "T", "F", "1", "0", "")))
    return(as.integer(xl %in% c("TRUE", "T", "1")))
  as.integer(x)
}
smd_pair <- function(x, w, group) {
  wm <- function(xx, ww) weighted.mean(xx, ww)
  wv <- function(xx, ww) { m <- wm(xx, ww); sum(ww * (xx - m)^2) / sum(ww) }
  if (is.numeric(x)) {
    m1 <- wm(x[group == 1], w[group == 1]); m0 <- wm(x[group == 0], w[group == 0])
    v1 <- wv(x[group == 1], w[group == 1]); v0 <- wv(x[group == 0], w[group == 0])
    s <- sqrt((v1 + v0) / 2)
    return(ifelse(s == 0, 0, (m1 - m0) / s))
  }
  NA
}

## ---------------------- organ-class mapping function -----------------------
icd10_bp_organ <- function(bp) {
  if (bp %in% c("0", "1", "2", "5"))            return("liver")
  if (bp == "4")                                return("biliary")
  if (bp %in% c("6", "7"))                      return("stomach")
  if (bp %in% c("8", "9", "A", "B"))            return("small_intestine")
  if (bp %in% c("C", "E", "F", "H", "K", "L", "M", "N", "P")) return("colorectal")
  if (bp %in% c("D", "G"))                      return("pancreas")
  NA_character_
}
icd9_prefix_organ <- function(code) {
  c <- gsub("\\.", "", trimws(as.character(code)))
  if (!grepl("^[0-9]+$", c)) return(NA_character_)
  blk <- suppressWarnings(as.integer(substr(c, 1, 2)))
  if (is.na(blk)) return(NA_character_)
  frac <- suppressWarnings(
    as.numeric(paste0(substr(c, 1, 2), ".", substr(c, 3, nchar(c)))))
  if (blk == 41) return("spleen")
  if (blk %in% c(42, 43, 44)) return("stomach")
  if (blk == 48) return("colorectal")
  if (blk == 50) return("liver")
  if (blk == 51) return("biliary")
  if (blk == 52) return("pancreas")
  if (!is.na(frac) && blk == 45) {
    if (frac >= 45.6 && frac < 45.7) return("small_intestine")
    if (frac >= 45.7 && frac < 45.9) return("colorectal")
    if (frac >= 45.9)                return("small_intestine")
  }
  if (!is.na(frac) && blk == 46) {
    if (frac >= 46.1 && frac < 46.2) return("colorectal")
    if (frac >= 46.2 && frac < 46.3) return("small_intestine")
    if (frac >= 46.4 && frac < 46.5) return("colorectal")
    if (frac >= 46.7)                return("small_intestine")
  }
  NA_character_
}
classify_organ <- function(code, version, title_map) {
  c <- trimws(as.character(code))
  v <- trimws(as.character(version))
  if (c %in% names(title_map)) {
    t <- tolower(title_map[[c]])
    organ <- if (grepl("liver|hepatic", t)) "liver"
      else if (grepl("gallbladder|bile duct|biliary", t)) "biliary"
      else if (grepl("pancrea", t)) "pancreas"
      else if (grepl("spleen|splenic", t)) "spleen"
      else if (grepl("stomach|gastric|pylor", t)) "stomach"
      else if (grepl("duoden|jejun|ileum|ileal|small intestine|small bowel", t)) "small_intestine"
      else if (grepl("colon|colonic|rect|cecum|cecal|sigmoid|large intestine|large bowel", t)) "colorectal"
      else NA_character_
    if (!is.na(organ)) return(organ)
  }
  if (identical(v, "10") && nchar(c) >= 4) {
    organ <- icd10_bp_organ(substr(c, 4, 4))
    if (!is.na(organ)) return(organ)
  }
  organ <- icd9_prefix_organ(c)
  if (!is.na(organ)) return(organ)
  "other"
}

## ----------------------------- weighted AUC --------------------------------
weighted_auc <- function(score, y, w, direction = ">") {
  ok <- is.finite(score) & is.finite(w) & w > 0
  score <- score[ok]; y <- y[ok]; w <- w[ok]
  ip <- which(y == 1); in0 <- which(y == 0)
  if (length(ip) == 0 || length(in0) == 0) return(NA_real_)
  diff <- outer(score[ip], score[in0], "-")
  ww <- outer(w[ip], w[in0])
  num <- sum(ww * (diff > 0)) + 0.5 * sum(ww * (diff == 0))
  u <- num / (sum(w[ip]) * sum(w[in0]))
  if (direction == ">") u else 1 - u
}
auc_binary <- function(score, y) {
  ok <- is.finite(score) & !is.na(y)
  score <- score[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(score, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## ------------------------------- load data ---------------------------------
roster <- read.csv(roster_f, stringsAsFactors = FALSE, check.names = FALSE)
index_f <- index_candidates[file.exists(index_candidates)][1]
if (is.na(index_f))
  stop("No available-case index table found; expected one of: ",
       paste(index_candidates, collapse = "; "))
idx <- read.csv(index_f, stringsAsFactors = FALSE, check.names = FALSE)

c_id    <- pick_col(roster, c("hadm_id", "hadm"))
c_age   <- pick_col(roster, c("age"))
c_charl <- pick_col(roster, c("charlson", "charlson_index"))
c_sex   <- pick_col(roster, c("sex"))
c_adm   <- pick_col(roster, c("admission_type", "admission"))
c_ver   <- pick_col(roster, c("last_icd_version", "icd_version"))
c_code  <- pick_col(roster, c("last_icd_code", "icd_code"))
c_hpanel<- pick_col(roster, c("has_any_panel"))
c_five  <- pick_col(roster, c("complete_five_index"))
c_death <- pick_col(roster, c("in_hospital_death", "hospital_expire_flag", "death", "hosp_death"))

roster$age   <- as.numeric(roster[[c_age]])
roster$charl <- as.numeric(roster[[c_charl]])
roster$sex   <- as.factor(roster[[c_sex]])
roster$adm   <- as.factor(roster[[c_adm]])
roster$ver   <- as.factor(roster[[c_ver]])
roster$has_panel <- to01(roster[[c_hpanel]])
roster$five  <- to01(roster[[c_five]])
roster$death <- to01(roster[[c_death]])
roster$id_key <- as.character(roster[[c_id]])

title_map <- list()
if (file.exists(icd10_cb)) {
  d <- read.csv(icd10_cb, stringsAsFactors = FALSE, check.names = FALSE)
  cc <- pick_col(d, c("icd_code", "code")); tt <- pick_col(d, c("long_title", "title"))
  for (i in seq_len(nrow(d))) title_map[[trimws(as.character(d[[cc]][i]))]] <- d[[tt]][i]
}
if (file.exists(icd9_cb)) {
  d <- read.csv(icd9_cb, stringsAsFactors = FALSE, check.names = FALSE)
  cc <- pick_col(d, c("icd_code", "code")); tt <- pick_col(d, c("long_title", "title"))
  for (i in seq_len(nrow(d))) title_map[[trimws(as.character(d[[cc]][i]))]] <- d[[tt]][i]
}
roster$organ <- as.factor(mapply(classify_organ,
                                 roster[[c_code]], roster[[c_ver]],
                                 MoreArgs = list(title_map = title_map)))

i_id <- pick_col(idx, c("hadm_id", "hadm"))
i_nlr <- pick_col(idx, c("nlr", "NLR")); i_sii <- pick_col(idx, c("sii", "SII"))
i_siri <- pick_col(idx, c("siri", "SIRI")); i_piv <- pick_col(idx, c("piv", "PIV"))
i_y <- pick_col(idx, c("hospital_expire_flag", "in_hospital_death", "hosp_death", "death"))
if (is.null(i_nlr) || is.null(i_sii) || is.null(i_siri) || is.null(i_piv))
  stop("Index table must contain nlr, sii, siri and piv columns")
idx$id_key <- as.character(idx[[i_id]])
idx$y <- to01(idx[[i_y]])
for (cc in c(i_nlr, i_sii, i_siri, i_piv)) idx[[cc]] <- as.numeric(idx[[cc]])

cat("Roster:", nrow(roster), "| has_panel:", sum(roster$has_panel),
    "| five:", sum(roster$five), "\n")
cat("Organ levels:", paste(levels(roster$organ), collapse = ", "), "\n")
cat("Index table:", basename(index_f), "rows:", nrow(idx),
    "| siri non-missing:", sum(is.finite(idx[[i_siri]]), na.rm = TRUE), "\n")

## --------------------- fit one availability PS + weights --------------------
fit_weights <- function(dat, avail, include_outcome = FALSE) {
  rhs <- "age + charl + sex + adm + ver + organ"
  if (include_outcome) rhs <- paste(rhs, "+ death")
  form <- as.formula(paste0(avail, " ~ ", rhs))
  fit <- glm(form, data = dat, family = binomial())
  ps <- predict(fit, type = "response")
  a <- dat[[avail]]
  p1 <- mean(a == 1)
  sw <- ifelse(a == 1, p1 / ps, (1 - p1) / (1 - ps))
  list(fit = fit, ps = ps, sw = sw)
}

balance_smd <- function(dat, fw) {
  out <- data.frame()
  num <- list(age = dat$age, charlson = dat$charl)
  for (nm in names(num)) {
    out <- rbind(out, data.frame(
      covariate = nm,
      smd_unweighted = smd_pair(num[[nm]], rep(1, nrow(dat)), dat$A),
      smd_weighted = smd_pair(num[[nm]], fw$sw, dat$A)))
  }
  facs <- list(sex = dat$sex, admission_type = dat$adm,
               icd_version = dat$ver, organ = dat$organ)
  for (nm in names(facs)) {
    for (L in levels(facs[[nm]])[-1]) {
      dummy <- as.integer(facs[[nm]] == L)
      out <- rbind(out, data.frame(
        covariate = paste0(nm, ": ", L),
        smd_unweighted = smd_pair(dummy, rep(1, nrow(dat)), dat$A),
        smd_weighted = smd_pair(dummy, fw$sw, dat$A)))
    }
  }
  out
}

run_one <- function(dat, avail, free_col, contain_col, label,
                    include_outcome = FALSE) {
  dat$A <- dat[[avail]]
  fw <- fit_weights(dat, avail, include_outcome)
  a <- dat$A

  diagnostics <- data.frame(
    comparison = label,
    availability_target = avail,
    n_total = nrow(dat),
    n_available = sum(a == 1),
    n_unavailable = sum(a == 0),
    availability_model_c = auc_binary(fw$ps, a),
    ps_min = min(fw$ps), ps_max = max(fw$ps),
    sw_min = min(fw$sw[a == 1]), sw_median = median(fw$sw[a == 1]),
    sw_mean = mean(fw$sw[a == 1]), sw_max = max(fw$sw[a == 1]),
    n_extreme_gt10 = sum(fw$sw[a == 1] > 10))

  bal <- balance_smd(dat, fw); bal$comparison <- label

  avail_rows <- dat[dat$A == 1, ]
  avail_rows$sw_col <- fw$sw[dat$A == 1]
  m <- merge(avail_rows[, c("id_key", "sw_col"), drop = FALSE], idx,
             by = "id_key", sort = FALSE)
  ww <- m$sw_col
  auc_free <- weighted_auc(m[[free_col]], m$y, ww, direction = ">")
  auc_cont <- weighted_auc(m[[contain_col]], m$y, ww, direction = ">")
  primary <- data.frame(
    comparison = label, n = nrow(m), events = sum(m$y),
    weighted_AUC_free = auc_free, weighted_AUC_containing = auc_cont,
    weighted_delta_AUC = auc_free - auc_cont)

  list(diagnostics = diagnostics, balance = bal, primary = primary)
}

resA <- run_one(roster, "has_panel", i_nlr, i_sii, "NLR_vs_SII")
resB <- run_one(roster, "five", i_siri, i_piv, "SIRI_vs_PIV")

bootstrap_ci <- function(dat, avail, free_col, contain_col, label,
                         include_outcome = FALSE) {
  deltas <- rep(NA_real_, B); fails <- 0
  for (b in seq_len(B)) {
    s <- sample(seq_len(nrow(dat)), nrow(dat), replace = TRUE)
    db <- dat[s, ]
    db$sex <- factor(db$sex, levels = levels(dat$sex))
    db$adm <- factor(db$adm, levels = levels(dat$adm))
    db$ver <- factor(db$ver, levels = levels(dat$ver))
    db$organ <- factor(db$organ, levels = levels(dat$organ))
    ok <- tryCatch({
      fwb <- fit_weights(db, avail, include_outcome)
      db$A <- db[[avail]]
      avail_db <- db[db$A == 1, ]
      avail_db$sw_col <- fwb$sw[db$A == 1]
      mb <- merge(avail_db[, c("id_key", "sw_col"), drop = FALSE], idx,
                  by = "id_key", sort = FALSE)
      ww <- mb$sw_col
      if (length(unique(mb$y)) < 2 || sum(ww) <= 0) stop("single class")
      af <- weighted_auc(mb[[free_col]], mb$y, ww, ">")
      ac <- weighted_auc(mb[[contain_col]], mb$y, ww, ">")
      if (!is.finite(af) || !is.finite(ac)) stop("non-finite AUC")
      deltas[b] <- af - ac; TRUE
    }, error = function(e) { fails <<- fails + 1; FALSE })
    if (b %% 250 == 0) cat("  ", label, "replicate", b, "\n")
  }
  d <- deltas[!is.na(deltas)]
  data.frame(comparison = label, bootstrap_replicates = B, failed = fails,
             n_used = length(d),
             CI95_lower = as.numeric(quantile(d, .025)),
             bootstrap_mean = mean(d),
             CI95_upper = as.numeric(quantile(d, .975)))
}

if (RUN_BOOTSTRAP) {
  cat("Bootstrapping PS_A (NLR vs SII) ...\n")
  bsA <- bootstrap_ci(roster, "has_panel", i_nlr, i_sii, "NLR_vs_SII")
  cat("Bootstrapping PS_B (SIRI vs PIV) ...\n")
  bsB <- bootstrap_ci(roster, "five", i_siri, i_piv, "SIRI_vs_PIV")
} else {
  bsA <- bsB <- data.frame()
}

resA_y <- run_one(roster, "has_panel", i_nlr, i_sii,
                  "NLR_vs_SII_outcome_conditioned", include_outcome = TRUE)
resB_y <- run_one(roster, "five", i_siri, i_piv,
                  "SIRI_vs_PIV_outcome_conditioned", include_outcome = TRUE)
if (RUN_BOOTSTRAP) {
  bsA_y <- bootstrap_ci(roster, "has_panel", i_nlr, i_sii,
                        "NLR_vs_SII_outcome_conditioned", include_outcome = TRUE)
  bsB_y <- bootstrap_ci(roster, "five", i_siri, i_piv,
                        "SIRI_vs_PIV_outcome_conditioned", include_outcome = TRUE)
} else {
  bsA_y <- bsB_y <- data.frame()
}

## ------------------------------- write out ---------------------------------
write.csv(rbind(resA$diagnostics, resB$diagnostics),
          file.path(out_dir, "ipw_diagnostics.csv"), row.names = FALSE)
write.csv(rbind(resA$balance, resB$balance),
          file.path(out_dir, "ipw_balance.csv"), row.names = FALSE)
write.csv(rbind(resA$primary, resB$primary),
          file.path(out_dir, "ipw_weighted_primary.csv"), row.names = FALSE)
if (RUN_BOOTSTRAP) write.csv(rbind(bsA, bsB),
          file.path(out_dir, "ipw_bootstrap_ci.csv"), row.names = FALSE)
write.csv(rbind(resA_y$primary, resB_y$primary),
          file.path(out_dir, "ipw_outcome_conditioned_primary.csv"), row.names = FALSE)
if (RUN_BOOTSTRAP) write.csv(rbind(bsA_y, bsB_y),
          file.path(out_dir, "ipw_outcome_conditioned_bootstrap_ci.csv"), row.names = FALSE)
write.csv(rbind(resA_y$diagnostics, resB_y$diagnostics),
          file.path(out_dir, "ipw_outcome_conditioned_diagnostics.csv"), row.names = FALSE)
write.csv(rbind(resA$diagnostics, resB$diagnostics)[,
            c("comparison", "n_total", "n_available", "n_unavailable", "availability_model_c")],
          file.path(out_dir, "ipw_availability_model_cstat.csv"), row.names = FALSE)
write.csv(rbind(resA_y$diagnostics, resB_y$diagnostics)[,
            c("comparison", "n_total", "n_available", "n_unavailable", "availability_model_c")],
          file.path(out_dir, "ipw_outcome_conditioned_model_cstat.csv"), row.names = FALSE)

cat("\n==== weighted primary ====\n"); print(rbind(resA$primary, resB$primary))
cat("\nMax |weighted SMD| across both models:",
    max(abs(rbind(resA$balance, resB$balance)$smd_weighted), na.rm = TRUE), "\n")
