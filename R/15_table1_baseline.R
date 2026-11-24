# =============================================================================
# 15_table1_baseline.R
# Build a standard baseline Table 1 (Overall / Survivors / Non-survivors) for the
# 879-stay complete five-index panel, adding race and admission_type and the
# descriptive statistics of the five inflammatory indices.
#
# Audit design (human-reviewed, locally executed):
#   * Reuses the frozen interim/cbc_export_v2.csv (879 rows / 87 deaths) produced
#     by 02_extract_index_panel.R; no heavy extraction chain is re-run.
#   * Opens the read-only MIMIC connection only to fetch race and admission_type
#     from mimiciv_hosp.admissions (the factory asserts read-only mode).
#   * The output is aggregate and contains no subject_id / hadm_id; still review
#     aggregate cells for disclosure risk before external release.
#
# Frozen assertions: 879 stays, 792 survivors, 87 non-survivors.
# =============================================================================
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
Sys.setenv(RA_ROOT = ra_root)
source(file.path(ra_root, "R", "00_common.R"))

panel_file <- file.path(P_INTERIM, "cbc_export_v2.csv")
if (!file.exists(panel_file))
  stop("interim/cbc_export_v2.csv not found; run 02_extract_index_panel.R first.")
panel <- as.data.table(fread(panel_file))

# cbc_export_v2.csv carries duplicated anc/alc/amc names (absolute counts and the
# percentage-derived columns share a name). Keep only the columns Table 1 needs;
# all kept names are unique, which also drops the duplicate-name conflict.
keep_cols <- c("hadm_id", "subject_id", "hospital_expire_flag",
               "nlr", "sii", "piv", "siri", "plr", "sex", "age", "charlson")
panel <- panel[, ..keep_cols]

# ---- frozen audit gate -------------------------------------------------------
if (nrow(panel) != 879L)
  stop("Panel does not match frozen row count 879.")
n_death <- sum(panel$hospital_expire_flag == 1L)
n_alive <- sum(panel$hospital_expire_flag == 0L)
if (n_death != 87L || n_alive != 792L)
  stop("Panel does not match frozen 792 survivors / 87 non-survivors.")
log_msg("Panel loaded:", nrow(panel), "stays |", n_alive, "survivors |", n_death, "deaths")

# ---- read-only query for race and admission_type ------------------------------
con <- connect_mimic()
on.exit(dbDisconnect(con), add = TRUE)

# hadm_id values are integers read from a trusted local file; safe to inline.
id_sql <- paste(panel$hadm_id, collapse = ",")
adm <- as.data.table(dbGetQuery(con, paste0(
  "SELECT hadm_id, race, admission_type FROM mimiciv_hosp.admissions ",
  "WHERE hadm_id IN (", id_sql, ")")))
log_msg("admissions rows fetched:", nrow(adm))
if (nrow(adm) != 879L)
  stop("Race/admission_type query did not return one row per panel stay.")

# ---- race grouping (explicit, human-reviewable keyword rules) -----------------
# MIMIC stores 30+ detailed labels. Collapse into five broad categories, checked
# in priority order:
#   Hispanic/Latino  : HISPANIC, LATINO, SOUTH AMERICAN
#   Black            : BLACK
#   Asian            : ASIAN
#   White            : WHITE, PORTUGUESE
#   Other/Unknown    : everything else (OTHER, UNKNOWN, UNABLE, DECLINED,
#                      MULTIPLE, AMERICAN INDIAN, PACIFIC ISLANDER)
classify_race <- function(x) {
  u <- toupper(trimws(x))
  out <- rep("Other/Unknown", length(u))
  out[grepl("WHITE|PORTUGUESE", u)] <- "White"
  out[grepl("ASIAN", u)] <- "Asian"
  out[grepl("BLACK", u)] <- "Black/African American"
  out[grepl("HISPANIC|LATINO|SOUTH AMERICAN", u)] <- "Hispanic/Latino"
  factor(out, levels = c("White", "Black/African American",
                        "Hispanic/Latino", "Asian", "Other/Unknown"))
}
adm[, race_group := classify_race(race)]

# ---- merge onto the panel (no identifiers retained afterwards) ---------------
dt <- merge(panel,
            adm[, .(hadm_id, race_group, admission_type)],
            by = "hadm_id", all.x = TRUE)
stopifnot(nrow(dt) == 879L, !anyNA(dt$race_group), !anyNA(dt$admission_type))
dt[, died := hospital_expire_flag == 1L]

# ---- formatting helpers -------------------------------------------------------
fmt_med_iqr <- function(v) {
  q <- quantile(v, probs = c(0.5, 0.25, 0.75), na.rm = TRUE, names = FALSE)
  sprintf("%s [%s, %s]",
          format(round(q[1], 1), nsmall = 1),
          format(round(q[2], 1), nsmall = 1),
          format(round(q[3], 1), nsmall = 1))
}
fmt_cat <- function(v, lev) {
  n_total <- sum(!is.na(v))
  n_lev <- sum(v == lev, na.rm = TRUE)
  sprintf("%d (%.1f)", n_lev, 100 * n_lev / n_total)
}

# Continuous variables: median [IQR]; Wilcoxon rank-sum p across outcome groups.
cont_vars <- c("age", "charlson", "nlr", "sii", "siri", "piv", "plr")
cont_rows <- rbindlist(lapply(cont_vars, function(v) {
  p_w <- suppressWarnings(
    wilcox.test(dt[[v]] ~ dt$died, exact = FALSE)$p.value)
  data.table(
    section = "Continuous", variable = v, level = "",
    overall = fmt_med_iqr(dt[[v]]),
    survivors = fmt_med_iqr(dt[[v]][!dt$died]),
    nonsurvivors = fmt_med_iqr(dt[[v]][dt$died]),
    p = ifelse(is.na(p_w), "", formatC(p_w, format = "g", digits = 3)))
}))

# Categorical variables: n (%); Fisher exact p (one test per variable).
set.seed(20261120)
build_cat_block <- function(v, levels) {
  p_f <- suppressWarnings(
    fisher.test(table(dt[[v]], dt$died), simulate.p.value = TRUE, B = 10000)$p.value)
  rbindlist(lapply(levels, function(lev) {
    data.table(
      section = "Categorical", variable = v, level = lev,
      overall = fmt_cat(dt[[v]], lev),
      survivors = fmt_cat(dt[[v]][!dt$died], lev),
      nonsurvivors = fmt_cat(dt[[v]][dt$died], lev),
      p = ifelse(lev == levels[1], formatC(p_f, format = "g", digits = 3), ""))
  }))
}
sex_block <- build_cat_block("sex", c("F", "M"))
race_block <- build_cat_block("race_group", levels(dt$race_group))
adm_levels <- sort(unique(dt$admission_type))
adm_block <- build_cat_block("admission_type", adm_levels)

# Order rows to match the planned Table 1.
out <- rbindlist(list(
  cont_rows[variable == "age"],
  sex_block,
  race_block,
  cont_rows[variable == "charlson"],
  adm_block,
  cont_rows[variable %in% c("nlr", "sii", "siri", "piv", "plr")]
), use.names = TRUE)

# Identifier-free, shareable table.
write_table(out, "table1_baseline_long")
log_msg("Written tables/table1_baseline_long.csv (", nrow(out), "rows; no subject/hadm IDs)")

cat("\n--- Table 1 preview (Overall n=879 / Survivors n=792 / Non-survivors n=87) ---\n")
print(out)
