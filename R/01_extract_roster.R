# =============================================================================
# 01_extract_roster.R
# Execute sql/cbc_available_case.sql through the read-only MIMIC connection and
# materialize row-level interim files with stable pseudo-identifiers:
#   Query A -> interim/cbc_available_case_panels.csv  (candidate panels)
#   Query B -> interim/cohort_roster_2376.csv         (full eligible roster)
# Frozen assertions (audit gate): 2376 rows / 181 deaths / 879 complete
# five-index cases / 87 complete-case deaths.
# These intermediates remain access-controlled under the PhysioNet DUA and must
# stay local; credentials come from config/db.env or pgpass, never this script.
# =============================================================================
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
Sys.setenv(RA_ROOT = ra_root)
source(file.path(ra_root, "R", "00_common.R"))

sql_file <- file.path(P_SQL, "cbc_available_case.sql")
lines <- readLines(sql_file, warn = FALSE, encoding = "UTF-8")
mark_a <- grep("^-- QUERY A:", lines)[1]
mark_b <- grep("^-- QUERY B:", lines)[1]
if (is.na(mark_a) || is.na(mark_b) || mark_b <= mark_a)
  stop("Could not find Query A / Query B markers in ", sql_file)

query_a <- paste(lines[seq_len(mark_b - 1L)], collapse = "\n")
shared_ctes <- paste(lines[seq_len(mark_a - 1L)], collapse = "\n")
query_b <- paste(shared_ctes, "
SELECT m.subject_id, m.hadm_id, m.hospital_expire_flag, m.admission_type,
       p.gender AS sex,
       p.anchor_age + (m.adm_year - p.anchor_year) AS age,
       COALESCE(ch.charlson_comorbidity_index,0) AS charlson,
       bl.last_icd_code, bl.last_icd_version,
       COALESCE(f.has_any_panel,FALSE) AS has_any_panel,
       COALESCE(f.complete_five_index,FALSE) AS complete_five_index
FROM main m
JOIN mimiciv_hosp.patients p ON p.subject_id = m.subject_id
LEFT JOIN abs_last bl ON bl.hadm_id = m.hadm_id
LEFT JOIN mimiciv_derived.charlson ch ON ch.hadm_id = m.hadm_id
LEFT JOIN panel_flags f ON f.hadm_id = m.hadm_id
ORDER BY m.hadm_id;", sep = "\n")

con <- connect_mimic()
on.exit(dbDisconnect(con), add = TRUE)

log_msg("Running Query A: candidate CBC panels")
panels <- dbGetQuery(con, query_a)
log_msg("Query A rows:", nrow(panels), "| unique hadm:", length(unique(panels$hadm_id)))
fwrite(panels, file.path(P_INTERIM, "cbc_available_case_panels.csv"))
rm(panels); gc()

log_msg("Running Query B: full eligible roster")
roster <- dbGetQuery(con, query_b)
n_row <- nrow(roster); n_death <- sum(roster$hospital_expire_flag)
n_complete <- sum(roster$complete_five_index)
n_complete_death <- sum(roster$hospital_expire_flag[roster$complete_five_index])
log_msg("Query B rows:", n_row, "| deaths:", n_death,
        "| complete five-index:", n_complete, "| complete deaths:", n_complete_death)

if (n_row != 2376L || n_death != 181L || n_complete != 879L || n_complete_death != 87L)
  stop("Roster does not match frozen counts (2376/181/879/87); inspect interim CSV")
fwrite(roster, file.path(P_INTERIM, "cohort_roster_2376.csv"))
log_msg("Roster written and frozen cohort counts verified.")

# ---------------------------------------------------------------------------
# Auxiliary codebooks and the 886-row available-case index, all derived from
# the same shared CTE chain on the existing connection (no extra table scan of
# the raw procedures table beyond what proc_all already materializes).
# ---------------------------------------------------------------------------

# (1) ICD-9 45/46 codebook: distinct codes with long_title and full-DB counts.
q_icd9_cb <- paste0(shared_ctes, "
SELECT p.icd_code, d.long_title, COUNT(*) AS n
FROM proc_all p
LEFT JOIN mimiciv_hosp.d_icd_procedures d
  ON d.icd_code = p.icd_code AND d.icd_version = 9
WHERE p.icd_code_version = 9
  AND (p.icd_code LIKE '45%' OR p.icd_code LIKE '46%')
GROUP BY p.icd_code, d.long_title
ORDER BY p.icd_code;")
icd9_cb <- as.data.table(dbGetQuery(con, q_icd9_cb))
stopifnot(nrow(icd9_cb) == 88L)
fwrite(icd9_cb, file.path(P_INTERIM, "D2_ICD9_45_46_codebook.csv"))

# (2) ICD-10 actual target codes: counts plus parsed root/body-part/approach.
q_icd10_cb <- paste0(shared_ctes, "
SELECT s.icd_code, d.long_title, COUNT(*) AS n_proc,
       SUBSTRING(s.icd_code FROM 3 FOR 1) AS root_operation,
       SUBSTRING(s.icd_code FROM 4 FOR 1) AS body_part,
       SUBSTRING(s.icd_code FROM 5 FOR 1) AS approach
FROM abs_surg s
JOIN mimiciv_hosp.d_icd_procedures d
  ON d.icd_code = s.icd_code AND d.icd_version = 10
WHERE s.icd_code_version = 10
GROUP BY s.icd_code, d.long_title
ORDER BY s.icd_code;")
icd10_cb <- as.data.table(dbGetQuery(con, q_icd10_cb))
stopifnot(nrow(icd10_cb) == 195L)
fwrite(icd10_cb, file.path(P_INTERIM, "ICD10_actual_target_codes.csv"))

# (3) 886-row index. Prefer the strict (amc_f>0) panel so the 879 frozen rows
#     match cbc_export_v2; the 7 stays without a strict panel contribute
#     panel-derived NLR/SII only (siri/piv are NA).
q_index <- paste0(shared_ctes, "
, strict_cand AS (
  SELECT r.*, ROW_NUMBER() OVER (PARTITION BY hadm_id
    ORDER BY all_abs DESC, panel_time DESC, specimen_id NULLS LAST) AS srn
  FROM resolved r
  WHERE anc_f > 0 AND alc_f > 0 AND amc_f > 0 AND plt > 0
),
strict_pick AS (SELECT * FROM strict_cand WHERE srn = 1),
loose_pick AS (SELECT * FROM best WHERE core_rank = 1),
combined AS (
  SELECT l.hadm_id, l.hospital_expire_flag,
    COALESCE(s.anc_f, l.anc_f) AS anc_f,
    COALESCE(s.alc_f, l.alc_f) AS alc_f,
    COALESCE(s.amc_f, l.amc_f) AS amc_f,
    COALESCE(s.plt, l.plt) AS plt,
    CASE WHEN s.hadm_id IS NOT NULL THEN 1 ELSE 0 END AS has_strict
  FROM loose_pick l
  LEFT JOIN strict_pick s ON s.hadm_id = l.hadm_id
)
SELECT hadm_id,
  anc_f / alc_f AS nlr,
  plt * anc_f / alc_f AS sii,
  CASE WHEN has_strict = 1 THEN anc_f * amc_f / alc_f END AS siri,
  CASE WHEN has_strict = 1 THEN plt * anc_f * amc_f / alc_f END AS piv,
  hospital_expire_flag,
  CASE WHEN has_strict = 1 THEN 'frozen_export' ELSE 'panel_derived' END AS source
FROM combined
ORDER BY hadm_id;")
case_index <- as.data.table(dbGetQuery(con, q_index))
stopifnot(nrow(case_index) == 886L,
          sum(case_index$source == "frozen_export") == 879L,
          sum(case_index$source == "panel_derived") == 7L,
          sum(is.na(case_index$siri)) == 7L,
          sum(is.na(case_index$piv)) == 7L)
fwrite(case_index, file.path(P_INTERIM, "cbc_available_case_index.csv"), na = "")
log_msg("Auxiliary codebooks and 886-row index written and verified.")
