# Load shared paths and the read-only connection factory (multi-line pgpass
# matching and verified read-only), replacing the former first-line-only read.
.this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
.ra_root <- normalizePath(file.path(dirname(.this_file), "..", ".."), winslash = "/")
Sys.setenv(RA_ROOT = .ra_root)
suppressMessages(source(file.path(.ra_root, "R", "00_common.R")))
# Exploration can generate row-level records and observed ROC thresholds.
# Keep every exploratory artifact inside the restricted, git-ignored tree.
P_OUT <- file.path(P_INTERIM, "exploratory")
dir.create(P_OUT, recursive = TRUE, showWarnings = FALSE)

con <- connect_eicu()
on.exit(dbDisconnect(con), add = TRUE)

schema <- as.data.table(dbGetQuery(con, paste0(
  "SELECT column_name,data_type FROM information_schema.columns ",
  "WHERE table_schema='eicu_crd' AND table_name='admissiondx' ",
  "ORDER BY ordinal_position"
)))
if (!all(c("patientunitstayid","admitdxpath","admitdxname") %in% schema$column_name)) {
  stop("admissiondx does not have the expected stay/path/name fields")
}
fwrite(schema, file.path(P_OUT, "eicu_admissiondx_schema.csv"))

units <- as.data.table(dbGetQuery(con, paste0(
  "SELECT unittype,COUNT(*) AS n_stay FROM eicu_crd.patient ",
  "GROUP BY unittype ORDER BY n_stay DESC"
)))
fwrite(units, file.path(P_OUT, "eicu_unittype_distribution.csv"))

gi_filter <- "admitdxpath LIKE 'admission diagnosis|All Diagnosis|Operative|Diagnosis|Gastrointestinal|%'"
gi_tags <- as.data.table(dbGetQuery(con, paste0(
  "SELECT admitdxpath,admitdxname,admitdxtext, ",
  "COUNT(DISTINCT patientunitstayid) AS n_stay ",
  "FROM eicu_crd.admissiondx WHERE ", gi_filter,
  " GROUP BY admitdxpath,admitdxname,admitdxtext ORDER BY n_stay DESC"
)))
fwrite(gi_tags, file.path(P_OUT, "eicu_gi_surgery_diagnosis_tags.csv"))

gi_units <- as.data.table(dbGetQuery(con, paste0(
  "WITH gi AS (SELECT DISTINCT patientunitstayid FROM eicu_crd.admissiondx WHERE ",
  gi_filter,
  ") SELECT p.unittype,COUNT(*) AS n_stay ",
  "FROM gi JOIN eicu_crd.patient p USING(patientunitstayid) ",
  "GROUP BY p.unittype ORDER BY n_stay DESC"
)))
fwrite(gi_units, file.path(P_OUT, "eicu_gi_surgery_by_unittype.csv"))

gi_offsets <- as.data.table(dbGetQuery(con, paste0(
  "WITH gi AS (SELECT DISTINCT patientunitstayid FROM eicu_crd.admissiondx WHERE ",
  gi_filter,
  ") SELECT COUNT(*) AS n_total, ",
  "COUNT(*) FILTER (WHERE p.hospitaldischargeoffset>0) AS n_hosp_offset_ok, ",
  "COUNT(*) FILTER (WHERE p.hospitaldischargeoffset>0 AND ",
  "p.hospitaldischargeoffset<100000) AS n_hosp_offset_sane, ",
  "COUNT(*) FILTER (WHERE p.hospitaldischargeoffset=0) AS n_offset_zero, ",
  "COUNT(*) FILTER (WHERE p.hospitaldischargeoffset>=100000) AS n_offset_ge_100000, ",
  "MIN(p.hospitaldischargeoffset) AS min_offset, ",
  "percentile_cont(0.5) WITHIN GROUP (ORDER BY p.hospitaldischargeoffset) AS median_offset, ",
  "MAX(p.hospitaldischargeoffset) AS max_offset ",
  "FROM gi JOIN eicu_crd.patient p USING(patientunitstayid)"
)))
fwrite(gi_offsets, file.path(P_OUT, "eicu_gi_surgery_offset_availability.csv"))

operative <- as.data.table(dbGetQuery(con, paste0(
  "SELECT admitdxpath,admitdxname,COUNT(DISTINCT patientunitstayid) AS n_stay ",
  "FROM eicu_crd.admissiondx WHERE admitdxpath LIKE 'admission diagnosis|Operative%' ",
  "GROUP BY admitdxpath,admitdxname ORDER BY n_stay DESC"
)))
from_or <- as.data.table(dbGetQuery(con, paste0(
  "SELECT admitdxpath,admitdxname,COUNT(DISTINCT patientunitstayid) AS n_stay ",
  "FROM eicu_crd.admissiondx WHERE admitdxpath ILIKE '%within 4 hours of admission%' ",
  "GROUP BY admitdxpath,admitdxname ORDER BY n_stay DESC"
)))
elective <- as.data.table(dbGetQuery(con, paste0(
  "SELECT admitdxpath,admitdxname,COUNT(DISTINCT patientunitstayid) AS n_stay ",
  "FROM eicu_crd.admissiondx WHERE admitdxpath LIKE 'admission diagnosis|Elective%' ",
  "GROUP BY admitdxpath,admitdxname ORDER BY n_stay DESC"
)))
fwrite(operative, file.path(P_OUT, "eicu_marker_operative_paths.csv"))
fwrite(from_or, file.path(P_OUT, "eicu_marker_from_or_paths.csv"))
fwrite(elective, file.path(P_OUT, "eicu_marker_elective_paths.csv"))

marker_summary <- as.data.table(dbGetQuery(con, paste0(
  "WITH flags AS (SELECT patientunitstayid,",
  "BOOL_OR(admitdxpath='admission diagnosis|Operative|Yes') AS op_yes,",
  "BOOL_OR(admitdxpath='admission diagnosis|Operative|No') AS op_no,",
  "BOOL_OR(admitdxpath='admission diagnosis|Was the patient admitted from the O.R. or went to the O.R. within 4 hours of admission?|Yes') AS or_yes,",
  "BOOL_OR(admitdxpath='admission diagnosis|Was the patient admitted from the O.R. or went to the O.R. within 4 hours of admission?|No') AS or_no,",
  "BOOL_OR(admitdxpath='admission diagnosis|Elective|Yes') AS elective_yes,",
  "BOOL_OR(admitdxpath='admission diagnosis|Elective|No') AS elective_no ",
  "FROM eicu_crd.admissiondx GROUP BY patientunitstayid),",
  "patient_flags AS (SELECT p.patientunitstayid, f.* FROM eicu_crd.patient p ",
  "LEFT JOIN flags f USING(patientunitstayid)),",
  "summary AS (SELECT 'Operative' AS marker,",
  "COUNT(*) FILTER (WHERE op_yes) AS n_yes,COUNT(*) FILTER (WHERE op_no) AS n_no,",
  "COUNT(*) FILTER (WHERE NOT COALESCE(op_yes,FALSE) AND NOT COALESCE(op_no,FALSE)) AS n_absent FROM patient_flags ",
  "UNION ALL SELECT 'FromOR',COUNT(*) FILTER (WHERE or_yes),COUNT(*) FILTER (WHERE or_no),",
  "COUNT(*) FILTER (WHERE NOT COALESCE(or_yes,FALSE) AND NOT COALESCE(or_no,FALSE)) FROM patient_flags ",
  "UNION ALL SELECT 'Elective',COUNT(*) FILTER (WHERE elective_yes),COUNT(*) FILTER (WHERE elective_no),",
  "COUNT(*) FILTER (WHERE NOT COALESCE(elective_yes,FALSE) AND NOT COALESCE(elective_no,FALSE)) FROM patient_flags) ",
  "SELECT * FROM summary"
)))
fwrite(marker_summary, file.path(P_OUT, "eicu_surgery_marker_summary.csv"))

all_gi <- as.data.table(dbGetQuery(con, paste0(
  "SELECT admitdxpath,admitdxname,COUNT(DISTINCT patientunitstayid) AS n_stay ",
  "FROM eicu_crd.admissiondx WHERE admitdxpath ILIKE '%gastrointestinal%' ",
  "GROUP BY admitdxpath,admitdxname ORDER BY n_stay DESC"
)))
fwrite(all_gi, file.path(P_OUT, "eicu_all_gastrointestinal_paths.csv"))

cat("admissiondx schema:\n"); print(schema)
cat("Full unit distribution:\n"); print(units)
cat("GI operative diagnosis labels:\n"); print(gi_tags)
cat("GI surgery by unit:\n"); print(gi_units)
cat("GI surgery discharge-offset availability:\n"); print(gi_offsets)
cat("Operative paths:\n"); print(operative)
cat("From OR paths:\n"); print(from_or)
cat("Elective paths:\n"); print(elective)
cat("Marker summary:\n"); print(marker_summary)
cat("All GI-containing paths:\n"); print(all_gi)
