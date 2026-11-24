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

schema <- dbGetQuery(con, paste0(
  "SELECT column_name,data_type FROM information_schema.columns ",
  "WHERE table_schema='eicu_crd' AND table_name IN ('lab','patient') ",
  "ORDER BY table_name,ordinal_position"
))
write.csv(schema, file.path(P_OUT, "eicu_lab_patient_schema.csv"), row.names = FALSE)

pattern_query <- paste0(
  "SELECT labname, labmeasurenamesystem AS unit, COUNT(*) AS n, ",
  "COUNT(DISTINCT patientunitstayid) AS n_stay, ",
  "MIN(labresult::numeric) AS min_v, MAX(labresult::numeric) AS max_v, ",
  "ROUND(AVG(labresult::numeric),2) AS mean_v ",
  "FROM eicu_crd.lab WHERE ",
  "lower(labname) LIKE '%neut%' OR lower(labname) LIKE '%lymph%' ",
  "OR lower(labname) LIKE '%mono%' OR lower(labname) LIKE '%platelet%' ",
  "OR lower(labname) LIKE '%plt%' OR lower(labname) LIKE '%wbc%' ",
  "OR lower(labname) LIKE '%leuko%' OR lower(labname) LIKE '%-eos%' ",
  "OR lower(labname) LIKE '%poly%' OR lower(labname) LIKE '%band%' ",
  "OR lower(labname) LIKE '%granul%' OR lower(labname) LIKE '%seg%' ",
  "GROUP BY labname,labmeasurenamesystem ORDER BY n DESC"
)
patterns <- as.data.table(dbGetQuery(con, pattern_query))
fwrite(patterns, file.path(P_OUT, "eicu_cbc_pattern_coverage.csv"))

negative <- as.data.table(dbGetQuery(con, paste0(
  "SELECT labname,labmeasurenamesystem AS unit,COUNT(*) AS n, ",
  "COUNT(DISTINCT patientunitstayid) AS n_stay ",
  "FROM eicu_crd.lab WHERE labname LIKE '-%' ",
  "GROUP BY labname,labmeasurenamesystem ORDER BY n DESC LIMIT 50"
)))
fwrite(negative, file.path(P_OUT, "eicu_negative_labnames.csv"))

units <- as.data.table(dbGetQuery(con, paste0(
  "SELECT unittype,COUNT(*) AS n,COUNT(DISTINCT patientunitstayid) AS n_stay ",
  "FROM eicu_crd.patient GROUP BY unittype ORDER BY n DESC"
)))
fwrite(units, file.path(P_OUT, "eicu_unittype_distribution.csv"))

cat("Lab schema:\n"); print(schema)
cat("CBC-related labnames:\n"); print(patterns)
cat("Negative-prefix names (top 50):\n"); print(negative)
cat("Unit types:\n"); print(units)
