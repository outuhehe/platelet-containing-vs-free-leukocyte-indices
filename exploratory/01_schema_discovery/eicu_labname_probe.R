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
print(dbGetQuery(con, paste0(
  "SELECT current_database() AS db, ",
  "current_setting('transaction_read_only') AS read_only, ",
  "to_regclass('eicu_crd.lab') AS lab"
)))
sample_sql <- paste0(
  "SELECT labname, labmeasurenamesystem, COUNT(*) AS n ",
  "FROM eicu_crd.lab WHERE ",
  "lower(labname) LIKE '%neut%' OR lower(labname) LIKE '%lymph%' ",
  "OR lower(labname) LIKE '%mono%' OR lower(labname) LIKE '%platelet%' ",
  "OR lower(labname) LIKE '%plt%' OR lower(labname) LIKE '%wbc%' ",
  "OR lower(labname) LIKE '%leuko%' OR lower(labname) LIKE '%-eos%' ",
  "GROUP BY labname,labmeasurenamesystem ORDER BY n DESC LIMIT 30"
)
print(dbGetQuery(con, sample_sql), row.names = FALSE)
