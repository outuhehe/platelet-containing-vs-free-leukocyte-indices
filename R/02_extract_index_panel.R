# =============================================================================
# 02_extract_index_panel.R
# Execute sql/cbc_export_v2.sql through the read-only MIMIC connection and
# materialize the 879-row index panel (all five inflammatory indices plus the
# single-marker ANC and platelets).
# This closes the original reproducibility gap: cbc_export_v2.sql previously had
# no scripted runner and was executed manually in a database GUI.
# Frozen assertion: 879 rows / 87 deaths, with the index columns used below.
# =============================================================================
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
Sys.setenv(RA_ROOT = ra_root)
source(file.path(ra_root, "R", "00_common.R"))

sql_file <- file.path(P_SQL, "cbc_export_v2.sql")
sql <- paste(readLines(sql_file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

required_cols <- c("nlr", "sii", "piv", "siri", "plr", "anc_single",
                   "plt", "hospital_expire_flag")

con <- connect_mimic()
on.exit(dbDisconnect(con), add = TRUE)

log_msg("Running cbc_export_v2.sql (full WITH chain)")
dt <- dbGetQuery(con, sql)
log_msg("Index panel rows:", nrow(dt), "| columns:", ncol(dt))

if (nrow(dt) != 879L)
  stop("Index panel does not match frozen row count 879; inspect sql/cbc_export_v2.sql")
if (sum(as.integer(dt$hospital_expire_flag)) != 87L)
  stop("Index panel does not match frozen death count 87")
absent <- setdiff(required_cols, names(dt))
if (length(absent) > 0L)
  stop("Index panel missing required columns: ", paste(absent, collapse = ", "))

fwrite(dt, file.path(P_INTERIM, "cbc_export_v2.csv"))
log_msg("Index panel written to interim/cbc_export_v2.csv (879 rows / 87 deaths verified).")
