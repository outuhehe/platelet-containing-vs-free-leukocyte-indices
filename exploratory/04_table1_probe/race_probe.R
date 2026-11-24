# Probe: confirm the read-only MIMIC connection works and that admissions has a
# race column, and list the distinct race labels actually present. Read-only;
# the connection factory already asserts default_transaction_read_only = on.
.this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
.ra_root <- normalizePath(file.path(dirname(.this_file), "..", ".."), winslash = "/")
Sys.setenv(RA_ROOT = .ra_root)
suppressMessages(source(file.path(.ra_root, "R", "00_common.R")))
# Exploration can generate row-level records and observed ROC thresholds.
# Keep every exploratory artifact inside the restricted, git-ignored tree.
P_OUT <- file.path(P_INTERIM, "exploratory")
dir.create(P_OUT, recursive = TRUE, showWarnings = FALSE)

con <- connect_mimic()
on.exit(dbDisconnect(con), add = TRUE)

cols <- dbGetQuery(con, paste0(
  "SELECT ordinal_position, column_name, data_type FROM information_schema.columns ",
  "WHERE table_schema='mimiciv_hosp' AND table_name='admissions' ORDER BY ordinal_position"))
cat("admissions columns:\n"); print(cols)

stopifnot("race" %in% cols$column_name)

race_levels <- dbGetQuery(con, paste0(
  "SELECT race, COUNT(*) AS n FROM mimiciv_hosp.admissions ",
  "GROUP BY race ORDER BY n DESC"))
cat("distinct race labels across admissions:\n"); print(race_levels)
