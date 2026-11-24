# Smoke test for 00_common.R: locate this file, source common, verify both DBs.
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
Sys.setenv(RA_ROOT = ra_root)
source(file.path(ra_root, "R", "00_common.R"))

cat("RA_ROOT:", RA_ROOT, "\n")
cat("user:", PGUSER, "| host:", PGHOST, "| port:", PGPORT,
    "| mimic:", MIMIC_DB, "| eicu:", EICU_DB, "| work_mem:", PG_WORK_MEM, "\n")

cm <- connect_mimic()
cat("MIMIC read_only:", dbGetQuery(cm, "SHOW default_transaction_read_only")[[1]],
    "| admissions:", dbGetQuery(cm, "SELECT count(*)::int n FROM mimiciv_hosp.admissions")$n, "\n")
dbDisconnect(cm)

ce <- connect_eicu()
cat("EICU read_only:", dbGetQuery(ce, "SHOW default_transaction_read_only")[[1]],
    "| admissiondx:", dbGetQuery(ce, "SELECT count(*)::int n FROM eicu_crd.admissiondx")$n, "\n")
dbDisconnect(ce)
cat("SMOKE OK\n")
