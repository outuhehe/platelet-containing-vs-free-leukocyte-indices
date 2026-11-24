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
  "WHERE table_schema='eicu_crd' AND table_name='patient' ",
  "AND column_name IN ('unittype','unitvisitnumber','age',",
  "'hospitaldischargelocation','hospitaldischargeoffset','hospitaldischargestatus') ",
  "ORDER BY ordinal_position"
)))
units <- as.data.table(dbGetQuery(con, paste0(
  "SELECT unittype,COUNT(*) AS n_stay FROM eicu_crd.patient ",
  "GROUP BY unittype ORDER BY n_stay DESC"
)))
visit <- as.data.table(dbGetQuery(con, paste0(
  "SELECT unittype,unitvisitnumber,COUNT(*) AS n_stay ",
  "FROM eicu_crd.patient GROUP BY unittype,unitvisitnumber ",
  "ORDER BY unittype,unitvisitnumber"
)))
locations <- as.data.table(dbGetQuery(con, paste0(
  "SELECT hospitaldischargelocation,COUNT(*) AS n_stay,",
  "COUNT(*) FILTER (WHERE hospitaldischargeoffset>0) AS n_offset_positive ",
  "FROM eicu_crd.patient WHERE unittype='SICU' ",
  "GROUP BY hospitaldischargelocation ORDER BY n_stay DESC NULLS LAST"
)))
panel <- as.data.table(dbGetQuery(con, paste0(
  "WITH panel AS MATERIALIZED (SELECT patientunitstayid,labresultoffset,",
  "COUNT(DISTINCT labname) AS n_items FROM eicu_crd.lab ",
  "WHERE labname IN ('-polys','-lymphs','-monos','WBC x 1000','platelets x 1000') ",
  "AND labresultoffset BETWEEN 0 AND 1440 ",
  "GROUP BY patientunitstayid,labresultoffset) ",
  "SELECT n_items,COUNT(*) AS n_panels,",
  "COUNT(DISTINCT patientunitstayid) AS n_stay ",
  "FROM panel GROUP BY n_items ORDER BY n_items"
)))

fwrite(schema, file.path(P_OUT, "eicu_sicu_patient_schema.csv"))
fwrite(units, file.path(P_OUT, "eicu_sicu_part0_unittype.csv"))
fwrite(visit, file.path(P_OUT, "eicu_sicu_unitvisitnumber.csv"))
fwrite(locations, file.path(P_OUT, "eicu_sicu_discharge_location.csv"))
fwrite(panel, file.path(P_OUT, "eicu_sicu_part0_panel_nitems.csv"))
cat("Schema:\n"); print(schema)
cat("Unit distribution:\n"); print(units)
cat("unitvisitnumber by unittype:\n"); print(visit)
cat("SICU discharge location:\n"); print(locations)
cat("All-stay panel n_items distribution, offsets 0..1440:\n"); print(panel)
