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

sql <- paste0(
  "WITH eligible AS MATERIALIZED (SELECT patientunitstayid,uniquepid,",
  "hospitaldischargeyear,hospitaldischargeoffset,hospitaldischargelocation ",
  "FROM eicu_crd.patient WHERE unittype='SICU' AND unitvisitnumber=1 ",
  "AND CASE WHEN age ~ '^[0-9]+$' THEN age::int>=18 WHEN age='> 89' THEN TRUE ELSE FALSE END),",
  "panels AS MATERIALIZED (SELECT l.patientunitstayid,l.labresultoffset,",
  "COUNT(DISTINCT l.labname) AS n_items FROM eicu_crd.lab l ",
  "JOIN eligible e USING(patientunitstayid) ",
  "WHERE l.labresultoffset BETWEEN 0 AND 1440 AND l.labname IN ",
  "('-polys','-lymphs','-monos','WBC x 1000','platelets x 1000') ",
  "GROUP BY l.patientunitstayid,l.labresultoffset),",
  "complete_stays AS MATERIALIZED (SELECT p.patientunitstayid FROM panels p WHERE n_items=5 GROUP BY p.patientunitstayid),",
  "rows AS (SELECT f.patientunitstayid,e.uniquepid,e.hospitaldischargeyear,",
  "e.hospitaldischargeoffset,e.hospitaldischargelocation FROM complete_stays f JOIN eligible e USING(patientunitstayid)),",
  "person_ranked AS (SELECT r.*,ROW_NUMBER() OVER(PARTITION BY uniquepid ",
  "ORDER BY hospitaldischargeyear NULLS LAST,patientunitstayid) AS rn FROM rows r),",
  "person_summary AS (SELECT uniquepid,BOOL_OR(lower(trim(hospitaldischargelocation))='death') AS any_death,",
  "BOOL_OR(hospitaldischargeoffset>0) AS any_offset_ok FROM rows GROUP BY uniquepid) ",
  "SELECT (SELECT COUNT(*) FROM panels WHERE n_items=5) AS n_five_item_panels,",
  "(SELECT COUNT(*) FROM complete_stays) AS n_eligible_stays,",
  "(SELECT COUNT(DISTINCT uniquepid) FROM rows) AS n_unique_people,",
  "(SELECT COUNT(*) FROM rows WHERE hospitaldischargeoffset>0) AS n_stays_outcome,",
  "(SELECT COUNT(*) FROM rows WHERE hospitaldischargeoffset>0 AND lower(trim(hospitaldischargelocation))='death') AS n_stay_deaths,",
  "(SELECT COUNT(*) FROM person_ranked WHERE rn=1) AS n_people_first_index,",
  "(SELECT COUNT(*) FROM person_ranked WHERE rn=1 AND hospitaldischargeoffset>0) AS n_people_first_outcome,",
  "(SELECT COUNT(*) FROM person_ranked WHERE rn=1 AND lower(trim(hospitaldischargelocation))='death') AS n_people_first_deaths,",
  "(SELECT COUNT(*) FROM person_summary WHERE any_death) AS n_people_any_death,",
  "(SELECT COUNT(*) FROM person_summary WHERE any_offset_ok) AS n_people_any_outcome_ok,",
  "(SELECT COUNT(*) FROM rows) - (SELECT COUNT(DISTINCT uniquepid) FROM rows) AS excess_repeat_stays"
)
res <- as.data.table(dbGetQuery(con, sql))
fwrite(res, file.path(P_OUT, "eicu_sicu_patient_panel_diagnostic.csv"))
cat("Patient vs panel diagnostic:\n"); print(res)
