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
  "WITH sicu_first AS MATERIALIZED (",
  " SELECT patientunitstayid,hospitaldischargeoffset,hospitaldischargelocation ",
  " FROM eicu_crd.patient WHERE unittype='SICU' AND unitvisitnumber=1 ",
  " AND CASE WHEN age ~ '^[0-9]+$' THEN age::int>=18 ",
  "          WHEN age='> 89' THEN TRUE ELSE FALSE END),",
  "lab_window AS MATERIALIZED (",
  " SELECT l.patientunitstayid,l.labresultoffset,l.labname ",
  " FROM eicu_crd.lab l JOIN sicu_first s USING(patientunitstayid) ",
  " WHERE l.labresultoffset BETWEEN -1440 AND 1440 ",
  " AND l.labname IN ('-polys','-lymphs','-monos','WBC x 1000','platelets x 1000')),",
  "panel AS MATERIALIZED (",
  " SELECT patientunitstayid,labresultoffset,COUNT(DISTINCT labname) AS n_items ",
  " FROM lab_window GROUP BY patientunitstayid,labresultoffset),",
  "main_ids AS MATERIALIZED (SELECT DISTINCT patientunitstayid FROM panel ",
  " WHERE labresultoffset BETWEEN 0 AND 1440 AND n_items=5),",
  "pre_ids AS MATERIALIZED (SELECT DISTINCT patientunitstayid FROM panel ",
  " WHERE labresultoffset BETWEEN -1440 AND 0 AND n_items=5),",
  "summary AS (SELECT ",
  " (SELECT COUNT(*) FROM sicu_first) AS n_sicu_first_adult,",
  " (SELECT COUNT(*) FROM main_ids) AS n_analyzable,",
  " (SELECT COUNT(*) FROM main_ids m JOIN sicu_first s USING(patientunitstayid) ",
  "  WHERE s.hospitaldischargeoffset>0) AS n_with_outcome,",
  " (SELECT COUNT(*) FROM main_ids m JOIN sicu_first s USING(patientunitstayid) ",
  "  WHERE s.hospitaldischargeoffset>0 AND lower(trim(s.hospitaldischargelocation))='death') AS n_death,",
  " (SELECT COUNT(*) FROM main_ids m JOIN sicu_first s USING(patientunitstayid) ",
  "  WHERE s.hospitaldischargeoffset>0 AND s.hospitaldischargeoffset<100000) AS n_outcome_sane,",
  " (SELECT COUNT(*) FROM pre_ids) AS n_analyzable_pre_window,",
  " (SELECT COUNT(*) FROM sicu_first s WHERE s.hospitaldischargeoffset>0) AS n_sicu_outcome,",
  " (SELECT COUNT(*) FROM sicu_first s WHERE s.hospitaldischargeoffset>0 ",
  "  AND lower(trim(s.hospitaldischargelocation))='death') AS n_sicu_death),",
  "panel_dist AS (SELECT n_items,COUNT(*) AS n_panels,",
  " COUNT(DISTINCT patientunitstayid) AS n_stay FROM panel ",
  " WHERE labresultoffset BETWEEN 0 AND 1440 GROUP BY n_items) ",
  "SELECT s.*,d.n_items,d.n_panels,d.n_stay FROM summary s CROSS JOIN panel_dist d ",
  "ORDER BY d.n_items"
)

res <- as.data.table(dbGetQuery(con, sql))
summary <- unique(res[, .(n_sicu_first_adult,n_analyzable,n_with_outcome,n_death,
                          n_outcome_sane,n_analyzable_pre_window,n_sicu_outcome,n_sicu_death)])
panel <- res[, .(n_items,n_panels,n_stay)]
if (nrow(summary) != 1L) stop("Expected one summary row")
summary[, `:=`(
  death_rate_by_location = n_death / n_with_outcome,
  n_panels_total = sum(panel$n_panels),
  n_five_item_panels = panel[n_items == 5L, n_panels],
  pct_five_item_panels = panel[n_items == 5L, n_panels] / sum(panel$n_panels)
)]
fwrite(summary, file.path(P_OUT, "eicu_sicu_main_counts.csv"))
fwrite(panel, file.path(P_OUT, "eicu_sicu_target_panel_distribution.csv"))
cat("Single lab-window scan summary:\n"); print(summary)
cat("SICU first adult panel distribution, main window [0,1440]:\n"); print(panel)
