# =============================================================================
# 13_eicu_extract_panels.R
# Extract SICU first-visit adult stays and their complete five-item, same-offset
# laboratory panels (offset window -1440..1440) through the read-only eICU
# connection. Writes interim/eicu_sicu_aligned_panels.csv.
# Credentials and work_mem come from 00_common; this script holds no password.
# =============================================================================
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
Sys.setenv(RA_ROOT = ra_root)
source(file.path(ra_root, "R", "00_common.R"))

con <- connect_eicu()
on.exit(dbDisconnect(con), add = TRUE)

sql <- paste0(
  "WITH eligible AS MATERIALIZED (",
  " SELECT patientunitstayid,uniquepid,hospitaldischargeyear,",
  " hospitaldischargeoffset,hospitaldischargelocation,",
  " CASE WHEN age ~ '^[0-9]+$' THEN age::int WHEN age='> 89' THEN 90 END AS age ",
  " FROM eicu_crd.patient WHERE unittype='SICU' AND unitvisitnumber=1 ",
  " AND CASE WHEN age ~ '^[0-9]+$' THEN age::int>=18 WHEN age='> 89' THEN TRUE ELSE FALSE END),",
  "panels AS MATERIALIZED (",
  " SELECT l.patientunitstayid,l.labresultoffset,COUNT(DISTINCT l.labname) AS n_items,",
  " COUNT(*) AS n_lab_rows,",
  " MAX(CASE WHEN l.labname='-polys' THEN l.labresult END) AS polys_pct,",
  " MAX(CASE WHEN l.labname='-lymphs' THEN l.labresult END) AS lymphs_pct,",
  " MAX(CASE WHEN l.labname='-monos' THEN l.labresult END) AS monos_pct,",
  " MAX(CASE WHEN l.labname='WBC x 1000' THEN l.labresult END) AS wbc_k,",
  " MAX(CASE WHEN l.labname='platelets x 1000' THEN l.labresult END) AS platelet_k ",
  " FROM eicu_crd.lab l JOIN eligible e USING(patientunitstayid) ",
  " WHERE l.labresultoffset BETWEEN -1440 AND 1440 ",
  " AND l.labname IN ('-polys','-lymphs','-monos','WBC x 1000','platelets x 1000') ",
  " GROUP BY l.patientunitstayid,l.labresultoffset ",
  " HAVING COUNT(DISTINCT l.labname)=5) ",
  "SELECT e.patientunitstayid,e.uniquepid,e.hospitaldischargeyear,e.age,",
  "e.hospitaldischargeoffset,e.hospitaldischargelocation,",
  "p.labresultoffset,p.n_items,p.n_lab_rows,p.polys_pct,p.lymphs_pct,",
  "p.monos_pct,p.wbc_k,p.platelet_k ",
  "FROM panels p JOIN eligible e USING(patientunitstayid) ",
  "ORDER BY e.uniquepid,e.hospitaldischargeyear NULLS LAST,e.patientunitstayid,p.labresultoffset"
)
panels <- as.data.table(dbGetQuery(con, sql))
fwrite(panels, file.path(P_INTERIM, "eicu_sicu_aligned_panels.csv"), na = "")
cat("Aligned panel rows =", nrow(panels),
    "; unique stays =", uniqueN(panels$patientunitstayid),
    "; unique patients =", uniqueN(panels$uniquepid),
    "; missing numeric components =",
    sum(!complete.cases(panels[,.(polys_pct,lymphs_pct,monos_pct,wbc_k,platelet_k)])), "\n")
