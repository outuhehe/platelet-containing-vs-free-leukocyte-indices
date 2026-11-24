-- ============================================================================
-- cbc_available_case.sql
-- Purpose: respond to reviewer Major Comments 2 & 3 (complete-case selection
--          bias and per-comparison available-case samples).
-- This script produces TWO result sets. Export them separately as:
--   Query A -> outputs/cbc_available_case_panels.csv
--   Query B -> outputs/cohort_roster_2376.csv
-- All upstream CTEs are identical to cbc_export_v2.sql; only the final
-- selection is relaxed (the monocyte requirement is removed for Query A).
-- ============================================================================
WITH proc_all AS MATERIALIZED (
    SELECT subject_id, hadm_id, seq_num, chartdate, icd_code,
           icd_version AS icd_code_version
    FROM mimiciv_hosp.procedures_icd
    WHERE icd_code IS NOT NULL
),
d2_45_46 AS (
    SELECT p.icd_code,
           CASE
             WHEN d.long_title ILIKE '%endoscop%' OR d.long_title ILIKE '%colonoscopy%'
               OR d.long_title ILIKE '%sigmoidoscopy%' OR d.long_title ILIKE '%biopsy%'
               OR d.long_title ILIKE '%diagnostic%' OR d.long_title ILIKE '%destruction%'
               OR d.long_title ILIKE '%insertion%' OR d.long_title ILIKE '%percutaneous%'
               OR d.long_title ILIKE '%stent%' OR d.long_title ILIKE '%incision%'
               OR d.long_title ILIKE '%isolation%' OR d.long_title ILIKE '%exteriorization%'
               OR d.long_title ILIKE '%colostomy%' OR d.long_title ILIKE '%ileostomy%'
               OR d.long_title ILIKE '%enterostomy%' OR d.long_title ILIKE '%fixation%'
               OR d.long_title ILIKE '%manipulation%' OR d.long_title ILIKE '%dilation%'
               OR d.long_title ILIKE '%perfusion%' OR d.long_title ILIKE '%pericolostomy hernia%'
               THEN 'exclude'
             WHEN d.long_title ~* '(repair|closure|revision|suture of laceration)'
               AND d.long_title NOT ILIKE '%pericolostomy hernia%' THEN 'review'
             WHEN d.long_title ~* '(resection|excision|anastomosis|enterectomy|colectomy|proctectomy|hemicolectomy|sigmoidectomy|cecectomy)'
               THEN 'keep'
             ELSE 'exclude'
           END AS tier
    FROM (SELECT DISTINCT icd_code FROM proc_all
          WHERE icd_code_version=9 AND (icd_code LIKE '45%' OR icd_code LIKE '46%')) p
    LEFT JOIN mimiciv_hosp.d_icd_procedures d
      ON d.icd_code=p.icd_code AND d.icd_version=9
),
abs_excl AS MATERIALIZED (
    SELECT DISTINCT p.hadm_id,p.chartdate
    FROM proc_all p
    LEFT JOIN mimiciv_hosp.d_icd_procedures d
      ON d.icd_code=p.icd_code AND d.icd_version=p.icd_code_version
    WHERE
      (p.icd_code_version=9 AND (
        d.long_title ILIKE '%endoscop%' OR d.long_title ILIKE '%colonoscopy%'
        OR d.long_title ILIKE '%sigmoidoscopy%' OR d.long_title ILIKE '%catheter%'
        OR d.long_title ILIKE '%stent%' OR p.icd_code IN ('4311','4432','4632','4686','4687')
        OR d.long_title ILIKE '%biopsy%' OR d.long_title ILIKE '%embolectomy%'
        OR d.long_title ILIKE '%thrombectomy%' OR d.long_title ILIKE '%embolization%'
        OR d.long_title ILIKE '%percutaneous%'
        OR p.icd_code IN ('543','544','5461','5462','5463','5464','5472','5473','5474','5475')
        OR p.icd_code LIKE '53%' OR d.long_title ILIKE '%hernia%'
        OR p.icd_code LIKE '49%' OR p.icd_code LIKE '487%' OR p.icd_code LIKE '488%'
        OR p.icd_code LIKE '489%'
        OR (d.long_title ILIKE '%anus%' AND d.long_title NOT ILIKE '%anastomosis%')
        OR d.long_title ILIKE '%transplant%' OR d.long_title ILIKE '%diagnostic%'
        OR d.long_title ILIKE '%exploration%' OR d.long_title ILIKE '%lysis%'
        OR d.long_title ILIKE '%drainage%' OR d.long_title ILIKE '%insertion%'
        OR d.long_title ILIKE '%removal%'))
      OR (p.icd_code_version=10 AND (
        (LEFT(p.icd_code,2) IN ('0D','0F') AND SUBSTRING(p.icd_code FROM 5 FOR 1) IN ('3','8'))
        OR (LEFT(p.icd_code,2) IN ('0D','0F') AND (RIGHT(p.icd_code,1)='X' OR d.long_title ILIKE '%diagnostic%'))
        OR (LEFT(p.icd_code,2) IN ('0D','0F') AND SUBSTRING(p.icd_code FROM 3 FOR 1) IN ('2','H','P','W') AND d.long_title ~* '(device|catheter|stent|tube)')
        OR (LEFT(p.icd_code,2) IN ('0D','0F') AND SUBSTRING(p.icd_code FROM 3 FOR 1)='J')
        OR (LEFT(p.icd_code,2)='04' AND SUBSTRING(p.icd_code FROM 3 FOR 1)='C')
        OR (LEFT(p.icd_code,3)='0J' AND SUBSTRING(p.icd_code FROM 4 FOR 1)='F')
        OR p.icd_code LIKE '0WQF%'
        OR (LEFT(p.icd_code,2)='0D' AND SUBSTRING(p.icd_code FROM 4 FOR 1)='Q')
        OR p.icd_code LIKE '0FY%' OR d.long_title ILIKE '%transplant%'))
),
abs_surg AS (
    SELECT p.subject_id,p.hadm_id,p.seq_num,p.chartdate,p.icd_code,p.icd_code_version
    FROM proc_all p
    JOIN mimiciv_hosp.d_icd_procedures d
      ON d.icd_code=p.icd_code AND d.icd_version=p.icd_code_version
    LEFT JOIN d2_45_46 w ON w.icd_code=p.icd_code AND p.icd_code_version=9
    WHERE
      (p.icd_code_version=9 AND (
        ((p.icd_code LIKE '45%' OR p.icd_code LIKE '46%') AND w.tier='keep')
        OR ((p.icd_code LIKE '43%' OR p.icd_code LIKE '44%' OR p.icd_code LIKE '48%'
          OR p.icd_code LIKE '50%' OR p.icd_code LIKE '510%' OR p.icd_code LIKE '51%'
          OR p.icd_code LIKE '52%' OR p.icd_code LIKE '47%')
          AND d.long_title ~* '(resection|excision|anastomosis|gastrectomy|bypass|colectomy|hemicolectomy|proctectomy|enterectomy|hepatectomy|pancreatectomy|cholecystectomy)')
        OR p.icd_code LIKE '47%'))
      OR (p.icd_code_version=10 AND LEFT(p.icd_code,2) IN ('0D','0F')
        AND SUBSTRING(p.icd_code FROM 3 FOR 1) IN ('1','B','T')
        AND ((LEFT(p.icd_code,2)='0D' AND SUBSTRING(p.icd_code FROM 4 FOR 1) IN ('6','7','8','9','A','B','C','D','E','F','G','H','K','L','M','N','P','Q'))
          OR (LEFT(p.icd_code,2)='0F' AND SUBSTRING(p.icd_code FROM 4 FOR 1) IN ('0','1','2','4','5','6','7','8','9','B','C','D','G')))
        AND SUBSTRING(p.icd_code FROM 5 FOR 1) NOT IN ('3','8')
        AND RIGHT(p.icd_code,1)<>'X' AND d.long_title !~* '(diagnostic|transplant)')
),
abs_clean AS (
    SELECT s.* FROM abs_surg s
    WHERE NOT EXISTS (SELECT 1 FROM abs_excl e WHERE e.hadm_id=s.hadm_id AND e.chartdate=s.chartdate)
),
abs_last AS MATERIALIZED (
    SELECT subject_id,hadm_id,chartdate,last_icd_code,last_icd_version
    FROM (
      SELECT c.subject_id,c.hadm_id,c.chartdate,c.icd_code AS last_icd_code,
             c.icd_code_version AS last_icd_version,
             ROW_NUMBER() OVER (PARTITION BY c.hadm_id ORDER BY c.chartdate DESC,c.seq_num DESC) AS rn
      FROM abs_clean c
    ) ranked
    WHERE rn=1
),
base AS MATERIALIZED (
    SELECT p.hadm_id,p.subject_id,p.chartdate,a.hospital_expire_flag,
           p.last_icd_code,p.last_icd_version
    FROM abs_last p
    JOIN mimiciv_hosp.admissions a USING(hadm_id)
    JOIN mimiciv_hosp.patients pt ON pt.subject_id=p.subject_id
    WHERE NOT EXISTS (
      SELECT 1 FROM mimiciv_hosp.diagnoses_icd d
      WHERE d.hadm_id=p.hadm_id
        AND ((d.icd_version=9 AND d.icd_code LIKE '571%')
          OR (d.icd_version=10 AND (d.icd_code LIKE 'K70%' OR d.icd_code LIKE 'K74%'
              OR d.icd_code IN ('K766','D731'))))
    )
      AND NOT (p.last_icd_version=9 AND (p.last_icd_code LIKE '47%' OR p.last_icd_code LIKE '510%'))
      AND NOT (p.last_icd_version=10 AND (
        (LEFT(p.last_icd_code,2)='0D' AND SUBSTRING(p.last_icd_code FROM 4 FOR 1)='D')
        OR (LEFT(p.last_icd_code,2)='0F' AND SUBSTRING(p.last_icd_code FROM 4 FOR 1)='4')))
),
main AS MATERIALIZED (
    SELECT b.hadm_id, b.subject_id, b.chartdate,
           b.chartdate - INTERVAL '3 days' AS w3,
           b.hospital_expire_flag, a.admittime, a.dischtime, a.deathtime,
           a.admission_type,
           EXTRACT(YEAR FROM a.admittime)::integer AS adm_year
    FROM base b
    JOIN mimiciv_hosp.admissions a USING (hadm_id)
    WHERE EXISTS (
      SELECT 1
      FROM mimiciv_hosp.labevents plt
      WHERE plt.subject_id = b.subject_id
        AND plt.hadm_id = b.hadm_id
        AND plt.charttime >= b.chartdate - INTERVAL '3 days'
        AND plt.charttime < b.chartdate
        AND plt.itemid = 51265
        AND plt.valuenum > 0
    )
),
lab_rows AS MATERIALIZED (
    SELECT m.hadm_id, m.subject_id, m.chartdate, m.w3,
           m.hospital_expire_flag, m.adm_year,
           l.specimen_id, l.charttime, l.itemid, l.valuenum
    FROM main m
    JOIN mimiciv_hosp.labevents l
      ON l.subject_id = m.subject_id
     AND l.hadm_id = m.hadm_id
     AND l.charttime >= m.w3
     AND l.charttime < m.chartdate
     AND l.itemid IN (52075,53159,51133,52769,52074,
                      51256,51244,51254,51301,51265)
    WHERE l.valuenum IS NOT NULL
),
cbc_by_specimen AS (
    SELECT hadm_id, subject_id, hospital_expire_flag, adm_year,
           specimen_id, MAX(charttime) AS panel_time,
           MAX(CASE WHEN itemid IN (52075,53159) AND valuenum > 0 THEN valuenum END) AS anc,
           MAX(CASE WHEN itemid=51133 AND valuenum > 0 THEN valuenum
                    WHEN itemid=52769 AND valuenum > 0 THEN valuenum/1000.0 END) AS alc,
           MAX(CASE WHEN itemid=52074 AND valuenum >= 0 THEN valuenum END) AS amc,
           MAX(CASE WHEN itemid=51256 AND valuenum BETWEEN 0 AND 100 THEN valuenum END) AS neup,
           MAX(CASE WHEN itemid=51244 AND valuenum BETWEEN 0 AND 100 THEN valuenum END) AS lymp,
           MAX(CASE WHEN itemid=51254 AND valuenum BETWEEN 0 AND 100 THEN valuenum END) AS monop,
           MAX(CASE WHEN itemid=51301 AND valuenum > 0 THEN valuenum END) AS wbc,
           MAX(CASE WHEN itemid=51265 AND valuenum > 0 THEN valuenum END) AS plt
    FROM lab_rows
    WHERE specimen_id IS NOT NULL
    GROUP BY hadm_id, subject_id, hospital_expire_flag, adm_year, specimen_id
),
cbc_by_time AS (
    SELECT hadm_id, subject_id, hospital_expire_flag, adm_year,
           NULL::integer AS specimen_id, charttime AS panel_time,
           MAX(CASE WHEN itemid IN (52075,53159) AND valuenum > 0 THEN valuenum END) AS anc,
           MAX(CASE WHEN itemid=51133 AND valuenum > 0 THEN valuenum
                    WHEN itemid=52769 AND valuenum > 0 THEN valuenum/1000.0 END) AS alc,
           MAX(CASE WHEN itemid=52074 AND valuenum >= 0 THEN valuenum END) AS amc,
           MAX(CASE WHEN itemid=51256 AND valuenum BETWEEN 0 AND 100 THEN valuenum END) AS neup,
           MAX(CASE WHEN itemid=51244 AND valuenum BETWEEN 0 AND 100 THEN valuenum END) AS lymp,
           MAX(CASE WHEN itemid=51254 AND valuenum BETWEEN 0 AND 100 THEN valuenum END) AS monop,
           MAX(CASE WHEN itemid=51301 AND valuenum > 0 THEN valuenum END) AS wbc,
           MAX(CASE WHEN itemid=51265 AND valuenum > 0 THEN valuenum END) AS plt
    FROM lab_rows
    GROUP BY hadm_id, subject_id, hospital_expire_flag, adm_year, charttime
),
panels AS MATERIALIZED (
    SELECT * FROM cbc_by_specimen
    UNION ALL
    SELECT * FROM cbc_by_time
),
resolved AS MATERIALIZED (
    SELECT p.*,
           COALESCE(p.anc, p.wbc*p.neup/100.0) AS anc_f,
           COALESCE(p.alc, p.wbc*p.lymp/100.0) AS alc_f,
           COALESCE(p.amc, p.wbc*p.monop/100.0) AS amc_f,
           (p.anc IS NOT NULL AND p.alc IS NOT NULL AND p.amc IS NOT NULL) AS all_abs
    FROM panels p
),
ranked_avail AS (
    SELECT r.*,
           ROW_NUMBER() OVER (
             PARTITION BY hadm_id
             ORDER BY all_abs DESC, panel_time DESC, specimen_id NULLS LAST
           ) AS core_rank
    FROM resolved r
    -- Relaxed filter: monocyte is NOT required here, so stays missing only the
    -- monocyte are retained for the NLR-vs-SII available-case comparison.
    WHERE r.anc_f > 0 AND r.alc_f > 0 AND r.plt > 0
),
best AS MATERIALIZED (
    SELECT * FROM ranked_avail
),
panel_flags AS MATERIALIZED (
    SELECT hadm_id,
           BOOL_OR(anc_f > 0 AND alc_f > 0 AND plt > 0) AS has_any_panel,
           BOOL_OR(anc_f > 0 AND alc_f > 0 AND amc_f > 0 AND plt > 0)
             AS complete_five_index
    FROM resolved
    GROUP BY hadm_id
)
-- ============================================================================
-- QUERY A: eligible candidate panels (multiple rows per hadm may be returned).
-- Export as: outputs/cbc_available_case_panels.csv
-- ============================================================================
SELECT b.subject_id, b.hadm_id, b.hospital_expire_flag,
       b.anc_f, b.alc_f, b.amc_f, b.plt,
       (b.amc_f IS NOT NULL AND b.amc_f > 0) AS has_amc,
       b.all_abs, b.specimen_id, b.panel_time,
       m.admission_type,
       a.admittime, a.dischtime, a.deathtime,
       p.gender AS sex,
       p.anchor_age + (b.adm_year - p.anchor_year) AS age,
       COALESCE(ch.charlson_comorbidity_index,0) AS charlson,
       bl.last_icd_code, bl.last_icd_version
FROM best b
JOIN main m ON m.hadm_id = b.hadm_id
JOIN mimiciv_hosp.admissions a ON a.hadm_id = b.hadm_id
JOIN mimiciv_hosp.patients p ON p.subject_id = b.subject_id
LEFT JOIN abs_last bl ON bl.hadm_id = b.hadm_id
LEFT JOIN mimiciv_derived.charlson ch ON ch.hadm_id = b.hadm_id
ORDER BY b.hadm_id, b.all_abs DESC, b.panel_time DESC, b.specimen_id NULLS LAST;

-- ============================================================================
-- QUERY B: full eligible roster (n=2,376) with the complete-case flag.
-- Run the same CTE chain, then execute the SELECT below instead of Query A.
-- Export as: outputs/cohort_roster_2376.csv
-- Hint: in psql, comment out Query A above and run this second SELECT against
-- the identical WITH chain; or wrap the WITH chain in a temporary view.
-- ============================================================================
-- SELECT m.subject_id, m.hadm_id, m.hospital_expire_flag, m.admission_type,
--        p.gender AS sex,
--        p.anchor_age + (m.adm_year - p.anchor_year) AS age,
--        COALESCE(ch.charlson_comorbidity_index,0) AS charlson,
--        bl.last_icd_code, bl.last_icd_version,
--        COALESCE(f.has_any_panel,FALSE) AS has_any_panel,
--        COALESCE(f.complete_five_index,FALSE) AS complete_five_index
-- FROM main m
-- JOIN mimiciv_hosp.patients p ON p.subject_id = m.subject_id
-- LEFT JOIN abs_last bl ON bl.hadm_id = m.hadm_id
-- LEFT JOIN mimiciv_derived.charlson ch ON ch.hadm_id = m.hadm_id
-- LEFT JOIN panel_flags f ON f.hadm_id = m.hadm_id
-- ORDER BY m.hadm_id;
