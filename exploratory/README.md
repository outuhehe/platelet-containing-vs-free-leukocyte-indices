# Exploratory eICU-CRD diagnostics (pre-finalization audit trail)

This directory preserves the **exploratory and diagnostic scripts** used while
defining the eICU-CRD analysis. They document how key design decisions were
reached from the data itself. They are **not part of `run_all.ps1`** and are not
needed to reproduce the reported results; the authoritative, reproducible eICU
pipeline is `R/13_eicu_extract_panels.R` and `R/14_eicu_analysis.R`.

## Why these are kept separate

The central decision evidenced here is that **eICU-CRD has no procedure-coded
definition of "major abdominal surgery"**. The surgery-identification gate
(§02) showed that only indirect APACHE admission-diagnosis labels and unit type
are available, so the eICU work is reported as **directional consistency**, not
external validation. Keeping the probe scripts makes that conclusion auditable
rather than asserted.

## Contents

### 01_schema_discovery — what the database actually contains

| Script | Purpose |
|---|---|
| `eicu_labname_probe.R` | Quick probe of CBC-related `lab` names (neutrophil/lymphocyte/monocyte/platelet/WBC) and read-only session check. |
| `eicu_labname_census.R` | Full census of CBC-related labnames, pattern coverage, negative-prefix names and unit-type distribution; writes the schema/coverage tables. |
| `eicu_sicu_part0.R` | SICU patient fields (`unittype`, `unitvisitnumber`, discharge location), and the distribution of panel item counts within 0–1440 min. |

### 02_surgery_gate — can abdominal surgery be identified?

| Script | Purpose |
|---|---|
| `eicu_abdominal_surgery_gate.R` | Tests `admissiondx` paths (Operative / Elective / from-O.R. / GI tags) as surgery proxies; establishes that no procedure-confirmed abdominal-surgery definition exists. |
| `eicu_sicu_counts_once.R` | Single lab-window scan of the SICU first-adult cohort: analyzable counts, outcome availability, deaths and five-item panel distribution. |
| `eicu_sicu_patient_panel_diagnostic.R` | Patient-vs-panel diagnostic quantifying repeat stays and the number of unique people, supporting patient-level collapsing. |

### 03_analysis_precursors — drafts that became stages 13/14

| Script | Purpose |
|---|---|
| `extract_eicu_sicu_aligned_panels.R` | Precursor of stage 13: extracts aligned five-item panels with main/pre-window counts. |
| `analyze_eicu_sicu_external.R` | Precursor of stage 14: value QC, patient-level index selection, single AUCs, paired DeLong, 2,000-replicate bootstrap and cross-database comparisons. |
| `export_eicu_roc_coords_raw.R` | Exports raw ROC coordinates per indicator for the ROC figure. |
| `run_eicu_prewindow_paired_bootstrap.R` | Paired DeLong + bootstrap for the alternative pre-admission window (final 1,388 patients / 131 events). |

## Running these (optional, audit only)

- Each script first sources `R/00_common.R` and opens its connection through the
  shared `connect_eicu()` factory. The factory resolves credentials from the
  local `pgpass.conf` with multi-line host/port/database/user matching, so no
  password is hard-coded, and it asserts the session is
  `default_transaction_read_only = on` before any query runs.
- Paths are resolved from the script's own location through the `RA_ROOT` and
  `P_OUT` constants in `00_common.R`, so a script runs directly with
  `Rscript <path-to-script>` from any working directory; it does not rely on a
  relative `outputs/` path.
- Their CSV outputs (which can carry `patientunitstayid` / `uniquepid`) are
  intentionally **not committed here**; only the scripts are version-controlled.
- The frozen, package-consistent implementations used by the main pipeline are
  the numbered stages 13/14, built on the same connection factory and path
  constants.

## Release-safety correction

All exploratory scripts now write to `outputs/interim/exploratory/`, including
aggregate diagnostics and ROC coordinates containing observed thresholds. This
whole directory is restricted and ignored by Git. The raw eICU figure input is
read from the local primary-pipeline row-level export under `outputs/tables/`;
that named file is also ignored. Shared primary aggregate tables are read from
`outputs/tables/`. Run the exploratory extraction before its precursor analysis.
No exploratory output is included automatically in a public archive.
