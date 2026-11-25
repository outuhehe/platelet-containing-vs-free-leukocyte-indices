# Reproducible analysis — platelet-free vs platelet-containing inflammatory indices after major abdominal surgery (MIMIC-IV / eICU-CRD)

Auditable, fully local, one-command reproduction of the analysis reported in the
manuscript submitted for peer review. The central finding is that, among patients with a complete
preoperative differential-CBC panel, the **platelet-free** indices (NLR, SIRI)
discriminate in-hospital death better than the **platelet-containing** products
(SII, PIV) on paired AUC comparisons. The eICU analysis is a **direction-consistency
analysis**, not an external validation; calibration is evaluated and decision-curve analysis is presented descriptively in current Figure 5 (file Figure5_calibration_dca) and Additional file 6; these analyses do not establish clinical utility.

This archive documents the analyses reported in the manuscript and does not establish
the timing of prospective registration or protocol approval.

## Design principles (following published MIMIC code repositories)

- **Database access is read-only.** The database-enabled scripts in this archive
  (01, 02, 13 and 15) use `00_common.R`, which sets and verifies PostgreSQL
  read-only mode. Analysis stages otherwise read local intermediates and contain
  no credentials.
- **Configuration is externalized.** Connection details live in `config/db.env`
  (local only, git-ignored); a template is provided in `config/db.env.example`.
- **Numbered, single-direction stages.** Each stage reads the previous stage's
  output and writes the next; hard assertions act as frozen-value gates.
- **Single orchestrator.** `run_all.ps1` runs every stage in order and stops at
  the first failure.
- **Row-level data remain local and access-controlled** under the PhysioNet DUA.
  The reproducibility archive contains code and aggregate result tables/figures
  only; it does not contain row-level extraction intermediates.

## Repository layout

```
reproducible_analysis/
├── config/
│   ├── db.env.example          # connection template (copy to db.env and edit)
│   └── db.env                  # local credentials (git-ignored; not for sharing)
├── sql/
│   ├── cbc_available_case.sql  # roster / Query A + Query B logic
│   └── cbc_export_v2.sql       # strict five-index panel extraction
├── R/
│   ├── 00_common.R             # paths, env loading, read-only connection factory, IO
│   ├── 01..15_*.R              # numbered pipeline stages (see table below)
├── figures/
│   └── 21_build_public_figures.R # ROC + aggregate Figures 2-5
├── exploratory/                # pre-finalization eICU probes (audit only; not in run_all)
│   └── README.md               # schema discovery, surgery gate, analysis precursors
├── Rlib/                       # local offline library (git-ignored; survival, survRM2)
│                               # install_pkgs.R installs the same versions from CRAN
├── tests/
│   ├── smoke_common.R          # environment / read-only / package smoke test
│   ├── test_pgpass.R           # pgpass multi-line matching unit checks
│   ├── parse_all.R             # syntax-parse every R script
│   └── make_session_info.R     # regenerate session_info.txt
├── outputs/
│   ├── interim/                # local row-level restricted data; never for release
│   ├── tables/                 # all result CSV tables
│   ├── figures/                # Figure 1 (and any later figures)
│   └── logs/                   # per-stage stdout/stderr
├── install_pkgs.R              # one-time CRAN dependency installer (idempotent)
├── session_info.txt            # R version and exact package versions
├── protocol.md                 # frozen analysis plan, decisions and revision log
├── run_all.ps1                 # single entry point
├── LICENSE                     # MIT License
└── README.md
```

## Requirements

- **R 4.5.1** (or compatible). Packages used: `data.table`, `pROC`, `DBI`,
  `RPostgres`, `survival`, and `survRM2`. Dependencies are **pinned to exact
  versions** and installed by `install_pkgs.R` (via `remotes::install_version`),
  so a later machine does not silently get newer CRAN releases. A local `Rlib/`
  (git-ignored) is used first for offline runs; `00_common.R` adds it to the
  library path when present. Exact installed versions are recorded in
  `session_info.txt` (regenerate with `Rscript tests/make_session_info.R`).
- A local **PostgreSQL** instance holding the `mimic` and `eicu` databases,
  reachable with a **read-only** role. The connection factory enforces
  `default_transaction_read_only=on` and refuses to proceed if the session is not
  read-only; `bigint="numeric"` avoids integer64 scalar-display issues.
- No internet connection is needed once the local databases and R packages are
  in place. No data leave the machine. Network use is confined to the one-time
  dependency install in `install_pkgs.R`; `run_all.ps1` and every numbered stage
  make no outbound request. When `Rlib/` already holds the pinned versions,
  `install_pkgs.R` returns without contacting CRAN.

## Install R packages (first time only)

```powershell
# Installs all six packages at their pinned versions (uses a local Rlib first
# if present, otherwise CRAN). Idempotent; already-pinned versions are skipped.
Rscript install_pkgs.R
```

## Configuration

1. Copy the template: `Copy-Item config/db.env.example config/db.env`
2. Edit `config/db.env` with `PGHOST`, `PGPORT`, `MIMIC_DB`, `EICU_DB`,
   `PGUSER`, `PGPASSWORD`. If `db.env` is absent, `00_common.R` falls back to
   the local `%APPDATA%\postgresql\pgpass.conf`.
3. First-time environment check:
   `Rscript tests\smoke_common.R`

## Run

```powershell
# Full pipeline (published configuration, including all bootstraps)
powershell -ExecutionPolicy Bypass -File .\run_all.ps1

# Fast logic check that skips the IPW bootstrap
powershell -ExecutionPolicy Bypass -File .\run_all.ps1 -SkipBootstrap
```

The orchestrator writes each stage's logs to `outputs/logs/` and stops at the
first non-zero exit code. To run stages manually, execute the numbered scripts in
order with `Rscript`; each script locates the project root from its own path.

## Stages

| # | Script | What it does | Key output |
|---|---|---|---|
| 01 | `extract_roster` | Query A candidate panels + Query B full roster (2,376); ICD-9/ICD-10 codebooks; 886-row index | `interim/cohort_roster_2376.csv`, `cbc_available_case_panels.csv`, `cbc_available_case_index.csv`, `D2_ICD9_45_46_codebook.csv`, `ICD10_actual_target_codes.csv` |
| 02 | `extract_index_panel` | Strict five-index panel via `cbc_export_v2.sql` (879) | `interim/cbc_export_v2.csv`, `cbc_export_with_survival.csv` |
| 03 | `paired_analysis` | Spearman correlations, single AUCs, paired DeLong, bootstrap, RMST | `tables/table2_single_auc.csv`, `table3_paired_delta*.csv`, `table3b_bootstrap_ci_v2.csv` |
| 04 | `valid_paired_comparisons` | Final 879-sample paired DeLong + 2,000 bootstrap, PLR | `tables/table3_paired_delta_v2.csv`, bootstrap CI |
| 05 | `lymphocyte_auc` | Lymphocyte / neutrophil single-AUC positive controls | `tables/lymphocyte_auc_positive_control.csv` |
| 06 | `lymphocyte_vs_nlr` | NLR vs lymphocyte paired DeLong | `tables/lymphocyte_vs_nlr_delong.csv` |
| 07 | `available_case_analysis` | Full-roster (2,376) available-case analysis | `tables/available_case_*.csv` |
| 08 | `review_sensitivity` | Panel-selection (absolute-first vs latest) sensitivity | `tables/review_sensitivity_*.csv` |
| 09 | `reviewer_additional` | PLR, platelet AUC, all-absolute subgroup pairs | `tables/reviewer_additional_*.csv` |
| 10 | `ipw_dual_ps` | Dual propensity models, overlap/IPW weighted AUC, diagnostics, balance, bootstrap | `tables/ipw_diagnostics.csv`, `ipw_balance.csv`, `ipw_weighted_primary.csv`, `ipw_bootstrap_ci.csv` |
| 11 | `complete_submission_tables` | Paper-format Table 2 / Table 3 with cross-checks, all-absolute subgroup | `tables/Table2_single_auc_paper.csv`, `Table3_paired_auc_paper.csv` |
| 12 | `figure1_flow` | Two-panel MIMIC/eICU flow diagram (300 dpi) | `figures/Figure1_flow.png` |
| 13 | `eicu_extract_panels` | SICU first-visit adult complete five-item panels | `interim/eicu_sicu_aligned_panels.csv` |
| 14 | `eicu_analysis` | QC, index derivation, single AUC, paired DeLong + bootstrap, pre-window sensitivity, eICU vs MIMIC | `tables/eicu_*.csv` |

Manuscript figure order and the file names, which now follow the manuscript numbering:

| Manuscript | Content | File stem |
|---|---|---|
| Figure 1 | two-panel cohort flow | `Figure1_flow` |
| Figure 2 | two-panel ROC curves | `Figure2_roc_two_panels` |
| Figure 3 | nine-row robustness forest | `Figure3_robustness_forest` |
| Figure 4 | adjusted platelet RCS | `Figure4_platelet_rcs` |
| Figure 5 | calibration and decision curves | `Figure5_calibration_dca` |

The platelet tertile and IPW love-plot figures are supplementary (Figure S1, S2).

## Frozen-value gates (reproduction must match)

- Roster: **2,376 stays / 181 deaths**; complete five-index panel **879 / 87 deaths**.
- Primary paired differences (879): SII vs NLR AUC **0.581 vs 0.690**, Δ **−0.109**
  (DeLong p ≈ 1.19×10⁻⁶); PIV vs SIRI **0.573 vs 0.661**, Δ **−0.087** (p ≈ 1.15×10⁻⁷).
- Single-marker AUCs (879): NLR 0.690, lymphocyte 0.687, platelet 0.669,
  SIRI 0.661, SII 0.581, PIV 0.573, PLR 0.544.
- Baseline IPW weighted Δ **+0.122 / +0.094**; max absolute SMD **0.019**;
  max stabilized weight **2.76**. Availability-model C-statistic ≈ 0.600 (does
  not prove or refute MAR; weighting rests on overlap and balance).
- eICU: **2,855 patients / 276 deaths**; paired Δ **−0.055 / −0.037**
  (same direction as MIMIC, smaller magnitude).

## Code provenance and audit

The pipeline is designed for a workflow in which a local model may **draft** code,
but every script is written to a file, **human-reviewed**, and then executed by a
local interpreter. The committed scripts are the reviewed, executed versions; the
intermediate CSVs and the gates above let any reader verify the numbers without
re-running the database extraction. The pre-finalization eICU probes that justify
the directional-consistency framing (schema discovery, the surgery-identification
gate, and the drafts that became stages 13/14) are preserved under
[`exploratory/`](exploratory/README.md) for audit; they are not part of `run_all.ps1`.

## Published MIMIC code examples

This archive follows two established patterns in published MIMIC research:

- The [MIMIC Code Repository](https://github.com/MIT-LCP/mimic-code), described
  by Johnson et al. in *The MIMIC Code Repository: enabling reproducibility in
  critical care research* ([JAMIA, 2018](https://doi.org/10.1093/jamia/ocx084)),
  organizes reusable SQL and analysis code by dataset and concept. Its code
  repository does not redistribute credentialed patient-level data.
- The paper-specific [mechanical ventilation liberation analysis repository](https://github.com/amanchoudhri/liberation-timing-paper)
  accompanies Shahn et al., *Effects of aggressive and conservative strategies
  for mechanical ventilation liberation* ([Journal of Critical Care, 2023](https://doi.org/10.1016/j.jcrc.2023.154275)).
  It links the code directly to the published article and keeps the analysis
  materials in a focused, inspectable repository.

Accordingly, this archive links its analysis stages to the manuscript outputs,
records the software environment and frozen-value checks, and excludes all
credentialed row-level data. The paper-specific repository is deliberately
distributed as a supplementary archive here; no public GitHub URL or publication
identifier is claimed for this manuscript.

## Data availability and PhysioNet DUA

Raw MIMIC-IV and eICU-CRD data are accessed under PhysioNet credentialed access
and are **not redistributed** here. Obtain them from [PhysioNet](https://physionet.org)
after completing the required training and data-use agreement.

The extraction stages create row-level files under `outputs/interim/`; these
retain stable PhysioNet pseudo-identifiers (`subject_id`/`hadm_id` for MIMIC-IV;
`patientunitstayid`/`uniquepid` for eICU-CRD) alongside clinical variables.
Treat these files as restricted data governed by the applicable PhysioNet DUA,
not as de-identified public data. Keep them and database credentials in the
access-controlled local research environment. Do not email, attach, commit,
publish, or upload them to cloud storage or external services.

The `.gitignore` excludes `outputs/interim/` and the two row-level eICU exports
(`eicu_sicu_patient_index_roster.csv` and
`eicu_sicu_external_validation_dataset.csv`). This is a guardrail, not a release
approval: inspect every archive and repository change before sharing. Removing
the pseudo-identifier columns alone does not establish that row-level data are
safe to disclose. Only release aggregate results after checking that they do not
reasonably permit identification; seek PhysioNet guidance for any proposed
release of row-level or otherwise potentially identifying derived data.

The submitted reproducibility archive contains analysis code and aggregate
tables/figures only. It does not include `outputs/interim/`, the row-level eICU
exports, database credentials, or query results at the patient/stay level.

## License

Code is released under the MIT License. The manuscript and analysis outputs are
governed by the journal's copyright on publication.

## Extension dependencies and release gate

Additional files 6/7 use `rms` 8.2-0 and `ggplot2` 4.0.3 in addition to the
primary dependencies listed above. The primary installer does not install these
optional extension packages. The current main adjusted spline is generated by
script 19; script 18 only produces a separately named historical unadjusted
curve. All exploratory products remain under `outputs/interim/exploratory/`.
Before distributing a ZIP, run `python tests/check_release.py <archive.zip>`.
This checks restricted paths, credential-file names, identifier/time columns,
and CSV paths/schemas/row bounds against `tests/release_csv_allowlist.json`; a passing scan is not a legal
certification or an approval to release new row-level derived data.


## Exploratory reclassification and winsorization extension

Script `R/21_reclassification_winsorization.R` reads only the restricted local
intermediates regenerated by primary stages 01-15. It does not access the
database and does not write patient-level predictions or records. Aggregate
results are provided in Additional file 7. This extension was added post hoc;
it does not change the original primary analyses or establish clinical utility.

Run after `run_all.ps1` has completed, from this extracted package directory:

```powershell
$env:RA_ROOT = (Get-Location).Path
Rscript .\R\21_reclassification_winsorization.R
```

Default settings: 20 repeated subject-grouped, outcome-stratified 10-fold
partitions, 500 full refit subject-cluster bootstrap samples for NRI/IDI,
2,000 subject-cluster bootstrap samples for paired-AUC winsorization, 4 workers.
Optional environment settings: `STUDY_NRI_OUT`, `STUDY_NRI_B`,
`STUDY_WINSOR_B`, `STUDY_NRI_REPEATS`, `STUDY_NRI_WORKERS`.
Reduced bootstrap settings are for local checks and do not reproduce the
released confidence intervals. No additional R packages are needed.

Prediction-model admission coding is EMER-labelled, ELECTIVE and other (URGENT
included in other); the M3 association models combine emergency and urgent. The
coding distinction is reported explicitly rather than silently changed.

Only continuous NRI is reported because no clinical risk categories were
prespecified. IDI in CSV files is on the absolute probability scale; the
manuscript supplementary table multiplies it by 100 (percentage points).
Component-wise winsorization rebuilds formulas after capping blood counts;
index-wise capping is a separate sensitivity. Both retain all common-sample stays.

## Supplementary material

`supplementary/` mirrors the files supplied with the manuscript.

| File | Content |
| --- | --- |
| `Supplementary_material.pdf` | Supplementary methods M1-M6, Figures S1-S2 and Tables S1-S6 |
| `coi_disclosure_QiYang.pdf` | Signed ICMJE conflict-of-interest disclosure |
| `Additional_file_1.zip` to `Additional_file_7.zip` | Code archive and aggregate result tables, matching the numbering in the manuscript |

These files contain analysis code, method documentation and aggregate results only. No patient-,
stay- or row-level records, identifiers, credentials or individual-level derived data from
MIMIC-IV or eICU-CRD are included. The source databases remain available only to credentialed
users on PhysioNet (https://physionet.org/) under the PhysioNet Credentialed Health Data Use
Agreement; redistribution of the source data is prohibited and nothing in this repository
provides access to it. `Additional_file_1.zip` duplicates the code at the repository root and is
supplied only so that the file numbering matches the manuscript.
