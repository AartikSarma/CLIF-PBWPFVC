# Predicted Forced Vital Capacity vs. Predicted Body Weight for Tidal Volume Dosing

## CLIF VERSION

2.1

## Objective

This project investigates whether **predicted forced vital capacity (PFVC)** is a
better scaling factor than **predicted body weight (PBW)** for tidal volume (VT)
dosing in mechanically ventilated ICU patients.

PBW is derived from height and sex alone (Devine formula) and does not account
for age- or race-related differences in lung size. Prior work
([PMC12313249](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12313249/)) showed
that PBW systematically overestimates lung size in older, female, and non-white
patients, which can translate into relatively higher delivered tidal volumes
(VT/PFVC) for those groups despite apparently lung-protective VT/PBW.

This repository replicates that analysis across the CLIF consortium: it derives
PBW and PFVC for each ventilated encounter, characterizes the PBW:PFVC
relationship across demographic strata, and relates VT scaled by each metric to
clinical outcomes.

For a high-level walkthrough of the pipeline and a catalog of every output file,
see [`EXECUTIVE_SUMMARY.md`](EXECUTIVE_SUMMARY.md).

## Required CLIF tables and fields

CLIF version 2.1. The following tables are required:

1. **patient**: `patient_id`, `race_category`, `ethnicity_category`, `sex_category`
2. **hospitalization**: `patient_id`, `hospitalization_id`, `admission_dttm`, `discharge_dttm`, `age_at_admission`, `discharge_category`
3. **adt**: `hospitalization_id`, `in_dttm`, `out_dttm`, `location_category`
4. **vitals**: `hospitalization_id`, `recorded_dttm`, `vital_category`, `vital_value`
   - `vital_category` = 'height_cm', 'weight_kg', 'spo2', 'map'
5. **labs**: `hospitalization_id`, `lab_result_dttm`, `lab_category`, `lab_value`
   - `lab_category` includes 'po2_arterial', 'pco2_arterial', 'creatinine', 'bilirubin_total', 'platelet_count'
6. **medication_admin_continuous**: `hospitalization_id`, `admin_dttm`, `med_category`, `med_dose`, `med_dose_unit`
   - `med_category` = 'norepinephrine', 'epinephrine', 'dopamine', 'dobutamine', 'phenylephrine', 'vasopressin' (vasopressors for SOFA cardiovascular scoring)
7. **patient_assessments**: `hospitalization_id`, `recorded_dttm`, `assessment_category`, `numerical_value`
   - `assessment_category` = 'gcs_total' (for the SOFA neurologic component)
8. **respiratory_support**: `hospitalization_id`, `recorded_dttm`, `device_category`, `mode_category`, `tracheostomy`, `fio2_set`, `peep_set`, `tidal_volume_set`, `tidal_volume_obs`, `resp_rate_set`, `resp_rate_obs`, `plateau_pressure_obs`, `mean_airway_pressure_obs`

See the [CLIF data dictionary](https://clif-icu.com/data-dictionary) for
guidance on constructing these tables.

## Cohort identification

Adult (age >= 18) ICU encounters receiving invasive mechanical ventilation, with
the height and ventilator data needed to compute PBW, PFVC, and delivered tidal
volume. Detailed inclusion/exclusion criteria and attrition are produced by
`01_cohort_identification.R` and logged to the cohort attrition table. Two control
cohorts are built the same way for figure 4: patients who never received advanced
respiratory support (`nosupport`) and patients whose first advanced support was
high-flow oxygen or noninvasive ventilation (`niv`).

## Key derived variables

- **PBW** — Devine formula (height, sex).
- **PFVC** — GLI-2012 predicted FVC via `rspiro::pred_GLI()` (height, age, sex,
  race/ethnicity; valid ages 3–95).
- **VT/PBW**, **VT/PFVC**, **PBW/PFVC** ratio.
- Driving pressure, compliance, elastance.
- SOFA (extremal aggregation; Severinghaus imputation of PaO2 from SpO2), SF/PF ratios.
- VFD-28 (ventilator-free days at 28 days).

Lung-protective ventilation is defined as VT/PBW between 6–8 mL/kg.

## The manuscript, and the scripts behind each figure

| Figure | Claim | Stage | Scripts |
|---|---|---|---|
| 1 | PBW overestimates lung size in older, shorter, female and non-white patients (replicating PMC12313249 in modern cohorts) | `cross_sectional` | `04_analysis.R` |
| 2 | That bias tracks otherwise unexplained differences in respiratory mechanics | `cross_sectional` | `04_analysis.R`, `05_normalization_analysis.R` |
| 3 | The bias is associated with mortality | `cross_sectional` | `04_analysis.R`, `05_normalization_analysis.R` |
| 4 | The bias is associated with rising markers of organ injury over time, in ventilated patients and not in the control cohorts | `injury`, `controls` | `20`–`29` |
| 5 | Causal inference: target trial emulation and a practice-preference instrument | `causal` | `30`–`38` |

[`code/README.md`](code/README.md) lists every script, what it reads and writes, and
the old script names.

## Outputs

Each site has one output folder, `output/<site_name>_output/`:

- `intermediate/` holds patient-level data and never leaves the site. The control
  cohorts' patient-level data sits under `intermediate/controls/<cohort>/`.
- `final/` holds aggregates only and is the folder a site returns. The control
  cohorts' aggregates all sit in `final/controls/`, and their file names carry
  `<site>_<cohort>`.

Re-running a stage updates `final/` in place, so a site can return the folder again
after any stage. All exports honor a minimum cell size of n >= 10. No patient-level
data leaves the site.

**Output file names are an interface between scripts.** Every pooling and figure
script finds its inputs by file-name prefix (`regression_results_long_`, `norm_`,
`jm_`, `injury_`, `tte_`, ...). A prefix may change until the code is distributed to
other sites, but change it together with its readers, in one commit: the pooling
scripts in `code/pooling/`, `24_biotrauma_figures.R`, `27_control_comparison.R`,
`28_biotrauma_summary.R` and `figures/*.py`. Once sites have returned `final/`
folders, a renamed prefix orphans their results.

## Running the project

### 1. Configure the site

Edit `config/config.json` with your site name, the path to your CLIF tables, and
the file format. See [config/README.md](config/README.md).

```json
{
  "site_name": "YOUR_SITE",
  "tables_path": "~/path/to/clif_tables",
  "file_type": "parquet"
}
```

### 2. Run the pipeline

One entry point restores the environment from `renv.lock` and runs the stages you
ask for, each script as a clean subprocess. From the repository root:

```bash
Rscript code/00_run_pipeline.R                          # prep + cross_sectional (figures 1-3)
Rscript code/00_run_pipeline.R --stages injury,controls # figure 4 and its control cohorts
Rscript code/00_run_pipeline.R --stages causal          # figure 5
Rscript code/00_run_pipeline.R --stages all
```

| Stage | Runs | Rough cost |
|---|---|---|
| `prep` | `01`–`03`: cohort, quality checks, derived variables | minutes |
| `cross_sectional` | `04`, `05` | minutes |
| `injury` | `29_run_biotrauma.sh`: panels, fixed-horizon comparators, joint models, figures | hours |
| `controls` | `29_run_controls.sh build` then `anchors`: builds the no-support and noninvasive cohorts together and writes the severity-anchor distributions | under an hour |
| `causal` | `32_tte_run_all.R`, `38_iv_preference.R` | hours |

If a step fails the runner stops and names it. The default, with no `--stages`, is
what this runner has always done, so existing site instructions still work.

The control cohorts need one decision that cannot be automated. The matched
no-support control is restricted to patients above a severity floor, and the floor
is chosen from the ventilated cohort's distribution, which the `controls` stage
writes to `final/jm_severity_anchor_*`. Then:

```bash
SEV_MIN="platelets=2,bilirubin=1" bash code/29_run_controls.sh fits
```

To run a single script while debugging, run it from the repository root, for
example `Rscript code/04_analysis.R`. Scripts read the previous step's outputs, so
they run in number order. Set `PBWPFVC_COHORT=nosupport` (or `niv`) to run a script
on a control cohort.

Scripts 04 and 05 report every exposure-to-outcome estimate both
**demographic-adjusted** (+ age/sex/race) and **unadjusted** (demographics dropped,
illness severity retained).

### Cross-site pooling (run centrally)

Pooling is **not** part of the per-site pipeline. The study coordinator runs the
scripts in `code/pooling/` after every site has returned its `final/` folder. They
expect a results root with one subfolder per site (each site's `final/` renamed to
the site name), by default the local `results/` folder, and write to an `All sites/`
subfolder there. Override the root with `PBWPFVC_RESULTS_ROOT`. They list each site
folder without recursing, so a site's `controls/` subfolder is never pooled with the
ventilated cohort by accident. `pooled_biotrauma.R` is in the repository;
`pooled_estimates.R` and `pooled_tte.R` are kept local and gitignored.

### Supplements and archive

`code/supplement/` holds tracked lead-site analyses that no runner calls:
sensitivity analyses, diagnostics and earlier lines of inquiry. Its
[README](code/supplement/README.md) is the inventory. `code/archive/` is gitignored
local scratch.

## Data safety

- Never commit patient data — only aggregated results belong in `output/`.
- Minimum cell size of n >= 10 for any reported group.
- `config/config.json` is site-specific and should not be committed with real paths.
