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
`01_cohort_identification.R` and logged to the cohort attrition table.

## Key derived variables

- **PBW** — Devine formula (height, sex).
- **PFVC** — GLI-2012 predicted FVC via `rspiro::pred_GLI()` (height, age, sex,
  race/ethnicity; valid ages 3–95).
- **VT/PBW**, **VT/PFVC**, **PBW/PFVC** ratio.
- Driving pressure, compliance, elastance.
- SOFA (extremal aggregation; Severinghaus imputation of PaO2 from SpO2), SF/PF ratios.
- VFD-28 (ventilator-free days at 28 days).

Lung-protective ventilation is defined as VT/PBW between 6–8 mL/kg.

## Expected results

Outputs are written to `output/<site_name>_output/`:

- `intermediate/` — filtered CLIF tables and the derived analysis datasets (Parquet).
- `final/` — Table 1, regression result tables (`regression_results_long_<site>.csv`),
  survival/Kaplan–Meier figures, conditional-bias diagnostic plots, and the
  federated PBW:PFVC distribution exports.

All exports honor a minimum cell size of n >= 10. No patient-level data leaves
the site; only aggregated results are written.

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

A single entry point restores the project environment and runs every script in
order. From the repository root:

```bash
Rscript code/00_run_pipeline.R
```

This restores packages from `renv.lock`, then runs scripts 01–05 sequentially,
each as a clean subprocess. If any step fails, the runner stops and reports which
script errored.

To run a single step (e.g. while debugging), run it directly from the repo root:

```bash
Rscript code/01_cohort_identification.R   # Filter CLIF tables to the eligible cohort
Rscript code/02_quality_checks.R          # Apply outlier thresholds, QC stats
Rscript code/03_variable_derivation.R     # PBW, PFVC, SOFA, SF/PF, VT metrics
Rscript code/04_analysis.R                # Regressions, survival, bias diagnostics
Rscript code/05_cross_cohort_forest.R     # Cross-cohort forest plots (after 04)
```

Scripts must be run in order — each reads the outputs of the previous step.
Script 05 aggregates the per-cohort regression tables produced by script 04 and
is site-agnostic: it discovers every `regression_results_long_*.csv` on disk and
stacks them into cross-cohort forest plots.

### Exploratory analyses (standalone, not in the runner)

These read the script 03 cross-sectional dataset and refit their own models;
they are not part of the consortium pipeline. Run individually after scripts 01–03.

```bash
Rscript code/explore_elastance_age_interaction.R       # Normalized elastance x age interaction
Rscript code/explore_mechanical_power_normalization.R  # Mechanical power scaled to PBW vs PFVC
Rscript code/explore_elastance_fingerprints.R          # Disentangling recoil vs size-surrogate error
Rscript code/explore_stress_strain_mortality.R         # Predicted-mortality stress-strain surfaces
Rscript code/explore_height_sensitivity.R              # Does VT/PBW+PFVC's fit survive height adjustment?
```

## Data safety

- Never commit patient data — only aggregated results belong in `output/`.
- Minimum cell size of n >= 10 for any reported group.
- `config/config.json` is site-specific and should not be committed with real paths.
