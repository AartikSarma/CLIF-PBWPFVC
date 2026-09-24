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

Figure 4 also reads three optional tables. Without them the pipeline still runs, and
the logs and panel summary say what was not available:

- **crrt_therapy**: `hospitalization_id`, `recorded_dttm` (the start of continuous renal
  replacement ends the creatinine trajectory)
- **patient_procedures**: `hospitalization_id`, `procedure_billed_dttm`, `procedure_code`
  (intermittent dialysis, CPT 90935/90937/90945/90947 or ICD-10-PCS 5A1D70Z/80Z/90Z,
  also ends it)
- **hospital_diagnosis**: `hospitalization_id`, `diagnosis_code` (ESRD, ICD-10
  N18.5/N18.6/Z99.2: no creatinine trajectory at all)

See the [CLIF data dictionary](https://clif-icu.com/data-dictionary) for
guidance on constructing these tables.

## Cohort identification

Adult (age >= 18) ICU encounters receiving invasive mechanical ventilation, with
the height and ventilator data needed to compute PBW, PFVC, and delivered tidal
volume. Detailed inclusion/exclusion criteria and attrition are produced by
`01_cohort_identification.R` and logged to the cohort attrition table. The control
cohort for figure 4 is built the same way: patients who never received advanced
respiratory support (`nosupport`), with its divergence read at the ventilated cohort's
severity. A
noninvasive cohort (`niv`) can be built on request, but it is not a control: NIPPV
delivers large, unlimited positive-pressure volumes.

## Key derived variables

- **PBW** — Devine formula (height, sex).
- **PFVC** — GLI-2012 predicted FVC via `rspiro::pred_GLI()` (height, age, sex,
  race/ethnicity; valid ages 3–95).
- **VT/PBW**, **VT/PFVC**, **PBW/PFVC** ratio.
- Driving pressure, compliance, elastance.
- SOFA from the worst values over the 24 hours from the index, scored by
  [clifR](https://github.com/AartikSarma/clifR)'s `compute_sofa()`; SF/PF ratios.
- VFD-28 (ventilator-free days at 28 days).

Lung-protective ventilation is defined as VT/PBW between 6–8 mL/kg.

## The manuscript, and the scripts behind each figure

| Figure | Claim | Stage | Scripts |
|---|---|---|---|
| 1 | PBW overestimates lung size in older, shorter, female and non-white patients (replicating PMC12313249 in modern cohorts) | `cross_sectional` | `04_analysis.R` |
| 2 | That bias tracks otherwise unexplained differences in respiratory mechanics | `cross_sectional` | `04_analysis.R`, `05_normalization_analysis.R` |
| 3 | The bias is associated with mortality | `cross_sectional` | `04_analysis.R`, `05_normalization_analysis.R` |
| 4 | The bias is associated with rising markers of organ injury over time, in ventilated patients and not in the control cohorts | `injury` | `20`–`29` |

A fifth figure, on causal inference, is not in the paper yet; its scripts are in the
tag `pre-prune-2026-09-19`. [`code/README.md`](code/README.md) lists every script and
what it writes, and says how to restore what was removed.

## Outputs

Each site has one output folder, `output/<site_name>_output/`:

- `intermediate/` holds patient-level data and never leaves the site. The control
  cohorts' patient-level data sits under `intermediate/controls/<cohort>/`.
- `final/` holds aggregates only and is the folder a site returns. It is sorted by
  manuscript block:

  | Subfolder | Holds | Written by |
  |---|---|---|
  | `final/cross_sectional/` | figures 1-3 | `03`-`05` |
  | `final/injury/` | figure 4 | `21`-`28` |
  | `final/controls/` | the control cohorts, all in one folder; file names carry `<site>_<cohort>` | `01`-`03`, `21`-`26` under `PBWPFVC_COHORT` |

  A script asks `utils/config.R` for its folder with `final_dir_for("<block>")`.

Re-running a stage updates `final/` in place, so a site can return the folder again
after any stage. No patient-level data is written to `final/`. Counts are not masked:
small-cell masking (censoring small cells, or a deterministic scheme agreed with CLIF)
will be a separate step applied to `final/` before a site shares it. Until then, do
not share `final/` outside the study team.

**Output file names are an interface between scripts.** Every pooling and figure
script finds its inputs by file-name prefix (`regression_results_long_`, `norm_`,
`jm_`, `injury_`, ...). A prefix may change until the code is distributed to other
sites, but change it together with its readers, in one commit: the pooling scripts in
`code/pooling/`, `24_biotrauma_figures.R` and `27_control_comparison.R`. Once sites
have returned `final/` folders, a renamed prefix orphans their results.

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

### 2. Install uvr

[uvr](https://github.com/nbafrank/uvr) manages the R packages. `uvr.toml` lists the
packages the scripts load and `uvr.lock` pins them and everything they depend on;
`uvr sync` installs them into `.uvr/library/`, and `uvr run` runs a script against
that library. Install uvr by following its
[instructions](https://github.com/nbafrank/uvr#installation). Run every script with
`uvr run`, never bare `Rscript`: bare `Rscript` does not see the project library.

### 3. Run the pipeline

One entry point installs the locked packages (`uvr sync`) and runs the stages you
ask for, each script as a clean subprocess. From the repository root:

```bash
uvr run code/00_run_pipeline.R                             # prep + cross_sectional (figures 1-3)
uvr run code/00_run_pipeline.R -- --stages injury          # figure 4, with its controls
uvr run code/00_run_pipeline.R -- --stages all
```

The `--` separates uvr's options from the runner's.

| Stage | Runs | Rough cost |
|---|---|---|
| `prep` | `01`–`03`: cohort, quality checks, derived variables | minutes |
| `cross_sectional` | `04`, `05` | minutes |
| `injury` | `29_run_figure4.R`: every analysis behind figure 4 and the figure itself (see below) | a few hours |

If a step fails the runner stops and names it. The default, with no `--stages`, is
what this runner has always done, so existing site instructions still work.

### Figure 4 in one command

`code/29_run_figure4.R` runs everything behind figure 4 for a site and draws it:

- **Markers:** platelets, bilirubin, creatinine, vasopressor dose on pressor days, and
  the oxygen saturation index, each over the first 7 days, adjusted and unadjusted.
  Creatinine ends at renal replacement of any kind, continuous or intermittent, which
  is modelled as a third competing event; patients with ESRD are censored at day 0.
- **Arms:** all ventilated patients, and the negative control, patients with no
  respiratory support. The ventilated cohort by baseline SF class is optional:
  `SF_BANDS="235,315 115,235 0,115"`.
- **Severity:** the control is standardised, not matched. Severity cannot confound a
  PFVC fixed by height, age, sex and race, but it could modify the divergence, so the
  control keeps every patient, its divergence varies with the marker's severity
  anchor, and it is read at the ventilated cohort's mean anchor. The severity x
  divergence term tests whether sicker controls diverge faster.

It builds the control cohort when missing, reuses finished fits, carries on past a
failed step and lists the failures at the end. With the figure it writes the
difference-in-differences (ventilated minus control divergence), and after it the
channel breakdown for the supplement (`CHANNEL_MARKERS`, empty skips it).
By default it runs 23 fits at 2,000
iterations, four at a time (`PAR=4`). An earlier estimate put one 7-day fit at a
7,000-patient site at 15-25 GB of memory, so four at once can need 60-100 GB: lower
`PAR` on a smaller machine. `ITER` and `BURNIN` lengthen the chains.
`uvr run code/29_run_figure4.R -- --dry-run` lists every step.

```bash
caffeinate -i nohup uvr run code/29_run_figure4.R > figure4.out 2>&1 &
```

The figure is `final/injury/biotrauma_fig_main_pfvc_7d_<site>.pdf`.

To run a single script while debugging, run it from the repository root, for
example `uvr run code/04_analysis.R`. Scripts read the previous step's outputs, so
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
subfolder there. Override the root with `PBWPFVC_RESULTS_ROOT`. Each script reads
only the block subfolder it pools (`cross_sectional/` or `injury/`) and never
`controls/`, so a control cohort cannot enter a pool of the ventilated cohort. `pooled_biotrauma.R` is in the repository;
`pooled_estimates.R` is kept local and gitignored.

### Archive

`code/archive/` is gitignored local scratch. Scripts removed from the repository are
recoverable from the tag `pre-prune-2026-09-19`.

## Data safety

- Never commit patient data — only aggregated results belong in `output/`.
- Counts in `final/` are unmasked until the masking step exists; see above.
- `config/config.json` is site-specific and should not be committed with real paths.
