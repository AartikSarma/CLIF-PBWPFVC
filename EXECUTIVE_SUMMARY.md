# CLIF-PBWPFVC — Executive Summary

## What the project asks

Whether **predicted forced vital capacity (PFVC)** is a better scaling factor than
**predicted body weight (PBW)** for setting tidal volume in mechanically
ventilated ICU adults. PBW (derived from height and sex only) ignores age- and
race-related differences in lung size, so it can systematically misjudge lung
size in older, female, and non-white patients. This pipeline replicates and
extends that analysis on any site's CLIF 2.1 data and emits only aggregated,
cell-suppressed results (n >= 10) designed to be pooled across the consortium.

## What the code does

The pipeline is five R scripts run in order by `code/00_run_pipeline.R`, which
first restores the `renv` environment, then runs 01 -> 05 each as a clean
subprocess.

- **01 — Cohort identification.** Loads the eight required CLIF tables; selects
  the index hospitalization for adults who received invasive ventilation with a
  set tidal volume in the ICU; extracts height (with imputation), weight, SpO2,
  MAP, the SOFA labs, vasopressors, and GCS; runs the respiratory-support
  waterfall; starts the attrition log.
- **02 — Quality checks.** Applies the `outlier-thresholds/` limits to labs,
  vitals, heights, and the waterfall, writing cleaned tables plus per-variable
  summary statistics.
- **03 — Variable derivation.** Computes the core quantities — PBW (Devine),
  PFVC (GLI-2012 via `rspiro`), VT/PBW, VT/PFVC, the **PBW/PFVC** ratio, driving
  pressure, compliance, elastance (and `ers_pbw` / `ers_pfvc`), SOFA, SF/PF
  ratios, VFD-28, and survival times — and builds the cross-sectional and
  longitudinal analysis datasets plus the federated PBW:PFVC distribution (fixed
  bins, n >= 10 suppression).
- **04 — Analysis (the main results).** Table 1 (stratified by PBW/PFVC tercile);
  logistic mortality regressions across the exposure specifications; the
  **elastance-normalized mortality models** (`ers_pbw` / `ers_pfvc`, z-scaled and
  VT/PBW-adjusted); linear models for elastance, compliance, VFD-28, and driving
  pressure; an AIC / evidence-ratio comparison; 60-day survival (Cox +
  Kaplan-Meier by PBW/PFVC tercile); demographic-bias models (each metric vs
  demographics); the H1 broad-cohort PFVC-vs-PBW relationship; conditional-bias
  diagnostic plots; and the unified long-format effect-estimate table that feeds
  the cross-cohort step.
- **05 — Cross-cohort aggregation.** Site-agnostic: discovers every site's
  `regression_results_long_*.csv`, stacks them, and renders forest plots (one per
  analysis + a combined PDF), covariate forests for the bias/H1 analyses, and
  pooled consort and PBW:PFVC distribution figures.

## What the output files are

Everything lands under `output/<site_name>_output/` (and a shared
`output/cross_cohort/`). The whole `output/` tree is gitignored — only
aggregated, n >= 10-suppressed results are written.

### `intermediate/` (scripts 01–03) — working data, not for sharing

- Cohort tables: `cohort_demographics`, `cohort_vitals(_clean)`,
  `cohort_labs(_clean)`, `cohort_meds`, `cohort_assessments`,
  `cohort_heights(_clean)`, `cohort_weights`, `resp_support_waterfall(_clean)`
  (Parquet); `cohort_hospitalization_ids.rds`; `summary_stats/` (lab & vital
  summaries); `attrition_log_partial.csv`
- Analysis datasets: `analysis_cross_sectional`, `analysis_all_timepoints`,
  `analysis_broad_pfvc`, `sofa_scores`, `sofa_daily` (Parquet)

### `final/` (scripts 03–04) — per-site results

- **Federated, meant to be pooled:** `regression_results_long_<site>.csv/.parquet`
  (every effect estimate — the input to script 05), `aic_comparison_all_<site>.csv`,
  `dist_histograms_<site>.csv` & `dist_quantiles_<site>.csv` (PBW:PFVC
  distribution), `attrition_log_<site>.csv`
- **Tables:** `table1_<site>`, `regression_mortality_<site>`,
  `regression_ers_mortality_<site>`, `regression_{ers,crs,vfd_28,dp}_<site>`,
  `table_demographic_bias_<site>`, `table_h1_pfvc_vs_pbw_<site>` (HTML/PDF)
- **Survival:** `km_curves_<site>.pdf`, `cox_model_summary_<site>.txt`
- **Figures:** `evidence_ratio_heatmap_all_<site>.pdf`,
  `distribution_pbwpfvc_<site>.pdf`, `consort_diagram_<site>.pdf`, and the
  conditional-bias plots `bias_pbw_vs_pfvc / bias_mortality / bias_elastance /
  bias_compliance / bias_vfd28 / bias_pbwpfvc_ratio_<site>.pdf`

### `output/cross_cohort/` (script 05) — pooled across sites

- `regression_results_all_cohorts.csv/.parquet` (all sites stacked)
- `forest_plots/forest_<analysis>.pdf` (one per analysis) +
  `forest_all_analyses.pdf`, plus `covforest_*` covariate forests for the
  demographic-bias and H1 analyses
- `consort_diagram_pooled.pdf`, `distribution_pbwpfvc_pooled.pdf`

---

The key cross-site artifact is **`regression_results_long_<site>.csv`** — the
file each site returns and that script 05 consumes to build the pooled forest
plots.
