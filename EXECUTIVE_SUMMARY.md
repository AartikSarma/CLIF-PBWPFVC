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
  VT/PBW-adjusted); linear models for elastance, compliance, and driving pressure;
  **28-day VFDs modeled as a competing-risks outcome** (extubation vs death, via
  Fine–Gray subdistribution models, per Yehya & Harhay 2019) rather than a
  continuous value; an AIC / evidence-ratio comparison; 60-day survival (Cox +
  Kaplan-Meier by PBW/PFVC tercile); demographic-bias models (each metric vs
  demographics); the predicted-FVC-vs-predicted-body-weight model (regressing
  PFVC on PBW plus age, sex, and race to show that PBW alone does not capture
  predicted lung size); conditional-bias diagnostic plots; and the unified
  long-format effect-estimate table that feeds the cross-cohort step.
- **05 — Cross-cohort aggregation.** Site-agnostic: discovers every site's
  `regression_results_long_*.csv`, stacks them, and renders forest plots (one per
  analysis + a combined PDF), covariate forests for the demographic-bias and
  PFVC-vs-PBW analyses, and pooled consort and PBW:PFVC distribution figures.

### Exploratory analyses (standalone — not in the pipeline runner)

These read the script 03 cross-sectional dataset and refit their own models. They
are run individually (not via `00_run_pipeline.R`) and are not consortium
deliverables. Mechanical power (`mp_pbw`/`mp_pfvc`) is still derived in script 03
with the mode-aware simplified power equation; the modeling lives here.

- **`explore_elastance_age_interaction.R`** — refits the mortality (and 28-day VFD
  competing-risks) models with an elastance × age interaction (Ers×PBW and Ers×PFVC,
  per SD), since recoil changes with age independently of injury; LRT + probability-
  scale slopes via `marginaleffects`.
- **`explore_mechanical_power_normalization.R`** — compares mechanical power scaled
  to PBW (J/min/kg) vs PFVC (J/min/L), z-scaled, as predictors of mortality and
  28-day VFDs, with an AIC evidence ratio and an MP × age interaction.
- **`explore_elastance_fingerprints.R`** — probes the statistical fingerprints that
  distinguish the size-surrogate error from non-constant recoil (size vs injury
  axes, predicted FEV1/FVC as a recoil surrogate, variance partition).
- **`explore_stress_strain_mortality.R`** — predicted-mortality surfaces over the
  strain (VT/PFVC) × stress/elastance plane, with a model-fit comparison.

## What the output files are

Everything lands under `output/<site_name>_output/` (and a shared
`output/cross_cohort/`). The whole `output/` tree is gitignored — only
aggregated, n >= 10-suppressed results are written.

### `intermediate/` (scripts 01–03) — working data, not for sharing

These are the building blocks consumed by later scripts; they contain row-level
data and should **not** leave the site.

- **`cohort_hospitalization_ids.rds`** — an R vector of the eligible index
  `hospitalization_id`s that define the analytic cohort; every downstream script
  uses it as the cohort filter.
- **`cohort_demographics.parquet`** — one row per index hospitalization:
  `patient_id`, `hospitalization_id`, age, sex, harmonized `race_category`
  (WHITE / BLACK / OTHER), ethnicity, admission/discharge times,
  `discharge_category`, and `deceased` (in-hospital mortality flag).
- **`cohort_vitals.parquet`** / **`cohort_vitals_clean.parquet`** — long format,
  one row per measurement, for SpO2 and MAP (`hospitalization_id`,
  `recorded_dttm`, `vital_category`, `vital_value`). The `_clean` version has
  out-of-range values set to `NA` per `outlier-thresholds/`.
- **`cohort_labs.parquet`** / **`cohort_labs_clean.parquet`** — long-format lab
  values used for SOFA and oxygenation (PaO2, PaCO2, creatinine,
  bilirubin_total, platelets); `_clean` is outlier-filtered.
- **`cohort_meds.parquet`** — continuous vasopressor/inotrope infusions
  (norepinephrine, epinephrine, dopamine, dobutamine, phenylephrine,
  vasopressin) for the SOFA cardiovascular component.
- **`cohort_assessments.parquet`** — GCS-total assessments for the SOFA
  neurologic component.
- **`cohort_heights.parquet`** / **`_clean`** — one row per hospitalization with
  imputed `height_cm` (measured height where available, otherwise the patient's
  median across admissions); drives PBW and PFVC.
- **`cohort_weights.parquet`** — one row per hospitalization with mean recorded
  weight (sanity-bounded), used to convert vasopressor doses to mcg/kg/min.
- **`resp_support_waterfall.parquet`** / **`_clean`** — the processed
  respiratory-support time series after the waterfall fill (device, mode, set
  parameters; plateau and mean airway pressure are **not** forward-filled);
  `_clean` is outlier-filtered.
- **`summary_stats/lab_summary_<site>.csv`**, **`vital_summary_<site>.csv`** —
  per-variable QC counts (n, missing, outliers removed, observed ranges).
- **`attrition_log_partial.csv`** — the first three cohort-funnel steps
  (adults → ICU → invasive ventilation with a set tidal volume), as distinct
  patient counts with exclusion reasons; script 03 appends the remaining steps.
- **`analysis_cross_sectional.parquet`** — the **primary analytic table**: one
  row per index hospitalization carrying every derived variable — PBW, PFVC,
  VT/PBW, VT/PFVC, PBW/PFVC, driving pressure, compliance, elastance,
  `ers_pbw`/`ers_pfvc`, mechanical power and its scalings (`mp_pbw`/`mp_pfvc`),
  SOFA, SF/PF ratios, VFD-28 (`vfd_28`) plus its competing-risks form
  (`vfd_time`/`vfd_status`: extubation vs death), BMI, survival time and event,
  mortality flags, and demographics. The index timepoint is the first lung-protective (VT/PBW 6–8),
  hypoxemic (SF < 315) timepoint with recorded pressures (so driving pressure /
  elastance are observed) within the first 6 h of ventilation, falling back to the
  first qualifying timepoint when no such early pressures exist.
- **`analysis_all_timepoints.parquet`** — the longitudinal counterpart: one row
  per hospitalization × ventilator timepoint, for time-varying analyses.
- **`analysis_broad_pfvc.parquet`** — a broader cohort (all eligible subjects with
  height, age, sex, race, and a computable PFVC — wider than the ventilated
  analytic cohort) used only for the predicted-FVC-vs-predicted-body-weight
  model.
- **`sofa_scores.parquet`** — per-encounter aggregated SOFA (worst component
  values and total). **`sofa_daily.parquet`** — the per-encounter-day components
  and totals.

### `final/` (scripts 03–04) — per-site results

These are the deliverables. The first group is aggregated and cell-suppressed
(n >= 10) and is what each site returns to the consortium; the rest are
human-readable tables and figures for local review.

**Federated tables (return these):**

- **`regression_results_long_<site>.csv` / `.parquet`** — the key cross-site
  artifact. One row per model term across **every** model in script 04, with
  columns: `site`, `term`, `estimate`, `conf_low`, `conf_high`, `std_error`,
  `statistic`, `p_value`, `estimate_type` (OR / HR / Beta), `analysis`,
  `model_spec` (the exposure specification), `model_family`, `formula`, and
  `n_obs`. Script 05 stacks these across sites to build the forest plots.
- **`aic_comparison_all_<site>.csv`** — model-fit comparison: per
  `outcome` × `exposure`, the `AIC`, `delta_AIC`, and `evidence_ratio` relative
  to the PBW-scaled reference within that outcome (VT/PBW, or Ers×PBW for the
  elastance-normalized mortality models), plus `is_reference` and the truncated
  ratio / display label.
- **`dist_histograms_<site>.csv`** — the federated PBW:PFVC distribution as
  histogram counts: one row per `group_type` (sex / race / age_bin / height_bin)
  × `group_value` × bin, with `bin_left`, `bin_right`, `count`. Fixed 0–50 bins
  (width 1) with tails clamped, so per-site histograms sum directly into a pooled
  distribution; groups with n < 10 are dropped.
- **`dist_quantiles_<site>.csv`** — the same groupings summarized as quantiles
  (`n`, `median`, `q1`, `q3`, and `p10`…`p90`) for exact per-site boxplots.
- **`attrition_log_<site>.csv`** — the complete CONSORT funnel (all steps): step
  label, patients remaining, patients excluded, and exclusion reason.

**Regression tables (HTML/PDF, for review):**

- **`table1_<site>.html` / `.pdf`** — cohort characteristics stratified by
  PBW/PFVC tercile (T1/T2/T3 plus an overall column), medians (IQR) for
  continuous variables and n (%) for categoricals, with group-comparison
  p-values.
- **`regression_mortality_<site>.html`** — merged logistic regression (odds
  ratios) for in-hospital mortality across the five exposure specifications
  (VT/PFVC; VT/PBW; VT/PFVC + VT/PBW; VT/PBW + PFVC; VT/PBW + PBW/PFVC).
- **`regression_ers_mortality_<site>.html`** — the elastance-normalized mortality
  models: in-hospital mortality on z-scaled Ers×PBW and Ers×PFVC (odds ratio per
  1 SD), adjusted for VT/PBW and the standard covariates.
- **`regression_ers_<site>.html`**, **`regression_crs_<site>.html`**,
  **`regression_dp_<site>.html`** — linear regressions of elastance, compliance,
  and static driving pressure, each across the same five exposure specifications.
- **`regression_vfd28_<site>.html`** — 28-day VFDs as a **competing-risks**
  outcome (Yehya & Harhay 2019): Fine–Gray subdistribution models for extubation
  (death as the competing risk), across the same five exposure specifications,
  reporting subdistribution hazard ratios (SHR > 1 = faster liberation). In the
  federated long table this analysis carries `estimate_type = HR`. Mortality (the
  other component of the composite) is reported separately above; VFD-28 is also
  shown descriptively in Table 1.
- **`table_demographic_bias_<site>.html` / `.pdf`** — the algorithmic-bias view:
  how each metric (VT/PBW, VT/PFVC, Ers×PBW, Ers×PFVC, Static DP, and mortality)
  varies by demographics (age per 10 yr, sex, race, height per 10 cm, SF ratio
  per 10, SOFA, BMI), as within-site standardized betas / odds ratios.
- **`table_pfvc_vs_pbw_<site>.html` / `.pdf`** — the predicted-FVC-vs-
  predicted-body-weight model: a broad-cohort linear regression
  `PFVC ~ PBW + age + sex + race`. The PBW coefficient shows how much of PFVC is
  explained by PBW, and the significant age / sex / race coefficients demonstrate
  that PBW alone does not capture predicted lung size — the core motivation for
  scaling tidal volume by PFVC instead.

**Elastance × age interaction (exploratory — `explore_elastance_age_interaction.R`):**

- **`regression_ers_age_interaction_<site>.html`** — the two mortality models with
  an elastance × age interaction (Ers×PBW and Ers×PFVC, per SD), showing the main
  effects and the interaction odds ratio.
- **`ers_age_interaction_<site>.csv`** — poolable one-row-per-model summary: the
  interaction odds ratio (per SD elastance × per 10 yr) with CI and p-value, plus
  the likelihood-ratio test against the no-interaction model (chi-square, df, p)
  and both models' AIC.
- **`ers_age_interaction_slopes_<site>.csv` / `.pdf`** — the elastance effect on
  mortality expressed on the probability scale (`marginaleffects`): the change in
  predicted mortality per 1 SD of elastance at ages 40, 60, and 80, as a table and
  a slope-by-age plot. A non-flat profile is the interaction in clinical terms.
- **`regression_ers_age_interaction_vfd_<site>.html`** and
  **`ers_age_interaction_vfd_<site>.csv`** — the same elastance × age interaction
  for the **28-day VFD competing-risks** outcome (Fine–Gray): the interaction
  subdistribution hazard ratio with CI/p and the likelihood-ratio test. (The
  unsuffixed files above are the mortality outcome.)

**Mechanical power normalization (exploratory — `explore_mechanical_power_normalization.R`):**

- **`regression_mp_normalization_<site>.html`** — the two mortality models with
  mechanical power scaled to PBW (J/min/kg) and to PFVC (J/min/L), z-scaled
  (odds ratio per 1 SD), adjusted for VT/PBW and the standard covariates.
- **`mp_normalization_<site>.csv`** — poolable one-row-per-model summary: the
  odds ratio per SD with CI and p-value, AIC, and the AIC evidence ratio relative
  to the PBW-scaled model (> 1 favors MP/PFVC), plus `n_obs` (the subset with a
  computable, mode-appropriate mechanical power).
- **`regression_mp_age_interaction_<site>.html`** — the two mortality models with
  an MP × age interaction (MP/PBW and MP/PFVC, per SD).
- **`mp_age_interaction_<site>.csv`** — poolable summary of the interaction odds
  ratio (per SD MP × per 10 yr) with CI and p-value, the likelihood-ratio test
  against the no-interaction model (chi-square, df, p), and both models' AIC.
- **`mp_age_interaction_slopes_<site>.csv` / `.pdf`** — the MP effect on mortality
  on the probability scale (`marginaleffects`): change in predicted mortality per
  1 SD of MP at ages 40, 60, and 80, as a table and a slope-by-age plot.
- **`regression_mp_normalization_vfd_<site>.html`**,
  **`mp_normalization_vfd_<site>.csv`**, and **`mp_age_interaction_vfd_<site>.csv`**
  — the MP/PBW-vs-MP/PFVC comparison and the MP × age interaction repeated for the
  **28-day VFD competing-risks** outcome (Fine–Gray): subdistribution hazard ratios
  per SD with the AIC evidence ratio, and the interaction SHR + LRT.

**Survival outputs:**

- **`km_curves_<site>.pdf`** — Kaplan–Meier 60-day survival curves by PBW/PFVC
  tercile (lowest / middle / highest) with a number-at-risk table.
- **`cox_model_summary_<site>.txt`** — the text summary of the Cox proportional-
  hazards model for 60-day survival (hazard ratios, CIs, p-values, concordance)
  on PBW/PFVC + VT/PBW + covariates.

**Figures:**

- **`evidence_ratio_heatmap_all_<site>.pdf`** — heatmap of the AIC evidence
  ratios (the visual form of `aic_comparison_all`).
- **`distribution_pbwpfvc_<site>.pdf`** — per-site PBW:PFVC distribution by
  demographic group (sex, race, age, height).
- **`consort_diagram_<site>.pdf`** — the cohort-attrition CONSORT flow diagram.
- **`bias_pbw_vs_pfvc_<site>.pdf`**, **`bias_mortality_<site>.pdf`**,
  **`bias_elastance_<site>.pdf`**, **`bias_compliance_<site>.pdf`**,
  **`bias_vfd28_<site>.pdf`**, **`bias_pbwpfvc_ratio_<site>.pdf`** — conditional-
  bias diagnostic plots (via `algorithmDiagnostics`) showing how each outcome and
  the PBW/PFVC ratio vary conditionally across demographic strata.

### `output/cross_cohort/` (script 05) — pooled across sites

Produced after each site's `regression_results_long_*.csv` files are collected
into the local `output/` tree.

- **`regression_results_all_cohorts.csv` / `.parquet`** — every site's long
  results table stacked into one file (same columns, with `site` distinguishing
  cohorts).
- **`forest_plots/forest_<analysis>.pdf`** — one forest plot per analysis
  (Mortality, Elastance, Compliance, 28-day VFDs, Static DP, Survival, …): rows
  are exposures, columns are model specifications, and each point is a cohort's
  estimate with 95% CI.
- **`forest_plots/forest_all_analyses.pdf`** — all analyses combined into a
  single multi-page PDF.
- **`forest_plots/covforest_*.pdf`** — covariate forests for the
  demographic-bias and PFVC-vs-PBW analyses (the demographic coefficients
  across cohorts, rather than the exposures).
- **`consort_diagram_pooled.pdf`** — the CONSORT funnel pooled across sites.
- **`distribution_pbwpfvc_pooled.pdf`** — the pooled PBW:PFVC distribution built
  by summing the per-site histogram exports.

---

The key cross-site artifact is **`regression_results_long_<site>.csv`** — the
file each site returns and that script 05 consumes to build the pooled forest
plots.
