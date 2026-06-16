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

The per-site pipeline is five R scripts run in order by `code/00_run_pipeline.R`,
which first restores the `renv` environment, then runs 01 -> 05 each as a clean
subprocess. Cross-cohort pooling (`code/pooled_estimates.R`) is a separate step run
centrally by the coordinator after every site returns its results — it is not part
of the per-site runner.

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
  logistic mortality regressions across the exposure specifications;
  linear models for the mechanics outcomes — elastance,
  compliance, driving pressure, normalized elastance (Ers×PBW, Ers×PFVC),
  mechanical power, and normalized mechanical power (MP/PBW, MP/PFVC) — each
  regressed across the exposure specifications;
  **28-day VFDs modeled as a competing-risks outcome** (extubation vs death, via
  Fine–Gray subdistribution models, per Yehya & Harhay 2019) rather than a
  continuous value; an AIC / evidence-ratio comparison; 60-day survival (Cox +
  Kaplan-Meier by PBW/PFVC tercile); demographic-bias models (each metric vs
  demographics); the predicted-FVC-vs-predicted-body-weight model (regressing
  PFVC on PBW plus age, sex, and race to show that PBW alone does not capture
  predicted lung size); conditional-bias diagnostic plots; and the unified
  long-format effect-estimate table that feeds the cross-cohort step. Every
  exposure-to-outcome model (mortality, mechanics, 28-day VFD, 60-day survival) is
  reported both **demographic-adjusted** (+ age/sex/race) and **unadjusted**
  (demographics dropped, illness severity retained), tagged by an `adjustment`
  column in the long table.
- **05 — Normalization analysis.** A systematic PBW-vs-PFVC comparison of the
  physiologic injury metrics whose landmark papers use PBW or no size reference —
  Goligher's normalized elastance (Ers x PBW), Gattinoni's MP/PBW, and Amato's
  driving pressure. It reports the physiologic **discordance** between the PBW- and
  PFVC-normalized metrics and injury-tertile **reclassification** (with direct
  age-interaction tests), the **prognostic head-to-head** with a form-vs-physiology
  decomposition, an **encompassing test** of whether PFVC adds information PBW
  misses, and age/lung-size **interaction ladders**. Every model is reported
  demographic-adjusted and unadjusted; a `dp <= 0` QC filter and centered
  predictors control collinearity. Writes per-site `norm_*` outputs.
- **Cross-cohort aggregation (`code/pooled_estimates.R`, run centrally).**
  Site-agnostic and run by the coordinator after every site returns its `final/`
  folder — not part of the per-site runner, and kept out of the repository. It
  discovers every site's `regression_results_long_*.csv` (script 04) and
  `norm_*.csv` (script 05) under a results root (one subfolder per site, each site's
  `final/` renamed to the site name; root set via `PBWPFVC_RESULTS_ROOT`),
  stacks/pools them, and renders forest plots (one per analysis + a combined PDF),
  the adjusted-vs-unadjusted comparison forest, covariate forests, pooled evidence
  ratios and E-values, the pooled normalization summaries, and pooled consort and
  PBW:PFVC distribution figures. Pooled outputs are written to an `All sites/`
  subfolder of the results root (excluded from discovery on re-run).

## What the output files are

Per-site output lands under `output/<site_name>_output/`. The whole `output/` tree
is gitignored — only aggregated, n >= 10-suppressed results are written. Pooled
cross-cohort output is written by `code/pooled_estimates.R` to an `All sites/`
subfolder of the central results root (by default the local, gitignored `results/`
folder into which each site's results are copied).

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

### `final/` (scripts 03–05) — per-site results

These are the deliverables. The first group is aggregated and cell-suppressed
(n >= 10) and is what each site returns to the consortium; the rest are
human-readable tables and figures for local review.

**Federated tables (return these):**

- **`regression_results_long_<site>.csv` / `.parquet`** — the key cross-site
  artifact. One row per model term across **every** model in script 04, with
  columns: `site`, `term`, `estimate`, `conf_low`, `conf_high`, `std_error`,
  `statistic`, `p_value`, `estimate_type` (OR / HR / Beta), `analysis`,
  `model_spec` (the exposure specification), `model_family`, `adjustment`
  (demographic-`adjusted` / `unadjusted`), `formula`, and `n_obs`.
  `pooled_estimates.R` stacks these across sites to build the forest plots (the
  adjusted estimates are primary, with an adjusted-vs-unadjusted comparison).
- **`norm_*_<site>.csv`** (script 05) — the normalization analysis outputs:
  discordance summary + demographics (`norm_discordance_*`), injury-tertile
  reclassification (`norm_reclassification_*`), prognostic fit and form-vs-physiology
  (`norm_prognostic_fit_*`, `norm_form_vs_physiology_*`), the encompassing test
  (`norm_encompassing_*`), and the age/size interaction ladder
  (`norm_age_interaction_*`). `pooled_estimates.R` pools these across cohorts.
- **`aic_comparison_all_<site>.csv`** — model-fit comparison: per
  `outcome` × `exposure`, the `AIC`, `delta_AIC`, and `evidence_ratio` relative
  to the VT/PBW reference within that outcome, plus `is_reference` and the
  truncated ratio / display label.
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
- **`regression_ers_<site>.html`**, **`regression_crs_<site>.html`**,
  **`regression_dp_<site>.html`**, **`regression_ers_pbw_<site>.html`**,
  **`regression_ers_pfvc_<site>.html`**, **`regression_mp_<site>.html`**,
  **`regression_mp_pbw_<site>.html`**, **`regression_mp_pfvc_<site>.html`** —
  linear regressions of the mechanics outcomes (elastance, compliance, static
  driving pressure, normalized elastance Ers×PBW / Ers×PFVC, mechanical power, and
  normalized mechanical power MP/PBW / MP/PFVC), each across the same five exposure
  specifications. The mechanical-power outcomes require a recorded peak inspiratory
  pressure, so their N is the (smaller) subset where MP is computable.
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

**Normalization analysis (script 05) — `norm_*_<site>.csv` / `.pdf`:**

- **`norm_discordance_summary_*`, `norm_discordance_demographics_*`,
  `norm_discordance_age_*`** — distribution of the PBW/PFVC size discordance and its
  demographic / age patterning.
- **`norm_reclassification_*`** — fraction of patients changing injury tertile when
  switching PBW->PFVC, by demographic group and by age.
- **`norm_prognostic_fit_*`, `norm_form_vs_physiology_*`** — the PBW-vs-PFVC
  prognostic head-to-head (AIC, C-statistic) and the decomposition into model-form
  vs physiology effects, each demographic-adjusted and unadjusted.
- **`norm_encompassing_*`** — the encompassing test: does PFVC add prognostic
  information beyond PBW (and vice versa), by AIC and LRT.
- **`norm_age_interaction_ladder_*`, `norm_age_interaction_slopes_*`** — whether each
  mechanic's (DP / Ers / MP) effect is modified by age (recoil) and lung volume.

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

### `All sites/` (`pooled_estimates.R`, run centrally) — pooled across sites

Produced after each site's `final/` folder (renamed to the site name) is copied
into the central results root (by default the local `results/` folder); written to
the `All sites/` subfolder there.

- **`regression_results_all_cohorts.csv` / `.parquet`** — every site's long
  results table stacked into one file (same columns, with `site` distinguishing
  cohorts).
- **`forest_plots/forest_<analysis>.pdf`** — one forest plot per analysis
  (Mortality, Elastance, Compliance, 28-day VFDs, Static DP, Survival, …), at the
  demographic-adjusted (primary) estimate: rows are exposures, columns are model
  specifications, and each point is a cohort's estimate with 95% CI.
- **`forest_plots/forest_all_analyses.pdf`** — all analyses combined into a
  single multi-page PDF.
- **`forest_plots/forest_adjusted_vs_unadjusted.pdf`** + **`adjusted_vs_unadjusted_pooled.csv`**
  — the headline exposures with both adjustment levels side by side.
- **`forest_plots/covforest_*.pdf`** — covariate forests for the
  demographic-bias and PFVC-vs-PBW analyses (the demographic coefficients
  across cohorts, rather than the exposures).
- **`norm_encompassing_pooled.{csv,pdf}`, `norm_form_vs_physiology_pooled.csv`,
  `norm_prognostic_fit_pooled.csv`, `norm_discordance_pooled.csv`** — the pooled
  normalization analyses (summed delta-AIC across cohorts).
- **`consort_diagram_pooled.pdf`** — the CONSORT funnel pooled across sites.
- **`distribution_pbwpfvc_pooled.pdf`** — the pooled PBW:PFVC distribution built
  by summing the per-site histogram exports.

---

The key cross-site artifacts are **`regression_results_long_<site>.csv`** (outcome
analyses, script 04) and **`norm_*_<site>.csv`** (normalization analyses, script
05) — the files each site returns and that `pooled_estimates.R` consumes to build
the pooled forest plots and summaries.
