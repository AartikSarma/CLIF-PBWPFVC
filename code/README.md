# Code

Scripts are numbered by the manuscript block they serve and run in number order.
`00_run_pipeline.R` is the entry point; the top-level [README](../README.md) says
how to run it. Everything here reads `config/config.json` through `utils/config.R`,
which also decides where files are written.

**Output file names are an interface between scripts.** The pooling and figure
scripts find their inputs by file-name prefix, so a prefix changes together with
its readers, in one commit. No other site has the code yet, so names are still free
to change; that stops once sites have returned `final/` folders.

## Pipeline

| Script | Block | What it does | Writes to `final/` |
|---|---|---|---|
| `00_run_pipeline.R` | | Entry point: restores `renv`, runs the requested stages | |
| `01_cohort_identification.R` | prep | Filters the CLIF tables to the eligible cohort (ventilated, or a control cohort under `PBWPFVC_COHORT`) | |
| `02_quality_checks.R` | prep | Outlier thresholds and QC summaries | `lab_summary_`, `vital_summary_` |
| `03_variable_derivation.R` | prep | PBW, PFVC, SOFA, SF and PF ratios, tidal-volume metrics | `attrition_log_`, `dist_` |
| `04_analysis.R` | figures 1–3 | Demographic bias of PBW, mechanics, mortality regressions and survival, negative controls, E-values | `regression_results_long_`, `table1_`, `bias_`, `negative_control_`, `evalues_`, ... |
| `05_normalization_analysis.R` | figures 2–3 | PBW versus PFVC normalization of the injury metrics: discordance, reclassification, prognostic head-to-head | `norm_` |
| `10_panel_common.R` | shared | The daily patient-day panel used by both the injury and the causal blocks. Sourced, never run; writes nothing | |
| `20_biotrauma_grid.R` | figure 4 | Time grid and the cohort-restriction knobs shared by 21–27. Sourced | |
| `21_biotrauma_panel.R` | figure 4 | Longitudinal and survival tables for the joint models | `jm_panel_summary_` |
| `22_biotrauma_fit.R` | figure 4 | One joint model per organ-injury marker | `jm_manifest_`, `jm_estimates_`, `jm_absorption_`, `jm_severity_anchor_` |
| `23_biotrauma_report.R` | figure 4 | Trajectory contrasts, hazard associations, marker movement | `jm_level_contrast_`, `jm_movement_`, `jm_association_hr_`, ... |
| `24_biotrauma_figures.R` | figure 4 | Figures from the aggregate tables only | `biotrauma_fig_` |
| `25_injury_at_horizon.R` | figure 4 | Fixed-horizon comparator (survivors only) | `injury_` |
| `26_quick_lme.R` | figure 4 | The longitudinal submodel alone, as a fast check | `quick_` |
| `27_control_comparison.R` | figure 4 | The divergence by lung size, arm by arm: ventilated, its SF strata, no support, matched no support, noninvasive | `jm_control_comparison_`, `jm_control_movement_` |
| `28_biotrauma_summary.R` | figure 4 | One table collecting a run's results | `overnight_summary_` |
| `29_run_biotrauma.sh` | figure 4 | Runner for 21–28 on one cohort | |
| `29_run_controls.sh` | figure 4 | Runner for the control cohorts: `build`, `anchors`, `fits` | |
| `29_run_sites.sh` | figure 4 | Runs `29_run_biotrauma.sh` over several sites in turn | |
| `30_tte_common.R` | figure 5 | Target trial emulation: shared setup and estimators. Sourced; writes nothing | |
| `31_tte_engine.R` | figure 5 | Cache in front of `30_tte_common.R`. Sourced by every TTE script | |
| `32_tte_run_all.R` | figure 5 | Builds the engine once, then runs 33–37 | |
| `33_tte_ceiling.R` | figure 5 | PRIMARY: PFVC-anchored versus PBW-anchored tidal-volume ceiling, bite-matched. Pieces: `33_tte_ceiling_{common,headtohead,cap,diagnostics,elastance}.R` | `tte_vt_` |
| `34_tte_titration.R` | figure 5 | Co-primary: one-step titration toward PFVC dosing | `vtpbw_titration_` |
| `35_tte_primary.R` | figure 5 | Strain-limiting policy: overall effect, bootstrap, E-value, subgroups | `tte_ccw_overall_`, `tte_ccw_subgroup_`, ... |
| `36_tte_diagnostics.R` | figure 5 | Positivity, weights, balance | `tte_ccw_diagnostics_`, `tte_ccw_balance_`, ... |
| `37_tte_discordance_benefit.R` | figure 5 | Does the policy benefit the most mis-dosed patients most? | `tte_ccw_disc_` |
| `38_iv_preference.R` | figure 5 | Practice-preference instrument for the strain-limiting strategy | `iv_preference_` |

## Folders

- `supplement/` holds tracked lead-site analyses that no runner calls. They are
  named by the block they belong to (`xsec_`, `injury_`, `tte_`, `iv_`) and are not
  numbered because they have no order. [Its README](supplement/README.md) is the
  inventory, including what was archived and why.
- `pooling/` holds the coordinator's cross-site pooling. It is never run at a site.
- `tools/` holds developer tools: `subset_synthetic.R` makes a small synthetic CLIF
  dataset for fast test loops, and `calc_external_pfvc.R` computes PFVC by arm for
  an external trial table.
- `archive/` is gitignored local scratch.

## Cohorts and where files go

`PBWPFVC_COHORT` selects the cohort: `imv` (default, the ventilated analytic cohort),
`nosupport` (room air or nasal cannula only: the negative control) or `niv`
(high-flow or noninvasive ventilation first: a point on the strain gradient, not a
clean control, because tidal volumes there are uncontrolled). `utils/config.R` sends
a control cohort's patient-level files to `intermediate/controls/<cohort>/` and its
aggregates to `final/controls/`, inside the site's one output folder, and tags the
file names `<site>_<cohort>`.

## Old names

The scripts were renumbered on 2026-09-18. `git log --follow` works across the move.

| Old | New |
|---|---|
| `05c_age_correction_bracket.R` | `supplement/xsec_age_correction_bracket.R` |
| `07_subgroup_harm_exposure.R` | `supplement/xsec_subgroup_harm_exposure.R` |
| `08_dose_heterogeneity.R` | `supplement/xsec_dose_heterogeneity.R` |
| `09_iv_policy.R` | `supplement/iv_height_policy.R` |
| `10_tte_common.R`, `10_tte_engine.R` | `30_tte_common.R`, `31_tte_engine.R` |
| `11_run_all.R` | `32_tte_run_all.R` |
| `11.M_ceiling_tte.R`, `11_ceiling_*.R` | `33_tte_ceiling.R`, `33_tte_ceiling_*.R` |
| `11_vtpbw_titration.R` | `34_tte_titration.R` |
| `11.A_primary.R` | `35_tte_primary.R` |
| `11.B_diagnostics.R` | `36_tte_diagnostics.R` |
| `11.X_discordance_benefit.R` | `37_tte_discordance_benefit.R` |
| `11.Y_mppbw_additive_hte.R` | `supplement/tte_mppbw_additive_hte.R` |
| `11.Z_within_demographic_hte.R` | `supplement/tte_within_demographic_hte.R` |
| `11_sensitivities.R`, `11_sens_*.R` | `supplement/tte_sensitivities.R`, `supplement/tte_sens_*.R` |
| `13_biotrauma_{grid,panel,fit,report,figures}.R` | `20`–`24_biotrauma_*.R` |
| `13_injury_at_horizon.R`, `13_quick_lme.R` | `25_injury_at_horizon.R`, `26_quick_lme.R` |
| `13_control_comparison.R`, `13_overnight_summary.R` | `27_control_comparison.R`, `28_biotrauma_summary.R` |
| `13_run_{biotrauma,controls,sites}.sh` | `29_run_{biotrauma,controls,sites}.sh` |
| `13_seven_day_feasibility.R`, `13_oi_diagnostics.R`, `13_dp_by_vtpbw.R`, `13_dp_vs_vtpfvc.R`, `13_vtpbw_band_scan.R` | `supplement/injury_*.R` |
| `14_preference_iv.R` | `38_iv_preference.R` |
| `cbias_federated_export.R` | `supplement/xsec_cbias_federated_export.R` (deferred) |
| `cbias_pooled_plots.R`, `pooled_biotrauma.R` | `pooling/` |
| `00_subset_synthetic.R`, `calc_external_pfvc.R` | `tools/` |
| `SENSITIVITY_INVENTORY.md` | `supplement/README.md` |
