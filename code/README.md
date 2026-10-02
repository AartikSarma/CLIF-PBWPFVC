# Code

Every script here feeds a manuscript figure, runs one, or pools the results.
Scripts are numbered by the block they serve and run in number order.
`00_run_pipeline.R` is the entry point; the top-level [README](../README.md) says
how to run it. Everything reads `config/config.json` through `utils/config.R`,
which also decides where files are written.

**Convergence.** Figure 4's estimates are gated on the lung-size terms the figure reads (for the pfvc form, `log_pfvc_sd` and `log_pfvc_sd:vent_day`): a fit counts as converged when their R-hat is at most 1.1 (`RHAT_GATE`, `20_biotrauma_grid.R`). The R-hat of the hazard links is reported beside it (`hazard_rhat`), and `jm_lme_check_` gives each estimate from the longitudinal model fitted alone, so a reader can see whether the joint model's correction for patients leaving the panel moves the answer.

**Vasopressors.** Vasopressors are a two-part (hurdle) outcome, reported as a pair: on/off (`any_pressor`, every patient-day, a logistic mixed model, in log-odds) and the dose on the days a pressor runs (`pressor_dose`, conditional on being on a pressor that day). Being on a pressor is itself an outcome, so the dose part is read only in the ventilated cohort; the controls, the difference-in-differences and the checks figure use the on/off part.

**The VT/PBW gate's sensitivity.** `PBWPFVC_VTPBW_GATE=0` runs scripts 03, 21, 22 and 23 on the ventilated cohort without the 6–8 mL/kg band. Script 03 then writes only `analysis_cross_sectional_ungated`, the panels are named `_ungated`, and every aggregate goes to `final/ungated/<block>/` under the paper's file names. `PBWPFVC_JM_NO_DOSE=1` drops the index VT/PBW and its daily change from the longitudinal model and keeps the main fit's rows (tag `nodose_`). With `PBWPFVC_JM_SHAPE_ONLY=1` the two give the longitudinal model alone, by maximum likelihood, in minutes. `PBWPFVC_JM_PLACEBO_N=500` then refits the unadjusted model with log PFVC replaced by GLI's other indices and by 500 placebo formulas, random weightings of the same age, sex, race and height inputs (`placebo_formulas()` in `20_biotrauma_grid.R`; `jm_placebo_` tables). The gate check and the reduction-timing supplement run the same placebos through their unadjusted models. The three ventilated-versus-control supplements (`xsec_pfvc_age_control.R`, `xsec_mortality_channel_equality.R`, `xsec_strain_invariance.R`) also take `PBWPFVC_VTPBW_GATE=0` for their ventilated arm, and the first two take `PBWPFVC_DID_VTPBW=0` to drop the arm's VT/PBW term, writing to a `novtpbw/` subfolder.

**Output file names are an interface between scripts.** The pooling and figure
scripts find their inputs by file-name prefix, so a prefix changes together with
its readers, in one commit. No other site has the code yet, so names are still free
to change; that stops once sites have returned `final/` folders.

## Pipeline

| Script | Block | What it does | Writes (prefix) |
|---|---|---|---|
| `00_run_pipeline.R` | | Entry point: installs the packages in `uvr.lock`, runs the requested stages | |
| `01_cohort_identification.R` | prep | Filters the CLIF tables to the eligible cohort (ventilated, or a control cohort under `PBWPFVC_COHORT`) | |
| `02_quality_checks.R` | prep | Outlier thresholds and QC summaries | nothing in `final/`; its summaries stay in `intermediate/summary_stats/` |
| `03_variable_derivation.R` | prep | PBW, PFVC, SOFA, SF and PF ratios, tidal-volume metrics | `attrition_log_`, `dist_` |
| `04_analysis.R` | figures 1–3, tables 1–2 | Demographic bias of PBW, mechanics, mortality regressions and survival (log PFVC the primary size term), negative controls, E-values, and figure 2's delivered-strain distribution and variance decomposition | `table1_`, `regression_results_long_`, `aic_comparison_all_`, `evalues_`, `consort_diagram_`, `negative_control_`, `size_`, `dose_` (figure 2) |
| `05_normalization_analysis.R` | figure 3, supplement | PBW versus PFVC normalization of elastance and mechanical power: discordance, tertile reclassification, prognostic fit | `norm_discordance_`, `norm_prognostic_` |
| `10_panel_common.R` | figure 4 | The daily patient-day panel. Sourced by 21, never run; writes nothing | |
| `20_biotrauma_grid.R` | figure 4 | Time grid, the cohort-restriction knobs and the convergence gate shared by 21–28 and the pooling. Sourced | |
| `21_biotrauma_panel.R` | figure 4 | Longitudinal and survival tables for the joint models | `jm_panel_summary_` |
| `22_biotrauma_fit.R` | figure 4 | One joint model per organ-injury marker | `jm_manifest_`, `jm_estimates_`, `jm_scale_`, `jm_severity_anchor_` |
| `23_biotrauma_report.R` | figure 4 | Level contrasts by horizon, marker movement, the size terms with and without the death correction, hazard associations | `jm_level_contrast_`, `jm_movement_`, `jm_lme_check_`, `jm_association_hr_` |
| `24_biotrauma_figures.R` | figure 4 | Figure 4 and its DiD check, from the aggregate tables only | `biotrauma_fig_main_`, `biotrauma_fig_checks_` |
| `26_measurement_check.R` | figure 4 | Whether a marker is missing in a way that biases the comparison: does predicted lung size predict being measured (entry and day by day), and does the platelet divergence survive on bilirubin's sampling frame | `measurement_model_`, `measurement_frame_` |
| `27_control_comparison.R` | figure 4 | The divergence by lung size, arm by arm: ventilated, its SF strata, and no support read at the ventilated severity, with the severity x divergence test, and the difference-in-differences (ventilated minus control, and minus the control hypoxemic on the index day) | `jm_control_comparison_`, `jm_control_did_`, `jm_hypoxemic_control_did_` |
| `28_height_fingerprint.R` | Claim 5c.2 | The height fingerprint: at a fixed VT/PBW, does the marker follow the ratio's sex-reversed height curve beyond a height function the sexes share? Platelets by default, on the 7-day panel; run for the no-support control first to get the DiD | `fingerprint_`, `fingerprint_ladder_`, `fingerprint_curves_`, `fingerprint_did_` |
| `29_run_figure4.R` | figure 4 | Runs every analysis behind figure 4 and draws it: both cohorts, all arms, the control standardised to the ventilated severity, SF as the positive control, and the channel breakdown (supplement). The comparisons with the control run both arms on the ICU-admission clock (`PBWPFVC_JM_CLOCK=icu`, panel files `jm_*_7d_icu`); every other arm counts from the index | |

## Folders

- `pooling/` holds the coordinator's cross-site pooling. It is never run at a site.
  `pooled_biotrauma.R` and `pooled_displays.R` are tracked; `pooled_estimates.R` is kept
  local and gitignored. `pooled_displays.R` draws figure 1C and the ratio's height
  channel from the GLI and Devine formulas alone, figure 2 from each site's `dose_`
  tables (the variance decomposition pooled exactly from site moments), and figure 3C's
  age-matched head-to-head from each site's `crs_channels_estimates_` and `_tests_`.
  Figure 4's estimates pool per 0.1 log units of PFVC, converted from each site's own
  SD with that site's `jm_scale_{h}_{site}.csv`, by common-effect inverse variance
  (two or three sites cannot support a random-effects variance), gated at rhat <= 1.1.
  Each restricted or sensitivity arm (the file tags `day0_`, `sevstd_`, `sf<lo>to<hi>_`,
  `nolag_`, `rrtcause_`, `offset_`) pools only with the same arm, and the script stops on
  a second copy of any table in a site folder:

  ```bash
  mkdir -p results/fig4/MIMIC results/fig4/UCSF     # results/ is gitignored
  cp -R output/MIMIC_output/final/injury output/MIMIC_output/final/supplement results/fig4/MIMIC/
  cp -R output/UCSF_output/final/injury  output/UCSF_output/final/supplement  results/fig4/UCSF/
  PBWPFVC_RESULTS_ROOT=results/fig4 uvr run code/pooling/pooled_biotrauma.R
  ```

  The same script pools the supplement's contrasts from each site's `supplement/`: the
  mortality control contrast and its GLI channel breakdown (converted to per 0.1 log
  units with each site's exported SD; the pieces' agreement tested by multivariate
  common-effect pooling of each site's covariance) and the compliance channels (the
  head-to-head AIC differences summed across sites).
- `supplement/` holds standalone supplementary analyses, run by hand and not by the
  pipeline; each is prefixed by the block it supports and writes to `final/supplement/`.
  `xsec_dp_vtpfvc_additive.R` asks whether driving pressure is a sufficient surrogate
  for strain: driving pressure and VT/PFVC as additive predictors of mortality, and
  whether age shifts the balance between them. `xsec_pfvc_age_control.R` asks whether
  log PFVC carries part of age's mortality gradient under ventilation only: the age
  curve with and without log PFVC in the ventilated arm on invasive ventilation at ICU
  admission (`icu_day0`) and the no-support cohort, and the cohort x log PFVC contrast,
  with the control censored at escalation (needs scripts 01-03 run for both cohorts).
  `xsec_intubation_overlap.R` asks whether a propensity score for intubation could
  weight the no-support controls to the ventilated arm: it fits the score from the 24
  hours before ICU admission (ventilated at ICU admission against each control) and
  reports overlap, balance and effective sample size, with no outcome model (run by
  29_run_figure4.R after the anchors).
  `xsec_crs_channels.R` asks whether measured compliance scales like predicted FVC,
  input by input: the Crs exponent through the height, age, sex and race pieces of
  log PFVC, the PFVC-against-PBW head-to-head (everyone, and short women), and the
  height elasticity of Crs by sex beside GLI's and Devine's (figure 3C).
  00_run_pipeline.R's cross_sectional stage runs `xsec_crs_channels.R`,
  `xsec_dp_vtpfvc_additive.R` and `xsec_mortality_prediction.R` after 05.
  `xsec_mortality_prediction.R` asks which dose or mechanics measure predicts death
  best, alone, given VT/PBW, and given VT/PBW, sex and race: VT/PBW, VT/PFVC, VT/PFVC at age 25, Ers scaled by
  each, mechanical power raw and scaled by Crs, PBW, PFVC and PFVC at age 25, and driving
  pressure, by cross-validated AUC on the patients who have every measure.
  `xsec_age_form_check.R` asks whether 04 and 05's linear age leaves curvature that
  log PFVC then carries: it refits every model with a size exposure with
  `ns(age10, 4)` in place of `age10`, nothing else changed, and reports how far each
  size coefficient moves and how much of each claimed AIC advantage survives. Not in
  the runner.
  `xsec_vtpbw_gate_collider.R` asks whether the VT/PBW 6-8 gate makes predicted lung
  size look sicker than it is. VT/PBW rises with the PBW/PFVC ratio and falls with
  illness, so selecting on it or adjusting for it can link the two. On an index chosen
  without the gate (03's `analysis_ungated_index`), it compares the size term's SOFA
  slope with and without the gate and the VT/PBW term, with bootstrap intervals for
  the change, split into respiratory and other SOFA, and then asks the same of
  in-hospital death, with and without 04's severity terms. Not in the runner.
  `xsec_vt_reduction_timing.R` follows it up: when clinicians turn the tidal volume
  down in the first 72 hours, which arrow from predicted size carries the reduction?
  The candidates are a high starting setting, a high driving pressure, or a fall in
  compliance (the only one that would mean small lungs caused the illness). It
  tabulates what preceded each first reduction by PFVC quartile, decomposes size's
  link to reductions in a discrete-time hazard, and asks whether size predicts a fall
  in compliance at 24-48 h. Not in the runner.
  `xsec_mortality_channel_equality.R` asks whether sex and race predict death as their
  share of predicted lung size says they should: with age held by `ns(age, 4)`, each
  group's coefficient divided by its GLI shift in log PFVC is tested against the
  height slope (sex, Black, Other, and jointly), in the ventilated arm, the no-support
  control and their difference (needs scripts 01-03 run for both cohorts). Not in the
  runner.
  `xsec_strain_invariance.R` asks whether ventilation's mortality gradient follows
  the strain error or any index of the same demographics, with severity-only models:
  the site's tidal-volume rule (slope of set VT on PBW before the 6-8 gate), the
  ventilated-minus-control contrast for log PBW/PFVC and log PFVC, and 500 random
  weightings of age, sex, race and height beside GLI's ratio, PFVC and pieces. The
  cross-site test (the contrast against the dosing rule) is a pooling step. Not in
  the runner.
- `tools/` holds developer tools: `subset_synthetic.R` makes a small synthetic CLIF
  dataset for fast test loops, `calc_external_pfvc.R` computes PFVC by arm for an
  external trial table, `creatinine_positive_control.R` checks that creatinine rises
  where kidney injury is expected (by ESRD, dialysis procedures and pressor dose; aggregates
  only).
- `archive/` is gitignored local scratch.

Every file each script writes to `final/`, its reader and the claim it supports are
listed in `docs/output_manifest.md` (kept outside the repository with the other
manuscript documents).

## Where aggregates go

`final/` is sorted by block: `cross_sectional/` (`03`-`05`), `injury/` (`21`-`28`),
`supplement/` and `controls/`. A script never builds the path itself: it calls
`final_dir_for("<block>")` from `utils/config.R`.

## Cohorts and where files go

`PBWPFVC_COHORT` selects the cohort: `imv` (default, the ventilated analytic cohort),
`nosupport` (room air or nasal cannula only: the negative control) or `niv`
(high-flow or noninvasive ventilation first). The noninvasive cohort is built only on
request and is never drawn as a control, because NIPPV delivers large, unlimited
positive-pressure volumes. `utils/config.R` sends
a control cohort's patient-level files to `intermediate/controls/<cohort>/` and all
its aggregates to `final/controls/`, whatever the block, inside the site's one output
folder, and tags the file names `<site>_<cohort>`.

## What was removed, and how to get it back

On 2026-09-24 the 48-hour runs went: `25_injury_at_horizon.R`, `26_quick_lme.R` and
`29_run_biotrauma.sh` (the fixed-horizon comparator, the longitudinal submodel alone
and their runner), with their pooling blocks. So did the VT/PFVC companion fits of
figure 4 (22-24 still accept the `vtpfvc` form). Bring them back from the commit
before their removal with `git log --diff-filter=D -- code/25_injury_at_horizon.R`.

On 2026-09-19 the repository was pruned to the scripts above. Figure 5 (causal
inference) is not in the paper yet, so the target trial emulation and the preference
instrument went, with the sensitivity suite, the one-off diagnostics, the deferred
conditional-bias export, the strain-figure material and the Python report builders.
All of it is in the tag `pre-prune-2026-09-19`, under the names that tag's
`code/README.md` lists:

```bash
git show pre-prune-2026-09-19:code/README.md                 # the index of what was there
git checkout pre-prune-2026-09-19 -- code/33_tte_ceiling.R   # bring one file back
```

Restoring the target trial emulation also needs its block in `utils/config.R`
(`FINAL_BLOCKS` gains `"causal"`) and its stage in `00_run_pipeline.R`; both are in
the tag.
