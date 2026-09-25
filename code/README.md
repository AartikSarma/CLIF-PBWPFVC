# Code

Every script here feeds a manuscript figure, runs one, or pools the results.
Scripts are numbered by the block they serve and run in number order.
`00_run_pipeline.R` is the entry point; the top-level [README](../README.md) says
how to run it. Everything reads `config/config.json` through `utils/config.R`,
which also decides where files are written.

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
| `20_biotrauma_grid.R` | figure 4 | Time grid and the cohort-restriction knobs shared by 21–27. Sourced | |
| `21_biotrauma_panel.R` | figure 4 | Longitudinal and survival tables for the joint models | `jm_panel_summary_` |
| `22_biotrauma_fit.R` | figure 4 | One joint model per organ-injury marker | `jm_manifest_`, `jm_estimates_`, `jm_scale_`, `jm_severity_anchor_` |
| `23_biotrauma_report.R` | figure 4 | Level contrasts by horizon, marker movement, the size terms with and without the death correction, hazard associations | `jm_level_contrast_`, `jm_movement_`, `jm_lme_check_`, `jm_association_hr_` |
| `24_biotrauma_figures.R` | figure 4 | Figure 4 and its DiD check, from the aggregate tables only | `biotrauma_fig_main_`, `biotrauma_fig_checks_` |
| `27_control_comparison.R` | figure 4 | The divergence by lung size, arm by arm: ventilated, its SF strata, and no support read at the ventilated severity, with the severity x divergence test, and the difference-in-differences (ventilated minus control, and minus the control hypoxemic on the index day) | `jm_control_comparison_`, `jm_control_did_`, `jm_hypoxemic_control_did_` |
| `28_height_fingerprint.R` | Claim 5c.2 | The height fingerprint: at a fixed VT/PBW, does the marker follow the ratio's sex-reversed height curve beyond a height function the sexes share? Platelets by default, on the 7-day panel; run for the no-support control first to get the DiD | `fingerprint_`, `fingerprint_ladder_`, `fingerprint_curves_`, `fingerprint_did_` |
| `29_run_figure4.R` | figure 4 | Runs every analysis behind figure 4 and draws it: both cohorts, all arms, the control standardised to the ventilated severity, SF as the positive control, and the channel breakdown (supplement) | |

## Folders

- `pooling/` holds the coordinator's cross-site pooling. It is never run at a site.
  `pooled_biotrauma.R` and `pooled_displays.R` are tracked; `pooled_estimates.R` is kept
  local and gitignored. `pooled_displays.R` draws figure 1C and the ratio's height
  channel from the GLI and Devine formulas alone, figure 2 from each site's `dose_`
  tables (the variance decomposition pooled exactly from site moments), and figure 3C's
  age-matched head-to-head from each site's `crs_channels_estimates_` and `_tests_`.
  Figure 4's estimates pool per 0.1 log units of PFVC, converted from each site's own
  SD with that site's `jm_scale_{h}_{site}.csv`, by common-effect inverse variance
  (two or three sites cannot support a random-effects variance), gated at rhat <= 1.1:

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
  curve with and without log PFVC in the ventilated and no-support cohorts, and the
  cohort x log PFVC contrast (needs scripts 01-03 run for both cohorts).
  `xsec_crs_channels.R` asks whether measured compliance scales like predicted FVC,
  input by input: the Crs exponent through the height, age, sex and race pieces of
  log PFVC, the PFVC-against-PBW head-to-head (everyone, and short women), and the
  height elasticity of Crs by sex beside GLI's and Devine's.
  `xsec_mortality_prediction.R` asks which dose or mechanics measure predicts death
  best, alone, given VT/PBW, and given VT/PBW, sex and race: VT/PBW, VT/PFVC, VT/PFVC at age 25, Ers scaled by
  each, mechanical power raw and scaled by Crs, PBW, PFVC and PFVC at age 25, and driving
  pressure, by cross-validated AUC on the patients who have every measure.
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
