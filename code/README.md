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
| `00_run_pipeline.R` | | Entry point: restores `renv`, runs the requested stages | |
| `01_cohort_identification.R` | prep | Filters the CLIF tables to the eligible cohort (ventilated, or a control cohort under `PBWPFVC_COHORT`) | |
| `02_quality_checks.R` | prep | Outlier thresholds and QC summaries | nothing in `final/`; its summaries stay in `intermediate/summary_stats/` |
| `03_variable_derivation.R` | prep | PBW, PFVC, SOFA, SF and PF ratios, tidal-volume metrics | `attrition_log_`, `dist_` |
| `04_analysis.R` | figures 1–3 | Demographic bias of PBW, mechanics, mortality regressions and survival, negative controls, E-values | `regression_results_long_`, `table1_`, `bias_`, `negative_control_`, `evalues_`, ... |
| `05_normalization_analysis.R` | figures 2–3 | PBW versus PFVC normalization of the injury metrics: discordance, reclassification, prognostic head-to-head | `norm_` |
| `10_panel_common.R` | figure 4 | The daily patient-day panel. Sourced by 21 and 25, never run; writes nothing | |
| `20_biotrauma_grid.R` | figure 4 | Time grid and the cohort-restriction knobs shared by 21–27. Sourced | |
| `21_biotrauma_panel.R` | figure 4 | Longitudinal and survival tables for the joint models | `jm_panel_summary_` |
| `22_biotrauma_fit.R` | figure 4 | One joint model per organ-injury marker | `jm_manifest_`, `jm_estimates_`, `jm_absorption_`, `jm_severity_anchor_` |
| `23_biotrauma_report.R` | figure 4 | Trajectory contrasts, hazard associations, marker movement | `jm_level_contrast_`, `jm_movement_`, `jm_association_hr_`, ... |
| `24_biotrauma_figures.R` | figure 4 | Figures from the aggregate tables only | `biotrauma_fig_` |
| `25_injury_at_horizon.R` | figure 4, robustness | Fixed-horizon comparator among survivors: what the joint model is compared against | `injury_` |
| `26_quick_lme.R` | figure 4, robustness | The longitudinal submodel alone, without the death correction | `quick_` |
| `27_control_comparison.R` | figure 4 | The divergence by lung size, arm by arm: ventilated, its SF strata, and no support read at the ventilated severity, with the severity x divergence test, and the difference-in-differences (ventilated minus control) | `jm_control_comparison_`, `jm_control_movement_`, `jm_control_did_` |
| `28_height_fingerprint.R` | Claim 5c.2 | The height fingerprint: at a fixed VT/PBW, does the marker follow the ratio's sex-reversed height curve beyond a height function the sexes share? Platelets by default, on the 7-day panel; run for the no-support control first to get the DiD | `fingerprint_`, `fingerprint_ladder_`, `fingerprint_curves_`, `fingerprint_did_` |
| `29_run_figure4.sh` | figure 4 | Runs every analysis behind figure 4 and draws it: both cohorts, all arms, the control standardised to the ventilated severity | |
| `29_run_biotrauma.sh` | figure 4, robustness | Runner for the 48-hour fits and the comparators 25 and 26; not in the pipeline | |

## Folders

- `pooling/` holds the coordinator's cross-site pooling. It is never run at a site.
  `pooled_biotrauma.R` is tracked; `pooled_estimates.R` is kept local and gitignored.
  Figure 4's estimates pool per 0.1 log units of PFVC, converted from each site's own
  SD with that site's `jm_scale_{h}_{site}.csv`, by common-effect inverse variance
  (two or three sites cannot support a random-effects variance), gated at rhat <= 1.1:

  ```bash
  mkdir -p results/fig4/MIMIC results/fig4/UCSF     # results/ is gitignored
  cp -R output/MIMIC_output/final/injury results/fig4/MIMIC/
  cp -R output/UCSF_output/final/injury  results/fig4/UCSF/
  PBWPFVC_RESULTS_ROOT=results/fig4 Rscript code/pooling/pooled_biotrauma.R
  ```
- `supplement/` holds standalone supplementary analyses, run by hand and not by the
  pipeline; each is prefixed by the block it supports and writes to `final/supplement/`.
  `xsec_dp_vtpfvc_additive.R` asks whether driving pressure is a sufficient surrogate
  for strain: driving pressure and VT/PFVC as additive predictors of mortality, and
  whether age shifts the balance between them.
- `tools/` holds developer tools: `subset_synthetic.R` makes a small synthetic CLIF
  dataset for fast test loops, `calc_external_pfvc.R` computes PFVC by arm for an
  external trial table, `creatinine_positive_control.R` checks that creatinine rises
  where kidney injury is expected (by ESRD, dialysis procedures and pressor dose; aggregates
  only, cells under 10 suppressed), and the two `migrate_*.sh` helpers move outputs written
  under older folder layouts. Delete the helpers once every site folder is migrated.
- `archive/` is gitignored local scratch.

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
