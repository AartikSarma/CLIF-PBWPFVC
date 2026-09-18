# Supplements

Tracked analyses that no runner calls. They are run by hand at the lead site, from
the repository root (`Rscript code/supplement/<name>.R`), after the stage they
depend on. Names carry the block they belong to and are not numbered, because they
have no order. Their output file names are unchanged from before the 2026-09-18
renumbering.

| Script | Needs | Question |
|---|---|---|
| `xsec_age_correction_bracket.R` | `prep` | Does PFVC's age correction add mortality discrimination beyond PBW? Brackets PBW, FVC at age 25, and PFVC (`norm_bracket_*`) |
| `xsec_subgroup_harm_exposure.R` | `prep` | Does PBW dosing concentrate harmful mechanical-stress exposure in the bias-prone subgroups? (`harm_*`, `pos_*`) |
| `xsec_dose_heterogeneity.R` | `prep` | How much dose heterogeneity hides inside the 6–8 mL/kg PBW band? (`dose_*`; read by `figures/make_strain_figures.py`) |
| `xsec_cbias_federated_export.R` | `prep` | DEFERRED, do not run or share yet: federated conditional-bias export (`cbias_export_*`) |
| `iv_height_policy.R` | `prep` | Height-instrumented PFVC-anchoring policy (`ivpolicy_*`; read by `figures/make_strain_figures.py`) |
| `injury_seven_day_feasibility.R` | 7-day panel | Can the joint-model window extend to 7 days, and with what time shape? (`sevenday_*`) |
| `injury_oi_diagnostics.R` | panel | Is an oxygenation-index signal the lung, or the ventilator setting in its own numerator? (`oi_*`) |
| `injury_dp_by_vtpbw.R`, `injury_dp_vs_vtpfvc.R` | `prep` | Driving pressure against VT/PFVC within and across VT/PBW bins (`dp_*`) |
| `injury_vtpbw_band_scan.R` | `prep` | Who enters the cohort if the VT/PBW gate is widened? (`vtpbw_band_scan_*`) |
| `tte_*.R` | `causal` | The target-trial sensitivity suite and two supplementary effect-modification analyses, inventoried below |

# TTE sensitivity / exploratory inventory (removed from the federated bundle)

The federated per-site deliverable is the PRIMARY analysis only:
`01–05` (+ `xsec_age_correction_bracket`) → `30_tte_common`/`31_tte_engine` → `35_tte_primary` (primary TTE), `36_tte_diagnostics`
(diagnostics), `37_tte_discordance_benefit` (discordance-HTE primary). The scripts below were moved to the
gitignored `code/archive/` (kept locally, not shipped to sites).

**Active lead-site supplements** (run locally at the lead site, NOT in `32_tte_run_all`, NOT
shipped to the federated sites):
- `tte_mppbw_additive_hte.R` — additive (normalizer-dependent) MP/PBW power-reduction, DR-LMTP
  CATE by discordance + PFVC; the mechanical-power sibling to `37_tte_discordance_benefit` (triangulation by exposure
  + estimator). Saves `mppbw_additive_{cate,slope}_*`.
- `33_tte_ceiling.R` — the normalizer head-to-head TTE, now the PRIMARY design and part of
  `32_tte_run_all` on the tidal-volume family (`PBWPFVC_TTE_EXPO_FAMILY=vt`); the mechanical-power
  family (`=mp`, the reframe of `tte_mppbw_additive_hte`) stays a lead-site supplement. Head-to-head clone-censor-weight
  TTE of a PFVC-anchored (X/PFVC ≤ τ) vs a PBW-anchored (X/PBW ≤ τ) ceiling, bite-matched (each binds
  the same share of post-grace days; `PBWPFVC_CEIL_BITE`, default 0.25), so positivity is symmetric by
  construction. RD = PFVC arm − PBW arm (negative favours PFVC), discordance HTE (tertiles +
  continuous CATE, SOFA-adjusted, as `37_tte_discordance_benefit` read A), per-tertile positivity (read B) and the per-lung
  correction each ceiling delivers (read C). Second design, **cap-on-top** (`tte_<fam>_cap_*`): usual
  PBW ceiling + a PFVC safety cap anchored to coincide with the PBW ceiling at the Concordant tertile's
  median discordance, vs the PBW ceiling alone — Concordant tertile is a negative control, all contrast
  accumulates in the discordant tail; positivity is NOT symmetric there, so it is the secondary,
  clinical-translation design. Pieces: `33_tte_ceiling_common.R` (guarded shared build + thresholds),
  `33_tte_ceiling_headtohead.R`, `33_tte_ceiling_cap.R`, `33_tte_ceiling_diagnostics.R` (weight-REFIT bootstrap
  for every design — written as `tte_<fam>_<design>_overall_refit_*`, the PRIMARY interval — plus the
  placebo no-differing-day contrast, `tte_<fam>_diag_*`). `PBWPFVC_CEIL_PIECES` selects pieces,
  `PBWPFVC_CEIL_REFIT_BOOT` the refit rep count.
- `tte_within_demographic_hte.R` — does the `37_tte_discordance_benefit` discordance slope survive demographics?
  Residualized (height-driven) slope vs raw + within-stratum slopes. **Result: discordance is ~99%
  demographics (R²≈0.987), residualized slope null at both sites** → repositioned from "linchpin" to
  a who-benefits characterization (the titration target is demographically patterned; physiology is
  in 05). Saves `tte_ccw_within_demo_{slope,strata,curve}_*`.
- `34_tte_titration.R` — the bedside-direct titration: an additive VT/PBW −0.5 mL/kg modified
  treatment policy (DR-LMTP), CATE/slope/gradient by discordance + PFVC. Identifiable where a static
  VT/PBW ceiling is not (feasible shift stays in support); normalizer-dependent → targets the
  misdosed; less severity-confounded than the MP sibling (`tte_mppbw_additive_hte`). Saves `vtpbw_titration_{cate,slope}_*`.
- **Split discordance-gradient sensitivity suite** — `tte_sens_common.R` (shared estimator +
  guarded engine build + lazy VR-panel builder) sourced by independently-runnable pieces:
  `tte_sens_thresholds.R` (C_LOW×C_HIGH), `tte_sens_censoring.R` (trim×cap), `tte_sens_weightmodel.R`
  (richer S/F + weight-timing), `tte_sens_deadspace.R` (VR), `tte_sens_severity_ladder.R` (daily
  SOFA / VR / DP ladder). `tte_sensitivities.R` is now a run-all wrapper over these. Output
  filenames unchanged from the old monolith (report + pooling unaffected). Run a single piece to
  add/rerun one test without rebuilding the others.

**Plan:** the `11.*` sensitivities below (still in `code/archive/`) are to be folded into the split
suite above as needed — each becomes (or extends) one `tte_sens_*.R` piece. This file is the
checklist of what that suite should eventually cover.

## TTE sensitivities to consolidate (`11.C`–`11.W`)

| script | what it varied / tested |
|---|---|
| 11.C_sens_weightcap | day-weight truncation cap (`DAYW_CAP`) |
| 11.D_sens_ceiling_grace | strain-ceiling value + grace-period length |
| 11.E_sens_rule | deviation rule: simple vs corrected (transient excursions) |
| 11.F_sens_aggregation | daily strain aggregation: median vs within-day peak |
| 11.G_sens_cumweight | cumulative-weight truncation bounds |
| 11.H_sens_ph | pH / acidosis added as a time-varying confounder |
| 11.I_sens_dp | driving pressure added as a time-varying confounder |
| 11.J_sens_pf | PaO2/FiO2 ratio added as a time-varying confounder |
| 11.K_sens_numerator | IPCW stabilizing numerator: time-only vs baseline-covariate |
| 11.L_sens_trim | common-support trim threshold (`TRIM_ALPHA`) |
| 11.M_sens_mtp | de-escalation modified-treatment-policy (dynamic decision) |
| 11.N_refit_boot | refit calibration / bootstrap reliability |
| 11.O_pfvc_overlap | PFVC positivity/overlap diagnostic by lung-size tertile |
| 11.P_sens_pfvc_overlap | PFVC overlap: adjust-vs-trim estimator |
| 11.Q_sens_dailysofa | daily (time-varying) SOFA as a confounder |
| 11.R_age_strata_balance | age-stratified weighted balance |
| 11.S_sofa_trajectory | SOFA trajectory by arm over follow-up |
| 11.T_demand_selection | deviation/adherence determinants (selection into arm) |
| 11.U_balance_reference | weighted SMD balance reference (both arms + survivorship floor) |
| 11.V_grace_compare | grace-period comparison (0/1/2/3 d) |
| 11.W_norm_positivity_compare | PFVC vs FVC_age25 normalizer positivity comparison |

**New to add:** strain-threshold sensitivity — sweep `C_LOW`/`C_HIGH` around the
ARMA-derived ceiling values; report the primary RD, the discordance gradient, and the
continuous CATE at each ceiling.

## Exploratory (MP-shift / CATE family, `12.*`) — retained locally, not for consolidation

Mostly null or superseded; kept for reference, not part of the deliverable.

| script | summary |
|---|---|
| 12.A_mp_overlap_diagnostic | MP/PFVC within-subject overlap diagnostic |
| 12.B_mp_shift_prototype | cross-sectional g-computation prototype (**superseded by 12.D**) |
| 12.C_mp_shift_ipw_longitudinal | longitudinal density-ratio IPW prototype (**superseded by 12.D**) |
| 12.D_mp_shift_lmtp | MP-shift doubly-robust LMTP, PBW vs PFVC head-to-head |
| 12.E_mp_hte_discordance | MP-reduction HTE by discordance tertile (null) |
| 12.F_mp_cate_continuous | continuous CATE of MP reduction vs discordance (null) |
| 12.G_mpcrs_cate | MP/Crs specific-power CATE by measured compliance |
| 12.H_mppbw_cate | PBW-anchored continuous CATE, *multiplicative* shift (normalizer-INVARIANT) → **promoted + reworked to active `tte_mppbw_additive_hte.R`**: switched to an *additive* MP/PBW reduction (normalizer-DEPENDENT), the LMTP/MP sibling to `37_tte_discordance_benefit`. Standalone supplement, NOT in `32_tte_run_all`. |

## Also archived
`05b_predicted_risk_grids.R` (secondary predicted-risk grids), `calc_external_pfvc.R`
(external-cohort PFVC helper).
