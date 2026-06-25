# TTE sensitivity / exploratory inventory (removed from the federated bundle)

The federated per-site deliverable is the PRIMARY analysis only:
`01–05` (+ `05c`) → `10_tte_common`/`10_tte_engine` → `11.A` (primary TTE), `11.B`
(diagnostics), `11.X` (discordance-HTE primary). The scripts below were moved to the
gitignored `code/archive/` (kept locally, not shipped to sites).

**Plan:** once the primary is finalized (ARMA-derived strain thresholds set, `11.X`
optimized), the `11.*` sensitivities below are to be **rebuilt into a single
consolidated sensitivity script**, plus a new **strain-threshold sensitivity** (sweep the
`C_LOW`/`C_HIGH` ceilings around the ARMA-derived values). This file is the checklist of
what that consolidated script must cover.

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
| 12.H_mppbw_cate | PBW-anchored continuous CATE |

## Also archived
`05b_predicted_risk_grids.R` (secondary predicted-risk grids), `calc_external_pfvc.R`
(external-cohort PFVC helper).
