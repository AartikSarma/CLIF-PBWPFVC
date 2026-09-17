# =============================================================================
# Script 10 (panel): the shared daily patient-day panel
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# Builds the baseline table and the day-level panel that both the longitudinal
# target trial emulation (10_tte_common.R and its 11.* leaves) and the biotrauma
# joint models (13_biotrauma_*.R) run on. Factored out of 10_tte_common.R on
# 2026-09-13 so the two analyses share one exposure grid, one set of daily
# confounders, and one pair of event definitions (death day, extubation day).
#
# This file performs NO file writes and applies NO analysis-specific exclusion.
# In particular the TTE's structural-positivity restriction (patients who cannot
# reach the strain ceiling even at the VT floor) stays in 10_tte_common.R: it is
# right for a ceiling emulation and wrong for a mechanism analysis, where those
# patients carry the dose-response.
#
# Contract. The caller defines, BEFORE sourcing:
#   output_dir    intermediate folder of the site (parquet inputs)
#   HORIZON       outcome horizon in days (death_day is NA past it)
#   MAX_VENT_DAY  last vent day carried in the panel
#   is_synthetic  TRUE at the synthetic site (simulated survival, plumbing only)
#   PANEL_NORM    "pfvc" or "pfvc_age25": the column that becomes `pfvc` in base
#                 and the denominator of the daily vtpfvc. The TTE passes its
#                 normalizer switch; the joint models pass "pfvc".
# Optional environment: PBWPFVC_TTE_VTPBW="lo,hi" re-selects the index over a
# wider VT/PBW band (positivity sensitivity); the panel then sets
# panel_cohort_tag = "vtpbw_lo_hi" so the caller can isolate its outputs.
#
# After sourcing, the caller has in .GlobalEnv:
#   base           one row per patient: t0, pfvc (switched), pfvc_gli (always the
#                  GLI value), pfvc_age25, pbw, death_day, ers, bmi, height_cm,
#                  age10, sex, race, sofa_total, age_grp, height_grp, disc_grp,
#                  imv_extub_day
#   daily          per patient-day ventilator settings: vtpfvc, vtpfvc_max, vt_ml,
#                  fio2, peep, rr
#   dp_daily       daily worst driving pressure on plateau-measured days
#   panel_full     daily + map + sf + on_pressor + ne_equiv_peak + labs + base
#   panel          panel_full restricted to complete core covariates
#   ph_daily, ph_daily_art, pao2_daily   indication-driven gases, kept separate
#   lab_daily      daily creatinine (max), platelets (min), bilirubin (max)
#   ne_daily       daily peak norepinephrine-equivalent dose
#   panel_drop_stats()   function: the missing-covariate drop summary for a
#                  given panel_full (the TTE recomputes it after its restriction)
#   panel_cohort_tag, age_breaks, rtrunc_lnorm
#
# Still to add (plan step 2): a daily RRT flag from the CLIF crrt_therapy table,
# which script 01 does not yet load; the creatinine joint model censors at RRT start.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)   # rolling SpO2->FiO2 join for the daily worst SF ratio
  library(tidyverse)    # loaded AFTER data.table so dplyr first()/last()/between() win
  library(arrow)
  library(here)
})

for (.need in c("output_dir", "HORIZON", "MAX_VENT_DAY", "is_synthetic", "PANEL_NORM"))
  if (!exists(.need)) stop("10_panel_common.R: caller must define `", .need, "` before sourcing.")
stopifnot(PANEL_NORM %in% c("pfvc", "pfvc_age25"))

# =============================================================================
# 10a. Baseline (PFVC, demographics, t0, death day) -- mortality handling
#      identical to scripts 06-08 incl. synthetic-only workaround.
# =============================================================================
# Cohort source. DEFAULT = the primary cross-sectional cohort (index selected at VT/PBW
# 6-8 mL/kg). Set PBWPFVC_TTE_VTPBW="lo,hi" to RE-SELECT the index over a wider VT/PBW band
# from the pre-gate per-timepoint dataset (analysis_all_eligible_timepoints) -- a positivity
# sensitivity that lets smaller lungs be dosed below 6 mL/kg and actually reach the strain
# ceiling. Same two-tier index rule and SF<315 hypoxemia gate as script 03; only the VT/PBW
# band changes. The caller isolates outputs under final/<panel_cohort_tag>/.
vtpbw_band <- Sys.getenv("PBWPFVC_TTE_VTPBW", "")
panel_cohort_tag <- ""
if (nzchar(vtpbw_band)) {
  bw <- suppressWarnings(as.numeric(strsplit(vtpbw_band, ",")[[1]]))
  stopifnot(length(bw) == 2, !any(is.na(bw)), bw[1] < bw[2])
  message("*** BROADER COHORT: re-selecting the index over VT/PBW [", bw[1], ", ", bw[2],
          "] from analysis_all_eligible_timepoints (positivity sensitivity) ***")
  ae <- read_parquet(file.path(output_dir, "analysis_all_eligible_timepoints.parquet"))
  qual <- ae %>% filter(has_all_data, vtpbw >= bw[1], vtpbw <= bw[2], sf_ratio < 315)
  imv_start <- ae %>% group_by(hospitalization_id) %>%
    summarise(imv_start_dttm = min(recorded_dttm), .groups = "drop")
  t1 <- qual %>% left_join(imv_start, by = "hospitalization_id") %>%
    filter(!is.na(dp), recorded_dttm <= imv_start_dttm + lubridate::hours(6)) %>%
    group_by(hospitalization_id) %>% slice_min(recorded_dttm, n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(-imv_start_dttm)
  t2 <- qual %>% filter(!hospitalization_id %in% t1$hospitalization_id) %>%
    group_by(hospitalization_id) %>% slice_min(recorded_dttm, n = 1, with_ties = FALSE) %>% ungroup()
  # The per-timepoint table lacks bmi (script 03 attaches weight and BMI to the
  # cross-sectional table only); derive it here exactly as 03 does so `base` is
  # the same shape on both cohort paths. Broken on main since 004c45a made bmi a
  # base column; fixed 2026-09-13.
  cs <- bind_rows(t1, t2) %>%
    left_join(read_parquet(file.path(output_dir, "cohort_weights.parquet")),
              by = "hospitalization_id") %>%
    mutate(bmi = if_else(!is.na(weight_kg) & height_cm > 0,
                         weight_kg / (height_cm / 100)^2, NA_real_))
  panel_cohort_tag <- paste0("vtpbw_", bw[1], "_", bw[2])
  message("  broader cohort: ", nrow(cs), " patients; output tag ", panel_cohort_tag)
} else {
  cs <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"))
}

if (!"pfvc_age25" %in% names(cs))
  stop("cross_sectional lacks pfvc_age25 -- re-run script 03 (it now derives the structural normalizer).")

rtrunc_lnorm <- function(n_needed, meanlog, sdlog, lo, hi) {
  acc <- numeric(0)
  while (length(acc) < n_needed) {
    cand <- rlnorm(max(n_needed * 2L, 1000L), meanlog, sdlog)
    cand <- cand[cand > lo & cand <= hi]; acc <- c(acc, cand) }
  acc[seq_len(n_needed)]
}
if (is_synthetic) {
  message("*** SYNTHETIC SITE: simulated survival (plumbing only; synthetic CLIF mortality is unreliable). ***")
  set.seed(20260615); n <- nrow(cs); died_h <- rbinom(n, 1L, 0.35)
  tte <- rep(NA_real_, n); tte[died_h == 1L] <- rtrunc_lnorm(sum(died_h), log(9), 0.95, 0.04, HORIZON)
  cs <- cs %>% mutate(death_day = if_else(died_h == 1L, floor(tte), NA_real_),
                      death_time_days = if_else(died_h == 1L, tte, NA_real_))
} else {
  cs <- cs %>% mutate(idx = as.numeric(difftime(death_dttm, recorded_dttm, units = "days")),
                      death_day = if_else(!is.na(idx) & idx >= 0 & idx <= HORIZON, floor(idx), NA_real_),
                      # unfloored death time, for the joint models' sub-daily grid
                      death_time_days = if_else(!is.na(idx) & idx >= 0 & idx <= HORIZON, idx, NA_real_))
}
# escalation to invasive ventilation (the never-intubated control only; NA otherwise)
if (!"escalation_dttm" %in% names(cs)) cs$escalation_dttm <- as.POSIXct(NA)
age_breaks <- quantile(cs$age_at_admission, c(1/3, 2/3), na.rm = TRUE)
base <- cs %>%
  filter(!is.na(pfvc), pfvc > 0, !is.na(pfvc_age25), pfvc_age25 > 0,
         !is.na(age_at_admission), !is.na(sex_category),
         !is.na(race_category), !is.na(sofa_total), !is.na(height_cm), !is.na(pbw), pbw > 0) %>%
  group_by(sex_category) %>% mutate(height_z = as.numeric(scale(height_cm))) %>% ungroup() %>%
  transmute(hospitalization_id, t0 = recorded_dttm,
            pfvc_gli = pfvc,            # ALWAYS the GLI-2012 PFVC (the joint models' normalizer)
            pfvc = .data[[PANEL_NORM]], # the CEILING normalizer (switch); downstream stays normalizer-agnostic
            pfvc_age25,                 # ALWAYS carry the structural normalizer (12 secondary CATE + positivity)
            pbw, death_day, death_time_days,
            escalation_time_days = as.numeric(difftime(escalation_dttm, recorded_dttm, units = "days")),
            # measured mechanics at the index timepoint (plateau subset only, so often NA).
            # ers (cmH2O/L) x the size normalizer is SPECIFIC elastance: near-constant across
            # lungs if the normalizer is right about this patient's aerated volume (Chiumello),
            # ABOVE that value when the aerated lung is smaller than predicted. Unlike PBW/PFVC
            # it is MEASURED, so it varies within demographic strata -- the one mis-sizing index
            # in this project that is not a deterministic function of age/sex/race/height.
            # bmi rides along because respiratory-system elastance includes the chest wall.
            ers, bmi, height_cm,
            age10 = age_at_admission / 10, sex_category, race_category, sofa_total,
            # non-respiratory SOFA: the joint models enter severity beside log SF, and
            # the respiratory component is computed from that same SF ratio, so the
            # two are collinear in the hazard (MIMIC: log SF x death R-hat 4.3)
            np_sofa = sofa_total - sofa_resp,
            age_grp = cut(age_at_admission, c(-Inf, age_breaks, Inf),
                          labels = c("Young", "Middle", "Old")),
            height_grp = cut(height_z, c(-Inf, quantile(height_z, c(1/3, 2/3), na.rm = TRUE), Inf),
                             labels = c("Short", "Middle", "Tall")),
            # PBW/PFVC discordance: at fixed VT/PBW the strain-limiting contrast lives ENTIRELY
            # in the high-discordance (small-lung) patients, so this tertile is the TTE's own
            # version of the MP discordance HTE (12.E) -- the effect should concentrate in the
            # Discordant band, where positivity is also thinnest. Cut on the switched pfvc, on
            # the full base (before any analysis-specific restriction).
            disc_grp = cut(pbw / pfvc, c(-Inf, quantile(pbw / pfvc, c(1/3, 2/3), na.rm = TRUE), Inf),
                           labels = c("Concordant", "Mid", "Discordant")))

# =============================================================================
# 10b. Daily exposure + time-varying confounder panel
# =============================================================================
# The panel spine: rows with a set tidal volume (the analytic cohort), or, for the
# never-intubated control, the HFNC / NIPPV / CPAP rows with a documented FiO2 and
# every ventilator setting blanked (no dose exists for them)
wf <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet")) %>%
  estimate_fio2_nosupport() %>%
  select(hospitalization_id, recorded_dttm, device_category, tidal_volume_set, fio2_set, peep_set,
         resp_rate_set, plateau_pressure_obs)
wf <- if (config$cohort != "imv") {
  spine_devices <- if (config$cohort == "niv") NIV_DEVICES else NOSUPPORT_DEVICES
  wf %>% filter(tolower(device_category) %in% spine_devices, !is.na(fio2_set)) %>%
    mutate(tidal_volume_set = NA_real_, peep_set = NA_real_, resp_rate_set = NA_real_, plateau_pressure_obs = NA_real_)
} else wf %>% filter(!is.na(tidal_volume_set), tidal_volume_set > 0)
wf <- wf %>% select(-device_category) %>%
  inner_join(base %>% select(hospitalization_id, t0, pfvc), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  mutate(vtpfvc = tidal_volume_set / pfvc * 0.1)
daily <- wf %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(vtpfvc = median(vtpfvc, na.rm = TRUE),
            vtpfvc_max = if (all(is.na(vtpfvc))) NA_real_ else max(vtpfvc, na.rm = TRUE),
            vt_ml = median(tidal_volume_set, na.rm = TRUE),   # absolute VT, for any other normalizer
            fio2 = median(fio2_set, na.rm = TRUE),
            peep = median(peep_set, na.rm = TRUE), rr = median(resp_rate_set, na.rm = TRUE),
            .groups = "drop")   # vtpfvc_max = within-day PEAK strain, for the Q9 aggregation sensitivity

# Daily WORST (max) driving pressure for the [T5c] sensitivity: DP = plateau - PEEP
# from RECORDED plateaus only (plateau_pressure_obs is never forward-filled), taking
# the worst value in each 24h vent-day. Clinicians titrate VT to plateau (ARMA) and
# driving pressure (Amato NEJM 2015), so lagged worst-DP is a behaviorally-real
# driver of the deviation decision -- a time-varying confounder, run as a sensitivity.
dp_daily <- wf %>%
  filter(!is.na(plateau_pressure_obs), !is.na(peep_set),
         plateau_pressure_obs - peep_set > 0) %>%
  mutate(dp = plateau_pressure_obs - peep_set) %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(dp = max(dp, na.rm = TRUE), .groups = "drop")
message("DP panel: ", nrow(dp_daily), " patient-days with a recorded plateau, ",
        n_distinct(dp_daily$hospitalization_id), " patients")
# MAP: daily median (typical) from vitals.
vit <- read_parquet(file.path(output_dir, "cohort_vitals_clean.parquet")) %>%
  filter(vital_category == "map") %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(map = median(vital_value, na.rm = TRUE), .groups = "drop")

# SF ratio computed PER SpO2 measurement -- each SpO2 matched to the most recent FiO2
# within 4h (the canonical rolling join from script 03, fio2 across ALL modes) -- then
# reduced to the daily WORST (lowest SF = worst oxygenation, the value most likely to
# drive a tidal-volume decision). SpO2 is CLAMPED to [80,97] (the linear part of the
# oxyhemoglobin dissociation curve) rather than filtered, so a fully-oxygenated day
# keeps a high (good) SF instead of being dropped, and off-curve readings are bounded.
fio2_dt <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet")) %>%
  estimate_fio2_nosupport() %>%
  filter(!is.na(fio2_set)) %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  transmute(hospitalization_id, fio2_set, t = as.numeric(recorded_dttm)) %>%
  as.data.table()
setkey(fio2_dt, hospitalization_id, t)
spo2_dt <- read_parquet(file.path(output_dir, "cohort_vitals_clean.parquet")) %>%
  filter(vital_category == "spo2", !is.na(vital_value)) %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  transmute(hospitalization_id, vent_day,
            spo2_clamped = pmin(pmax(vital_value, 80), 97), t = as.numeric(recorded_dttm)) %>%
  as.data.table()
setkey(spo2_dt, hospitalization_id, t)
sf_daily <- fio2_dt[spo2_dt, roll = 4 * 3600, on = .(hospitalization_id, t)] %>%
  as_tibble() %>%
  filter(!is.na(fio2_set)) %>%
  mutate(fio2_frac = if_else(fio2_set > 1.5, fio2_set / 100, fio2_set),
         sf_pt = spo2_clamped / fio2_frac) %>%
  filter(is.finite(sf_pt)) %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(sf = min(sf_pt), .groups = "drop")   # daily WORST (lowest) SF
message("Worst-SF panel: ", nrow(sf_daily), " patient-days, ",
        n_distinct(sf_daily$hospitalization_id), " patients")
med <- read_parquet(file.path(output_dir, "cohort_meds.parquet")) %>%
  filter(med_group == "vasoactives") %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(admin_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  distinct(hospitalization_id, vent_day) %>% mutate(on_pressor = 1L)

# Daily PEAK norepinephrine-equivalent dose (mcg/kg/min) from the per-administration
# table script 03 writes (catecholamines standardized to mcg/kg/min, vasopressin at
# 2.5 per unit/min, dobutamine 0). The binary on_pressor flag above is the TTE's
# confounder; the dose is the joint models' hemodynamic injury marker. A day with no
# vasoactive administration is a true zero, not a missing value.
ne_daily <- read_parquet(file.path(output_dir, "ne_equiv_admin.parquet")) %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(admin_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(ne_equiv_peak = max(ne_equiv_total, na.rm = TRUE), .groups = "drop")
message("NE-equivalent panel: ", nrow(ne_daily), " patient-days with a vasoactive dose, ",
        n_distinct(ne_daily$hospitalization_id), " patients")

# Daily organ-injury labs for the joint models: the worst value of the day in the
# direction of injury (creatinine and bilirubin rise, platelets fall). Labs are
# drawn once or twice a day, so the panel carries NA on days without a draw; the
# joint models treat those as unobserved, never as unchanged.
lab_daily <- read_parquet(file.path(output_dir, "cohort_labs_clean.parquet")) %>%
  filter(lab_category %in% c("creatinine", "platelet_count", "bilirubin_total"),
         !is.na(lab_value_numeric)) %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(lab_result_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(creatinine = { v <- lab_value_numeric[lab_category == "creatinine"];      if (length(v)) max(v) else NA_real_ },
            platelets  = { v <- lab_value_numeric[lab_category == "platelet_count"];  if (length(v)) min(v) else NA_real_ },
            bilirubin  = { v <- lab_value_numeric[lab_category == "bilirubin_total"]; if (length(v)) max(v) else NA_real_ },
            .groups = "drop")
message("Lab panel: creatinine on ", sum(!is.na(lab_daily$creatinine)), ", platelets on ",
        sum(!is.na(lab_daily$platelets)), ", bilirubin on ", sum(!is.na(lab_daily$bilirubin)),
        " patient-days (", n_distinct(lab_daily$hospitalization_id), " patients)")

# True extubation from the FULL IMV course (any ventilator mode), independent of
# tidal_volume_set. The old proxy (last volume-targeted day + 1) fires at the END of
# volume-control ventilation, but patients are routinely weaned onto pressure support
# before extubation -- at our sites ~28-38% of IMV patients have their last IMV day on a
# non-volume mode, so the old proxy truncated ventilation/liberation early. The IPCW
# weight still freezes at the last VOLUME-TARGETED day (no exposure info accrues on non-
# volume days); only the competing-risk LIBERATION time uses this true extubation.
imv_extub <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet")) %>%
  select(hospitalization_id, recorded_dttm, device_category) %>%
  filter(tolower(device_category) == "imv") %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  group_by(hospitalization_id) %>%
  summarise(imv_extub_day = max(vent_day) + 1L, .groups = "drop")
base <- base %>% left_join(imv_extub, by = "hospitalization_id")
if (config$cohort != "imv") {
  message("Control cohort (", config$cohort, "): escalation within the window for ",
          sum(!is.na(base$escalation_time_days)), " of ", nrow(base), " patients")
} else {
  message("IMV-course extubation derived for ", sum(!is.na(base$imv_extub_day)), " of ",
          nrow(base), " patients (every index-IMV patient should resolve).")
}

panel_full <- daily %>%
  left_join(vit, by = c("hospitalization_id", "vent_day")) %>%          # map (daily median)
  left_join(sf_daily, by = c("hospitalization_id", "vent_day")) %>%      # sf (daily worst)
  left_join(med, by = c("hospitalization_id", "vent_day")) %>%
  left_join(ne_daily, by = c("hospitalization_id", "vent_day")) %>%      # ne_equiv_peak (daily peak dose)
  left_join(lab_daily, by = c("hospitalization_id", "vent_day")) %>%     # creatinine / platelets / bilirubin
  left_join(base, by = "hospitalization_id") %>%
  mutate(on_pressor = coalesce(on_pressor, 0L),
         ne_equiv_peak = coalesce(ne_equiv_peak, 0),
         keep = (config$cohort != "imv" | (is.finite(vtpfvc) & is.finite(peep) & is.finite(rr))) &
                is.finite(fio2) & is.finite(sf) & is.finite(map))
panel <- panel_full %>% filter(keep) %>% select(-keep)
# Extubation for the liberation endpoint = the IMV-course day (imv_extub_day, derived
# above, carried in via base); the IPCW weight-freeze point stays at the last volume-
# targeted day (arm idsum last_vent). The old panel-based extub proxy is retired ([T12]).
message("Panel: ", nrow(panel), " patient-days, ", n_distinct(panel$hospitalization_id), " patients")

# [T11] Missing-covariate diagnostic (the drop above must not be silent). ICU charting
# should make vtpfvc/fio2/peep/rr/sf/map dense, but a dropped vent-day (a) breaks the
# lag-1 spacing of the time-varying confounders (lag() is row-based, so a gap makes
# "yesterday" actually 2-3 days back) and (b) if it falls at the END of a stay, pulls
# the extubation proxy (= last observed vent-day + 1) earlier, biasing ventilation
# duration and the liberation competing-risk endpoint. Report the magnitude so it is
# auditable; a non-trivial share here is a signal to LOCF-fill rather than list-delete.
# A function, because the TTE recomputes it after its structural restriction; the
# CSV write is performed by 11.B_diagnostics.R.
panel_drop_stats <- function(pf) {
  drop_diag <- pf %>% group_by(hospitalization_id) %>%
    summarise(days_total = n(), days_kept = sum(keep),
              last_full = max(vent_day),
              last_kept = { k <- vent_day[keep]; if (length(k)) max(k) else NA_real_ },
              .groups = "drop")
  tibble(
    patient_days_total     = sum(drop_diag$days_total),
    patient_days_dropped   = sum(drop_diag$days_total - drop_diag$days_kept),
    frac_days_dropped      = sum(drop_diag$days_total - drop_diag$days_kept) / sum(drop_diag$days_total),
    patients_any_drop      = sum(drop_diag$days_kept < drop_diag$days_total),
    patients_extub_shifted = sum(!is.na(drop_diag$last_kept) & drop_diag$last_kept < drop_diag$last_full),
    patients_lost_entirely = sum(drop_diag$days_kept == 0))
}
panel_drop_summary <- panel_drop_stats(panel_full)
message(sprintf(paste0("Missing-covariate drops: %.1f%% of patient-days (%d of %d); ",
        "%d patients lose >=1 day, %d have extubation pulled earlier, %d lost entirely."),
        100 * panel_drop_summary$frac_days_dropped, panel_drop_summary$patient_days_dropped,
        panel_drop_summary$patient_days_total, panel_drop_summary$patients_any_drop,
        panel_drop_summary$patients_extub_shifted, panel_drop_summary$patients_lost_entirely))

# Daily WORST (lowest) pH for the [T5b] sensitivity -- the acidosis nadir is the gas
# most likely to drive a tidal-volume increase (permissive hypercapnia -> deviation
# from a low-VT arm). Read both gases once, then derive two daily worst-pH series:
# (1) POOLED arterial-equivalent = arterial + venous imputed as venous +
# 0.05 (venous pH runs ~0.03-0.05 below arterial); (2) ARTERIAL-ONLY (drops the
# imputed venous values, to confirm the +0.05 imputation isn't doing the work).
# Kept SEPARATE from the core panel: gas sampling is indication-driven (MNAR).
gas <- read_parquet(file.path(output_dir, "cohort_labs_clean.parquet")) %>%
  filter(lab_category %in% c("ph_arterial", "ph_venous"), !is.na(lab_value_numeric)) %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(lab_result_dttm, t0, units = "days"))),
         ph_art_eq = if_else(lab_category == "ph_venous",
                             lab_value_numeric + 0.05, lab_value_numeric)) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY)
ph_daily <- gas %>% group_by(hospitalization_id, vent_day) %>%        # pooled (art + venous+0.05)
  summarise(ph = min(ph_art_eq, na.rm = TRUE), .groups = "drop")      # worst (lowest) pH
ph_daily_art <- gas %>% filter(lab_category == "ph_arterial") %>%     # arterial only
  group_by(hospitalization_id, vent_day) %>%
  summarise(ph = min(lab_value_numeric, na.rm = TRUE), .groups = "drop")  # worst (lowest) pH
message("pH panel: ", nrow(ph_daily), " patient-days (pooled), ",
        n_distinct(ph_daily$hospitalization_id), " patients; ",
        n_distinct(ph_daily_art$hospitalization_id), " with an arterial gas")

# Daily worst (lowest) PaO2 for the PF-ratio sensitivity (10i2). Oxygenation drives the
# FiO2/PEEP titration that co-determines whether a low-VT arm can be held, so lagged
# P/F is a behaviorally-real confounder of the deviation decision -- run as a SENSITIVITY
# (arterial gas is indication-driven / MNAR, like pH). PF is formed in the sensitivity
# block by joining this PaO2 to the panel's daily FiO2 (same fraction/percent handling).
pao2_daily <- read_parquet(file.path(output_dir, "cohort_labs_clean.parquet")) %>%
  filter(lab_category == "po2_arterial", !is.na(lab_value_numeric)) %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(lab_result_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(pao2 = min(lab_value_numeric, na.rm = TRUE), .groups = "drop")   # worst-of-day
message("PaO2 panel: ", nrow(pao2_daily), " patient-days, ",
        n_distinct(pao2_daily$hospitalization_id), " patients with an arterial PaO2")
