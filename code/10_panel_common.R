# =============================================================================
# Script 10 (panel): the shared daily patient-day panel
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# Builds, for the cohort in config.json (the ventilated cohort or a control), one
# baseline row per patient from analysis_cross_sectional (script 03) and a daily
# panel from the index onward: ventilator settings, the worst SF ratio, mean
# arterial pressure, vasopressors and the organ-injury labs. Sourced by
# 21_biotrauma_panel.R, which turns it into the joint-model tables; nothing else
# sources it. It writes no files.
#
# One clock. Every time is in days from the index (index_dttm, script 03). Each
# patient's follow-up ends at the first of death, the competing event (the last
# invasive-ventilation record of the stay in the ventilated cohort, so a
# reintubation counts as continuous ventilation; the first advanced-support record
# after the index in a control) and FOLLOWUP_END_D. No measurement recorded after
# that time enters any daily aggregate, so no marker value or lagged covariate
# comes from after the patient's event.
#
# Contract. The caller sources utils/config.R (config, estimate_fio2_nosupport,
# NIV_DEVICES, NOSUPPORT_DEVICES) and defines, BEFORE sourcing:
#   output_dir      the site's folder of parquet inputs (config$output_dir)
#   HORIZON         death window in days (death_day is NA past it)
#   MAX_VENT_DAY    last day since the index carried in the panel
#   FOLLOWUP_END_D  the administrative end of follow-up, days from the index
#                   (the joint-model horizon)
#   is_synthetic    TRUE at a synthetic site (simulated survival, plumbing only)
#
# After sourcing, the caller has in .GlobalEnv:
#   base           one row per patient: t0 (the index time), pfvc and pfvc_gli
#                  (both the GLI-2012 PFVC), pfvc_age25, pbw, death_day,
#                  death_time_days, extub_time_days, escalation_time_days,
#                  competing_time_days, followup_end_days, sf_index, icu_day0,
#                  ers, bmi, height_cm, age10, sex, race, sofa_total, np_sofa,
#                  SOFA components, age_grp, height_grp, disc_grp
#   wf             the spine rows (ventilator or support records) with vent_day
#   daily          per patient-day ventilator settings: vtpfvc, vtpfvc_max, vt_ml,
#                  fio2, peep, rr
#   dp_daily       daily worst driving pressure on plateau-measured days
#   maw_daily      daily median mean airway pressure (recorded values only)
#   lab_daily      daily creatinine (max), platelets (min), bilirubin (max)
#   ne_daily       daily peak norepinephrine-equivalent dose and on_pressor
#   panel_full     daily + map + sf + on_pressor + ne_equiv_peak + labs + base
#   fio2_dt, spo2_dt, maw_dt, pao2_dt   keyed data.tables for rolling joins
#   before_followup_end(df, dttm_col)   drops rows recorded after follow-up ends
#   age_breaks, rtrunc_lnorm
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)   # rolling SpO2->FiO2 join for the daily worst SF ratio
  library(tidyverse)    # loaded AFTER data.table so dplyr first()/last()/between() win
  library(arrow)
  library(here)
})

for (.need in c("output_dir", "HORIZON", "MAX_VENT_DAY", "FOLLOWUP_END_D", "is_synthetic"))
  if (!exists(.need)) stop("10_panel_common.R: caller must define `", .need, "` before sourcing.")

FIO2_PERCENT_THRESHOLD <- 1.5   # fio2_set values above are percent, below are fractions
FIO2_LOOKBACK_H        <- 4     # hours an FiO2 is carried forward to an SpO2 (chosen to capture early ventilation)

# =============================================================================
# 10a. Baseline: one row per patient (PFVC, demographics, index time, death day)
# =============================================================================
cs <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"))
if (!"pfvc_age25" %in% names(cs))
  stop("cross_sectional lacks pfvc_age25 -- re-run script 03.")

# Death within HORIZON days of the index, all causes, in continuous days
# (death_time_days) and whole days (death_day). Script 03's death_day is already days
# from the index: an expired patient with no death time died at discharge, and a
# stamp between admission and the index counts on the index day (03 logs how many).
# Synthetic CLIF mortality is unreliable, so the synthetic site draws a survival time
# instead (35% die, log-normal time to death with median 9 days, truncated to the
# horizon); it tests the plumbing only and never runs at a real site.
rtrunc_lnorm <- function(n_needed, meanlog, sdlog, lo, hi) {
  acc <- numeric(0)
  while (length(acc) < n_needed) {
    cand <- rlnorm(max(n_needed * 2L, 1000L), meanlog, sdlog)
    cand <- cand[cand > lo & cand <= hi]; acc <- c(acc, cand) }
  acc[seq_len(n_needed)]
}
stopifnot("death_day" %in% names(cs))   # days from the index, script 03
if (any(cs$death_day < 0, na.rm = TRUE))
  stop("script 03's death_day is negative for ", sum(cs$death_day < 0, na.rm = TRUE),
       " patients: it is days from the index and never before it; re-run script 03")
if (is_synthetic) {
  message("*** SYNTHETIC SITE: simulated survival (plumbing only; synthetic CLIF mortality is unreliable). ***")
  set.seed(20260615); n <- nrow(cs); died_h <- rbinom(n, 1L, 0.35)
  tte <- rep(NA_real_, n); tte[died_h == 1L] <- rtrunc_lnorm(sum(died_h), log(9), 0.95, 0.04, HORIZON)
  cs <- cs %>% mutate(death_day = if_else(died_h == 1L, floor(tte), NA_real_),
                      death_time_days = if_else(died_h == 1L, tte, NA_real_))
} else {
  cs <- cs %>% mutate(death_time_days = if_else(!is.na(death_day) & death_day <= HORIZON, death_day, NA_real_),
                      death_day = floor(death_time_days))
  message("Deaths within ", HORIZON, " days of the index: ", sum(!is.na(cs$death_time_days)), ", ",
          sum(cs$death_day == 0, na.rm = TRUE), " of them on the index day")
}
# escalation to the next level of support: script 03 writes it for the control
# cohorts; the ventilated cohort has none
if (config$cohort == "imv") cs$escalation_dttm <- as.POSIXct(NA)
age_breaks <- quantile(cs$age_at_admission, c(1/3, 2/3), na.rm = TRUE)
# Patients without a usable PFVC, age, sex, race, SOFA, height or PBW cannot enter
# any model; count each reason (a patient can fail more than one) before dropping them.
base_missing <- c(pfvc        = sum(is.na(cs$pfvc) | cs$pfvc <= 0),
                  pfvc_age25  = sum(is.na(cs$pfvc_age25) | cs$pfvc_age25 <= 0),
                  age         = sum(is.na(cs$age_at_admission)),
                  sex         = sum(is.na(cs$sex_category)),
                  race        = sum(is.na(cs$race_category)),
                  sofa        = sum(is.na(cs$sofa_total)),
                  height      = sum(is.na(cs$height_cm)),
                  pbw         = sum(is.na(cs$pbw) | cs$pbw <= 0))
base <- cs %>%
  filter(!is.na(pfvc), pfvc > 0, !is.na(pfvc_age25), pfvc_age25 > 0,
         !is.na(age_at_admission), !is.na(sex_category),
         !is.na(race_category), !is.na(sofa_total), !is.na(height_cm), !is.na(pbw), pbw > 0) %>%
  group_by(sex_category) %>% mutate(height_z = as.numeric(scale(height_cm))) %>% ungroup() %>%
  transmute(hospitalization_id, t0 = index_dttm,
            pfvc_gli = pfvc,           # the GLI-2012 PFVC, under the name 21 reads
            pfvc,                       # the same value, the denominator of the daily VT/PFVC
            pfvc_age25,                 # PFVC at age 25 (script 03), carried to the survival table
            pbw, death_day, death_time_days,
            # extubation: the last invasive-ventilation record of the stay (script 03;
            # a reintubation counts as continuous ventilation), ventilated cohort only
            extub_time_days      = as.numeric(difftime(last_imv_dttm, index_dttm, units = "days")),
            escalation_time_days = as.numeric(difftime(escalation_dttm, index_dttm, units = "days")),
            # the SF ratio at the index timepoint (script 03), which gates every hypoxemic
            # arm, and the patient's status at ICU admission (the comparison arms)
            sf_index = sf_ratio, icu_day0,
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
            # two are collinear in the hazard
            np_sofa = sofa_total - sofa_resp,
            # SOFA components for the per-marker severity anchor (20_biotrauma_grid.R)
            sofa_cv_97, sofa_coag, sofa_liver, sofa_renal,
            age_grp = cut(age_at_admission, c(-Inf, age_breaks, Inf),
                          labels = c("Young", "Middle", "Old")),
            height_grp = cut(height_z, c(-Inf, quantile(height_z, c(1/3, 2/3), na.rm = TRUE), Inf),
                             labels = c("Short", "Middle", "Tall")),
            # PBW/PFVC discordance tertile, cut on the full base before any analysis restriction
            disc_grp = cut(pbw / pfvc, c(-Inf, quantile(pbw / pfvc, c(1/3, 2/3), na.rm = TRUE), Inf),
                           labels = c("Concordant", "Mid", "Discordant")))
message("Baseline: ", nrow(base), " of ", nrow(cs), " patients; dropped for missing or non-positive ",
        paste(sprintf("%s %d", names(base_missing), base_missing), collapse = ", "))

# End of follow-up: the first of death, the competing event and FOLLOWUP_END_D.
base <- base %>%
  # a death within DEATH_ON_VENT_TOL_H of the last IMV record is a death on the
  # ventilator (utils/config.R), as in 03's ventilator-free days: no extubation
  mutate(extub_time_days = if_else(!is.na(death_time_days) & !is.na(extub_time_days) &
                                     death_time_days <= extub_time_days + DEATH_ON_VENT_TOL_H / 24,
                                   NA_real_, extub_time_days),
         competing_time_days = if (config$cohort == "imv") extub_time_days else escalation_time_days,
         followup_end_days   = pmin(death_time_days, competing_time_days, FOLLOWUP_END_D, na.rm = TRUE))
if (config$cohort == "imv") {
  message("Extubation (last IMV record of the stay) resolved for ", sum(!is.na(base$extub_time_days)), " of ",
          nrow(base), " patients; within ", FOLLOWUP_END_D, " days of the index for ",
          sum(base$extub_time_days <= FOLLOWUP_END_D, na.rm = TRUE))
} else {
  message("Control cohort (", config$cohort, "): escalation within ", FOLLOWUP_END_D, " days of the index for ",
          sum(base$escalation_time_days <= FOLLOWUP_END_D, na.rm = TRUE), " of ", nrow(base), " patients")
}
# Rows recorded after the patient's follow-up ends are dropped from every source
# before its daily reduction.
followup_end <- base %>% transmute(hospitalization_id, end_t = as.numeric(t0) + followup_end_days * 86400)
before_followup_end <- function(df, dttm_col) df %>%
  inner_join(followup_end, by = "hospitalization_id") %>%
  filter(as.numeric(.data[[dttm_col]]) <= end_t) %>%
  select(-end_t)

# =============================================================================
# 10b. Daily ventilator settings, markers and time-varying covariates
# =============================================================================
# The panel spine: rows with a set tidal volume (the ventilated cohort), or, for a
# control cohort, the rows on that cohort's devices (room air or nasal cannula for
# nosupport; HFNC, NIPPV or CPAP for niv) with a documented FiO2, every ventilator
# setting blanked (no dose exists for them).
wf <-read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet")) %>%
  estimate_fio2_nosupport() %>%
  select(hospitalization_id, recorded_dttm, device_category, tidal_volume_set, fio2_set, peep_set,
         resp_rate_set, plateau_pressure_obs, mean_airway_pressure_obs)
wf <- if (config$cohort != "imv") {
  spine_devices <- if (config$cohort == "niv") NIV_DEVICES else NOSUPPORT_DEVICES
  wf %>% filter(tolower(device_category) %in% spine_devices, !is.na(fio2_set)) %>%
    mutate(tidal_volume_set = NA_real_, peep_set = NA_real_, resp_rate_set = NA_real_,
           plateau_pressure_obs = NA_real_, mean_airway_pressure_obs = NA_real_)
} else wf %>% filter(!is.na(tidal_volume_set), tidal_volume_set > 0)
# vent_day is whole days since the index time t0, in every cohort (in a control
# cohort it counts days since the index, not days of ventilation). Rows for patients
# outside `base` drop here.
wf <- wf %>% select(-device_category) %>%
  before_followup_end("recorded_dttm") %>%
  inner_join(base %>% select(hospitalization_id, t0, pfvc), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  mutate(vtpfvc = tidal_volume_set / pfvc * 0.1)   # VT (mL) / PFVC (L) x 0.1 = percent of predicted FVC
message("Spine: ", nrow(wf), " rows on days 0-", MAX_VENT_DAY, " for ", n_distinct(wf$hospitalization_id),
        " of ", nrow(base), " baseline patients")
daily <- wf %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(vtpfvc = median(vtpfvc, na.rm = TRUE),
            vtpfvc_max = if (all(is.na(vtpfvc))) NA_real_ else max(vtpfvc, na.rm = TRUE),   # the day's peak VT/PFVC
            vt_ml = median(tidal_volume_set, na.rm = TRUE),   # absolute VT, for any other normalizer
            fio2 = median(fio2_set, na.rm = TRUE),
            peep = median(peep_set, na.rm = TRUE), rr = median(resp_rate_set, na.rm = TRUE),
            .groups = "drop")

# Daily WORST (max) driving pressure, a joint-model marker: DP = plateau - PEEP from
# RECORDED plateaus only (plateau_pressure_obs is never forward-filled), the worst
# value in each day since the index.
dp_daily <- wf %>%
  filter(!is.na(plateau_pressure_obs), !is.na(peep_set),
         plateau_pressure_obs - peep_set > 0) %>%
  mutate(dp = plateau_pressure_obs - peep_set) %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(dp = max(dp, na.rm = TRUE), .groups = "drop")
message("DP panel: ", nrow(dp_daily), " patient-days with a recorded plateau, ",
        n_distinct(dp_daily$hospitalization_id), " patients")

# Daily mean AIRWAY pressure, the numerator of the oxygenation indices. Like the
# plateau it is never forward-filled, so this is recorded values only, and its
# coverage limits the OI and OSI markers. Daily median (the typical support that
# day), reported beside the day's worst oxygenation index.
maw_daily <- wf %>%
  filter(!is.na(mean_airway_pressure_obs), mean_airway_pressure_obs > 0) %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(map_aw = median(mean_airway_pressure_obs, na.rm = TRUE), .groups = "drop")
message("Mean airway pressure panel: ", nrow(maw_daily), " patient-days with a recorded value, ",
        n_distinct(maw_daily$hospitalization_id), " patients")
# MAP: daily median (typical) from vitals.
vit <- read_parquet(file.path(output_dir, "cohort_vitals_clean.parquet")) %>%
  filter(vital_category == "map") %>%
  before_followup_end("recorded_dttm") %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(map = median(vital_value, na.rm = TRUE), .groups = "drop")

# SF ratio computed PER SpO2 measurement -- each SpO2 matched to the most recent FiO2
# within FIO2_LOOKBACK_H hours (the rolling join of script 03, FiO2 across all
# devices) -- then reduced to the daily WORST (lowest SF = worst oxygenation). SpO2 is CLAMPED to [80,97] (the linear part of the
# oxyhemoglobin dissociation curve) rather than filtered, so a fully-oxygenated day
# keeps a high (good) SF instead of being dropped, and off-curve readings are bounded.
fio2_dt <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet")) %>%
  estimate_fio2_nosupport() %>%
  filter(!is.na(fio2_set)) %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  transmute(hospitalization_id, fio2_set, t = as.numeric(recorded_dttm)) %>%
  as.data.table()
setkey(fio2_dt, hospitalization_id, t)
# Mean AIRWAY pressure and arterial PaO2 as keyed tables, for rolling either one
# onto a measurement time (the oxygenation indices built in 21_biotrauma_panel.R).
# Mean airway pressure is never forward-filled, so these are recorded values only.
maw_dt <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet")) %>%
  filter(!is.na(mean_airway_pressure_obs), mean_airway_pressure_obs > 0) %>%
  before_followup_end("recorded_dttm") %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  transmute(hospitalization_id, map_aw = mean_airway_pressure_obs, t = as.numeric(recorded_dttm)) %>%
  as.data.table()
setkey(maw_dt, hospitalization_id, t)
pao2_dt <- read_parquet(file.path(output_dir, "cohort_labs_clean.parquet")) %>%
  filter(lab_category == "po2_arterial", !is.na(lab_value_numeric), lab_value_numeric > 0) %>%
  before_followup_end("lab_result_dttm") %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  transmute(hospitalization_id, pao2 = lab_value_numeric, t = as.numeric(lab_result_dttm)) %>%
  as.data.table()
setkey(pao2_dt, hospitalization_id, t)
spo2_dt <- read_parquet(file.path(output_dir, "cohort_vitals_clean.parquet")) %>%
  filter(vital_category == "spo2", !is.na(vital_value)) %>%
  before_followup_end("recorded_dttm") %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  transmute(hospitalization_id, vent_day,
            spo2_clamped = pmin(pmax(vital_value, 80), 97), t = as.numeric(recorded_dttm)) %>%
  as.data.table()
setkey(spo2_dt, hospitalization_id, t)
sf_daily <- fio2_dt[spo2_dt, roll = FIO2_LOOKBACK_H * 3600, on = .(hospitalization_id, t)] %>%
  as_tibble() %>%
  filter(!is.na(fio2_set)) %>%
  mutate(fio2_frac = if_else(fio2_set > FIO2_PERCENT_THRESHOLD, fio2_set / 100, fio2_set),
         sf_pt = spo2_clamped / fio2_frac) %>%
  filter(is.finite(sf_pt)) %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(sf = min(sf_pt), .groups = "drop")   # daily WORST (lowest) SF
message("Worst-SF panel: ", nrow(sf_daily), " patient-days, ",
        n_distinct(sf_daily$hospitalization_id), " patients")
# Daily PEAK norepinephrine-equivalent dose (mcg/kg/min) from the series script 03
# writes: the summed dose in force (each vasopressor converted with published
# norepinephrine-equivalence factors, citation in the Methods), with a row at every
# change, every clock hour while positive and a zero row where it ends, so the
# maximum of a day's rows is that day's peak. The dose is the hemodynamic marker;
# on_pressor (a positive dose in force at some time in the day) enters the joint
# models as the previous day's pressor covariate, so a stopped or held infusion does
# not count. A day with no dose is a true zero, not a missing value.
ne_daily <- read_parquet(file.path(output_dir, "ne_equiv_admin.parquet")) %>%
  before_followup_end("admin_dttm") %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(admin_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(ne_equiv_peak = max(ne_equiv_total, na.rm = TRUE), .groups = "drop") %>%
  mutate(on_pressor = as.integer(ne_equiv_peak > 0))
message("NE-equivalent panel: ", sum(ne_daily$on_pressor), " patient-days with a positive dose, ",
        n_distinct(ne_daily$hospitalization_id[ne_daily$on_pressor == 1L]), " patients")

# Daily organ-injury labs for the joint models: the worst value of the day in the
# direction of injury (creatinine and bilirubin rise, platelets fall). Labs are
# drawn once or twice a day, so the panel carries NA on days without a draw; the
# joint models treat those as unobserved, never as unchanged.
lab_daily <- read_parquet(file.path(output_dir, "cohort_labs_clean.parquet")) %>%
  filter(lab_category %in% c("creatinine", "platelet_count", "bilirubin_total"),
         !is.na(lab_value_numeric)) %>%
  before_followup_end("lab_result_dttm") %>%
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

panel_full <- daily %>%
  left_join(vit, by = c("hospitalization_id", "vent_day")) %>%          # map (daily median)
  left_join(sf_daily, by = c("hospitalization_id", "vent_day")) %>%      # sf (daily worst)
  left_join(ne_daily, by = c("hospitalization_id", "vent_day")) %>%      # ne_equiv_peak (daily peak dose), on_pressor
  left_join(lab_daily, by = c("hospitalization_id", "vent_day")) %>%     # creatinine / platelets / bilirubin
  left_join(base, by = "hospitalization_id") %>%
  mutate(on_pressor = coalesce(on_pressor, 0L),
         ne_equiv_peak = coalesce(ne_equiv_peak, 0))
message("Daily panel: ", nrow(panel_full), " patient-days, ", n_distinct(panel_full$hospitalization_id), " patients")
