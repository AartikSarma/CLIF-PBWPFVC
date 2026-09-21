# =============================================================================
# Script 21 (panel): Biotrauma joint models -- the longitudinal and survival tables
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# Builds the two tables the biotrauma joint models (22_biotrauma_fit.R) consume,
# from the shared daily panel of 10_panel_common.R:
#
#   jm_long_{48h|72h|24h|7d}.parquet  one row per patient-period (index period 0
#                          to the horizon) with the organ-injury markers observed
#                          in that period, the previous period's strain and
#                          confounders, and the patient's baseline covariates
#   jm_surv_{tag}.parquet  one row per patient: competing-risk coding of death vs
#                          extubation within the horizon (tie counts as death),
#                          the RRT start, every baseline covariate, and the GLI
#                          channel pieces of log PFVC
#   jm_meta_{tag}.rds      horizon, cohort tag, counts
#
# and one aggregate table for the deliverable:
#
#   final/jm_panel_summary_{tag}_{site}.csv  patients, patient-periods and events
#                                       per marker; RRT censoring counts; the size
#                                       of the plateau-measured subset
# The tag is set by 20_biotrauma_grid.R (PBWPFVC_JM_GRID, PBWPFVC_JM_HORIZON_H);
# panels for different horizons sit side by side.
#
# Design (docs/joint_model_plan_2026-09.md, sections 3 and 4):
#   * every index-IMV patient enters at day 0; no survival-based restriction, and
#     no structural-positivity exclusion
#   * the exposure is the PREVIOUS day's median VT/PBW (the clinician's dose) as a
#     deviation from the patient's mean over the course, taken by joining on
#     vent_day - 1 so a missing day gives a missing lag rather than a two-day-old
#     one; PBW/PFVC discordance (log, centred) is the effect modifier; VT/PFVC,
#     its running mean and the count of days above 11% are carried for the
#     sensitivity forms
#   * markers on the day of observation: creatinine (daily max), platelets (daily
#     min), bilirubin (daily max), SF ratio (daily worst), driving pressure (daily
#     max, plateau-measured days only), NE-equivalent dose (daily peak)
#   * creatinine is censored at the first CRRT record: days on or after RRT start
#     are set to missing, and a patient already on CRRT at the index has no
#     creatinine trajectory at all
#   * rows are truncated at the event day, as a joint model requires
#
# Horizon: PBWPFVC_JM_HORIZON days (default 7; 14 is the sensitivity). The shared
# panel is built with the TTE's 28-day death window so death_day and the
# extubation day are identical objects in both analyses.
#
# Usage: Rscript code/21_biotrauma_panel.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(arrow)
  library(here)
})
rm(list = ls())
source("utils/config.R")

site_name  <- config$site_name
output_dir <- config$output_dir
final_dir  <- final_dir_for("injury")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

# --- shared-panel contract (identical windows to the TTE, so the event objects match)
HORIZON      <- 28L
MAX_VENT_DAY <- 27L
is_synthetic <- grepl("^synthetic_clif", site_name)   # any synthetic site (synthetic_clif, synthetic_clif_b, ...)
PANEL_NORM   <- "pfvc"          # the joint models always normalize to GLI PFVC
source(here("code", "10_panel_common.R"))

# --- time grid. PRIMARY = "6h": six-hour periods over the first 48 hours, where
# the ventilator acts on the injured lung most, a linear trend is plausible, and
# death and extubation have not yet removed a third of the cohort. The dense
# markers (SF, driving pressure, NE-equivalent dose) have a value in most
# periods; labs are drawn about daily and carry one or two values. SENSITIVITY =
# "daily": one row per ventilator day over PBWPFVC_JM_HORIZON days (7), with the
# day spline in the fit. Both grids share the 22_biotrauma_fit.R / _report.R code
# through the `period` index and the numeric time `vent_day` (days).
source(here("code", "20_biotrauma_grid.R"))   # JM_GRID, STEP_H, STEP, JM_HORIZON, N_PERIODS, h_suffix
STRAIN_CEILING <- 11   # VT/PFVC % above which a period counts toward the cumulative-strain exposure
message("=== 21_biotrauma_panel: grid ", JM_GRID, ", horizon ", JM_HORIZON, " days (",
        N_PERIODS, " periods), site ", site_name, " ===")

# =============================================================================
# The period panel: one row per patient-period, either the shared daily panel
# (grid = daily) or a six-hour reduction of the same sources (grid = 6h)
# =============================================================================
# ---- oxygenation index, on whichever grid is in force
# OI  = FiO2(%) x mean AIRWAY pressure / PaO2      (arterial gas; indication-driven, sparse)
# OSI = FiO2(%) x mean airway pressure / SpO2      (the saturation analogue, dense)
#     = 100 x mean airway pressure / SF, exactly, since SF = SpO2 / FiO2 as a fraction
# Higher is worse for both. Each index is computed AT A MEASUREMENT from values taken
# at the same time (user, 2026-09-21: an index whose numerator and denominator come
# from different hours is uninterpretable):
#   OSI at each SpO2: SpO2 clamped to 80-97 as for SF, FiO2 from the last 4 h (as for
#                     SF), mean airway pressure the nearest RECORDED value within
#                     OXY_MATCH_H hours either side (never forward-filled)
#   OI at each PaO2:  FiO2 from the last 4 h, mean airway pressure as above
# and a period keeps its WORST (highest) index. A measurement with no mean airway
# pressure within the window has no index, so coverage is the thing to read first.
#
# CAUTION, and it is the reason OI was parked in the first place: mean airway
# pressure is a ventilator setting that the exposure moves arithmetically. A
# bigger tidal volume at the same PEEP and compliance raises mean airway
# pressure, so OI can worsen with VT/PFVC through its own numerator, with
# nothing happening in the lung. Read OI beside its components (this script
# writes the mean airway pressure by day), and treat PEEP-per-PFVC as the
# mediator it is, not as a nuisance.
OXY_MATCH_H <- as.numeric(Sys.getenv("PBWPFVC_OXY_MATCH_H", "1"))
maw_near <- copy(maw_dt)[, t_maw := t]                     # keeps the matched reading's own time
near_maw <- function(pts) {                                # nearest mean airway pressure, within the window
  as_tibble(maw_near[pts, roll = "nearest", on = .(hospitalization_id, t)]) %>%
    filter(!is.na(map_aw), abs(t - t_maw) <= OXY_MATCH_H * 3600)
}
osi_pts <- fio2_dt[spo2_dt, roll = 4 * 3600, on = .(hospitalization_id, t)][!is.na(fio2_set)] %>%
  near_maw() %>%
  mutate(fio2_frac = if_else(fio2_set > 1.5, fio2_set / 100, fio2_set),
         osi = 100 * map_aw / (spo2_clamped / fio2_frac)) %>%
  filter(is.finite(osi), osi > 0) %>% select(hospitalization_id, t, osi)
oi_pts <- fio2_dt[copy(pao2_dt), roll = 4 * 3600, on = .(hospitalization_id, t)][!is.na(fio2_set)] %>%
  near_maw() %>%
  mutate(fio2_pct = if_else(fio2_set > 1.5, fio2_set, fio2_set * 100), oi = fio2_pct * map_aw / pao2) %>%
  filter(is.finite(oi), oi > 0) %>% select(hospitalization_id, t, oi)
message("Oxygenation indices from matched measurements (mean airway pressure within ", OXY_MATCH_H, " h): ",
        nrow(osi_pts), " OSI and ", nrow(oi_pts), " OI values")
t0_secs <- base %>% transmute(hospitalization_id, t0n = as.numeric(t0))
worst_index <- function(pts, col, period_hours, max_period) pts %>%
  inner_join(t0_secs, by = "hospitalization_id") %>%
  mutate(period = as.integer(floor((t - t0n) / 3600 / period_hours))) %>%
  filter(period >= 0L, period <= max_period) %>%
  group_by(hospitalization_id, period) %>%
  summarise(!!col := max(.data[[col]]), .groups = "drop")
if (JM_GRID == "daily") {
  pf  <- panel_full %>% mutate(period = as.integer(vent_day))
  dpp <- dp_daily   %>% mutate(period = as.integer(vent_day))
  # the day's median mean airway pressure (descriptive) beside the day's worst matched indices
  oxy <- maw_daily %>% transmute(hospitalization_id, period = as.integer(vent_day), map_aw) %>%
    full_join(worst_index(osi_pts, "osi", 24, MAX_VENT_DAY), by = c("hospitalization_id", "period")) %>%
    full_join(worst_index(oi_pts,  "oi",  24, MAX_VENT_DAY), by = c("hospitalization_id", "period"))
  # the competing event: extubation (analytic cohort) or escalation to invasive ventilation (control)
  extub_time <- base %>% transmute(hospitalization_id,
                                   extub_time = if (config$cohort != "imv") escalation_time_days else as.numeric(imv_extub_day))
} else {
  per <- function(dttm, t0) as.integer(floor(as.numeric(difftime(dttm, t0, units = "hours")) / STEP_H))
  b0  <- base %>% select(hospitalization_id, t0)
  MAXP <- as.integer(ceiling(MAX_VENT_DAY * 24 / STEP_H))
  # ventilator settings per period from the same waterfall rows the daily panel uses (wf)
  set_p <- wf %>% mutate(period = per(recorded_dttm, t0)) %>%
    filter(period >= 0L, period <= MAXP) %>%
    group_by(hospitalization_id, period) %>%
    summarise(vtpfvc = median(vtpfvc, na.rm = TRUE),
              vtpfvc_max = if (all(is.na(vtpfvc))) NA_real_ else max(vtpfvc, na.rm = TRUE),
              vt_ml = median(tidal_volume_set, na.rm = TRUE),
              fio2 = median(fio2_set, na.rm = TRUE), peep = median(peep_set, na.rm = TRUE),
              rr = median(resp_rate_set, na.rm = TRUE), .groups = "drop")
  dpp <- wf %>% mutate(period = per(recorded_dttm, t0)) %>%
    filter(period >= 0L, period <= MAXP, !is.na(plateau_pressure_obs), !is.na(peep_set),
           plateau_pressure_obs - peep_set > 0) %>%
    mutate(dp = plateau_pressure_obs - peep_set) %>%
    group_by(hospitalization_id, period) %>% summarise(dp = max(dp), .groups = "drop")
  map_p <- read_parquet(file.path(output_dir, "cohort_vitals_clean.parquet")) %>%
    filter(vital_category == "map") %>% inner_join(b0, by = "hospitalization_id") %>%
    mutate(period = per(recorded_dttm, t0)) %>% filter(period >= 0L, period <= MAXP) %>%
    group_by(hospitalization_id, period) %>% summarise(map = median(vital_value, na.rm = TRUE), .groups = "drop")
  # SF per SpO2 measurement (FiO2 rolled back within 4 h, SpO2 clamped to 80-97, as in the daily panel), worst per period
  t0_num <- b0 %>% transmute(hospitalization_id, t0n = as.numeric(t0))
  sf_p <- fio2_dt[spo2_dt, roll = 4 * 3600, on = .(hospitalization_id, t)] %>% as_tibble() %>%
    filter(!is.na(fio2_set)) %>%
    mutate(fio2_frac = if_else(fio2_set > 1.5, fio2_set / 100, fio2_set), sf_pt = spo2_clamped / fio2_frac) %>%
    filter(is.finite(sf_pt)) %>% inner_join(t0_num, by = "hospitalization_id") %>%
    mutate(period = as.integer(floor((t - t0n) / 3600 / STEP_H))) %>% filter(period >= 0L, period <= MAXP) %>%
    group_by(hospitalization_id, period) %>% summarise(sf = min(sf_pt), .groups = "drop")
  # the period's median recorded mean airway pressure (descriptive; the indices are matched above)
  maw_p <- wf %>% mutate(period = per(recorded_dttm, t0)) %>%
    filter(period >= 0L, period <= MAXP, !is.na(mean_airway_pressure_obs), mean_airway_pressure_obs > 0) %>%
    group_by(hospitalization_id, period) %>%
    summarise(map_aw = median(mean_airway_pressure_obs, na.rm = TRUE), .groups = "drop")
  pao2_p <- read_parquet(file.path(output_dir, "cohort_labs_clean.parquet")) %>%
    filter(lab_category == "po2_arterial", !is.na(lab_value_numeric)) %>%
    inner_join(b0, by = "hospitalization_id") %>%
    mutate(period = per(lab_result_dttm, t0)) %>% filter(period >= 0L, period <= MAXP) %>%
    group_by(hospitalization_id, period) %>%
    summarise(pao2 = min(lab_value_numeric, na.rm = TRUE), .groups = "drop")
  oxy <- maw_p %>%
    full_join(worst_index(osi_pts, "osi", STEP_H, MAXP), by = c("hospitalization_id", "period")) %>%
    full_join(worst_index(oi_pts,  "oi",  STEP_H, MAXP), by = c("hospitalization_id", "period"))
  med_p <- read_parquet(file.path(output_dir, "cohort_meds.parquet")) %>%
    filter(med_group == "vasoactives") %>% inner_join(b0, by = "hospitalization_id") %>%
    mutate(period = per(admin_dttm, t0)) %>% filter(period >= 0L, period <= MAXP) %>%
    distinct(hospitalization_id, period) %>% mutate(on_pressor = 1L)
  ne_p <- read_parquet(file.path(output_dir, "ne_equiv_admin.parquet")) %>%
    inner_join(b0, by = "hospitalization_id") %>%
    mutate(period = per(admin_dttm, t0)) %>% filter(period >= 0L, period <= MAXP) %>%
    group_by(hospitalization_id, period) %>% summarise(ne_equiv_peak = max(ne_equiv_total), .groups = "drop")
  lab_p <- read_parquet(file.path(output_dir, "cohort_labs_clean.parquet")) %>%
    filter(lab_category %in% c("creatinine", "platelet_count", "bilirubin_total"), !is.na(lab_value_numeric)) %>%
    inner_join(b0, by = "hospitalization_id") %>%
    mutate(period = per(lab_result_dttm, t0)) %>% filter(period >= 0L, period <= MAXP) %>%
    group_by(hospitalization_id, period) %>%
    summarise(creatinine = { v <- lab_value_numeric[lab_category == "creatinine"];      if (length(v)) max(v) else NA_real_ },
              platelets  = { v <- lab_value_numeric[lab_category == "platelet_count"];  if (length(v)) min(v) else NA_real_ },
              bilirubin  = { v <- lab_value_numeric[lab_category == "bilirubin_total"]; if (length(v)) max(v) else NA_real_ },
              .groups = "drop")
  # the competing event at period resolution: extubation = last IMV period + 1
  # (analytic cohort); escalation to invasive ventilation (the control)
  extub_time <- if (config$cohort != "imv") {
    base %>% transmute(hospitalization_id, extub_time = escalation_time_days)
  } else read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet")) %>%
    select(hospitalization_id, recorded_dttm, device_category) %>%
    filter(tolower(device_category) == "imv") %>% inner_join(b0, by = "hospitalization_id") %>%
    mutate(period = per(recorded_dttm, t0)) %>% filter(period >= 0L, period <= MAXP) %>%
    group_by(hospitalization_id) %>% summarise(extub_time = (max(period) + 1L) * STEP, .groups = "drop")
  pf <- set_p %>%
    left_join(map_p, by = c("hospitalization_id", "period")) %>%
    left_join(sf_p,  by = c("hospitalization_id", "period")) %>%
    left_join(med_p, by = c("hospitalization_id", "period")) %>%
    left_join(ne_p,  by = c("hospitalization_id", "period")) %>%
    left_join(lab_p, by = c("hospitalization_id", "period")) %>%
    left_join(base,  by = "hospitalization_id") %>%
    mutate(on_pressor = coalesce(on_pressor, 0L), ne_equiv_peak = coalesce(ne_equiv_peak, 0),
           vent_day = period * STEP)
  dpp <- dpp %>% mutate(vent_day = period * STEP)
  message("Six-hour panel: ", nrow(pf), " patient-periods, ", n_distinct(pf$hospitalization_id),
          " patients; SF on ", sum(!is.na(pf$sf)), ", creatinine on ", sum(!is.na(pf$creatinine)), " periods")
}
COMPETING_EVENT <- switch(config$cohort, imv = "extubation", niv = "intubation", nosupport = "escalation")
death_time <- base %>%
  transmute(hospitalization_id,
            death_time = if (JM_GRID == "daily") death_day else death_time_days)

# =============================================================================
# SYNTHETIC SITE ONLY: give the markers a patient-level structure
# =============================================================================
# The synthetic CLIF labs are independent draws per timestamp: within a patient
# the day-to-day autocorrelation of log creatinine and log SF is 0.02 to 0.03,
# and the fitted random-intercept variances are zero. A joint model is identified
# through those random effects, so on the raw synthetic data the sampler has
# nothing to move (random-effects acceptance 0.004 for SF) and every fit fails
# the R-hat gate for a reason that has nothing to do with the model. This guard,
# the marker analogue of the simulated survival in 10_panel_common.R, multiplies
# each marker by a patient-specific log-normal intercept and slope so the
# machinery can be exercised end to end. It changes nothing about which patients
# or days are present and never runs at a real site. Requested 2026-09-13.
if (is_synthetic) {
  message("*** SYNTHETIC SITE: adding a patient-level random intercept and slope to every marker (plumbing only). ***")
  marker_sd <- c(creatinine = 0.6, platelets = 0.4, bilirubin = 0.6, sf = 0.2, dp = 0.15,
                 ne_equiv_peak = 0.8, map_aw = 0.15)
  set.seed(20260913)
  synth_re <- base %>% select(hospitalization_id) %>%
    bind_cols(map_dfc(names(marker_sd), function(m) {
      n <- nrow(base)
      tibble(!!paste0("u_", m) := rnorm(n, 0, marker_sd[[m]]),      # intercept, log scale
             !!paste0("v_", m) := rnorm(n, 0, marker_sd[[m]] / 8))  # slope per day, log scale
    }))
  perturb <- function(df, m) {
    df %>% left_join(synth_re %>% select(hospitalization_id, u = all_of(paste0("u_", m)),
                                         v = all_of(paste0("v_", m))), by = "hospitalization_id") %>%
      mutate(!!m := .data[[m]] * exp(u + v * vent_day)) %>% select(-u, -v)
  }
  for (m in setdiff(names(marker_sd), c("dp", "map_aw"))) pf <- perturb(pf, m)
  dpp <- perturb(dpp, "dp")
  # The oxygenation indices are computed from matched measurements, so they are
  # perturbed through their numerator: both move with the patient's mean airway
  # pressure perturbation.
  oxy <- oxy %>% mutate(vent_day = period * STEP, map_aw_raw = map_aw)
  oxy <- perturb(oxy, "map_aw") %>%
    mutate(ratio = if_else(is.finite(map_aw / map_aw_raw), map_aw / map_aw_raw, 1),
           oi = oi * ratio, osi = osi * ratio) %>%
    select(-vent_day, -map_aw_raw, -ratio)
}

# =============================================================================
# Non-positive marker values
# =============================================================================
# The markers are modelled on the log scale. The lab outlier thresholds (script
# 02) admit zero for creatinine, platelets and bilirubin, and a zero is charted
# at real sites (MIMIC has creatinine rows of 0), which is not a measurement and
# is -Inf on the log scale: nlme then fails with "NA/NaN/Inf in foreign function
# call". Such values are set to missing here and counted per marker in the
# summary table. NE-equivalent dose keeps its true zeros (the fit adds an offset).
LOG_MARKERS <- c("creatinine", "platelets", "bilirubin", "sf", "dp", "oi", "osi")
nonpositive_counts <- setNames(integer(length(LOG_MARKERS)), LOG_MARKERS)
for (m in LOG_MARKERS) {
  df <- if (m == "dp") dpp else if (m %in% c("oi", "osi")) oxy else pf
  bad <- !is.na(df[[m]]) & df[[m]] <= 0
  nonpositive_counts[[m]] <- sum(bad)
  df[[m]][bad] <- NA_real_
  if (m == "dp") dpp <- df else if (m %in% c("oi", "osi")) oxy <- df else pf <- df
}
if (any(nonpositive_counts > 0))
  message("Non-positive marker values set to missing (patient-periods): ",
          paste(sprintf("%s %d", names(nonpositive_counts), nonpositive_counts), collapse = ", "))

# =============================================================================
# 13a. RRT start day: the first CRRT record or dialysis procedure (script 01)
# =============================================================================
# Renal replacement of any kind, continuous or intermittent, ends the creatinine
# trajectory and is its competing event. A patient with end-stage renal disease is
# on renal replacement before the index: censored at day 0, so no creatinine
# trajectory at all (user, 2026-09-21).
crrt_available <- readRDS(file.path(output_dir, "crrt_available.rds"))
rrt_sources    <- readRDS(file.path(output_dir, "rrt_sources_available.rds"))
if (!crrt_available)
  message("*** This site has no crrt_therapy table: CRRT does not censor the creatinine trajectory. ***")
if (!rrt_sources[["dialysis"]])
  message("*** This site has no patient_procedures table: intermittent dialysis does not censor the creatinine trajectory. ***")
if (!rrt_sources[["esrd"]])
  message("*** This site has no hospital_diagnosis table: ESRD patients are not identified. ***")
esrd_ids <- read_parquet(file.path(output_dir, "cohort_esrd.parquet"))$hospitalization_id
rrt <- bind_rows(read_parquet(file.path(output_dir, "cohort_crrt.parquet")) %>% select(hospitalization_id, recorded_dttm),
                 read_parquet(file.path(output_dir, "cohort_dialysis.parquet")) %>% select(hospitalization_id, recorded_dttm)) %>%
  group_by(hospitalization_id) %>%
  summarise(rrt_start_dttm = min(recorded_dttm), .groups = "drop") %>%
  right_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  transmute(hospitalization_id,
            rrt_period = as.integer(floor(as.numeric(difftime(rrt_start_dttm, t0, units = "hours")) / STEP_H)),
            esrd = hospitalization_id %in% esrd_ids,
            rrt_period = if_else(esrd, -1L, rrt_period),       # ESRD: on RRT before the index
            rrt_day = rrt_period * STEP) %>%
  filter(!is.na(rrt_period))
message("RRT (CRRT or dialysis procedure, or ESRD): ", nrow(rrt), " patients; ",
        sum(rrt$esrd), " with ESRD (censored at day 0), ",
        sum(rrt$rrt_day < 0 & !rrt$esrd), " others already on RRT at the index, ",
        sum(rrt$rrt_day >= 0 & rrt$rrt_day <= JM_HORIZON), " start within the horizon")
rrt <- rrt %>% select(-esrd)

# =============================================================================
# 13b. Survival table: death vs extubation within JM_HORIZON, tie = death
# =============================================================================
# death_day (index-anchored, all-cause, NA past 28 days; synthetic site simulated)
# and imv_extub_day (last IMV day + 1) come from the shared panel. Censoring at
# the horizon otherwise. JMbayes2 needs strictly positive times, so an event on
# day 0 is placed at day 1 (the TTE's pmax(., 1) convention).
# Marker baselines are the FIRST OBSERVED daily value of each marker within the
# horizon (the same reduction as the trajectory: creatinine and bilirubin max,
# platelets min, SF worst, DP max, NE-equivalent peak), with the day it was
# observed recorded as {marker}_0_day; the trajectory for that marker starts the
# day after. Day 0 is the baseline for everyone with a day-0 value (SF, DP and
# the NE-equivalent dose always; creatinine for most). Labs are not drawn every
# day, so requiring a day-0 value cost a third of the creatinine cohort at MIMIC
# (3,427 of 5,382 patients with two or more observations); the summary table
# reports how many baselines fall on day 0 and how many later. The
# cross-sectional table's index-timepoint labs are not used: they are matched
# inside a narrow window, missing for most patients, and carry no bilirubin.
marker_cols <- c(creatinine = "creatinine", platelets = "platelets", bilirubin = "bilirubin",
                 sf = "sf", dp = "dp", ne_equiv_peak = "ne_equiv_peak", oi = "oi", osi = "osi")
baseline_names <- c(creatinine = "creatinine_0", platelets = "platelet_0", bilirubin = "bilirubin_0",
                    sf = "sf_0", dp = "dp_0", ne_equiv_peak = "ne_equiv_0", oi = "oi_0", osi = "osi_0")
with_dp <- pf %>% filter(period <= N_PERIODS) %>%
  left_join(dpp %>% select(hospitalization_id, period, dp), by = c("hospitalization_id", "period")) %>%
  left_join(oxy %>% select(hospitalization_id, period, oi, osi), by = c("hospitalization_id", "period"))
first_obs <- map(names(marker_cols), function(m) {
  with_dp %>% filter(!is.na(.data[[m]])) %>%
    group_by(hospitalization_id) %>% slice_min(period, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(hospitalization_id, !!baseline_names[[m]] := .data[[m]], !!paste0(baseline_names[[m]], "_day") := period)
}) %>% reduce(full_join, by = "hospitalization_id")
day0 <- pf %>%
  filter(period == 0L) %>%
  select(hospitalization_id, vt_ml_0 = vt_ml, vtpfvc_0 = vtpfvc) %>%
  left_join(first_obs, by = "hospitalization_id")
# Patient-level strain: the mean of the daily VT/PFVC over the observed course
# within the horizon. It is the BETWEEN-patient term of the within-between
# decomposition in the longitudinal submodel (dose level plus PBW/PFVC
# discordance, which within a patient is a constant); the WITHIN term is each
# day's deviation from it, the clinician's dose change rescaled by that constant.
pt_strain <- pf %>%
  filter(period <= N_PERIODS, !is.na(vtpfvc)) %>%
  group_by(hospitalization_id) %>%
  summarise(vtpfvc_pt_mean = mean(vtpfvc), vtpbw_pt_mean = mean(vt_ml / pbw),
            vtpfvc_pt_n = n(), .groups = "drop")
surv <- base %>%
  left_join(rrt, by = "hospitalization_id") %>%
  left_join(day0, by = "hospitalization_id") %>%
  left_join(pt_strain, by = "hospitalization_id") %>%
  left_join(death_time, by = "hospitalization_id") %>%
  left_join(extub_time, by = "hospitalization_id") %>%
  mutate(
    # hazard-model exposures on the paper's primary scale: the clinician's dose
    # (VT/PBW at the index) and the size term (log PFVC); VT/PFVC at the index is
    # kept for reference
    vtpbw_idx  = vt_ml_0 / pbw,
    log_pfvc   = log(pfvc_gli),
    log_pbw    = log(pbw),
    # log PBW/PFVC discordance, centred at the cohort median: the effect modifier
    # of the dose slope in the primary longitudinal model
    ldisc_c    = log(pbw / pfvc_gli) - median(log(pbw / pfvc_gli), na.rm = TRUE),
    # per-SD versions for the PFVC-level question (the paper reports PFVC per SD)
    log_pfvc_sd = as.numeric(scale(log_pfvc)),
    ldisc_sd    = as.numeric(scale(log(pbw / pfvc_gli))),
    # VT/PFVC as the reader meets it: the patient's mean VT/PFVC over the window in
    # percent of predicted FVC (the project's unit; the 11% ceiling), centred at the
    # cohort median, per point; the index value for the hazard on the same scale.
    # At a given VT/PBW this is the PBW/PFVC discordance contrast scaled by the dose.
    vtpfvc_c   = vtpfvc_pt_mean - median(vtpfvc_pt_mean, na.rm = TRUE),
    vtpfvc_idx = vtpfvc_0,
    death_in  = !is.na(death_time) & death_time <= JM_HORIZON,
    extub_in  = !is.na(extub_time) & extub_time <= JM_HORIZON,
    event = case_when(
      death_in & (!extub_in | death_time <= extub_time) ~ 1L,   # death (tie counts as death)
      extub_in                                          ~ 2L,   # extubation
      TRUE                                              ~ 0L),  # censored at the horizon
    event_day = case_when(event == 1L ~ death_time,
                          event == 2L ~ extub_time,
                          TRUE        ~ as.numeric(JM_HORIZON)),
    event_time = pmax(event_day, STEP),   # JMbayes2 needs strictly positive times
    event_factor = factor(c("censored", "death", COMPETING_EVENT)[event + 1L],
                          levels = c("censored", "death", COMPETING_EVENT)),
    # A THIRD cause for the creatinine model: renal replacement. Dialysis does not
    # end the patient's course, but it ends the creatinine trajectory, and it is
    # started BECAUSE the creatinine is rising, so dropping those patient-days (as
    # the panel does) censors informatively on the very signal being measured. A
    # shared-parameter joint model corrects that only for events it models, so
    # these columns let RRT enter as a cause and share the random effects. Death
    # wins a tie with RRT, RRT wins a tie with extubation, because the question is
    # which one ends the marker series. Used only when PBWPFVC_JM_RRT_EVENT=1, and
    # only for creatinine: with RRT in, "death" means death before dialysis, so
    # every other marker's hazards would silently change meaning.
    t_death = if_else(death_in, death_time, Inf),
    t_extub = if_else(extub_in, extub_time, Inf),
    t_rrt   = if_else(!is.na(rrt_day) & rrt_day >= 0 & rrt_day <= JM_HORIZON, as.numeric(rrt_day), Inf),
    t_first = pmin(t_death, t_rrt, t_extub),
    event_rrt = case_when(!is.finite(t_first) ~ 0L, t_death <= t_first ~ 1L,
                          t_rrt <= t_first ~ 3L, TRUE ~ 2L),
    event_day_rrt  = if_else(is.finite(t_first), t_first, as.numeric(JM_HORIZON)),
    event_time_rrt = pmax(event_day_rrt, STEP),
    event_factor_rrt = factor(c("censored", "death", COMPETING_EVENT, "rrt")[event_rrt + 1L],
                              levels = c("censored", "death", COMPETING_EVENT, "rrt")),
    ers_pfvc_0 = ers * pfvc_gli,                 # specific elastance at the index (plateau subset)
    disc       = pbw / pfvc_gli,                 # PBW/PFVC discordance
    rrt_before_index = !is.na(rrt_day) & rrt_day < 0,
    creatinine_0 = if_else(rrt_before_index, NA_real_, creatinine_0)   # no creatinine trajectory on CRRT at the index
  ) %>%
  select(hospitalization_id, t0, event, event_day, event_time, event_factor,
         event_rrt, event_day_rrt, event_time_rrt, event_factor_rrt,
         death_day, imv_extub_day, death_time, extub_time, rrt_day, rrt_period, rrt_before_index,
         pfvc_gli, pfvc_age25, pbw, disc, disc_grp, age_grp, height_grp,
         age10, sex_category, race_category, sofa_total, np_sofa, sofa_cv_97, sofa_coag, sofa_liver, sofa_renal, bmi, height_cm,
         vtpbw_idx, log_pfvc, log_pbw, ldisc_c, log_pfvc_sd, ldisc_sd, vtpfvc_c, vtpfvc_idx,
         vtpfvc_0, vtpfvc_pt_mean, vtpbw_pt_mean, vtpfvc_pt_n,
         ers, ers_pfvc_0, creatinine_0, platelet_0, bilirubin_0, sf_0, dp_0, ne_equiv_0, oi_0, osi_0,
         ends_with("_0_day"))
# channel pieces of log PFVC (20_biotrauma_grid.R): the size term of the "channels" joint-model form
surv <- bind_cols(surv, pfvc_channels(surv, "log_pfvc"))
message("Survival table: ", nrow(surv), " patients; deaths ", sum(surv$event == 1L),
        ", ", COMPETING_EVENT, "s ", sum(surv$event == 2L), ", censored ", sum(surv$event == 0L))

# =============================================================================
# 13c. Longitudinal table: markers by day with the previous day's exposure
# =============================================================================
# One row per patient-day with a set tidal volume (the `daily` grid), day 0 to
# JM_HORIZON. The previous day's values are attached by joining on vent_day - 1,
# so a gap in charting yields a missing lag (reported below) instead of a stale one.
prev <- pf %>%
  transmute(hospitalization_id, period = period + 1L,
            l_vtpfvc = vtpfvc, l_vtpbw = vt_ml / pbw,
            l_sf = sf, l_pressor = on_pressor, l_fio2 = fio2, l_peep = peep)
# Cumulative strain through the previous day, two forms. mean_prior_vtpfvc (the
# fit's default) is the mean of the daily VT/PFVC over days 0..t-1: it carries the
# dose history without growing with time. cum_days_above (days above 11%) is the
# sensitivity form: for a patient above the ceiling throughout it equals the
# ventilator day exactly, so it is an interaction of day with a patient indicator
# and competes with the day spline and the random slope (R-hat 3 to 4 on it in
# every fit, synthetic and MIMIC).
cum_above <- pf %>%
  group_by(hospitalization_id) %>% arrange(period, .by_group = TRUE) %>%
  transmute(hospitalization_id, period = period + 1L,
            cum_days_above = cumsum(vtpfvc > STRAIN_CEILING) * STEP,   # days above, through the previous period
            mean_prior_vtpfvc = cummean(vtpfvc)) %>%
  ungroup()
long <- pf %>%
  filter(period <= N_PERIODS) %>%
  select(hospitalization_id, period, vent_day, vtpfvc, vt_ml, fio2, peep, rr, map, sf, on_pressor,
         ne_equiv_peak, creatinine, platelets, bilirubin) %>%
  left_join(dpp %>% select(hospitalization_id, period, dp), by = c("hospitalization_id", "period")) %>%
  left_join(oxy %>% select(hospitalization_id, period, map_aw, oi, osi), by = c("hospitalization_id", "period")) %>%
  left_join(prev, by = c("hospitalization_id", "period")) %>%
  left_join(cum_above, by = c("hospitalization_id", "period")) %>%
  mutate(cum_days_above = if_else(period == 0L, 0, cum_days_above)) %>%   # mean_prior_vtpfvc stays NA at period 0
  inner_join(surv %>% select(hospitalization_id, event_day, rrt_period, rrt_before_index,
                             vtpfvc_pt_mean, vtpbw_pt_mean, ends_with("_0_day")),
             by = "hospitalization_id") %>%
  filter(vent_day <= event_day) %>%
  mutate(
    # each marker's trajectory starts the period after its baseline observation
    creatinine    = if_else(!is.na(creatinine_0_day) & period <= creatinine_0_day, NA_real_, creatinine),
    platelets     = if_else(!is.na(platelet_0_day)   & period <= platelet_0_day,   NA_real_, platelets),
    bilirubin     = if_else(!is.na(bilirubin_0_day)  & period <= bilirubin_0_day,  NA_real_, bilirubin),
    sf            = if_else(!is.na(sf_0_day)         & period <= sf_0_day,         NA_real_, sf),
    dp            = if_else(!is.na(dp_0_day)         & period <= dp_0_day,         NA_real_, dp),
    ne_equiv_peak = if_else(!is.na(ne_equiv_0_day)   & period <= ne_equiv_0_day,   NA_real_, ne_equiv_peak),
    oi            = if_else(!is.na(oi_0_day)         & period <= oi_0_day,         NA_real_, oi),
    osi           = if_else(!is.na(osi_0_day)        & period <= osi_0_day,        NA_real_, osi),
    # within-patient strain: yesterday's VT/PFVC relative to the patient's own mean
    l_vtpfvc_within = l_vtpfvc - vtpfvc_pt_mean,
    l_vtpbw_within  = l_vtpbw  - vtpbw_pt_mean,    # the clinician's dose change (mL/kg PBW)
    # creatinine censored at RRT start; no trajectory if on CRRT at the index
    creat_censored_rrt = rrt_before_index | (!is.na(rrt_period) & period >= rrt_period),
    creatinine = if_else(creat_censored_rrt, NA_real_, creatinine)
  ) %>%
  select(-event_day, -rrt_period, -rrt_before_index, -vtpfvc_pt_mean, -vtpbw_pt_mean, -ends_with("_0_day")) %>%
  arrange(hospitalization_id, period)
message("Longitudinal table: ", nrow(long), " patient-periods, ",
        n_distinct(long$hospitalization_id), " patients; lag missing on ",
        sum(is.na(long$l_vtpfvc) & long$period > 0L), " post-index rows; creatinine periods removed for RRT: ",
        sum(long$creat_censored_rrt))

# =============================================================================
# 13d. Aggregate summary (deliverable) and persistence
# =============================================================================
markers <- c("creatinine", "platelets", "bilirubin", "sf", "dp", "ne_equiv_peak", "oi", "osi")
per_marker <- map_dfr(markers, function(m) {
  obs <- long %>% filter(!is.na(.data[[m]]))
  per_pt <- obs %>% count(hospitalization_id)
  ids2 <- per_pt$hospitalization_id[per_pt$n >= 2L]
  ev <- surv %>% filter(hospitalization_id %in% ids2)
  y0 <- baseline_names[[m]]
  tibble(marker = m,
         patient_days = nrow(obs),
         nonpositive_set_missing = if (m %in% LOG_MARKERS) nonpositive_counts[[m]] else 0L,
         patients_any = nrow(per_pt),
         patients_with_baseline = sum(!is.na(surv[[y0]])),
         patients_day0_baseline = sum(surv[[paste0(y0, "_day")]] == 0L, na.rm = TRUE),
         patients_ge2_obs = length(ids2),
         median_obs_per_patient = if (nrow(per_pt)) median(per_pt$n) else NA_real_,
         deaths_ge2 = sum(ev$event == 1L), extubations_ge2 = sum(ev$event == 2L),
         plateau_subset_ge2 = sum(!is.na(ev$ers_pfvc_0)))
})
summary_tbl <- bind_rows(
  per_marker,
  tibble(marker = "cohort",
         patient_days = nrow(long), nonpositive_set_missing = NA_integer_, patients_any = nrow(surv),
         patients_with_baseline = NA_integer_, patients_day0_baseline = NA_integer_,
         patients_ge2_obs = NA_integer_, median_obs_per_patient = NA_real_,
         deaths_ge2 = sum(surv$event == 1L), extubations_ge2 = sum(surv$event == 2L),
         plateau_subset_ge2 = sum(!is.na(surv$ers_pfvc_0)))) %>%
  mutate(grid = JM_GRID, step_hours = STEP_H, horizon_days = JM_HORIZON,
         crrt_available = crrt_available, dialysis_procedures_available = rrt_sources[["dialysis"]],
         esrd_diagnoses_available = rrt_sources[["esrd"]], esrd_censored_day0 = sum(surv$hospitalization_id %in% esrd_ids),
         rrt_before_index = sum(surv$rrt_before_index),
         rrt_within_horizon = sum(!is.na(surv$rrt_day) & surv$rrt_day >= 0 & surv$rrt_day <= JM_HORIZON),
         creatinine_days_removed_rrt = sum(long$creat_censored_rrt),
         lag_missing_rows = sum(is.na(long$l_vtpfvc) & long$period > 0L),
         site = site_name)
print(as.data.frame(summary_tbl), row.names = FALSE)
write_csv(mask_small_counts(summary_tbl), file.path(final_dir, paste0("jm_panel_summary_", h_suffix, "_", site_name, ".csv")))

write_parquet(long, file.path(output_dir, paste0("jm_long_", h_suffix, ".parquet")))
write_parquet(surv, file.path(output_dir, paste0("jm_surv_", h_suffix, ".parquet")))
saveRDS(list(grid = JM_GRID, step_hours = STEP_H, step = STEP, horizon = JM_HORIZON, n_periods = N_PERIODS,
             h_suffix = h_suffix, strain_ceiling = STRAIN_CEILING,
             panel_cohort_tag = panel_cohort_tag, site_name = site_name,
             n_patients = nrow(surv), n_days = nrow(long), built_at = as.character(Sys.time())),
        file.path(output_dir, paste0("jm_meta_", h_suffix, ".rds")))
message("21_biotrauma_panel complete (", h_suffix, "): tables in ", output_dir,
        "; summary in ", final_dir)
