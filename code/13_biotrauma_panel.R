# =============================================================================
# Script 13 (panel): Biotrauma joint models -- the longitudinal and survival tables
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# Builds the two tables the biotrauma joint models (13_biotrauma_fit.R) consume,
# from the shared daily panel of 10_panel_common.R:
#
#   jm_long_{H}d.parquet   one row per patient-day (index day 0 to day H) with the
#                          organ-injury markers observed that day, the previous
#                          day's strain and confounders, and the patient's
#                          baseline covariates
#   jm_surv_{H}d.parquet   one row per patient: competing-risk coding of death vs
#                          extubation within H days (same-day tie counts as death),
#                          the RRT start day, and every baseline covariate
#   jm_meta_{H}d.rds       horizon, cohort tag, counts
#
# and one aggregate table for the deliverable:
#
#   final/jm_panel_summary_{site}.csv   patients, patient-days and events per
#                                       marker; RRT censoring counts; the size of
#                                       the plateau-measured subset
#
# Design (docs/joint_model_plan_2026-09.md, sections 3 and 4):
#   * every index-IMV patient enters at day 0; no survival-based restriction, and
#     no structural-positivity exclusion (that is TTE-only)
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
# Usage: Rscript code/13_biotrauma_panel.R
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
output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

# --- shared-panel contract (identical windows to the TTE, so the event objects match)
HORIZON      <- 28L
MAX_VENT_DAY <- 27L
is_synthetic <- identical(site_name, "synthetic_clif")
PANEL_NORM   <- "pfvc"          # the joint models always normalize to GLI PFVC
source(here("code", "10_panel_common.R"))

JM_HORIZON <- as.integer(Sys.getenv("PBWPFVC_JM_HORIZON", "7"))
stopifnot(is.finite(JM_HORIZON), JM_HORIZON >= 2L, JM_HORIZON <= HORIZON)
h_suffix <- paste0(JM_HORIZON, "d")
STRAIN_CEILING <- 11   # VT/PFVC % above which a day counts toward the cumulative-strain exposure
message("=== 13_biotrauma_panel: horizon ", JM_HORIZON, " days, site ", site_name, " ===")

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
  marker_sd <- c(creatinine = 0.6, platelets = 0.4, bilirubin = 0.6, sf = 0.2, dp = 0.15, ne_equiv_peak = 0.8)
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
  for (m in setdiff(names(marker_sd), "dp")) panel_full <- perturb(panel_full, m)
  dp_daily <- perturb(dp_daily, "dp")
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
LOG_MARKERS <- c("creatinine", "platelets", "bilirubin", "sf", "dp")
nonpositive_counts <- setNames(integer(length(LOG_MARKERS)), LOG_MARKERS)
for (m in LOG_MARKERS) {
  df <- if (m == "dp") dp_daily else panel_full
  bad <- !is.na(df[[m]]) & df[[m]] <= 0
  nonpositive_counts[[m]] <- sum(bad)
  df[[m]][bad] <- NA_real_
  if (m == "dp") dp_daily <- df else panel_full <- df
}
if (any(nonpositive_counts > 0))
  message("Non-positive marker values set to missing (patient-days): ",
          paste(sprintf("%s %d", names(nonpositive_counts), nonpositive_counts), collapse = ", "))

# =============================================================================
# 13a. RRT start day (CRRT table from script 01)
# =============================================================================
crrt_available <- readRDS(file.path(output_dir, "crrt_available.rds"))
if (!crrt_available)
  message("*** This site has no crrt_therapy table: the creatinine trajectory is NOT censored at RRT start. ***")
rrt <- read_parquet(file.path(output_dir, "cohort_crrt.parquet")) %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  group_by(hospitalization_id) %>%
  summarise(rrt_start_dttm = min(recorded_dttm), .groups = "drop") %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  transmute(hospitalization_id,
            rrt_day = floor(as.numeric(difftime(rrt_start_dttm, t0, units = "days"))))
message("CRRT: ", nrow(rrt), " patients with any record; ",
        sum(rrt$rrt_day < 0), " already on CRRT at the index, ",
        sum(rrt$rrt_day >= 0 & rrt$rrt_day <= JM_HORIZON), " start within the horizon")

# =============================================================================
# 13b. Survival table: death vs extubation within JM_HORIZON, tie = death
# =============================================================================
# death_day (index-anchored, all-cause, NA past 28 days; synthetic site simulated)
# and imv_extub_day (last IMV day + 1) come from the shared panel. Censoring at
# the horizon otherwise. JMbayes2 needs strictly positive times, so an event on
# day 0 is placed at day 1 (the TTE's pmax(., 1) convention).
# Marker baselines are the DAY-0 daily values (the same reduction as the trajectory:
# creatinine and bilirubin max, platelets min, SF worst, DP max, NE-equivalent peak).
# The cross-sectional table's index-timepoint labs are matched inside a narrow
# window and are missing for most patients (synthetic: 218 of 676 for creatinine
# against 428 with a day-0 value), and it carries no bilirubin at all. A patient
# without a day-0 value has no baseline for that marker and leaves that marker's
# model; the count is in the summary table.
day0 <- panel_full %>%
  filter(vent_day == 0L) %>%
  left_join(dp_daily, by = c("hospitalization_id", "vent_day")) %>%
  select(hospitalization_id, creatinine_0 = creatinine, platelet_0 = platelets,
         bilirubin_0 = bilirubin, sf_0 = sf, dp_0 = dp, ne_equiv_0 = ne_equiv_peak,
         vt_ml_0 = vt_ml, vtpfvc_0 = vtpfvc)
# Patient-level strain: the mean of the daily VT/PFVC over the observed course
# within the horizon. It is the BETWEEN-patient term of the within-between
# decomposition in the longitudinal submodel (dose level plus PBW/PFVC
# discordance, which within a patient is a constant); the WITHIN term is each
# day's deviation from it, the clinician's dose change rescaled by that constant.
pt_strain <- panel_full %>%
  filter(vent_day <= JM_HORIZON, !is.na(vtpfvc)) %>%
  group_by(hospitalization_id) %>%
  summarise(vtpfvc_pt_mean = mean(vtpfvc), vtpbw_pt_mean = mean(vt_ml / pbw),
            vtpfvc_pt_n = n(), .groups = "drop")
surv <- base %>%
  left_join(rrt, by = "hospitalization_id") %>%
  left_join(day0, by = "hospitalization_id") %>%
  left_join(pt_strain, by = "hospitalization_id") %>%
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
    death_in  = !is.na(death_day) & death_day <= JM_HORIZON,
    extub_in  = !is.na(imv_extub_day) & imv_extub_day <= JM_HORIZON,
    event = case_when(
      death_in & (!extub_in | death_day <= imv_extub_day) ~ 1L,   # death (tie counts as death)
      extub_in                                            ~ 2L,   # extubation
      TRUE                                                ~ 0L),  # censored at the horizon
    event_day = case_when(event == 1L ~ death_day,
                          event == 2L ~ as.numeric(imv_extub_day),
                          TRUE        ~ as.numeric(JM_HORIZON)),
    event_time = pmax(event_day, 1),
    event_factor = factor(c("censored", "death", "extubation")[event + 1L],
                          levels = c("censored", "death", "extubation")),
    ers_pfvc_0 = ers * pfvc_gli,                 # specific elastance at the index (plateau subset)
    disc       = pbw / pfvc_gli,                 # PBW/PFVC discordance
    rrt_before_index = !is.na(rrt_day) & rrt_day < 0,
    creatinine_0 = if_else(rrt_before_index, NA_real_, creatinine_0)   # no creatinine trajectory on CRRT at the index
  ) %>%
  select(hospitalization_id, t0, event, event_day, event_time, event_factor,
         death_day, imv_extub_day, rrt_day, rrt_before_index,
         pfvc_gli, pfvc_age25, pbw, disc, disc_grp, age_grp, height_grp,
         age10, sex_category, race_category, sofa_total, np_sofa, bmi, height_cm,
         vtpbw_idx, log_pfvc, log_pbw, ldisc_c, vtpfvc_0, vtpfvc_pt_mean, vtpbw_pt_mean, vtpfvc_pt_n,
         ers, ers_pfvc_0, creatinine_0, platelet_0, bilirubin_0, sf_0, dp_0, ne_equiv_0)
message("Survival table: ", nrow(surv), " patients; deaths ", sum(surv$event == 1L),
        ", extubations ", sum(surv$event == 2L), ", censored ", sum(surv$event == 0L))

# =============================================================================
# 13c. Longitudinal table: markers by day with the previous day's exposure
# =============================================================================
# One row per patient-day with a set tidal volume (the `daily` grid), day 0 to
# JM_HORIZON. The previous day's values are attached by joining on vent_day - 1,
# so a gap in charting yields a missing lag (reported below) instead of a stale one.
prev <- panel_full %>%
  transmute(hospitalization_id, vent_day = vent_day + 1L,
            l_vtpfvc = vtpfvc, l_vtpbw = vt_ml / pbw,
            l_sf = sf, l_pressor = on_pressor, l_fio2 = fio2, l_peep = peep)
# Cumulative strain through the previous day, two forms. mean_prior_vtpfvc (the
# fit's default) is the mean of the daily VT/PFVC over days 0..t-1: it carries the
# dose history without growing with time. cum_days_above (days above 11%) is the
# sensitivity form: for a patient above the ceiling throughout it equals the
# ventilator day exactly, so it is an interaction of day with a patient indicator
# and competes with the day spline and the random slope (R-hat 3 to 4 on it in
# every fit, synthetic and MIMIC).
cum_above <- panel_full %>%
  group_by(hospitalization_id) %>% arrange(vent_day, .by_group = TRUE) %>%
  transmute(hospitalization_id, vent_day = vent_day + 1L,
            cum_days_above = cumsum(vtpfvc > STRAIN_CEILING),
            mean_prior_vtpfvc = cummean(vtpfvc)) %>%
  ungroup()
long <- panel_full %>%
  filter(vent_day <= JM_HORIZON) %>%
  select(hospitalization_id, vent_day, vtpfvc, vt_ml, fio2, peep, rr, map, sf, on_pressor,
         ne_equiv_peak, creatinine, platelets, bilirubin) %>%
  left_join(dp_daily, by = c("hospitalization_id", "vent_day")) %>%
  left_join(prev, by = c("hospitalization_id", "vent_day")) %>%
  left_join(cum_above, by = c("hospitalization_id", "vent_day")) %>%
  mutate(cum_days_above = if_else(vent_day == 0L, 0L, cum_days_above)) %>%   # mean_prior_vtpfvc stays NA on day 0
  inner_join(surv %>% select(hospitalization_id, event_day, rrt_day, rrt_before_index,
                             vtpfvc_pt_mean, vtpbw_pt_mean),
             by = "hospitalization_id") %>%
  filter(vent_day <= event_day) %>%
  mutate(
    # within-patient strain: yesterday's VT/PFVC relative to the patient's own mean
    l_vtpfvc_within = l_vtpfvc - vtpfvc_pt_mean,
    l_vtpbw_within  = l_vtpbw  - vtpbw_pt_mean,    # the clinician's dose change (mL/kg PBW)
    # creatinine censored at RRT start; no trajectory if on CRRT at the index
    creat_censored_rrt = rrt_before_index | (!is.na(rrt_day) & vent_day >= rrt_day),
    creatinine = if_else(creat_censored_rrt, NA_real_, creatinine)
  ) %>%
  select(-event_day, -rrt_day, -rrt_before_index, -vtpfvc_pt_mean, -vtpbw_pt_mean) %>%
  arrange(hospitalization_id, vent_day)
message("Longitudinal table: ", nrow(long), " patient-days, ",
        n_distinct(long$hospitalization_id), " patients; lag missing on ",
        sum(is.na(long$l_vtpfvc) & long$vent_day > 0L), " post-index rows; creatinine days removed for RRT: ",
        sum(long$creat_censored_rrt))

# =============================================================================
# 13d. Aggregate summary (deliverable) and persistence
# =============================================================================
markers <- c("creatinine", "platelets", "bilirubin", "sf", "dp", "ne_equiv_peak")
per_marker <- map_dfr(markers, function(m) {
  obs <- long %>% filter(!is.na(.data[[m]]))
  per_pt <- obs %>% count(hospitalization_id)
  ids2 <- per_pt$hospitalization_id[per_pt$n >= 2L]
  ev <- surv %>% filter(hospitalization_id %in% ids2)
  y0 <- c(creatinine = "creatinine_0", platelets = "platelet_0", bilirubin = "bilirubin_0",
          sf = "sf_0", dp = "dp_0", ne_equiv_peak = "ne_equiv_0")[[m]]
  tibble(marker = m,
         patient_days = nrow(obs),
         nonpositive_set_missing = if (m %in% LOG_MARKERS) nonpositive_counts[[m]] else 0L,
         patients_any = nrow(per_pt),
         patients_day0_baseline = sum(!is.na(surv[[y0]])),
         patients_ge2_obs = length(ids2),
         median_obs_per_patient = if (nrow(per_pt)) median(per_pt$n) else NA_real_,
         deaths_ge2 = sum(ev$event == 1L), extubations_ge2 = sum(ev$event == 2L),
         plateau_subset_ge2 = sum(!is.na(ev$ers_pfvc_0)))
})
summary_tbl <- bind_rows(
  per_marker,
  tibble(marker = "cohort",
         patient_days = nrow(long), nonpositive_set_missing = NA_integer_, patients_any = nrow(surv),
         patients_day0_baseline = NA_integer_,
         patients_ge2_obs = NA_integer_, median_obs_per_patient = NA_real_,
         deaths_ge2 = sum(surv$event == 1L), extubations_ge2 = sum(surv$event == 2L),
         plateau_subset_ge2 = sum(!is.na(surv$ers_pfvc_0)))) %>%
  mutate(horizon_days = JM_HORIZON,
         crrt_available = crrt_available,
         rrt_before_index = sum(surv$rrt_before_index),
         rrt_within_horizon = sum(!is.na(surv$rrt_day) & surv$rrt_day >= 0 & surv$rrt_day <= JM_HORIZON),
         creatinine_days_removed_rrt = sum(long$creat_censored_rrt),
         lag_missing_rows = sum(is.na(long$l_vtpfvc) & long$vent_day > 0L),
         site = site_name)
print(as.data.frame(summary_tbl), row.names = FALSE)
write_csv(summary_tbl, file.path(final_dir, paste0("jm_panel_summary_", h_suffix, "_", site_name, ".csv")))

write_parquet(long, file.path(output_dir, paste0("jm_long_", h_suffix, ".parquet")))
write_parquet(surv, file.path(output_dir, paste0("jm_surv_", h_suffix, ".parquet")))
saveRDS(list(horizon = JM_HORIZON, h_suffix = h_suffix, strain_ceiling = STRAIN_CEILING,
             panel_cohort_tag = panel_cohort_tag, site_name = site_name,
             n_patients = nrow(surv), n_days = nrow(long), built_at = as.character(Sys.time())),
        file.path(output_dir, paste0("jm_meta_", h_suffix, ".rds")))
message("13_biotrauma_panel complete (", h_suffix, "): tables in ", output_dir,
        "; summary in ", final_dir)
