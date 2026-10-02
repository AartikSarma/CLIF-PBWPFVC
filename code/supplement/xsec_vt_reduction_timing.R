# =============================================================================
# Supplement (cross-sectional): when clinicians turn the tidal volume down, is it
# because a small lung got sicker?
# =============================================================================
# xsec_vtpbw_gate_collider.R found that the VT/PBW gate and the VT/PBW term make a
# small predicted lung look more lethal: among patients at the same VT/PBW, the
# small-PFVC ones carry more of the illness that leads clinicians to lower the volume.
# That reading treats the illness as independent of lung size. The alternative is
# that the small lung CAUSED it: a smaller lung takes more strain at the same VT/PBW,
# deteriorates, and is turned down for that reason. This script asks which arrow
# from predicted size into the set tidal volume carries the reductions. Three can:
#
#   (i)   the starting setting. Round-number volumes put a short patient high in
#         mL/kg PBW, and a high start is corrected down. Dosing habit, not illness.
#   (ii)  mechanics. At the same VT/PBW a small lung has a higher driving pressure,
#         and a high driving pressure prompts a lower volume. Size acting through
#         the pressure it produces, not illness.
#   (iii) deterioration. Compliance falls over the first days (lung injury, possibly
#         ventilator-induced), and the volume is lowered in response. Only this arrow
#         is "small lungs cause the illness".
#
# Stated before the data:
#   collider reading     size's link to reductions runs through (i) and (ii); size
#                        does not predict a fall in compliance; small-PFVC reductions
#                        are no more often preceded by one
#   size-caused illness  size predicts a fall in compliance (Link 1 below), and
#                        small-PFVC reductions are preceded by one more often
#
# Population: every patient at the ungated index (03, 3e3), on the clock from the
# first IMV row, over the first HORIZON_H hours. A reduction is a fall of at least
# REDUCTION_MIN_ML in the set tidal volume between consecutive rows that are both in
# a volume-targeted mode. The threshold is in mL, not mL/kg PBW: a mL/kg threshold
# would make a given change in volume count more often in small patients and build
# the result in. A volume set in a volume-targeted mode and carried into a pressure
# mode is not a setting, so a switch of mode never counts as a reduction.
#
# Triggers are read only at measurement times, never carried forward as a value.
# Plateau pressure is charted only when someone measures it, so for each moment
# (a reduction, or the start of an hour at risk) the script asks whether a plateau
# was measured in the TRIGGER_WINDOWS_H hours before it, and if so:
#   DP high      driving pressure above DP_HIGH (arrow ii)
#   Crs fall     log compliance at least CRS_FALL below the patient's first measured
#                compliance (arrow iii); "first measurement" when the trigger
#                measurement is itself the first
# "No plateau measured" is its own category and an informative one: a reduction with
# no pressure measured before it is a protocol or practice correction. An SF fall
# (SF_FALL below the patient's baseline) is reported beside it, but early SF is not
# a clean marker of injury: more volume per unit of lung recruits more and raises SF
# in the first days, so its early sign can run opposite to injury.
#
# 1. Reductions by size (no model). Patients by quartile of log PFVC and by where
#    their first volume-targeted VT/PBW sat (below, in or above 6-8): how many had a
#    reduction, how soon, and what preceded the first one. The "above 8" rows of the
#    paper's cohort are the patients the gate indexed only once a reduction brought
#    them into the band.
# 2. Decomposition (discrete-time hazard of the FIRST reduction, hourly). The size
#    term's coefficient alone, then with the current VT/PBW (arrow i), then with the
#    DP trigger (arrow ii), then with the Crs and SF triggers (arrow iii). How far it
#    moves at each step says which arrow carries it. This conditions on the current
#    VT/PBW, the collider of the companion script, so it is a descriptive
#    decomposition of who gets turned down, not a causal estimate.
# 3. Link 1 (does size predict deterioration at all). Fixed-horizon ANCOVA: log Crs
#    over FOLLOWUP_H hours on the size term and log Crs in the first BASELINE_H
#    hours; SF the same way. Fitted without the initial VT/PBW (primary: the total
#    effect of size, which the companion script found is the identifiable one) and
#    with it (the standing "at a given VT/PBW" framing). The two answer different
#    questions, so both are reported. Patients extubated or dead before the
#    follow-up window drop out of it, and the table counts them.
#
# Exposures: log PFVC and the PBW/PFVC ratio, per SD of the ungated sample (as in
# the companion script); unadjusted and adjusted for ns(age10, 4), sex and race.
#
# Thresholds are fixed here, with one sensitivity (a 2-hour trigger window, a
# response to a pressure just measured):
#   REDUCTION_MIN_ML 25   the smallest step clinicians use (0.5 mL/kg on a 50 kg PBW)
#   TRIGGER_WINDOWS_H 6, 2
#   DP_HIGH 15 cmH2O      the driving pressure above which mortality rises (Amato 2015)
#   CRS_FALL 10%          a fall in compliance beyond charting noise
#   SF_FALL 25            a quarter of the way from SF 315 (hypoxemia) to 215
#
# Writes to final/supplement/:
#   vt_reduction_timing_events_{site}.csv     section 1, per trigger window
#   vt_reduction_timing_hazard_{site}.csv     section 2: the size term at each step,
#                                             and the triggers in the full model
#   vt_reduction_timing_ancova_{site}.csv     section 3
#   vt_reduction_timing_{site}.pdf            the size term across the decomposition
#
# Needs scripts 01-03 (the ventilated cohort; 03 writes analysis_ungated_index).
# Usage: uvr run code/supplement/xsec_vt_reduction_timing.R
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(splines); library(data.table) })
rm(list = ls())
source("utils/config.R")
if (config$cohort != "imv") stop("Tidal-volume reductions exist only in the ventilated cohort; unset PBWPFVC_COHORT")
site_name  <- config$site_name
output_dir <- config$output_dir
final_dir  <- final_dir_for("supplement")

HORIZON_H             <- 72
VOLUME_TARGETED_MODES <- c("assist control-volume control", "pressure-regulated volume control", "simv")
REDUCTION_MIN_ML      <- 25
TRIGGER_WINDOWS_H     <- c(6, 2)    # primary first
DP_HIGH               <- 15
CRS_FALL              <- 0.10
SF_FALL               <- 25
BASELINE_H            <- 6
FOLLOWUP_H            <- c(24, 48)
VTPBW_BAND            <- c(6, 8)

# --- patients: the ungated index, one row each
patients <- read_parquet(file.path(output_dir, "analysis_ungated_index.parquet"),
                         col_select = c("hospitalization_id", "pbw", "pfvc", "pbwpfvc", "age_at_admission",
                                        "sex_category", "race_category", "in_paper_cohort")) %>%
  mutate(sex_category  = factor(sex_category,  levels = c("Male", "Female")),
         race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
         age10 = age_at_admission / 10,
         log_pfvc = log(pfvc),
         pfvc_quartile = paste0("Q", ntile(log_pfvc, 4)))
EXPOSURES   <- c("PBW/PFVC" = "pbwpfvc", "log PFVC" = "log_pfvc")
ADJUSTMENTS <- c(unadjusted = "", adjusted = "+ ns(age10, 4) + sex_category + race_category")
exposure_sd <- map_dbl(EXPOSURES, ~ sd(patients[[.x]]))
for (exposure in EXPOSURES) patients[[paste0(exposure, "_sd")]] <- patients[[exposure]] / sd(patients[[exposure]])

# --- every ventilator row of those patients over the first HORIZON_H hours of IMV
timepoints <- open_dataset(file.path(output_dir, "analysis_all_eligible_timepoints.parquet")) %>%
  filter(hospitalization_id %in% patients$hospitalization_id) %>%
  select(hospitalization_id, recorded_dttm, mode_category, tidal_volume_set, dp, crs, sf_ratio) %>%
  collect() %>%
  as.data.table()
setorder(timepoints, hospitalization_id, recorded_dttm)
timepoints[, t_h := as.numeric(difftime(recorded_dttm, min(recorded_dttm), units = "hours")), by = hospitalization_id]
timepoints <- timepoints[t_h < HORIZON_H]
timepoints[, volume_mode := tolower(mode_category) %in% VOLUME_TARGETED_MODES & !is.na(tidal_volume_set)]
timepoints[, `:=`(previous_vt = shift(tidal_volume_set), previous_volume_mode = shift(volume_mode)), by = hospitalization_id]
timepoints[, reduction := volume_mode & previous_volume_mode %in% TRUE &
             (previous_vt - tidal_volume_set) >= REDUCTION_MIN_ML]
timepoints <- merge(timepoints, as.data.table(patients)[, .(hospitalization_id, pbw)], by = "hospitalization_id")
message("=== xsec_vt_reduction_timing: ", site_name, ", ", nrow(patients), " patients, ", nrow(timepoints),
        " ventilator rows in the first ", HORIZON_H, " h, ", sum(timepoints$reduction), " reductions of >= ",
        REDUCTION_MIN_ML, " mL")

# --- measurement events: plateau (DP and compliance) and SF, with each patient's baseline
plateau_events <- timepoints[!is.na(dp) & !is.na(crs) & crs > 0, .(hospitalization_id, t_h, dp, log_crs = log(crs))]
plateau_events[, first_log_crs := log_crs[1], by = hospitalization_id]
plateau_events[, is_first := seq_len(.N) == 1L, by = hospitalization_id]
sf_events <- timepoints[!is.na(sf_ratio), .(hospitalization_id, t_h, sf_ratio)]
sf_baseline <- sf_events[t_h < BASELINE_H, .(baseline_sf = median(sf_ratio)), by = hospitalization_id]

# the triggers at a set of moments (hospitalization_id, t_h): the most recent plateau
# and SF measured in the window_h hours strictly before each moment
triggers_at <- function(moments, window_h) {
  query <- copy(moments)[, query_t := t_h - 1e-6]
  setkey(plateau_events, hospitalization_id, t_h)
  last_plateau <- plateau_events[query, on = .(hospitalization_id, t_h = query_t), roll = window_h,
                                 .(dp, log_crs, first_log_crs, is_first)]
  setkey(sf_events, hospitalization_id, t_h)
  last_sf <- sf_events[query, on = .(hospitalization_id, t_h = query_t), roll = window_h, .(sf_ratio)]
  moments[, `:=`(
    dp_trigger = factor(fcase(is.na(last_plateau$dp), "no plateau measured",
                              last_plateau$dp > DP_HIGH, "DP high",
                              default = "DP not high"),
                        levels = c("no plateau measured", "DP not high", "DP high")),
    crs_trigger = factor(fcase(is.na(last_plateau$dp), "no plateau measured",
                               last_plateau$is_first, "first measurement",
                               last_plateau$log_crs - last_plateau$first_log_crs <= log(1 - CRS_FALL), "Crs fall",
                               default = "no Crs fall"),
                         levels = c("no plateau measured", "first measurement", "no Crs fall", "Crs fall")))]
  # For the hazard model: whether a plateau was measured is already in dp_trigger, so
  # the compliance factor contrasts only among measurements (otherwise the two
  # factors share their "no plateau measured" column)
  moments[, crs_change := factor(fcase(crs_trigger == "first measurement", "first measurement",
                                       crs_trigger == "Crs fall", "Crs fall",
                                       default = "no Crs fall, or no plateau"),
                                 levels = c("no Crs fall, or no plateau", "first measurement", "Crs fall"))]
  baseline_sf <- sf_baseline$baseline_sf[match(moments$hospitalization_id, sf_baseline$hospitalization_id)]
  moments[, sf_trigger := factor(fcase(is.na(last_sf$sf_ratio) | is.na(baseline_sf), "SF not measured",
                                       baseline_sf - last_sf$sf_ratio >= SF_FALL, "SF fall",
                                       default = "no SF fall"),
                                 levels = c("SF not measured", "no SF fall", "SF fall"))]
  moments
}

# --- each patient's first volume-targeted setting and first reduction
first_setting <- timepoints[volume_mode == TRUE, .(initial_vtpbw = tidal_volume_set[1] / pbw[1]), by = hospitalization_id]
first_reduction <- timepoints[reduction == TRUE, .(reduction_t_h = t_h[1]), by = hospitalization_id]
patient_course <- as.data.table(patients) %>%
  merge(first_setting, by = "hospitalization_id", all.x = TRUE) %>%
  merge(first_reduction, by = "hospitalization_id", all.x = TRUE)
patient_course[, initial_band := factor(fcase(is.na(initial_vtpbw), "never volume-targeted",
                                              initial_vtpbw < VTPBW_BAND[1], "below 6",
                                              initial_vtpbw > VTPBW_BAND[2], "above 8",
                                              default = "6-8"),
                                        levels = c("below 6", "6-8", "above 8", "never volume-targeted"))]

# =============================================================================
# 1. Reductions by size, and what preceded the first one
# =============================================================================
events <- map_dfr(TRIGGER_WINDOWS_H, function(window_h) {
  reduced <- patient_course[!is.na(reduction_t_h), .(hospitalization_id, t_h = reduction_t_h)]
  reduced <- triggers_at(reduced, window_h)
  course <- merge(patient_course, reduced[, .(hospitalization_id, dp_trigger, crs_trigger, sf_trigger)],
                  by = "hospitalization_id", all.x = TRUE)
  summarise_group <- function(frame) {
    frame %>% summarise(
      n_patients = n(), n_in_paper_cohort = sum(in_paper_cohort),
      initial_vtpbw_mean = mean(initial_vtpbw, na.rm = TRUE),
      n_reduced = sum(!is.na(reduction_t_h)),
      hours_to_first_reduction_median = median(reduction_t_h, na.rm = TRUE),
      n_reduced_no_plateau = sum(dp_trigger %in% "no plateau measured"),
      n_reduced_dp_high = sum(dp_trigger %in% "DP high"),
      n_reduced_crs_fall = sum(crs_trigger %in% "Crs fall"),
      n_reduced_plateau_neither = sum(dp_trigger %in% "DP not high" & !crs_trigger %in% "Crs fall"),
      n_reduced_sf_fall = sum(sf_trigger %in% "SF fall"),
      .groups = "drop")
  }
  bind_rows(
    course %>% group_by(pfvc_quartile) %>% summarise_group() %>% mutate(initial_band = "all"),
    course %>% group_by(pfvc_quartile, initial_band) %>% summarise_group(),
    course %>% filter(in_paper_cohort) %>% group_by(pfvc_quartile, initial_band) %>% summarise_group() %>%
      mutate(population = "paper cohort")) %>%
    mutate(population = replace_na(population, "ungated"), trigger_window_h = window_h,
           share_reduced = n_reduced / n_patients,
           share_of_reductions_dp_high = n_reduced_dp_high / n_reduced,
           share_of_reductions_crs_fall = n_reduced_crs_fall / n_reduced,
           share_of_reductions_no_plateau = n_reduced_no_plateau / n_reduced)
}) %>%
  mutate(site = site_name, .before = 1) %>%
  relocate(trigger_window_h, population, pfvc_quartile, initial_band, .after = site)
write_csv(events, file.path(final_dir, paste0("vt_reduction_timing_events_", site_name, ".csv")))

# =============================================================================
# 2. Decomposition: hourly hazard of the first reduction
# =============================================================================
# At risk in hour h: the patient has a row at or before h, is in a volume-targeted
# mode at the start of h, and has had no reduction before h. The state at the start
# of h is the last row at or before it.
last_row_time <- timepoints[, .(last_t_h = max(t_h)), by = hospitalization_id]
hours <- last_row_time[, .(hour = seq(0, floor(last_t_h))), by = hospitalization_id]
setkey(timepoints, hospitalization_id, t_h)
state <- timepoints[hours[, .(hospitalization_id, t_h = as.numeric(hour))], on = .(hospitalization_id, t_h),
                    roll = Inf, mult = "last",
                    .(hospitalization_id, hour = i.t_h, volume_mode, tidal_volume_set, pbw)]
state <- merge(state, patient_course[, .(hospitalization_id, reduction_t_h)], by = "hospitalization_id", all.x = TRUE)
at_risk <- state[volume_mode %in% TRUE & (is.na(reduction_t_h) | hour <= floor(reduction_t_h)),
                 .(hospitalization_id, t_h = hour, vtpbw_now = tidal_volume_set / pbw,
                   event = as.integer(!is.na(reduction_t_h) & hour == floor(reduction_t_h)))]
message("Hazard: ", nrow(at_risk), " patient-hours at risk, ", sum(at_risk$event), " first reductions")

STEPS <- c("size alone"                = "",
           "+ current VT/PBW (i)"      = "+ ns(vtpbw_now, 3)",
           "+ DP trigger (ii)"         = "+ ns(vtpbw_now, 3) + dp_trigger",
           "+ Crs and SF triggers (iii)" = "+ ns(vtpbw_now, 3) + dp_trigger + crs_change + sf_trigger")
TRIGGER_TERMS <- "^(dp_trigger|crs_change|sf_trigger)"
hazard <- map_dfr(TRIGGER_WINDOWS_H, function(window_h) {
  frame <- triggers_at(copy(at_risk), window_h) %>%
    merge(as.data.table(patients), by = "hospitalization_id") %>%
    as_tibble()
  expand_grid(exposure = names(EXPOSURES), adjustment = names(ADJUSTMENTS), step = names(STEPS)) %>%
    pmap_dfr(function(exposure, adjustment, step) {
      size_term <- paste0(EXPOSURES[[exposure]], "_sd")
      formula <- paste("event ~ ns(t_h, 3) +", size_term, STEPS[[step]], ADJUSTMENTS[[adjustment]])
      model <- glm(as.formula(formula), data = frame, family = binomial)
      if (!model$converged) stop("Did not converge: ", formula)
      if (anyNA(coef(model))) stop("Aliased terms (", paste(names(coef(model))[is.na(coef(model))], collapse = ", "),
                                   "): ", formula)
      coefficients <- summary(model)$coefficients
      keep <- rownames(coefficients) == size_term | grepl(TRIGGER_TERMS, rownames(coefficients))
      tibble(trigger_window_h = window_h, exposure, adjustment, step,
             term = rownames(coefficients)[keep], estimate = coefficients[keep, 1], se = coefficients[keep, 2],
             patient_hours = stats::nobs(model), events = sum(frame$event))
    })
}) %>%
  mutate(site = site_name, exposure_sd = exposure_sd[exposure], .before = 1) %>%
  mutate(step = factor(step, levels = names(STEPS)))
write_csv(hazard, file.path(final_dir, paste0("vt_reduction_timing_hazard_", site_name, ".csv")))

# =============================================================================
# 3. Link 1: does size predict a fall in compliance (or SF) at all?
# =============================================================================
marker_windows <- function(events_table, value) {
  events_table[, .(baseline = median(get(value)[t_h < BASELINE_H]),
                   followup = median(get(value)[t_h >= FOLLOWUP_H[1] & t_h < FOLLOWUP_H[2]])),
               by = hospitalization_id]
}
MARKERS <- list("log Crs" = marker_windows(plateau_events, "log_crs"), "SF" = marker_windows(sf_events, "sf_ratio"))
ancova <- imap_dfr(MARKERS, function(marker_values, marker) {
  frame <- patient_course %>% as_tibble() %>% inner_join(as_tibble(marker_values), by = "hospitalization_id")
  n_with_baseline <- sum(!is.na(frame$baseline))
  frame <- frame %>% filter(!is.na(baseline), !is.na(followup))
  expand_grid(exposure = names(EXPOSURES), adjustment = names(ADJUSTMENTS),
              initial_vtpbw_term = c("without initial VT/PBW", "with initial VT/PBW")) %>%
    pmap_dfr(function(exposure, adjustment, initial_vtpbw_term) {
      size_term <- paste0(EXPOSURES[[exposure]], "_sd")
      formula <- paste("followup ~", size_term, "+ baseline",
                       if (initial_vtpbw_term == "with initial VT/PBW") "+ initial_vtpbw" else "", ADJUSTMENTS[[adjustment]])
      coefficients <- summary(lm(as.formula(formula), data = frame))$coefficients
      tibble(marker, exposure, adjustment, initial_vtpbw_term,
             estimate = coefficients[size_term, 1], se = coefficients[size_term, 2],
             n = nrow(frame), n_with_baseline)
    })
}) %>%
  mutate(site = site_name, exposure_sd = exposure_sd[exposure], .before = 1)
write_csv(ancova, file.path(final_dir, paste0("vt_reduction_timing_ancova_", site_name, ".csv")))

# --- figure: the size term across the decomposition (primary window)
ADJUSTMENT_COLOURS <- c(unadjusted = "#E69F00", adjusted = "#0072B2")    # Okabe-Ito
decomposition_plot <- hazard %>%
  filter(trigger_window_h == TRIGGER_WINDOWS_H[1], !grepl(TRIGGER_TERMS, term)) %>%
  mutate(step = fct_rev(step)) %>%
  ggplot(aes(estimate, step, colour = adjustment)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_pointrange(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se),
                  position = position_dodge(width = 0.5), size = 0.25) +
  facet_wrap(~ exposure, scales = "free_x") +
  scale_colour_manual(values = ADJUSTMENT_COLOURS, name = NULL) +
  labs(x = "Log-odds of a first VT reduction in the hour, per SD of the size term", y = NULL,
       title = paste0("Which arrow carries size into tidal-volume reductions? (", site_name, ", ",
                      TRIGGER_WINDOWS_H[1], "-h trigger window)")) +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(final_dir, paste0("vt_reduction_timing_", site_name, ".pdf")), decomposition_plot, width = 8, height = 3.5)

message("Wrote vt_reduction_timing_* to ", final_dir)
print(events %>% filter(trigger_window_h == TRIGGER_WINDOWS_H[1], population == "ungated", initial_band == "all") %>%
        select(pfvc_quartile, n_patients, initial_vtpbw_mean, share_reduced, share_of_reductions_no_plateau,
               share_of_reductions_dp_high, share_of_reductions_crs_fall))
print(hazard %>% filter(trigger_window_h == TRIGGER_WINDOWS_H[1], exposure == "log PFVC", !grepl(TRIGGER_TERMS, term)) %>%
        select(adjustment, step, estimate, se))
print(ancova %>% filter(exposure == "log PFVC") %>% select(marker, adjustment, initial_vtpbw_term, estimate, se, n))
