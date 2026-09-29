# =============================================================================
# Script 26 (measurement check): is a marker missing in a way that biases figure 4?
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# Platelets and creatinine are drawn as a matter of routine; bilirubin is drawn
# when someone suspects a liver. In the 7-day panels bilirubin reaches about
# two thirds of ventilated patients but only a third to a half of the no-support
# controls, and the patients who get one are sicker in both arms -- more so in
# the control arm. A difference-in-differences survives selection that is equally
# strong on both sides; it does not survive selection that differs between them.
#
# Two tests, neither of which needs a joint model:
#
#   frame check (the falsification).  Refit the PLATELET divergence twice with one
#     model: on every platelet day, and on the platelet days that also carry a
#     bilirubin. Platelets are measured on nearly every panel day, so the second
#     fit changes the sampling frame and nothing else. If the platelet divergence
#     survives on bilirubin's frame, that frame is not what distorts bilirubin; if
#     it moves, the distortion has been measured directly, in a marker whose
#     answer we already trust. The restricted sample is nested in the full one, so
#     the difference between the two carries the Hausman variance, var(restricted)
#     minus var(full), and not the restricted fit's own standard error: a subset is
#     not an independent study, and treating it as one reads a shift as bigger than
#     the smaller sample can support. Written to measurement_frame_{tag}.csv.
#
#   measurement model (the mechanism).  Does predicted lung size predict being
#     measured at all? Two parts per marker:
#       entry        one row per patient: does the patient enter that marker's
#                    analysis (a day-0 baseline and two or more later days, the
#                    joint model's own requirement)? Logistic regression.
#       measurement  one row per patient-day: was the marker resulted that day?
#                    Logistic mixed model with a patient random intercept, and
#                    log PFVC x day, so that a measurement rate which DIVERGES by
#                    lung size over the week is visible and not just a level
#                    difference.
#     Run for bilirubin and, as the reference, for platelets and creatinine: the
#     question is not whether sicker patients are measured more, which they are,
#     but whether bilirubin depends on lung size in a way the routine labs do not.
#     Written to measurement_model_{tag}.csv.
#
# Both are reported demographic-adjusted and unadjusted, as everything here is.
#
# What the denominator is, and is not. The panel holds one row per patient-day
# that carries a respiratory-support record, so a day after ICU discharge is not
# in it and is not counted as an unmeasured day. These tests therefore ask
# whether measurement depends on lung size AMONG THE DAYS THAT ENTER FIGURE 4,
# which is the selection figure 4 is exposed to. Days lost to discharge are a
# different question, and the no-support control arm is what addresses it.
#
# The trajectory of each marker starts the day after its baseline (21), so the
# day-level model is about measurement after the baseline draw.
#
# Usage (per cohort, as 29_run_figure4.R runs it):
#   Rscript code/26_measurement_check.R                      # the ventilated cohort
#   PBWPFVC_COHORT=nosupport Rscript code/26_measurement_check.R
# Switches: PBWPFVC_JM_CLOCK (index or icu, as for the arm being checked),
#   PBWPFVC_MEAS_MARKERS (default bilirubin,platelets,creatinine).
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(arrow); library(splines)
  library(nlme); library(GLMMadaptive)
})
rm(list = ls())
source("utils/config.R")

site_name  <- config$site_name
output_dir <- config$output_dir
final_dir  <- final_dir_for("injury")
source(here("code", "20_biotrauma_grid.R"))   # JM_GRID, JM_HORIZON, h_suffix, panel_path, clock_tag

MARKERS <- strsplit(Sys.getenv("PBWPFVC_MEAS_MARKERS", "bilirubin,platelets,creatinine"), ",")[[1]]
BASELINE <- c(bilirubin = "bilirubin_0", platelets = "platelet_0", creatinine = "creatinine_0")
# the sparsely measured marker whose sampling frame the frame check imposes on
# platelets; a knob so the check can be exercised on a dataset without bilirubin
FRAME_MARKER <- Sys.getenv("PBWPFVC_MEAS_FRAME", "bilirubin")
out_tag <- paste0(clock_tag, h_suffix, "_", site_name)
message("=== 26_measurement_check: ", site_name, ", cohort ", config$cohort, ", clock ", JM_CLOCK,
        ", markers ", paste(MARKERS, collapse = ", "))

long <- read_parquet(panel_path("long"))
surv <- read_parquet(panel_path("surv"))

# One model frame for every test: the panel's patient-days with the covariates the
# joint model adjusts for, plus a measured/not indicator per marker.
frame <- long %>%
  select(hospitalization_id, vent_day, any_of(MARKERS)) %>%
  inner_join(surv %>% select(hospitalization_id, np_sofa, sf_0, age10, sex_category, race_category,
                             log_pfvc_sd, any_of(unname(BASELINE))),
             by = "hospitalization_id") %>%
  filter(!is.na(np_sofa), !is.na(sf_0), !is.na(log_pfvc_sd)) %>%
  # the logistic mixed model below diverges when its covariates sit on scales as
  # different as a SOFA score, a log SF ratio and a day count, so the continuous
  # adjusters are standardised and the day is centred. log_pfvc_sd is already a
  # z-score, and the reported terms are its own, so nothing reported changes scale.
  mutate(log_sf_0 = as.numeric(scale(log(sf_0))), np_sofa = as.numeric(scale(np_sofa)),
         vent_day_c = vent_day - mean(vent_day))
for (m in MARKERS) frame[[paste0("measured_", m)]] <- as.integer(!is.na(frame[[m]]))
message("Panel days in the frame: ", nrow(frame), " for ", n_distinct(frame$hospitalization_id), " patients")

DEMO <- "ns(age10, 4) + sex_category + race_category"
BASE <- "np_sofa + log_sf_0"
rhs  <- function(adjusted) paste(BASE, if (adjusted) paste("+", DEMO) else "")
# the estimates are per SD of log PFVC, the unit every figure-4 table carries;
# the pooling converts it to per 0.1 log units with each site's own SD
EXPOSURE <- "log_pfvc_sd"
tidy_row <- function(est, se, ...) tibble(estimate = est, se = se, lo = est - 1.96 * se,
                                          hi = est + 1.96 * se, p = 2 * pnorm(-abs(est / se)), ...)

# --- 1. the measurement model ------------------------------------------------
# entry: does this patient contribute a trajectory for this marker at all?
entry_rows <- list(); day_rows <- list()
for (m in MARKERS) {
  y0 <- BASELINE[[m]]
  if (!y0 %in% names(frame)) { message("  ", m, ": no baseline column ", y0, ", skipped"); next }
  measured <- paste0("measured_", m)
  per_patient <- frame %>%
    group_by(hospitalization_id) %>%
    summarise(in_sample = as.integer(!is.na(first(.data[[y0]])) & sum(.data[[measured]]) >= 2),
              across(c(np_sofa, log_sf_0, age10, sex_category, race_category, log_pfvc_sd), first),
              .groups = "drop")
  message("  ", m, ": ", sum(per_patient$in_sample), " of ", nrow(per_patient),
          " patients enter its analysis (", round(100 * mean(per_patient$in_sample)), "%); ",
          sum(frame[[measured]]), " of ", nrow(frame), " panel days carry a value")
  for (adjusted in c(TRUE, FALSE)) {
    # a coefficient this large is separation, not an association: it is flagged in
    # the table rather than left to be read as a huge odds ratio. R's own warning is
    # left to print.
    fit <- glm(as.formula(paste("in_sample ~", EXPOSURE, "+", rhs(adjusted))),
               family = binomial(), data = per_patient)
    s <- summary(fit)$coefficients[EXPOSURE, ]
    separated <- !fit$converged || abs(s[["Estimate"]]) > 10 || s[["Std. Error"]] > 10
    if (separated)
      message("  ", m, " (", if (adjusted) "adjusted" else "unadjusted",
              "): the entry model separated; its row is flagged, not interpreted")
    entry_rows[[length(entry_rows) + 1L]] <- tidy_row(
      s[["Estimate"]], s[["Std. Error"]], test = "entry", marker = m,
      adjustment = if (adjusted) "adjusted" else "unadjusted", term = EXPOSURE,
      n_patients = nrow(per_patient), n_in_sample = sum(per_patient$in_sample), n_obs = NA_integer_,
      separated = separated)
  }
  # measurement: was it resulted on this day, and does that DIVERGE by lung size?
  for (adjusted in c(TRUE, FALSE)) {
    fixed <- as.formula(paste(measured, "~", EXPOSURE, "* vent_day_c +", rhs(adjusted)))
    # the EM start diverges on some panels; the message GLMMadaptive prints for that
    # names iter_EM = 0 as the remedy, so that is the second attempt, not a fallback
    # to a different model
    fit <- tryCatch(mixed_model(fixed = fixed, random = ~ 1 | hospitalization_id,
                                data = frame, family = binomial()),
                    error = function(e) {
                      message("  ", m, " (", if (adjusted) "adjusted" else "unadjusted",
                              "): the measurement model diverged from the EM start; refitting with iter_EM = 0")
                      tryCatch(mixed_model(fixed = fixed, random = ~ 1 | hospitalization_id,
                                           data = frame, family = binomial(),
                                           control = list(iter_EM = 0)),
                               error = function(e2) {
                                 message("  ", m, ": measurement model failed: ", conditionMessage(e2)); NULL })
                    })
    if (is.null(fit)) next
    fe <- fixef(fit); V <- vcov(fit, parm = "fixed-effects")
    for (tm in intersect(c(EXPOSURE, paste0(EXPOSURE, ":vent_day_c"), paste0("vent_day_c:", EXPOSURE)), names(fe)))
      day_rows[[length(day_rows) + 1L]] <- tidy_row(
        unname(fe[[tm]]), sqrt(V[tm, tm]), test = "measurement", marker = m,
        adjustment = if (adjusted) "adjusted" else "unadjusted",
        term = if (tm == EXPOSURE) "level" else "divergence per day",
        n_patients = n_distinct(frame$hospitalization_id), n_in_sample = NA_integer_, n_obs = nrow(frame))
  }
}
measurement <- bind_rows(entry_rows, day_rows) %>%
  mutate(odds_ratio = exp(estimate), or_lo = exp(lo), or_hi = exp(hi),
         unit = "per SD of log PFVC (log-odds; odds ratio in odds_ratio)",
         scale = if_else(test == "entry", "odds of entering the marker's analysis",
                         "odds that the marker is resulted on a panel day"),
         cohort = config$cohort, clock = JM_CLOCK, panel = h_suffix, site = site_name)
write_csv(measurement, file.path(final_dir, paste0("measurement_model_", out_tag, ".csv")))
message("measurement model -> ", file.path(final_dir, paste0("measurement_model_", out_tag, ".csv")))

# --- 2. the frame check ------------------------------------------------------
# The same platelet model on two sampling frames. Nothing but the rows changes, so
# a difference between the two is the frame's doing -- or the smaller sample's noise,
# which is what the nested variance below separates.
frame_rows <- list()
if (all(c("platelets", FRAME_MARKER) %in% MARKERS) && FRAME_MARKER != "platelets") {
  platelet_fit <- function(d, adjusted, label) {
    d <- d %>% filter(!is.na(platelets), !is.na(platelet_0), platelet_0 > 0) %>%
      mutate(log_y = log(platelets), log_y0 = log(platelet_0)) %>%
      group_by(hospitalization_id) %>% filter(n() >= 2) %>% ungroup()
    message("  frame check, ", label, ": ", n_distinct(d$hospitalization_id), " patients, ", nrow(d), " days")
    if (n_distinct(d$hospitalization_id) < 20) {
      message("  frame check, ", label, ": too few patients, not fitted"); return(NULL)
    }
    f <- as.formula(paste("log_y ~ ns(vent_day, 3) +", EXPOSURE, "+", EXPOSURE, ": vent_day + log_y0 +", rhs(adjusted)))
    fit <- tryCatch(lme(f, random = ~ vent_day | hospitalization_id, data = d,
                        control = lmeControl(opt = "optim", msMaxIter = 200)),
                    error = function(e) { message("  frame check failed: ", conditionMessage(e)); NULL })
    if (is.null(fit)) return(NULL)
    fe <- fixef(fit); V <- vcov(fit)
    tm <- intersect(c(paste0(EXPOSURE, ":vent_day"), paste0("vent_day:", EXPOSURE)), names(fe))[1]
    list(estimate = unname(fe[[tm]]), se = sqrt(V[tm, tm]),
         n_patients = n_distinct(d$hospitalization_id), n_obs = nrow(d))
  }
  for (adjusted in c(TRUE, FALSE)) {
    adj_note <- if (adjusted) "adjusted" else "unadjusted"
    restricted_label <- paste0("platelet days carrying a ", FRAME_MARKER)
    all_days <- platelet_fit(frame, adjusted, paste(adj_note, "every platelet day"))
    on_bili  <- platelet_fit(frame %>% filter(.data[[paste0("measured_", FRAME_MARKER)]] == 1),
                             adjusted, paste(adj_note, restricted_label))
    if (is.null(all_days) || is.null(on_bili)) next
    adj_label <- if (adjusted) "adjusted" else "unadjusted"
    for (s in list(c(list(sample = "every platelet day"), all_days),
                   c(list(sample = restricted_label), on_bili)))
      frame_rows[[length(frame_rows) + 1L]] <- tidy_row(
        s$estimate, s$se, marker = "platelets", adjustment = adj_label, sample = s$sample,
        n_patients = s$n_patients, n_obs = s$n_obs)
    # The restricted sample is NESTED in the full one, so the two estimates are
    # correlated and the difference does not carry the restricted fit's standard
    # error: under the full model it is var(restricted) - var(full), the Hausman
    # form. Using the restricted SE instead treats a subset as an independent study
    # and reads a shift as larger than the smaller sample can support. Where the
    # difference of variances is not positive the assumption behind it has failed,
    # and no z is reported rather than a fabricated one.
    diff <- on_bili$estimate - all_days$estimate
    nested_var <- on_bili$se^2 - all_days$se^2
    nested_se  <- if (nested_var > 0) sqrt(nested_var) else NA_real_
    frame_rows[[length(frame_rows) + 1L]] <- tibble(
      estimate = diff, se = nested_se,
      lo = if (is.na(nested_se)) NA_real_ else diff - 1.96 * nested_se,
      hi = if (is.na(nested_se)) NA_real_ else diff + 1.96 * nested_se,
      p = if (is.na(nested_se)) NA_real_ else 2 * pnorm(-abs(diff / nested_se)),
      marker = "platelets", adjustment = adj_label,
      sample = "difference (restricted minus every day)",
      n_patients = on_bili$n_patients, n_obs = on_bili$n_obs,
      z_nested = if (is.na(nested_se)) NA_real_ else diff / nested_se)
  }
}
if (length(frame_rows)) {
  frame_check <- bind_rows(frame_rows) %>%
    mutate(term = "divergence per day", unit = "log platelets per day per SD of log PFVC",
           cohort = config$cohort, clock = JM_CLOCK, panel = h_suffix, site = site_name)
  write_csv(frame_check, file.path(final_dir, paste0("measurement_frame_", out_tag, ".csv")))
  message("frame check -> ", file.path(final_dir, paste0("measurement_frame_", out_tag, ".csv")))
  print(as.data.frame(frame_check %>% filter(adjustment == "adjusted") %>%
                        transmute(sample, estimate = signif(estimate, 3), lo = signif(lo, 3),
                                  hi = signif(hi, 3), n_patients, n_obs)), row.names = FALSE)
} else if (!all(c("platelets", FRAME_MARKER) %in% MARKERS)) {
  message("frame check: not run, it needs platelets and ", FRAME_MARKER, " in PBWPFVC_MEAS_MARKERS")
} else {
  message("frame check: no fit produced a divergence term; see the messages above")
}
message("=== 26_measurement_check complete")
