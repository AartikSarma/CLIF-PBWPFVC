# =============================================================================
# Supplement (causal): does tidal volume, as it changes within a patient, move the marker?
# =============================================================================
#
# Strain acts over time: a marker's course should follow the tidal volume the lung
# receives day by day, not only the volume set at intubation. Every design in
# causal_triangulation.R compares patients; this one compares a patient with
# themselves. Each ventilated patient gets their own intercept and their own slope in
# time, so anything fixed about the patient (age, sex, race, height, PBW, PFVC, so
# lambda, and baseline illness) drops out. What is left to identify the effect is when,
# and by how much, the patient's set tidal volume was changed.
#
# What it can and cannot answer. A fixed lung size cancels within a patient, so this
# asks whether TIDAL VOLUME matters for the marker, not whether PFVC is the right
# denominator. Under strain the effect of a log unit of tidal volume equals the effect
# of a log unit of lambda (strain = log VT - log lung size), so the estimate is compared
# with the lambda gradient of the unadjusted design on the same patients.
#
# Pre-specified analysis (fixed before any real-data run of this script):
#   Sample      ventilated patients (the cohort of config: gated unless
#               PBWPFVC_VTPBW_GATE=0) with a volume-targeted set tidal volume for at least
#               LONG_MIN_COVER (half) of their ventilator days, and at least
#               LONG_MIN_OBS (3) marker values on the ventilator. Ventilator days run
#               from t0 to liberation or ON_DAYS (7), as in causal_triangulation.R.
#   Exposure    log set tidal volume on volume-targeted modes (assist control-volume
#               control, pressure-regulated volume control, SIMV: VOLUME_TARGETED_MODES),
#               as a step function in time: each charted value holds until the next,
#               for at most LONG_HOLD_H (12) hours. Time on other modes, or past the
#               hold, contributes nothing (K, the volume-targeted time, is a covariate).
#               C(t) is the time integral from t0 to t of log VT about the patient's own
#               time-weighted mean.
#   Model       log marker = patient intercept + patient slope x t + common natural
#               spline in t (3 df) + b C(t) + K(t) + the integrals of FiO2, PEEP and
#               norepinephrine-equivalent dose, each LONG_LAG_H (6) hours earlier,
#               least squares, standard errors clustered on the patient
#               (fit_longitudinal() in causal_common.R).
#   Primary     b for platelets: the change in the marker's slope (log units per day) per
#               log unit of tidal volume. Strain predicts b < 0 for platelets, > 0 for
#               creatinine, and b equal to the lambda gradient (p_equal_lambda).
#   Secondary   creatinine; the model without the severity terms; the falsification row,
#               which adds tidal volume LONG_LEAD_H (24) hours in the FUTURE (C_lead) and
#               should find nothing; and C x lambda, whether a change in tidal volume
#               acts more where lambda is high (a linear, not logarithmic, dose).
#
# Read the description table first. If tidal volume rarely changes, or changes mostly at
# mode switches or in the last day before liberation (weaning), the estimate is a
# contrast of those moments, not of strain; the table says which.
#
# Limits. Tidal volume set lower because the lung got worse is the main threat; the
# lagged FiO2, PEEP and pressor terms adjust for the severity the chart shows, not for
# all of it, and they are themselves partly the result of earlier tidal volumes (a
# g-method would be needed to separate that). Platelet transfusions are not removed.
# Pressure-targeted time has no set tidal volume and is left out, not imputed.
#
# Writes to final/supplement/ (final/ungated/supplement/ with PBWPFVC_VTPBW_GATE=0):
#   causal_longitudinal_description_<site>.csv   how tidal volume changed within patients
#   causal_longitudinal_modes_<site>.csv         ventilator records by mode, VT set or not
#   causal_longitudinal_estimates_<site>.csv     marker x model x term
# No patient-level row is written.
#
# Runtime: a minute or two. Needs scripts 01-03 (and PBWPFVC_COHORT=nosupport, which
# load_causal_patients() reads, although the control is not used here).
#
# Usage: uvr run code/supplement/causal_longitudinal.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(here)
})

rm(list = ls())
source("utils/config.R")
site_name <- config$site_name
final_dir <- final_dir_for("supplement")
source(here::here("code", "20_biotrauma_grid.R"))
source(here::here("code", "supplement", "causal_common.R"))

MARKERS <- strsplit(Sys.getenv("PBWPFVC_CAUSAL_MARKERS", "platelets,creatinine"), ",")[[1]]
stopifnot(all(MARKERS %in% CAUSAL_MARKERS$marker[CAUSAL_MARKERS$source == "lab"]))
cohort <- if (config$vtpbw_gate == "0") "ungated" else "gated"
patients <- load_causal_patients(config)
records <- load_vent_records(config, patients)
vp <- patients[patients$vent == 1, ]
t_end <- pmin(coalesce(vp$t_off_d, ON_DAYS), ON_DAYS)
message("Ventilated patients: ", nrow(vp), "; ventilator records: ", nrow(records$vent))

modes <- records$vent %>% filter(s >= 0, s < t_end[match(id, vp$id)]) %>%
  count(mode, vt_set = !is.na(vt) & vt > 0, name = "records") %>%
  mutate(volume_targeted = mode %in% VOLUME_TARGETED_MODES, site = site_name, cohort = cohort, .before = 1)
description <- describe_vt_changes(vp, records$vent, t_end) %>% mutate(site = site_name, cohort = cohort, .before = 1)

MODELS <- list(primary = c("C", "K", "Z_fio2", "Z_peep", "Z_pressor"),
               no_severity = c("C", "K"),
               falsification_lead = c("C", "C_lead", "K", "Z_fio2", "Z_peep", "Z_pressor"),
               interaction_lambda = c("C", "CL", "K", "Z_fio2", "Z_peep", "Z_pressor"))
REPORT <- c("C", "C_lead", "CL", "Z_fio2", "Z_peep", "Z_pressor")

estimates <- list()
for (marker in MARKERS) {
  m <- CAUSAL_MARKERS[CAUSAL_MARKERS$marker == marker, ]
  panel <- load_causal_outcome(config, patients, marker, complete = FALSE)
  frame <- longitudinal_frame(patients, panel, records)
  if (length(unique(frame$pid)) < MIN_PATIENTS) { message("  ", marker, ": too few patients"); next }
  ids <- unique(frame$id)
  # the lambda gradient of the unadjusted design on the same patients, with its influences
  design <- causal_designs(TRUE) %>% filter(design == "unadjusted")
  pats <- patients[patients$id %in% ids, ]; pan <- panel[panel$id %in% ids, ]
  prep <- prepare_design(design, pats, pan)
  lam <- if (prep$ok) fit_coefs(prep, align_outcome(prep, pan)) else NULL
  lam_est <- if (is.null(lam)) NA_real_ else lam$beta[[prep$k]]
  lam_infl <- if (is.null(lam)) NULL else setNames(lam$influence[, prep$k], prep$patient_ids[lam$patient])
  for (model in names(MODELS)) {
    fit <- fit_longitudinal(frame, MODELS[[model]])
    for (term in intersect(REPORT, names(fit$beta))) {
      j <- which(names(fit$beta) == term)
      row <- tibble(marker = marker, model = model, term = term, estimate = fit$beta[[j]], se = sqrt(fit$V[j, j]),
                    n_patients = fit$n_patients, n_obs = fit$n_obs, lambda_gradient = lam_est,
                    lambda_gradient_se = if (is.null(lam)) NA_real_ else sqrt(lam$V[prep$k, prep$k]), p_equal_lambda = NA_real_)
      if (term == "C" && !is.null(lam_infl)) {
        # the two estimators share patients: the covariance from their per-patient influences
        a <- setNames(fit$influence[, j], frame$id[match(fit$pid, frame$pid)])
        common <- intersect(names(a), names(lam_infl))
        G <- length(union(names(a), names(lam_infl)))
        cov_ab <- sum(a[common] * lam_infl[common]) * G / (G - 1)
        row$p_equal_lambda <- 2 * pnorm(-abs(row$estimate - lam_est) / sqrt(row$se^2 + row$lambda_gradient_se^2 - 2 * cov_ab))
      }
      estimates[[paste(marker, model, term)]] <- row
    }
  }
  message("  ", marker, ": ", length(ids), " patients, ", nrow(frame), " values (of ", attr(frame, "patients_on"),
          " with a value on the ventilator)")
}
estimates <- bind_rows(estimates) %>%
  mutate(site = site_name, cohort = cohort, .before = 1) %>%
  mutate(conf_low = estimate - qnorm(0.975) * se, conf_high = estimate + qnorm(0.975) * se,
         p_value = 2 * pnorm(-abs(estimate / se)), .after = se)

out_file <- function(name) file.path(final_dir, paste0("causal_longitudinal_", name, site_name, ".csv"))
write_csv(description, out_file("description_"))
write_csv(modes, out_file("modes_"))
write_csv(estimates, out_file("estimates_"))

message("\nHow tidal volume changed within patients (read this first):")
print(as.data.frame(description %>% select(-site, -cohort) %>% mutate(across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)
message("\nVentilator records on the ventilator days, by mode:")
print(as.data.frame(modes %>% select(mode, volume_targeted, vt_set, records)), row.names = FALSE)
message("\nEstimates (C: log VT integrated over time; strain predicts C = the lambda gradient; C_lead should be 0):")
print(as.data.frame(estimates %>% transmute(marker, model, term, estimate = signif(estimate, 3), se = signif(se, 3),
                                            p_value = signif(p_value, 2), lambda_gradient = signif(lambda_gradient, 3),
                                            p_equal_lambda = signif(p_equal_lambda, 2), n_patients, n_obs)), row.names = FALSE)
message("\ncausal_longitudinal complete: ", final_dir)
