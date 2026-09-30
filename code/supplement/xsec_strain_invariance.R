# =============================================================================
# Supplement (cross-sectional): does ventilation's mortality gradient follow the
# strain error, or any index of the same demographics?
# =============================================================================
# PBW/PFVC is a fixed function of age, sex, race and height, so adjusting for the
# demographics removes the variation it is made of. This script uses no demographic
# adjustment. Ventilation's effect is isolated by comparison with the no-support
# control: patients with the same formulas, the same demographic paths to death and
# no tidal volume. The weak point of that design is age: the ratio is mostly age, and
# age could matter more under ventilation for reasons other than strain. Each section
# below gives the data a way to separate the two.
#
# "Unadjusted" here means severity only: SF and SOFA, standardised within cohort. No
# age, sex, race or height term. (xsec_pfvc_age_control.R's "unadjusted" arm keeps a
# 4-df age spline; this one does not.)
#
# 1. The site's dosing rule. How set tidal volume scales with PBW decides how strain
#    depends on the demographics. Under strict per-kg dosing, VT/PFVC moves with
#    PBW/PFVC (mostly age, sex and race); under a fixed volume in mL, it moves with
#    1/PFVC (height too). Frailty does not know the site's rule. The rule is the slope
#    of log set VT on log PBW at the first complete IMV row of every ventilated
#    patient, BEFORE the lung-protective VT/PBW 6-8 gate (a slope measured inside the
#    gate is pushed toward 1 by the gate itself): 1 = per-kg dosing, 0 = a fixed
#    volume. Reported pooled and within sex (a site that fixes different volumes for
#    men and women looks per-kg on the pooled slope), with the SDs of set VT (mL) and
#    VT/PBW. The site-level test (does the ratio's ventilation-specific effect grow
#    with the slope, and PFVC's shrink?) is a pooling step across sites.
#
# 2. The ventilation contrast. Per cohort, in-hospital death (logistic) and 60-day
#    death before invasive ventilation (cause-specific Cox):
#      death ~ SF (z) + SOFA (z) + exposure       [+ VT/PBW, ventilated, "at fixed VT/PBW"]
#    for log PBW/PFVC and log PFVC, per 0.1 log units. The difference, ventilated
#    minus control, is ventilation's share. Two dose versions of the ventilated side:
#      no dose           the association as the site's practice delivers it: under
#                        strain, the ratio's difference grows with the site's per-kg
#                        slope and PFVC's shrinks
#      at fixed VT/PBW   at a given VT/PBW, log PBW/PFVC is log VT/PFVC less a
#                        constant: the delivered strain. The control has no dose term,
#                        so this difference is asymmetric by construction.
#    Also the ventilated cohort's delivered strain on its own (log VT/PFVC, severity
#    only), the quantity whose coefficient strain predicts to be the same at every
#    site. Confounding by indication could also be the same at every site, so the
#    discriminating prediction is the proxies moving with practice, not the strain
#    coefficient holding still.
#    Populations: everyone, and hypoxemic at the index (SF < 315, the ventilated
#    cohort's own gate), the control matched on the ventilated cohort's hypoxemia.
#
# 3. Placebo formulas. If the ratio's ventilation-specific gradient is strain, GLI's
#    weighting of the inputs should matter; if it is any demographic index, arbitrary
#    weightings should do as well. The inputs: a 4-df natural spline in age, sex, race
#    (Black, Other) and log height, whitened over both cohorts together so every
#    direction has SD 1 on one scale. PLACEBO_N random directions on the unit sphere
#    are each run through the in-hospital contrast (everyone, no dose). Beside them,
#    named reference indices, each standardised the same way: log PBW/PFVC, log PFVC,
#    the ratio's projection onto the inputs (R2 reported), and GLI's four pieces
#    (height, age, sex, race; pfvc_channels() in 20_biotrauma_grid.R). Read where each
#    reference falls in the placebo cloud, not GLI's percentile alone: the ratio is
#    mostly age, so if age alone lands in the tail the ratio will too.
#    The statistic is |z| of the difference; random directions have no sign.
#
# The cohorts are xsec_mortality_channel_equality.R's (copied from it, same
# synthetic seed): the ventilated arm on IMV at ICU admission (icu_day0), gated at
# VT/PBW 6-8 and SF < 315 by script 03, and the no-support control.
#
# Inputs : intermediate/analysis_cross_sectional.parquet, analysis_all_eligible_timepoints.parquet
#          (script 03, ventilated), intermediate/controls/nosupport/analysis_cross_sectional.parquet
#          and resp_support_waterfall_clean.parquet (scripts 01-03, PBWPFVC_COHORT=nosupport)
# Outputs: final/supplement/
#   strain_invariance_dosing_{site}.csv     the dosing rule: slopes, R2, SDs, before the
#                                           VT/PBW gate and (for reference) inside it
#   strain_invariance_contrast_{site}.csv   per cohort, the difference, and the delivered-
#                                           strain row, by population, outcome and dose
#   strain_invariance_placebo_{site}.csv    every index: the two cohorts' coefficients,
#                                           the difference, its z, and (references) the
#                                           percentile of |z| among the placebos
# Usage: uvr run code/supplement/xsec_strain_invariance.R   (PBWPFVC_COHORT unset)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(survival)
  library(splines)
})

options(width = 220)
source("utils/config.R")
if (config$cohort != "imv") stop("xsec_strain_invariance.R reads both cohorts itself: unset PBWPFVC_COHORT")
site_name <- config$site_name
final_dir <- final_dir_for("supplement")
source(here::here("code", "20_biotrauma_grid.R"))   # pfvc_channels()

MIN_EVENTS <- 10L               # minimum deaths per model: the CLIF minimum-count standard
SAME_DAY_DEATH_D <- 0.5         # script 03's day for a death stamped between admission and the index
HORIZON_DAYS <- 60
SF_HYPOXEMIA_THRESHOLD <- 315   # script 03's index gate for the ventilated cohort
PLACEBO_N <- 500L
PLACEBO_SEED <- 20260930
COHORTS <- c("Ventilated", "No support")
DIFFERENCE <- "ventilated minus no support"

# =============================================================================
# Data: both cross-sectional cohorts, as xsec_mortality_channel_equality.R builds them
# =============================================================================
cohort_columns <- c("hospitalization_id", "index_dttm", "age_at_admission", "sex_category",
                    "race_category", "pfvc", "pbw", "sf_ratio", "sofa_total", "deceased", "admission_dttm", "death_dttm",
                    "discharge_dttm", "height_cm")
control_file <- file.path(config$output_dir, "controls", "nosupport", "analysis_cross_sectional.parquet")
if (!file.exists(control_file))
  stop("no no-support cohort: run scripts 01-03 with PBWPFVC_COHORT=nosupport first")
# the ventilated arm: patients on invasive ventilation at ICU admission (icu_day0)
ventilated_cohort <- read_parquet(file.path(config$output_dir, "analysis_cross_sectional.parquet")) %>%
  select(all_of(cohort_columns), vtpbw, icu_day0)
ventilated <- ventilated_cohort %>% filter(icu_day0) %>% select(-icu_day0) %>%
  mutate(cohort = "Ventilated", escalation_dttm = as.POSIXct(NA))
message("Ventilated arm: ", nrow(ventilated), " of ", nrow(ventilated_cohort),
        " ventilated-cohort patients on invasive ventilation at ICU admission (icu_day0)")
no_support <- read_parquet(control_file) %>%
  select(all_of(cohort_columns), escalation_dttm) %>%
  mutate(cohort = "No support", vtpbw = NA_real_)
# the control's first invasive ventilation after the index (device "imv" or a set
# tidal volume, as 03's escalation rule)
control_imv <- read_parquet(file.path(config$output_dir, "controls", "nosupport", "resp_support_waterfall_clean.parquet"),
                            col_select = c("hospitalization_id", "recorded_dttm", "device_category", "tidal_volume_set")) %>%
  filter(tolower(device_category) == "imv" | (!is.na(tidal_volume_set) & tidal_volume_set > 0)) %>%
  inner_join(no_support %>% select(hospitalization_id, index_dttm), by = "hospitalization_id") %>%
  filter(recorded_dttm >= index_dttm) %>%
  group_by(hospitalization_id) %>% summarise(imv_dttm = min(recorded_dttm), .groups = "drop")
no_support <- no_support %>% left_join(control_imv, by = "hospitalization_id")
both_cohorts <- bind_rows(ventilated, no_support)

# SYNTHETIC SITE ONLY: synthetic CLIF mortality is unreliable, so death is simulated
# independently of every exposure (35% by day 60, time from the index log-normal
# with median 9 days; a simulated death is in hospital), with the seed and row order
# of xsec_pfvc_age_control.R. The run exercises the machinery and can show no real
# effect. Never runs at a real site.
if (grepl("^synthetic_clif", site_name)) {
  message("*** SYNTHETIC SITE: simulated mortality (plumbing only; synthetic CLIF mortality is unreliable). ***")
  set.seed(20260615)
  simulated_death <- rbinom(nrow(both_cohorts), 1L, 0.35)
  simulated_day   <- pmin(pmax(rlnorm(nrow(both_cohorts), log(9), 0.95), 0.04), 60)
  both_cohorts <- both_cohorts %>%
    mutate(deceased = simulated_death,
           death_dttm = if_else(simulated_death == 1L, index_dttm + simulated_day * 86400, as.POSIXct(NA)),
           discharge_dttm = if_else(simulated_death == 1L, death_dttm, pmax(discharge_dttm, index_dttm)))
}
# script 03 dates an expired patient with no death time at discharge
if (any(both_cohorts$deceased == 1 & is.na(both_cohorts$death_dttm), na.rm = TRUE))
  stop("expired patients without a death time: rebuild the cohorts with script 03, which dates them at discharge")

index_day <- function(dttm, index) as.numeric(difftime(dttm, index, units = "days"))
# Patients missing PFVC, PBW, SF, SOFA, the outcome, a demographic, height or the
# discharge time are excluded; the counts before and after are printed per cohort.
n_before_missing_data <- both_cohorts %>% count(cohort, name = "n_before")
both_cohorts <- both_cohorts %>%
  filter(!is.na(pfvc), pfvc > 0, !is.na(pbw), pbw > 0, !is.na(sf_ratio), !is.na(sofa_total), !is.na(deceased),
         !is.na(age_at_admission), !is.na(sex_category), !is.na(race_category), !is.na(discharge_dttm),
         !is.na(height_cm), height_cm > 0) %>%
  group_by(cohort) %>%
  mutate(sf_z = as.numeric(scale(sf_ratio)), sofa_z = as.numeric(scale(sofa_total))) %>%
  ungroup() %>%
  mutate(cohort        = factor(cohort, levels = COHORTS),
         sex_category  = factor(sex_category, levels = c("Male", "Female")),
         race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
         log_pfvc      = log(pfvc),
         log_ratio     = log(pbw / pfvc),                  # log PBW/PFVC: the strain error at a given VT/PBW
         log_vtpfvc    = log(vtpbw) + log_ratio,           # log VT/PFVC up to a constant (ventilated only)
         death_stamped_before_index = !is.na(death_dttm) & death_dttm < index_dttm,
         death_index_day     = if_else(death_stamped_before_index, SAME_DAY_DEATH_D, index_day(death_dttm, index_dttm)),
         discharge_index_day = index_day(discharge_dttm, index_dttm),
         imv_day             = index_day(imv_dttm, index_dttm),   # control only; NA in the ventilated cohort
         age10 = age_at_admission / 10,
         hypoxemic = sf_ratio < SF_HYPOXEMIA_THRESHOLD)
missing_data_counts <- n_before_missing_data %>%
  left_join(both_cohorts %>% count(cohort = as.character(cohort), name = "n_after"), by = "cohort") %>%
  mutate(n_after = coalesce(n_after, 0L), n_patients_excluded_missing_data = n_before - n_after)
message("Excluded for missing PFVC, PBW, SF, SOFA, outcome, demographics, height or discharge time:")
print(as.data.frame(missing_data_counts), row.names = FALSE)
if (anyNA(both_cohorts$vtpbw[both_cohorts$cohort == "Ventilated"]))
  stop("ventilated patients without VT/PBW in the cross-sectional table")

# A model that warns stops the script; separation (a coefficient running to infinity)
# is caught and reported as not estimable.
SEPARATION_PATTERN <- "coefficient may be infinite|fitted probabilities numerically 0 or 1"
fit_strict <- function(expr) withCallingHandlers(expr, warning = function(w) {
  message_text <- conditionMessage(w)
  if (grepl(SEPARATION_PATTERN, message_text))
    stop(structure(class = c("separation", "error", "condition"),
                   list(message = paste("not estimable (separation):", trimws(message_text)), call = NULL)))
  stop("model warning, stopping: ", message_text, call. = FALSE)
})
fit_or_separation <- function(expr) tryCatch(expr, separation = function(condition) condition)
separated <- function(fit) inherits(fit, "separation")

# time (days from the index) and event for 60-day death before invasive ventilation,
# as xsec_pfvc_age_control.R defines it (identical to 60-day death in the ventilated
# arm, which is intubated at the index)
outcome_data <- function(cohort_data) cohort_data %>% mutate(
  censor_day = pmin(HORIZON_DAYS, coalesce(imv_day, Inf)),
  event = as.integer(!is.na(death_index_day) & death_index_day <= censor_day),
  end_day = pmax(if_else(event == 1L, death_index_day, censor_day), 0.01))
OUTCOMES <- c(inhosp_logistic = "in-hospital death (logistic)",
              day60_before_imv = "60-day death, before invasive ventilation (Cox)")

# One coefficient: `term` in the model event ~ SF + SOFA + [extra] + term
fit_term <- function(dat, outcome_key, term, extra = character(0)) {
  dat <- if (outcome_key == "inhosp_logistic") dat %>% mutate(event = deceased) else outcome_data(dat)
  row_head <- tibble(n_patients = nrow(dat), n_deaths = sum(dat$event == 1))
  if (row_head$n_deaths < MIN_EVENTS)
    return(row_head %>% mutate(log_ratio = NA_real_, se = NA_real_, note = paste("skipped: fewer than", MIN_EVENTS, "deaths")))
  rhs <- paste(c("sf_z", "sofa_z", extra, term), collapse = " + ")
  fit <- fit_or_separation(if (outcome_key == "inhosp_logistic")
    fit_strict(glm(as.formula(paste("event ~", rhs)), family = binomial, data = dat)) else
    fit_strict(coxph(as.formula(paste("Surv(end_day, event) ~", rhs)), data = dat)))
  if (separated(fit)) return(row_head %>% mutate(log_ratio = NA_real_, se = NA_real_, note = conditionMessage(fit)))
  row_head %>% mutate(log_ratio = unname(coef(fit)[term]), se = unname(sqrt(vcov(fit)[term, term])), note = NA_character_)
}

# =============================================================================
# 1. The site's dosing rule, before the VT/PBW gate
# =============================================================================
pre_gate <- read_parquet(file.path(config$output_dir, "analysis_all_eligible_timepoints.parquet"),
                         col_select = c("hospitalization_id", "recorded_dttm", "has_all_data", "tidal_volume_set",
                                        "pbw", "pfvc", "vtpbw", "sex_category")) %>%
  filter(has_all_data, !is.na(tidal_volume_set), tidal_volume_set > 0, !is.na(pbw), pbw > 0, !is.na(pfvc), pfvc > 0,
         sex_category %in% c("Male", "Female")) %>%
  group_by(hospitalization_id) %>% slice_min(recorded_dttm, n = 1, with_ties = FALSE) %>% ungroup() %>%
  mutate(sex_category = factor(sex_category, levels = c("Male", "Female")))
dosing_rule <- function(dat, population) {
  pooled <- lm(log(tidal_volume_set) ~ log(pbw), data = dat)
  within_sex <- lm(log(tidal_volume_set) ~ log(pbw) + sex_category, data = dat)
  on_pfvc <- lm(log(tidal_volume_set) ~ log(pfvc), data = dat)
  slope_row <- function(fit, term) c(summary(fit)$coefficients[term, 1:2], r2 = summary(fit)$r.squared)
  s_pooled <- slope_row(pooled, "log(pbw)"); s_sex <- slope_row(within_sex, "log(pbw)"); s_pfvc <- slope_row(on_pfvc, "log(pfvc)")
  tibble(population = population, n_patients = nrow(dat),
         slope_vt_on_pbw = s_pooled[1], slope_vt_on_pbw_se = s_pooled[2], r2_vt_on_pbw = s_pooled[3],
         slope_vt_on_pbw_within_sex = s_sex[1], slope_vt_on_pbw_within_sex_se = s_sex[2],
         slope_vt_on_pfvc = s_pfvc[1], slope_vt_on_pfvc_se = s_pfvc[2], r2_vt_on_pfvc = s_pfvc[3],
         median_vt_ml = median(dat$tidal_volume_set), sd_vt_ml = sd(dat$tidal_volume_set),
         median_vtpbw = median(dat$vtpbw), sd_vtpbw = sd(dat$vtpbw),
         share_vtpbw_6_to_8 = mean(dat$vtpbw >= 6 & dat$vtpbw <= 8))
}
dosing <- bind_rows(
  dosing_rule(pre_gate, "first complete IMV row, before the VT/PBW 6-8 gate (the site's rule)"),
  dosing_rule(pre_gate %>% filter(vtpbw >= 6, vtpbw <= 8),
              "the same rows inside the 6-8 gate (for reference: the gate pushes the slope toward 1)")) %>%
  mutate(site = site_name)
message("\nThe site's dosing rule (slope of log set VT on log PBW: 1 = per-kg dosing, 0 = a fixed volume):")
print(as.data.frame(dosing %>% select(population, n_patients, slope_vt_on_pbw, slope_vt_on_pbw_within_sex, r2_vt_on_pbw,
                                      sd_vt_ml, sd_vtpbw, share_vtpbw_6_to_8) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)

# =============================================================================
# 2. The ventilation contrast
# =============================================================================
POPULATIONS <- c(everyone = "everyone", hypoxemic = "hypoxemic at the index (SF < 315)")
population_data <- function(population, cohort_now) both_cohorts %>%
  filter(cohort == cohort_now, if (population == "hypoxemic") hypoxemic else TRUE)
EXPOSURES <- c(log_ratio = "log PBW/PFVC", log_pfvc = "log PFVC")
DOSES <- c(no_dose = "no dose", fixed_vtpbw = "at fixed VT/PBW")
contrast_grid <- expand_grid(population = names(POPULATIONS), outcome_key = names(OUTCOMES), exposure = names(EXPOSURES))
per_cohort <- contrast_grid %>%
  mutate(fits = pmap(list(population, outcome_key, exposure), function(population, outcome_key, exposure) {
    control <- fit_term(population_data(population, "No support"), outcome_key, exposure) %>% mutate(quantity = "No support")
    ventilated <- map_dfr(names(DOSES), function(dose)
      fit_term(population_data(population, "Ventilated"), outcome_key, exposure,
               extra = if (dose == "fixed_vtpbw") "vtpbw" else character(0)) %>%
        mutate(quantity = "Ventilated", dose = dose))
    differences <- ventilated %>% transmute(dose, quantity = DIFFERENCE, log_ratio = log_ratio - control$log_ratio,
                                            se = sqrt(se^2 + control$se^2), n_patients = NA_integer_, n_deaths = NA_integer_,
                                            note = coalesce(note, control$note))
    bind_rows(ventilated, control %>% mutate(dose = "no dose (a control has no tidal volume)"), differences)
  })) %>% unnest(fits)
# the delivered strain, ventilated only: log VT/PFVC, severity only
delivered <- expand_grid(population = names(POPULATIONS), outcome_key = names(OUTCOMES)) %>%
  mutate(fits = map2(population, outcome_key, ~ fit_term(population_data(.x, "Ventilated"), .y, "log_vtpfvc"))) %>%
  unnest(fits) %>% mutate(exposure = "log_vtpfvc", quantity = "Ventilated", dose = "delivered strain")
contrast <- bind_rows(per_cohort, delivered) %>%
  mutate(population = POPULATIONS[population], outcome = OUTCOMES[outcome_key],
         exposure = c(EXPOSURES, log_vtpfvc = "log VT/PFVC (delivered strain)")[exposure],
         dose = coalesce(DOSES[dose], dose),
         ratio_per_0.1 = exp(0.1 * log_ratio), lo_per_0.1 = exp(0.1 * (log_ratio - 1.96 * se)),
         hi_per_0.1 = exp(0.1 * (log_ratio + 1.96 * se)), p = 2 * pnorm(-abs(log_ratio / se)),
         scale = "log OR or log HR per log unit of the exposure; ratio_per_0.1 is per 0.1 log units",
         site = site_name) %>%
  select(population, outcome, exposure, dose, quantity, log_ratio, se, ratio_per_0.1, lo_per_0.1, hi_per_0.1, p,
         n_patients, n_deaths, note, scale, site)
message("\nThe ventilation contrast, in-hospital death, everyone (per 0.1 log units; severity only):")
print(as.data.frame(contrast %>% filter(population == "everyone", grepl("^in-hospital", outcome)) %>%
                      select(exposure, dose, quantity, ratio_per_0.1, lo_per_0.1, hi_per_0.1, p, n_deaths) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)

# =============================================================================
# 3. Placebo formulas
# =============================================================================
# The inputs, whitened over both cohorts together: X centred, then rotated and scaled
# by its covariance's eigen-decomposition, so Z has identity covariance and any unit
# vector w gives an index Z w with SD 1.
input_matrix <- cbind(ns(both_cohorts$age_at_admission, df = 4),
                      female = as.numeric(both_cohorts$sex_category == "Female"),
                      black  = as.numeric(both_cohorts$race_category == "BLACK"),
                      other  = as.numeric(both_cohorts$race_category == "OTHER"),
                      log_height = log(both_cohorts$height_cm))
colnames(input_matrix)[1:4] <- paste0("age_spline_", 1:4)
input_centred <- scale(input_matrix, center = TRUE, scale = FALSE)
eigen_inputs <- eigen(cov(input_centred), symmetric = TRUE)
if (min(eigen_inputs$values) < 1e-10 * max(eigen_inputs$values))
  stop("the placebo inputs are collinear at this site (a race or sex group is empty?): eigenvalues ",
       paste(signif(eigen_inputs$values, 3), collapse = ", "))
whitened <- input_centred %*% eigen_inputs$vectors %*% diag(1 / sqrt(eigen_inputs$values))
standardise <- function(x) (x - mean(x)) / sd(x)
gli_pieces <- pfvc_channels(both_cohorts, "log_pfvc")
ratio_projection <- lm(both_cohorts$log_ratio ~ input_matrix)
reference_indices <- list(
  `log PBW/PFVC (GLI and Devine)`       = both_cohorts$log_ratio,
  `log PFVC (GLI)`                      = both_cohorts$log_pfvc,
  `log PBW/PFVC projected on the inputs` = fitted(ratio_projection),
  `GLI height piece`                    = gli_pieces$ch_height,
  `GLI age piece`                       = gli_pieces$ch_age,
  `GLI sex piece`                       = gli_pieces$ch_sex,
  `GLI race piece`                      = gli_pieces$ch_race) %>% map(standardise)
set.seed(PLACEBO_SEED)
placebo_directions <- matrix(rnorm(PLACEBO_N * ncol(whitened)), nrow = ncol(whitened))
placebo_directions <- sweep(placebo_directions, 2, sqrt(colSums(placebo_directions^2)), "/")   # unit vectors
placebo_indices <- map(seq_len(PLACEBO_N), ~ as.vector(whitened %*% placebo_directions[, .x]))
names(placebo_indices) <- sprintf("placebo %03d", seq_len(PLACEBO_N))
# the in-hospital contrast (everyone, no dose) of one index
index_contrast <- function(index_values) {
  dat <- both_cohorts %>% mutate(index = index_values)
  fits <- map(set_names(COHORTS), ~ fit_term(dat %>% filter(cohort == .x), "inhosp_logistic", "index"))
  tibble(ventilated = fits$Ventilated$log_ratio, ventilated_se = fits$Ventilated$se,
         control = fits$`No support`$log_ratio, control_se = fits$`No support`$se,
         note = coalesce(fits$Ventilated$note, fits$`No support`$note)) %>%
    mutate(difference = ventilated - control, difference_se = sqrt(ventilated_se^2 + control_se^2),
           z = difference / difference_se)
}
message("\nPlacebo formulas: ", PLACEBO_N, " random directions and ", length(reference_indices), " references ...")
placebo_results <- imap_dfr(c(reference_indices, placebo_indices), ~ index_contrast(.x) %>% mutate(index = .y)) %>%
  mutate(kind = if_else(index %in% names(reference_indices), "reference", "placebo"))
placebo_abs_z <- abs(placebo_results$z[placebo_results$kind == "placebo"])
placebo_results <- placebo_results %>%
  mutate(percentile_abs_z_among_placebos = if_else(kind == "reference",
                                                   map_dbl(abs(z), ~ mean(placebo_abs_z < .x, na.rm = TRUE)), NA_real_),
         ratio_projection_r2 = summary(ratio_projection)$r.squared,
         scale = "log OR of in-hospital death per SD of the index (SD over both cohorts); severity only",
         site = site_name) %>%
  relocate(kind, index)
message("\nReference indices against ", PLACEBO_N, " placebo directions (in-hospital, everyone, severity only;",
        " per SD of the index; the ratio's projection on the inputs has R2 ", signif(summary(ratio_projection)$r.squared, 3), "):")
print(as.data.frame(placebo_results %>% filter(kind == "reference") %>%
                      select(index, ventilated, control, difference, difference_se, z, percentile_abs_z_among_placebos) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
message("Placebo |z|: median ", signif(median(placebo_abs_z, na.rm = TRUE), 3), ", 95th percentile ",
        signif(quantile(placebo_abs_z, 0.95, na.rm = TRUE), 3))

write_csv(dosing, file.path(final_dir, paste0("strain_invariance_dosing_", site_name, ".csv")))
write_csv(contrast, file.path(final_dir, paste0("strain_invariance_contrast_", site_name, ".csv")))
write_csv(placebo_results, file.path(final_dir, paste0("strain_invariance_placebo_", site_name, ".csv")))
message("\nWrote strain_invariance_{dosing,contrast,placebo}_", site_name, ".csv to ", final_dir)
