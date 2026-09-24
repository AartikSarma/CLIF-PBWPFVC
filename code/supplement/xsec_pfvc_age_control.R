# =============================================================================
# Supplement (cross-sectional): does predicted lung size carry age's mortality
# gradient under ventilation only?
# =============================================================================
# In the ventilated cohort, adding log PFVC to a model with a 4-df age spline, sex
# and race flattens the age curve of in-hospital death by about half between 40
# and 90 (MIMIC, 2026-09-23). That split between "through predicted lung size"
# and "not" rests on an assumption the ventilated cohort cannot check. With the
# demographics in, the PFVC coefficient is identified mostly by height, and
# crediting part of age's gradient to lung size assumes a log unit of PFVC does
# the same whether height or age moved it. A direct effect of height (body size,
# socioeconomic position) or of predicted lung reserve would produce the same
# flattening without any ventilator.
#
# The no-support control separates the two. Its patients have the same formulas,
# the same demographic paths to death and no PBW-scaled tidal volume, so strain
# cannot act. Strain predicts a protective PFVC association in the ventilated
# cohort and an attenuated or absent one in the control; lung reserve or a direct
# height effect predicts the same association in both. The absorption of the age
# gradient is the PFVC coefficient times GLI's age slope, which the formula fixes,
# so the formal test is the difference in the PFVC coefficient between cohorts:
# the cohort x log PFVC interaction in one model where every covariate is free by
# cohort.
#
# Model, per cohort (in-hospital death, logistic):
#   deceased ~ [VT/PBW] + SF (z) + SOFA (z) + ns(age, 4) + sex + race  [+ log PFVC]
# VT/PBW, the delivered dose, only in the ventilated cohort. SF and SOFA are
# standardised within cohort: the control's SF uses an FiO2 estimated from device
# and flow, and the ventilated cohort is gated at SF < 315, so severity overlaps
# little and adjusts within each population rather than matching across them.
# Age stays in every model, because its curve is the object. The unadjusted arm
# drops sex and race only.
#
# The age curve is the model's log-odds of death for a white man at the cohort's
# median dose, severity and log PFVC, relative to age 40, with and without log
# PFVC. The share of the 40-to-90 rise that log PFVC absorbs is reported per
# cohort, as a description; the interaction is the inference. The share is a ratio
# of two rises, so it is unstable when the rise without PFVC is small (read the two
# rises beside it).
#
# The contrast is read across an outcome grid and two severity forms (described
# at the contrast section below). The grid separates two explanations of MIMIC's
# first run, where the control was null for in-hospital death but protective for
# 60-day death before escalation: deaths after discharge (a reserve channel outside
# the ventilator) and informative censoring at escalation. MIMIC's second run
# pointed at escalation (the control turned protective only when escalation
# censored), so two reads follow: the control's deaths counted up to intubation
# rather than up to any support, and the PFVC association with escalation itself,
# which says whether censoring at escalation selects on PFVC. Script 03 removes
# control patients escalated within 24 hours of the index; those escalated later
# stay in, and their deaths can follow ventilation. The severity form reads the
# control at the ventilated cohort's severity, as figure 4 does, because the
# control is much less sick and a smaller lung may show only under stress.
#
# Log PFVC is per SD of the ventilated cohort, so both cohorts share one unit. A
# model that warns (non-convergence, separation) stops the script: no model is
# silently replaced by a simpler one.
#
# Inputs : intermediate/analysis_cross_sectional.parquet (script 03, ventilated)
#          intermediate/controls/nosupport/analysis_cross_sectional.parquet
#          (script 03 run with PBWPFVC_COHORT=nosupport)
#          intermediate/controls/nosupport/resp_support_waterfall_clean.parquet
#          clif_code_status and clif_hospitalization (config$tables_path; optional)
# Outputs: final/supplement/
#   pfvc_age_control_estimates_{site}.csv  per cohort and adjustment: PFVC OR per SD,
#                                          likelihood-ratio p, share of the age
#                                          gradient absorbed, counts
#   pfvc_age_control_curves_{site}.csv     the age curves with and without log PFVC
#   pfvc_age_control_contrast_{site}.csv   PFVC ratio per cohort and their difference,
#                                          by outcome, severity form and adjustment,
#                                          with the PFVC x anchor terms
#   pfvc_age_control_anchor_{site}.csv     the severity anchor by cohort: how many
#                                          controls reach the ventilated mean
#   pfvc_age_control_escalation_paths_{site}.csv   the control's patients and deaths by
#                                          escalation path (never, noninvasive only,
#                                          invasive ventilation)
#   pfvc_age_control_escalation_hazard_{site}.csv  the control's escalation hazard per
#                                          SD of log PFVC (any support; invasive)
#   pfvc_age_control_channels_{site}.csv   the channel control contrast: each GLI piece's
#                                          coefficient per cohort and ventilated minus
#                                          no support; pfvc_age_control_channel_tests_
#                                          {site}.csv, whether the differences agree;
#                                          pfvc_age_control_channel_conversion_{site}.csv,
#                                          how each piece rescales between the two
#                                          exposure scales. Every channel table carries
#                                          both scales (column exposure): log PFVC
#                                          (predicted size) and log PBW/PFVC (strain
#                                          error), with each piece's identifying SD;
#                                          pfvc_age_control_channel_vcov_{site}.csv, the
#                                          differences' covariance, for the pooled test
#   pfvc_age_control_code_status_{site}.csv  patients and deaths by code status at the
#                                          index and later limitation, per cohort
#   pfvc_age_control_hypoxemia_{site}.csv  the hypoxemia pathway in the control: onset
#                                          of hypoxemia by PFVC, and PFVC's death HR
#                                          before and after hypoxemia (needs the 7-day
#                                          control panel of 21_biotrauma_panel.R)
# The contrast and the escalation hazard are fitted in five populations (column
# population): everyone; full code at the index; full code throughout; hypoxemic at
# the index (SF < 315, the ventilated cohort's own gate); hypoxemic and full code.
# The full-code ones need the optional CLIF code_status table (section "Code status").
#   pfvc_age_control_{site}.pdf
# Usage: uvr run code/supplement/xsec_pfvc_age_control.R   (PBWPFVC_COHORT unset)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(survival)
  library(splines)
  library(patchwork)
})

options(width = 220)
source("utils/config.R")
if (config$cohort != "imv") stop("xsec_pfvc_age_control.R reads both cohorts itself: unset PBWPFVC_COHORT")
site_name <- config$site_name
final_dir <- final_dir_for("supplement")

AGE_CURVE_GRID <- seq(20, 90, by = 10)
MIN_DEATHS <- 10L
SF_HYPOXEMIA_THRESHOLD <- 315   # script 03's index gate for the ventilated cohort
COHORTS <- c("Ventilated", "No support")
OKABE_ITO <- c(without_pfvc = "#E69F00", with_pfvc = "#0072B2")

# =============================================================================
# Data: both cross-sectional cohorts, on shared columns
# =============================================================================
cohort_columns <- c("hospitalization_id", "recorded_dttm", "age_at_admission", "sex_category",
                    "race_category", "pfvc", "sf_ratio", "sofa_total",
                    "sofa_cv_97", "sofa_coag", "sofa_liver", "sofa_renal",
                    "deceased", "death_dttm", "discharge_dttm", "height_cm")
control_file <- file.path(config$output_dir, "controls", "nosupport", "analysis_cross_sectional.parquet")
if (!file.exists(control_file))
  stop("no no-support cohort: run scripts 01-03 with PBWPFVC_COHORT=nosupport first")
ventilated <- read_parquet(file.path(config$output_dir, "analysis_cross_sectional.parquet")) %>%
  select(all_of(cohort_columns), vtpbw) %>%
  mutate(cohort = "Ventilated", escalation_dttm = as.POSIXct(NA))
no_support <- read_parquet(control_file) %>%
  select(all_of(cohort_columns), escalation_dttm) %>%
  mutate(cohort = "No support", vtpbw = NA_real_)
# The control's first invasive ventilation after the index, from its own respiratory
# support table: escalation_dttm (script 03) is the first advanced support of any
# kind, and only invasive ventilation delivers a PBW-scaled tidal volume. A patient
# escalated to high-flow or noninvasive ventilation who is intubated later counts
# from the intubation. The same rule as 03's escalation: device "imv" or a set tidal
# volume.
control_imv <- read_parquet(file.path(config$output_dir, "controls", "nosupport", "resp_support_waterfall_clean.parquet"),
                            col_select = c("hospitalization_id", "recorded_dttm", "device_category", "tidal_volume_set")) %>%
  filter(tolower(device_category) == "imv" | (!is.na(tidal_volume_set) & tidal_volume_set > 0)) %>%
  inner_join(no_support %>% select(hospitalization_id, index_dttm = recorded_dttm), by = "hospitalization_id") %>%
  filter(recorded_dttm >= index_dttm) %>%
  group_by(hospitalization_id) %>% summarise(imv_dttm = min(recorded_dttm), .groups = "drop")
no_support <- no_support %>% left_join(control_imv, by = "hospitalization_id")
both_cohorts <- bind_rows(ventilated, no_support)

# SYNTHETIC SITE ONLY: synthetic CLIF mortality is unreliable, so death is simulated
# independently of every exposure (35% by day 60, time from the index log-normal
# with median 9 days; a simulated death is in hospital), as the other supplement
# scripts do. The run exercises the machinery and can show no real effect. Never
# runs at a real site.
if (grepl("^synthetic_clif", site_name)) {
  message("*** SYNTHETIC SITE: simulated mortality (plumbing only; synthetic CLIF mortality is unreliable). ***")
  set.seed(20260615)
  simulated_death <- rbinom(nrow(both_cohorts), 1L, 0.35)
  simulated_day   <- pmin(pmax(rlnorm(nrow(both_cohorts), log(9), 0.95), 0.04), 60)
  both_cohorts <- both_cohorts %>%
    mutate(deceased = simulated_death,
           death_dttm = if_else(simulated_death == 1L, recorded_dttm + simulated_day * 86400, as.POSIXct(NA)),
           discharge_dttm = if_else(simulated_death == 1L, death_dttm, pmax(discharge_dttm, recorded_dttm)))
}

HORIZON_DAYS <- 60
index_day <- function(dttm, index) as.numeric(difftime(dttm, index, units = "days"))
ventilated_log_pfvc_sd <- sd(log(ventilated$pfvc[ventilated$pfvc > 0]), na.rm = TRUE)
both_cohorts <- both_cohorts %>%
  filter(!is.na(pfvc), pfvc > 0, !is.na(sf_ratio), !is.na(sofa_total), !is.na(deceased),
         !is.na(age_at_admission), !is.na(sex_category), !is.na(race_category), !is.na(discharge_dttm),
         !is.na(height_cm), height_cm > 0) %>%
  group_by(cohort) %>%
  mutate(sf_z = as.numeric(scale(sf_ratio)), sofa_z = as.numeric(scale(sofa_total))) %>%
  ungroup() %>%
  mutate(cohort        = factor(cohort, levels = COHORTS),
         sex_category  = factor(sex_category, levels = c("Male", "Female")),
         race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
         log_pfvc_z    = log(pfvc) / ventilated_log_pfvc_sd,
         # every clock starts at the index; 03's surv_time runs from hospital admission
         death_index_day     = index_day(death_dttm, recorded_dttm),
         discharge_index_day = index_day(discharge_dttm, recorded_dttm),
         escalation_day      = index_day(escalation_dttm, recorded_dttm),
         imv_day             = index_day(imv_dttm, recorded_dttm),   # control only; NA in the ventilated cohort
         # figure 4's severity anchor: SOFA without its respiratory and neurological parts
         anchor = sofa_cv_97 + sofa_coag + sofa_liver + sofa_renal,
         age10 = age_at_admission / 10)
if (anyNA(both_cohorts$vtpbw[both_cohorts$cohort == "Ventilated"]))
  stop("ventilated patients without VT/PBW in the cross-sectional table")
if (anyNA(both_cohorts$anchor)) stop("patients without the SOFA components of the severity anchor")
# Deaths timestamped before the index are excluded (user, 2026-09-23). A death
# cannot precede the index, so either the death time or the index is wrong for
# these patients (MIMIC: 14 of about 19,000, most of them in-hospital deaths within
# a day before a ventilated index), and neither can be repaired from here. The
# breakdown is printed and the counts excluded per cohort are written with the
# estimates.
death_before_index <- both_cohorts %>% filter(!is.na(death_index_day), death_index_day < 0) %>%
  mutate(how_far = if_else(death_index_day >= -1, "within 1 day before the index", "more than 1 day before the index")) %>%
  count(cohort, how_far, in_hospital_death = deceased == 1, name = "n_patients")
if (nrow(death_before_index)) {
  message("Excluded, death timestamped before the index:")
  print(as.data.frame(death_before_index), row.names = FALSE)
}
excluded_by_cohort <- both_cohorts %>% group_by(cohort) %>%
  summarise(n_patients_excluded_death_before_index = sum(!is.na(death_index_day) & death_index_day < 0), .groups = "drop")
both_cohorts <- both_cohorts %>% filter(is.na(death_index_day) | death_index_day >= 0)
ventilated_mean_anchor <- mean(both_cohorts$anchor[both_cohorts$cohort == "Ventilated"])
both_cohorts <- both_cohorts %>% mutate(anchor_c = anchor - ventilated_mean_anchor)

cohort_counts <- both_cohorts %>% group_by(cohort) %>%
  summarise(n_patients = n(), n_deaths = sum(deceased == 1),
            n_escalated = sum(!is.na(escalation_day)), median_age = median(age_at_admission),
            median_sf = median(sf_ratio), median_sofa = median(sofa_total),
            n_deaths_60d = sum(!is.na(death_index_day) & death_index_day <= HORIZON_DAYS),
            n_deaths_after_discharge_60d = sum(deceased == 0 & !is.na(death_index_day) & death_index_day <= HORIZON_DAYS),
            .groups = "drop")
print(as.data.frame(cohort_counts), row.names = FALSE)
if (any(cohort_counts$n_deaths < MIN_DEATHS) || nrow(cohort_counts) < 2)
  stop("each cohort needs at least ", MIN_DEATHS, " in-hospital deaths")

# =============================================================================
# Code status: the full-code populations
# =============================================================================
# MIMIC's third run left one explanation for the control's protective PFVC
# association before intubation: goals of care. A control patient with a smaller
# predicted lung (at fixed age, sex and race: shorter, perhaps frailer) who
# deteriorates may die without intubation under a do-not-intubate order, while a
# larger one is intubated. Restricting both cohorts to full-code patients removes
# that path. Two populations, beside everyone:
#   full code at the index   the last code status recorded up to 24 hours after
#                            the index is Full or Presume Full (orders are often
#                            written hours after ICU admission). A limitation
#                            entered later, during deterioration, is not removed.
#   full code throughout     full code at the index AND no other status recorded
#                            before death, discharge or day 60. This removes the
#                            later limitations too, but it conditions on the future:
#                            limitations are often written as death approaches, so
#                            it drops many deaths and selects survivors. A reading
#                            aid for the first population, not a replacement.
# Code status is a patient-level CLIF table (clif_code_status: patient_id,
# start_dttm, code_status_category), mapped to hospitalizations through
# clif_hospitalization. It is optional in CLIF: without it the full-code
# populations are skipped, announced, and only everyone is analysed.
CODE_STATUS_WINDOW_H <- 24
FULL_CODE_CATEGORIES <- c("full", "presume full")
read_clif_table <- function(table_name, columns) {
  path <- file.path(path.expand(config$tables_path), paste0("clif_", table_name, ".", config$file_type))
  switch(config$file_type,
         parquet = read_parquet(path, col_select = all_of(columns)),
         csv     = readr::read_csv(path, col_select = all_of(columns), show_col_types = FALSE),
         fst     = fst::read_fst(path, columns = columns))
}
code_status_file <- file.path(path.expand(config$tables_path), paste0("clif_code_status.", config$file_type))
HAS_CODE_STATUS <- file.exists(code_status_file)
if (HAS_CODE_STATUS) {
  code_status <- read_clif_table("code_status", c("patient_id", "start_dttm", "code_status_category")) %>%
    inner_join(read_clif_table("hospitalization", c("patient_id", "hospitalization_id")) %>%
                 filter(hospitalization_id %in% both_cohorts$hospitalization_id),
               by = "patient_id", relationship = "many-to-many") %>%
    # keyed by cohort too: a hospitalization can hold a no-support index and, later,
    # a ventilated one
    inner_join(both_cohorts %>% transmute(cohort, hospitalization_id, index_dttm = recorded_dttm,
                                          end_dttm = pmin(death_dttm, discharge_dttm,
                                                          recorded_dttm + HORIZON_DAYS * 86400, na.rm = TRUE)),
               by = "hospitalization_id", relationship = "many-to-many") %>%
    mutate(is_full = tolower(code_status_category) %in% FULL_CODE_CATEGORIES)
  baseline_status <- code_status %>%
    filter(start_dttm <= index_dttm + CODE_STATUS_WINDOW_H * 3600) %>%
    group_by(cohort, hospitalization_id) %>% slice_max(start_dttm, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(cohort, hospitalization_id, code_status_at_index = if_else(is_full, "full code", "limited or other"))
  limited_later <- code_status %>%
    filter(start_dttm > index_dttm + CODE_STATUS_WINDOW_H * 3600, start_dttm <= end_dttm, !is_full) %>%
    distinct(cohort, hospitalization_id) %>% mutate(limited_later = TRUE)
  both_cohorts <- both_cohorts %>%
    left_join(baseline_status, by = c("cohort", "hospitalization_id")) %>%
    left_join(limited_later, by = c("cohort", "hospitalization_id")) %>%
    mutate(code_status_at_index = coalesce(code_status_at_index, "no record"),
           limited_later = coalesce(limited_later, FALSE))
} else {
  message("*** No code_status table at ", config$tables_path, ": the full-code populations are skipped. ***")
  both_cohorts <- both_cohorts %>% mutate(code_status_at_index = "no record", limited_later = FALSE)
}
both_cohorts <- both_cohorts %>%
  mutate(population_everyone = TRUE,
         population_full_code_at_index = code_status_at_index == "full code",
         population_full_code_throughout = code_status_at_index == "full code" & !limited_later,
         # Hypoxemic at the index: the ventilated cohort's own gate (script 03 selects its
         # index at SF < 315), applied to the control, so the two arms differ in ventilation
         # and not in hypoxemia (MIMIC and UCSF, 2026-09-23: the control's median SF was 323,
         # so its contrast with the ventilated cohort also contrasted hypoxemia). The
         # control's SF uses an estimated FiO2, so its gate is not measured quite as the
         # ventilated cohort's is.
         population_hypoxemic_at_index = sf_ratio < SF_HYPOXEMIA_THRESHOLD,
         population_hypoxemic_full_code_at_index = population_hypoxemic_at_index & population_full_code_at_index)
POPULATIONS <- c(everyone = "everyone",
                 full_code_at_index = "full code at the index",
                 full_code_throughout = "full code throughout (conditions on the future)",
                 hypoxemic_at_index = "hypoxemic at the index (SF < 315)",
                 hypoxemic_full_code_at_index = "hypoxemic and full code at the index")
if (!HAS_CODE_STATUS) POPULATIONS <- POPULATIONS[c("everyone", "hypoxemic_at_index")]
n_ventilated_not_hypoxemic <- sum(both_cohorts$cohort == "Ventilated" & !both_cohorts$population_hypoxemic_at_index)
if (n_ventilated_not_hypoxemic > 0)
  message("  ", n_ventilated_not_hypoxemic, " ventilated patients have an index SF of 315 or more ",
          "(script 03 gates at SF < 315): they leave the hypoxemic populations")
code_status_counts <- both_cohorts %>%
  group_by(cohort, code_status_at_index, limited_later) %>%
  summarise(n_patients = n(), n_deaths = sum(deceased == 1), .groups = "drop") %>%
  mutate(code_status_table = HAS_CODE_STATUS, site = site_name)
message("\nCode status at the index (last status up to ", CODE_STATUS_WINDOW_H, " h after it) and later limitation:")
print(as.data.frame(code_status_counts), row.names = FALSE)

# The control's escalation paths within the horizon, and where its deaths fall: the
# deaths after escalation are the ones strain could reach, but only on the path
# through invasive ventilation
escalation_paths <- both_cohorts %>% filter(cohort == "No support") %>%
  mutate(path = case_when(
           !is.na(imv_day) & imv_day <= HORIZON_DAYS               ~ "invasive ventilation (after any noninvasive support)",
           !is.na(escalation_day) & escalation_day <= HORIZON_DAYS ~ "high-flow or noninvasive ventilation only",
           TRUE                                                    ~ "never escalated"),
         died_in_hospital = deceased == 1,
         died_after_escalation = died_in_hospital & !is.na(escalation_day) &
           coalesce(death_index_day, discharge_index_day) > escalation_day) %>%
  group_by(path) %>%
  summarise(n_patients = n(), n_deaths = sum(died_in_hospital), n_deaths_after_escalation = sum(died_after_escalation),
            median_pfvc_litres = median(pfvc), .groups = "drop") %>%
  mutate(site = site_name)
message("\nThe control's escalation paths (in-hospital deaths):")
print(as.data.frame(escalation_paths), row.names = FALSE)

# A model that warns stops the script (no silent fallback), with one exception.
# Separation (a coefficient that runs to infinity, as when a category of a small
# subgroup has no deaths; UCSF, 2026-09-23: hypoxemic controls, death before
# escalation) raises a "separation" condition instead. The subgroup fits below catch
# it and report the fit as not estimable, with the model's own message, in their
# output row; no model is replaced by a simpler one, and every other warning stops.
SEPARATION_PATTERN <- "coefficient may be infinite|fitted probabilities numerically 0 or 1"
fit_strict <- function(expr) withCallingHandlers(expr, warning = function(w) {
  message_text <- conditionMessage(w)
  if (grepl(SEPARATION_PATTERN, message_text))
    stop(structure(class = c("separation", "error", "condition"),
                   list(message = paste("not estimable (separation):", trimws(message_text)), call = NULL)))
  stop("model warning, stopping: ", message_text, call. = FALSE)
})
# the fit, or the separation condition itself (the caller writes its message as the note)
fit_or_separation <- function(expr) tryCatch(expr, separation = function(condition) condition)
separated <- function(fit) inherits(fit, "separation")

# =============================================================================
# Per cohort: the age curve with and without log PFVC
# =============================================================================
ADJUSTMENTS <- c(adjusted = "sex_category + race_category", unadjusted = NA_character_)
fit_cohort <- function(cohort_data, adjustment) {
  is_ventilated <- cohort_data$cohort[1] == "Ventilated"
  rhs <- c(if (is_ventilated) "vtpbw", "sf_z", "sofa_z", "ns(age_at_admission, 4)",
           if (!is.na(ADJUSTMENTS[[adjustment]])) ADJUSTMENTS[[adjustment]])
  without_pfvc_fit <- fit_strict(glm(as.formula(paste("deceased ~", paste(rhs, collapse = " + "))),
                                     family = binomial, data = cohort_data))
  with_pfvc_fit <- fit_strict(update(without_pfvc_fit, . ~ . + log_pfvc_z))
  age_grid <- tibble(age_at_admission = AGE_CURVE_GRID,
                     sex_category  = factor("Male", levels = levels(cohort_data$sex_category)),
                     race_category = factor("WHITE", levels = levels(cohort_data$race_category)),
                     vtpbw = median(cohort_data$vtpbw), sf_z = 0, sofa_z = 0,
                     log_pfvc_z = median(cohort_data$log_pfvc_z))
  age_curve <- age_grid %>%
    mutate(without_pfvc = predict(without_pfvc_fit, newdata = age_grid),
           with_pfvc    = predict(with_pfvc_fit, newdata = age_grid)) %>%
    transmute(age_at_admission,
              without_pfvc = without_pfvc - without_pfvc[age_at_admission == 40],
              with_pfvc    = with_pfvc - with_pfvc[age_at_admission == 40])
  coefficient <- summary(with_pfvc_fit)$coefficients["log_pfvc_z", ]
  rise_40_to_90 <- age_curve %>% filter(age_at_admission == 90)
  estimate_row <- tibble(
    log_or_per_sd = unname(coefficient["Estimate"]), se = unname(coefficient["Std. Error"]),
    or_per_sd = exp(log_or_per_sd), or_lo = exp(log_or_per_sd - 1.96 * se), or_hi = exp(log_or_per_sd + 1.96 * se),
    p = unname(coefficient["Pr(>|z|)"]),
    lr_p = anova(without_pfvc_fit, with_pfvc_fit, test = "LRT")$`Pr(>Chi)`[2],
    rise_40_to_90_without_pfvc = rise_40_to_90$without_pfvc, rise_40_to_90_with_pfvc = rise_40_to_90$with_pfvc,
    share_of_age_gradient_absorbed = 1 - rise_40_to_90$with_pfvc / rise_40_to_90$without_pfvc,
    n_patients = nrow(cohort_data), n_deaths = sum(cohort_data$deceased == 1))
  list(estimates = estimate_row, curves = age_curve)
}
cohort_fits <- expand_grid(cohort = COHORTS, adjustment = names(ADJUSTMENTS)) %>%
  mutate(fit = map2(cohort, adjustment, ~ fit_cohort(filter(both_cohorts, cohort == .x), .y)))
estimates <- cohort_fits %>% mutate(estimates = map(fit, "estimates")) %>% select(-fit) %>% unnest(estimates) %>%
  left_join(excluded_by_cohort %>% mutate(cohort = as.character(cohort)), by = "cohort") %>%
  mutate(scale = "OR per SD of log PFVC (ventilated cohort SD)", site = site_name)
curves <- cohort_fits %>% mutate(curves = map(fit, "curves")) %>% select(-fit) %>% unnest(curves) %>%
  mutate(scale = "log-odds of in-hospital death relative to age 40", site = site_name)

# =============================================================================
# The contrast: the PFVC coefficient, no support minus ventilated
# =============================================================================
# Each cohort is fitted on its own, which for the logistic model is the same as
# one model with every covariate free by cohort, and for the Cox model also lets
# each cohort keep its own baseline hazard. The difference is the no-support log
# ratio minus the ventilated one, with the variances of independent samples.
#
# The outcome grid. Every clock starts at the index (intubation for the ventilated
# cohort, ICU admission for the control); 03's surv_time runs from hospital
# admission and is not used here. Four cause-specific Cox models, 60-day horizon:
#   in-hospital death, all          censored at discharge
#   in-hospital death, before escalation   censored at discharge and at escalation
#   60-day death, all               includes deaths after discharge
#   60-day death, before escalation  censored at escalation
# The in-hospital logistic model of the age curves is kept beside them. If the
# control's PFVC association moves with the outcome (in-hospital against 60-day),
# deaths after discharge carry it, a reserve channel outside the ventilator; if it
# moves with escalation censoring, the censoring is informative.
#
# The severity forms. "within cohort": SF and SOFA standardised inside each cohort,
# adjustment only. "standardised to ventilated severity": as figure 4 does, the
# anchor (SOFA cardiovascular + coagulation + liver + renal, leaving out the
# respiratory component, which is computed from SF, and the neurological one, which
# on the day of intubation scores sedation) is centred at the ventilated cohort's
# mean and allowed to modify the PFVC term. The PFVC coefficient is then the
# association at the ventilated mean severity, and the PFVC x anchor term tests
# whether sicker patients show a stronger association. SF stays standardised within
# cohort, because the control's FiO2 is estimated. The control's value at the
# ventilated severity is an extrapolation where few controls are that sick: the
# anchor table says how few.
severity_rhs <- c(within = "sf_z + sofa_z", standardised = "sf_z + anchor_c + log_pfvc_z:anchor_c")
SEVERITY_LABELS <- c(within = "within cohort", standardised = "standardised to ventilated severity")
OUTCOMES <- tribble(
  ~outcome_key,          ~outcome,                                    ~model,
  "inhosp_logistic",     "in-hospital death (logistic)",               "logistic",
  "inhosp_all",          "in-hospital death, all",                     "cox",
  "inhosp_before_esc",   "in-hospital death, before escalation",       "cox",
  "inhosp_before_imv",   "in-hospital death, before invasive ventilation", "cox",
  "day60_all",           "60-day death, all",                          "cox",
  "day60_before_esc",    "60-day death, before escalation",            "cox",
  "day60_before_imv",    "60-day death, before invasive ventilation",  "cox")
outcome_data <- function(cohort_data, outcome_key) {
  # time (days from the index) and event for each Cox outcome; "before invasive
  # ventilation" keeps deaths after high-flow or noninvasive support and censors
  # only at intubation, where a PBW-scaled volume begins (identical to "all" in the
  # ventilated cohort, which is intubated at the index)
  censor_escalation <- grepl("before_esc", outcome_key)
  censor_imv <- grepl("before_imv", outcome_key)
  in_hospital <- grepl("^inhosp", outcome_key)
  # an in-hospital death is dated by death_dttm, or by discharge if that is missing;
  # discharge censors survivors only (a death's timestamp can trail its discharge)
  cohort_data %>% mutate(
    death_day = if (in_hospital) if_else(deceased == 1, coalesce(death_index_day, discharge_index_day), NA_real_)
                else death_index_day,
    censor_day = pmin(HORIZON_DAYS,
                      if (in_hospital) if_else(deceased == 1, Inf, discharge_index_day) else Inf,
                      if (censor_escalation) coalesce(escalation_day, Inf) else Inf,
                      if (censor_imv) coalesce(imv_day, Inf) else Inf),
    event = as.integer(!is.na(death_day) & death_day <= censor_day),
    end_day = pmax(if_else(event == 1L, death_day, censor_day), 0.01))
}
fit_pfvc <- function(cohort_data, outcome_key, model, adjustment, severity) {
  is_ventilated <- cohort_data$cohort[1] == "Ventilated"
  rhs <- paste(c(if (is_ventilated) "vtpbw", severity_rhs[[severity]], "ns(age_at_admission, 4)",
                 if (!is.na(ADJUSTMENTS[[adjustment]])) ADJUSTMENTS[[adjustment]], "log_pfvc_z"), collapse = " + ")
  dat <- if (model == "logistic") cohort_data %>% mutate(event = deceased) else outcome_data(cohort_data, outcome_key)
  events <- sum(dat$event == 1)
  # a restricted population can run short of deaths: the row says so and has no estimate
  if (events < MIN_DEATHS)
    return(tibble(term = "pfvc", log_ratio = NA_real_, se = NA_real_, n_patients = nrow(cohort_data),
                  n_deaths = events, note = paste("skipped: fewer than", MIN_DEATHS, "deaths")))
  fit <- fit_or_separation(if (model == "logistic")
    fit_strict(glm(as.formula(paste("deceased ~", rhs)), family = binomial, data = dat)) else
    fit_strict(coxph(as.formula(paste("Surv(end_day, event) ~", rhs)), data = dat)))
  if (separated(fit))
    return(tibble(term = "pfvc", log_ratio = NA_real_, se = NA_real_, n_patients = nrow(cohort_data),
                  n_deaths = events, note = conditionMessage(fit)))
  b <- coef(fit); V <- vcov(fit)
  interaction_term <- names(b)[sapply(strsplit(names(b), ":"), setequal, c("log_pfvc_z", "anchor_c"))]
  terms <- c(pfvc = "log_pfvc_z", pfvc_x_anchor = if (length(interaction_term)) interaction_term)
  tibble(term = names(terms), log_ratio = unname(b[terms]), se = unname(sqrt(diag(V)[terms])),
         n_patients = nrow(cohort_data), n_deaths = events, note = NA_character_)
}
population_data <- function(population, cohort_now)
  both_cohorts %>% filter(cohort == cohort_now, .data[[paste0("population_", population)]])
per_cohort <- expand_grid(population = names(POPULATIONS), cohort = COHORTS, OUTCOMES,
                          adjustment = names(ADJUSTMENTS), severity = names(severity_rhs)) %>%
  mutate(fit = pmap(list(population, cohort, outcome_key, model, adjustment, severity),
                    function(population, cohort_now, outcome_key, model, adjustment, severity)
                      fit_pfvc(population_data(population, cohort_now), outcome_key, model, adjustment, severity))) %>%
  unnest(fit)
difference <- per_cohort %>% filter(term == "pfvc") %>%
  select(population, cohort, outcome_key, adjustment, severity, log_ratio, se) %>%
  pivot_wider(names_from = cohort, values_from = c(log_ratio, se)) %>%
  transmute(population, outcome_key, adjustment, severity, cohort = "no support minus ventilated", term = "pfvc",
            log_ratio = `log_ratio_No support` - log_ratio_Ventilated,
            se = sqrt(`se_No support`^2 + se_Ventilated^2))
contrast <- bind_rows(per_cohort %>% select(-outcome, -model), difference) %>%
  left_join(OUTCOMES, by = "outcome_key") %>%
  mutate(ratio_type = if_else(model == "logistic", "OR", "HR"),
         quantity = case_when(term == "pfvc_x_anchor" ~ paste0(cohort, ": PFVC x anchor (per SOFA point)"),
                              TRUE ~ cohort),
         ratio = exp(log_ratio), ratio_lo = exp(log_ratio - 1.96 * se), ratio_hi = exp(log_ratio + 1.96 * se),
         p = 2 * pnorm(-abs(log_ratio / se)), severity = SEVERITY_LABELS[severity],
         population = POPULATIONS[population],
         scale = "per SD of log PFVC (ventilated cohort SD)",
         ventilated_log_pfvc_sd = ventilated_log_pfvc_sd, site = site_name) %>%
  select(population, outcome, ratio_type, adjustment, severity, quantity, ratio, ratio_lo, ratio_hi, p, log_ratio, se,
         n_patients, n_deaths, note, scale, ventilated_log_pfvc_sd, site)

# =============================================================================
# The control's escalation hazard by PFVC
# =============================================================================
# Censoring at escalation is uninformative about PFVC only if PFVC does not predict
# escalation. Cause-specific Cox, control only, for two escalations: any advanced
# support, and invasive ventilation. Death and discharge censor; the horizon is 60
# days. An HR below 1 says larger predicted lungs escalate less, so censoring at
# escalation removes smaller-lung patients preferentially, and the before-escalation
# death estimate is read on a population that PFVC itself selected.
fit_escalation <- function(population, escalation_column, escalation_label, adjustment, severity) {
  dat <- population_data(population, "No support") %>%
    mutate(escalation_time = .data[[escalation_column]],
           death_or_discharge = pmin(coalesce(death_index_day, Inf), discharge_index_day),
           censor_day = pmin(HORIZON_DAYS, death_or_discharge),
           event = as.integer(!is.na(escalation_time) & escalation_time <= censor_day),
           end_day = pmax(if_else(event == 1L, escalation_time, censor_day), 0.01))
  rhs <- paste(c(severity_rhs[[severity]], "ns(age_at_admission, 4)",
                 if (!is.na(ADJUSTMENTS[[adjustment]])) ADJUSTMENTS[[adjustment]], "log_pfvc_z"), collapse = " + ")
  row_head <- tibble(population = POPULATIONS[[population]], escalation = escalation_label, adjustment = adjustment,
                     severity = SEVERITY_LABELS[[severity]])
  if (sum(dat$event) < MIN_DEATHS)
    return(row_head %>% mutate(term = "pfvc", log_hr = NA_real_, se = NA_real_, n_patients = nrow(dat),
                               n_patients_escalated = sum(dat$event),
                               note = paste("skipped: fewer than", MIN_DEATHS, "escalations")))
  fit <- fit_or_separation(fit_strict(coxph(as.formula(paste("Surv(end_day, event) ~", rhs)), data = dat)))
  if (separated(fit))
    return(row_head %>% mutate(term = "pfvc", log_hr = NA_real_, se = NA_real_, n_patients = nrow(dat),
                               n_patients_escalated = sum(dat$event), note = conditionMessage(fit)))
  b <- coef(fit); se <- sqrt(diag(vcov(fit)))
  interaction_term <- names(b)[sapply(strsplit(names(b), ":"), setequal, c("log_pfvc_z", "anchor_c"))]
  terms <- c(pfvc = "log_pfvc_z", pfvc_x_anchor = if (length(interaction_term)) interaction_term)
  row_head %>% cross_join(tibble(term = names(terms), log_hr = unname(b[terms]), se = unname(se[terms]))) %>%
    mutate(n_patients = nrow(dat), n_patients_escalated = sum(dat$event), note = NA_character_)
}
escalation_hazard <- expand_grid(
  population = names(POPULATIONS),
  tibble(escalation_column = c("escalation_day", "imv_day"),
         escalation_label = c("any advanced support", "invasive ventilation")),
  adjustment = names(ADJUSTMENTS), severity = names(severity_rhs)) %>%
  pmap_dfr(fit_escalation) %>%
  mutate(hr = exp(log_hr), hr_lo = exp(log_hr - 1.96 * se), hr_hi = exp(log_hr + 1.96 * se),
         p = 2 * pnorm(-abs(log_hr / se)), cohort = "No support",
         scale = "cause-specific HR of escalation per SD of log PFVC", site = site_name)

# =============================================================================
# The hypoxemia pathway in the control
# =============================================================================
# Among control patients not hypoxemic on the index day (day-0 worst SF >= 315),
# does a smaller predicted lung lead to hypoxemia, and does its association with
# death sit after hypoxemia develops? Two fits per population (everyone, full code at
# the index), adjusted and unadjusted:
#   onset   cause-specific Cox for the first day with a worst SF < 315 (days 1-7);
#           escalation, death and discharge censor
#   death   cause-specific Cox for death before escalation, with hypoxemia as a
#           time-varying state and its interaction with log PFVC: the PFVC ratio
#           before hypoxemia, after it, and their ratio. Onset is dated at the start
#           of the day whose worst SF first falls below 315.
# Hypoxemia is observed only while the patient is unsupported and only to day 7 (the
# figure-4 panel of 21_biotrauma_panel.R, daily worst SF on an estimated FiO2): a
# patient escalated before any low SF, or hypoxemic after day 7, counts as never
# hypoxemic. The death model is therefore read at 7 days, where the state is fully
# observed, and at 60 days as a companion. Formal mediation is not identifiable:
# PFVC is fixed by the demographics, and hypoxemia is also driven by the illness.
# The panel is optional here: without it the section is skipped and announced.
HYPOXEMIA_HORIZONS <- c(7, 60)
control_panel_dir <- file.path(config$output_dir, "controls", "nosupport")
HAS_CONTROL_PANEL <- all(file.exists(file.path(control_panel_dir, c("jm_long_7d.parquet", "jm_surv_7d.parquet"))))
hypoxemia_pathway <- NULL
if (HAS_CONTROL_PANEL) {
  day0_sf <- read_parquet(file.path(control_panel_dir, "jm_surv_7d.parquet"), col_select = c("hospitalization_id", "sf_0"))
  onset <- read_parquet(file.path(control_panel_dir, "jm_long_7d.parquet"), col_select = c("hospitalization_id", "vent_day", "sf")) %>%
    filter(vent_day >= 1, !is.na(sf), sf < SF_HYPOXEMIA_THRESHOLD) %>%
    group_by(hospitalization_id) %>% summarise(hypoxemia_day = min(vent_day), .groups = "drop")
  hypoxemia_population <- function(population) population_data(population, "No support") %>%
    inner_join(day0_sf, by = "hospitalization_id") %>%
    filter(!is.na(sf_0), sf_0 >= SF_HYPOXEMIA_THRESHOLD) %>%
    left_join(onset, by = "hospitalization_id") %>%
    # onset counts only while the patient is still unsupported
    mutate(hypoxemia_day = if_else(!is.na(escalation_day) & hypoxemia_day > escalation_day, NA_real_, hypoxemia_day))
  covariate_rhs <- function(adjustment) paste(c("sf_z", "sofa_z", "ns(age_at_admission, 4)",
                                                if (!is.na(ADJUSTMENTS[[adjustment]])) ADJUSTMENTS[[adjustment]]), collapse = " + ")
  fit_onset <- function(population, adjustment) {
    dat <- hypoxemia_population(population) %>%
      mutate(censor_day = pmin(7, coalesce(escalation_day, Inf), coalesce(death_index_day, Inf), discharge_index_day),
             event = as.integer(!is.na(hypoxemia_day) & hypoxemia_day <= censor_day),
             end_day = pmax(if_else(event == 1L, hypoxemia_day, censor_day), 0.01))
    row_head <- tibble(population = POPULATIONS[[population]], adjustment, analysis = "onset of hypoxemia (days 1-7)",
                       horizon_days = 7, n_patients = nrow(dat), n_patients_hypoxemic = sum(dat$event))
    if (sum(dat$event) < MIN_DEATHS) return(row_head %>% mutate(term = "PFVC", note = paste("skipped: fewer than", MIN_DEATHS, "onsets")))
    fit <- fit_or_separation(fit_strict(coxph(as.formula(paste("Surv(end_day, event) ~", covariate_rhs(adjustment), "+ log_pfvc_z")),
                                              data = dat)))
    if (separated(fit)) return(row_head %>% mutate(term = "PFVC", note = conditionMessage(fit)))
    row_head %>% mutate(term = "PFVC", log_hr = unname(coef(fit)["log_pfvc_z"]),
                        se = unname(sqrt(vcov(fit)["log_pfvc_z", "log_pfvc_z"])), note = NA_character_)
  }
  fit_death_by_state <- function(population, adjustment, horizon) {
    dat <- hypoxemia_population(population) %>%
      mutate(censor_day = pmin(horizon, coalesce(escalation_day, Inf), discharge_index_day),
             death_day = death_index_day,
             event = as.integer(!is.na(death_day) & death_day <= censor_day),
             end_day = pmax(if_else(event == 1L, death_day, censor_day), 0.01),
             hypoxemia_day = if_else(!is.na(hypoxemia_day) & hypoxemia_day < end_day, hypoxemia_day, NA_real_))
    row_head <- tibble(population = POPULATIONS[[population]], adjustment,
                       analysis = "death before escalation, by hypoxemic state", horizon_days = horizon,
                       n_patients = nrow(dat), n_patients_hypoxemic = sum(!is.na(dat$hypoxemia_day)),
                       n_deaths = sum(dat$event),
                       n_deaths_after_hypoxemia = sum(dat$event == 1 & !is.na(dat$hypoxemia_day)))
    if (row_head$n_deaths_after_hypoxemia < MIN_DEATHS || row_head$n_deaths - row_head$n_deaths_after_hypoxemia < MIN_DEATHS)
      return(row_head %>% mutate(term = "PFVC before hypoxemia",
                                 note = paste("skipped: fewer than", MIN_DEATHS, "deaths in a hypoxemic state")))
    # split each patient's follow-up at hypoxemia onset (counting-process form)
    split <- tmerge(dat %>% select(-death_day, -event, -censor_day),
                    dat %>% select(hospitalization_id, follow_up_end = end_day, died = event),
                    id = hospitalization_id, death = event(follow_up_end, died))
    split <- tmerge(split, dat %>% filter(!is.na(hypoxemia_day)) %>% select(hospitalization_id, hypoxemia_day),
                    id = hospitalization_id, hypoxemic = tdc(hypoxemia_day))
    fit <- fit_or_separation(fit_strict(coxph(as.formula(paste("Surv(tstart, tstop, death) ~", covariate_rhs(adjustment),
                                                               "+ hypoxemic + log_pfvc_z + log_pfvc_z:hypoxemic")),
                                              data = split, cluster = hospitalization_id)))
    if (separated(fit)) return(row_head %>% mutate(term = "PFVC before hypoxemia", note = conditionMessage(fit)))
    b <- coef(fit); V <- vcov(fit)
    interaction_term <- names(b)[sapply(strsplit(names(b), ":"), setequal, c("log_pfvc_z", "hypoxemic"))]
    after_b <- b[["log_pfvc_z"]] + b[[interaction_term]]
    after_se <- sqrt(V["log_pfvc_z", "log_pfvc_z"] + V[interaction_term, interaction_term] + 2 * V["log_pfvc_z", interaction_term])
    row_head %>% cross_join(tibble(
      term = c("PFVC before hypoxemia", "PFVC after hypoxemia", "after / before (interaction)", "hypoxemic state"),
      log_hr = c(b[["log_pfvc_z"]], after_b, b[[interaction_term]], b[["hypoxemic"]]),
      se = c(sqrt(V["log_pfvc_z", "log_pfvc_z"]), after_se, sqrt(V[interaction_term, interaction_term]),
             sqrt(V["hypoxemic", "hypoxemic"])))) %>%
      mutate(note = NA_character_)
  }
  hypoxemia_populations <- intersect(c("everyone", "full_code_at_index"), names(POPULATIONS))
  hypoxemia_pathway <- bind_rows(
    expand_grid(population = hypoxemia_populations, adjustment = names(ADJUSTMENTS)) %>% pmap_dfr(fit_onset),
    expand_grid(population = hypoxemia_populations, adjustment = names(ADJUSTMENTS), horizon = HYPOXEMIA_HORIZONS) %>%
      pmap_dfr(fit_death_by_state)) %>%
    mutate(hr = exp(log_hr), hr_lo = exp(log_hr - 1.96 * se), hr_hi = exp(log_hr + 1.96 * se),
           p = 2 * pnorm(-abs(log_hr / se)), cohort = "No support, not hypoxemic on the index day",
           scale = "cause-specific HR per SD of log PFVC (the hypoxemic-state row: HR for being hypoxemic)",
           site = site_name)
} else message("*** No 7-day control panel (jm_long_7d, jm_surv_7d) in ", control_panel_dir,
               ": the hypoxemia pathway is skipped. Build it with PBWPFVC_COHORT=nosupport PBWPFVC_JM_GRID=daily ",
               "PBWPFVC_JM_HORIZON=7 uvr run code/21_biotrauma_panel.R ***")

# =============================================================================
# The channel control contrast
# =============================================================================
# log PFVC split into its GLI pieces (height, age, sex, race; pfvc_channels() in
# 20_biotrauma_grid.R, each in log-PFVC units), each with its own coefficient, in each
# cohort. The pieces replace the demographics, so a demographic direct path (age's
# frailty, race's selection into the ICU) loads onto its piece in both cohorts; the
# difference ventilated minus no support removes it, and strain, which acts only where
# a PBW-scaled volume is delivered, remains. Two readings:
#   each piece's difference   the ventilator-specific association through that input
#   equal differences         strain's prediction: one lung-size effect, so the four
#                             differences agree (Wald, 3 df; an overidentification
#                             test). Also height = sex = race (2 df), leaving out age,
#                             where the two candidate strain denominators part.
# Assumption, as for the whole contrast: each direct path is the same in both cohorts
# (at the within-cohort severity). The pieces are computed once on both cohorts, so
# they share one reference patient. Outcomes: the in-hospital logistic model, 60-day
# death, and 60-day death before intubation; populations as elsewhere. The one-beta
# model (the four pieces' sum, which is log PFVC to a small remainder) is fitted beside
# it for reference. No demographics: the pieces are the demographics in GLI's shape.
#
# Two exposure scales (2026-09-24). log PFVC asks about predicted lung size; log
# PBW/PFVC asks about the strain error, which at a fixed VT/PBW is the delivered
# strain and so the paper's clinical quantity. Devine has no age or race term, so the
# ratio's age and race pieces are GLI's negated exactly, and its sex piece is GLI's
# rescaled (+0.102 against -0.173 at a 170 cm reference): for those three the ratio
# scale is a re-expression, and the conversion table below gives the factor and the R2
# of the proportionality. Height is the exception and is reported, not converted: the
# two formulas' height functions nearly coincide, so the ratio's height piece spans
# about 0.017 log units against GLI's 0.43, and the ratio model's height column
# therefore carries almost no height adjustment. Its coefficient is weakly identified
# and the other three pieces on that scale can absorb height confounding; the
# identifying SD of each piece is reported beside every estimate so this is visible.
source(here::here("code", "20_biotrauma_grid.R"))   # pfvc_channels(), CHANNELS, channels_equal_p()
CHANNEL_EXPOSURES <- c(`log PFVC` = "log_pfvc", `log PBW/PFVC (strain error)` = "ldisc")
channel_frames <- map(CHANNEL_EXPOSURES, ~ bind_cols(both_cohorts, pfvc_channels(both_cohorts, .x)))
# how each piece rescales between the two exposure scales, from the formulas alone: the
# slope through the origin of the ratio piece on the PFVC piece, and the R2 of that
# proportionality. R2 = 1 means the ratio scale is the PFVC scale times a constant
# (beta on the ratio scale = beta on the PFVC scale divided by the slope); below 1
# means the two scales are different functions of the input and no factor exists.
channel_conversion <- map_dfr(CHANNELS, function(piece) {
  x <- channel_frames[["log PFVC"]][[piece]]; y <- channel_frames[["log PBW/PFVC (strain error)"]][[piece]]
  slope <- sum(x * y) / sum(x * x)
  tibble(piece = sub("^ch_", "", piece), sd_pfvc_piece = sd(x), sd_ratio_piece = sd(y),
         ratio_per_pfvc_unit = slope, r2_proportional = 1 - sum((y - slope * x)^2) / sum(y^2),
         convertible = r2_proportional > 0.999, site = site_name)
})
# the SD of each piece left after the other three: what identifies its coefficient
identifying_sd <- function(frame, piece) {
  others <- setdiff(CHANNELS, piece)
  sd(resid(lm(as.formula(paste(piece, "~", paste(others, collapse = " + "))), data = frame)))
}
CHANNEL_OUTCOMES <- OUTCOMES %>% filter(outcome_key %in% c("inhosp_logistic", "day60_all", "day60_before_imv"))
fit_channels <- function(frame, population, cohort_now, outcome_key, model) {
  cohort_data <- frame %>% filter(cohort == cohort_now, .data[[paste0("population_", population)]])
  dat <- if (model == "logistic") cohort_data %>% mutate(event = deceased) else outcome_data(cohort_data, outcome_key)
  if (sum(dat$event == 1) < MIN_DEATHS) return(NULL)
  base_terms <- c(if (cohort_now == "Ventilated") "vtpbw", "sf_z", "sofa_z")
  fit_one <- function(terms) {
    rhs <- paste(c(terms, base_terms), collapse = " + ")
    if (model == "logistic") fit_strict(glm(as.formula(paste("deceased ~", rhs)), family = binomial, data = dat)) else
      fit_strict(coxph(as.formula(paste("Surv(end_day, event) ~", rhs)), data = dat))
  }
  pieces_fit <- fit_or_separation(fit_one(CHANNELS)); one_fit <- fit_or_separation(fit_one("ch_sum"))
  if (separated(pieces_fit) || separated(one_fit)) {
    message("  channel contrast not estimable for ", population, ", ", cohort_now, ", ", outcome_key, ": ",
            conditionMessage(if (separated(pieces_fit)) pieces_fit else one_fit))
    return(NULL)
  }
  list(b = coef(pieces_fit)[CHANNELS], V = vcov(pieces_fit)[CHANNELS, CHANNELS],
       b_one = coef(one_fit)[["ch_sum"]], se_one = sqrt(vcov(one_fit)["ch_sum", "ch_sum"]),
       n_patients = nrow(dat), n_deaths = sum(dat$event == 1))
}
piece_rows <- function(b, V, b_one, se_one, quantity, n_patients, n_deaths) tibble(
  quantity = quantity, piece = c(sub("^ch_", "", CHANNELS), "all four (one beta)"),
  log_ratio = c(unname(b), b_one), se = c(sqrt(diag(V)), se_one), n_patients = n_patients, n_deaths = n_deaths)
channel_rows <- list(); channel_tests <- list(); channel_vcov <- list()
for (exposure in names(CHANNEL_EXPOSURES)) {
frame <- channel_frames[[exposure]]
piece_identifying_sd <- setNames(map_dbl(CHANNELS, ~ identifying_sd(frame, .x)), sub("^ch_", "", CHANNELS))
for (population in names(POPULATIONS)) for (k in seq_len(nrow(CHANNEL_OUTCOMES))) {
  outcome_key <- CHANNEL_OUTCOMES$outcome_key[k]; model <- CHANNEL_OUTCOMES$model[k]
  fits <- map(set_names(COHORTS), ~ fit_channels(frame, population, .x, outcome_key, model))
  if (any(map_lgl(fits, is.null))) next
  ventilated <- fits[["Ventilated"]]; control <- fits[["No support"]]
  difference_b <- ventilated$b - control$b
  difference_V <- ventilated$V + control$V          # independent cohorts
  # exported so the pool can test the pieces' agreement across sites (it needs the
  # covariance between pieces, not only their standard errors)
  channel_vcov[[length(channel_vcov) + 1]] <- as_tibble(as.table(difference_V), .name_repair = "minimal") %>%
    set_names(c("piece_row", "piece_col", "covariance")) %>%
    mutate(across(c(piece_row, piece_col), ~ sub("^ch_", "", .x)),
           exposure = exposure, population = POPULATIONS[[population]], outcome = CHANNEL_OUTCOMES$outcome[k],
           quantity = "ventilated minus no support", .before = 1)
  channel_rows[[length(channel_rows) + 1]] <- bind_rows(
    piece_rows(ventilated$b, ventilated$V, ventilated$b_one, ventilated$se_one, "Ventilated",
               ventilated$n_patients, ventilated$n_deaths),
    piece_rows(control$b, control$V, control$b_one, control$se_one, "No support", control$n_patients, control$n_deaths),
    piece_rows(difference_b, difference_V, ventilated$b_one - control$b_one,
               sqrt(ventilated$se_one^2 + control$se_one^2), "ventilated minus no support", NA_integer_, NA_integer_)) %>%
    mutate(exposure = exposure, population = POPULATIONS[[population]], outcome = CHANNEL_OUTCOMES$outcome[k],
           ratio_type = if_else(model == "logistic", "OR", "HR"),
           identifying_sd = unname(piece_identifying_sd[piece]))
  size_pieces <- c("ch_height", "ch_sex", "ch_race")
  size_contrast <- rbind(c(1, -1, 0), c(1, 0, -1))
  size_d <- size_contrast %*% difference_b[size_pieces]
  size_stat <- as.numeric(t(size_d) %*% solve(size_contrast %*% difference_V[size_pieces, size_pieces] %*% t(size_contrast)) %*% size_d)
  channel_tests[[length(channel_tests) + 1]] <- tibble(
    exposure = exposure, population = POPULATIONS[[population]], outcome = CHANNEL_OUTCOMES$outcome[k],
    test = c("the four differences are equal (3 df)", "height = sex = race differences (2 df)",
             "the four ventilated pieces are equal (3 df, for reference)",
             "the four no-support pieces are equal (3 df, for reference)"),
    p = c(channels_equal_p(difference_b, difference_V), pchisq(size_stat, 2, lower.tail = FALSE),
          channels_equal_p(ventilated$b, ventilated$V), channels_equal_p(control$b, control$V)))
}
}
channel_contrast <- bind_rows(channel_rows) %>%
  mutate(ratio_per_0.1 = exp(0.1 * log_ratio), lo_per_0.1 = exp(0.1 * (log_ratio - 1.96 * se)),
         hi_per_0.1 = exp(0.1 * (log_ratio + 1.96 * se)), p = 2 * pnorm(-abs(log_ratio / se)),
         scale = if_else(exposure == "log PFVC",
                         "per log unit of the piece (log-PFVC units); ratio_per_0.1 is per 0.1 log units (about 10% of PFVC)",
                         "per log unit of the piece (log PBW/PFVC units); ratio_per_0.1 is per 0.1 log units (about 10% more strain than the protocol intends)"),
         site = site_name) %>%
  relocate(exposure, population, outcome, ratio_type, quantity, piece)
channel_tests <- bind_rows(channel_tests) %>% mutate(site = site_name)
channel_vcov <- bind_rows(channel_vcov) %>% mutate(site = site_name)

# how far the control's standardised estimate extrapolates
anchor_overlap <- both_cohorts %>% group_by(cohort) %>%
  summarise(anchor_mean = mean(anchor), anchor_sd = sd(anchor),
            n_patients = n(), n_patients_at_or_above_ventilated_mean = sum(anchor_c >= 0), .groups = "drop") %>%
  mutate(ventilated_mean_anchor = ventilated_mean_anchor, anchor = "SOFA cardiovascular + coagulation + liver + renal",
         site = site_name)

message("\nPFVC per SD of log PFVC, per cohort, and the share of the 40-to-90 age gradient it absorbs:")
print(as.data.frame(estimates %>% select(cohort, adjustment, or_per_sd, or_lo, or_hi, lr_p,
                                         rise_40_to_90_without_pfvc, rise_40_to_90_with_pfvc,
                                         share_of_age_gradient_absorbed, n_patients, n_deaths) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
message("\nSeverity anchor by cohort (the standardised control estimate extrapolates where few controls reach the ventilated mean):")
print(as.data.frame(anchor_overlap %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
message("\nThe contrast by population and outcome, adjusted, severity within cohort (strain predicts the no-support ratio nearer 1;",
        " the severity-standardised rows and the PFVC x anchor terms are in the CSV):")
print(as.data.frame(contrast %>% filter(adjustment == "adjusted", severity == SEVERITY_LABELS[["within"]],
                                        quantity %in% c(COHORTS, "no support minus ventilated")) %>%
                      arrange(factor(population, levels = POPULATIONS), factor(outcome, levels = OUTCOMES$outcome), quantity) %>%
                      select(population, outcome, quantity, ratio, ratio_lo, ratio_hi, p, n_deaths, note) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
message("\nThe channel control contrast (everyone, 60-day death before intubation; per 0.1 log units of each piece;",
        " strain predicts equal ventilated-minus-control differences):")
print(as.data.frame(channel_contrast %>%
                      filter(population == POPULATIONS[["everyone"]], outcome == "60-day death, before invasive ventilation") %>%
                      select(exposure, quantity, piece, ratio_per_0.1, lo_per_0.1, hi_per_0.1, p, identifying_sd) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
print(as.data.frame(channel_tests %>% mutate(p = signif(p, 3))), row.names = FALSE)
write_csv(channel_contrast, file.path(final_dir, paste0("pfvc_age_control_channels_", site_name, ".csv")))
write_csv(channel_tests, file.path(final_dir, paste0("pfvc_age_control_channel_tests_", site_name, ".csv")))
write_csv(channel_vcov, file.path(final_dir, paste0("pfvc_age_control_channel_vcov_", site_name, ".csv")))
message("\nThe two exposure scales: how each piece rescales, from the formulas (R2 = 1: a constant factor exists)")
print(as.data.frame(channel_conversion %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
write_csv(channel_conversion, file.path(final_dir, paste0("pfvc_age_control_channel_conversion_", site_name, ".csv")))
write_csv(code_status_counts,
          file.path(final_dir, paste0("pfvc_age_control_code_status_", site_name, ".csv")))

write_csv(estimates, file.path(final_dir, paste0("pfvc_age_control_estimates_", site_name, ".csv")))
write_csv(curves, file.path(final_dir, paste0("pfvc_age_control_curves_", site_name, ".csv")))
write_csv(contrast, file.path(final_dir, paste0("pfvc_age_control_contrast_", site_name, ".csv")))
write_csv(anchor_overlap, file.path(final_dir, paste0("pfvc_age_control_anchor_", site_name, ".csv")))
message("\nThe control's escalation hazard per SD of log PFVC (below 1: larger predicted lungs escalate less):")
print(as.data.frame(escalation_hazard %>% filter(adjustment == "adjusted") %>%
                      select(population, escalation, severity, term, hr, hr_lo, hr_hi, p, n_patients_escalated, note) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
write_csv(escalation_paths,
          file.path(final_dir, paste0("pfvc_age_control_escalation_paths_", site_name, ".csv")))
write_csv(escalation_hazard,
          file.path(final_dir, paste0("pfvc_age_control_escalation_hazard_", site_name, ".csv")))
if (!is.null(hypoxemia_pathway)) {
  message("\nThe hypoxemia pathway in the control (not hypoxemic on the index day), adjusted: ",
          "onset of hypoxemia by PFVC, and PFVC's death HR before and after hypoxemia develops")
  print(as.data.frame(hypoxemia_pathway %>% filter(adjustment == "adjusted") %>%
                        select(population, analysis, horizon_days, term, hr, hr_lo, hr_hi, p,
                               n_patients, n_patients_hypoxemic, any_of(c("n_deaths", "n_deaths_after_hypoxemia")), note) %>%
                        mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
  write_csv(hypoxemia_pathway,
            file.path(final_dir, paste0("pfvc_age_control_hypoxemia_", site_name, ".csv")))
}

# =============================================================================
# Figure: the age curves by cohort, and the PFVC ratios across the outcome grid
# =============================================================================
curve_panel <- curves %>% filter(adjustment == "adjusted") %>%
  pivot_longer(c(without_pfvc, with_pfvc), names_to = "model", values_to = "log_odds") %>%
  mutate(cohort = factor(cohort, levels = COHORTS)) %>%
  ggplot(aes(age_at_admission, log_odds, colour = model)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_line(linewidth = 0.9) + geom_point(size = 1.5) + facet_wrap(~ cohort, ncol = 1) +
  scale_colour_manual(values = OKABE_ITO, labels = c(without_pfvc = "without log PFVC", with_pfvc = "with log PFVC"),
                      name = NULL) +
  labs(title = "A. Age curve of in-hospital death",
       subtitle = "adjusted; white man at the cohort's median\ndose, severity and PFVC; relative to age 40",
       x = "Age", y = "log-odds of death") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ratio_panel <- contrast %>%
  filter(adjustment == "adjusted", severity == SEVERITY_LABELS[["within"]], !is.na(ratio),
         quantity %in% c(COHORTS, "no support minus ventilated")) %>%
  mutate(quantity = factor(quantity, levels = rev(c("Ventilated", "No support", "no support minus ventilated"))),
         outcome = factor(outcome, levels = OUTCOMES$outcome),
         population = factor(population, levels = POPULATIONS)) %>%
  ggplot(aes(ratio, quantity, colour = population)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = ratio_lo, xmax = ratio_hi), position = position_dodge(width = 0.7)) +
  facet_wrap(~ outcome, ncol = 1) + scale_x_log10() +
  scale_colour_manual(values = c("#009E73", "#CC79A7", "#56B4E9")[seq_along(POPULATIONS)], name = NULL) +
  guides(colour = guide_legend(ncol = 1)) +
  labs(title = "B. Per SD of log PFVC, by outcome and population",
       subtitle = "adjusted, severity within cohort; OR for the logistic row, cause-specific HR otherwise.\nStrain predicts a protective ventilated ratio and a no-support ratio nearer 1",
       x = "ratio per SD of log PFVC, log scale", y = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(final_dir, paste0("pfvc_age_control_", site_name, ".pdf")),
       curve_panel + ratio_panel + plot_layout(widths = c(1, 1.5)), width = 13, height = 14)
message("xsec_pfvc_age_control complete -> ", final_dir)
