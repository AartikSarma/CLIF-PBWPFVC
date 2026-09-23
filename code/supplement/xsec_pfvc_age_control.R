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
#   pfvc_age_control_{site}.pdf
# Usage: Rscript code/supplement/xsec_pfvc_age_control.R   (PBWPFVC_COHORT unset)
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
COHORTS <- c("Ventilated", "No support")
OKABE_ITO <- c(without_pfvc = "#E69F00", with_pfvc = "#0072B2")

# =============================================================================
# Data: both cross-sectional cohorts, on shared columns
# =============================================================================
cohort_columns <- c("hospitalization_id", "recorded_dttm", "age_at_admission", "sex_category",
                    "race_category", "pfvc", "sf_ratio", "sofa_total",
                    "sofa_cv_97", "sofa_coag", "sofa_liver", "sofa_renal",
                    "deceased", "death_dttm", "discharge_dttm")
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
         !is.na(age_at_admission), !is.na(sex_category), !is.na(race_category), !is.na(discharge_dttm)) %>%
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
         anchor = sofa_cv_97 + sofa_coag + sofa_liver + sofa_renal)
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

# a model that warns stops the script (no silent fallback)
fit_strict <- function(expr) withCallingHandlers(expr, warning = function(w)
  stop("model warning, stopping: ", conditionMessage(w), call. = FALSE))

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
  if (model == "logistic") {
    fit <- fit_strict(glm(as.formula(paste("deceased ~", rhs)), family = binomial, data = cohort_data))
    events <- sum(cohort_data$deceased == 1)
  } else {
    dat <- outcome_data(cohort_data, outcome_key)
    fit <- fit_strict(coxph(as.formula(paste("Surv(end_day, event) ~", rhs)), data = dat))
    events <- sum(dat$event)
  }
  b <- coef(fit); V <- vcov(fit)
  interaction_term <- names(b)[sapply(strsplit(names(b), ":"), setequal, c("log_pfvc_z", "anchor_c"))]
  terms <- c(pfvc = "log_pfvc_z", pfvc_x_anchor = if (length(interaction_term)) interaction_term)
  tibble(term = names(terms), log_ratio = unname(b[terms]), se = unname(sqrt(diag(V)[terms])),
         n_patients = nrow(cohort_data), n_deaths = events)
}
per_cohort <- expand_grid(cohort = COHORTS, OUTCOMES, adjustment = names(ADJUSTMENTS), severity = names(severity_rhs)) %>%
  mutate(fit = pmap(list(cohort, outcome_key, model, adjustment, severity),
                    function(cohort_now, outcome_key, model, adjustment, severity)
                      fit_pfvc(filter(both_cohorts, cohort == cohort_now), outcome_key, model, adjustment, severity))) %>%
  unnest(fit)
difference <- per_cohort %>% filter(term == "pfvc") %>%
  select(cohort, outcome_key, adjustment, severity, log_ratio, se) %>%
  pivot_wider(names_from = cohort, values_from = c(log_ratio, se)) %>%
  transmute(outcome_key, adjustment, severity, cohort = "no support minus ventilated", term = "pfvc",
            log_ratio = `log_ratio_No support` - log_ratio_Ventilated,
            se = sqrt(`se_No support`^2 + se_Ventilated^2))
contrast <- bind_rows(per_cohort %>% select(-outcome, -model), difference) %>%
  left_join(OUTCOMES, by = "outcome_key") %>%
  mutate(ratio_type = if_else(model == "logistic", "OR", "HR"),
         quantity = case_when(term == "pfvc_x_anchor" ~ paste0(cohort, ": PFVC x anchor (per SOFA point)"),
                              TRUE ~ cohort),
         ratio = exp(log_ratio), ratio_lo = exp(log_ratio - 1.96 * se), ratio_hi = exp(log_ratio + 1.96 * se),
         p = 2 * pnorm(-abs(log_ratio / se)), severity = SEVERITY_LABELS[severity],
         scale = "per SD of log PFVC (ventilated cohort SD)", site = site_name) %>%
  select(outcome, ratio_type, adjustment, severity, quantity, ratio, ratio_lo, ratio_hi, p, log_ratio, se,
         n_patients, n_deaths, scale, site)

# =============================================================================
# The control's escalation hazard by PFVC
# =============================================================================
# Censoring at escalation is uninformative about PFVC only if PFVC does not predict
# escalation. Cause-specific Cox, control only, for two escalations: any advanced
# support, and invasive ventilation. Death and discharge censor; the horizon is 60
# days. An HR below 1 says larger predicted lungs escalate less, so censoring at
# escalation removes smaller-lung patients preferentially, and the before-escalation
# death estimate is read on a population that PFVC itself selected.
fit_escalation <- function(escalation_column, escalation_label, adjustment, severity) {
  dat <- both_cohorts %>% filter(cohort == "No support") %>%
    mutate(escalation_time = .data[[escalation_column]],
           death_or_discharge = pmin(coalesce(death_index_day, Inf), discharge_index_day),
           censor_day = pmin(HORIZON_DAYS, death_or_discharge),
           event = as.integer(!is.na(escalation_time) & escalation_time <= censor_day),
           end_day = pmax(if_else(event == 1L, escalation_time, censor_day), 0.01))
  rhs <- paste(c(severity_rhs[[severity]], "ns(age_at_admission, 4)",
                 if (!is.na(ADJUSTMENTS[[adjustment]])) ADJUSTMENTS[[adjustment]], "log_pfvc_z"), collapse = " + ")
  fit <- fit_strict(coxph(as.formula(paste("Surv(end_day, event) ~", rhs)), data = dat))
  b <- coef(fit); se <- sqrt(diag(vcov(fit)))
  interaction_term <- names(b)[sapply(strsplit(names(b), ":"), setequal, c("log_pfvc_z", "anchor_c"))]
  terms <- c(pfvc = "log_pfvc_z", pfvc_x_anchor = if (length(interaction_term)) interaction_term)
  tibble(escalation = escalation_label, adjustment = adjustment, severity = SEVERITY_LABELS[[severity]],
         term = names(terms), log_hr = unname(b[terms]), se = unname(se[terms]),
         n_patients = nrow(dat), n_patients_escalated = sum(dat$event))
}
escalation_hazard <- expand_grid(
  tibble(escalation_column = c("escalation_day", "imv_day"),
         escalation_label = c("any advanced support", "invasive ventilation")),
  adjustment = names(ADJUSTMENTS), severity = names(severity_rhs)) %>%
  pmap_dfr(fit_escalation) %>%
  mutate(hr = exp(log_hr), hr_lo = exp(log_hr - 1.96 * se), hr_hi = exp(log_hr + 1.96 * se),
         p = 2 * pnorm(-abs(log_hr / se)), cohort = "No support",
         scale = "cause-specific HR of escalation per SD of log PFVC", site = site_name)

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
message("\nThe contrast by outcome and severity form, adjusted (strain predicts the no-support ratio nearer 1):")
print(as.data.frame(contrast %>% filter(adjustment == "adjusted") %>%
                      arrange(factor(outcome, levels = OUTCOMES$outcome), severity, quantity) %>%
                      select(outcome, severity, quantity, ratio, ratio_lo, ratio_hi, p, n_deaths) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)

write_csv(mask_small_counts(estimates), file.path(final_dir, paste0("pfvc_age_control_estimates_", site_name, ".csv")))
write_csv(curves, file.path(final_dir, paste0("pfvc_age_control_curves_", site_name, ".csv")))
write_csv(mask_small_counts(contrast), file.path(final_dir, paste0("pfvc_age_control_contrast_", site_name, ".csv")))
write_csv(mask_small_counts(anchor_overlap), file.path(final_dir, paste0("pfvc_age_control_anchor_", site_name, ".csv")))
message("\nThe control's escalation hazard per SD of log PFVC (below 1: larger predicted lungs escalate less):")
print(as.data.frame(escalation_hazard %>% filter(adjustment == "adjusted") %>%
                      select(escalation, severity, term, hr, hr_lo, hr_hi, p, n_patients_escalated) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
write_csv(mask_small_counts(escalation_paths),
          file.path(final_dir, paste0("pfvc_age_control_escalation_paths_", site_name, ".csv")))
write_csv(mask_small_counts(escalation_hazard),
          file.path(final_dir, paste0("pfvc_age_control_escalation_hazard_", site_name, ".csv")))

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
  filter(adjustment == "adjusted", quantity %in% c(COHORTS, "no support minus ventilated")) %>%
  mutate(quantity = factor(quantity, levels = rev(c("Ventilated", "No support", "no support minus ventilated"))),
         outcome = factor(outcome, levels = OUTCOMES$outcome),
         severity = factor(severity, levels = SEVERITY_LABELS)) %>%
  ggplot(aes(ratio, quantity, colour = severity)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = ratio_lo, xmax = ratio_hi), position = position_dodge(width = 0.6)) +
  facet_wrap(~ outcome, ncol = 1) + scale_x_log10() +
  scale_colour_manual(values = c("#009E73", "#CC79A7"), name = NULL) +
  labs(title = "B. Per SD of log PFVC, by outcome and severity form",
       subtitle = "adjusted; OR for the logistic row, cause-specific HR otherwise. Strain predicts\na protective ventilated ratio and a no-support ratio nearer 1",
       x = "ratio per SD of log PFVC, log scale", y = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(final_dir, paste0("pfvc_age_control_", site_name, ".pdf")),
       curve_panel + ratio_panel + plot_layout(widths = c(1, 1.5)), width = 13, height = 14)
message("xsec_pfvc_age_control complete -> ", final_dir)
