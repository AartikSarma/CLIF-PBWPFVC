# =============================================================================
# Supplement (cross-sectional, exploratory): does the mortality gradient of GLI's
# age decline change with the ventilator's rate and PEEP?
# =============================================================================
# The strain denominator is plausibly the lung volume that can take part in tidal
# ventilation, which shrinks with age as GLI's FVC does: closing volume rises, units
# close above end-expiratory volume, and slow units cannot fill at ventilator rates.
# If so, GLI's age decline is part of the strain error that PBW (which has no age
# term) delivers, and two settings should move its mortality gradient:
#   set respiratory rate   more units drop out as the rate rises (time constants):
#                          the age gradient should be STEEPER at higher rates
#   PEEP                   PEEP holds closed units open: the age gradient should be
#                          FLATTER at higher PEEP
# These are exploratory. The clinician chooses both settings for acidosis and
# severity, and acidosis and severity can interact with frailty, so an age x setting
# interaction is a pattern consistent with the mechanism, not a test that separates
# strain from frailty.
#
# The age exposure is GLI's age decline: log PFVC at age 25 minus log PFVC (script
# 03's pfvc_age25 and pfvc), the fraction of predicted vital capacity lost to age for
# the patient's height, sex and race, per 0.1 log units. Linear age per decade is
# reported beside it.
#
# Model, ventilated cohort (script 03's cross-sectional table, every patient):
#   death ~ SF (z) + SOFA (z) + VT/PBW + setting (z) + age decline + age decline x setting
#           [+ sex + race, the adjusted fit]
# for each setting alone and for both together; in-hospital death (logistic) and
# 60-day death (Cox, script 03's surv_time and mortality_event_60, from the index).
# The setting is the index row's set rate or PEEP, standardised over the cohort. The
# age gradient within three groups of the setting (the rate's tertiles with ties kept
# together; PEEP at 5 or less, 6 to 9 and 10 or more cmH2O) is written beside the
# interaction, for reading. The unadjusted fit keeps severity and the dose and drops
# sex and race; age cannot be adjusted away because it is the exposure.
#
# A model that warns stops the script; no model is replaced by a simpler one.
#
# Inputs : intermediate/analysis_cross_sectional.parquet (script 03, ventilated)
# Outputs: final/supplement/
#   age_setting_interaction_{site}.csv  the age-decline (or age) gradient, the setting's
#                                       main term and their interaction, per outcome,
#                                       setting and adjustment
#   age_setting_tertiles_{site}.csv     the age gradient within each tertile of the setting
#   age_setting_settings_{site}.csv     the settings' distribution and missingness
# Usage: uvr run code/supplement/xsec_age_ventilator_settings.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(survival)
})

options(width = 220)
source("utils/config.R")
if (config$cohort != "imv") stop("xsec_age_ventilator_settings.R reads the ventilated cohort: unset PBWPFVC_COHORT")
site_name <- config$site_name
final_dir <- final_dir_for("supplement")

MIN_EVENTS <- 10L   # minimum deaths per model: the CLIF minimum-count standard

cohort <- read_parquet(file.path(config$output_dir, "analysis_cross_sectional.parquet"),
                       col_select = c("hospitalization_id", "age_at_admission", "sex_category", "race_category",
                                      "pfvc", "pfvc_age25", "vtpbw", "sf_ratio", "sofa_total", "resp_rate_set", "peep_set",
                                      "deceased", "surv_time", "mortality_event_60"))
# SYNTHETIC SITE ONLY: synthetic CLIF mortality is unreliable, so death is simulated
# independently of every exposure (35% by day 60, time log-normal with median 9 days;
# a simulated death is in hospital). The run exercises the machinery and can show no
# real effect. Never runs at a real site.
if (grepl("^synthetic_clif", site_name)) {
  message("*** SYNTHETIC SITE: simulated mortality (plumbing only; synthetic CLIF mortality is unreliable). ***")
  set.seed(20260930)
  simulated_death <- rbinom(nrow(cohort), 1L, 0.35)
  simulated_day   <- pmin(pmax(rlnorm(nrow(cohort), log(9), 0.95), 0.04), 60)
  cohort <- cohort %>% mutate(deceased = simulated_death, mortality_event_60 = simulated_death,
                              surv_time = if_else(simulated_death == 1L, simulated_day, 60))
}

settings_summary <- cohort %>%
  summarise(n_patients = n(),
            across(c(resp_rate_set, peep_set),
                   list(n_missing = ~ sum(is.na(.x)), median = ~ median(.x, na.rm = TRUE),
                        p25 = ~ quantile(.x, 0.25, na.rm = TRUE), p75 = ~ quantile(.x, 0.75, na.rm = TRUE)))) %>%
  mutate(site = site_name)
message("Index settings (set rate, breaths/min; PEEP, cmH2O):")
print(as.data.frame(settings_summary), row.names = FALSE)

# Patients missing an input, a setting or an outcome are excluded; the count is printed.
n_before <- nrow(cohort)
cohort <- cohort %>%
  filter(!is.na(pfvc), pfvc > 0, !is.na(pfvc_age25), pfvc_age25 > 0, !is.na(vtpbw), !is.na(sf_ratio),
         !is.na(sofa_total), !is.na(resp_rate_set), resp_rate_set > 0, !is.na(peep_set), peep_set >= 0,
         !is.na(deceased), !is.na(surv_time), !is.na(mortality_event_60),
         !is.na(sex_category), !is.na(race_category)) %>%
  mutate(sex_category  = factor(sex_category, levels = c("Male", "Female")),
         race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
         # GLI's age decline in log units: what the patient's age has taken off the
         # vital capacity predicted for their height, sex and race at age 25
         age_decline = log(pfvc_age25) - log(pfvc),
         age10 = age_at_admission / 10,
         sf_z = as.numeric(scale(sf_ratio)), sofa_z = as.numeric(scale(sofa_total)),
         rate_z = as.numeric(scale(resp_rate_set)), peep_z = as.numeric(scale(peep_set)),
         surv_time = pmax(surv_time, 0.01))
message("Excluded for a missing input, setting or outcome: ", n_before - nrow(cohort), " of ", n_before)
if (sum(cohort$deceased == 1) < MIN_EVENTS) stop("fewer than ", MIN_EVENTS, " in-hospital deaths")

fit_strict <- function(expr) withCallingHandlers(expr, warning = function(w)
  stop("model warning, stopping: ", conditionMessage(w), call. = FALSE))
OUTCOMES <- c(inhosp_logistic = "in-hospital death (logistic)", day60_cox = "60-day death (Cox)")
AGE_TERMS <- c(age_decline = "GLI age decline (log PFVC at 25 minus log PFVC)", age10 = "age, per decade")
SETTINGS <- list(rate = "rate_z", peep = "peep_z", both = c("rate_z", "peep_z"))
SETTING_LABELS <- c(rate_z = "set respiratory rate", peep_z = "PEEP")
ADJUSTMENTS <- c(adjusted = "sex_category + race_category", unadjusted = NA_character_)
fit_outcome <- function(outcome_key, rhs, data = cohort) {
  if (outcome_key == "inhosp_logistic")
    fit_strict(glm(as.formula(paste("deceased ~", rhs)), family = binomial, data = data)) else
    fit_strict(coxph(as.formula(paste("Surv(surv_time, mortality_event_60) ~", rhs)), data = data))
}
# the coefficient of a term named by its components, in either order
coefficient_named <- function(b, parts) names(b)[vapply(strsplit(names(b), ":"), setequal, logical(1), parts)]

interaction_rows <- expand_grid(outcome_key = names(OUTCOMES), age_term = names(AGE_TERMS),
                                setting_set = names(SETTINGS), adjustment = names(ADJUSTMENTS)) %>%
  pmap_dfr(function(outcome_key, age_term, setting_set, adjustment) {
    settings <- SETTINGS[[setting_set]]
    rhs <- paste(c("sf_z", "sofa_z", "vtpbw", settings, age_term, paste0(age_term, ":", settings),
                   if (!is.na(ADJUSTMENTS[[adjustment]])) ADJUSTMENTS[[adjustment]]), collapse = " + ")
    fit <- fit_outcome(outcome_key, rhs)
    b <- coef(fit); V <- vcov(fit)
    terms <- c(setNames(age_term, "age gradient (at the mean setting)"),
               setNames(settings, paste("main term:", SETTING_LABELS[settings])),
               setNames(map_chr(settings, ~ coefficient_named(b, c(age_term, .x))),
                        paste("interaction: age x", SETTING_LABELS[settings])))
    tibble(outcome = OUTCOMES[[outcome_key]], age_term = AGE_TERMS[[age_term]],
           settings_in_model = paste(SETTING_LABELS[settings], collapse = " and "), adjustment = adjustment,
           quantity = names(terms), estimate = unname(b[terms]), se = unname(sqrt(diag(V)[terms])),
           n_patients = nrow(cohort),
           n_deaths = if (outcome_key == "inhosp_logistic") sum(cohort$deceased == 1) else sum(cohort$mortality_event_60 == 1))
  }) %>%
  mutate(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se, p = 2 * pnorm(-abs(estimate / se)),
         scale = case_when(grepl("^age gradient", quantity) & grepl("decline", age_term) ~ "log ratio per log unit of age decline (x 0.1 for 10% of vital capacity)",
                           grepl("^age gradient", quantity) ~ "log ratio per decade",
                           grepl("^main term", quantity) ~ "log ratio per SD of the setting",
                           grepl("decline", age_term) ~ "change in the age-decline gradient (per log unit) per SD of the setting",
                           TRUE ~ "change in the per-decade gradient per SD of the setting"),
         site = site_name)

# the age gradient within each tertile of each setting, for reading the interaction
tertile_rows <- expand_grid(outcome_key = names(OUTCOMES), age_term = names(AGE_TERMS),
                            setting = c("resp_rate_set", "peep_set"), adjustment = names(ADJUSTMENTS)) %>%
  pmap_dfr(function(outcome_key, age_term, setting, adjustment) {
    # PEEP clusters at 5 cmH2O, so equal-count tertiles would split one value across
    # groups: PEEP is grouped at clinical cut points (5 or less, 6 to 9, 10 or more),
    # and the rate at its tertiles with tied values kept together
    tertile_data <- cohort %>% mutate(setting_tertile = if (setting == "peep_set")
      cut(peep_set, breaks = c(-Inf, 5, 9, Inf), labels = c("low", "middle", "high")) else
      cut(resp_rate_set, breaks = unique(quantile(resp_rate_set, c(0, 1 / 3, 2 / 3, 1))), include.lowest = TRUE,
          labels = c("low", "middle", "high")))
    rhs <- paste(c("sf_z", "sofa_z", "vtpbw", "setting_tertile", paste0(age_term, ":setting_tertile"),
                   if (!is.na(ADJUSTMENTS[[adjustment]])) ADJUSTMENTS[[adjustment]]), collapse = " + ")
    fit <- fit_outcome(outcome_key, rhs, tertile_data)
    b <- coef(fit); V <- vcov(fit)
    map_dfr(c("low", "middle", "high"), function(tertile) {
      term <- coefficient_named(b, c(age_term, paste0("setting_tertile", tertile)))
      # computed before tibble(), whose own `setting` column would mask the argument
      tertile_range <- paste(range(tertile_data[[setting]][tertile_data$setting_tertile == tertile]), collapse = " to ")
      tibble(outcome = OUTCOMES[[outcome_key]], age_term = AGE_TERMS[[age_term]],
             setting = if (setting == "resp_rate_set") "set respiratory rate" else "PEEP", adjustment = adjustment,
             tertile = tertile,
             setting_range = tertile_range,
             n_patients = sum(tertile_data$setting_tertile == tertile),
             estimate = unname(b[term]), se = unname(sqrt(V[term, term])))
    })
  }) %>%
  mutate(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se, site = site_name)

message("\nAge x setting, in-hospital death, adjusted, each setting alone (strain predicts rate +, PEEP -):")
print(as.data.frame(interaction_rows %>%
                      filter(grepl("^in-hospital", outcome), adjustment == "adjusted", settings_in_model != "set respiratory rate and PEEP") %>%
                      select(age_term, settings_in_model, quantity, estimate, lo, hi, p) %>%
                      mutate(age_term = substr(age_term, 1, 16), across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)
message("\nThe age-decline gradient by tertile of the setting (in-hospital, adjusted; log odds per log unit):")
print(as.data.frame(tertile_rows %>% filter(grepl("^in-hospital", outcome), adjustment == "adjusted", grepl("decline", age_term)) %>%
                      select(setting, tertile, setting_range, n_patients, estimate, lo, hi) %>%
                      mutate(across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)

write_csv(interaction_rows, file.path(final_dir, paste0("age_setting_interaction_", site_name, ".csv")))
write_csv(tertile_rows, file.path(final_dir, paste0("age_setting_tertiles_", site_name, ".csv")))
write_csv(settings_summary, file.path(final_dir, paste0("age_setting_settings_", site_name, ".csv")))
message("\nWrote age_setting_{interaction,tertiles,settings}_", site_name, ".csv to ", final_dir)
