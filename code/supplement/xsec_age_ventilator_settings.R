# =============================================================================
# Supplement (cross-sectional, exploratory): do the mortality gradients of GLI's age
# decline and of predicted lung size change with the ventilator's rate and PEEP?
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
# The same interactions for predicted lung size and the strain error: log PFVC, log
# PFVC at age 25 (height, sex and race only) and log PBW/PFVC. Under the usable-range
# model (a breath must fit between end-expiratory volume and total lung capacity),
# PEEP raises end-expiratory volume and shrinks the room left, most in the smallest
# lungs, so a larger PFVC should protect MORE at higher PEEP (log PFVC x PEEP below 0;
# log PBW/PFVC x PEEP above 0). This prediction was written after the age x PEEP
# result (steeper at higher PEEP, against the closure prediction above), so it is
# post hoc.
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
#   age_setting_interaction_{site}.csv  each exposure's gradient, the setting's main term
#                                       and their interaction, per outcome, setting and
#                                       adjustment (column exposure)
#   age_setting_tertiles_{site}.csv     each exposure's gradient within groups of the setting
#   age_setting_settings_{site}.csv     the settings' distribution and missingness
#   age_setting_peep_volume_{site}.csv  PEEP as volume (PEEP / specific elastance) in the
#                                       plateau-measured patients: the interactions, and
#                                       dynamic against total strain (section at the end)
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
                                      "pfvc", "pfvc_age25", "pbw", "vtpbw", "sf_ratio", "sofa_total", "resp_rate_set", "peep_set",
                                      "deceased", "surv_time", "mortality_event_60", "tidal_volume_set", "dp", "ers_pfvc"))
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
  filter(!is.na(pfvc), pfvc > 0, !is.na(pfvc_age25), pfvc_age25 > 0, !is.na(pbw), pbw > 0, !is.na(vtpbw), !is.na(sf_ratio),
         !is.na(sofa_total), !is.na(resp_rate_set), resp_rate_set > 0, !is.na(peep_set), peep_set >= 0,
         !is.na(deceased), !is.na(surv_time), !is.na(mortality_event_60),
         !is.na(sex_category), !is.na(race_category)) %>%
  mutate(sex_category  = factor(sex_category, levels = c("Male", "Female")),
         race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
         # GLI's age decline in log units: what the patient's age has taken off the
         # vital capacity predicted for their height, sex and race at age 25
         age_decline = log(pfvc_age25) - log(pfvc),
         log_pfvc = log(pfvc), log_pfvc25 = log(pfvc_age25), log_ratio = log(pbw / pfvc),
         age10 = age_at_admission / 10,
         sf_z = as.numeric(scale(sf_ratio)), sofa_z = as.numeric(scale(sofa_total)),
         rate_z = as.numeric(scale(resp_rate_set)), peep_z = as.numeric(scale(peep_set)),
         surv_time = pmax(surv_time, 0.01)) %>%
  # every exposure centred at the cohort mean, so a setting's main term is its
  # association at the average patient, not at an exposure of zero
  mutate(across(c(age_decline, age10, log_pfvc, log_pfvc25, log_ratio), ~ .x - mean(.x)))
message("Excluded for a missing input, setting or outcome: ", n_before - nrow(cohort), " of ", n_before)
if (sum(cohort$deceased == 1) < MIN_EVENTS) stop("fewer than ", MIN_EVENTS, " in-hospital deaths")

fit_strict <- function(expr) withCallingHandlers(expr, warning = function(w)
  stop("model warning, stopping: ", conditionMessage(w), call. = FALSE))
OUTCOMES <- c(inhosp_logistic = "in-hospital death (logistic)", day60_cox = "60-day death (Cox)")
EXPOSURES <- c(age_decline = "GLI age decline (log PFVC at 25 minus log PFVC)", age10 = "age, per decade",
               log_pfvc = "log PFVC", log_pfvc25 = "log PFVC at age 25", log_ratio = "log PBW/PFVC")
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

interaction_rows <- expand_grid(outcome_key = names(OUTCOMES), exposure = names(EXPOSURES),
                                setting_set = names(SETTINGS), adjustment = names(ADJUSTMENTS)) %>%
  pmap_dfr(function(outcome_key, exposure, setting_set, adjustment) {
    settings <- SETTINGS[[setting_set]]
    rhs <- paste(c("sf_z", "sofa_z", "vtpbw", settings, exposure, paste0(exposure, ":", settings),
                   if (!is.na(ADJUSTMENTS[[adjustment]])) ADJUSTMENTS[[adjustment]]), collapse = " + ")
    fit <- fit_outcome(outcome_key, rhs)
    b <- coef(fit); V <- vcov(fit)
    terms <- c(setNames(exposure, "exposure gradient (at the mean setting)"),
               setNames(settings, paste("main term:", SETTING_LABELS[settings])),
               setNames(map_chr(settings, ~ coefficient_named(b, c(exposure, .x))),
                        paste("interaction: exposure x", SETTING_LABELS[settings])))
    tibble(outcome = OUTCOMES[[outcome_key]], exposure = EXPOSURES[[exposure]],
           settings_in_model = paste(SETTING_LABELS[settings], collapse = " and "), adjustment = adjustment,
           quantity = names(terms), estimate = unname(b[terms]), se = unname(sqrt(diag(V)[terms])),
           n_patients = nrow(cohort),
           n_deaths = if (outcome_key == "inhosp_logistic") sum(cohort$deceased == 1) else sum(cohort$mortality_event_60 == 1))
  }) %>%
  mutate(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se, p = 2 * pnorm(-abs(estimate / se)),
         scale = case_when(grepl("^main term", quantity) ~ "log ratio per SD of the setting",
                           grepl("^exposure gradient", quantity) & grepl("decade", exposure) ~ "log ratio per decade",
                           grepl("^exposure gradient", quantity) ~ "log ratio per log unit of the exposure (x 0.1 for about 10%)",
                           grepl("decade", exposure) ~ "change in the per-decade gradient per SD of the setting",
                           TRUE ~ "change in the gradient per log unit per SD of the setting"),
         site = site_name)

# the age gradient within each tertile of each setting, for reading the interaction
tertile_rows <- expand_grid(outcome_key = names(OUTCOMES), exposure = names(EXPOSURES),
                            setting = c("resp_rate_set", "peep_set"), adjustment = names(ADJUSTMENTS)) %>%
  pmap_dfr(function(outcome_key, exposure, setting, adjustment) {
    # PEEP clusters at 5 cmH2O, so equal-count tertiles would split one value across
    # groups: PEEP is grouped at clinical cut points (5 or less, 6 to 9, 10 or more),
    # and the rate at its tertiles with tied values kept together
    tertile_data <- cohort %>% mutate(setting_tertile = if (setting == "peep_set")
      cut(peep_set, breaks = c(-Inf, 5, 9, Inf), labels = c("low", "middle", "high")) else
      cut(resp_rate_set, breaks = unique(quantile(resp_rate_set, c(0, 1 / 3, 2 / 3, 1))), include.lowest = TRUE,
          labels = c("low", "middle", "high")))
    rhs <- paste(c("sf_z", "sofa_z", "vtpbw", "setting_tertile", paste0(exposure, ":setting_tertile"),
                   if (!is.na(ADJUSTMENTS[[adjustment]])) ADJUSTMENTS[[adjustment]]), collapse = " + ")
    fit <- fit_outcome(outcome_key, rhs, tertile_data)
    b <- coef(fit); V <- vcov(fit)
    map_dfr(c("low", "middle", "high"), function(tertile) {
      term <- coefficient_named(b, c(exposure, paste0("setting_tertile", tertile)))
      # computed before tibble(), whose own `setting` column would mask the argument
      tertile_range <- paste(range(tertile_data[[setting]][tertile_data$setting_tertile == tertile]), collapse = " to ")
      tibble(outcome = OUTCOMES[[outcome_key]], exposure = EXPOSURES[[exposure]],
             setting = if (setting == "resp_rate_set") "set respiratory rate" else "PEEP", adjustment = adjustment,
             tertile = tertile,
             setting_range = tertile_range,
             n_patients = sum(tertile_data$setting_tertile == tertile),
             estimate = unname(b[term]), se = unname(sqrt(V[term, term])))
    })
  }) %>%
  mutate(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se, site = site_name)

message("\nExposure x setting, in-hospital death, adjusted, each setting alone:")
print(as.data.frame(interaction_rows %>%
                      filter(grepl("^in-hospital", outcome), adjustment == "adjusted", settings_in_model != "set respiratory rate and PEEP") %>%
                      select(exposure, settings_in_model, quantity, estimate, lo, hi, p) %>%
                      mutate(exposure = substr(exposure, 1, 16), across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)
message("\nThe age-decline gradient by tertile of the setting (in-hospital, adjusted; log odds per log unit):")
print(as.data.frame(tertile_rows %>% filter(grepl("^in-hospital", outcome), adjustment == "adjusted", grepl("decline", exposure)) %>%
                      select(setting, tertile, setting_range, n_patients, estimate, lo, hi) %>%
                      mutate(across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)

# =============================================================================
# PEEP as volume, in the plateau-measured patients
# =============================================================================
# The volume PEEP adds is PEEP x Crs, and as a fraction of predicted lung size it is
# PEEP / (Ers x PFVC): PEEP over specific elastance. Specific elastance falls steeply
# with age (xsec_crs_channels.R), so the same PEEP in cmH2O raises an older patient's
# end-expiratory volume by a larger fraction, and leaves less of the usable range for
# the breath. If that is the mechanism, the age x PEEP interaction above should be
# carried by the PEEP volume fraction rather than by PEEP in cmH2O. In the patients
# with a plateau pressure (driving pressure at least DP_FLOOR cmH2O, as compliance
# from a smaller driving pressure is not physiologic), per exposure:
#   PEEP in cmH2O          the interaction above, in this subset
#   PEEP volume fraction   PEEP / specific elastance, the fraction of PFVC PEEP adds
#   both                   which of the two carries the interaction
# and, for the whole breath, a head-to-head of dynamic strain (VT / PFVC) against total
# strain ((VT + PEEP x Crs) / PFVC), each on the log scale, by AIC (below 0 favours
# total strain).
# Caveats: Crs is measured with the disease, so it carries severity; the volume assumes
# a linear pressure-volume relation at the set PEEP (no recruitment); and PFVC is in
# both the volume fraction and the size exposures, so their interaction is partly
# arithmetic (the ratio problem). The subset is selected on a recorded plateau.
DP_FLOOR <- 5
peep_data <- cohort %>%
  filter(!is.na(dp), dp >= DP_FLOOR, !is.na(ers_pfvc), ers_pfvc > 0, !is.na(tidal_volume_set), tidal_volume_set > 0) %>%
  mutate(peep_fraction = peep_set / ers_pfvc,                            # fraction of PFVC added by PEEP
         dynamic_strain = tidal_volume_set / (1000 * pfvc),
         log_dynamic_strain = log(dynamic_strain),
         log_total_strain = log(dynamic_strain + peep_fraction),
         peep_z = as.numeric(scale(peep_set)), peep_fraction_z = as.numeric(scale(peep_fraction)),
         sf_z = as.numeric(scale(sf_ratio)), sofa_z = as.numeric(scale(sofa_total)))
message("\nPEEP as volume: ", nrow(peep_data), " of ", nrow(cohort), " patients with a plateau (driving pressure >= ", DP_FLOOR,
        "); PEEP volume fraction median ", signif(median(peep_data$peep_fraction), 3), " of PFVC (IQR ",
        paste(signif(quantile(peep_data$peep_fraction, c(0.25, 0.75)), 3), collapse = " to "), ")")
PEEP_FORMS <- list(`PEEP (cmH2O)` = "peep_z", `PEEP volume fraction` = "peep_fraction_z", both = c("peep_z", "peep_fraction_z"))
PEEP_LABELS <- c(peep_z = "PEEP (cmH2O)", peep_fraction_z = "PEEP volume fraction")
peep_volume_rows <- if (sum(peep_data$deceased == 1) < MIN_EVENTS) {
  message("  fewer than ", MIN_EVENTS, " deaths among the plateau-measured patients: PEEP as volume skipped")
  tibble()
} else bind_rows(
  expand_grid(outcome_key = names(OUTCOMES), exposure = c("age_decline", "log_pfvc", "log_pfvc25"),
              peep_form = names(PEEP_FORMS), adjustment = names(ADJUSTMENTS)) %>%
    pmap_dfr(function(outcome_key, exposure, peep_form, adjustment) {
      settings <- PEEP_FORMS[[peep_form]]
      # fitted on the exposure in SD units (a log exposure's small SD beside two
      # correlated PEEP terms stalls the Cox fit), reported per log unit
      exposure_sd <- sd(peep_data[[exposure]])
      fit_data <- peep_data %>% mutate(exposure_scaled = .data[[exposure]] / exposure_sd)
      rhs <- paste(c("sf_z", "sofa_z", "vtpbw", settings, "exposure_scaled", paste0("exposure_scaled:", settings),
                     if (!is.na(ADJUSTMENTS[[adjustment]])) ADJUSTMENTS[[adjustment]]), collapse = " + ")
      fit <- fit_outcome(outcome_key, rhs, fit_data)
      b <- coef(fit); V <- vcov(fit)
      terms <- setNames(map_chr(settings, ~ coefficient_named(b, c("exposure_scaled", .x))),
                        paste("interaction: exposure x", PEEP_LABELS[settings]))
      tibble(analysis = "interaction", outcome = OUTCOMES[[outcome_key]], exposure = EXPOSURES[[exposure]],
             peep_in_model = peep_form, adjustment = adjustment, quantity = names(terms),
             estimate = unname(b[terms]) / exposure_sd, se = unname(sqrt(diag(V)[terms])) / exposure_sd)
    }),
  expand_grid(outcome_key = names(OUTCOMES), adjustment = names(ADJUSTMENTS)) %>%
    pmap_dfr(function(outcome_key, adjustment) {
      fit_for <- function(strain) fit_outcome(outcome_key, paste(c("sf_z", "sofa_z", strain,
        if (!is.na(ADJUSTMENTS[[adjustment]])) ADJUSTMENTS[[adjustment]]), collapse = " + "), peep_data)
      dynamic_fit <- fit_for("log_dynamic_strain"); total_fit <- fit_for("log_total_strain")
      tibble(analysis = "strain head-to-head", outcome = OUTCOMES[[outcome_key]], adjustment = adjustment,
             quantity = c("log dynamic strain (VT / PFVC)", "log total strain ((VT + PEEP x Crs) / PFVC)",
                          "AIC, total minus dynamic (below 0 favours total strain)"),
             estimate = c(coef(dynamic_fit)[["log_dynamic_strain"]], coef(total_fit)[["log_total_strain"]],
                          AIC(total_fit) - AIC(dynamic_fit)),
             se = c(sqrt(vcov(dynamic_fit)["log_dynamic_strain", "log_dynamic_strain"]),
                    sqrt(vcov(total_fit)["log_total_strain", "log_total_strain"]), NA_real_))
    })) %>%
  mutate(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se, p = 2 * pnorm(-abs(estimate / se)),
         n_patients = nrow(peep_data), n_deaths_in_hospital = sum(peep_data$deceased == 1),
         peep_fraction_median = median(peep_data$peep_fraction), dp_floor = DP_FLOOR, site = site_name)
if (nrow(peep_volume_rows)) {
  message("PEEP as volume, in-hospital, adjusted (interactions per SD of each PEEP form; strain per log unit):")
  print(as.data.frame(peep_volume_rows %>% filter(grepl("^in-hospital", outcome), adjustment == "adjusted") %>%
                        select(analysis, exposure, peep_in_model, quantity, estimate, lo, hi, p) %>%
                        mutate(exposure = substr(exposure, 1, 16), across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)
}

write_csv(interaction_rows, file.path(final_dir, paste0("age_setting_interaction_", site_name, ".csv")))
write_csv(tertile_rows, file.path(final_dir, paste0("age_setting_tertiles_", site_name, ".csv")))
write_csv(settings_summary, file.path(final_dir, paste0("age_setting_settings_", site_name, ".csv")))
if (nrow(peep_volume_rows))
  write_csv(peep_volume_rows, file.path(final_dir, paste0("age_setting_peep_volume_", site_name, ".csv")))
message("\nWrote age_setting_{interaction,tertiles,settings,peep_volume}_", site_name, ".csv to ", final_dir)
