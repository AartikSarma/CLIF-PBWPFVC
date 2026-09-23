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
# Sensitivity, 60-day death before escalation (cause-specific Cox). Script 03
# removes control patients escalated to any advanced support within 24 hours of
# the index; those escalated later stay in, and their in-hospital deaths can follow
# ventilation, which lets strain into the control. Here death counts only before
# escalation, and escalation censors.
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
#   pfvc_age_control_contrast_{site}.csv   cohort x log PFVC: logistic and
#                                          cause-specific Cox
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
                    "race_category", "pfvc", "sf_ratio", "sofa_total", "deceased",
                    "mortality_event_60", "surv_time")
control_file <- file.path(config$output_dir, "controls", "nosupport", "analysis_cross_sectional.parquet")
if (!file.exists(control_file))
  stop("no no-support cohort: run scripts 01-03 with PBWPFVC_COHORT=nosupport first")
ventilated <- read_parquet(file.path(config$output_dir, "analysis_cross_sectional.parquet")) %>%
  select(all_of(cohort_columns), vtpbw) %>%
  mutate(cohort = "Ventilated", escalation_dttm = as.POSIXct(NA))
no_support <- read_parquet(control_file) %>%
  select(all_of(cohort_columns), escalation_dttm) %>%
  mutate(cohort = "No support", vtpbw = NA_real_)
both_cohorts <- bind_rows(ventilated, no_support)

# SYNTHETIC SITE ONLY: synthetic CLIF mortality is unreliable, so death is simulated
# independently of every exposure (35% by day 60, time to death log-normal with
# median 9 days), as the other supplement scripts do. The run exercises the machinery
# and can show no real effect. Never runs at a real site.
if (grepl("^synthetic_clif", site_name)) {
  message("*** SYNTHETIC SITE: simulated mortality (plumbing only; synthetic CLIF mortality is unreliable). ***")
  set.seed(20260615)
  simulated_death <- rbinom(nrow(both_cohorts), 1L, 0.35)
  simulated_day   <- pmin(pmax(rlnorm(nrow(both_cohorts), log(9), 0.95), 0.04), 60)
  both_cohorts <- both_cohorts %>%
    mutate(deceased = simulated_death, mortality_event_60 = simulated_death,
           surv_time = if_else(simulated_death == 1L, simulated_day, 60))
}

ventilated_log_pfvc_sd <- sd(log(ventilated$pfvc[ventilated$pfvc > 0]), na.rm = TRUE)
both_cohorts <- both_cohorts %>%
  filter(!is.na(pfvc), pfvc > 0, !is.na(sf_ratio), !is.na(sofa_total), !is.na(deceased),
         !is.na(age_at_admission), !is.na(sex_category), !is.na(race_category)) %>%
  group_by(cohort) %>%
  mutate(sf_z = as.numeric(scale(sf_ratio)), sofa_z = as.numeric(scale(sofa_total))) %>%
  ungroup() %>%
  mutate(cohort        = factor(cohort, levels = COHORTS),
         sex_category  = factor(sex_category, levels = c("Male", "Female")),
         race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
         log_pfvc_z    = log(pfvc) / ventilated_log_pfvc_sd,
         # The dose exists only under ventilation: (ventilated) x (VT/PBW - ventilated median),
         # zero in the control. With the cohort main effect in the model, its slope is
         # estimated from ventilated patients alone, and the zero never places a control
         # patient on the dose scale; the centring moves only the ventilated intercept
         # (to the median dose), not the PFVC terms or their cohort contrast.
         vtpbw_ventilated = if_else(cohort == "Ventilated", vtpbw - median(vtpbw[cohort == "Ventilated"], na.rm = TRUE), 0),
         escalation_day   = as.numeric(difftime(escalation_dttm, recorded_dttm, units = "days")))
if (anyNA(both_cohorts$vtpbw_ventilated)) stop("ventilated patients without VT/PBW in the cross-sectional table")

cohort_counts <- both_cohorts %>% group_by(cohort) %>%
  summarise(n_patients = n(), n_deaths = sum(deceased == 1),
            n_escalated = sum(!is.na(escalation_day)), median_age = median(age_at_admission),
            median_sf = median(sf_ratio), median_sofa = median(sofa_total), .groups = "drop")
print(as.data.frame(cohort_counts), row.names = FALSE)
if (any(cohort_counts$n_deaths < MIN_DEATHS) || nrow(cohort_counts) < 2)
  stop("each cohort needs at least ", MIN_DEATHS, " in-hospital deaths")

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
  mutate(scale = "OR per SD of log PFVC (ventilated cohort SD)", site = site_name)
curves <- cohort_fits %>% mutate(curves = map(fit, "curves")) %>% select(-fit) %>% unnest(curves) %>%
  mutate(scale = "log-odds of in-hospital death relative to age 40", site = site_name)

# =============================================================================
# The contrast: the PFVC coefficient, no support minus ventilated
# =============================================================================
# Every covariate free by cohort, so each cohort keeps its own age curve and
# severity slopes; the interaction row is the difference in the PFVC log-OR.
contrast_rhs <- function(adjustment)
  paste("vtpbw_ventilated + cohort * (sf_z + sofa_z + ns(age_at_admission, 4) +",
        if (!is.na(ADJUSTMENTS[[adjustment]])) paste(ADJUSTMENTS[[adjustment]], "+") else "", "log_pfvc_z)")
both_cohorts_cause_specific <- both_cohorts %>%
  mutate(time_cause_specific  = pmax(if_else(!is.na(escalation_day), pmin(surv_time, escalation_day), surv_time), 0.01),
         death_cause_specific = as.integer(mortality_event_60 == 1 & (is.na(escalation_day) | surv_time <= escalation_day)))
contrast_rows <- function(fit, outcome, adjustment, ratio_label) {
  co <- summary(fit)$coefficients
  se_col <- if ("Std. Error" %in% colnames(co)) "Std. Error" else "se(coef)"
  estimate_col <- if ("Estimate" %in% colnames(co)) "Estimate" else "coef"
  terms <- c(ventilated = "log_pfvc_z", difference = "cohortNo support:log_pfvc_z")
  b <- coef(fit)[terms]; V <- vcov(fit)[terms, terms]
  # the control's own coefficient: ventilated + difference
  control_b <- sum(b); control_se <- sqrt(sum(V))
  tibble(quantity = c("ventilated", "no support minus ventilated", "no support"),
         log_ratio = c(co[terms, estimate_col], control_b),
         se = c(co[terms, se_col], control_se)) %>%
    mutate(ratio = exp(log_ratio), ratio_lo = exp(log_ratio - 1.96 * se), ratio_hi = exp(log_ratio + 1.96 * se),
           p = 2 * pnorm(-abs(log_ratio / se)), ratio_type = ratio_label,
           outcome = outcome, adjustment = adjustment)
}
contrast <- map_dfr(names(ADJUSTMENTS), function(adjustment) bind_rows(
  contrast_rows(fit_strict(glm(as.formula(paste("deceased ~", contrast_rhs(adjustment))),
                               family = binomial, data = both_cohorts)),
                "in-hospital death", adjustment, "OR per SD of log PFVC"),
  contrast_rows(fit_strict(coxph(as.formula(paste("Surv(time_cause_specific, death_cause_specific) ~", contrast_rhs(adjustment))),
                                 data = both_cohorts_cause_specific)),
                "60-day death before escalation (cause-specific)", adjustment, "HR per SD of log PFVC"))) %>%
  cross_join(both_cohorts_cause_specific %>% summarise(n_patients = n(), n_deaths = sum(deceased == 1),
                                                       n_deaths_before_escalation = sum(death_cause_specific))) %>%
  mutate(site = site_name)

message("\nPFVC per SD of log PFVC, per cohort, and the share of the 40-to-90 age gradient it absorbs:")
print(as.data.frame(estimates %>% select(cohort, adjustment, or_per_sd, or_lo, or_hi, lr_p,
                                         share_of_age_gradient_absorbed, n_patients, n_deaths) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
message("\nThe contrast (strain predicts the no-support ratio nearer 1, so a difference above 1 for a protective PFVC):")
print(as.data.frame(contrast %>% select(outcome, adjustment, quantity, ratio, ratio_lo, ratio_hi, p) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)

write_csv(mask_small_counts(estimates), file.path(final_dir, paste0("pfvc_age_control_estimates_", site_name, ".csv")))
write_csv(curves, file.path(final_dir, paste0("pfvc_age_control_curves_", site_name, ".csv")))
write_csv(mask_small_counts(contrast), file.path(final_dir, paste0("pfvc_age_control_contrast_", site_name, ".csv")))

# =============================================================================
# Figure: the age curves by cohort, and the PFVC ratios with the contrast
# =============================================================================
curve_panel <- curves %>% filter(adjustment == "adjusted") %>%
  pivot_longer(c(without_pfvc, with_pfvc), names_to = "model", values_to = "log_odds") %>%
  mutate(cohort = factor(cohort, levels = COHORTS)) %>%
  ggplot(aes(age_at_admission, log_odds, colour = model)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_line(linewidth = 0.9) + geom_point(size = 1.5) + facet_wrap(~ cohort) +
  scale_colour_manual(values = OKABE_ITO, labels = c(without_pfvc = "without log PFVC", with_pfvc = "with log PFVC"),
                      name = NULL) +
  labs(title = "A. Age curve of in-hospital death, with and without predicted lung size",
       subtitle = "adjusted; white man at the cohort's median dose, severity and PFVC; relative to age 40",
       x = "Age", y = "log-odds of death") +
  theme_minimal(base_size = 10)
ratio_panel <- contrast %>% filter(adjustment == "adjusted") %>%
  mutate(quantity = factor(quantity, levels = rev(c("ventilated", "no support", "no support minus ventilated")))) %>%
  ggplot(aes(ratio, quantity)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = ratio_lo, xmax = ratio_hi), colour = OKABE_ITO[["with_pfvc"]]) +
  facet_wrap(~ outcome, ncol = 1) + scale_x_log10() +
  labs(title = "B. Per SD of log PFVC, and the cohort contrast",
       subtitle = "strain predicts a protective ventilated ratio\nand a no-support ratio nearer 1",
       x = "OR (in-hospital) or HR (60-day, before escalation), log scale", y = NULL) +
  theme_minimal(base_size = 10)
ggsave(file.path(final_dir, paste0("pfvc_age_control_", site_name, ".pdf")),
       curve_panel + ratio_panel + plot_layout(widths = c(1.6, 1)), width = 14, height = 5)
message("xsec_pfvc_age_control complete -> ", final_dir)
