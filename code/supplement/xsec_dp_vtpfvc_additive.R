# =============================================================================
# Supplement (cross-sectional): is driving pressure a sufficient surrogate for strain?
# Driving pressure and VT/PFVC as additive predictors of mortality, modified by age
# =============================================================================
# A skeptic may argue that driving pressure already measures what VT/PFVC is meant
# to measure: DP = VT x Ers, and elastance scales with the aerated lung, so DP
# should carry the strain whatever the patient's size. That holds only if the
# lung's elastic recoil per unit volume is constant. It is not: recoil falls with
# age, so the same driving pressure distends an older lung further. If so, DP and
# VT/PFVC should each predict death with the other held fixed, and the balance
# between them should shift with age.
#
# The identity behind the test. With DP in cmH2O and VT/PFVC in percent of
# predicted FVC, define the PFVC-normalised (specific) elastance
#     Espec = DP / (VT/PFVC)            cmH2O per percent of predicted FVC.
# Then log DP = log Espec + log VT/PFVC,
# and a model with both logs,
#     b_dp log DP + b_vt log VT/PFVC  =  b_dp log Espec + (b_dp + b_vt) log VT/PFVC,
# weights recoil (Espec) and strain (VT/PFVC) separately. "DP is a sufficient
# surrogate" is b_vt = 0: recoil and strain carry equal weight, so only their
# product matters (Gattinoni's stress = specific elastance x strain). A nonzero
# b_vt says strain matters beyond the pressure it takes to deliver it. The
# likelihood-ratio test of VT/PFVC beyond DP is that test; the reparametrised
# model (log Espec + log VT/PFVC) is reported beside it so the coefficients read
# physiologically.
#
# The rival reading. A single plateau pressure is a noisy measure of the true
# driving pressure (and set PEEP is not total PEEP), so its coefficient is
# attenuated and VT/PFVC, computed from a formula without that noise, can absorb
# signal DP lost. The spline rung does not remove this. An age gradient in the
# balance between the two is harder to manufacture that way, unless the plateau's
# error itself changes with age.
#
# Every mortality model adjusts for VT/PBW (the delivered dose; DP models must), BMI
# (DP includes the chest wall), SOFA and SF ratio, the covariates of scripts 04-05.
# The premise model (specific elastance across age, panel D) is not a mortality
# model and holds no VT/PBW.
# With VT/PBW held, the remaining variation in VT/PFVC is the PBW/PFVC discordance,
# so the VT/PFVC terms here are the discordance at a fixed dose. Each model is fitted
# demographic-adjusted (+ age, sex, race) and unadjusted (demographics dropped).
# The age-modification rungs need age as a main effect, so their unadjusted arm
# keeps age (the modifier) and drops sex and race.
#
# Age as a main effect enters two ways, and every model is fitted both ways
# (column age_adjustment): linear, and a 4-df natural spline. The spline matters
# because VT/PFVC carries the GLI equations' curved age terms; with age linear, a
# VT/PFVC x age product can absorb curvature in mortality by age rather than a
# change in the VT/PFVC effect. The modification itself stays linear in age
# (dp_c:age_c, vt_c:age_c), so a VT/PFVC x age term that survives the spline is a
# change in slope, not unmodelled age. Rungs without age are identical in the two
# forms in the unadjusted arm.
#
# VT/PBW is held fixed in every mortality model, because VT/PFVC only escapes confounding by
# indication at a given VT/PBW: clinicians set the dose from PBW (and lower it for
# sicker patients), while the PBW/PFVC discordance left over is assigned by the
# formulas. So that a straight line in VT/PBW cannot leave dose-severity confounding
# behind, every model is also fitted with VT/PBW as a 3-df natural spline (column
# vtpbw_adjustment; the figure uses the spline).
#
# Predictors are centred logs (dp_c = log DP - mean, vt_c = log VT/PFVC - mean)
# and age_c = age/10 - 6, so a main effect in a model with interactions is the
# effect at the mean DP, the mean VT/PFVC and age 60. Ratios are reported per SD
# of the log predictor, interactions with age per SD per decade.
#
# Rungs (in-hospital death, logistic; 60-day death, Cox, as the companion):
#   dp          log DP
#   vt          log VT/PFVC
#   additive    log DP + log VT/PFVC
#   product     log DP * log VT/PFVC          departure from additivity (logit scale)
#   espec       log Espec + log VT/PFVC       the additive model, reparametrised
#   dp_ns       ns(log DP, 4)                 DP allowed any shape ...
#   dp_ns_vt    ns(log DP, 4) + log VT/PFVC   ... does VT/PFVC still add?
#   age_main    additive + age
#   dp_x_age    DP modified by age
#   vt_x_age    VT/PFVC modified by age
#   both_x_age  both modified by age
#   full        both modified by age, plus DP x VT/PFVC
#   age_band    cell means: DP and VT/PFVC slopes within each age band
#
# Restricted to patients with a measured plateau (dp > 0 at the index timepoint;
# pressures are never forward-filled). A model that warns (non-convergence,
# separation) stops the script: no rung is silently replaced by a simpler one.
#
# Inputs : intermediate/analysis_cross_sectional.parquet (script 03)
# Outputs: final/supplement/
#   dp_vtpfvc_additive_estimates_{site}.csv   exposure terms of every rung: log-scale
#                                             coefficient, ratio per SD, CI, p
#   dp_vtpfvc_additive_lrt_{site}.csv         the likelihood-ratio tests and AICs
#   dp_vtpfvc_additive_age_slopes_{site}.csv  per-SD ratio for DP and VT/PFVC by age:
#                                             along continuous age and within age bands
#   dp_vtpfvc_additive_espec_age_{site}.csv   specific elastance across age (the premise)
#   dp_vtpfvc_additive_{site}.pdf             slopes by age, band estimates, predicted
#                                             mortality, specific elastance by age
# Usage: uvr run code/supplement/xsec_dp_vtpfvc_additive.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(survival)
  library(splines)
  library(patchwork)
})

source("utils/config.R")
# An age band is drawn only with enough data for a stable estimate (not masking):
# panel B needs this many deaths, panel D this many patients. The CLIF minimum-count
# standard.
MIN_DEATHS_TO_DRAW <- 10L
MIN_PATIENTS_TO_DRAW <- 10L
site_name <- config$site_name
final_dir <- final_dir_for("supplement")

AGE_BAND_BREAKS <- c(18, 50, 65, 80, Inf)
AGE_BAND_LABELS <- c("18-49", "50-64", "65-79", "80+")
AGE_CURVE_GRID  <- seq(30, 90, by = 5)
OKABE_ITO <- c(dp = "#D55E00", vt = "#0072B2", espec = "#009E73",
               q25 = "#E69F00", q50 = "#009E73", q75 = "#0072B2")

# =============================================================================
# Data: the mechanics subset of the cross-sectional cohort
# =============================================================================

cross_sectional <- read_parquet(file.path(config$output_dir, "analysis_cross_sectional.parquet"))

# SYNTHETIC SITE ONLY: synthetic CLIF mortality is unreliable (a handful of deaths),
# so death is simulated independently of every exposure (35% by day 60, time to
# death log-normal with median 9 days), as 10_panel_common.R does. The run then
# exercises the machinery and can show no real effect. Never runs at a real site.
if (grepl("^synthetic_clif", site_name)) {
  message("*** SYNTHETIC SITE: simulated mortality (plumbing only; synthetic CLIF mortality is unreliable). ***")
  set.seed(20260615)
  simulated_death <- rbinom(nrow(cross_sectional), 1L, 0.35)
  simulated_day   <- pmin(pmax(rlnorm(nrow(cross_sectional), log(9), 0.95), 0.04), 60)
  cross_sectional <- cross_sectional %>%
    mutate(deceased = simulated_death, mortality_event_60 = simulated_death,
           surv_time = if_else(simulated_death == 1L, simulated_day, 60))
}

mechanics <- cross_sectional %>%
  filter(!is.na(dp), dp > 0) %>%
  mutate(
    sex_category  = factor(sex_category,  levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
    age10 = age_at_admission / 10,
    sf10  = sf_ratio / 10,
    espec = dp / vtpfvc,
    age_band = cut(age_at_admission, AGE_BAND_BREAKS, labels = AGE_BAND_LABELS, right = FALSE)
  )

model_variables <- c("deceased", "dp", "vtpfvc", "vtpbw", "bmi", "sofa_total", "sf10",
                     "age10", "sex_category", "race_category")
complete_mechanics <- mechanics %>% filter(if_all(all_of(model_variables), ~ !is.na(.x)))
message("Patients with a measured plateau: ", nrow(mechanics), "; complete for the models: ",
        nrow(complete_mechanics), " (", nrow(mechanics) - nrow(complete_mechanics), " dropped for missing covariates)")

# Centre on this cohort's means; the SDs set the per-SD scale of the reported ratios.
mean_log_dp     <- mean(log(complete_mechanics$dp))
mean_log_vtpfvc <- mean(log(complete_mechanics$vtpfvc))
mean_log_espec  <- mean(log(complete_mechanics$espec))
sd_log <- c(dp_c = sd(log(complete_mechanics$dp)), vt_c = sd(log(complete_mechanics$vtpfvc)),
            espec_c = sd(log(complete_mechanics$espec)))

centre <- function(df) df %>% mutate(
  dp_c    = log(dp)     - mean_log_dp,
  vt_c    = log(vtpfvc) - mean_log_vtpfvc,
  espec_c = log(espec)  - mean_log_espec,
  age_c   = age10 - 6)
complete_mechanics <- centre(complete_mechanics)

survival_mechanics <- complete_mechanics %>%
  filter(!is.na(surv_time), surv_time > 0, !is.na(mortality_event_60)) %>%
  mutate(event = as.integer(mortality_event_60))

message("Logistic (in-hospital death): ", nrow(complete_mechanics), " patients, ",
        sum(complete_mechanics$deceased == 1), " deaths; Cox (60-day death): ",
        nrow(survival_mechanics), " patients, ", sum(survival_mechanics$event), " deaths")
message("Correlation of log DP with log VT/PFVC: ",
        round(cor(complete_mechanics$dp_c, complete_mechanics$vt_c), 3))

# =============================================================================
# Model rungs
# =============================================================================

# VT/PBW, the delivered dose, enters linearly and as a 3-df natural spline
# (column vtpbw_adjustment). "VTPBW" in the covariate set is replaced by it.
VTPBW_TERM <- c(linear = "vtpbw", spline = "ns(vtpbw, 3)")
base_covariates <- "VTPBW + bmi + sofa_total + sf10"
# The age main effect, by age_adjustment. "AGE" in a rung (and in the demographic
# set) is replaced by it; the age modification terms stay linear in age_c.
AGE_MAIN_EFFECT <- c(linear = "age_c", spline = "ns(age_c, 4)")
demographic_covariates <- "AGE + sex_category + race_category"

rungs <- c(
  dp         = "dp_c",
  vt         = "vt_c",
  additive   = "dp_c + vt_c",
  product    = "dp_c * vt_c",
  espec      = "espec_c + vt_c",
  dp_ns      = "ns(dp_c, 4)",
  dp_ns_vt   = "ns(dp_c, 4) + vt_c",
  age_main   = "dp_c + vt_c + AGE",
  dp_x_age   = "dp_c + vt_c + AGE + dp_c:age_c",
  vt_x_age   = "dp_c + vt_c + AGE + vt_c:age_c",
  both_x_age = "dp_c + vt_c + AGE + dp_c:age_c + vt_c:age_c",
  full       = "dp_c + vt_c + AGE + dp_c:age_c + vt_c:age_c + dp_c:vt_c",
  age_band   = "age_band + dp_c:age_band + vt_c:age_band"
)

# Likelihood-ratio tests: (smaller rung, larger rung, what the test asks).
lrt_specs <- tribble(
  ~test,                          ~smaller,     ~larger,      ~question,
  "vt_beyond_dp",                 "dp",         "additive",   "VT/PFVC predicts death beyond DP (DP is not a sufficient surrogate)",
  "dp_beyond_vt",                 "vt",         "additive",   "DP predicts death beyond VT/PFVC",
  "nonadditivity",                "additive",   "product",    "DP x VT/PFVC product term (departure from additivity, logit/log-hazard scale)",
  "vt_beyond_flexible_dp",        "dp_ns",      "dp_ns_vt",   "VT/PFVC predicts death beyond a 4-df spline in DP",
  "dp_x_age",                     "age_main",   "dp_x_age",   "Age modifies the DP effect",
  "vt_x_age",                     "age_main",   "vt_x_age",   "Age modifies the VT/PFVC effect",
  "dp_x_age_given_vt_x_age",      "vt_x_age",   "both_x_age", "Age modifies DP, with VT/PFVC x age in",
  "vt_x_age_given_dp_x_age",      "dp_x_age",   "both_x_age", "Age modifies VT/PFVC, with DP x age in",
  "age_modification_joint",       "age_main",   "both_x_age", "Age modifies either effect (2 df)",
  "nonadditivity_given_age",      "both_x_age", "full",       "DP x VT/PFVC product, with both age interactions in"
)

# Age enters the unadjusted arm only where a rung needs it as the modifier. The
# band rung's linear form holds age within bands only; its spline form adds the
# spline within them.
rung_formula <- function(outcome_lhs, rung, adjusted, age_adjustment, vtpbw_adjustment) {
  covariates <- if (adjusted) paste(base_covariates, "+", demographic_covariates) else base_covariates
  rhs <- rungs[[rung]]
  if (rung == "age_band" && age_adjustment == "spline") rhs <- paste(rhs, "+ AGE")
  formula_text <- str_replace_all(paste(outcome_lhs, "~", rhs, "+", covariates),
                                  fixed("AGE"), AGE_MAIN_EFFECT[[age_adjustment]]) %>%
    str_replace_all(fixed("VTPBW"), VTPBW_TERM[[vtpbw_adjustment]])
  as.formula(formula_text)
}

# Any warning (non-convergence, separation, an infinite coefficient) stops the run,
# naming the rung: no fallback to a simpler model.
fit_without_warnings <- function(expr, label) {
  withCallingHandlers(expr, warning = function(w)
    stop("Model '", label, "' warned: ", conditionMessage(w), call. = FALSE))
}

fit_rungs <- function(outcome_model, adjusted, age_adjustment, vtpbw_adjustment) {
  adjustment <- if (adjusted) "adjusted" else "unadjusted"
  imap(rungs, function(rhs, rung) {
    label <- paste(outcome_model, adjustment, paste(age_adjustment, "age"), paste(vtpbw_adjustment, "VT/PBW"), rung, sep = " / ")
    if (outcome_model == "logistic") {
      model <- fit_without_warnings(glm(rung_formula("deceased", rung, adjusted, age_adjustment, vtpbw_adjustment),
                                        data = complete_mechanics, family = binomial), label)
      if (!model$converged) stop("Model '", label, "' did not converge", call. = FALSE)
    } else {
      model <- fit_without_warnings(coxph(rung_formula("Surv(surv_time, event)", rung, adjusted, age_adjustment, vtpbw_adjustment),
                                          data = survival_mechanics), label)
    }
    if (anyNA(coef(model))) stop("Model '", label, "' has aliased coefficients", call. = FALSE)
    model
  })
}

model_grid <- crossing(outcome_model = c("logistic", "cox"), adjusted = c(TRUE, FALSE),
                       age_adjustment = names(AGE_MAIN_EFFECT),
                       vtpbw_adjustment = names(VTPBW_TERM)) %>%
  mutate(adjustment = if_else(adjusted, "adjusted", "unadjusted"),
         key = paste(outcome_model, adjustment, age_adjustment, vtpbw_adjustment, sep = "_"))
fitted_rungs <- pmap(select(model_grid, outcome_model, adjusted, age_adjustment, vtpbw_adjustment), fit_rungs) %>%
  set_names(model_grid$key)
grid_row <- function(key) as.list(model_grid[model_grid$key == key, ])

# =============================================================================
# Estimates: exposure terms only (the adjustment covariates are not reported)
# =============================================================================

exposure_term <- function(term) {
  str_detect(term, "dp_c|vt_c|espec_c") & !str_detect(term, "^ns\\(")
}
# Per-SD scale of a term: SD of each log predictor in it; age contributes per decade.
term_scale <- function(term) {
  parts <- str_split(str_remove_all(term, "age_band[^:]*"), ":")[[1]]
  prod(vapply(parts, function(p) if (p %in% names(sd_log)) sd_log[[p]] else 1, numeric(1)))
}

cohort_counts <- function(outcome_model) {
  if (outcome_model == "logistic")
    c(n_patients = nrow(complete_mechanics), n_deaths = sum(complete_mechanics$deceased == 1))
  else c(n_patients = nrow(survival_mechanics), n_deaths = sum(survival_mechanics$event))
}

estimates <- imap_dfr(fitted_rungs, function(models, key) {
  meta <- grid_row(key)
  outcome_model <- meta$outcome_model
  counts <- cohort_counts(outcome_model)
  imap_dfr(models, function(model, rung) {
    coefficient_table <- summary(model)$coefficients
    se_column <- if (outcome_model == "logistic") "Std. Error" else "se(coef)"
    estimate_column <- if (outcome_model == "logistic") "Estimate" else "coef"
    tibble(term = rownames(coefficient_table),
           estimate = coefficient_table[, estimate_column],
           se = coefficient_table[, se_column],
           p = coefficient_table[, "Pr(>|z|)"]) %>%
      filter(exposure_term(term)) %>%
      mutate(scale_per_sd = map_dbl(term, term_scale),
             ratio_per_sd = exp(estimate * scale_per_sd),
             ratio_lo = exp((estimate - 1.96 * se) * scale_per_sd),
             ratio_hi = exp((estimate + 1.96 * se) * scale_per_sd),
             ratio = if (outcome_model == "logistic") "odds ratio" else "hazard ratio",
             rung = rung, .before = 1)
  }) %>%
    mutate(site = site_name,
           outcome = if (outcome_model == "logistic") "in-hospital death" else "60-day death",
           outcome_model = outcome_model, adjustment = meta$adjustment,
           age_adjustment = meta$age_adjustment, vtpbw_adjustment = meta$vtpbw_adjustment,
           n_patients = counts[["n_patients"]], n_deaths = counts[["n_deaths"]], .before = 1)
})

# =============================================================================
# Likelihood-ratio tests and fit
# =============================================================================

lrt_p <- function(smaller, larger) {
  if (inherits(larger, "coxph")) {
    comparison <- anova(smaller, larger)
    c(chisq = comparison[2, "Chisq"], df = comparison[2, "Df"], p = comparison[2, "Pr(>|Chi|)"])
  } else {
    comparison <- anova(smaller, larger, test = "LRT")
    c(chisq = comparison[2, "Deviance"], df = comparison[2, "Df"], p = comparison[2, "Pr(>Chi)"])
  }
}

lrt_table <- imap_dfr(fitted_rungs, function(models, key) {
  meta <- grid_row(key)
  counts <- cohort_counts(meta$outcome_model)
  lrt_specs %>%
    mutate(result = map2(smaller, larger, ~ lrt_p(models[[.x]], models[[.y]])),
           chisq = map_dbl(result, "chisq"), df = map_dbl(result, "df"), p = map_dbl(result, "p"),
           aic_smaller = map_dbl(smaller, ~ AIC(models[[.x]])),
           aic_larger  = map_dbl(larger,  ~ AIC(models[[.x]]))) %>%
    select(-result) %>%
    mutate(site = site_name, outcome_model = meta$outcome_model,
           adjustment = meta$adjustment, age_adjustment = meta$age_adjustment,
           vtpbw_adjustment = meta$vtpbw_adjustment,
           n_patients = counts[["n_patients"]], n_deaths = counts[["n_deaths"]], .before = 1)
})

# The reparametrised model must reproduce the additive model's fit exactly.
walk(fitted_rungs, function(models) {
  if (abs(logLik(models$additive) - logLik(models$espec)) > 1e-6)
    stop("log Espec + log VT/PFVC does not reproduce the additive model's likelihood")
})

# =============================================================================
# Slopes by age: continuous (both_x_age) and within bands (age_band)
# =============================================================================

# Per-SD ratio for a predictor at a given age from the both_x_age rung:
# log ratio = (b_main + b_int * age_c) * SD, delta-method SE from the covariance.
slope_at_age <- function(model, predictor, ages) {
  interaction <- paste0(predictor, ":age_c")
  beta <- coef(model); covariance <- vcov(model)
  age_c <- ages / 10 - 6
  log_ratio <- beta[[predictor]] + beta[[interaction]] * age_c
  se <- sqrt(covariance[predictor, predictor] + age_c^2 * covariance[interaction, interaction] +
             2 * age_c * covariance[predictor, interaction])
  per_sd <- sd_log[[predictor]]
  tibble(age = ages, predictor = predictor,
         ratio_per_sd = exp(log_ratio * per_sd),
         ratio_lo = exp((log_ratio - 1.96 * se) * per_sd),
         ratio_hi = exp((log_ratio + 1.96 * se) * per_sd))
}

band_counts <- bind_rows(
  complete_mechanics %>% group_by(age_band) %>%
    summarise(n_patients = n(), n_deaths = sum(deceased == 1),
              cor_log_dp_log_vtpfvc = cor(dp_c, vt_c), .groups = "drop") %>%
    mutate(outcome_model = "logistic"),
  survival_mechanics %>% group_by(age_band) %>%
    summarise(n_patients = n(), n_deaths = sum(event),
              cor_log_dp_log_vtpfvc = cor(dp_c, vt_c), .groups = "drop") %>%
    mutate(outcome_model = "cox"))

age_slopes <- imap_dfr(fitted_rungs, function(models, key) {
  meta <- grid_row(key)
  outcome_model <- meta$outcome_model
  adjustment <- meta$adjustment
  age_adjustment <- meta$age_adjustment
  vtpbw_adjustment <- meta$vtpbw_adjustment
  continuous <- bind_rows(slope_at_age(models$both_x_age, "dp_c", AGE_CURVE_GRID),
                          slope_at_age(models$both_x_age, "vt_c", AGE_CURVE_GRID)) %>%
    mutate(age_form = "continuous (linear interaction)", age_band = NA_character_)
  band <- estimates %>%
    filter(outcome_model == !!outcome_model, adjustment == !!adjustment,
           age_adjustment == !!age_adjustment, vtpbw_adjustment == !!vtpbw_adjustment,
           rung == "age_band") %>%
    transmute(predictor = str_extract(term, "dp_c|vt_c"),
              age_band = str_remove(str_remove(term, ":?(dp_c|vt_c):?"), "^age_band"),
              ratio_per_sd, ratio_lo, ratio_hi, p) %>%
    left_join(filter(band_counts, outcome_model == !!outcome_model) %>%
                mutate(age_band = as.character(age_band)) %>% select(-outcome_model), by = "age_band") %>%
    mutate(age_form = "age band (cell means)")
  bind_rows(continuous, band) %>%
    mutate(site = site_name, outcome_model = outcome_model, adjustment = adjustment,
           age_adjustment = age_adjustment, vtpbw_adjustment = vtpbw_adjustment,
           predictor = recode(predictor, dp_c = "Driving pressure", vt_c = "VT/PFVC"), .before = 1)
})

# =============================================================================
# The premise: specific elastance (DP per percent of predicted FVC) across age
# =============================================================================
# Not a mortality model: log Espec on an age spline with sex, race, BMI, SOFA and
# SF, with no VT/PBW term.

espec_model <- fit_without_warnings(
  lm(log(espec) ~ ns(age_at_admission, 4) + sex_category + race_category + bmi + sofa_total + sf10,
     data = complete_mechanics), "specific elastance by age")
espec_reference <- tibble(
  age_at_admission = AGE_CURVE_GRID,
  sex_category  = factor("Male",  levels = levels(complete_mechanics$sex_category)),
  race_category = factor("WHITE", levels = levels(complete_mechanics$race_category)),
  bmi = median(complete_mechanics$bmi), sofa_total = median(complete_mechanics$sofa_total),
  sf10 = median(complete_mechanics$sf10))
espec_prediction <- predict(espec_model, newdata = espec_reference, interval = "confidence")
espec_by_age <- bind_rows(
  espec_reference %>% transmute(age = age_at_admission,
                                espec = exp(espec_prediction[, "fit"]),
                                espec_lo = exp(espec_prediction[, "lwr"]),
                                espec_hi = exp(espec_prediction[, "upr"]),
                                summary = "adjusted fit (male; White; median BMI SOFA and SF)"),
  complete_mechanics %>% group_by(age_band) %>%
    summarise(age = median(age_at_admission),
              espec_lo = quantile(espec, 0.25), espec_hi = quantile(espec, 0.75),
              espec = median(espec), n_patients = n(), .groups = "drop") %>%
    mutate(summary = "observed median and IQR by age band", age_band = as.character(age_band))) %>%
  mutate(site = site_name, unit = "cmH2O per percent of predicted FVC", .before = 1)

# =============================================================================
# Write the tables
# =============================================================================

write_csv(estimates,
          file.path(final_dir, paste0("dp_vtpfvc_additive_estimates_", site_name, ".csv")))
write_csv(lrt_table,
          file.path(final_dir, paste0("dp_vtpfvc_additive_lrt_", site_name, ".csv")))
write_csv(age_slopes,
          file.path(final_dir, paste0("dp_vtpfvc_additive_age_slopes_", site_name, ".csv")))
write_csv(espec_by_age,
          file.path(final_dir, paste0("dp_vtpfvc_additive_espec_age_", site_name, ".csv")))

lrt_table %>% filter(outcome_model == "logistic") %>%
  select(adjustment, age_adjustment, vtpbw_adjustment, test, chisq, df, p) %>%
  mutate(across(c(chisq, p), ~ signif(.x, 3))) %>%
  arrange(adjustment, age_adjustment, vtpbw_adjustment) %>%
  print(n = Inf)

# =============================================================================
# Figure (logistic models, VT/PBW as a spline)
# =============================================================================

predictor_colours <- c("Driving pressure" = OKABE_ITO[["dp"]], "VT/PFVC" = OKABE_ITO[["vt"]])
age_adjustment_labels <- c(linear = "age linear", spline = "age spline (4 df)")

slopes_panel <- age_slopes %>%
  filter(outcome_model == "logistic", vtpbw_adjustment == "spline",
         age_form != "age band (cell means)") %>%
  mutate(age_adjustment = age_adjustment_labels[age_adjustment]) %>%
  ggplot(aes(age, ratio_per_sd, colour = predictor, fill = predictor, linetype = age_adjustment)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = ratio_lo, ymax = ratio_hi, group = interaction(predictor, age_adjustment)),
              alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ adjustment) +
  scale_colour_manual(values = predictor_colours) + scale_fill_manual(values = predictor_colours) +
  scale_linetype_manual(values = c("age linear" = "dashed", "age spline (4 df)" = "solid")) +
  scale_y_log10() +
  labs(title = "A. Each predictor, the other held fixed, by age",
       x = "Age (years)", y = "Odds ratio per SD (log scale)", colour = NULL, fill = NULL, linetype = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

band_panel <- age_slopes %>%
  filter(outcome_model == "logistic", vtpbw_adjustment == "spline",
         age_form == "age band (cell means)",
         n_deaths >= MIN_DEATHS_TO_DRAW) %>%
  mutate(age_band = factor(age_band, levels = AGE_BAND_LABELS),
         age_adjustment = age_adjustment_labels[age_adjustment]) %>%
  ggplot(aes(ratio_per_sd, age_band, colour = predictor)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = ratio_lo, xmax = ratio_hi), position = position_dodge(width = 0.5)) +
  facet_grid(~ adjustment + age_adjustment) +
  scale_colour_manual(values = predictor_colours) + scale_x_log10() +
  labs(title = "B. Within age bands", x = "Odds ratio per SD (log scale)", y = "Age band",
       colour = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

# Predicted in-hospital mortality across DP at three VT/PFVC levels, one facet per
# age band at that band's median age, from the full adjusted model with age as a
# spline and VT/PBW as a spline. Other covariates at the reference profile (male, White, median VT/PBW,
# BMI, SOFA, SF). DP is drawn only over its 5th-95th percentile within the band.
full_adjusted <- fitted_rungs$logistic_adjusted_spline_spline$full
vtpfvc_levels <- quantile(complete_mechanics$vtpfvc, c(0.25, 0.5, 0.75))
reference_profile <- tibble(
  vtpbw = median(complete_mechanics$vtpbw), bmi = median(complete_mechanics$bmi),
  sofa_total = median(complete_mechanics$sofa_total), sf10 = median(complete_mechanics$sf10),
  sex_category  = factor("Male",  levels = levels(complete_mechanics$sex_category)),
  race_category = factor("WHITE", levels = levels(complete_mechanics$race_category)))
prediction_grid <- complete_mechanics %>% group_by(age_band) %>%
  summarise(age = median(age_at_admission), dp_low = quantile(dp, 0.05), dp_high = quantile(dp, 0.95),
            .groups = "drop") %>%
  mutate(dp = map2(dp_low, dp_high, ~ seq(.x, .y, length.out = 40))) %>% unnest(dp) %>%
  crossing(vtpfvc_quantile = names(vtpfvc_levels)) %>%
  mutate(vtpfvc = vtpfvc_levels[vtpfvc_quantile], age10 = age / 10) %>%
  crossing(reference_profile) %>%
  mutate(espec = dp / vtpfvc) %>% centre()
link_prediction <- predict(full_adjusted, newdata = prediction_grid, type = "link", se.fit = TRUE)
prediction_grid <- prediction_grid %>%
  mutate(risk = plogis(link_prediction$fit),
         risk_lo = plogis(link_prediction$fit - 1.96 * link_prediction$se.fit),
         risk_hi = plogis(link_prediction$fit + 1.96 * link_prediction$se.fit),
         vtpfvc_label = paste0("VT/PFVC ", vtpfvc_quantile, " (", round(vtpfvc, 1), "%)"),
         band_label = paste0(age_band, " (age ", round(age), ")"))
vtpfvc_colours <- setNames(OKABE_ITO[c("q25", "q50", "q75")], sort(unique(prediction_grid$vtpfvc_label)))

risk_panel <- ggplot(prediction_grid, aes(dp, risk, colour = vtpfvc_label, fill = vtpfvc_label)) +
  geom_ribbon(aes(ymin = risk_lo, ymax = risk_hi), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ band_label, nrow = 1) +
  scale_colour_manual(values = vtpfvc_colours) + scale_fill_manual(values = vtpfvc_colours) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "C. Predicted in-hospital mortality: same DP, different strain",
       x = "Driving pressure (cmH2O)", y = "Predicted mortality", colour = NULL, fill = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

espec_panel <- ggplot() +
  geom_ribbon(data = filter(espec_by_age, str_starts(summary, "adjusted")),
              aes(age, ymin = espec_lo, ymax = espec_hi), fill = OKABE_ITO[["espec"]], alpha = 0.15) +
  geom_line(data = filter(espec_by_age, str_starts(summary, "adjusted")),
            aes(age, espec), colour = OKABE_ITO[["espec"]], linewidth = 0.9) +
  geom_pointrange(data = filter(espec_by_age, str_starts(summary, "observed"), n_patients >= MIN_PATIENTS_TO_DRAW),
                  aes(age, espec, ymin = espec_lo, ymax = espec_hi), colour = "grey30") +
  labs(title = "D. Pressure per unit of strain, by age",
       x = "Age (years)", y = "DP / (VT/PFVC)\n(cmH2O per % predicted FVC)") +
  theme_minimal(base_size = 11)

figure <- slopes_panel / band_panel / risk_panel / espec_panel +
  plot_layout(heights = c(1, 1, 1, 0.8)) +
  plot_annotation(
    title = paste0("Driving pressure and VT/PFVC as additive predictors of mortality (", site_name, ")"),
    subtitle = paste("Logistic models; VT/PBW (3-df spline), BMI, SOFA and SF ratio held fixed. A-B: both adjustment sets, age linear and spline.",
                     "C: demographic-adjusted, age spline. D: no VT/PBW."))
ggsave(file.path(final_dir, paste0("dp_vtpfvc_additive_", site_name, ".pdf")),
       figure, width = 12, height = 15)

message("Wrote dp_vtpfvc_additive_* to ", final_dir)
