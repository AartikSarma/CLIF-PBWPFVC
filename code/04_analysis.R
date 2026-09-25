# =============================================================================
# Script 04: Cross-sectional analysis -- demographic bias, mechanics, mortality
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Pipeline script (run after 03). Written for the ventilated analytic cohort:
# run it with PBWPFVC_COHORT unset (imv). The control cohorts reach the paper
# through scripts 21-29, not through this script.
#
# Question (figures 1-3, Tables 1-2): predicted body weight (PBW) sets the tidal
# volume, but predicted FVC (PFVC, race-specific GLI-2012) is the better estimate of
# lung size. Does PBW mis-size the lung differently by age, sex and race, and, at a
# fixed delivered dose (VT/PBW), do the mis-sizing ratio PBW/PFVC and predicted lung
# size itself go with lung mechanics, ventilator liberation and death?
#
# Populations
#   analysis_cross_sectional  hypoxemic, ventilated adults at VT/PBW 6-8 mL/kg, one
#                             row per patient (script 03); every section except:
#   analysis_broad_pfvc       every eligible patient with height, age, sex, race and
#                             PFVC, ventilated or not (4f3)
#   analysis_negative_control ventilated non-hypoxemic and non-ventilated cohorts (4j)
#
# Model families (section)
#   4c  logistic, in-hospital death              4d  linear, lung mechanics (Ers, Crs,
#   4d2 Fine-Gray, liberation by day 28              DP, Ers x PBW/PFVC, MP and MP/size)
#   4f  Cox, 60-day all-cause death              4f2 demographic bias of each metric
#   4f3 PFVC on PBW + demographics               4g2 E-values, spline-age sensitivity
#   4j  negative-control cohorts                 4k  saturated log model (VT, PBW, PFVC;
#                                                    unadjusted) and size-term forms
#   4l  delivered strain inside the band (figure 2)
# Covariates in 4c-4f: unadjusted = SOFA + SF ratio; adjusted = unadjusted + age10 +
# sex + race; either + BMI when the outcome or exposure is driving-pressure derived
# (DP_DERIVED). Every 4c, 4d, 4d2 and 4f model is reported both ways (adjustment
# column). Sections 4f2, 4f3, 4j and 4k state their own covariate sets.
# Exposure specs (4c, 4d, 4d2): VT/PFVC; VT/PBW; VT/PFVC + VT/PBW; VT/PBW + PFVC;
# VT/PBW + log PFVC; VT/PBW + PBW/PFVC; VT/PBW + VT excess (mL). The Cox models (4f)
# carry three: VT/PBW with PBW/PFVC, with PFVC, and with log PFVC.
#
# Environment: PBWPFVC_NC_MIN_EVENTS (default 10) is the minimum number of deaths
# for a 4j or 4k model; lower it only to test the plumbing on a small site.
#
# Inputs : analysis_cross_sectional, analysis_broad_pfvc, analysis_negative_control
#          (.parquet, script 03); final/cross_sectional/attrition_log_<site>.csv (03)
# Outputs: final/cross_sectional/, each file suffixed _<site>, pooled across sites
#          by the coordinator's pooling script (not distributed with this code):
#   table1_                          Table 1 (4a)
#   regression_results_long_         every model's coefficients, adjusted and
#                                    unadjusted (Table 2 and the forest plots)
#   aic_comparison_all_              evidence ratios vs VT/PBW alone (4e)
#   evalues_                         E-values and spline-age estimates (4g2)
#   consort_diagram_ (pdf)           inclusion flow (4i, site QC)
#   negative_control_, negative_control_{counts, identifying_variation,
#     interaction}_                  the negative-control cohorts (4j)
#   size_saturated_log_model_, size_functional_form_   (4k)
#   dose_vtpfvc_histograms_, dose_variance_decomposition_   figure 2 (4l;
#                                    code/pooling/pooled_displays.R)
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(survival)
library(splines)
library(broom)
library(EValue)

source("utils/config.R")
source("utils/consort_diagram.R")   # supplies render_consort() (section 4i)
site_name <- config$site_name

output_dir <- config$output_dir
final_dir <- final_dir_for("cross_sectional")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

# Minimum counts per model (the CLIF consortium minimum-count standard)
MIN_DP_OBS             <- 30   # driving-pressure observations, Static DP demographic-bias model (4f2)
NC_MIN_PATIENTS        <- 50   # patients per cohort, negative-control models (4j)
NC_MIN_EVENTS          <- as.integer(Sys.getenv("PBWPFVC_NC_MIN_EVENTS", "10"))   # deaths per model (4j, 4k)
SATURATED_MIN_PATIENTS <- 100  # complete patients, saturated log model (4k)

# Dose thresholds summarised for each negative-control cohort (4j)
VTPBW_BAND_MAX  <- 8    # mL/kg, upper edge of the 6-8 mL/kg lung-protective band
VTPFVC_ARMA_P75 <- 11   # % of predicted FVC, about the 75th percentile of VT/PFVC in the ARMA low tidal volume arm

# Age and height strata of figure 2 (4l), left-closed. Reason for these cut
# points: author to supply.
DOSE_AGE_BREAKS    <- c(18, 40, 50, 60, 70, 80, Inf)
DOSE_HEIGHT_BREAKS <- c(150, 160, 170, 180, 190, 210)

# =============================================================================
# Load data
# =============================================================================

cross_sectional <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"))

message("Cross-sectional: ", nrow(cross_sectional), " observations")

# Reference categories: Male and White. Set in script 03, but re-applied here so
# the reference is explicit and robust to the parquet round-trip; all regressions
# below report effects relative to male / white patients.
cross_sectional <- cross_sectional %>%
  mutate(
    sex_category  = factor(sex_category,  levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
    # Per-10-unit covariates so the adjusted age and SF-ratio coefficients are
    # reported per 10 years / 10 SF units (Table 1 keeps the raw scales).
    age10 = age_at_admission / 10,
    sf10  = sf_ratio / 10,
    height10 = height_cm / 10,
    # the primary size term of the mortality models (Table 2): log PFVC fits better
    # than PFVC in litres (4k), and its coefficient is per proportional change in
    # predicted lung size, the scale the GLI pieces are read on
    log_pfvc = log(pfvc)
  )

# =============================================================================
# 4a. Table 1
# =============================================================================

# The whole cohort, one row per characteristic (continuous: median and quartiles) or
# per level (categorical: count and percent). The pooled Table 1 is the cohort by
# site, so the table is not stratified. Small cells are not masked yet: a
# deterministic masking step will be applied once the outputs are locked.
TABLE1_CONTINUOUS <- c(age_at_admission = "Age (years)", height_cm = "Height (cm)",
                       pbw = "PBW (kg)", pfvc = "PFVC (L)", pbwpfvc = "PBW/PFVC (kg/L)",
                       vtpbw = "VT/PBW (mL/kg)", vtpfvc = "VT/PFVC (%)",
                       sofa_total = "SOFA", sf_ratio = "SF ratio", crs = "Compliance (mL/cmH2O)",
                       ers = "Elastance (cmH2O/L)", vfd_28 = "28-day VFDs")
TABLE1_CATEGORICAL <- c(sex_category = "Sex", race_category = "Race",
                        deceased = "In-hospital death", mortality_event_60 = "Death by day 60")
n_cohort <- nrow(cross_sectional)
table1_continuous <- imap_dfr(TABLE1_CONTINUOUS, function(label, v) {
  x <- cross_sectional[[v]]
  tibble(characteristic = label, level = NA_character_,
         n_patients = sum(!is.na(x)), percent = NA_real_,
         median = median(x, na.rm = TRUE), q1 = quantile(x, 0.25, na.rm = TRUE, names = FALSE),
         q3 = quantile(x, 0.75, na.rm = TRUE, names = FALSE))
})
table1_categorical <- imap_dfr(TABLE1_CATEGORICAL, function(label, v) {
  cross_sectional %>% filter(!is.na(.data[[v]])) %>%
    count(level = as.character(.data[[v]]), name = "n_patients") %>%
    transmute(characteristic = label, level, n_patients, percent = 100 * n_patients / n_cohort,
              median = NA_real_, q1 = NA_real_, q3 = NA_real_)
})
table1_long <- bind_rows(tibble(characteristic = "Patients", level = NA_character_, n_patients = n_cohort),
                         table1_continuous, table1_categorical) %>%
  mutate(site = site_name, .before = 1)
write_csv(table1_long, file.path(final_dir, paste0("table1_", site_name, ".csv")))
message("Table 1 written (", n_cohort, " patients)")

# =============================================================================
# 4b. Model specifications
# =============================================================================

# Common covariates for all models. The unadjusted set drops the demographics
# (age/sex/race) -- deterministic parents of the PBW/PFVC exposures -- and keeps
# illness severity; both adjusted and unadjusted estimates are reported per the
# adjusted+unadjusted convention.
covariates       <- "race_category + age10 + sex_category + sofa_total + sf10"
covariates_unadj <- "sofa_total + sf10"

# As in the original paper (PMC12313249), any model whose outcome OR exposure is
# derived from driving pressure (static DP, elastance, compliance, the
# elastance-normalized Ers x PBW / Ers x PFVC, and mechanical power with its
# PBW-, PFVC- and Crs-normalized and elastic-component variants) is additionally
# adjusted for BMI. Models without a driving-pressure component use the standard
# covariate set.
DP_DERIVED <- c("dp", "ers", "crs", "ers_pbw", "ers_pfvc",
                "mechanical_power", "mp_pbw", "mp_pfvc", "mp_crs",
                "mp_elastic", "mp_elastic_pbw", "mp_elastic_pfvc")
uses_dp <- function(...) {
  vars <- trimws(unlist(strsplit(paste(c(...), collapse = " + "), "\\+")))
  any(vars %in% DP_DERIVED)
}
model_covariates <- function(..., adjusted = TRUE) {
  base <- if (adjusted) covariates else covariates_unadj
  if (uses_dp(...)) paste(base, "+ bmi") else base
}

# Define all exposure specifications
# vt_excess_ml (script 03) is the DIFFERENCE parameterization of the same PBW-vs-PFVC
# disagreement the pbwpfvc ratio carries: the millilitres of tidal volume the PBW rule
# prescribes over the PFVC rule at the protective dose, from two EXTERNAL (ARMA)
# constants, so nothing is scaled within the cohort. It is not a restatement of the
# ratio -- the ratio treats a 10% mis-sizing alike in a 40 kg and a 90 kg patient, the
# difference does not -- and it is the quantity a clinician can act on (mL on the vent).
exposure_specs <- list(
  vtpfvc        = "vtpfvc",
  vtpbw         = "vtpbw",
  vtpfvc_vtpbw  = "vtpfvc + vtpbw",
  vtpbw_pfvc    = "vtpbw + pfvc",
  vtpbw_logpfvc = "vtpbw + log_pfvc",
  vtpbw_pbwpfvc = "vtpbw + pbwpfvc",
  vtpbw_excess  = "vtpbw + vt_excess_ml"
)

# PFVC enters twice: in litres, the scale of the mechanics betas (Crs per litre of
# PFVC, Claim 3), and as log PFVC, the primary size term of the mortality models.
exposure_labels <- c(
  vtpfvc        = "VT/PFVC",
  vtpbw         = "VT/PBW",
  vtpfvc_vtpbw  = "VT/PFVC + VT/PBW",
  vtpbw_pfvc    = "VT/PBW + PFVC",
  vtpbw_logpfvc = "VT/PBW + log PFVC",
  vtpbw_pbwpfvc = "VT/PBW + PBW/PFVC",
  vtpbw_excess  = "VT/PBW + VT excess (mL)"
)

# =============================================================================
# 4c. Logistic regression — mortality
# =============================================================================

has_mortality_variation <- length(unique(na.omit(cross_sectional$deceased))) > 1

if (has_mortality_variation) {
  mortality_models <- map(exposure_specs, ~ {
    formula_str <- paste("deceased ~", .x, "+", model_covariates(.x))
    glm(as.formula(formula_str), data = cross_sectional, family = binomial)
  })

  message("Logistic regression (mortality) — AIC:")
  iwalk(mortality_models, ~ message("  ", exposure_labels[.y], ": ", round(AIC(.x), 1)))
} else {
  mortality_models <- NULL
  message("Skipping mortality regression: no variation in outcome (all deceased = ",
          unique(na.omit(cross_sectional$deceased)), ")")
}

# =============================================================================
# 4d. Linear regression — elastance, compliance, driving pressure, and the
#     PBW- and PFVC-normalized Ers and mechanical power
# =============================================================================

continuous_outcomes <- list(
  ers      = list(var = "ers",              label = "Elastance"),
  crs      = list(var = "crs",              label = "Compliance"),
  dp       = list(var = "dp",               label = "Static DP"),
  ers_pbw  = list(var = "ers_pbw",          label = "Ers x PBW"),
  ers_pfvc = list(var = "ers_pfvc",         label = "Ers x PFVC"),
  mp       = list(var = "mechanical_power", label = "Mechanical power"),
  mp_pbw   = list(var = "mp_pbw",           label = "MP / PBW"),
  mp_pfvc  = list(var = "mp_pfvc",          label = "MP / PFVC"),
  mp_crs   = list(var = "mp_crs",           label = "MP / Crs"),
  # elastic tidal power per PFVC: the "specific elastic power" sensitivity to total MP
  mp_el_pfvc = list(var = "mp_elastic_pfvc", label = "Elastic MP / PFVC")
)

continuous_models <- list()

for (outcome_name in names(continuous_outcomes)) {
  outcome_var <- continuous_outcomes[[outcome_name]]$var
  outcome_label <- continuous_outcomes[[outcome_name]]$label

  covars_used <- model_covariates(outcome_var)
  models_for_outcome <- map(exposure_specs, ~ {
    formula_str <- paste(outcome_var, "~", .x, "+", covars_used)
    lm(as.formula(formula_str), data = cross_sectional)
  })

  continuous_models[[outcome_name]] <- models_for_outcome

  message(outcome_label, " models — AIC:")
  iwalk(models_for_outcome, ~ message("  ", exposure_labels[.y], ": ", round(AIC(.x), 1)))
}

# =============================================================================
# 4d2. 28-day VFDs — competing-risks (Fine-Gray) regression
# =============================================================================
# Per Yehya & Harhay (AJRCCM 2019), VFDs are analyzed as a competing-risks
# outcome on script 03's clock, which starts at the index: event of interest =
# liberation (vfd_status 1: the last IMV record of the stay falls within 28 days
# and the patient is alive then), competing risk = death while still ventilated
# (vfd_status 2), censored at day 28 if still ventilated (vfd_status 0). A death
# after liberation leaves the patient liberated. Each exposure specification is fit
# with a Fine-Gray subdistribution hazard model (survival::finegray weights +
# coxph), so effects are subdistribution hazard ratios for liberation (SHR > 1 =
# faster liberation). The mortality component is modeled in section 4c.
vfd_cr_covariates <- covariates  # VFD + these exposures are never DP-derived (no BMI)

# The fit carries the patients and liberations it used as attributes: finegray()
# repeats patients across risk-set rows, so nobs() of the coxph fit counts rows.
fit_vfd_finegray <- function(exposure_spec, covariate_rhs = vfd_cr_covariates) {
  model_rhs <- paste(exposure_spec, "+", covariate_rhs)
  rhs_vars  <- all.vars(as.formula(paste("~", model_rhs)))
  df <- cross_sectional %>%
    mutate(vfd_status_f = factor(vfd_status, levels = c(0, 1, 2),
                                 labels = c("censored", "extubation", "death"))) %>%
    select(vfd_time, vfd_status_f, all_of(rhs_vars)) %>%
    filter(!is.na(vfd_time), vfd_time > 0, !is.na(vfd_status_f)) %>%
    drop_na(all_of(rhs_vars))
  fg <- survival::finegray(survival::Surv(vfd_time, vfd_status_f) ~ ., data = df,
                           etype = "extubation")
  fit <- survival::coxph(
    as.formula(paste0("survival::Surv(fgstart, fgstop, fgstatus) ~ ", model_rhs)),
    weights = fgwt, data = fg
  )
  attr(fit, "n_patients") <- nrow(df)
  attr(fit, "n_events")   <- sum(df$vfd_status_f == "extubation")
  fit
}

# Patients and events a fitted model used: the Fine-Gray attributes above, the
# subjects and deaths of a Cox fit (nobs() of a coxph fit is its event count),
# the rows and deaths of a logistic fit, the rows of a linear fit.
model_patients <- function(model) {
  if (!is.null(attr(model, "n_patients"))) attr(model, "n_patients")
  else if (inherits(model, "coxph")) model$n
  else stats::nobs(model)
}
model_events <- function(model) {
  if (!is.null(attr(model, "n_events"))) attr(model, "n_events")
  else if (inherits(model, "coxph")) model$nevent
  else if (inherits(model, "glm") && family(model)$family == "binomial") sum(model$y)
  else NA_real_
}

vfd_cr_models <- map(exposure_specs, fit_vfd_finegray)

message("28-day VFDs (Fine-Gray, extubation SHR) — AIC:")
iwalk(vfd_cr_models, ~ message("  ", exposure_labels[.y], ": ", round(AIC(.x), 1)))

# =============================================================================
# 4e. AIC comparison across all models and outcomes
# =============================================================================

aic_results <- list()

# Two DIAGNOSTIC rungs, fit here only (not in the regression tables), to read the
# gap between "VT/PBW + PFVC" and the ratio models. log PFVC = log PBW - log(PBW/PFVC)
# and log PBW = f(height, sex), so the PFVC model is the ratio model PLUS a height
# main effect, which the adjustment set deliberately omits. Rung 1 adds height to the
# ratio model: if its AIC matches the PFVC model, the whole gap is height. Rung 2
# swaps PFVC for FVC_age25 (age pinned to 25): how much of the PFVC fit is the
# structural (height/sex/race) leg versus the age slope.
aic_extra_specs <- c(
  "VT/PBW + PBW/PFVC + height" = "vtpbw + pbwpfvc + height10",
  "VT/PBW + FVC_age25"         = "vtpbw + pfvc_age25"
)
# n_obs is the patients each model used (a missing covariate, BMI in the DP-derived
# models, drops a patient). AIC compares models only on the same patients, so an
# outcome whose models used different numbers is announced.
aic_row <- function(models, extra_models, outcome_label) {
  all_models <- c(models, extra_models)
  n_used <- map_dbl(all_models, model_patients)
  if (n_distinct(n_used) > 1)
    message("AIC comparison, ", outcome_label, ": the models used different numbers of patients (",
            paste(sort(unique(n_used)), collapse = ", "), "); their AICs are not on one sample")
  tibble(exposure = c(exposure_labels, names(aic_extra_specs)),
         AIC = map_dbl(all_models, AIC),
         n_obs = unname(n_used)) %>%
    mutate(is_reference = exposure == "VT/PBW")
}

# Mortality
if (!is.null(mortality_models)) {
  mort_extra <- map(aic_extra_specs, ~ glm(as.formula(paste("deceased ~", .x, "+", model_covariates(.x))),
                                          data = cross_sectional, family = binomial))
  aic_results[["Mortality"]] <- aic_row(mortality_models, mort_extra, "Mortality")
}

# Continuous outcomes. Exclude the normalized-mechanics outcomes (Ers x PBW/PFVC,
# MP/PBW/PFVC) from the AIC comparison/heatmap: the outcome shares a predicted-size
# variable with the size-containing exposures, so their evidence ratios are inflated
# by the shared denominator rather than a dosing-outcome relationship. Their
# standalone regression tables and long-format results are still produced.
AIC_EXCLUDE <- c("ers_pbw", "ers_pfvc", "mp_pbw", "mp_pfvc", "mp_el_pfvc")
for (outcome_name in setdiff(names(continuous_outcomes), AIC_EXCLUDE)) {
  outcome_var <- continuous_outcomes[[outcome_name]]$var
  cont_extra <- map(aic_extra_specs, ~ lm(as.formula(paste(outcome_var, "~", .x, "+", model_covariates(outcome_var))),
                                         data = cross_sectional))
  aic_results[[continuous_outcomes[[outcome_name]]$label]] <-
    aic_row(continuous_models[[outcome_name]], cont_extra, continuous_outcomes[[outcome_name]]$label)
}

# 28-day VFDs (competing-risks Fine-Gray models). AICs are comparable within the
# outcome and referenced to the VT/PBW model, as elsewhere.
vfd_extra <- map(aic_extra_specs, fit_vfd_finegray)
aic_results[["28-day VFDs"]] <- aic_row(vfd_cr_models, vfd_extra, "28-day VFDs")

# Evidence ratios are all referenced to the VT/PBW-alone model WITHIN each
# outcome: ER = exp(-0.5 * (AIC_model - AIC_VT/PBW)). ER > 1 means more support
# than VT/PBW alone, ER = 1 for VT/PBW itself. Truncated to [0.001, 1000].
ER_FLOOR <- 0.001
ER_CEIL  <- 1000
aic_all <- bind_rows(aic_results, .id = "outcome") %>%
  group_by(outcome) %>%
  mutate(
    aic_ref = AIC[is_reference][1],
    delta_AIC = AIC - aic_ref,
    evidence_ratio = exp(-0.5 * delta_AIC),
    evidence_ratio_trunc = pmin(pmax(evidence_ratio, ER_FLOOR), ER_CEIL),
    er_label = case_when(
      evidence_ratio > ER_CEIL  ~ ">1000",
      evidence_ratio < ER_FLOOR ~ "<0.001",
      TRUE                      ~ formatC(evidence_ratio, format = "g", digits = 2)
    ),
    # Sample-size-standardized: delta_AIC per 1000 patients and its evidence ratio,
    # so the figure separates models the raw ratio saturates at 1000 and cohorts of
    # different size are comparable (the coordinator's pooling script, not
    # distributed with this code, uses the same scale).
    delta_AIC_per_1k = delta_AIC / n_obs * 1000,
    er_per_1k        = exp(-0.5 * delta_AIC_per_1k),
    er_per_1k_trunc  = pmin(pmax(er_per_1k, ER_FLOOR), ER_CEIL),
    er_per_1k_label  = case_when(
      er_per_1k > ER_CEIL  ~ ">1000",
      er_per_1k < ER_FLOOR ~ "<0.001",
      TRUE                 ~ formatC(er_per_1k, format = "g", digits = 2))
  ) %>%
  ungroup() %>%
  # Column order: clinical outcomes first, then compliance/elastance, then static
  # driving pressure and mechanical power; any other outcome label keeps its place
  # at the end instead of becoming NA.
  mutate(outcome = factor(outcome, levels = {
    known <- c("Mortality", "28-day VFDs", "Compliance", "Elastance", "Static DP", "Mechanical power")
    c(intersect(known, unique(outcome)), setdiff(unique(outcome), known))
  }))

# Row order: exposure specs by the strongest evidence ratio they reach in any
# outcome (best at the top). Raw AIC is not comparable across outcomes, so the
# AIC-based evidence ratio (truncated) is used; ties are broken by the summed
# log10 ratio.
exposure_order <- aic_all %>%
  group_by(exposure) %>%
  summarise(max_er = max(evidence_ratio_trunc, na.rm = TRUE),
            sum_lr = sum(log10(evidence_ratio_trunc), na.rm = TRUE),
            .groups = "drop") %>%
  arrange(max_er, sum_lr) %>%
  pull(exposure)
aic_all <- aic_all %>% mutate(exposure = factor(exposure, levels = exposure_order))

message("AIC comparison (evidence ratios vs VT/PBW-alone within each outcome):")
print(aic_all)

# The evidence-ratio heatmap across sites is drawn from this table by the
# coordinator's pooling script (not distributed with this code).
write_csv(aic_all, file.path(final_dir, paste0("aic_comparison_all_", site_name, ".csv")))

# =============================================================================
# 4f. Survival analysis
# =============================================================================

# Event = all-cause death within 60 days of the index, in- or out-of-hospital,
# derived in script 03 (mortality_event_60, surv_time in days from the index) from
# the patient-level death_dttm, dated at discharge for an expired patient with no
# death time. Using this
# instead of the in-hospital-only `deceased` flag stops survivors from being
# censored at hospital discharge and counts post-discharge deaths as events.
surv_data <- cross_sectional %>%
  filter(!is.na(surv_time), surv_time > 0, !is.na(mortality_event_60)) %>%
  mutate(event = as.integer(mortality_event_60))

n_deaths <- sum(surv_data$event, na.rm = TRUE)
message("Survival analysis: ", nrow(surv_data), " patients, ", n_deaths,
        " deaths within 60 days (in- and out-of-hospital), ",
        nrow(surv_data) - n_deaths, " censored at day 60")

if (n_deaths > 0 && length(unique(surv_data$event)) > 1) {
  cox_model <- coxph(
    Surv(surv_time, event) ~ pbwpfvc + vtpbw + age10 +
      sex_category + race_category + sf10 + sofa_total,
    data = surv_data
  )

  # Companion model with PFVC as the scaling exposure (mirrors the VT/PBW + PFVC
  # mortality model), so the survival analysis carries both headline exposures.
  cox_model_pfvc <- coxph(
    Surv(surv_time, event) ~ pfvc + vtpbw + age10 +
      sex_category + race_category + sf10 + sofa_total,
    data = surv_data
  )

  # the primary size term, log PFVC (Table 2)
  cox_model_logpfvc <- coxph(
    Surv(surv_time, event) ~ log_pfvc + vtpbw + age10 +
      sex_category + race_category + sf10 + sofa_total,
    data = surv_data
  )

  message("Cox model (PBW/PFVC):")
  print(summary(cox_model))
  message("Cox model (PFVC):")
  print(summary(cox_model_pfvc))
  message("Cox model (log PFVC):")
  print(summary(cox_model_logpfvc))
} else {
  message("Skipping survival analysis: no mortality events in data")
}

# =============================================================================
# 4f2. Demographic-bias models (outcome ~ age + sex + race + height + SF + SOFA)
# =============================================================================
# How each dosing / driving-pressure / elastance metric varies by demographics
# (the algorithmic-bias view). Every model adjusts for age, sex, race, height, SF
# ratio and SOFA; driving-pressure-derived outcomes add BMI. These models are
# descriptive, not causal: the coefficients of interest are age, sex and race, and
# height is held fixed so that they compare patients of the same height rather than
# proxy for height differences between groups. (4g2's over-adjustment argument
# concerns the exposure-mortality models, where height is part of the exposure's
# identifying variation.) Predictors are scaled per 10 units (age, height,
# SF ratio). Outcomes for the dosing and elastance-normalized metrics are
# z-scored WITHIN SITE (so standardized betas are comparable across cohorts in SD
# units, and no row-level data is needed to pool). Static DP is on the raw cmH2O
# scale and additionally adjusts for BMI; Mortality is logistic (OR). Reference
# categories are Male and White.

demo_data <- cross_sectional %>%
  mutate(age10 = age_at_admission / 10,
         height10 = height_cm / 10,
         sf10 = sf_ratio / 10)

demo_covars     <- "age10 + sex_category + race_category + height10 + sf10 + sofa_total"
demo_covars_bmi <- paste(demo_covars, "+ bmi")

# Each entry: fitted model + metadata for the unified long table.
demo_models <- list()

# z-scored linear outcomes. The elastance- and power-normalized outcomes are
# driving-pressure derivatives, so they are adjusted for BMI (matching the
# original paper, PMC12313249); VT/PBW and VT/PFVC are not. MP/PBW vs MP/PFVC is the
# mechanical-power parallel of the Ers x PBW vs Ers x PFVC pair: the sign of the
# demographic term should flip with the denominator in the same way.
demo_z_outcomes <- c(vtpbw = "VT/PBW", vtpfvc = "VT/PFVC (%)",
                     ers_pbw = "Ers x PBW", ers_pfvc = "Ers x PFVC",
                     mp_pbw = "MP / PBW", mp_pfvc = "MP / PFVC", mp_crs = "MP / Crs")
for (v in names(demo_z_outcomes)) {
  cov_v <- if (uses_dp(v)) demo_covars_bmi else demo_covars
  fstr <- paste0("scale(", v, ") ~ ", cov_v)
  demo_models[[demo_z_outcomes[[v]]]] <- list(
    model = lm(as.formula(fstr), data = demo_data),
    type = "Beta", family = "linear", formula = fstr,
    outcome_scale = "z-scored within site (SD units)"
  )
}

# Mortality logistic (no BMI)
if (has_mortality_variation) {
  fstr <- paste("deceased ~", demo_covars)
  demo_models[["Mortality"]] <- list(
    model = glm(as.formula(fstr), data = demo_data, family = binomial),
    type = "OR", family = "logistic", formula = fstr,
    outcome_scale = "natural units"
  )
}

# Static DP raw linear + BMI (guard on minimum driving-pressure N)
if (sum(!is.na(demo_data$dp)) >= MIN_DP_OBS) {
  fstr <- paste("dp ~", demo_covars_bmi)
  demo_models[["Static DP"]] <- list(
    model = lm(as.formula(fstr), data = demo_data),
    type = "Beta", family = "linear", formula = fstr,
    outcome_scale = "natural units"
  )
} else {
  message("Skipping Static DP demographic-bias model: < ", MIN_DP_OBS, " driving-pressure observations")
}

message("Demographic-bias models fitted (", length(demo_models), " outcomes); rows in regression_results_long")

# =============================================================================
# 4f3. Predicted FVC vs predicted body weight
# =============================================================================
# Does PBW capture predicted lung size? PFVC regressed on PBW + demographics on
# the BROAD cohort (all eligible patients with height/age/sex/race/PFVC, not just
# the ventilated cross-sectional cohort). PFVC is a deterministic function of
# height, age, sex and race, and PBW of height and sex, so the residual is the
# misfit of a linear approximation, not noise: standard errors, intervals and
# p-values have no sampling meaning here and are left empty. The read is the size
# of the age, sex and race coefficients (litres of PFVC at a fixed PBW) and the R^2,
# which is in the rows' note column.
broad_pfvc <- read_parquet(file.path(output_dir, "analysis_broad_pfvc.parquet")) %>%
  mutate(
    sex_category  = factor(sex_category,  levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
    age10 = age_at_admission / 10
  )

pfvc_vs_pbw_model <- lm(pfvc ~ pbw + age10 + sex_category + race_category, data = broad_pfvc)
pfvc_vs_pbw_formula <- "pfvc ~ pbw + age10 + sex_category + race_category"
message("PFVC-vs-PBW model fitted (N = ", nrow(broad_pfvc), "); rows in regression_results_long")

# =============================================================================
# 4g. Unified long-format regression results table
# =============================================================================
# One row per (model, term). Every estimate is reported on its natural scale —
# OR for logistic regression, HR for the Cox model, Beta for linear regression —
# alongside its 95% CI, p-value, estimate type, variable name, and the full model
# specification. The coordinator's pooling script (not distributed with this code)
# stacks this table across sites for the cross-site forest plots, so the column
# schema must stay stable across sites. outcome_scale marks the rows whose outcome
# was z-scored within site (4f2): their Betas are in SD units of the outcome.

extract_model_results <- function(model, estimate_type, analysis,
                                  model_spec, model_family, formula_str,
                                  adjustment = "adjusted",
                                  outcome_scale = "natural units") {
  # OR / HR are reported on the exponentiated (ratio) scale; Beta is the raw
  # linear coefficient. broom::tidy applies the matching transform to estimate
  # and CI together so they stay internally consistent. `adjustment` flags whether
  # the model is demographic-adjusted (age/sex/race) or unadjusted; both are
  # reported, the adjusted as primary.
  exponentiate <- estimate_type %in% c("OR", "HR")
  broom::tidy(model, conf.int = TRUE, exponentiate = exponentiate) %>%
    transmute(
      term,
      estimate,
      conf_low      = conf.low,
      conf_high     = conf.high,
      std_error     = std.error,
      statistic,
      p_value       = p.value,
      estimate_type = estimate_type,
      analysis      = analysis,
      model_spec    = model_spec,
      model_family  = model_family,
      adjustment    = adjustment,
      formula       = formula_str,
      n_obs         = model_patients(model),   # patients, for every family
      n_events      = model_events(model),     # deaths, liberations (Fine-Gray); NA for linear
      outcome_scale = outcome_scale
    )
}

results_long <- list()

# Mortality (logistic regression -> odds ratios)
if (!is.null(mortality_models)) {
  results_long <- c(results_long, imap(mortality_models, ~ {
    formula_str <- paste("deceased ~", exposure_specs[[.y]], "+",
                         model_covariates(exposure_specs[[.y]]))
    extract_model_results(.x, "OR", "Mortality",
                          exposure_labels[[.y]], "logistic", formula_str)
  }))
}

# Continuous outcomes (linear regression -> beta coefficients)
for (outcome_name in names(continuous_outcomes)) {
  outcome_var   <- continuous_outcomes[[outcome_name]]$var
  outcome_label <- continuous_outcomes[[outcome_name]]$label
  results_long <- c(results_long, imap(continuous_models[[outcome_name]], ~ {
    formula_str <- paste(outcome_var, "~", exposure_specs[[.y]], "+",
                         model_covariates(outcome_var))
    extract_model_results(.x, "Beta", outcome_label,
                          exposure_labels[[.y]], "linear", formula_str)
  }))
}

# 28-day VFDs (competing-risks Fine-Gray -> subdistribution hazard ratios for
# extubation; estimate_type "HR" so the cross-cohort forest plots it on the log
# scale, like the other ratio outcomes).
results_long <- c(results_long, imap(vfd_cr_models, ~ {
  formula_str <- paste("finegray(Surv(vfd_time, vfd_status) [extubation vs death]) ~",
                       exposure_specs[[.y]], "+", vfd_cr_covariates)
  extract_model_results(.x, "HR", "28-day VFDs",
                        exposure_labels[[.y]], "finegray", formula_str)
}))

# 60-day death hazard (Cox proportional hazards -> hazard ratios). Three exposure
# specs (PBW/PFVC, PFVC, log PFVC, each with VT/PBW), mirroring the mortality models; HR > 1 = higher death hazard (worse), the
# same direction as the mortality OR.
if (exists("cox_model")) {
  cox_covars  <- "vtpbw + age10 + sex_category + race_category + sf10 + sofa_total"
  results_long <- c(results_long, list(
    # Same vtpbw + pbwpfvc exposure spec as the mortality model — label it with
    # the shared convention so it collapses into one column cross-cohort.
    extract_model_results(cox_model, "HR", "Survival", "VT/PBW + PBW/PFVC", "cox",
                          paste("Surv(surv_time, event) ~ pbwpfvc +", cox_covars)),
    extract_model_results(cox_model_pfvc, "HR", "Survival", "VT/PBW + PFVC", "cox",
                          paste("Surv(surv_time, event) ~ pfvc +", cox_covars)),
    extract_model_results(cox_model_logpfvc, "HR", "Survival", "VT/PBW + log PFVC", "cox",
                          paste("Surv(surv_time, event) ~ log_pfvc +", cox_covars))
  ))
}

# --- Demographic-UNADJUSTED variants of the exposure->outcome models -----------
# Drop age/sex/race; keep illness severity (SOFA + SF) and BMI for DP-derived.
# Reported alongside the adjusted estimates (adjustment = "unadjusted"); the
# adjusted are primary, and the pair shows what the demographics change. The
# rendered tables, E-values, and evidence ratios above use the adjusted models.
make_formula_unadj <- function(lhs, exposure, determinant) {
  cov <- model_covariates(determinant, adjusted = FALSE)
  rhs <- if (cov == "") exposure else paste(exposure, "+", cov)
  as.formula(paste(lhs, "~", rhs))
}
fit_vfd_finegray_unadj <- function(exposure_spec) fit_vfd_finegray(exposure_spec, covariates_unadj)

if (!is.null(mortality_models)) {
  results_long <- c(results_long, imap(exposure_specs, ~ {
    m <- glm(make_formula_unadj("deceased", .x, .x), data = cross_sectional,
             family = binomial)
    fstr <- paste("deceased ~", .x, "+", model_covariates(.x, adjusted = FALSE))
    extract_model_results(m, "OR", "Mortality", exposure_labels[[.y]], "logistic",
                          fstr, adjustment = "unadjusted")
  }))
}

for (outcome_name in names(continuous_outcomes)) {
  outcome_var   <- continuous_outcomes[[outcome_name]]$var
  outcome_label <- continuous_outcomes[[outcome_name]]$label
  results_long <- c(results_long, imap(exposure_specs, ~ {
    m <- lm(make_formula_unadj(outcome_var, .x, outcome_var), data = cross_sectional)
    fstr <- paste(outcome_var, "~", .x, "+", model_covariates(outcome_var, adjusted = FALSE))
    extract_model_results(m, "Beta", outcome_label, exposure_labels[[.y]], "linear",
                          fstr, adjustment = "unadjusted")
  }))
}

results_long <- c(results_long, imap(exposure_specs, ~ {
  m <- fit_vfd_finegray_unadj(.x)
  fstr <- paste("finegray(Surv(vfd_time, vfd_status) [extubation vs death]) ~",
                .x, "+", covariates_unadj)
  extract_model_results(m, "HR", "28-day VFDs", exposure_labels[[.y]], "finegray",
                        fstr, adjustment = "unadjusted")
}))

if (exists("cox_model")) {
  cox_covars_unadj <- "vtpbw + sf10 + sofa_total"   # demographics dropped
  cox_unadj      <- coxph(Surv(surv_time, event) ~ pbwpfvc + vtpbw + sf10 + sofa_total,
                          data = surv_data)
  cox_pfvc_unadj <- coxph(Surv(surv_time, event) ~ pfvc + vtpbw + sf10 + sofa_total,
                          data = surv_data)
  cox_logpfvc_unadj <- coxph(Surv(surv_time, event) ~ log_pfvc + vtpbw + sf10 + sofa_total,
                             data = surv_data)
  results_long <- c(results_long, list(
    extract_model_results(cox_logpfvc_unadj, "HR", "Survival", "VT/PBW + log PFVC", "cox",
                          paste("Surv(surv_time, event) ~ log_pfvc +", cox_covars_unadj),
                          adjustment = "unadjusted"),
    extract_model_results(cox_unadj, "HR", "Survival", "VT/PBW + PBW/PFVC", "cox",
                          paste("Surv(surv_time, event) ~ pbwpfvc +", cox_covars_unadj),
                          adjustment = "unadjusted"),
    extract_model_results(cox_pfvc_unadj, "HR", "Survival", "VT/PBW + PFVC", "cox",
                          paste("Surv(surv_time, event) ~ pfvc +", cox_covars_unadj),
                          adjustment = "unadjusted")
  ))
}

# Demographic-bias family (one analysis name per outcome, single "Demographics" spec)
results_long <- c(results_long, imap(demo_models, ~ {
  extract_model_results(.x$model, if (.x$type == "OR") "OR" else "Beta",
                        paste0("Demo bias: ", .y), "Demographics",
                        .x$family, .x$formula,
                        outcome_scale = .x$outcome_scale)
}))

# Predicted FVC vs PBW (broad cohort)
# (4f3: coefficient sizes and R^2 only; the inferential columns are emptied)
results_long <- c(results_long, list(
  extract_model_results(pfvc_vs_pbw_model, "Beta", "PFVC vs PBW",
                        "PFVC ~ PBW", "linear", pfvc_vs_pbw_formula) %>%
    mutate(across(c(conf_low, conf_high, std_error, statistic, p_value), ~ NA_real_),
           note = sprintf(paste("PFVC is a deterministic function of the regressors: coefficient sizes",
                                "and R^2 only (R^2 = %.4f); no standard errors or p-values"),
                          summary(pfvc_vs_pbw_model)$r.squared))
))

# Drop intercepts (not a reportable effect) and stamp the site name so the
# table is self-describing once pooled across sites.
regression_results_long <- bind_rows(results_long) %>%
  filter(term != "(Intercept)") %>%
  mutate(site = site_name, .before = 1)

write_csv(regression_results_long,
          file.path(final_dir, paste0("regression_results_long_", site_name, ".csv")))

message("Unified regression results table: ", nrow(regression_results_long),
        " rows across ", n_distinct(regression_results_long$analysis), " analyses; saved to ",
        file.path(final_dir, paste0("regression_results_long_", site_name, ".csv")))

# =============================================================================
# 4g2. Residual confounding: E-values + age functional-form (spline) sensitivity
# =============================================================================
# Because PBW and PFVC are DETERMINISTIC functions of {height, age, sex, race},
# those four variables are the complete parent set of every PBW/PFVC-derived
# exposure, and the parents of treatment are a sufficient backdoor adjustment set.
# The confounder space is therefore closed (not open-ended), and residual
# confounding can only enter through the gap between that full parent set and what
# the models actually adjust for. Height is NOT in that gap: in the DAG it reaches
# mortality only through predicted lung size (height -> lung size -> strain ->
# mortality), so it lies on the causal pathway / is the identifying variation in
# the exposure, not a backdoor — conditioning on it is over-adjustment. The single
# enumerable residual confounder is therefore the FUNCTIONAL FORM of age (a genuine
# confounder, entered linearly in the main models). Anything beyond that would have
# to be a truly unmeasured common cause acting outside this deterministic structure
# — exactly what the E-value bounds.
#
# This section reports, for each exposure term in EXPOSURE_TERMS in the mortality
# OR, 60-day survival HR and 28-day liberation Fine-Gray SHR models:
#   1. the per-1-SD estimate under LINEAR age (the main-model specification),
#   2. the per-1-SD estimate under SPLINE age (race + ns(age, 4) + sex + SOFA + SF):
#      the direct bound on nonlinear-age confounding — an estimate unchanged under
#      flexible age is robust to it; one that collapses was carrying age, and
#   3. the E-value (point and CI-limit) for the linear-age estimate: the minimum
#      association, on the risk-ratio scale, an unmeasured confounder would need
#      with BOTH the exposure and the outcome to explain the estimate (or, for the
#      CI E-value, to move the CI to include the null).
#
# Estimates are rescaled to a 1-SD increase in the exposure (E-values on the raw
# per-unit scale would sit artificially close to 1). Mortality and liberation are
# common outcomes, so each OR/HR is converted to an approximate risk ratio before
# the E-value (rare = FALSE; VanderWeele & Ding, Ann Intern Med 2017); the
# Fine-Gray subdistribution HR uses the same common-outcome HR conversion.

# vt_excess_ml is left out: it is the difference parameterization of the same
# PBW-vs-PFVC disagreement that pbwpfvc carries (see exposure_specs), so the
# "VT/PBW + VT excess (mL)" models are refitted with spline age but contribute only
# their VT/PBW row to this table.
EXPOSURE_TERMS <- c("vtpfvc", "vtpbw", "pfvc", "log_pfvc", "pbwpfvc")
exposure_term_labels_ev <- c(vtpfvc = "VT/PFVC", vtpbw = "VT/PBW",
                             pfvc = "PFVC", log_pfvc = "log PFVC", pbwpfvc = "PBW/PFVC")

# Spline-age covariate set (mirrors the linear-age set, age10 -> ns(age, 4)).
# Reused for the mortality and Fine-Gray refits below; the Cox refit is spelled
# out separately because it carries its own exposure terms.
spline_covars <- "race_category + ns(age_at_admission, 4) + sex_category + sofa_total + sf10"

# Per-SD-rescaled ratio estimate + 95% CI from a fitted model. beta and its SE
# are the log-OR / log-HR (link scale); multiplying by the exposure SD gives the
# per-SD log-ratio, exponentiated back to the ratio scale.
persd_ratio <- function(model, term, exposure_sd) {
  beta <- coef(model)[[term]]
  se   <- sqrt(diag(vcov(model)))[[term]]
  list(
    estimate  = exp(beta * exposure_sd),
    conf_low  = exp((beta - 1.96 * se) * exposure_sd),
    conf_high = exp((beta + 1.96 * se) * exposure_sd)
  )
}

# Per-SD estimate + CI for every exposure term in a model (no E-value). The SD is
# the marginal, patient-level SD of the exposure in sd_data, so the linear- and
# spline-age fits are compared on an identical contrast.
persd_terms <- function(model, sd_data, terms = EXPOSURE_TERMS) {
  present <- intersect(terms, names(coef(model)))
  map_dfr(present, function(tm) {
    exposure_sd <- sd(sd_data[[tm]], na.rm = TRUE)
    pr <- persd_ratio(model, tm, exposure_sd)
    tibble(term = tm, exposure_sd = exposure_sd,
           estimate = pr$estimate, conf_low = pr$conf_low, conf_high = pr$conf_high)
  })
}

# Point and CI E-values for one per-SD ratio estimate. type is "OR" or "HR";
# both use rare = FALSE (common outcomes). The CI E-value is the non-NA bound
# the package returns (the limit nearest the null; 1 if the CI crosses it).
evalue_for <- function(est, lo, hi, type) {
  ev_mat <- switch(type,
    OR = EValue::evalues.OR(est, lo, hi, rare = FALSE),
    HR = EValue::evalues.HR(est, lo, hi, rare = FALSE)
  )
  ev    <- ev_mat["E-values", ]
  ci_ev <- ev[c("lower", "upper")]
  ci_ev <- ci_ev[!is.na(ci_ev)]
  list(evalue_point = unname(ev[["point"]]),
       evalue_ci    = if (length(ci_ev)) unname(ci_ev[1]) else NA_real_)
}

# One model's exposure rows: per-SD estimate under linear age (model_lin) and
# under spline age (model_spl) side by side, plus the E-value for the linear-age
# estimate. Both fits use the same sd_data so the per-SD contrast is identical.
residual_conf_rows <- function(model_lin, model_spl, model_spec, analysis, type, sd_data) {
  lin <- persd_terms(model_lin, sd_data)
  spl <- persd_terms(model_spl, sd_data)
  if (nrow(lin) == 0) return(tibble())
  ev  <- pmap_dfr(list(lin$estimate, lin$conf_low, lin$conf_high),
                  function(e, l, h) as_tibble(evalue_for(e, l, h, type)))
  spl_i <- match(lin$term, spl$term)
  tibble(
    analysis      = analysis,
    model_spec    = model_spec,
    term          = lin$term,
    term_label    = unname(exposure_term_labels_ev[lin$term]),
    estimate_type = type,
    exposure_sd   = lin$exposure_sd,
    est_linage    = lin$estimate,
    lo_linage     = lin$conf_low,
    hi_linage     = lin$conf_high,
    est_splineage = spl$estimate[spl_i],
    lo_splineage  = spl$conf_low[spl_i],
    hi_splineage  = spl$conf_high[spl_i],
    evalue_point  = ev$evalue_point,
    evalue_ci     = ev$evalue_ci
  )
}

# --- Spline-age refits of each ratio model (age10 -> ns(age_at_admission, 4)) ---
# The exposure specs are dosing ratios, never DP-derived, so no BMI is added.

mortality_models_spline <- if (!is.null(mortality_models)) {
  map(exposure_specs, ~ glm(
    as.formula(paste("deceased ~", .x, "+", spline_covars)),
    data = cross_sectional, family = binomial))
} else NULL

cox_model_spline <- if (exists("cox_model")) {
  coxph(Surv(surv_time, event) ~ pbwpfvc + vtpbw + ns(age_at_admission, 4) +
          sex_category + race_category + sf10 + sofa_total, data = surv_data)
} else NULL
cox_model_pfvc_spline <- if (exists("cox_model_pfvc")) {
  coxph(Surv(surv_time, event) ~ pfvc + vtpbw + ns(age_at_admission, 4) +
          sex_category + race_category + sf10 + sofa_total, data = surv_data)
} else NULL
cox_model_logpfvc_spline <- if (exists("cox_model_logpfvc")) {
  coxph(Surv(surv_time, event) ~ log_pfvc + vtpbw + ns(age_at_admission, 4) +
          sex_category + race_category + sf10 + sofa_total, data = surv_data)
} else NULL

# Fine-Gray refit with spline age (fit_vfd_finegray of section 4d2).
fit_vfd_finegray_spline <- function(exposure_spec) fit_vfd_finegray(exposure_spec, spline_covars)
vfd_cr_models_spline <- map(exposure_specs, fit_vfd_finegray_spline)

# --- Assemble residual-confounding rows across all ratio models ----------------
residual_list <- list()

# In-hospital mortality (logistic, OR) — one model per exposure specification.
if (!is.null(mortality_models)) {
  residual_list <- c(residual_list, imap(mortality_models,
    ~ residual_conf_rows(.x, mortality_models_spline[[.y]], exposure_labels[[.y]],
                         "Mortality (in-hospital)", "OR", cross_sectional)))
}

# 60-day death hazard (Cox, HR) — PBW/PFVC, PFVC and log PFVC exposure specs.
if (!is.null(cox_model_spline)) {
  residual_list <- c(residual_list, list(
    residual_conf_rows(cox_model, cox_model_spline, "VT/PBW + PBW/PFVC",
                       "Survival (60-day)", "HR", surv_data),
    residual_conf_rows(cox_model_pfvc, cox_model_pfvc_spline, "VT/PBW + PFVC",
                       "Survival (60-day)", "HR", surv_data),
    residual_conf_rows(cox_model_logpfvc, cox_model_logpfvc_spline, "VT/PBW + log PFVC",
                       "Survival (60-day)", "HR", surv_data)))
}

# 28-day ventilator liberation (Fine-Gray subdistribution HR). The finegray()
# expansion duplicates rows per subject, so the per-SD contrast uses the
# patient-level SD from the cross-sectional cohort (rows with a valid VFD time),
# not the expanded risk set.
vfd_sd_data <- cross_sectional %>% filter(!is.na(vfd_time), vfd_time > 0)
residual_list <- c(residual_list, imap(vfd_cr_models,
  ~ residual_conf_rows(.x, vfd_cr_models_spline[[.y]], exposure_labels[[.y]],
                       "28-day VFDs (liberation)", "HR", vfd_sd_data)))

residual_confounding <- bind_rows(residual_list) %>%
  mutate(analysis = factor(analysis, levels = c(
    "Survival (60-day)", "Mortality (in-hospital)", "28-day VFDs (liberation)"))) %>%
  arrange(analysis) %>%               # stable: preserves model & term order within analysis
  mutate(analysis = as.character(analysis)) %>%
  mutate(site = site_name, .before = 1)

# evalues_<site>: per-SD estimates under linear and spline age, with E-values for
# the linear-age estimate.
write_csv(residual_confounding, file.path(final_dir, paste0("evalues_", site_name, ".csv")))
message("Residual-confounding table written (", nrow(residual_confounding),
        " ratio-scale exposure estimates across ",
        n_distinct(residual_confounding$analysis),
        " analyses; linear/spline age + E-values)")

# =============================================================================
# 4i. Inclusion CONSORT diagram (site QC)
# =============================================================================
# CONSORT flow from the 7-step attrition log written in script 03.
attrition <- read_csv(file.path(final_dir, paste0("attrition_log_", site_name, ".csv")),
                      show_col_types = FALSE)
consort_fig <- render_consort(attrition, title = paste0("Cohort inclusion - ", site_name))
ggsave(file.path(final_dir, paste0("consort_diagram_", site_name, ".pdf")),
       consort_fig, width = 9, height = 11)
message("CONSORT diagram saved")

# =============================================================================
# 4j. Negative-control cohorts: PBW/PFVC, PFVC, height (and dose) vs mortality
# =============================================================================
# The analytic cohort is hypoxemic AND ventilated with PBW-dosed tidal volumes. Two
# non-hypoxemic cohorts (script 01 / 03k) separate the pathways:
#   * Ventilated, non-hypoxemic (dosed, uninjured lung): still receives PBW-dosed
#     tidal volumes (negative_control_counts reports its delivered VT/PBW). A
#     NEGATIVE control for hypoxemia-specific mechanisms and a POSITIVE control for
#     the dosing pathway: PBW/PFVC harm should persist here. The strict definition
#     (never hypoxemic during ventilation) leaves few patients per site, so this
#     cohort is read only when pooled.
#   * Not ventilated (no tidal volume): the true no-dose control. A ventilatory
#     pathway predicts attenuation of PBW/PFVC, PFVC and height here.
# The three cohorts share ONE adjustment set (age, sex, race) because SOFA and the
# SF ratio exist only for the analytic cohort; the analytic cohort's fully adjusted
# estimates live in 4c/4f. Each exposure is standardized by the ANALYTIC cohort's
# SD so estimates are on one scale across cohorts. VT/PBW and VT/PFVC are fit in
# the two ventilated cohorts (the ventilated control's delivered dose comes from
# 03k). Every model is fit with linear age AND with ns(age, 4): what remains of the
# ratio after linear age, sex and race is height plus the convex part of the age
# curve, so a residual ratio effect where height is null is read against the spline.
# Outcomes: in-hospital death (logistic) and 60-day all-cause death (Cox), each
# cohort timed from its own origin (script 03: the index for the analytic cohort,
# the first IMV record for the ventilated control, ICU admission for the
# non-ventilated one). A cell with fewer than NC_MIN_EVENTS deaths is not fitted.
# The delivered dose is defined differently in the two ventilated cohorts: the
# analytic cohort's VT/PBW and VT/PFVC are the values at its index timepoint, the
# ventilated control's the per-patient median over all its volume-targeted IMV
# timepoints (03k). The dose_definition column of the tables says which.
nc_file <- file.path(output_dir, "analysis_negative_control.parquet")
nc_data <- read_parquet(nc_file) %>%
  select(hospitalization_id, nc_cohort, age_at_admission, sex_category, race_category,
         height_cm, pbw, pfvc, pfvc_age25, pbwpfvc, vt_excess_ml, vtpbw, vtpfvc,
         deceased, mortality_event_60, surv_time)
nc_cohort_levels <- c("Hypoxemic, ventilated (analytic)",
                      "Ventilated, non-hypoxemic (dosed, uninjured lung)",
                      "Not ventilated (no tidal volume)")
nc_frames <- bind_rows(
  cross_sectional %>%
    transmute(hospitalization_id, nc_cohort = nc_cohort_levels[1],
              age_at_admission, sex_category, race_category, height_cm, pbw, pfvc, pfvc_age25,
              pbwpfvc, vt_excess_ml, vtpbw, vtpfvc, deceased, mortality_event_60, surv_time),
  nc_data) %>%
  mutate(age10 = age_at_admission / 10,
         sex_category  = factor(sex_category,  levels = c("Male", "Female")),
         race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")))
nc_sd <- cross_sectional %>%
  summarise(across(c(pbwpfvc, pfvc, height_cm, vt_excess_ml, vtpbw, vtpfvc), ~ sd(.x, na.rm = TRUE)))
nc_exposures <- c(pbwpfvc = "PBW/PFVC", pfvc = "PFVC", height_cm = "Height",
                  vt_excess_ml = "VT excess (mL)", vtpbw = "VT/PBW", vtpfvc = "VT/PFVC")
nc_age_forms <- c(linear = "age10", spline = "splines::ns(age_at_admission, 4)")
NC_DOSE_DEFINITION <- setNames(c("value at the index timepoint",
                                 "per-patient median over all volume-targeted IMV timepoints",
                                 NA_character_), nc_cohort_levels)

nc_fit_one <- function(df, expo, cohort_lab) {
  d <- df %>% filter(nc_cohort == cohort_lab) %>%
    mutate(z = .data[[expo]] / nc_sd[[expo]]) %>%
    filter(is.finite(z))
  n_death <- sum(d$deceased == 1, na.rm = TRUE); n_death60 <- sum(d$mortality_event_60 == 1, na.rm = TRUE)
  out <- tibble()
  for (af in names(nc_age_forms)) {
    rhs <- paste("z +", nc_age_forms[[af]], "+ sex_category + race_category")
    if (nrow(d) >= NC_MIN_PATIENTS && n_death >= NC_MIN_EVENTS) {
      m <- glm(as.formula(paste("deceased ~", rhs)), data = d, family = binomial)
      cf <- summary(m)$coefficients["z", ]
      out <- bind_rows(out, tibble(age_form = af, outcome = "In-hospital mortality", estimate_type = "OR",
        estimate = exp(cf["Estimate"]), conf_low = exp(cf["Estimate"] - 1.96 * cf["Std. Error"]),
        conf_high = exp(cf["Estimate"] + 1.96 * cf["Std. Error"]), std_error = cf["Std. Error"],
        p_value = cf["Pr(>|z|)"], n = nrow(d), events = n_death))
    }
    if (nrow(d) >= NC_MIN_PATIENTS && n_death60 >= NC_MIN_EVENTS) {
      d60 <- d %>% filter(!is.na(surv_time), surv_time > 0)
      m <- survival::coxph(as.formula(paste("survival::Surv(surv_time, mortality_event_60) ~", rhs)), data = d60)
      cf <- summary(m)$coefficients["z", ]
      out <- bind_rows(out, tibble(age_form = af, outcome = "60-day mortality", estimate_type = "HR",
        estimate = cf["exp(coef)"], conf_low = exp(cf["coef"] - 1.96 * cf["se(coef)"]),
        conf_high = exp(cf["coef"] + 1.96 * cf["se(coef)"]), std_error = cf["se(coef)"],
        p_value = cf["Pr(>|z|)"], n = nrow(d60), events = n_death60))
    }
  }
  if (nrow(out)) out %>% mutate(cohort = cohort_lab, exposure = nc_exposures[[expo]],
                                scale = "per analytic-cohort SD", .before = 1) else out
}
nc_results <- map_dfr(nc_cohort_levels, function(cl)
  map_dfr(names(nc_exposures), function(e) nc_fit_one(nc_frames, e, cl)))
if (nrow(nc_results)) nc_results <- nc_results %>%
  mutate(adjustment = "sex + race + age (linear or ns4)",
         dose_definition = if_else(exposure %in% c("VT/PBW", "VT/PFVC"), unname(NC_DOSE_DEFINITION[cohort]), NA_character_),
         site = site_name)
nc_counts <- nc_frames %>% group_by(cohort = nc_cohort) %>%
  summarise(n = n(), deaths_inhosp = sum(deceased == 1, na.rm = TRUE),
            deaths_60d = sum(mortality_event_60 == 1, na.rm = TRUE),
            median_age = median(age_at_admission), pct_female = mean(sex_category == "Female"),
            median_pbwpfvc = median(pbwpfvc, na.rm = TRUE),
            # delivered dose (ventilated cohorts only): the dosing-pathway read
            n_with_vt = sum(!is.na(vtpbw)),
            median_vtpbw = median(vtpbw, na.rm = TRUE), q25_vtpbw = quantile(vtpbw, .25, na.rm = TRUE),
            q75_vtpbw = quantile(vtpbw, .75, na.rm = TRUE),
            pct_vtpbw_over_8 = mean(vtpbw > VTPBW_BAND_MAX, na.rm = TRUE),
            median_vtpfvc = median(vtpfvc, na.rm = TRUE),
            pct_vtpfvc_over_11 = mean(vtpfvc > VTPFVC_ARMA_P75, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(dose_definition = unname(NC_DOSE_DEFINITION[cohort]),
         cohort = factor(cohort, nc_cohort_levels)) %>% arrange(cohort)
write_csv(nc_results, file.path(final_dir, paste0("negative_control_", site_name, ".csv")))
write_csv(nc_counts,  file.path(final_dir, paste0("negative_control_counts_", site_name, ".csv")))

# --- Identifying variation: what is left of each exposure after the adjusters -------
# The ratio is almost entirely a function of the demographics (the share is reported
# in residual_variance_share), so after age, sex and race its coefficient is
# estimated from a sliver of residual variation; PFVC keeps far more. The share of
# each exposure's variance that survives the adjustment set explains the CI widths
# in one number, and the share of THAT residual explained by height says what the
# surviving variation is. Height enters as a spline within sex: GLI FVC is a power
# law in height but Devine PBW is a line with an intercept, so at fixed age, sex and
# race log PBW/PFVC is curved in height and runs opposite ways by sex (a hump near
# 167 cm in men, rising to ~185 cm in women). A linear height term finds almost none
# of that; PFVC's height channel is near log-linear and either form captures it.
# Reported per cohort, under linear and spline age, with and without VT/PBW.
nc_idvar <- map_dfr(nc_cohort_levels, function(cl) {
  d <- nc_frames %>% filter(nc_cohort == cl)
  map_dfr(names(nc_exposures), function(e) {
    z <- d[[e]]; ok <- is.finite(z)
    if (sum(ok) < NC_MIN_PATIENTS) return(tibble())
    map_dfr(names(nc_age_forms), function(af) {
      map_dfr(c("demographics", "demographics + VT/PBW"), function(adj) {
        if (adj == "demographics + VT/PBW" && (e == "vtpbw" || all(is.na(d$vtpbw)))) return(tibble())
        rhs <- paste(nc_age_forms[[af]], "+ sex_category + race_category",
                     if (adj == "demographics + VT/PBW") "+ vtpbw" else "")
        dd <- d[ok, ]; if (adj == "demographics + VT/PBW") dd <- dd %>% filter(is.finite(vtpbw))
        if (nrow(dd) < NC_MIN_PATIENTS) return(tibble())
        fit <- lm(as.formula(paste(e, "~", rhs)), data = dd)
        tibble(cohort = cl, exposure = nc_exposures[[e]], age_form = af, adjustment = adj,
               n = nrow(dd), r2_on_adjusters = summary(fit)$r.squared,
               residual_variance_share = 1 - summary(fit)$r.squared,
               residual_sd_in_analytic_sd = sd(resid(fit)) / nc_sd[[e]],
               # how much of what survives the adjusters is height: R2 of the residual on a
               # sex-specific natural spline in height (the ratio's height channel is curved)
               height_share_of_residual = summary(lm(resid(fit) ~ dd$sex_category * ns(dd$height_cm, 3)))$r.squared)
      })
    })
  })
}) %>% mutate(site = site_name)
write_csv(nc_idvar, file.path(final_dir, paste0("negative_control_identifying_variation_", site_name, ".csv")))
cat("--- identifying variation (share of exposure variance left after the adjusters, and how much of that is height; linear age) ---\n")
print(as.data.frame(nc_idvar %>% filter(age_form == "linear") %>%
        transmute(cohort, exposure, adjustment, n, resid_share = round(residual_variance_share, 3),
                  height_share_of_resid = round(height_share_of_residual, 3))), row.names = FALSE)

# --- Formal cohort contrast: does the exposure's effect differ where tidal volume is set? --
# One model over the stacked cohorts with cohort-specific effects of every adjuster
# (equivalent to fitting each cohort separately) and an exposure x cohort interaction;
# the LRT against the no-interaction model tests heterogeneity, and the pairwise
# contrasts (each control minus the reference, log scale) say where it lies. The
# reference is the analytic cohort. When the analytic cohort falls below
# NC_MIN_PATIENTS or NC_MIN_EVENTS for an outcome, the first remaining cohort
# becomes the reference: the script says so, and reference_cohort names it in every
# row. Cox models stratify the baseline hazard by cohort. Linear age. The contrasts
# are poolable across sites (random effects on the log difference); the p-values by
# Fisher.
nc_interaction <- map_dfr(names(nc_exposures), function(e) {
  d <- nc_frames %>% mutate(z = .data[[e]] / nc_sd[[e]]) %>% filter(is.finite(z))
  keep <- d %>% count(nc_cohort) %>% filter(n >= NC_MIN_PATIENTS) %>% pull(nc_cohort)
  d <- d %>% filter(nc_cohort %in% keep) %>%
    mutate(cohort = factor(nc_cohort, levels = intersect(nc_cohort_levels, keep)))
  if (n_distinct(d$cohort) < 2) return(tibble())
  one <- function(outcome_lab) {
    if (outcome_lab == "In-hospital mortality") {
      ev_ok <- d %>% group_by(cohort) %>% summarise(ev = sum(deceased == 1), .groups = "drop")
      dd <- d %>% filter(cohort %in% ev_ok$cohort[ev_ok$ev >= NC_MIN_EVENTS]) %>% droplevels()
      if (n_distinct(dd$cohort) < 2) return(tibble())
      f0 <- glm(deceased ~ z + cohort * (age10 + sex_category + race_category), data = dd, family = binomial)
      f1 <- glm(deceased ~ z * cohort + cohort * (age10 + sex_category + race_category), data = dd, family = binomial)
      V <- vcov(f1); b <- coef(f1); est_type <- "OR"
    } else {
      dd <- d %>% filter(!is.na(surv_time), surv_time > 0)
      ev_ok <- dd %>% group_by(cohort) %>% summarise(ev = sum(mortality_event_60 == 1), .groups = "drop")
      dd <- dd %>% filter(cohort %in% ev_ok$cohort[ev_ok$ev >= NC_MIN_EVENTS]) %>% droplevels()
      if (n_distinct(dd$cohort) < 2) return(tibble())
      f0 <- survival::coxph(survival::Surv(surv_time, mortality_event_60) ~ z + cohort:(age10 + sex_category + race_category) +
                              survival::strata(cohort), data = dd)
      f1 <- survival::coxph(survival::Surv(surv_time, mortality_event_60) ~ z * cohort + cohort:(age10 + sex_category + race_category) +
                              survival::strata(cohort), data = dd)
      V <- vcov(f1); b <- coef(f1); est_type <- "HR"
    }
    # the reference is the first cohort left after the patient and death gates
    ref <- levels(dd$cohort)[1]
    if (ref != nc_cohort_levels[1])
      message("*** 4j contrast, ", nc_exposures[[e]], ", ", outcome_lab, ": the analytic cohort is below ",
              NC_MIN_PATIENTS, " patients or ", NC_MIN_EVENTS, " deaths, so the reference is '", ref,
              "', not the analytic cohort ***")
    lrt <- 2 * (as.numeric(logLik(f1)) - as.numeric(logLik(f0)))
    df  <- n_distinct(dd$cohort) - 1
    int_terms <- grep("^z:cohort", names(b), value = TRUE)
    contrasts <- map_dfr(int_terms, function(tm) {
      cl <- sub("^z:cohort", "", tm)
      tibble(contrast = paste0(cl, " minus ", ref), log_diff = unname(b[tm]),
             se = sqrt(V[tm, tm]), ratio_of_effects = exp(unname(b[tm])),
             lo = exp(unname(b[tm]) - 1.96 * sqrt(V[tm, tm])), hi = exp(unname(b[tm]) + 1.96 * sqrt(V[tm, tm])))
    })
    bind_rows(tibble(contrast = "LRT: exposure x cohort", log_diff = NA_real_, se = NA_real_,
                     ratio_of_effects = NA_real_, lo = NA_real_, hi = NA_real_),
              contrasts) %>%
      mutate(outcome = outcome_lab, estimate_type = est_type, lrt_chi2 = lrt, lrt_df = df,
             lrt_p = pchisq(lrt, df, lower.tail = FALSE), reference_cohort = ref,
             cohorts = paste(levels(dd$cohort), collapse = " | "), n = nrow(dd), .before = 1)
  }
  res <- bind_rows(one("In-hospital mortality"), one("60-day mortality"))
  if (nrow(res) == 0) return(tibble())   # mutate() on an empty tibble would manufacture a row
  res %>% mutate(exposure = nc_exposures[[e]], .before = 1)
})
if (nrow(nc_interaction)) nc_interaction <- nc_interaction %>%
  mutate(scale = paste("effect per analytic-cohort SD; ratio_of_effects = cohort / reference_cohort",
                       "(< 1 = attenuated relative to the reference, the analytic cohort unless reference_cohort says otherwise)"),
         site = site_name)
write_csv(nc_interaction, file.path(final_dir, paste0("negative_control_interaction_", site_name, ".csv")))
if (nrow(nc_interaction)) {
  cat("--- exposure x cohort heterogeneity (LRT) and contrasts vs the analytic cohort ---\n")
  print(as.data.frame(nc_interaction %>% transmute(exposure, outcome, contrast,
          ratio = ifelse(is.na(ratio_of_effects), NA, sprintf("%.2f [%.2f, %.2f]", ratio_of_effects, lo, hi)),
          lrt_p = signif(lrt_p, 2))), row.names = FALSE)
}
message("Negative-control models: ", nrow(nc_results), " estimates across ",
        n_distinct(nc_results$cohort), " cohort(s)")
print(as.data.frame(nc_counts %>% transmute(cohort, n, deaths_inhosp, deaths_60d, n_with_vt,
        vtpbw = ifelse(is.na(median_vtpbw), NA, sprintf("%.1f [%.1f-%.1f]", median_vtpbw, q25_vtpbw, q75_vtpbw)),
        pct_over_8 = round(100 * pct_vtpbw_over_8), vtpfvc = round(median_vtpfvc, 1))), row.names = FALSE)
# no cohort reached NC_MIN_EVENTS deaths (a small site, or the synthetic subset): nothing to print
if (nrow(nc_results)) print(as.data.frame(nc_results %>% filter(age_form == "linear") %>%
        transmute(cohort, exposure, outcome, n, events,
                  est = sprintf("%.2f [%.2f, %.2f]", estimate, conf_low, conf_high))), row.names = FALSE)
# =============================================================================
# 4k. The saturated log model: log VT, log PBW, log PFVC (demographic-unadjusted)
# =============================================================================
# In a log-linear model every combination of VT/PBW, PFVC, PBW/PFVC and VT/PFVC lives in
# the span of three columns, log VT, log PBW and log PFVC. The saturated model with all
# three is the reference; each two-term model is that model plus one linear constraint,
# and each "parameterization" is a change of basis with the same likelihood:
#   log PFVC at fixed VT and PBW: same breath, same PBW label, larger predicted lung --
#       lower strain and lower ratio together (the strain-error effect)
#   log PBW  at fixed VT and PFVC: same breath, same lung, higher PBW -- the mL/kg label
#       falls and the ratio rises with no change in strain (label and height, no dose)
#   log VT   at fixed PBW and PFVC: pure dose, confounded by indication within 6-8 mL/kg
# The saturated model is fitted WITHOUT the demographics (SOFA and SF only). With age,
# sex and race in the model, log PBW is a function of height within sex and log PFVC
# nearly so, so the two columns are almost collinear and their separate coefficients
# are not identified; log PBW and log PFVC are never fitted together with the
# demographics. Unadjusted, the log PFVC coefficient at fixed PBW carries the age,
# sex and race content of PFVC, including age's own path to death: read it as an
# association, beside the adjusted two-term models of 4c and 4f.
# The two-term models become tests: "VT/PBW + PFVC" imposes b_VT + b_PBW = 0; the ratio
# model "VT/PBW + PBW/PFVC" imposes b_VT + b_PBW + b_PFVC = 0, which is SCALE INVARIANCE
# (scale breath, body and lung together and nothing changes -- only dimensionless ratios
# carry information), a physiologic hypothesis tested as a one-df Wald test, plus the
# LRT of each constraint. Coefficients per log unit (x1.1 = a 10% increase). (B) is
# the functional-form ladder for the size term (PFVC linear / log / 1/x; ratio linear /
# log), with the demographics and linear or spline age; it holds one size term at a
# time, so log PBW and log PFVC never meet there. In-hospital (logistic) and 60-day
# (Cox, from the index) death; each outcome only where it has at least NC_MIN_EVENTS
# deaths, and only with SATURATED_MIN_PATIENTS complete patients.
sz <- cross_sectional %>%
  filter(vtpbw > 0, pbwpfvc > 0, pbw > 0, pfvc > 0, !is.na(sofa_total), !is.na(sf10)) %>%
  mutate(l_vt = log(tidal_volume_set), l_vtpbw = log(vtpbw), l_ratio = log(pbwpfvc), l_pbw = log(pbw),
         l_pfvc = log(pfvc), l_vtpfvc = log(vtpfvc), inv_pfvc = 1 / pfvc)
message("4k frame (positive VT/PBW, PBW/PFVC, PBW and PFVC; SOFA and SF recorded): ",
        nrow(cross_sectional), " -> ", nrow(sz), " patients")
sz_age <- c(linear = "age10", spline = "splines::ns(age_at_admission, 4)")
sz_base <- "sex_category + race_category + sofa_total + sf10"
sz_saturated_covariates <- "sofa_total + sf10"   # (A): no demographics (see the header)
sz_fit <- function(rhs, outcome) {
  if (outcome == "In-hospital mortality")
    glm(as.formula(paste("deceased ~", rhs)), data = sz, family = binomial)
  else survival::coxph(as.formula(paste("survival::Surv(surv_time, mortality_event_60) ~", rhs)),
                       data = sz %>% filter(!is.na(surv_time), surv_time > 0))
}
# each outcome runs only where it has at least NC_MIN_EVENTS deaths
sz_events <- c("In-hospital mortality" = sum(sz$deceased == 1, na.rm = TRUE),
               "60-day mortality" = sum(sz$mortality_event_60 == 1, na.rm = TRUE))
sz_outcomes <- names(sz_events)[sz_events >= NC_MIN_EVENTS]
sz_ok <- length(sz_outcomes) >= 1 && nrow(sz) >= SATURATED_MIN_PATIENTS

if (sz_ok) {
  # --- (A) the saturated log model + constraint tests ---------------------------------
  sz_wald <- function(b, V, w) { est <- sum(w * b); se <- sqrt(as.numeric(t(w) %*% V %*% w)); c(est = est, se = se, p = 2 * pnorm(-abs(est / se))) }
  size_saturated <- map_dfr(sz_outcomes, function(oc) {
    cov <- sz_saturated_covariates
    f_sat   <- sz_fit(paste("l_vt + l_pbw + l_pfvc +", cov), oc)
    f_pfvc  <- sz_fit(paste("l_vtpbw + l_pfvc +", cov), oc)          # b_VT + b_PBW = 0
    f_ratio <- sz_fit(paste("l_vtpbw + l_ratio +", cov), oc)         # b_VT + b_PBW + b_PFVC = 0 (scale invariance)
    f_strain<- sz_fit(paste("l_vtpfvc +", cov), oc)                  # b_VT + b_PFVC = 0 and b_PBW = 0
    terms <- c("l_vt", "l_pbw", "l_pfvc"); b <- coef(f_sat)[terms]; V <- vcov(f_sat)[terms, terms]
    cf <- summary(f_sat)$coefficients[terms, ]
    est <- if (inherits(f_sat, "coxph")) cf[, "coef"] else cf[, "Estimate"]
    se  <- if (inherits(f_sat, "coxph")) cf[, "se(coef)"] else cf[, "Std. Error"]
    lrt <- function(f0) { x <- 2 * (as.numeric(logLik(f_sat)) - as.numeric(logLik(f0))); c(x, attr(logLik(f_sat), "df") - attr(logLik(f0), "df")) }
    w_scale <- c(1, 1, 1); w_pfvc <- c(1, 1, 0)
    ws <- sz_wald(b, V, w_scale); wp <- sz_wald(b, V, w_pfvc)
    l_p <- lrt(f_pfvc); l_r <- lrt(f_ratio); l_s <- lrt(f_strain)
    bind_rows(
      tibble(row_type = "coefficient", term = terms,
             term_label = c("log VT at fixed PBW and PFVC (dose)", "log PBW at fixed VT and PFVC (label and height, no strain change)",
                            "log PFVC at fixed VT and PBW (strain error)"),
             estimate = exp(est), conf_low = exp(est - 1.96 * se), conf_high = exp(est + 1.96 * se), std_error = se,
             p_value = cf[, ncol(cf)]),
      tibble(row_type = "constraint", term = c("b_VT + b_PBW + b_PFVC = 0", "b_VT + b_PBW = 0"),
             term_label = c("scale invariance (only dimensionless ratios matter; the ratio model)",
                            "VT enters only as VT/PBW (the PFVC model)"),
             estimate = exp(c(ws["est"], wp["est"])), conf_low = exp(c(ws["est"] - 1.96 * ws["se"], wp["est"] - 1.96 * wp["se"])),
             conf_high = exp(c(ws["est"] + 1.96 * ws["se"], wp["est"] + 1.96 * wp["se"])), std_error = c(ws["se"], wp["se"]),
             p_value = c(ws["p"], wp["p"]),
             lrt_chi2 = c(l_r[1], l_p[1]), lrt_df = c(l_r[2], l_p[2]), lrt_p = pchisq(c(l_r[1], l_p[1]), c(l_r[2], l_p[2]), lower.tail = FALSE)),
      tibble(row_type = "constraint", term = "b_PBW = 0 and b_VT + b_PFVC = 0", term_label = "only strain VT/PFVC matters (the strain model)",
             lrt_chi2 = l_s[1], lrt_df = l_s[2], lrt_p = pchisq(l_s[1], l_s[2], lower.tail = FALSE))) %>%
      mutate(outcome = oc, adjustment = paste("unadjusted:", sz_saturated_covariates),
             estimate_type = if (inherits(f_sat, "coxph")) "HR" else "OR",
             n = model_patients(f_sat), n_events = model_events(f_sat), aic_saturated = AIC(f_sat), .before = 1)
  }) %>% mutate(scale = "per log unit (x1.1 = +10%); constraint rows: exp(sum of coefficients), 1 = constraint holds", site = site_name)
  write_csv(size_saturated, file.path(final_dir, paste0("size_saturated_log_model_", site_name, ".csv")))

  # --- (B) functional-form ladder for the size term -------------------------------------
  # One size term per model: a form holding log PBW beside log PBW/PFVC would put log PBW
  # and log PFVC together with the demographics, which is not identified (see (A)).
  sz_forms <- c("PFVC (linear)" = "pfvc", "log PFVC" = "l_pfvc", "1/PFVC" = "inv_pfvc",
                "PBW/PFVC (linear)" = "pbwpfvc", "log PBW/PFVC" = "l_ratio")
  size_form <- map_dfr(sz_outcomes, function(oc) map_dfr(names(sz_age), function(af) {
    cov <- paste(sz_age[[af]], "+", sz_base)
    fits <- map(sz_forms, ~ sz_fit(paste("vtpbw +", .x, "+", cov), oc))
    aics <- map_dbl(fits, AIC)
    tibble(outcome = oc, age_form = af, size_term = names(sz_forms), rhs = paste("vtpbw +", unname(sz_forms)),
           AIC = aics, delta_AIC_vs_linear_pfvc = aics - aics[["PFVC (linear)"]],
           n = map_dbl(fits, model_patients))
  })) %>% mutate(site = site_name)
  write_csv(size_form, file.path(final_dir, paste0("size_functional_form_", site_name, ".csv")))

  cat("\n--- 4k(A) saturated log model (per log unit), demographic-unadjusted ---\n")
  print(as.data.frame(size_saturated %>% filter(row_type == "coefficient") %>%
          transmute(outcome, term_label, est = sprintf("%.2f [%.2f, %.2f]", estimate, conf_low, conf_high))), row.names = FALSE)
  cat("--- constraint tests (exp(sum) with CI; Wald p; LRT p) ---\n")
  print(as.data.frame(size_saturated %>% filter(row_type == "constraint") %>%
          transmute(outcome, term_label,
                    exp_sum = ifelse(is.na(estimate), NA, sprintf("%.2f [%.2f, %.2f]", estimate, conf_low, conf_high)),
                    wald_p = signif(p_value, 2), lrt_p = signif(lrt_p, 2))), row.names = FALSE)
  cat("--- 4k(B) functional form of the size term (delta AIC vs linear PFVC; negative = better) ---\n")
  print(as.data.frame(size_form %>% transmute(outcome, age_form, size_term, dAIC = round(delta_AIC_vs_linear_pfvc, 1))),
        row.names = FALSE)
} else {
  message("4k skipped: no outcome with >= ", NC_MIN_EVENTS, " deaths, or fewer than ",
          SATURATED_MIN_PATIENTS, " complete patients in the analytic cohort")
}

# =============================================================================
# 4l. Figure 2: the strain the protocol delivers inside the band
# =============================================================================
# VT/PFVC (% of predicted FVC) = VT/PBW (mL/kg) x PBW/PFVC (kg/L) / 10, so in logs the
# delivered strain is the clinician's dose plus the label's mis-sizing:
#   log VT/PFVC = log VT/PBW + log PBW/PFVC - log 10
#   Var(log VT/PFVC) = Var(log VT/PBW) + Var(log PBW/PFVC) + 2 Cov
# (A) the distribution of VT/PFVC inside the 6-8 mL/kg band, by sex, race, age band
#     and height band, as histograms on fixed 0.5-point bins (tails clamped into the
#     edge bins) so they sum across sites; the 11% line (ARMA's low-VT arm near its
#     75th percentile) is drawn centrally.
# (B) the variance decomposition, exported as moments (n, means, variances,
#     covariance) so the pooled decomposition is exact: pooled variance is the
#     within-site variance plus the spread of the site means.
VTPFVC_BIN_EDGES <- seq(4, 25, by = 0.5)
# The 6-8 mL/kg band itself is applied in script 03 (the analytic cohort's
# inclusion); this filter only drops rows whose logs are undefined.
dose_band <- cross_sectional %>%
  filter(vtpbw > 0, pbwpfvc > 0, vtpfvc > 0) %>%
  mutate(l_vtpbw = log(vtpbw), l_ratio = log(pbwpfvc), l_vtpfvc = log(vtpfvc))
message("4l frame (positive VT/PBW, PBW/PFVC and VT/PFVC): ",
        nrow(cross_sectional), " -> ", nrow(dose_band), " patients")
stopifnot(max(abs(dose_band$l_vtpfvc - (dose_band$l_vtpbw + dose_band$l_ratio - log(10)))) < 1e-6)
dose_groups <- bind_rows(
  dose_band %>% transmute(group_type = "overall", group_value = "all", value = vtpfvc),
  dose_band %>% transmute(group_type = "sex", group_value = as.character(sex_category), value = vtpfvc),
  dose_band %>% transmute(group_type = "race", group_value = as.character(race_category), value = vtpfvc),
  dose_band %>% transmute(group_type = "age_bin",
                          group_value = as.character(cut(age_at_admission, DOSE_AGE_BREAKS, right = FALSE)),
                          value = vtpfvc),
  dose_band %>% transmute(group_type = "height_bin",
                          group_value = as.character(cut(height_cm, DOSE_HEIGHT_BREAKS, right = FALSE)),
                          value = vtpfvc)) %>%
  filter(!is.na(group_value))
dose_histograms <- dose_groups %>%
  mutate(value = pmin(pmax(value, min(VTPFVC_BIN_EDGES)), max(VTPFVC_BIN_EDGES) - 1e-9),
         bin_i = findInterval(value, VTPFVC_BIN_EDGES, rightmost.closed = TRUE)) %>%
  count(group_type, group_value, bin_i, name = "count") %>%
  transmute(site = site_name, group_type, group_value, bin_left = VTPFVC_BIN_EDGES[bin_i],
            bin_right = VTPFVC_BIN_EDGES[bin_i + 1], count)
write_csv(dose_histograms, file.path(final_dir, paste0("dose_vtpfvc_histograms_", site_name, ".csv")))

dose_moments <- function(d, group_type, group_value) tibble(
  group_type = group_type, group_value = group_value, n = nrow(d),
  mean_log_vtpbw = mean(d$l_vtpbw), mean_log_ratio = mean(d$l_ratio),
  var_log_vtpbw = var(d$l_vtpbw), var_log_ratio = var(d$l_ratio), cov_log_vtpbw_ratio = cov(d$l_vtpbw, d$l_ratio))
dose_decomposition <- bind_rows(
  dose_moments(dose_band, "overall", "all"),
  map_dfr(c("sex_category", "race_category"), function(g)
    map_dfr(sort(unique(as.character(dose_band[[g]]))), function(v)
      dose_moments(dose_band[as.character(dose_band[[g]]) == v, ], sub("_category$", "", g), v)))) %>%
  mutate(var_log_vtpfvc = var_log_vtpbw + var_log_ratio + 2 * cov_log_vtpbw_ratio,
         share_clinician = var_log_vtpbw / var_log_vtpfvc,
         share_missizing = var_log_ratio / var_log_vtpfvc,
         share_covariance = 2 * cov_log_vtpbw_ratio / var_log_vtpfvc,
         site = site_name, .before = 1)
write_csv(dose_decomposition, file.path(final_dir, paste0("dose_variance_decomposition_", site_name, ".csv")))
with(dose_decomposition %>% filter(group_type == "overall"),
     message(sprintf("4l: VT/PFVC variance inside the band: %.0f%% mis-sizing, %.0f%% clinician, %.0f%% covariance (n = %d)",
                     100 * share_missizing, 100 * share_clinician, 100 * share_covariance, n)))

message("All outputs saved to: ", final_dir)
message("Script 04 complete.")
