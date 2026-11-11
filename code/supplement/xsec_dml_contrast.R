# =============================================================================
# Supplement (cross-sectional and marker): the ventilated-minus-control contrast with
# severity entered flexibly (double machine learning) and with the two cohorts
# weighted to their common severity (overlap weights)
# =============================================================================
# xsec_strain_invariance.R and the figure-4 marker models compare the strain error's
# association in ventilated patients with the same association in no-support
# controls, with no age, sex, race or height term (they are the strain error's own
# inputs) and severity entered linearly: SF and SOFA standardised within cohort for
# mortality, a SOFA anchor projected to the ventilated mean for the markers. This
# script asks whether that contrast depends on how severity enters. It adds no new
# test of strain: same cohorts, same exposures, same comparison, and still no
# demographic term anywhere.
#
# 1. Double machine learning (the partialling-out estimator of Chernozhukov et al.,
#    Econometrics J 2018). In each cohort separately,
#      outcome = theta x strain error + g(severity) + error,
#    with g learned by a flexible learner. The outcome and the strain error are each
#    predicted from severity on held-out folds (K-fold cross-fitting, repeated with
#    new folds), and theta is the slope of one residual on the other. The contrast is
#    theta(ventilated) - theta(control), the cohorts' variances summed, as in the
#    strain-invariance contrast. Two learners: a random forest (ranger) and a lasso on
#    a spline-and-interaction expansion of severity (glmnet). Beside them, the same
#    model with severity linear (ordinary least squares on the same patients and the
#    same outcome scale), so a change from the linear fit to the learners is the
#    flexibility and not the scale.
#    The out-of-fold R2 of the strain error on severity is reported: it is how much of
#    the strain error severity can stand in for. A flexible g can absorb demographic
#    signal that reaches severity (creatinine carries sex and muscle mass, for
#    example), so an attenuation under the learners can be fit or leakage; the R2 says
#    how much room there is for the second.
# 2. Overlap weights. Ventilated and control patients differ in severity, and the
#    marker models project the control onto the ventilated severity, which
#    extrapolates. A logistic model of P(ventilated | severity) gives each ventilated
#    patient the weight 1 - p and each control the weight p, which weights both
#    cohorts toward the patients whose severity both cohorts contain. The existing
#    linear models are refitted with those weights (robust SEs). The propensity
#    quantiles and effective sample sizes are written as the overlap diagnostic.
#
# Outcomes, one number per patient (the learners need one):
#   inhosp_death        in-hospital death (0/1); theta is a risk difference
#   day60_censor_imv    60-day death, controls censored at escalation to IMV: the
#                       pseudo-observation of 60-day death probability (1 - the
#                       jackknife Kaplan-Meier pseudo-value, per cohort), which carries
#                       the censoring without a proportional-hazards model
#   death_on_support    death during the exposure (ventilated: on IMV or within one day
#                       of the last IMV record, censored at liberation; control: in
#                       hospital before IMV, censored at escalation and discharge), as
#                       a pseudo-observation in the same way
#   platelets_slope,    the patient's least-squares slope of the log marker per day
#   creatinine_slope    over days 1-7 of the figure-4 panel (two or more values), the
#                       marker models' trajectory reduced to one number; creatinine
#                       is already censored at renal replacement in the panel
# Mortality uses xsec_strain_invariance.R's cohorts (copied from it, same synthetic
# seed and row order) and its two populations, everyone and intubation-eligible. The
# markers use the icu-clock panels of 21_biotrauma_panel.R (ventilated: icu_day0, the
# ungated panel under PBWPFVC_VTPBW_GATE=0; control: the no-support panel), everyone.
#
# Severity, measured at the index in both cohorts:
#   mortality  SF ratio and the six SOFA components
#   markers    SF and the marker's leave-one-out SOFA anchor components (20's
#              anchor_components(), as the marker models use), plus the log baseline
#              value of the marker itself
# Exposure: the strain error along the lambda axis of xsec_strain_invariance.R 2f,
#   log PBW/V(lambda) = log PBW/PFVC at 25 - lambda (log PFVC - log PFVC at 25),
# at PBWPFVC_DML_LAMBDAS (default 0, 0.5, 1). Estimates are per log unit; the _per_0.1
# columns are per 0.1 log units.
#
# Knobs: PBWPFVC_DML_LAMBDAS ("0,0.5,1"), PBWPFVC_DML_FOLDS (5), PBWPFVC_DML_REPS (3),
#        PBWPFVC_DML_TREES (500)
# Inputs : as xsec_strain_invariance.R, plus the icu-clock panels
#          jm_{long,surv}_<horizon>_icu[_ungated].parquet (ventilated) and
#          controls/nosupport/jm_{long,surv}_<horizon>_icu.parquet (21_biotrauma_panel.R)
# Outputs: final/supplement/ (final/ungated/supplement/ when ungated)
#   dml_contrast_{site}.csv   theta per cohort and the difference, by method (linear,
#                             forest, lasso, overlap-weighted linear), outcome,
#                             population and lambda, with nuisance R2 and patient counts
#   dml_overlap_{site}.csv    the propensity model's quantiles by cohort and the
#                             effective sample size under the overlap weights
# Usage: PBWPFVC_VTPBW_GATE=0 uvr run code/supplement/xsec_dml_contrast.R   (PBWPFVC_COHORT unset)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(survival)
  library(splines)
})

options(width = 220)
source("utils/config.R")
if (config$cohort != "imv") stop("xsec_dml_contrast.R reads both cohorts itself: unset PBWPFVC_COHORT")
site_name <- config$site_name
final_dir <- final_dir_for("supplement")
source(here::here("code", "20_biotrauma_grid.R"))   # anchor_components(), h_suffix

parse_number_list <- function(text) as.numeric(strsplit(text, ",")[[1]])
DML_LAMBDAS <- parse_number_list(Sys.getenv("PBWPFVC_DML_LAMBDAS", "0,0.5,1"))
DML_FOLDS   <- as.integer(Sys.getenv("PBWPFVC_DML_FOLDS", "5"))
DML_REPS    <- as.integer(Sys.getenv("PBWPFVC_DML_REPS", "3"))
DML_TREES   <- as.integer(Sys.getenv("PBWPFVC_DML_TREES", "500"))
DML_SEED    <- 20261028
FOREST_MIN_NODE <- 20L          # ranger minimum node size: smooths the outcome surface a little
MIN_PATIENTS <- 50L             # a cohort-outcome cell below this is skipped, not fitted
MIN_EVENTS <- 10L               # deaths needed for a mortality outcome (xsec_strain_invariance.R)
SAME_DAY_DEATH_D <- 0.5
HORIZON_DAYS <- 60
DEATH_ON_SUPPORT_GRACE_D <- 1
SOFA_COMPONENTS <- c("sofa_cv_97", "sofa_coag", "sofa_liver", "sofa_resp", "sofa_cns", "sofa_renal")
COHORTS <- c("Ventilated", "No support")
DIFFERENCE <- "ventilated minus no support"
if (anyNA(DML_LAMBDAS) || any(DML_LAMBDAS < 0 | DML_LAMBDAS > 1)) stop("PBWPFVC_DML_LAMBDAS must be numbers in [0, 1]")
message("=== xsec_dml_contrast (", site_name, if (config$vtpbw_gate == "0") ", ungated" else ", gated",
        "): lambda ", paste(DML_LAMBDAS, collapse = ", "), "; ", DML_FOLDS, "-fold cross-fitting x ", DML_REPS,
        " repetitions; forest ", DML_TREES, " trees ===")

# =============================================================================
# Data, mortality: both cross-sectional cohorts, as xsec_strain_invariance.R builds them
# =============================================================================
cohort_columns <- c("hospitalization_id", "index_dttm", "age_at_admission", "sex_category",
                    "race_category", "pfvc", "pfvc_age25", "pva", "ptlc", "pfrc", "pbw", "sf_ratio", "sofa_total", "deceased", "admission_dttm", "death_dttm",
                    "discharge_dttm", "height_cm")
control_file <- file.path(config$output_dir, "controls", "nosupport", "analysis_cross_sectional.parquet")
if (!file.exists(control_file))
  stop("no no-support cohort: run scripts 01-03 with PBWPFVC_COHORT=nosupport first")
ventilated_cohort <- read_parquet(file.path(config$output_dir, paste0("analysis_cross_sectional", config$cs_suffix, ".parquet"))) %>%
  select(all_of(cohort_columns), all_of(SOFA_COMPONENTS), vtpbw, icu_day0, patient_id, last_imv_dttm)
ventilated <- ventilated_cohort %>% filter(icu_day0) %>% select(-icu_day0) %>%
  mutate(cohort = "Ventilated", escalation_dttm = as.POSIXct(NA))
no_support <- read_parquet(control_file) %>%
  select(all_of(cohort_columns), all_of(SOFA_COMPONENTS), escalation_dttm, patient_id) %>%
  mutate(cohort = "No support", vtpbw = NA_real_, last_imv_dttm = as.POSIXct(NA))
# no patient in both arms (by patient_id, any hospitalization), as xsec_strain_invariance.R
shared_patients <- intersect(no_support$patient_id, ventilated$patient_id)
message("Controls also in the ventilated arm, dropped: ", length(shared_patients))
# the control hospitalizations kept, for the marker panels below
control_hospitalizations <- no_support %>% filter(!patient_id %in% shared_patients) %>% pull(hospitalization_id)
no_support <- no_support %>% filter(!patient_id %in% shared_patients) %>% select(-patient_id)
ventilated <- ventilated %>% select(-patient_id)
control_imv <- read_parquet(file.path(config$output_dir, "controls", "nosupport", "resp_support_waterfall_clean.parquet"),
                            col_select = c("hospitalization_id", "recorded_dttm", "device_category", "tidal_volume_set")) %>%
  filter(tolower(device_category) == "imv" | (!is.na(tidal_volume_set) & tidal_volume_set > 0)) %>%
  inner_join(no_support %>% select(hospitalization_id, index_dttm), by = "hospitalization_id") %>%
  filter(recorded_dttm >= index_dttm) %>%
  group_by(hospitalization_id) %>% summarise(imv_dttm = min(recorded_dttm), .groups = "drop")
no_support <- no_support %>% left_join(control_imv, by = "hospitalization_id")
both_cohorts <- bind_rows(ventilated, no_support)

# SYNTHETIC SITE ONLY: simulated mortality with xsec_strain_invariance.R's seed and row
# order (synthetic CLIF mortality is unreliable). Plumbing only; never runs at a real site.
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
if (any(both_cohorts$deceased == 1 & is.na(both_cohorts$death_dttm), na.rm = TRUE))
  stop("expired patients without a death time: rebuild the cohorts with script 03, which dates them at discharge")

index_day <- function(dttm, index) as.numeric(difftime(dttm, index, units = "days"))
# the strain error along the lambda axis, from PBW, PFVC and PFVC at 25
add_lambda_exposures <- function(dat, pfvc_column) {
  for (lambda in DML_LAMBDAS)
    dat[[lambda_column(lambda)]] <- log(dat$pbw / dat$pfvc_age25) - lambda * (log(dat[[pfvc_column]]) - log(dat$pfvc_age25))
  dat
}
lambda_column <- function(lambda) sprintf("strain_error_lambda%03d", round(100 * lambda))

n_before_missing_data <- both_cohorts %>% count(cohort, name = "n_before")
both_cohorts <- both_cohorts %>%
  filter(!is.na(pfvc), pfvc > 0, !is.na(pfvc_age25), pfvc_age25 > 0, !is.na(pbw), pbw > 0, !is.na(sf_ratio), !is.na(sofa_total),
         !is.na(deceased), !is.na(discharge_dttm), if_all(all_of(SOFA_COMPONENTS), ~ !is.na(.x)),
         # the strain-invariance cohort's exclusions, so the patients match it
         !is.na(age_at_admission), !is.na(sex_category), !is.na(race_category), !is.na(height_cm), height_cm > 0) %>%
  group_by(cohort) %>%
  mutate(sf_z = as.numeric(scale(sf_ratio)), sofa_z = as.numeric(scale(sofa_total))) %>%
  ungroup() %>%
  mutate(cohort = factor(cohort, levels = COHORTS),
         death_stamped_before_index = !is.na(death_dttm) & death_dttm < index_dttm,
         death_index_day     = if_else(death_stamped_before_index, SAME_DAY_DEATH_D, index_day(death_dttm, index_dttm)),
         discharge_index_day = index_day(discharge_dttm, index_dttm),
         imv_day             = index_day(imv_dttm, index_dttm),
         last_imv_day        = index_day(last_imv_dttm, index_dttm)) %>%
  add_lambda_exposures("pfvc")
missing_data_counts <- n_before_missing_data %>%
  left_join(both_cohorts %>% count(cohort = as.character(cohort), name = "n_after"), by = "cohort")
message("Excluded for missing PFVC, PBW, SF, SOFA (or a component), outcome, demographics, height or discharge time:")
print(as.data.frame(missing_data_counts), row.names = FALSE)

# Goals of care: the intubation-eligible population keeps every ventilated patient and
# drops the controls with a documented limitation (xsec_strain_invariance.R, copied).
CODE_STATUS_WINDOW_H <- 24
FULL_CODE_CATEGORIES <- c("full", "presume full")
code_status_file <- file.path(path.expand(config$tables_path), paste0("clif_code_status.", config$file_type))
HAS_CODE_STATUS <- file.exists(code_status_file)
if (HAS_CODE_STATUS) {
  read_clif_table <- function(table_name, columns) {
    path <- file.path(path.expand(config$tables_path), paste0("clif_", table_name, ".", config$file_type))
    switch(config$file_type,
           parquet = read_parquet(path, col_select = all_of(columns)),
           csv     = readr::read_csv(path, col_select = all_of(columns), show_col_types = FALSE),
           fst     = fst::read_fst(path, columns = columns))
  }
  site_hospitalizations <- read_clif_table("hospitalization", c("patient_id", "hospitalization_id")) %>%
    filter(hospitalization_id %in% both_cohorts$hospitalization_id)
  if (nrow(site_hospitalizations) == 0)
    stop("clif_hospitalization at ", config$tables_path, " holds none of this site's cohort hospitalizations: ",
         "set PBWPFVC_TABLES_PATH (or config.json) to ", site_name, "'s CLIF tables.")
  status_at_index <- read_clif_table("code_status", c("patient_id", "start_dttm", "code_status_category")) %>%
    inner_join(site_hospitalizations, by = "patient_id", relationship = "many-to-many") %>%
    inner_join(both_cohorts %>% transmute(cohort, hospitalization_id, admission_dttm, index_dttm),
               by = "hospitalization_id", relationship = "many-to-many") %>%
    filter(start_dttm >= admission_dttm, start_dttm <= index_dttm + CODE_STATUS_WINDOW_H * 3600) %>%
    group_by(cohort, hospitalization_id) %>% slice_max(start_dttm, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(cohort, hospitalization_id,
              code_status_at_index = if_else(tolower(code_status_category) %in% FULL_CODE_CATEGORIES, "full code", "limited or other"))
  both_cohorts <- both_cohorts %>% left_join(status_at_index, by = c("cohort", "hospitalization_id")) %>%
    mutate(code_status_at_index = coalesce(code_status_at_index, "no record"))
} else {
  message("\n*** No code_status table at ", config$tables_path, ": the intubation-eligible population is skipped. ***")
  both_cohorts <- both_cohorts %>% mutate(code_status_at_index = "no record")
}
both_cohorts <- both_cohorts %>% mutate(intubation_eligible = cohort == "Ventilated" | code_status_at_index != "limited or other")
POPULATIONS <- c(everyone = "everyone",
                 if (HAS_CODE_STATUS) c(eligible = "intubation-eligible (every ventilated patient; controls without a documented limitation)"))

# Each mortality outcome's follow-up, as xsec_strain_invariance.R defines it, then the
# per-patient value: in-hospital death as 0/1, the two time-to-event outcomes as
# pseudo-observations of the probability of death by day 60, computed within cohort.
MORTALITY_OUTCOMES <- c(inhosp_death = "in-hospital death (risk difference)",
                        day60_censor_imv = "60-day death, controls censored at escalation to IMV (pseudo-observation, probability)",
                        death_on_support = "death during the exposure (pseudo-observation, probability)")
survival_followup <- function(dat, outcome_key) dat %>% mutate(
  censor_day = case_when(outcome_key == "day60_censor_imv" ~ pmin(HORIZON_DAYS, coalesce(imv_day, Inf)),
                         cohort == "Ventilated" ~ pmin(HORIZON_DAYS, coalesce(last_imv_day, 0) + DEATH_ON_SUPPORT_GRACE_D),
                         TRUE ~ pmin(HORIZON_DAYS, coalesce(imv_day, Inf), coalesce(discharge_index_day, Inf))),
  event = as.integer(!is.na(death_index_day) & death_index_day <= censor_day),
  end_day = pmax(if_else(event == 1L, death_index_day, censor_day), 0.01))
mortality_outcome <- function(dat, outcome_key) {
  if (outcome_key == "inhosp_death") return(dat %>% mutate(outcome_value = deceased, n_events = sum(deceased == 1)))
  dat <- survival_followup(dat, outcome_key)
  km <- survfit(Surv(end_day, event) ~ 1, data = dat)
  # pseudo() rebuilds the model frame by re-evaluating the call in another
  # environment, where `dat` is not this one; hand it the frame instead
  km$model <- model.frame(Surv(end_day, event) ~ 1, data = dat)
  survival_pseudo <- as.numeric(pseudo(km, times = HORIZON_DAYS, type = "surv"))
  dat %>% mutate(outcome_value = 1 - survival_pseudo, n_events = sum(event))
}

# =============================================================================
# Data, markers: the icu-clock panels, one slope per patient
# =============================================================================
MARKERS <- c(platelets = "platelets", creatinine = "creatinine")
MARKER_BASELINE <- c(platelets = "platelet_0", creatinine = "creatinine_0")
MARKER_OUTCOMES <- c(platelets_slope = "platelets, log per day over days 1-7",
                     creatinine_slope = "creatinine, log per day over days 1-7 (censored at renal replacement)")
ventilated_panel_sfx <- paste0("_icu", if (config$vtpbw_gate == "0") "_ungated" else "")
panel_file <- function(kind, dir, sfx) file.path(dir, paste0("jm_", kind, "_", h_suffix, sfx, ".parquet"))
panel_dirs <- list(Ventilated = list(dir = config$output_dir, sfx = ventilated_panel_sfx),
                   `No support` = list(dir = file.path(config$output_dir, "controls", "nosupport"), sfx = "_icu"))
for (cohort_now in COHORTS) for (kind in c("long", "surv"))
  if (!file.exists(panel_file(kind, panel_dirs[[cohort_now]]$dir, panel_dirs[[cohort_now]]$sfx)))
    stop("no ", cohort_now, " icu-clock panel (", basename(panel_file(kind, panel_dirs[[cohort_now]]$dir, panel_dirs[[cohort_now]]$sfx)),
         "): run 21_biotrauma_panel.R with PBWPFVC_JM_CLOCK=icu for this cohort first (29_run_figure4.R does)")
marker_patients <- map(set_names(COHORTS), function(cohort_now) {
  location <- panel_dirs[[cohort_now]]
  surv <- read_parquet(panel_file("surv", location$dir, location$sfx))
  long <- read_parquet(panel_file("long", location$dir, location$sfx))
  surv <- if (cohort_now == "Ventilated") surv %>% filter(icu_day0)
          else surv %>% filter(hospitalization_id %in% control_hospitalizations)   # the shared-patient rule above
  # the slope of the log marker on day, days 1 to the horizon, two or more values
  slopes <- map(names(MARKERS), function(marker) long %>%
                  filter(period >= 1L, !is.na(.data[[MARKERS[[marker]]]]), .data[[MARKERS[[marker]]]] > 0) %>%
                  transmute(hospitalization_id, day = vent_day, log_y = log(.data[[MARKERS[[marker]]]])) %>%
                  group_by(hospitalization_id) %>%
                  filter(n() >= 2L, var(day) > 0) %>%
                  summarise(!!paste0(marker, "_slope") := sum((day - mean(day)) * (log_y - mean(log_y))) / sum((day - mean(day))^2),
                            .groups = "drop")) %>%
    reduce(full_join, by = "hospitalization_id")
  surv %>% inner_join(slopes, by = "hospitalization_id") %>%
    mutate(cohort = cohort_now) %>%
    add_lambda_exposures("pfvc_gli")
})
message("Marker panels: ", paste(sprintf("%s %d patients with a slope", names(marker_patients), map_int(marker_patients, nrow)),
                                collapse = "; "))

# =============================================================================
# Estimators
# =============================================================================
# Severity design for the lasso: each continuous column as a 4-df natural spline,
# each column with five or fewer distinct values (a SOFA component) as indicators,
# and every pairwise product of the raw columns.
lasso_design <- function(severity) {
  pieces <- map(names(severity), function(column) {
    x <- severity[[column]]
    if (n_distinct(x) <= 5L) { m <- model.matrix(~ factor(x))[, -1, drop = FALSE] }
    else m <- ns(x, df = 4)
    colnames(m) <- paste0(column, "_", seq_len(ncol(m))); m
  })
  pairwise <- model.matrix(~ .^2 - 1, data = as.data.frame(severity))
  pairwise <- pairwise[, grepl(":", colnames(pairwise)), drop = FALSE]
  cbind(do.call(cbind, pieces), pairwise)
}
# Out-of-fold predictions of `target` from severity, for one fold assignment
predict_out_of_fold <- function(target, severity, design, folds, learner) {
  predictions <- numeric(length(target))
  for (fold in sort(unique(folds))) {
    train <- folds != fold
    predictions[!train] <- if (learner == "forest") {
      fit <- ranger::ranger(x = severity[train, , drop = FALSE], y = target[train], num.trees = DML_TREES,
                            min.node.size = FOREST_MIN_NODE, num.threads = 1, seed = DML_SEED + fold)
      predict(fit, data = severity[!train, , drop = FALSE])$predictions
    } else {
      fit <- glmnet::cv.glmnet(design[train, , drop = FALSE], target[train], alpha = 1, nfolds = 5)
      as.numeric(predict(fit, newx = design[!train, , drop = FALSE], s = "lambda.min"))
    }
  }
  predictions
}
out_of_fold_r2 <- function(target, predictions) 1 - sum((target - predictions)^2) / sum((target - mean(target))^2)
# Partialling-out estimate of theta for each exposure: the outcome's residual on each
# exposure's residual, per repetition, aggregated by the median with the
# between-repetition spread added to the variance (Chernozhukov et al. 2018, sec 3.4).
# The outcome nuisance is fitted once per repetition and shared by the exposures.
dml_thetas <- function(outcome, exposures, severity, learner) {
  severity <- as.data.frame(severity)
  design <- if (learner == "lasso") lasso_design(severity) else NULL
  n <- length(outcome)
  per_rep <- map_dfr(seq_len(DML_REPS), function(rep) {
    set.seed(DML_SEED + rep)
    folds <- sample(rep(seq_len(DML_FOLDS), length.out = n))
    outcome_hat <- predict_out_of_fold(outcome, severity, design, folds, learner)
    outcome_residual <- outcome - outcome_hat
    map_dfr(names(exposures), function(exposure) {
      exposure_hat <- predict_out_of_fold(exposures[[exposure]], severity, design, folds, learner)
      exposure_residual <- exposures[[exposure]] - exposure_hat
      theta <- sum(exposure_residual * outcome_residual) / sum(exposure_residual^2)
      score <- (outcome_residual - theta * exposure_residual) * exposure_residual
      tibble(exposure, rep, theta, se = sqrt(mean(score^2) / mean(exposure_residual^2)^2 / n),
             r2_outcome_on_severity = out_of_fold_r2(outcome, outcome_hat),
             r2_exposure_on_severity = out_of_fold_r2(exposures[[exposure]], exposure_hat))
    })
  })
  per_rep %>% group_by(exposure) %>%
    summarise(estimate = median(theta), se = sqrt(median(se^2 + (theta - median(theta))^2)),
              r2_outcome_on_severity = median(r2_outcome_on_severity),
              r2_exposure_on_severity = median(r2_exposure_on_severity), .groups = "drop")
}
# Least squares with severity linear, optionally weighted; HC0 sandwich SEs, which are
# robust for the weighted fit and for the 0/1 and pseudo-observation outcomes alike
linear_theta <- function(outcome, exposure, severity, weights = rep(1, length(outcome))) {
  design <- cbind(intercept = 1, exposure = exposure, as.matrix(severity))
  weighted_cross <- crossprod(design, design * weights)
  bread <- solve(weighted_cross)
  coefficients <- bread %*% crossprod(design, weights * outcome)
  residuals <- as.numeric(outcome - design %*% coefficients)
  meat <- crossprod(design * (weights * residuals))
  tibble(estimate = coefficients["exposure", 1], se = sqrt((bread %*% meat %*% bread)["exposure", "exposure"]))
}
# Overlap weights from P(ventilated | severity): logistic, splines in the continuous
# columns, SOFA components linear
overlap_weights <- function(dat, severity_columns) {
  continuous <- severity_columns[map_lgl(severity_columns, ~ n_distinct(dat[[.x]]) > 5L)]
  ordinal <- setdiff(severity_columns, continuous)
  rhs <- paste(c(sprintf("ns(%s, 4)", continuous), ordinal), collapse = " + ")
  # Fitted probabilities of 0 or 1 are allowed: they are patients outside the common
  # severity, whom the overlap weights set to about zero, and the overlap table reports
  # how many there are (share_below_0.05, share_above_0.95). Any other warning stops.
  propensity_fit <- withCallingHandlers(
    glm(as.formula(paste("ventilated ~", rhs)), family = binomial, data = dat %>% mutate(ventilated = as.integer(cohort == "Ventilated"))),
    warning = function(w) {
      if (grepl("fitted probabilities numerically 0 or 1", conditionMessage(w))) invokeRestart("muffleWarning")
      stop("propensity model warning, stopping: ", conditionMessage(w), call. = FALSE)
    })
  propensity <- fitted(propensity_fit)
  dat %>% mutate(propensity = propensity, overlap_weight = if_else(cohort == "Ventilated", 1 - propensity, propensity))
}

# One outcome in one population: every method and lambda, both cohorts and the difference
METHODS <- c(linear = "severity linear (least squares)", forest = "double ML, random forest",
             lasso = "double ML, lasso on splines and interactions",
             overlap = "severity linear, overlap-weighted to the common severity")
contrast_rows <- function(dat, outcome_key, outcome_label, population_label, severity_columns) {
  dat <- dat %>% filter(!is.na(outcome_value), if_all(all_of(severity_columns), ~ !is.na(.x)))
  if (any(table(factor(dat$cohort, levels = COHORTS)) < MIN_PATIENTS)) {
    message("  ", outcome_key, " (", population_label, "): fewer than ", MIN_PATIENTS, " patients in a cohort, skipped")
    return(NULL)
  }
  exposure_columns <- set_names(map_chr(DML_LAMBDAS, lambda_column), DML_LAMBDAS)
  weighted <- overlap_weights(dat, severity_columns)
  per_cohort <- map_dfr(COHORTS, function(cohort_now) {
    cohort_data <- dat %>% filter(cohort == cohort_now)
    cohort_weighted <- weighted %>% filter(cohort == cohort_now)
    severity <- cohort_data %>% select(all_of(severity_columns))
    exposures <- map(exposure_columns, ~ cohort_data[[.x]])
    linear <- imap_dfr(exposures, ~ linear_theta(cohort_data$outcome_value, .x, severity) %>% mutate(exposure = .y)) %>%
      mutate(method = "linear")
    overlap <- imap_dfr(exposure_columns, ~ linear_theta(cohort_weighted$outcome_value, cohort_weighted[[.x]],
                                                         cohort_weighted %>% select(all_of(severity_columns)),
                                                         cohort_weighted$overlap_weight) %>% mutate(exposure = .y)) %>%
      mutate(method = "overlap")
    learners <- map_dfr(c("forest", "lasso"), ~ dml_thetas(cohort_data$outcome_value, exposures, severity, .x) %>% mutate(method = .x))
    bind_rows(linear, overlap, learners) %>%
      mutate(quantity = cohort_now, n_patients = nrow(cohort_data),
             n_events = if ("n_events" %in% names(cohort_data)) cohort_data$n_events[1] else NA_integer_,
             effective_n_overlap = sum(cohort_weighted$overlap_weight)^2 / sum(cohort_weighted$overlap_weight^2))
  })
  difference <- per_cohort %>%
    select(method, exposure, quantity, estimate, se) %>%
    pivot_wider(names_from = quantity, values_from = c(estimate, se)) %>%
    transmute(method, exposure, quantity = DIFFERENCE,
              estimate = `estimate_Ventilated` - `estimate_No support`,
              se = sqrt(`se_Ventilated`^2 + `se_No support`^2))
  bind_rows(per_cohort, difference) %>%
    mutate(outcome = outcome_label, outcome_key = outcome_key, population = population_label,
           lambda = as.numeric(exposure), severity = paste(severity_columns, collapse = " + "))
}

# =============================================================================
# Fits
# =============================================================================
results <- list()
overlap_diagnostics <- list()
mortality_severity <- c("sf_ratio", SOFA_COMPONENTS)
for (population in names(POPULATIONS)) for (outcome_key in names(MORTALITY_OUTCOMES)) {
  message("Mortality: ", outcome_key, ", ", POPULATIONS[[population]])
  dat <- both_cohorts %>% filter(if (population == "eligible") intubation_eligible else TRUE)
  dat <- map_dfr(COHORTS, ~ mortality_outcome(dat %>% filter(cohort == .x), outcome_key))
  if (any(dat %>% group_by(cohort) %>% summarise(e = first(n_events)) %>% pull(e) < MIN_EVENTS)) {
    message("  fewer than ", MIN_EVENTS, " deaths in a cohort, skipped"); next
  }
  results[[length(results) + 1L]] <- contrast_rows(dat, outcome_key, MORTALITY_OUTCOMES[[outcome_key]], POPULATIONS[[population]],
                                                   mortality_severity)
  if (outcome_key == "inhosp_death")
    overlap_diagnostics[[length(overlap_diagnostics) + 1L]] <- overlap_weights(dat, mortality_severity) %>%
      mutate(analysis = paste("mortality,", POPULATIONS[[population]]))
}
for (marker in names(MARKERS)) {
  outcome_key <- paste0(marker, "_slope")
  severity_columns <- c("sf_0", anchor_components(marker), "log_baseline")
  message("Marker: ", outcome_key, " (severity: ", paste(severity_columns, collapse = " + "), ")")
  dat <- map_dfr(marker_patients, ~ .x %>%
                   transmute(hospitalization_id, cohort, outcome_value = .data[[outcome_key]], sf_0,
                             across(all_of(anchor_components(marker))),
                             log_baseline = log(.data[[MARKER_BASELINE[[marker]]]]),
                             across(all_of(map_chr(DML_LAMBDAS, lambda_column))))) %>%
    filter(is.finite(log_baseline)) %>%
    mutate(cohort = factor(cohort, levels = COHORTS))
  results[[length(results) + 1L]] <- contrast_rows(dat, outcome_key, MARKER_OUTCOMES[[outcome_key]], "everyone", severity_columns)
  overlap_diagnostics[[length(overlap_diagnostics) + 1L]] <-
    overlap_weights(dat %>% filter(!is.na(outcome_value), if_all(all_of(severity_columns), ~ !is.na(.x))), severity_columns) %>%
    mutate(analysis = paste("marker,", marker))
}

dml_contrast <- bind_rows(results) %>%
  mutate(method_label = METHODS[method],
         estimate_per_0.1 = 0.1 * estimate, lo_per_0.1 = 0.1 * (estimate - 1.96 * se), hi_per_0.1 = 0.1 * (estimate + 1.96 * se),
         z = estimate / se, p = 2 * pnorm(-abs(z)),
         scale = if_else(grepl("slope", outcome_key), "change in log marker per day per log unit of strain error",
                         "change in probability of death per log unit of strain error"),
         folds = DML_FOLDS, repetitions = DML_REPS, trees = DML_TREES, vtpbw_gate = config$vtpbw_gate != "0", site = site_name) %>%
  select(outcome, outcome_key, population, lambda, method, method_label, quantity, estimate, se, estimate_per_0.1,
         lo_per_0.1, hi_per_0.1, z, p, r2_outcome_on_severity, r2_exposure_on_severity, n_patients, n_events,
         effective_n_overlap, severity, scale, folds, repetitions, trees, vtpbw_gate, site)
overlap_summary <- bind_rows(overlap_diagnostics) %>%
  group_by(analysis, cohort) %>%
  summarise(n_patients = n(),
            propensity_p01 = quantile(propensity, 0.01), propensity_p10 = quantile(propensity, 0.10),
            propensity_median = median(propensity), propensity_p90 = quantile(propensity, 0.90),
            propensity_p99 = quantile(propensity, 0.99),
            share_below_0.05 = mean(propensity < 0.05), share_above_0.95 = mean(propensity > 0.95),
            effective_n_overlap = sum(overlap_weight)^2 / sum(overlap_weight^2), .groups = "drop") %>%
  mutate(site = site_name)

message("\nVentilated minus control, per 0.1 log units of the strain error:")
print(as.data.frame(dml_contrast %>% filter(quantity == DIFFERENCE) %>%
                      transmute(outcome_key, population = substr(population, 1, 12), lambda, method,
                                estimate_per_0.1 = signif(estimate_per_0.1, 3), lo = signif(lo_per_0.1, 3),
                                hi = signif(hi_per_0.1, 3), z = round(z, 2))), row.names = FALSE)
message("\nOverlap of severity between the cohorts (propensity of ventilation):")
print(as.data.frame(overlap_summary %>% mutate(across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)

write_csv(dml_contrast, file.path(final_dir, paste0("dml_contrast_", site_name, ".csv")))
write_csv(overlap_summary, file.path(final_dir, paste0("dml_overlap_", site_name, ".csv")))
message("\nWrote dml_contrast_", site_name, ".csv and dml_overlap_", site_name, ".csv to ", final_dir)
