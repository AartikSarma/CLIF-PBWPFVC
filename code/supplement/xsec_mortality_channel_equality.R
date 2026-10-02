# =============================================================================
# Supplement (cross-sectional): do sex and race predict death as their share of
# predicted lung size says they should?
# =============================================================================
# GLI-2012 log PFVC is, to a small remainder, a sum of four pieces (pfvc_channels()
# in 20_biotrauma_grid.R): a height term, an age curve, a sex shift and a race shift.
# The age piece is a smooth curve that a 4-df age spline reproduces, so the spline
# holds it here, with age's direct path to death. What is left is height, sex and race,
# and the question this script asks of them is whether one log unit of predicted lung
# size carries the same mortality association whichever of the three moved it.
#
# Model, per cohort:
#   death ~ [VT/PBW] + SF (z) + SOFA (z) + ns(age, 4) + ch_height + sex + race
# Sex and race enter as factors, so each group's coefficient is free: race is not
# forced into GLI's proportions between its Black and Other shifts. GLI fixes each
# group's shift in log PFVC (s, a constant: female against male, and each race against
# white, at the reference height and age), so a group's coefficient divided by its
# shift is its implied slope per log unit of PFVC. If predicted lung size is what
# carries the association, every implied slope equals the height slope:
#   sex     gamma_female / s_female = beta_height         1 df
#   Black   gamma_black  / s_black  = beta_height         1 df
#   Other   gamma_other  / s_other  = beta_height         1 df
#   race    Black and Other together                      2 df
#   all     sex and both race levels                      3 df (Wald, and an exact
#           likelihood-ratio test against the model with one column, ch_height +
#           ch_sex + ch_race, which lies in the free model's span; the two should
#           agree, and a gap between them points to misspecification)
#
# What the tests mean. A pass says the group's mortality association is what its PFVC
# shift predicts at the height slope. It would mislead in two ways: height's own
# non-lung paths to death can contaminate the height slope, the yardstick, or a direct
# effect of sex or race can match that contamination by coincidence. A failure says the
# group has a path to death other than predicted lung size, or that the height slope is
# contaminated. The no-support control, fitted with the same model, separates the two:
# its patients carry the same formulas and the same demographic paths to death, but no
# PBW-scaled tidal volume. The ventilated-minus-control difference of each contrast is
# reported beside the cohorts' own.
#
# The unadjusted counterpart. This model is the demographic-adjusted form: age by its
# spline, sex and race by their factors. The form without demographics, all four
# pieces with their own coefficients, is the channel section of xsec_pfvc_age_control.R.
#
# Details:
#   - ch_height uses GLI's male height exponent for everyone (the reference patient is
#     a man; the female exponent is about 6% smaller). The project's decomposition keeps
#     it; the difference falls in the remainder, which this model leaves out.
#   - The cohorts and their clocks are xsec_pfvc_age_control.R's, copied from it: the
#     ventilated arm is the ventilated cohort on invasive ventilation at ICU admission
#     (script 03's icu_day0), the control the no-support cohort; every clock starts at
#     the index; a death stamped between admission and the index counts on the index
#     day, at SAME_DAY_DEATH_D days.
#   - Outcomes: in-hospital death (logistic); 60-day death (cause-specific Cox); 60-day
#     death before invasive ventilation (Cox; identical to 60-day death in the
#     ventilated arm, which is intubated at the index, and the control outcome that
#     keeps strain out).
#   - A model that warns stops the script; separation is reported as not estimable in
#     the output row, and a model with fewer than MIN_EVENTS deaths is skipped with a note.
#     No model is replaced by a simpler one.
#
# Inputs : intermediate/analysis_cross_sectional.parquet (script 03, ventilated)
#          intermediate/controls/nosupport/analysis_cross_sectional.parquet
#          intermediate/controls/nosupport/resp_support_waterfall_clean.parquet
#          (scripts 01-03 run with PBWPFVC_COHORT=nosupport)
# Outputs: final/supplement/
#   mortality_channel_equality_slopes_{site}.csv  per cohort and outcome: the height
#                                                 slope, each group's coefficient, shift
#                                                 and implied slope, and each implied
#                                                 slope minus the height slope; the
#                                                 ventilated-minus-control rows; the
#                                                 one-column and log-PFVC slopes for
#                                                 reference; patients and deaths per group
#   mortality_channel_equality_tests_{site}.csv   the Wald tests (and the LR test) per
#                                                 cohort, outcome and difference
#   mortality_channel_equality_vcov_{site}.csv    covariance of (height slope, implied
#                                                 slopes) per cohort and outcome, for
#                                                 pooling the tests across sites
# Usage: uvr run code/supplement/xsec_mortality_channel_equality.R   (PBWPFVC_COHORT unset)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(survival)
  library(splines)
})

options(width = 220)
source("utils/config.R")
if (config$cohort != "imv") stop("xsec_mortality_channel_equality.R reads both cohorts itself: unset PBWPFVC_COHORT")
site_name <- config$site_name
final_dir <- final_dir_for("supplement")
source(here::here("code", "20_biotrauma_grid.R"))   # pfvc_channels()

MIN_EVENTS <- 10L               # minimum deaths per model: the CLIF minimum-count standard
SAME_DAY_DEATH_D <- 0.5         # script 03's day for a death stamped between admission and the index
HORIZON_DAYS <- 60
COHORTS <- c("Ventilated", "No support")
DIFFERENCE <- "ventilated minus no support"

# =============================================================================
# Data: both cross-sectional cohorts, as xsec_pfvc_age_control.R builds them
# =============================================================================
cohort_columns <- c("hospitalization_id", "index_dttm", "age_at_admission", "sex_category",
                    "race_category", "pfvc", "sf_ratio", "sofa_total", "deceased", "admission_dttm", "death_dttm", "discharge_dttm", "height_cm")
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
# tidal volume, as 03's escalation rule): only invasive ventilation delivers a
# PBW-scaled tidal volume
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
# Patients missing PFVC, SF, SOFA, the outcome, a demographic, height or the discharge
# time are excluded; the counts before and after are printed per cohort.
n_before_missing_data <- both_cohorts %>% count(cohort, name = "n_before")
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
         log_pfvc      = log(pfvc),
         death_stamped_before_index = !is.na(death_dttm) & death_dttm < index_dttm,
         death_index_day     = if_else(death_stamped_before_index, SAME_DAY_DEATH_D, index_day(death_dttm, index_dttm)),
         discharge_index_day = index_day(discharge_dttm, index_dttm),
         escalation_day      = index_day(escalation_dttm, index_dttm),
         imv_day             = index_day(imv_dttm, index_dttm),   # control only; NA in the ventilated cohort
         age10 = age_at_admission / 10)
missing_data_counts <- n_before_missing_data %>%
  left_join(both_cohorts %>% count(cohort = as.character(cohort), name = "n_after"), by = "cohort") %>%
  mutate(n_after = coalesce(n_after, 0L), n_patients_excluded_missing_data = n_before - n_after)
message("Excluded for missing PFVC, SF, SOFA, outcome, demographics, height or discharge time:")
print(as.data.frame(missing_data_counts), row.names = FALSE)
if (anyNA(both_cohorts$vtpbw[both_cohorts$cohort == "Ventilated"]))
  stop("ventilated patients without VT/PBW in the cross-sectional table")

# The GLI pieces, computed once on both cohorts so they share one reference patient
# (the pooled median height and age, a white man). ch_size is the height, sex and race
# pieces together: log PFVC less its age piece and the remainder.
both_cohorts <- bind_cols(both_cohorts, pfvc_channels(both_cohorts, "log_pfvc")) %>%
  mutate(ch_size = ch_height + ch_sex + ch_race)

# Each group's shift in log PFVC, from the pieces themselves so it follows
# pfvc_channels()'s reference. The sex and race pieces are evaluated at the reference
# height and age, so each group has exactly one value.
group_shift <- function(piece, factor_column, level) {
  shift <- unique(round(both_cohorts[[piece]][both_cohorts[[factor_column]] == level], 12))
  if (length(shift) != 1) stop(piece, " takes ", length(shift), " values in group ", level, "; expected one")
  shift
}
GROUPS <- tribble(
  ~group,   ~term,                 ~piece,    ~factor_column,  ~level,
  "female", "sex_categoryFemale",  "ch_sex",  "sex_category",  "Female",
  "Black",  "race_categoryBLACK",  "ch_race", "race_category", "BLACK",
  "Other",  "race_categoryOTHER",  "ch_race", "race_category", "OTHER") %>%
  mutate(shift = pmap_dbl(list(piece, factor_column, level), group_shift))
message("\nGLI shifts in log PFVC at the reference height and age (female against male; each race against white):")
print(as.data.frame(GROUPS %>% select(group, shift) %>% mutate(shift = signif(shift, 4))), row.names = FALSE)
if (any(abs(GROUPS$shift) < 1e-6)) stop("a group has no GLI shift: its implied slope is undefined")

# =============================================================================
# The model and its tests
# =============================================================================
# The parameters the tests read, theta = (height slope, female, Black and Other implied
# slopes), are a linear map of the free model's coefficients b = (beta_height,
# gamma_female, gamma_black, gamma_other): theta = A b, A = diag(1, 1/s_female,
# 1/s_black, 1/s_other). Every test is a contrast D theta = 0, with D's rows the
# implied slope minus the height slope.
THETA <- c("height", GROUPS$group)
B_TERMS <- c("ch_height", GROUPS$term)
A <- diag(c(1, 1 / GROUPS$shift)); dimnames(A) <- list(THETA, B_TERMS)
D <- cbind(-1, diag(nrow(GROUPS))); dimnames(D) <- list(paste(GROUPS$group, "- height"), THETA)
TESTS <- list(
  `sex: female = height (1 df)`                 = "female - height",
  `race: Black = height (1 df)`                 = "Black - height",
  `race: Other = height (1 df)`                 = "Other - height",
  `race: Black = Other = height (2 df)`         = c("Black - height", "Other - height"),
  `all: female = Black = Other = height (3 df)` = rownames(D))
wald_test <- function(theta, V_theta, rows) {
  d <- D[rows, , drop = FALSE] %*% theta
  statistic <- as.numeric(t(d) %*% solve(D[rows, , drop = FALSE] %*% V_theta %*% t(D[rows, , drop = FALSE])) %*% d)
  tibble(statistic = statistic, df = length(rows), p = pchisq(statistic, length(rows), lower.tail = FALSE))
}

OUTCOMES <- tribble(
  ~outcome_key,          ~outcome,                                     ~model,
  "inhosp_logistic",     "in-hospital death (logistic)",               "logistic",
  "day60_all",           "60-day death, all",                          "cox",
  "day60_before_imv",    "60-day death, before invasive ventilation",  "cox")
# time (days from the index) and event for a Cox outcome, as xsec_pfvc_age_control.R
# defines them: "before invasive ventilation" censors at intubation, where a
# PBW-scaled volume begins
outcome_data <- function(cohort_data, outcome_key) {
  censor_imv <- grepl("before_imv", outcome_key)
  cohort_data %>% mutate(
    death_day = death_index_day,
    censor_day = pmin(HORIZON_DAYS, if (censor_imv) coalesce(imv_day, Inf) else Inf),
    event = as.integer(!is.na(death_day) & death_day <= censor_day),
    end_day = pmax(if_else(event == 1L, death_day, censor_day), 0.01))
}

# A model that warns stops the script; separation (a coefficient running to infinity,
# as when a small group has no deaths) is caught and reported as not estimable.
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

fit_cohort <- function(cohort_now, outcome_key, model) {
  cohort_data <- both_cohorts %>% filter(cohort == cohort_now)
  dat <- if (model == "logistic") cohort_data %>% mutate(event = deceased) else outcome_data(cohort_data, outcome_key)
  events <- sum(dat$event == 1)
  row_head <- tibble(cohort = cohort_now, outcome_key = outcome_key, n_patients = nrow(dat), n_deaths = events)
  if (events < MIN_EVENTS) return(list(head = row_head, note = paste("skipped: fewer than", MIN_EVENTS, "deaths")))
  base_rhs <- c(if (cohort_now == "Ventilated") "vtpbw", "sf_z", "sofa_z", "ns(age_at_admission, 4)")
  fit_one <- function(terms) {
    rhs <- paste(c(base_rhs, terms), collapse = " + ")
    fit_or_separation(if (model == "logistic")
      fit_strict(glm(as.formula(paste("event ~", rhs)), family = binomial, data = dat)) else
      fit_strict(coxph(as.formula(paste("Surv(end_day, event) ~", rhs)), data = dat)))
  }
  free_fit        <- fit_one(c("ch_height", "sex_category", "race_category"))
  one_column_fit  <- fit_one("ch_size")      # every implied slope equal to the height slope
  log_pfvc_fit    <- fit_one("log_pfvc")     # the age-spline-adjusted PFVC slope, for reference
  separated_fit <- keep(list(free_fit, one_column_fit, log_pfvc_fit), separated)
  if (length(separated_fit)) return(list(head = row_head, note = conditionMessage(separated_fit[[1]])))
  b <- coef(free_fit)[B_TERMS]; V <- vcov(free_fit)[B_TERMS, B_TERMS]
  if (anyNA(b)) stop("a group coefficient is not estimable (aliased or empty) for ", cohort_now, ", ", outcome_key)
  # what identifies the height slope: ch_height's SD once the spline, sex and race are removed
  height_identifying_sd <- sd(resid(lm(ch_height ~ ns(age_at_admission, 4) + sex_category + race_category, data = dat)))
  group_counts <- GROUPS %>% mutate(n_group = map2_int(factor_column, level, ~ sum(dat[[.x]] == .y)),
                                    n_group_deaths = map2_int(factor_column, level, ~ sum(dat$event[dat[[.x]] == .y] == 1)))
  list(head = row_head, note = NA_character_,
       theta = as.vector(A %*% b) %>% set_names(THETA), V_theta = A %*% V %*% t(A),
       gamma = b[GROUPS$term], gamma_se = sqrt(diag(V))[GROUPS$term], group_counts = group_counts,
       lr_statistic = 2 * (as.numeric(logLik(free_fit)) - as.numeric(logLik(one_column_fit))),
       one_column = c(estimate = unname(coef(one_column_fit)["ch_size"]), se = unname(sqrt(vcov(one_column_fit)["ch_size", "ch_size"]))),
       log_pfvc = c(estimate = unname(coef(log_pfvc_fit)["log_pfvc"]), se = unname(sqrt(vcov(log_pfvc_fit)["log_pfvc", "log_pfvc"]))),
       height_identifying_sd = height_identifying_sd)
}

# the rows of the slopes table for one set of theta (a cohort, or the difference)
slope_rows <- function(theta, V_theta, quantity) {
  contrast <- as.vector(D %*% theta); contrast_se <- sqrt(diag(D %*% V_theta %*% t(D)))
  bind_rows(
    tibble(quantity = quantity, parameter = "height slope", group = "height",
           log_ratio = theta[["height"]], se = sqrt(V_theta["height", "height"])),
    tibble(quantity = quantity, parameter = "implied slope", group = GROUPS$group,
           log_ratio = unname(theta[GROUPS$group]), se = sqrt(diag(V_theta)[GROUPS$group])),
    tibble(quantity = quantity, parameter = "implied slope minus height slope", group = GROUPS$group,
           log_ratio = contrast, se = contrast_se))
}

slopes <- list(); tests <- list(); vcov_rows <- list(); skipped <- list()
for (k in seq_len(nrow(OUTCOMES))) {
  outcome_key <- OUTCOMES$outcome_key[k]; model <- OUTCOMES$model[k]
  fits <- map(set_names(COHORTS), ~ fit_cohort(.x, outcome_key, model))
  for (cohort_now in COHORTS) {
    fit <- fits[[cohort_now]]
    if (!is.na(fit$note)) { skipped[[length(skipped) + 1]] <- fit$head %>% mutate(note = fit$note); next }
    counts <- fit$group_counts %>% select(group, n_group, n_group_deaths)
    slopes[[length(slopes) + 1]] <- bind_rows(
      slope_rows(fit$theta, fit$V_theta, cohort_now) %>%
        left_join(tibble(group = GROUPS$group, parameter = "implied slope", shift = GROUPS$shift,
                         group_coefficient = unname(fit$gamma), group_coefficient_se = unname(fit$gamma_se)),
                  by = c("group", "parameter")),
      tibble(quantity = cohort_now, group = "height, sex and race together",
             parameter = c("one slope (ch_height + ch_sex + ch_race)", "log PFVC slope (age spline)"),
             log_ratio = c(fit$one_column[["estimate"]], fit$log_pfvc[["estimate"]]),
             se = c(fit$one_column[["se"]], fit$log_pfvc[["se"]]))) %>%
      left_join(counts, by = "group") %>%
      mutate(outcome_key = outcome_key, n_patients = fit$head$n_patients, n_deaths = fit$head$n_deaths,
             height_identifying_sd = fit$height_identifying_sd)
    tests[[length(tests) + 1]] <- bind_rows(
      imap_dfr(TESTS, ~ wald_test(fit$theta, fit$V_theta, .x) %>% mutate(test = .y, method = "Wald")),
      tibble(test = "all: female = Black = Other = height (3 df)", method = "likelihood ratio",
             statistic = fit$lr_statistic, df = nrow(GROUPS), p = pchisq(fit$lr_statistic, nrow(GROUPS), lower.tail = FALSE))) %>%
      mutate(quantity = cohort_now, outcome_key = outcome_key)
    vcov_rows[[length(vcov_rows) + 1]] <- as_tibble(as.table(fit$V_theta), .name_repair = "minimal") %>%
      set_names(c("parameter_row", "parameter_col", "covariance")) %>%
      mutate(quantity = cohort_now, outcome_key = outcome_key, .before = 1)
  }
  # the difference, ventilated minus no support: the demographic paths that act alike
  # in both cohorts cancel (independent cohorts, so the covariances add)
  if (all(map_lgl(fits, ~ is.na(.x$note)))) {
    difference_theta <- fits[["Ventilated"]]$theta - fits[["No support"]]$theta
    difference_V <- fits[["Ventilated"]]$V_theta + fits[["No support"]]$V_theta
    slopes[[length(slopes) + 1]] <- slope_rows(difference_theta, difference_V, DIFFERENCE) %>% mutate(outcome_key = outcome_key)
    tests[[length(tests) + 1]] <- imap_dfr(TESTS, ~ wald_test(difference_theta, difference_V, .x) %>% mutate(test = .y, method = "Wald")) %>%
      mutate(quantity = DIFFERENCE, outcome_key = outcome_key)
    vcov_rows[[length(vcov_rows) + 1]] <- as_tibble(as.table(difference_V), .name_repair = "minimal") %>%
      set_names(c("parameter_row", "parameter_col", "covariance")) %>%
      mutate(quantity = DIFFERENCE, outcome_key = outcome_key, .before = 1)
  }
}
if (length(skipped)) {
  message("\nModels not fitted:")
  print(as.data.frame(bind_rows(skipped)), row.names = FALSE)
}
if (!length(slopes)) stop("no model was fitted: every cohort and outcome was skipped")

slopes <- bind_rows(slopes) %>%
  left_join(OUTCOMES, by = "outcome_key") %>%
  mutate(ratio_type = if_else(model == "logistic", "OR", "HR"),
         ratio_per_0.1 = exp(0.1 * log_ratio), lo_per_0.1 = exp(0.1 * (log_ratio - 1.96 * se)),
         hi_per_0.1 = exp(0.1 * (log_ratio + 1.96 * se)), p = 2 * pnorm(-abs(log_ratio / se)),
         scale = "log ratio per log unit of PFVC; ratio_per_0.1 is per 0.1 log units (about 10% of PFVC)",
         site = site_name) %>%
  select(outcome, ratio_type, quantity, parameter, group, log_ratio, se, ratio_per_0.1, lo_per_0.1, hi_per_0.1, p,
         shift, group_coefficient, group_coefficient_se, n_patients, n_deaths, n_group, n_group_deaths,
         height_identifying_sd, scale, site)
tests <- bind_rows(tests) %>% left_join(OUTCOMES %>% select(outcome_key, outcome), by = "outcome_key") %>%
  mutate(site = site_name) %>% select(outcome, quantity, test, method, statistic, df, p, site)
vcov_rows <- bind_rows(vcov_rows) %>% left_join(OUTCOMES %>% select(outcome_key, outcome), by = "outcome_key") %>%
  mutate(site = site_name) %>% select(outcome, quantity, parameter_row, parameter_col, covariance, site)

message("\nHeight slope and each group's implied slope, per 0.1 log units of PFVC (equality says the group",
        " predicts death as its share of predicted lung size does):")
print(as.data.frame(slopes %>% filter(parameter %in% c("height slope", "implied slope")) %>%
                      select(outcome, quantity, group, ratio_per_0.1, lo_per_0.1, hi_per_0.1, n_group, n_group_deaths) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
message("\nEquality tests:")
print(as.data.frame(tests %>% mutate(across(c(statistic, p), ~ signif(.x, 3)))), row.names = FALSE)

write_csv(slopes, file.path(final_dir, paste0("mortality_channel_equality_slopes_", site_name, ".csv")))
write_csv(tests, file.path(final_dir, paste0("mortality_channel_equality_tests_", site_name, ".csv")))
write_csv(vcov_rows, file.path(final_dir, paste0("mortality_channel_equality_vcov_", site_name, ".csv")))
message("\nWrote mortality_channel_equality_{slopes,tests,vcov}_", site_name, ".csv to ", final_dir)
