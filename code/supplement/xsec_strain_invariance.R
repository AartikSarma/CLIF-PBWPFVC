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
#    for log PBW/PFVC and log PFVC, and their age-standardised versions (GLI at age 25:
#    height, sex and race only, the structural lung size and strain error), per 0.1
#    log units. The difference, ventilated
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
#    both at age 25, the ratio's projection onto the inputs (R2 reported), and GLI's four pieces
#    (height, age, sex, race; pfvc_channels() in 20_biotrauma_grid.R). Read where each
#    reference falls in the placebo cloud, not GLI's percentile alone: the ratio is
#    mostly age, so if age alone lands in the tail the ratio will too.
#    The statistic is |z| of the difference; random directions have no sign.
#    The ratio carries almost no height (the two formulas' height functions nearly
#    cancel), so it cannot compete on the height axis. A second cloud drops log height
#    from the inputs, and the height-free references (the ratio, its projection, the
#    age, sex and race pieces) are ranked in it as well: the strain-error question.
#    The GLI age piece fixes age's shape, so one more reference frees it: the 4-df
#    age spline in each cohort, its four ventilated-minus-control differences tested
#    jointly (Wald, 4 df). It asks whether ventilation changes age's mortality curve in
#    any shape; its p is reported, with the equivalent |z| only for rough placement.
#
# 4. Flexible age with its curve (section 4 below), in everyone and in the
#    intubation-eligible population. Goals of care act on the control only: a
#    ventilated patient's intubation shows intubation was within their goals, while a
#    control under a do-not-intubate order can only die unintubated. That population
#    keeps every ventilated patient and drops the controls with a documented
#    limitation (no record stays in); sections 2 and 4 are repeated in it. Needs the
#    optional clif_code_status table.
#
# The cohorts are xsec_mortality_channel_equality.R's (copied from it, same
# synthetic seed): the ventilated arm on IMV at ICU admission (icu_day0), gated at
# VT/PBW 6-8 and SF < 315 by script 03, and the no-support control.
#
# Ungated (PBWPFVC_VTPBW_GATE=0): the ventilated arm comes from script 03's
# analysis_cross_sectional_ungated, the cohort without the VT/PBW 6-8 gate, and every
# table goes to final/ungated/supplement/ under the same name.
#
# Inputs : intermediate/analysis_cross_sectional.parquet, analysis_all_eligible_timepoints.parquet
#          (script 03, ventilated), intermediate/controls/nosupport/analysis_cross_sectional.parquet
#          and resp_support_waterfall_clean.parquet (scripts 01-03, PBWPFVC_COHORT=nosupport)
#          clif_code_status and clif_hospitalization (config$tables_path; optional)
# Outputs: final/supplement/
#   strain_invariance_dosing_{site}.csv     the dosing rule: slopes, R2, SDs, before the
#                                           VT/PBW gate and (for reference) inside it
#   strain_invariance_contrast_{site}.csv   per cohort, the difference, and the delivered-
#                                           strain row, by population, outcome and dose
#                                           (outcomes include 60-day death with no
#                                           censoring at intubation, and 60-day death or
#                                           intubation, which bracket section 2b)
#   strain_invariance_censoring_{site}.csv  section 2b: each exposure's hazard of
#                                           intubation in the control (informative
#                                           censoring?)
#   strain_invariance_age25_pieces_{site}.csv  section 2c: the age-25 strain error split
#                                           into sex, Black, Other and the height curve,
#                                           per cohort and the difference, each on the
#                                           ratio's scale
#   strain_invariance_placebo_{site}.csv    every index: the two cohorts' coefficients,
#                                           the difference, its z, and (references) the
#                                           percentile of |z| among the placebos; the
#                                           flexible-age tests; the references in the
#                                           intubation-eligible population (column population)
#   strain_invariance_age_curve_{site}.csv  each cohort's age curve and their
#                                           difference, relative to age 60, by population
#   strain_invariance_code_status_{site}.csv  patients and deaths by code status at the
#                                           index, per cohort (with the code_status table)
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
                    "race_category", "pfvc", "pfvc_age25", "pbw", "sf_ratio", "sofa_total", "deceased", "admission_dttm", "death_dttm",
                    "discharge_dttm", "height_cm")
control_file <- file.path(config$output_dir, "controls", "nosupport", "analysis_cross_sectional.parquet")
if (!file.exists(control_file))
  stop("no no-support cohort: run scripts 01-03 with PBWPFVC_COHORT=nosupport first")
# the ventilated arm: patients on invasive ventilation at ICU admission (icu_day0)
# the ungated cohort's table with PBWPFVC_VTPBW_GATE=0 (utils/config.R, script 03)
ventilated_cohort <- read_parquet(file.path(config$output_dir, paste0("analysis_cross_sectional", config$cs_suffix, ".parquet"))) %>%
  select(all_of(cohort_columns), vtpbw, icu_day0, patient_id)
ventilated <- ventilated_cohort %>% filter(icu_day0) %>% select(-icu_day0) %>%
  mutate(cohort = "Ventilated", escalation_dttm = as.POSIXct(NA))
message("Ventilated arm: ", nrow(ventilated), " of ", nrow(ventilated_cohort),
        " ventilated-cohort patients on invasive ventilation at ICU admission (icu_day0)")
no_support <- read_parquet(control_file) %>%
  select(all_of(cohort_columns), escalation_dttm, patient_id) %>%
  mutate(cohort = "No support", vtpbw = NA_real_)
# No patient in both arms. Script 03 removed the gated ventilated arm's patients from the
# control; the ungated arm (PBWPFVC_VTPBW_GATE=0) is larger, so the same rule (by
# patient_id, any hospitalization) is applied again here. On the gated arm it drops no one.
shared_patients <- intersect(no_support$patient_id, ventilated$patient_id)
message("Controls also in the ventilated arm, dropped: ", length(shared_patients))
no_support <- no_support %>% filter(!patient_id %in% shared_patients) %>% select(-patient_id)
ventilated <- ventilated %>% select(-patient_id)
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
  filter(!is.na(pfvc), pfvc > 0, !is.na(pfvc_age25), pfvc_age25 > 0, !is.na(pbw), pbw > 0, !is.na(sf_ratio), !is.na(sofa_total), !is.na(deceased),
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
         # age-standardised (GLI at age 25): the structural lung size and strain error,
         # height, sex and race only
         log_pfvc_age25  = log(pfvc_age25),
         log_ratio_age25 = log(pbw / pfvc_age25),
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
#
# Two companions bracket what censoring a control at intubation can do (section 2b
# asks whether that censoring is informative): every 60-day death, the control's
# deaths after intubation included, and 60-day death or invasive ventilation, in
# which a control's intubation is itself the event, so nothing is censored before day
# 60. In the ventilated arm all three are the same 60-day death.
outcome_data <- function(cohort_data, outcome_key) cohort_data %>% mutate(
  censor_day = if (outcome_key == "day60_before_imv") pmin(HORIZON_DAYS, coalesce(imv_day, Inf)) else HORIZON_DAYS,
  failure_day = if (outcome_key == "day60_death_or_imv") pmin(coalesce(death_index_day, Inf), coalesce(imv_day, Inf))
                else death_index_day,
  failure_day = if_else(is.infinite(failure_day), NA_real_, failure_day),
  event = as.integer(!is.na(failure_day) & failure_day <= censor_day),
  end_day = pmax(if_else(event == 1L, failure_day, censor_day), 0.01))
OUTCOMES <- c(inhosp_logistic = "in-hospital death (logistic)",
              day60_before_imv = "60-day death, before invasive ventilation (Cox)",
              day60_all = "60-day death, all (Cox; control deaths after intubation counted)",
              day60_death_or_imv = "60-day death or invasive ventilation (Cox; control intubation is the event)")

# One coefficient: `term` in the model event ~ SF + SOFA + [extra] + term
fit_term <- function(dat, outcome_key, term, extra = character(0)) {
  dat <- if (outcome_key == "inhosp_logistic") dat %>% mutate(event = deceased) else outcome_data(dat, outcome_key)
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
# Code status, and the intubation-eligible population
# =============================================================================
# Goals of care act on one side of the contrast only. A ventilated patient was
# intubated, which shows that intubation was within their goals at the index, whatever
# the chart records. A control patient under a do-not-intubate order can only die
# unintubated, and such patients are older, so they steepen the control's age and
# PBW/PFVC gradients without any ventilator. The comparison population is therefore
# every ventilated patient against the control patients with no documented limitation
# (a status other than Full or Presume Full, recorded between hospital admission and
# CODE_STATUS_WINDOW_H hours after the index). A control with no record stays in, as
# the default in practice is full code. Restricting the ventilated side instead would
# drop the patients whose limitation is written in the hours after intubation, most of
# whom die, and would select on the outcome.
# Code status is an optional CLIF table (clif_code_status, patient-level, mapped
# through clif_hospitalization, both read from config$tables_path, as
# xsec_pfvc_age_control.R reads them); without it the population is skipped and
# announced.
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
  # the raw tables come from config$tables_path, which a PBWPFVC_SITE_NAME override
  # alone does not change: if they hold none of these hospitalizations they are
  # another site's tables, and the population would be silently wrong
  site_hospitalizations <- read_clif_table("hospitalization", c("patient_id", "hospitalization_id")) %>%
    filter(hospitalization_id %in% both_cohorts$hospitalization_id)
  if (nrow(site_hospitalizations) == 0)
    stop("clif_hospitalization at ", config$tables_path, " holds none of this site's ", nrow(both_cohorts),
         " cohort hospitalizations: config$tables_path points at another site's tables. ",
         "Set PBWPFVC_TABLES_PATH (or config.json) to ", site_name, "'s CLIF tables.")
  # a status counts only from this hospitalization's admission to CODE_STATUS_WINDOW_H
  # after the index (the table is patient-level; orders are often written hours after
  # ICU admission); the last one in that window is the status at the index. Keyed by
  # cohort too: a hospitalization can hold a no-support index and a ventilated one.
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
  code_status_counts <- both_cohorts %>% group_by(cohort, code_status_at_index) %>%
    summarise(n_patients = n(), n_deaths = sum(deceased == 1), .groups = "drop") %>%
    mutate(in_eligible_population = cohort == "Ventilated" | code_status_at_index != "limited or other", site = site_name)
  message("\nCode status at the index (last status up to ", CODE_STATUS_WINDOW_H,
          " h after it); the intubation-eligible population drops only controls with a documented limitation:")
  print(as.data.frame(code_status_counts), row.names = FALSE)
} else {
  message("\n*** No code_status table at ", config$tables_path, ": the intubation-eligible population is skipped. ***")
  both_cohorts <- both_cohorts %>% mutate(code_status_at_index = "no record")
  code_status_counts <- NULL
}
both_cohorts <- both_cohorts %>% mutate(intubation_eligible = cohort == "Ventilated" | code_status_at_index != "limited or other")
ELIGIBLE_LABEL <- "intubation-eligible (every ventilated patient; controls without a documented limitation)"

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
POPULATIONS <- c(everyone = "everyone", hypoxemic = "hypoxemic at the index (SF < 315)",
                 if (HAS_CODE_STATUS) c(eligible = ELIGIBLE_LABEL))
population_data <- function(population, cohort_now) both_cohorts %>%
  filter(cohort == cohort_now,
         if (population == "hypoxemic") hypoxemic else if (population == "eligible") intubation_eligible else TRUE)
EXPOSURES <- c(log_ratio = "log PBW/PFVC", log_pfvc = "log PFVC",
               log_ratio_age25 = "log PBW/PFVC at age 25", log_pfvc_age25 = "log PFVC at age 25")
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
# 2b. Is censoring the control at intubation informative?
# =============================================================================
# "60-day death before invasive ventilation" censors a control at intubation. If an
# exposure predicts which controls are intubated, the censoring is informative: the
# control's coefficient is estimated on the patients the exposure kept unventilated,
# and the difference shifts with it. Per control population, the cause-specific
# hazard of invasive ventilation within HORIZON_DAYS (death is the competing event,
# censored), on SF, SOFA and each exposure, as the mortality models. A hazard ratio
# near 1 means the censoring does not select on the exposure; the two companion
# outcomes of section 2 show how far the contrast moves when nothing is censored.
censoring_check <- expand_grid(population = names(POPULATIONS), exposure = names(EXPOSURES)) %>%
  pmap_dfr(function(population, exposure) {
    dat <- population_data(population, "No support") %>%
      mutate(event = as.integer(!is.na(imv_day) & imv_day <= HORIZON_DAYS &
                                  (is.na(death_index_day) | imv_day <= death_index_day)),
             end_day = pmax(pmin(HORIZON_DAYS, coalesce(imv_day, Inf), coalesce(death_index_day, Inf)), 0.01))
    head_row <- tibble(population = POPULATIONS[[population]], exposure = EXPOSURES[[exposure]],
                       n_patients = nrow(dat), n_intubated = sum(dat$event))
    if (head_row$n_intubated < MIN_EVENTS)
      return(head_row %>% mutate(log_hr = NA_real_, se = NA_real_, note = paste("skipped: fewer than", MIN_EVENTS, "intubations")))
    term <- exposure   # the column name; inside mutate(), `exposure` is the label column
    fit <- fit_or_separation(fit_strict(coxph(as.formula(paste("Surv(end_day, event) ~ sf_z + sofa_z +", term)), data = dat)))
    if (separated(fit)) return(head_row %>% mutate(log_hr = NA_real_, se = NA_real_, note = conditionMessage(fit)))
    head_row %>% mutate(log_hr = unname(coef(fit)[term]), se = unname(sqrt(vcov(fit)[term, term])), note = NA_character_)
  }) %>%
  mutate(hr_per_0.1 = exp(0.1 * log_hr), lo_per_0.1 = exp(0.1 * (log_hr - 1.96 * se)), hi_per_0.1 = exp(0.1 * (log_hr + 1.96 * se)),
         p = 2 * pnorm(-abs(log_hr / se)), outcome = paste0("invasive ventilation within ", HORIZON_DAYS, " days, control only (cause-specific Cox)"),
         scale = "log HR per log unit of the exposure; hr_per_0.1 is per 0.1 log units", site = site_name)
message("\nIs the control's censoring at intubation informative? Hazard of intubation per 0.1 log units:")
print(as.data.frame(censoring_check %>% select(population, exposure, n_patients, n_intubated, hr_per_0.1, lo_per_0.1, hi_per_0.1, p) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)

# =============================================================================
# 2c. The age-free strain error by its inputs: sex, race and height
# =============================================================================
# log PBW/PFVC at age 25 is Devine's PBW over GLI's FVC at age 25: a fixed function of
# sex, race and height. Its ventilator-specific contrast is split here by input, so
# that one input (at MIMIC, the OTHER category, which holds the unknown race) cannot
# carry it unseen. Per cohort, with no age term (the age-free index has none):
#   death ~ SF + SOFA + female + Black + Other + height piece
# The height piece is the ratio's own height curve at each patient's sex
# (height_fingerprint(), 20_biotrauma_grid.R), in log PBW/PFVC units. Each indicator
# is also put on the ratio's scale: its coefficient over the ratio-at-25 shift that
# GLI and Devine assign it at the cohort's median height (implied_per_0.1). If the
# contrast is strain error, the pieces' implied slopes agree with one another and
# with the one-beta contrast of section 2; a piece that runs alone is that input's
# own ventilator-specific effect, not strain.
ratio25_at <- function(height_cm, sex, race) {
  sex_code <- if (sex == "Female") 2L else 1L
  race_code <- switch(race, WHITE = 1L, BLACK = 2L, OTHER = 5L)
  devine <- if (sex == "Female") 45.5 + 2.3 * (height_cm / 2.54 - 60) else 50 + 2.3 * (height_cm / 2.54 - 60)
  log(devine) - log(rspiro::pred_GLI(age = 25, height = height_cm / 100, gender = sex_code, ethnicity = race_code, param = "FVC"))
}
median_height <- median(both_cohorts$height_cm)
PIECE_SHIFTS <- c(female = ratio25_at(median_height, "Female", "WHITE") - ratio25_at(median_height, "Male", "WHITE"),
                  black  = ratio25_at(median_height, "Male", "BLACK") - ratio25_at(median_height, "Male", "WHITE"),
                  other  = ratio25_at(median_height, "Male", "OTHER") - ratio25_at(median_height, "Male", "WHITE"),
                  height_piece = 1)   # already in log PBW/PFVC units
PIECE_LABELS <- c(female = "sex (female)", black = "race (Black)", other = "race (Other)", height_piece = "height (own-sex curve)")
both_cohorts <- both_cohorts %>%
  mutate(female = as.numeric(sex_category == "Female"), black = as.numeric(race_category == "BLACK"),
         other = as.numeric(race_category == "OTHER"),
         height_piece = height_fingerprint(height_cm, as.character(sex_category)))
fit_pieces <- function(dat, outcome_key) {
  dat <- if (outcome_key == "inhosp_logistic") dat %>% mutate(event = deceased) else outcome_data(dat, outcome_key)
  pieces <- names(PIECE_SHIFTS)
  head_row <- tibble(n_patients = nrow(dat), n_deaths = sum(dat$event == 1))
  if (head_row$n_deaths < MIN_EVENTS)
    return(head_row %>% mutate(piece = pieces, log_ratio = NA_real_, se = NA_real_, note = paste("skipped: fewer than", MIN_EVENTS, "deaths")))
  rhs <- paste(c("sf_z", "sofa_z", pieces), collapse = " + ")
  fit <- fit_or_separation(if (outcome_key == "inhosp_logistic")
    fit_strict(glm(as.formula(paste("event ~", rhs)), family = binomial, data = dat)) else
    fit_strict(coxph(as.formula(paste("Surv(end_day, event) ~", rhs)), data = dat)))
  if (separated(fit)) return(head_row %>% mutate(piece = pieces, log_ratio = NA_real_, se = NA_real_, note = conditionMessage(fit)))
  head_row %>% slice(rep(1, length(pieces))) %>%
    mutate(piece = pieces, log_ratio = unname(coef(fit)[pieces]), se = unname(sqrt(diag(vcov(fit))[pieces])), note = NA_character_)
}
age25_pieces <- expand_grid(population = names(POPULATIONS), outcome_key = names(OUTCOMES)) %>%
  mutate(fits = map2(population, outcome_key, function(population, outcome_key) {
    ventilated <- fit_pieces(population_data(population, "Ventilated"), outcome_key) %>% mutate(quantity = "Ventilated")
    control <- fit_pieces(population_data(population, "No support"), outcome_key) %>% mutate(quantity = "No support")
    difference <- ventilated %>% transmute(piece, quantity = DIFFERENCE, log_ratio = log_ratio - control$log_ratio,
                                           se = sqrt(se^2 + control$se^2), n_patients = NA_integer_, n_deaths = NA_integer_,
                                           note = coalesce(note, control$note))
    bind_rows(ventilated, control, difference)
  })) %>% unnest(fits) %>%
  mutate(population = POPULATIONS[population], outcome = OUTCOMES[outcome_key],
         shift_in_log_ratio25 = PIECE_SHIFTS[piece], input = PIECE_LABELS[piece],
         implied_log_ratio = log_ratio / shift_in_log_ratio25, implied_se = se / abs(shift_in_log_ratio25),
         implied_per_0.1 = exp(0.1 * implied_log_ratio),
         implied_lo_per_0.1 = exp(0.1 * (implied_log_ratio - 1.96 * implied_se)),
         implied_hi_per_0.1 = exp(0.1 * (implied_log_ratio + 1.96 * implied_se)),
         p = 2 * pnorm(-abs(log_ratio / se)),
         scale = "log_ratio: log OR or log HR per unit of the piece (an indicator, or log PBW/PFVC units for height); implied_*: per 0.1 log units of PBW/PFVC at age 25",
         site = site_name) %>%
  select(population, outcome, quantity, input, piece, log_ratio, se, p, shift_in_log_ratio25,
         implied_per_0.1, implied_lo_per_0.1, implied_hi_per_0.1, n_patients, n_deaths, note, scale, site)
message("\nThe age-25 strain error by input, ventilated minus control, everyone (implied per 0.1 log units of PBW/PFVC at 25):")
print(as.data.frame(age25_pieces %>% filter(population == "everyone", quantity == DIFFERENCE) %>%
                      select(outcome, input, implied_per_0.1, implied_lo_per_0.1, implied_hi_per_0.1, p) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)

# =============================================================================
# 3. Placebo formulas
# =============================================================================
# The inputs, whitened over both cohorts together: X centred, then rotated and scaled
# by its covariance's eigen-decomposition, so Z has identity covariance and any unit
# vector w gives an index Z w with SD 1.
# the age spline's basis, kept so the flexible-age curve can be evaluated at any age
age_basis <- ns(both_cohorts$age_at_admission, df = 4)
input_matrix <- cbind(age_basis,
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
  `log PBW/PFVC at age 25`            = both_cohorts$log_ratio_age25,
  `log PFVC at age 25`                  = both_cohorts$log_pfvc_age25,
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
# A second cloud without height. log PBW - log PFVC nearly cancels height (Devine's
# and GLI's height functions almost coincide), so the ratio has almost no height in it,
# and it does not compete on the height axis, where lung size carries a strong
# ventilation-specific signal. The strain-error question is whether the ratio stands
# out among indices of the same height-free inputs (age spline, sex, race), so the
# whitening, the directions and the ranking are repeated without log height, for the
# references that carry no height by construction.
height_free_matrix <- input_matrix[, colnames(input_matrix) != "log_height"]
height_free_centred <- scale(height_free_matrix, center = TRUE, scale = FALSE)
eigen_height_free <- eigen(cov(height_free_centred), symmetric = TRUE)
whitened_height_free <- height_free_centred %*% eigen_height_free$vectors %*% diag(1 / sqrt(eigen_height_free$values))
set.seed(PLACEBO_SEED + 1)
height_free_directions <- matrix(rnorm(PLACEBO_N * ncol(whitened_height_free)), nrow = ncol(whitened_height_free))
height_free_directions <- sweep(height_free_directions, 2, sqrt(colSums(height_free_directions^2)), "/")
height_free_indices <- map(seq_len(PLACEBO_N), ~ as.vector(whitened_height_free %*% height_free_directions[, .x]))
names(height_free_indices) <- sprintf("height-free placebo %03d", seq_len(PLACEBO_N))
HEIGHT_FREE_REFERENCES <- c("log PBW/PFVC (GLI and Devine)", "log PBW/PFVC projected on the inputs", "log PBW/PFVC at age 25",
                            "GLI age piece", "GLI sex piece", "GLI race piece")
message("\nThe ratio's height content: correlation of log PBW/PFVC with log height ",
        signif(cor(both_cohorts$log_ratio, log(both_cohorts$height_cm)), 3), "; R2 of log PBW/PFVC on the height-free inputs ",
        signif(summary(lm(both_cohorts$log_ratio ~ height_free_matrix))$r.squared, 3), " and on all inputs ",
        signif(summary(lm(both_cohorts$log_ratio ~ input_matrix))$r.squared, 3))
# the in-hospital contrast (no dose) of one index, in the patients `rows` selects
index_contrast <- function(index_values, rows = rep(TRUE, nrow(both_cohorts))) {
  dat <- both_cohorts %>% mutate(index = index_values) %>% filter(rows)
  fits <- map(set_names(COHORTS), ~ fit_term(dat %>% filter(cohort == .x), "inhosp_logistic", "index"))
  tibble(ventilated = fits$Ventilated$log_ratio, ventilated_se = fits$Ventilated$se,
         control = fits$`No support`$log_ratio, control_se = fits$`No support`$se,
         note = coalesce(fits$Ventilated$note, fits$`No support`$note),
         n_ventilated = fits$Ventilated$n_patients, deaths_ventilated = fits$Ventilated$n_deaths,
         n_control = fits$`No support`$n_patients, deaths_control = fits$`No support`$n_deaths) %>%
    mutate(difference = ventilated - control, difference_se = sqrt(ventilated_se^2 + control_se^2),
           z = difference / difference_se)
}
# The placebos and references are run in everyone and, with the code_status table, in
# the intubation-eligible population (section "Code status"); each reference is
# placed among its own population's placebos. The indices keep everyone's
# standardisation in both.
placebo_populations <- c(everyone = "everyone", if (HAS_CODE_STATUS) c(eligible = ELIGIBLE_LABEL))
placebo_abs_z <- list()
placebo_results <- imap_dfr(placebo_populations, function(label, key) {
  rows <- if (key == "eligible") both_cohorts$intubation_eligible else rep(TRUE, nrow(both_cohorts))
  message("\nPlacebo formulas, ", label, ": ", PLACEBO_N, " random directions (all inputs), ", PLACEBO_N,
          " without height, and ", length(reference_indices), " references ...")
  reference_rows <- imap_dfr(reference_indices, ~ index_contrast(.x, rows) %>% mutate(index = .y))
  cloud_rows <- imap_dfr(placebo_indices, ~ index_contrast(.x, rows) %>% mutate(index = .y))
  height_free_rows <- imap_dfr(height_free_indices, ~ index_contrast(.x, rows) %>% mutate(index = .y))
  abs_z <- abs(cloud_rows$z); abs_z_height_free <- abs(height_free_rows$z)
  placebo_abs_z[[label]] <<- abs_z
  results <- bind_rows(
    reference_rows %>% mutate(kind = "reference", cloud = "all inputs",
                              percentile_abs_z_among_placebos = map_dbl(abs(z), ~ mean(abs_z < .x, na.rm = TRUE))),
    reference_rows %>% filter(index %in% HEIGHT_FREE_REFERENCES) %>%
      mutate(kind = "reference", cloud = "without height",
             percentile_abs_z_among_placebos = map_dbl(abs(z), ~ mean(abs_z_height_free < .x, na.rm = TRUE))),
    cloud_rows %>% mutate(kind = "placebo", cloud = "all inputs"),
    height_free_rows %>% mutate(kind = "placebo", cloud = "without height")) %>%
    mutate(population = label)
  message("Reference indices against the placebos (in-hospital, severity only, per SD of the index; the ratio's",
          " projection on the inputs has R2 ", signif(summary(ratio_projection)$r.squared, 3), "):")
  print(as.data.frame(results %>% filter(kind == "reference") %>%
                        select(index, cloud, ventilated, control, difference, difference_se, z, percentile_abs_z_among_placebos,
                               n_ventilated, deaths_ventilated, n_control, deaths_control) %>%
                        mutate(across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)
  message("Placebo |z|, all inputs: median ", signif(median(abs_z, na.rm = TRUE), 3), ", 95th percentile ",
          signif(quantile(abs_z, 0.95, na.rm = TRUE), 3), "; without height: median ",
          signif(median(abs_z_height_free, na.rm = TRUE), 3), ", 95th percentile ", signif(quantile(abs_z_height_free, 0.95, na.rm = TRUE), 3))
  results
}) %>%
  mutate(ratio_projection_r2 = summary(ratio_projection)$r.squared,
         scale = "log OR of in-hospital death per SD of the index (SD over both cohorts, everyone); severity only",
         site = site_name) %>%
  relocate(kind, index, population, cloud)

# =============================================================================
# 4. Flexible age, in everyone and in the intubation-eligible population
# =============================================================================
# The GLI age piece gives age one shape (GLI's curve) and one coefficient; ventilation
# could change age's mortality curve in another shape (flat, then steep in the
# oldest). Here age enters each cohort as the placebo inputs' 4-df natural spline
# (knots from both cohorts together, so the two fits share one basis), and the four
# ventilated-minus-control differences are tested jointly (Wald, 4 df): does
# ventilation change age's curve in ANY shape? Its statistic is a chi-square on 4 df,
# not a z on 1 df; the row carries the equivalent |z| (the normal quantile of its p)
# and, in everyone, that |z|'s percentile among the placebos, only as a rough
# placement. The curve behind the test is written too: each cohort's log-odds of
# in-hospital death by age relative to age CURVE_REFERENCE_AGE, and their difference,
# at the cohort's mean severity.
#
# Goals of care can reshape the control's age curve without any ventilator (the
# section "Code status" above), so the test, the curve and the named reference
# indices are repeated in the intubation-eligible population: every ventilated
# patient, and the controls without a documented limitation. A difference that
# collapses there is do-not-intubate selection, not ventilation.
CURVE_AGES <- seq(30, 90, by = 5)
CURVE_REFERENCE_AGE <- 60
age_spline_terms <- paste0("age_spline_", 1:4)
age_spline_data <- both_cohorts %>% bind_cols(as_tibble(input_matrix[, age_spline_terms]))
# the spline's rows at each curve age, less its row at the reference age
curve_contrast_rows <- predict(age_basis, CURVE_AGES) -
  predict(age_basis, rep(CURVE_REFERENCE_AGE, length(CURVE_AGES)))
colnames(curve_contrast_rows) <- age_spline_terms
flexible_age_for <- function(rows, population) {
  fits <- map(set_names(COHORTS), function(cohort_now) {
    dat <- age_spline_data %>% filter(rows, cohort == cohort_now)
    if (sum(dat$deceased == 1) < MIN_EVENTS) return(NULL)
    fit <- fit_strict(glm(as.formula(paste("deceased ~ sf_z + sofa_z +", paste(age_spline_terms, collapse = " + "))),
                          family = binomial, data = dat))
    list(b = coef(fit)[age_spline_terms], V = vcov(fit)[age_spline_terms, age_spline_terms],
         n_patients = nrow(dat), n_deaths = sum(dat$deceased == 1))
  })
  if (any(map_lgl(fits, is.null))) {
    message("  flexible age not estimable in ", population, ": fewer than ", MIN_EVENTS, " deaths in a cohort")
    return(NULL)
  }
  difference <- fits$Ventilated$b - fits$`No support`$b
  difference_V <- fits$Ventilated$V + fits$`No support`$V          # independent cohorts
  chi2 <- as.numeric(t(difference) %*% solve(difference_V) %*% difference)
  p <- pchisq(chi2, df = length(age_spline_terms), lower.tail = FALSE)
  test_row <- tibble(kind = "reference", index = "age, 4-df spline (any shape; Wald on 4 df)", population = population,
                     chi2 = chi2, df = length(age_spline_terms), p = p, z = qnorm(p / 2, lower.tail = FALSE),
                     n_ventilated = fits$Ventilated$n_patients, deaths_ventilated = fits$Ventilated$n_deaths,
                     n_control = fits$`No support`$n_patients, deaths_control = fits$`No support`$n_deaths)
  curve_for <- function(b, V, quantity) tibble(
    population = population, quantity = quantity, age = CURVE_AGES, reference_age = CURVE_REFERENCE_AGE,
    log_odds = as.vector(curve_contrast_rows %*% b),
    se = sqrt(rowSums((curve_contrast_rows %*% V) * curve_contrast_rows)))
  curve <- bind_rows(curve_for(fits$Ventilated$b, fits$Ventilated$V, "Ventilated"),
                     curve_for(fits$`No support`$b, fits$`No support`$V, "No support"),
                     curve_for(difference, difference_V, DIFFERENCE))
  list(test = test_row, curve = curve)
}

age_populations <- c(everyone = "everyone", if (HAS_CODE_STATUS) c(eligible = ELIGIBLE_LABEL))
flexible_age <- imap(age_populations, function(label, key)
  flexible_age_for(if (key == "eligible") age_spline_data$intubation_eligible else rep(TRUE, nrow(age_spline_data)), label)) %>%
  compact()
age_tests <- map_dfr(flexible_age, "test") %>%
  mutate(percentile_abs_z_among_placebos = map2_dbl(population, z, ~ mean(placebo_abs_z[[.x]] < .y, na.rm = TRUE)),
         scale = "Wald chi-square (4 df) on the ventilated-minus-control differences of the age spline; z is the normal equivalent of p")
age_curve <- map_dfr(flexible_age, "curve") %>%
  mutate(lo = log_odds - 1.96 * se, hi = log_odds + 1.96 * se,
         scale = paste0("log-odds of in-hospital death relative to age ", CURVE_REFERENCE_AGE, ", at the cohort's mean severity"),
         site = site_name)
message("\nFlexible age (4-df spline, ventilated minus control, joint Wald on 4 df):")
print(as.data.frame(age_tests %>% select(population, chi2, p, z, percentile_abs_z_among_placebos, n_ventilated, deaths_ventilated,
                                         n_control, deaths_control) %>%
                      mutate(across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)
message("Ventilated minus control, log-odds by age relative to ", CURVE_REFERENCE_AGE, ":")
print(as.data.frame(age_curve %>% filter(quantity == DIFFERENCE) %>%
                      transmute(population, age, difference = signif(log_odds, 3), lo = signif(lo, 3), hi = signif(hi, 3)) %>%
                      pivot_wider(names_from = population, values_from = c(difference, lo, hi))), row.names = FALSE)

placebo_results <- bind_rows(placebo_results, age_tests) %>%
  mutate(ratio_projection_r2 = summary(ratio_projection)$r.squared, site = site_name) %>%
  relocate(kind, index, population)

write_csv(dosing, file.path(final_dir, paste0("strain_invariance_dosing_", site_name, ".csv")))
write_csv(contrast, file.path(final_dir, paste0("strain_invariance_contrast_", site_name, ".csv")))
write_csv(censoring_check, file.path(final_dir, paste0("strain_invariance_censoring_", site_name, ".csv")))
write_csv(age25_pieces, file.path(final_dir, paste0("strain_invariance_age25_pieces_", site_name, ".csv")))
write_csv(placebo_results, file.path(final_dir, paste0("strain_invariance_placebo_", site_name, ".csv")))
write_csv(age_curve, file.path(final_dir, paste0("strain_invariance_age_curve_", site_name, ".csv")))
if (!is.null(code_status_counts))
  write_csv(code_status_counts, file.path(final_dir, paste0("strain_invariance_code_status_", site_name, ".csv")))
message("\nWrote strain_invariance_{dosing,contrast,placebo,age_curve", if (!is.null(code_status_counts)) ",code_status" else "",
        "}_", site_name, ".csv to ", final_dir)
