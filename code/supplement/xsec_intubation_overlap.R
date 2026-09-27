# =============================================================================
# Supplement (cross-sectional): could a propensity score for intubation weight the
# no-support controls to the ventilated arm? An overlap diagnostic.
# =============================================================================
# Figure 4 compares the lung-size divergence of ventilated patients with that of
# no-support controls, reading the controls at the ventilated severity through a
# linear severity-anchor term. A propensity score for intubation could instead
# weight the controls to the ventilated arm, but only where the arms overlap. This
# script fits that score and reports overlap, balance and effective sample size. It
# fits no outcome model.
#
# Populations: two comparisons, each "ventilated at ICU admission" against a control.
#   ventilated arm   the ventilated cohort's patients on invasive ventilation at ICU
#                    admission (script 03's icu_day0)
#   control A        the no-support cohort, every patient
#   control B        the no-support cohort with index SF < 315 (the hypoxemic control)
# Script 03 has already removed from the control every patient of the ventilated
# arm, so no patient is in both.
#
# Time zero (the treatment decision) is the same event in both arms: the first ICU
# admission of the stay (the earliest in_dttm in cohort_icu_stays). Covariates come
# from [T0 - PRE_WINDOW_H, T0). A ventilated patient intubated before ICU admission
# (in the ED or the operating room) has every measurement taken at or after the first
# IMV record removed from that window, so no covariate is recorded on the ventilator.
# The first IMV record is the first respiratory-support row with device imv or a set
# tidal volume above zero; the same cut is applied to the control, which has no IMV
# before T0 by construction. The patients intubated before ICU admission, and those
# who lost a measurement to the cut, are counted per arm. A control's index is its
# first qualifying row within 6 hours of an ICU admission (script 03); a control
# indexed after a later ICU admission is still read at its first one, and the count
# of those is printed and written.
#
# Covariates (PS_COVARIATES below; edit there):
#   structural     ns(age, 3), sex, race (PFVC's inputs besides height), ns(height, 3)
#                  (so the arms' PFVC distributions are balanced); PS_INCLUDE_STRUCTURAL
#                  = FALSE drops all four from the score (their balance is still
#                  reported)
#   oxygenation    lowest SpO2; highest FiO2; lowest SF ratio; highest support level
#                  (room air < nasal cannula < face mask or other < high-flow nasal
#                  cannula < NIPPV or CPAP)
#   circulation    lowest MAP; any vasopressor; highest NE-equivalent dose
#   neurological   lowest GCS
#   labs           lowest arterial pH, highest PaCO2, lowest PaO2/FiO2, highest
#                  creatinine, highest bilirubin, lowest platelets
#   timing         hours from hospital admission to T0 (ICU admission)
# FiO2 on room air and nasal cannula is estimated by one rule in both arms: the
# documented fio2_set where present, else 0.21 on room air and 0.21 + 0.03 per L/min
# on a cannula capped at 0.60 (utils/config.R, estimate_fio2_nosupport, the rule the
# control's own index SF uses). The ventilated cohort's waterfall is not estimated by
# script 03, so without the shared rule its pre-intubation SF would be missing where
# the control's is not, and missingness would predict intubation by measurement alone.
# SpO2 and PaO2 each take the most recent FiO2 recorded up to FIO2_LOOKBACK_H hours
# before them, as in script 03, and the SF ratio uses SpO2 80-97 only (03's range);
# the lowest SpO2 uses every SpO2. The NE-equivalent dose is script 03's series in
# force (ne_equiv_admin), read over the window including the dose already running
# when it opens; a patient with no dose in the window has 0 and no vasopressor.
# Labs are timed by result (lab_result_dttm), as in scripts 03 and 10, so a lab
# resulted before T0 was drawn before it.
# Missing values: each covariate with any missing value in a comparison gets a
# missing indicator, and the value is set to that comparison's pooled median (the
# median level for the support level). A covariate never observed, or constant, in a
# comparison cannot enter its score: it is left out and listed in the summary
# (covariates_not_in_score). A continuous covariate whose spline knots (its 0, 1/3,
# 2/3 and 1 quantiles, after the median fill) are not all distinct enters linearly
# (covariates_linear): sparse data filled at the median pile up there. The share missing per arm is
# written: pre-ICU data are often sparse, and missingness itself predicts intubation.
# Continuous covariates enter as ns(x, 3); binary and categorical ones as factors.
#
# Model: logistic regression, P(ventilated arm | covariates), fitted separately for
# comparisons A and B. The control is built with no high-flow, NIPPV, CPAP or IMV
# before its index, so the top support levels can occur only in the ventilated arm:
# the fit may separate. Separation is the overlap finding, not a failure: its warning
# is recorded in the summary's note and every diagnostic is still computed (a score
# of 1 is a patient no control resembles). Any other warning stops the script.
#
# Diagnostics, aggregates only:
#   AUC (c-statistic) and patients per arm
#   the score's distribution per arm: counts in bins of width PS_BIN_WIDTH, and its
#   minimum, median and maximum
#   common support: [larger of the two minima, smaller of the two maxima] and the
#   share of each arm inside it
#   weights: ATT odds (ventilated 1, control ps / (1 - ps)) and overlap (ATO:
#   ventilated 1 - ps, control ps); per arm, Kish effective sample size
#   (sum w)^2 / sum w^2, the 50th, 90th and 99th percentiles and the maximum of the
#   weights, and the share of the arm's total weight carried by its top 1% of patients
#   balance: standardised mean difference per covariate term (each level of a
#   categorical, each missing indicator), unweighted, ATT-weighted and
#   overlap-weighted, all over the unweighted pooled SD sqrt((s1^2 + s0^2) / 2);
#   |SMD| > SMD_FLAG flagged
#   context: the mean log PFVC (the exposure) and the mean figure-4 severity anchor
#   (SOFA cardiovascular + coagulation + liver + renal, from analysis_cross_sectional)
#   per arm and weighting, beside the current standardisation
#
# Inputs (patient-level, never shared), for the ventilated cohort
# (intermediate/) and the no-support cohort (intermediate/controls/nosupport/):
#   analysis_cross_sectional, cohort_icu_stays, resp_support_waterfall_clean,
#   cohort_vitals_clean, cohort_labs_clean, cohort_assessments, ne_equiv_admin
# Outputs: final/supplement/
#   intubation_overlap_summary_{site}.csv   per comparison, weighting and arm: n, AUC,
#                                           common support, score range, ESS, weight
#                                           quantiles, top-1% share, mean log PFVC and
#                                           anchor, the separation note, the covariates
#                                           left out or entered linearly
#   intubation_overlap_balance_{site}.csv   per comparison, covariate term and
#                                           weighting: SMD and the |SMD| > 0.1 flag
#   intubation_overlap_ps_{site}.csv        the binned score per comparison and arm
#   intubation_overlap_missing_{site}.csv   share missing per covariate, comparison
#                                           and arm
#   intubation_overlap_time_zero_{site}.csv per arm: patients intubated before ICU
#                                           admission, patients with a measurement
#                                           removed by the IMV cut, controls indexed
#                                           more than 6 h after their first ICU admission
#   intubation_overlap_{site}.pdf           per comparison: the score mirrored by arm
#                                           (ventilated up, control down) and a love plot
# Usage: uvr run code/supplement/xsec_intubation_overlap.R   (PBWPFVC_COHORT unset)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(splines)
  library(data.table)
  library(patchwork)
})

options(width = 220)
source("utils/config.R")
if (config$cohort != "imv") stop("xsec_intubation_overlap.R reads both cohorts itself: unset PBWPFVC_COHORT")
site_name <- config$site_name
final_dir <- final_dir_for("supplement")

PRE_WINDOW_H          <- 24     # the pre-decision window, hours before T0
PS_INCLUDE_STRUCTURAL <- TRUE   # FALSE drops age, sex, race and height from the score
SF_HYPOXEMIA_THRESHOLD <- 315   # script 03's index gate for the ventilated cohort; control B
SPO2_SF_RANGE         <- c(80, 97)   # script 03's SpO2 range for the SF ratio
CONTROL_ICU_WINDOW_H  <- 6      # script 03: the control's index lies within this long of an ICU admission
PS_BIN_WIDTH          <- 0.02
SMD_FLAG              <- 0.1
SEPARATION_PATTERN    <- "fitted probabilities numerically 0 or 1|coefficient may be infinite"

# The score's covariates. type: continuous (ns(x, 3)), binary or categorical.
# structural: dropped from the score when PS_INCLUDE_STRUCTURAL is FALSE.
PS_COVARIATES <- tribble(
  ~covariate,              ~label,                                    ~type,         ~structural,
  "age_at_admission",      "Age",                                     "continuous",  TRUE,
  "sex_category",          "Sex",                                     "categorical", TRUE,
  "race_category",         "Race",                                    "categorical", TRUE,
  "height_cm",             "Height",                                  "continuous",  TRUE,
  "lowest_spo2",           "Lowest SpO2",                             "continuous",  FALSE,
  "highest_fio2",          "Highest FiO2",                            "continuous",  FALSE,
  "lowest_sf",             "Lowest SF ratio",                         "continuous",  FALSE,
  "support_level",         "Highest support level",                   "categorical", FALSE,
  "lowest_map",            "Lowest MAP",                              "continuous",  FALSE,
  "any_vasopressor",       "Any vasopressor",                         "binary",      FALSE,
  "highest_ne_equiv",      "Highest NE-equivalent dose",              "continuous",  FALSE,
  "lowest_gcs",            "Lowest GCS",                              "continuous",  FALSE,
  "lowest_ph",             "Lowest arterial pH",                      "continuous",  FALSE,
  "highest_paco2",         "Highest PaCO2",                           "continuous",  FALSE,
  "lowest_pf",             "Lowest PaO2/FiO2",                        "continuous",  FALSE,
  "highest_creatinine",    "Highest creatinine",                      "continuous",  FALSE,
  "highest_bilirubin",     "Highest bilirubin",                       "continuous",  FALSE,
  "lowest_platelets",      "Lowest platelets",                        "continuous",  FALSE,
  "hours_admission_to_t0", "Hours from hospital admission to ICU",    "continuous",  FALSE)
SUPPORT_LEVELS <- c("room air", "nasal cannula", "face mask or other", "high-flow nasal cannula", "NIPPV or CPAP")
SUPPORT_LEVEL_OF_DEVICE <- c("room air" = "room air", "nasal cannula" = "nasal cannula",
                             "face mask" = "face mask or other", "trach collar" = "face mask or other",
                             "other" = "face mask or other", "high flow nc" = "high-flow nasal cannula",
                             "nippv" = "NIPPV or CPAP", "cpap" = "NIPPV or CPAP")
ARMS <- c("Ventilated", "Control")
ARM_COLOURS <- c(Ventilated = "#0072B2", Control = "#E69F00")    # Okabe-Ito
WEIGHTINGS <- c(unweighted = "unweighted", att = "ATT odds", overlap = "overlap (ATO)")
WEIGHTING_COLOURS <- c(unweighted = "#D55E00", att = "#009E73", overlap = "#CC79A7")   # Okabe-Ito

# =============================================================================
# Data: both cohorts' tables
# =============================================================================
ventilated_dir <- config$output_dir
control_dir    <- file.path(config$output_dir, "controls", "nosupport")
INPUT_TABLES <- c("analysis_cross_sectional", "cohort_icu_stays", "resp_support_waterfall_clean",
                  "cohort_vitals_clean", "cohort_labs_clean", "cohort_assessments", "ne_equiv_admin")
for (cohort_dir in c(ventilated_dir, control_dir)) {
  missing_tables <- INPUT_TABLES[!file.exists(file.path(cohort_dir, paste0(INPUT_TABLES, ".parquet")))]
  if (length(missing_tables))
    stop("missing in ", cohort_dir, ": ", paste(missing_tables, collapse = ", "),
         ". Run scripts 01-03 for the ventilated cohort and then with PBWPFVC_COHORT=nosupport.")
}
read_table <- function(cohort_dir, table_name, columns) read_parquet(file.path(cohort_dir, paste0(table_name, ".parquet")),
                                                                      col_select = all_of(columns))

cross_sectional_columns <- c("hospitalization_id", "admission_dttm", "index_dttm", "age_at_admission", "sex_category",
                             "race_category", "height_cm", "pfvc", "sf_ratio",
                             "sofa_cv_97", "sofa_coag", "sofa_liver", "sofa_renal")
ventilated_cohort <- read_table(ventilated_dir, "analysis_cross_sectional", c(cross_sectional_columns, "icu_day0"))
ventilated_arm <- ventilated_cohort %>% filter(icu_day0) %>% select(-icu_day0) %>% mutate(arm = "Ventilated")
message("Ventilated arm: ", nrow(ventilated_arm), " of ", nrow(ventilated_cohort),
        " ventilated-cohort patients on invasive ventilation at ICU admission (icu_day0)")
control_arm <- read_table(control_dir, "analysis_cross_sectional", cross_sectional_columns) %>% mutate(arm = "Control")
message("Control arm: ", nrow(control_arm), " no-support patients; ",
        sum(control_arm$sf_ratio < SF_HYPOXEMIA_THRESHOLD), " with index SF < ", SF_HYPOXEMIA_THRESHOLD)

# T0 = the first ICU admission of the stay; the first IMV record cuts the window
first_icu_admission <- function(cohort_dir) read_table(cohort_dir, "cohort_icu_stays", c("hospitalization_id", "in_dttm")) %>%
  summarise(t0 = min(in_dttm), .by = hospitalization_id)
first_imv_record <- function(waterfall) waterfall %>%
  filter(tolower(device_category) == "imv" | (!is.na(tidal_volume_set) & tidal_volume_set > 0)) %>%
  summarise(first_imv_dttm = min(recorded_dttm), .by = hospitalization_id)
waterfall_columns <- c("hospitalization_id", "recorded_dttm", "device_category", "fio2_set", "lpm_set", "tidal_volume_set")
# FiO2 on room air and nasal cannula, the one rule applied to both arms
estimate_fio2 <- function(waterfall) {
  device <- tolower(waterfall$device_category)
  waterfall %>% mutate(fio2_set = case_when(!is.na(fio2_set) ~ as.numeric(fio2_set),
                                            device == "room air" ~ 0.21,
                                            device == "nasal cannula" & !is.na(lpm_set) ~ pmin(0.21 + 0.03 * as.numeric(lpm_set), 0.60),
                                            TRUE ~ NA_real_))
}

# =============================================================================
# Pre-decision covariates, per cohort
# =============================================================================
# The window per patient: [T0 - PRE_WINDOW_H, window_end), window_end = T0, or the
# first IMV record if that comes first.
pre_decision_covariates <- function(cohort_dir, arm_patients) {
  waterfall <- read_table(cohort_dir, "resp_support_waterfall_clean", waterfall_columns) %>%
    filter(hospitalization_id %in% arm_patients$hospitalization_id)
  windows <- arm_patients %>% select(hospitalization_id, admission_dttm) %>%
    left_join(first_icu_admission(cohort_dir), by = "hospitalization_id") %>%
    left_join(first_imv_record(waterfall), by = "hospitalization_id")
  if (anyNA(windows$t0)) stop(cohort_dir, ": ", sum(is.na(windows$t0)), " patients without an ICU admission in cohort_icu_stays")
  windows <- windows %>%
    mutate(window_start = t0 - PRE_WINDOW_H * 3600,
           intubated_before_icu = !is.na(first_imv_dttm) & first_imv_dttm < t0,
           window_end = if_else(intubated_before_icu, first_imv_dttm, t0))
  waterfall <- estimate_fio2(waterfall)
  # measurements in [window_start, t0): kept if before window_end, removed (and counted) otherwise
  in_window <- function(measurements, time_column) {
    measurements %>%
      inner_join(windows %>% select(hospitalization_id, window_start, window_end, t0), by = "hospitalization_id") %>%
      filter(.data[[time_column]] >= window_start, .data[[time_column]] < t0) %>%
      mutate(removed_by_imv_cut = .data[[time_column]] >= window_end)
  }
  removed_patients <- character(0)
  keep_before_imv <- function(windowed) {
    removed_patients <<- union(removed_patients, as.character(windowed$hospitalization_id[windowed$removed_by_imv_cut]))
    windowed %>% filter(!removed_by_imv_cut)
  }

  vitals <- read_table(cohort_dir, "cohort_vitals_clean", c("hospitalization_id", "recorded_dttm", "vital_category", "vital_value")) %>%
    filter(hospitalization_id %in% arm_patients$hospitalization_id, vital_category %in% c("spo2", "map"), !is.na(vital_value)) %>%
    mutate(vital_value = as.numeric(vital_value))
  labs <- read_table(cohort_dir, "cohort_labs_clean", c("hospitalization_id", "lab_result_dttm", "lab_category", "lab_value_numeric")) %>%
    filter(hospitalization_id %in% arm_patients$hospitalization_id, !is.na(lab_value_numeric),
           lab_category %in% c("ph_arterial", "pco2_arterial", "po2_arterial", "creatinine", "bilirubin_total", "platelet_count"))
  gcs <- read_table(cohort_dir, "cohort_assessments", c("hospitalization_id", "recorded_dttm", "assessment_category", "numerical_value")) %>%
    filter(hospitalization_id %in% arm_patients$hospitalization_id, assessment_category == "gcs_total") %>%
    mutate(gcs = as.numeric(numerical_value)) %>% filter(!is.na(gcs))

  # SpO2 and PaO2 paired with the most recent FiO2 up to FIO2_LOOKBACK_H before them (script 03)
  fio2_rows <- as.data.table(waterfall %>% filter(!is.na(fio2_set)) %>%
                               transmute(hospitalization_id, t = as.numeric(recorded_dttm), fio2_set))
  setkey(fio2_rows, hospitalization_id, t)
  with_fio2 <- function(measurements, time_column) {
    query <- as.data.table(measurements %>% mutate(t = as.numeric(.data[[time_column]])))
    as_tibble(fio2_rows[query, roll = FIO2_LOOKBACK_H * 3600, on = .(hospitalization_id, t)])
  }

  waterfall_window <- keep_before_imv(in_window(waterfall, "recorded_dttm"))
  vitals_window    <- keep_before_imv(in_window(vitals, "recorded_dttm"))
  labs_window      <- keep_before_imv(in_window(labs, "lab_result_dttm"))
  gcs_window       <- keep_before_imv(in_window(gcs, "recorded_dttm"))

  oxygenation <- waterfall_window %>%
    mutate(support_level = unname(SUPPORT_LEVEL_OF_DEVICE[tolower(device_category)]),
           support_rank = match(support_level, SUPPORT_LEVELS)) %>%
    summarise(highest_fio2 = if (any(!is.na(fio2_set))) max(fio2_set, na.rm = TRUE) else NA_real_,
              support_rank = if (any(!is.na(support_rank))) max(support_rank, na.rm = TRUE) else NA_integer_,
              .by = hospitalization_id)
  spo2 <- vitals_window %>% filter(vital_category == "spo2")
  sf <- with_fio2(spo2 %>% filter(vital_value >= SPO2_SF_RANGE[1], vital_value <= SPO2_SF_RANGE[2]) %>%
                    select(hospitalization_id, recorded_dttm, vital_value), "recorded_dttm") %>%
    filter(!is.na(fio2_set)) %>%
    summarise(lowest_sf = min(vital_value / fio2_set), .by = hospitalization_id)
  pf <- with_fio2(labs_window %>% filter(lab_category == "po2_arterial") %>%
                    select(hospitalization_id, lab_result_dttm, lab_value_numeric), "lab_result_dttm") %>%
    filter(!is.na(fio2_set)) %>%
    summarise(lowest_pf = min(lab_value_numeric / fio2_set), .by = hospitalization_id)
  vitals_summary <- vitals_window %>%
    summarise(lowest_spo2 = if (any(vital_category == "spo2")) min(vital_value[vital_category == "spo2"]) else NA_real_,
              lowest_map  = if (any(vital_category == "map"))  min(vital_value[vital_category == "map"])  else NA_real_,
              .by = hospitalization_id)
  lab_extreme <- function(category, direction, name) {
    labs_window %>% filter(lab_category == category) %>%
      summarise(value = if (direction == "lowest") min(lab_value_numeric) else max(lab_value_numeric), .by = hospitalization_id) %>%
      rename(!!name := value)
  }
  lab_summary <- list(lab_extreme("ph_arterial", "lowest", "lowest_ph"),
                      lab_extreme("pco2_arterial", "highest", "highest_paco2"),
                      lab_extreme("creatinine", "highest", "highest_creatinine"),
                      lab_extreme("bilirubin_total", "highest", "highest_bilirubin"),
                      lab_extreme("platelet_count", "lowest", "lowest_platelets")) %>%
    reduce(full_join, by = "hospitalization_id")
  gcs_summary <- gcs_window %>% summarise(lowest_gcs = min(gcs), .by = hospitalization_id)

  # NE-equivalent dose: rows in the window, plus the dose in force when it opens (the
  # series' last row at or before the window start; it holds an explicit zero row where
  # a dose ends)
  ne_equiv <- read_table(cohort_dir, "ne_equiv_admin", c("hospitalization_id", "admin_dttm", "ne_equiv_total")) %>%
    filter(hospitalization_id %in% arm_patients$hospitalization_id)
  ne_in_window <- keep_before_imv(in_window(ne_equiv, "admin_dttm")) %>% select(hospitalization_id, ne_equiv_total)
  ne_at_open <- ne_equiv %>%
    inner_join(windows %>% select(hospitalization_id, window_start, window_end), by = "hospitalization_id") %>%
    filter(admin_dttm <= window_start, window_start < window_end) %>%
    slice_max(admin_dttm, n = 1, with_ties = FALSE, by = hospitalization_id) %>%
    select(hospitalization_id, ne_equiv_total)
  ne_summary <- bind_rows(ne_in_window, ne_at_open) %>%
    summarise(highest_ne_equiv = max(ne_equiv_total), .by = hospitalization_id)

  covariates <- arm_patients %>%
    left_join(windows %>% select(hospitalization_id, t0, intubated_before_icu), by = "hospitalization_id") %>%
    left_join(oxygenation, by = "hospitalization_id") %>%
    left_join(sf, by = "hospitalization_id") %>%
    left_join(pf, by = "hospitalization_id") %>%
    left_join(vitals_summary, by = "hospitalization_id") %>%
    left_join(lab_summary, by = "hospitalization_id") %>%
    left_join(gcs_summary, by = "hospitalization_id") %>%
    left_join(ne_summary, by = "hospitalization_id") %>%
    mutate(highest_ne_equiv = coalesce(highest_ne_equiv, 0),
           any_vasopressor = as.integer(highest_ne_equiv > 0),
           support_level = factor(SUPPORT_LEVELS[support_rank], levels = SUPPORT_LEVELS, ordered = TRUE),
           hours_admission_to_t0 = as.numeric(difftime(t0, admission_dttm, units = "hours")),
           measurement_removed_by_imv_cut = as.character(hospitalization_id) %in% removed_patients,
           index_hours_after_first_icu = as.numeric(difftime(index_dttm, t0, units = "hours"))) %>%
    select(-support_rank)
  covariates
}
ventilated_arm <- pre_decision_covariates(ventilated_dir, ventilated_arm)
control_arm    <- pre_decision_covariates(control_dir, control_arm)
both_arms <- bind_rows(ventilated_arm, control_arm) %>%
  mutate(arm = factor(arm, levels = ARMS),
         sex_category  = factor(sex_category, levels = c("Male", "Female")),
         race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
         log_pfvc = log(pfvc),
         anchor = sofa_cv_97 + sofa_coag + sofa_liver + sofa_renal)
structural_missing <- both_arms %>% select(age_at_admission, sex_category, race_category, height_cm, pfvc, anchor) %>%
  summarise(across(everything(), ~ sum(is.na(.x))))
if (any(structural_missing > 0)) {
  print(as.data.frame(structural_missing))
  stop("patients missing age, sex, race, height, PFVC or the severity anchor in analysis_cross_sectional: rebuild with script 03")
}

time_zero <- both_arms %>% group_by(arm) %>%
  summarise(n_patients = n(),
            n_intubated_before_icu_admission = sum(intubated_before_icu),
            n_with_measurement_removed_by_imv_cut = sum(measurement_removed_by_imv_cut),
            n_index_more_than_6h_after_first_icu_admission = sum(index_hours_after_first_icu > CONTROL_ICU_WINDOW_H),
            median_hours_admission_to_t0 = median(hours_admission_to_t0), .groups = "drop") %>%
  mutate(pre_window_h = PRE_WINDOW_H, site = site_name)
message("\nTime zero (first ICU admission) and the IMV cut, per arm:")
print(as.data.frame(time_zero), row.names = FALSE)

# =============================================================================
# The two comparisons: the score and its diagnostics
# =============================================================================
COMPARISONS <- c(all_controls = "ventilated vs no support (all)",
                 hypoxemic_controls = "ventilated vs no support, index SF < 315")
comparison_data <- function(comparison) {
  if (comparison == "all_controls") both_arms else
    both_arms %>% filter(arm == "Ventilated" | sf_ratio < SF_HYPOXEMIA_THRESHOLD)
}
score_covariates <- PS_COVARIATES %>% filter(PS_INCLUDE_STRUCTURAL | !structural)

# share missing per covariate and arm
missing_share <- map_dfr(names(COMPARISONS), function(comparison) {
  comparison_data(comparison) %>%
    group_by(arm) %>%
    summarise(across(all_of(PS_COVARIATES$covariate), ~ mean(is.na(.x))), n_patients = n(), .groups = "drop") %>%
    pivot_longer(all_of(PS_COVARIATES$covariate), names_to = "covariate", values_to = "share_missing") %>%
    mutate(comparison = COMPARISONS[[comparison]])
}) %>%
  left_join(PS_COVARIATES %>% select(covariate, label), by = "covariate") %>%
  mutate(site = site_name) %>%
  select(comparison, arm, covariate, label, n_patients, share_missing, site)

# Missing indicators and pooled-median values; the design and the model formula
prepare_design <- function(data) {
  not_in_score <- character(0); linear <- character(0); terms <- character(0)
  for (row in seq_len(nrow(PS_COVARIATES))) {
    covariate <- PS_COVARIATES$covariate[row]; type <- PS_COVARIATES$type[row]
    in_score <- covariate %in% score_covariates$covariate
    values <- data[[covariate]]
    n_missing <- sum(is.na(values))
    if (n_missing == nrow(data)) {
      if (in_score) not_in_score <- c(not_in_score, covariate)
      next
    }
    if (n_missing > 0) {
      indicator <- paste0(covariate, "_missing")
      data[[indicator]] <- as.integer(is.na(values))
      data[[covariate]][is.na(values)] <- if (type == "categorical")
        levels(values)[quantile(as.integer(values), 0.5, type = 1, na.rm = TRUE)] else median(values, na.rm = TRUE)
      if (in_score) terms <- c(terms, indicator)
    }
    if (!in_score) next
    # a covariate with one observed value carries nothing beyond its missing indicator
    if (n_distinct(data[[covariate]]) < 2) {
      not_in_score <- c(not_in_score, covariate)
      next
    }
    if (type == "continuous") {
      if (n_distinct(quantile(data[[covariate]], c(0, 1 / 3, 2 / 3, 1), names = FALSE)) < 4) {
        linear <- c(linear, covariate); terms <- c(terms, covariate)
      } else terms <- c(terms, paste0("ns(", covariate, ", 3)"))
    } else if (type == "categorical") {
      terms <- c(terms, paste0("factor(", covariate, ", ordered = FALSE)"))
    } else terms <- c(terms, covariate)
  }
  list(data = data, formula = as.formula(paste("treated ~", paste(terms, collapse = " + "))),
       not_in_score = not_in_score, linear = linear)
}

# the fit: separation is recorded and the fit kept; any other warning stops
fit_score <- function(formula, data) {
  separation_note <- NA_character_
  fit <- withCallingHandlers(glm(formula, family = binomial, data = data), warning = function(w) {
    if (grepl(SEPARATION_PATTERN, conditionMessage(w))) {
      separation_note <<- paste("separation:", trimws(conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
    stop("propensity model warning, stopping: ", conditionMessage(w), call. = FALSE)
  })
  list(fit = fit, separation_note = separation_note)
}

auc_of <- function(score, treated) {
  ranks <- rank(score)
  n_treated <- sum(treated); n_control <- sum(!treated)
  (sum(ranks[treated]) - n_treated * (n_treated + 1) / 2) / (n_treated * n_control)
}
kish_ess <- function(w) sum(w)^2 / sum(w^2)
top_share <- function(w, fraction = 0.01) {
  n_top <- max(1L, ceiling(fraction * length(w)))
  sum(sort(w, decreasing = TRUE)[seq_len(n_top)]) / sum(w)
}
weighted_var <- function(x, w) sum(w * (x - weighted.mean(x, w))^2) / sum(w)

# the balance terms: each continuous or binary covariate, each level of a categorical,
# each missing indicator
balance_terms <- function(data) {
  columns <- list()
  for (row in seq_len(nrow(PS_COVARIATES))) {
    covariate <- PS_COVARIATES$covariate[row]; label <- PS_COVARIATES$label[row]
    if (!covariate %in% names(data) || all(is.na(data[[covariate]]))) next
    if (PS_COVARIATES$type[row] == "categorical") {
      for (level in levels(data[[covariate]]))
        columns[[paste0(label, ": ", level)]] <- as.numeric(data[[covariate]] == level)
    } else columns[[label]] <- as.numeric(data[[covariate]])
    indicator <- paste0(covariate, "_missing")
    if (indicator %in% names(data)) columns[[paste0(label, " (missing)")]] <- data[[indicator]]
  }
  as_tibble(columns)
}
smd_table <- function(terms, treated, weights) {
  map_dfr(names(terms), function(term) {
    x <- terms[[term]]
    pooled_sd <- sqrt((var(x[treated]) + var(x[!treated])) / 2)
    map_dfr(names(weights), function(weighting) {
      w <- weights[[weighting]]
      difference <- weighted.mean(x[treated], w[treated]) - weighted.mean(x[!treated], w[!treated])
      tibble(term = term, weighting = weighting, mean_ventilated = weighted.mean(x[treated], w[treated]),
             mean_control = weighted.mean(x[!treated], w[!treated]),
             smd = if (pooled_sd > 0) difference / pooled_sd else NA_real_)
    })
  })
}

run_comparison <- function(comparison) {
  design <- prepare_design(comparison_data(comparison) %>% mutate(treated = as.integer(arm == "Ventilated")))
  data <- design$data
  if (length(design$not_in_score))
    message(COMPARISONS[[comparison]], ": left out of the score (never observed or constant): ",
            paste(design$not_in_score, collapse = ", "))
  if (length(design$linear))
    message(COMPARISONS[[comparison]], ": entered linearly (spline knots not distinct): ",
            paste(design$linear, collapse = ", "))
  score_fit <- fit_score(design$formula, data)
  if (!is.na(score_fit$separation_note)) message(COMPARISONS[[comparison]], ": ", score_fit$separation_note)
  ps <- fitted(score_fit$fit)
  treated <- data$treated == 1
  weights <- list(unweighted = rep(1, nrow(data)),
                  att = if_else(treated, 1, ps / (1 - ps)),
                  overlap = if_else(treated, 1 - ps, ps))
  common_lo <- max(min(ps[treated]), min(ps[!treated]))
  common_hi <- min(max(ps[treated]), max(ps[!treated]))
  summary_rows <- expand_grid(weighting = names(WEIGHTINGS), arm = ARMS) %>%
    pmap_dfr(function(weighting, arm) {
      in_arm <- if (arm == "Ventilated") treated else !treated
      w <- weights[[weighting]][in_arm]
      tibble(weighting = WEIGHTINGS[[weighting]], arm = arm, n_patients = sum(in_arm),
             ps_min = min(ps[in_arm]), ps_median = median(ps[in_arm]), ps_max = max(ps[in_arm]),
             share_in_common_support = mean(ps[in_arm] >= common_lo & ps[in_arm] <= common_hi),
             ess = kish_ess(w), weight_p50 = quantile(w, 0.5, names = FALSE), weight_p90 = quantile(w, 0.9, names = FALSE),
             weight_p99 = quantile(w, 0.99, names = FALSE), weight_max = max(w),
             share_of_arm_weight_top_1pct = top_share(w),
             mean_log_pfvc = weighted.mean(data$log_pfvc[in_arm], w),
             sd_log_pfvc = sqrt(weighted_var(data$log_pfvc[in_arm], w)),
             mean_anchor = weighted.mean(data$anchor[in_arm], w))
    }) %>%
    mutate(comparison = COMPARISONS[[comparison]], auc = auc_of(ps, treated),
           n_ventilated = sum(treated), n_control = sum(!treated),
           common_support_lo = common_lo, common_support_hi = common_hi,
           n_score_terms = length(coef(score_fit$fit)) - 1L,
           n_score_terms_aliased = sum(is.na(coef(score_fit$fit))),
           separation_note = score_fit$separation_note,
           covariates_not_in_score = paste(design$not_in_score, collapse = ";"),
           covariates_linear = paste(design$linear, collapse = ";"),
           ps_include_structural = PS_INCLUDE_STRUCTURAL, pre_window_h = PRE_WINDOW_H, site = site_name, .before = 1)
  balance <- smd_table(balance_terms(data), treated, weights) %>%
    mutate(abs_smd = abs(smd), imbalanced = abs_smd > SMD_FLAG, weighting = WEIGHTINGS[weighting],
           comparison = COMPARISONS[[comparison]], site = site_name, .before = 1)
  bins <- tibble(bin_lo = seq(0, 1 - PS_BIN_WIDTH, by = PS_BIN_WIDTH)) %>% mutate(bin_hi = bin_lo + PS_BIN_WIDTH)
  bin_index <- pmin(floor(ps / PS_BIN_WIDTH) + 1, nrow(bins))
  ps_bins <- expand_grid(arm = ARMS, bin = seq_len(nrow(bins))) %>%
    mutate(n_patients = map2_int(arm, bin, ~ sum(bin_index == .y & (if (.x == "Ventilated") treated else !treated))),
           bin_lo = bins$bin_lo[bin], bin_hi = bins$bin_hi[bin]) %>%
    select(-bin) %>%
    mutate(comparison = COMPARISONS[[comparison]], site = site_name, .before = 1)
  list(summary = summary_rows, balance = balance, ps_bins = ps_bins)
}
results <- map(set_names(names(COMPARISONS)), run_comparison)
overlap_summary <- map_dfr(results, "summary")
overlap_balance <- map_dfr(results, "balance")
overlap_ps      <- map_dfr(results, "ps_bins")

message("\nSummary (ESS and weights per arm and weighting):")
print(as.data.frame(overlap_summary %>% select(comparison, weighting, arm, n_patients, auc, common_support_lo, common_support_hi,
                                               share_in_common_support, ess, weight_p99, weight_max,
                                               share_of_arm_weight_top_1pct, mean_log_pfvc, mean_anchor)),
      row.names = FALSE, digits = 3)
message("\nCovariate terms with |SMD| > ", SMD_FLAG, ":")
print(as.data.frame(overlap_balance %>% count(comparison, weighting, n_imbalanced = imbalanced) %>%
                      filter(n_imbalanced) %>% select(-n_imbalanced)), row.names = FALSE)

write_csv(overlap_summary, file.path(final_dir, paste0("intubation_overlap_summary_", site_name, ".csv")))
write_csv(overlap_balance, file.path(final_dir, paste0("intubation_overlap_balance_", site_name, ".csv")))
write_csv(overlap_ps,      file.path(final_dir, paste0("intubation_overlap_ps_", site_name, ".csv")))
write_csv(missing_share,   file.path(final_dir, paste0("intubation_overlap_missing_", site_name, ".csv")))
write_csv(time_zero,       file.path(final_dir, paste0("intubation_overlap_time_zero_", site_name, ".csv")))

# =============================================================================
# Figure: per comparison, the mirrored score and the love plot
# =============================================================================
theme_set(theme_minimal(base_size = 10))
comparison_page <- function(comparison_label) {
  mirrored <- overlap_ps %>% filter(comparison == comparison_label) %>%
    mutate(height = if_else(arm == "Ventilated", n_patients, -n_patients), arm = factor(arm, levels = ARMS))
  support <- overlap_summary %>% filter(comparison == comparison_label) %>% slice(1)
  histogram <- ggplot(mirrored, aes(x = (bin_lo + bin_hi) / 2, y = height, fill = arm)) +
    geom_col(width = PS_BIN_WIDTH) +
    geom_hline(yintercept = 0, colour = "grey30") +
    annotate("rect", xmin = support$common_support_lo, xmax = support$common_support_hi, ymin = -Inf, ymax = Inf,
             alpha = 0.08, fill = "grey20") +
    scale_fill_manual(values = ARM_COLOURS, name = NULL) +
    scale_y_continuous(labels = abs) +
    labs(x = "Propensity score for intubation", y = "Patients (ventilated up, control down)",
         subtitle = sprintf("AUC %.3f; common support %.3f to %.3f (shaded)", support$auc,
                            support$common_support_lo, support$common_support_hi)) +
    theme(legend.position = "top")
  # terms ordered by their unweighted |SMD|; a term constant in both arms has no SMD and is not drawn
  love <- overlap_balance %>% filter(comparison == comparison_label, !is.na(abs_smd))
  term_order <- love %>% filter(weighting == WEIGHTINGS[["unweighted"]]) %>% arrange(abs_smd) %>% pull(term)
  love <- love %>% mutate(weighting = factor(weighting, levels = WEIGHTINGS), term = factor(term, levels = term_order))
  love_plot <- ggplot(love, aes(x = abs_smd, y = term, colour = weighting)) +
    geom_vline(xintercept = SMD_FLAG, linetype = "dashed", colour = "grey40") +
    geom_point(size = 1.8) +
    scale_colour_manual(values = setNames(WEIGHTING_COLOURS, WEIGHTINGS), name = NULL) +
    labs(x = "|Standardised mean difference|", y = NULL) +
    theme(legend.position = "top")
  (histogram | love_plot) + plot_layout(widths = c(1, 1.2)) +
    plot_annotation(title = comparison_label,
                    caption = paste0("Covariates from ", PRE_WINDOW_H, " h before ICU admission; ",
                                     if (PS_INCLUDE_STRUCTURAL) "age, sex, race and height in the score" else
                                       "age, sex, race and height not in the score", ". ", site_name))
}
pdf_path <- file.path(final_dir, paste0("intubation_overlap_", site_name, ".pdf"))
pdf(pdf_path, width = 13, height = 8)
for (comparison_label in COMPARISONS) print(comparison_page(comparison_label))
invisible(dev.off())
message("\nWrote intubation_overlap_{summary,balance,ps,missing,time_zero}_", site_name, ".csv and the PDF to ", final_dir)
