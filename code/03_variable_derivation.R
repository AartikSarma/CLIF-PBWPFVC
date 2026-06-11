# =============================================================================
# Script 03: Variable Derivation
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(data.table)
library(lubridate)
library(rspiro)

source("utils/config.R")
source("utils/process_sofa_scores.R")
source("utils/standardize_pressor_dose.R")
source("utils/attrition_log.R")

site_name <- config$site_name
output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

# Fixed bin edges for the federated PBW:PFVC distribution exports (must match
# across sites so per-site histograms can be summed into a pooled distribution).
PBWPFVC_BIN_EDGES <- seq(0, 50, by = 1)
AGE_BIN_EDGES     <- c(18, 30, 40, 50, 60, 70, 80, Inf)
HEIGHT_BIN_EDGES  <- c(150, 155, 160, 165, 170, 175, 180, 185, 190, 210)
DIST_MIN_CELL     <- 10   # group-level small-cell suppression (CLAUDE.md)

# =============================================================================
# Load cleaned intermediate data
# =============================================================================

cohort_ids <- readRDS(file.path(output_dir, "cohort_hospitalization_ids.rds"))
resp_waterfall <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet"))
cohort_demographics <- read_parquet(file.path(output_dir, "cohort_demographics.parquet"))
cohort_vitals <- read_parquet(file.path(output_dir, "cohort_vitals_clean.parquet"))
cohort_labs <- read_parquet(file.path(output_dir, "cohort_labs_clean.parquet"))
cohort_heights <- read_parquet(file.path(output_dir, "cohort_heights_clean.parquet"))
cohort_weights <- read_parquet(file.path(output_dir, "cohort_weights.parquet"))
cohort_meds <- read_parquet(file.path(output_dir, "cohort_meds.parquet"))
cohort_assessments <- read_parquet(file.path(output_dir, "cohort_assessments.parquet"))

message("Loaded all cleaned intermediate data")

# =============================================================================
# 3a. PBW and PFVC calculation
# =============================================================================

# Join demographics with heights
pbw_pfvc_data <- cohort_demographics %>%
  select(hospitalization_id, patient_id, age_at_admission,
         sex_category, race_category, deceased) %>%
  inner_join(cohort_heights, by = "hospitalization_id") %>%
  # Height inclusion range from the original analysis (PBWvsFVC): 150-210 cm.
  # The GLI-2012 reference equations are only valid within this range, and the
  # published cohort excluded heights outside it. Restricting here cascades to all
  # downstream analyses, since they inner-join pbw_pfvc_data.
  filter(!is.na(height_cm), height_cm >= 150, height_cm <= 210) %>%
  mutate(
    # Map sex to numeric for pred_GLI: 1=male, 2=female
    sex_numeric = case_when(
      sex_category == "Male"   ~ 1L,
      sex_category == "Female" ~ 2L,
      TRUE                     ~ NA_integer_
    ),
    # GLI ethnicity: 1=Caucasian, 2=African-American, 5=Other/Mixed
    race_numeric = case_when(
      race_category == "WHITE" ~ 1L,
      race_category == "BLACK" ~ 2L,
      TRUE                     ~ 5L
    ),
    # PBW via Devine formula (kg)
    pbw = case_when(
      sex_numeric == 1 ~ 50.0 + 2.3 * (height_cm / 2.54 - 60),
      sex_numeric == 2 ~ 45.5 + 2.3 * (height_cm / 2.54 - 60)
    ),
    # PFVC via GLI-2012 (litres)
    pfvc = pred_GLI(
      age       = age_at_admission,
      height    = height_cm / 100,
      gender    = sex_numeric,
      ethnicity = race_numeric,
      param     = "FVC"
    )
  ) %>%
  filter(!is.na(pbw), !is.na(pfvc), pfvc > 0) %>%
  # Reference categories: Male and White (first factor level). Releveling here
  # cascades to every downstream model (analysis_data inner-joins these columns),
  # so all regressions report effects relative to male / white patients.
  mutate(
    sex_category  = factor(sex_category,  levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER"))
  )

message("PBW/PFVC calculated: ", nrow(pbw_pfvc_data), " hospitalizations")

# Broad PFVC cohort for the H1 (PFVC vs PBW) table: all eligible patients with
# height/age/sex/race/PFVC, before any ventilation-based inclusion. Saved for
# script 04.
write_parquet(pbw_pfvc_data, file.path(output_dir, "analysis_broad_pfvc.parquet"))

# =============================================================================
# 3b. Compute SF and PF ratios at measurement time using concurrent FiO2
# =============================================================================
# Ratios are computed at the time SpO2/PaO2 were measured, using the FiO2
# that was active at that moment (most recent prior FiO2 within 4h).
# These pre-computed ratios are then carried forward to IMV timepoints.

# --- FiO2 from waterfall (for matching to SpO2/PaO2 measurement times) ---
fio2_data <- resp_waterfall %>%
  filter(!is.na(fio2_set)) %>%
  select(hospitalization_id, recorded_dttm, fio2_set) %>%
  mutate(fio2_dttm = as.numeric(recorded_dttm))

fio2_dt <- as.data.table(fio2_data)
setkey(fio2_dt, hospitalization_id, fio2_dttm)

# --- SF ratio: SpO2 / FiO2 at time of SpO2 measurement ---
spo2_data <- cohort_vitals %>%
  filter(vital_category == "spo2", !is.na(vital_value)) %>%
  # SpO2 inclusion range from the original analysis (PBWvsFVC): 80-97. The SF
  # ratio is only valid within the (near-)linear part of the oxyhemoglobin
  # dissociation curve; SpO2 >= 98 saturates and SpO2 < 80 falls off the curve.
  # Restricting the SpO2 that defines the SF ratio enforces the original's
  # `o2sat <= 97 & o2sat >= 80` hypoxemia criterion (the <= 97 bound was already
  # applied in script 01; this adds the missing >= 80 floor).
  filter(vital_value >= 80, vital_value <= 97) %>%
  select(hospitalization_id, recorded_dttm, vital_value) %>%
  rename(spo2_value = vital_value, spo2_dttm = recorded_dttm) %>%
  mutate(spo2_dttm_num = as.numeric(spo2_dttm))

spo2_dt <- as.data.table(spo2_data)
setkey(spo2_dt, hospitalization_id, spo2_dttm_num)

# Rolling join: for each SpO2, find most recent prior FiO2 within 4h
sf_joined <- fio2_dt[spo2_dt,
                      roll = 4 * 3600,
                      on = .(hospitalization_id, fio2_dttm = spo2_dttm_num)]

sf_joined <- as_tibble(sf_joined) %>%
  filter(!is.na(fio2_set), !is.na(spo2_value)) %>%
  mutate(sf_ratio = spo2_value / fio2_set,
         sf_dttm = fio2_dttm)  # timestamp = when SpO2 was measured

message("SF ratio computed at ", nrow(sf_joined), " SpO2 measurement times")

# --- PF ratio: PaO2 / FiO2 at time of PaO2 measurement ---
pao2_data <- cohort_labs %>%
  filter(lab_category == "po2_arterial", !is.na(lab_value_numeric)) %>%
  select(hospitalization_id, lab_result_dttm, lab_value_numeric) %>%
  rename(pao2_value = lab_value_numeric, pao2_dttm = lab_result_dttm) %>%
  mutate(pao2_dttm_num = as.numeric(pao2_dttm))

pao2_dt <- as.data.table(pao2_data)
setkey(pao2_dt, hospitalization_id, pao2_dttm_num)

pf_joined <- fio2_dt[pao2_dt,
                      roll = 4 * 3600,
                      on = .(hospitalization_id, fio2_dttm = pao2_dttm_num)]

pf_joined <- as_tibble(pf_joined) %>%
  filter(!is.na(fio2_set), !is.na(pao2_value)) %>%
  mutate(pf_ratio = pao2_value / fio2_set,
         pf_dttm = fio2_dttm)  # timestamp = when PaO2 was measured

message("PF ratio computed at ", nrow(pf_joined), " PaO2 measurement times")

# =============================================================================
# 3c. SOFA scores — daily (worst values per calendar day)
# =============================================================================

# Aggregate each data source to worst-per-day BEFORE joining,
# so the join produces compact rows instead of millions of sparse ones.

admission_times <- cohort_demographics %>%
  select(hospitalization_id, admission_dttm)

# Helper: add sofa_day and filter to the ventilation window.
#
# Daily SOFA is required by the cross-sectional inclusion gate (has_all_data needs
# a same-day sofa_total). Previously the window was the first 72h (days 0-2), which
# silently excluded any patient whose first qualifying IMV timepoint fell later than
# 72h after admission. The original analysis required only a non-missing stay-level
# SOFA, so late-intubation patients were retained. We extend the window to span the
# 28-day study/ventilation window so a same-day SOFA exists on whichever day the
# index timepoint lands. (Pre-aggregation to patient-days keeps this inexpensive.)
SOFA_MAX_DAY <- 28L
add_day <- function(df, dttm_col = "recorded_dttm") {
  df %>%
    inner_join(admission_times, by = "hospitalization_id") %>%
    mutate(sofa_day = as.integer(floor(as.numeric(
      difftime(.data[[dttm_col]], admission_dttm, units = "days")
    )))) %>%
    filter(sofa_day >= 0L, sofa_day <= SOFA_MAX_DAY) %>%
    select(-admission_dttm)
}

# Labs: worst per day (max for creatinine/bilirubin, min for platelets/PaO2)
# Worst lab per day in long format, then pivot once
labs_daily <- cohort_labs %>%
  filter(!is.na(lab_value_numeric),
         lab_category %in% c("creatinine", "bilirubin_total", "platelet_count", "po2_arterial")) %>%
  mutate(lab_category = case_when(
    lab_category == "platelets" ~ "platelet_count",
    TRUE ~ lab_category
  )) %>%
  select(hospitalization_id, lab_result_dttm, lab_category, lab_value_numeric) %>%
  add_day("lab_result_dttm") %>%
  group_by(hospitalization_id, sofa_day, lab_category) %>%
  summarise(
    # max for worse-when-higher, min for worse-when-lower
    lab_value_numeric = case_when(
      first(lab_category) %in% c("creatinine", "bilirubin_total") ~ max(lab_value_numeric, na.rm = TRUE),
      TRUE ~ min(lab_value_numeric, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  pivot_wider(
    id_cols = c(hospitalization_id, sofa_day),
    names_from = lab_category,
    values_from = lab_value_numeric
  )

# Vitals: worst per day (min MAP, min SpO2)
vitals_daily <- cohort_vitals %>%
  filter(!is.na(vital_value), vital_category %in% c("map", "spo2")) %>%
  select(hospitalization_id, recorded_dttm, vital_category, vital_value) %>%
  add_day() %>%
  group_by(hospitalization_id, sofa_day, vital_category) %>%
  summarise(vital_value = min(vital_value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    id_cols = c(hospitalization_id, sofa_day),
    names_from = vital_category,
    values_from = vital_value
  )

# GCS: worst (min) per day
gcs_daily <- cohort_assessments %>%
  filter(assessment_category == "gcs_total") %>%
  mutate(gcs_total = as.numeric(numerical_value)) %>%
  select(hospitalization_id, recorded_dttm, gcs_total) %>%
  add_day() %>%
  group_by(hospitalization_id, sofa_day) %>%
  summarise(gcs_total = min(gcs_total, na.rm = TRUE), .groups = "drop") %>%
  mutate(gcs_total = if_else(is.infinite(gcs_total), NA_real_, gcs_total))

# Vasopressors: max dose per day
# Norepinephrine doses are standardized to mcg/kg/min using patient weight so
# that mcg/min and mcg/kg/min entries are comparable.
vaso_daily <- cohort_meds %>%
  filter(med_category == "norepinephrine") %>%
  standardize_pressor_dose(
    weights = cohort_weights,
    out_col = "norepinephrine_mcg_kg_min",
    label = "norepinephrine (SOFA)"
  ) %>%
  filter(!is.na(norepinephrine_mcg_kg_min)) %>%
  select(hospitalization_id, admin_dttm, norepinephrine_mcg_kg_min) %>%
  add_day("admin_dttm") %>%
  group_by(hospitalization_id, sofa_day) %>%
  summarise(norepinephrine_mcg_kg_min = max(norepinephrine_mcg_kg_min, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(norepinephrine_mcg_kg_min = if_else(is.infinite(norepinephrine_mcg_kg_min),
                                              NA_real_, norepinephrine_mcg_kg_min))

# Respiratory: worst device + max FiO2 per day
resp_daily <- resp_waterfall %>%
  select(hospitalization_id, recorded_dttm, device_category, fio2_set) %>%
  add_day() %>%
  mutate(device_category = str_to_title(device_category)) %>%
  group_by(hospitalization_id, sofa_day) %>%
  summarise(
    fio2_set = max(fio2_set, na.rm = TRUE),
    device_category = if (all(is.na(device_category))) NA_character_ else first(na.omit(device_category)),
    .groups = "drop"
  ) %>%
  mutate(fio2_set = if_else(is.infinite(fio2_set), NA_real_, fio2_set))

# Join pre-aggregated daily data (small: ~N_patients * N_days rows)
sofa_collapsed <- labs_daily %>%
  full_join(vitals_daily, by = c("hospitalization_id", "sofa_day")) %>%
  full_join(gcs_daily, by = c("hospitalization_id", "sofa_day")) %>%
  full_join(vaso_daily, by = c("hospitalization_id", "sofa_day")) %>%
  full_join(resp_daily, by = c("hospitalization_id", "sofa_day"))

# Add missing vasopressor columns
for (col in c("epinephrine_mcg_kg_min", "dopamine_mcg_kg_min", "dobutamine_mcg_kg_min")) {
  if (!col %in% names(sofa_collapsed)) {
    sofa_collapsed[[col]] <- NA_real_
  }
}

message("Daily SOFA input: ", nrow(sofa_collapsed), " patient-days, ",
        n_distinct(sofa_collapsed$hospitalization_id), " patients")

# Data is already pre-collapsed to one row per patient-day.
# Score directly: impute PaO2, join device ranks, compute components.
# No group_by aggregation needed — each row is already one observation unit.
sofa_daily <- sofa_collapsed %>%
  mutate(
    po2_arterial = if_else(between(po2_arterial, 0, 700), po2_arterial, NA_real_),
    fio2_set = if_else(between(fio2_set, 0.21, 1), fio2_set, NA_real_),
    spo2 = if_else(between(spo2, 50, 100), spo2, NA_real_)
  ) %>%
  .impute_pao2_from_spo2() %>%
  left_join(DEVICE_RANK_MAPPING, by = "device_category") %>%
  mutate(
    p_f = po2_arterial / fio2_set,
    p_f_imputed = pao2_imputed / fio2_set,
    sofa_cv_97 = case_when(
      norepinephrine_mcg_kg_min > 0.1 ~ 4L,
      norepinephrine_mcg_kg_min > 0   ~ 3L,
      map < 70  ~ 1L,
      map >= 70 ~ 0L
    ),
    sofa_coag = case_when(
      platelet_count < 20   ~ 4L, platelet_count < 50  ~ 3L,
      platelet_count < 100  ~ 2L, platelet_count < 150 ~ 1L,
      platelet_count >= 150 ~ 0L
    ),
    sofa_liver = case_when(
      bilirubin_total >= 12  ~ 4L, bilirubin_total >= 6   ~ 3L,
      bilirubin_total >= 2   ~ 2L, bilirubin_total >= 1.2 ~ 1L,
      bilirubin_total < 1.2  ~ 0L
    ),
    sofa_resp = case_when(
      p_f < 100 & device_category %in% c("Imv", "Nippv", "Cpap") ~ 4L,
      p_f < 200 & device_category %in% c("Imv", "Nippv", "Cpap") ~ 3L,
      p_f < 300 ~ 2L, p_f < 400 ~ 1L, p_f >= 400 ~ 0L
    ),
    sofa_cns = case_when(
      gcs_total < 6  ~ 4L, gcs_total <= 9  ~ 3L,
      gcs_total <= 12 ~ 2L, gcs_total <= 14 ~ 1L,
      gcs_total == 15 ~ 0L
    ),
    sofa_renal = case_when(
      creatinine >= 5   ~ 4L, creatinine >= 3.5 ~ 3L,
      creatinine >= 2   ~ 2L, creatinine >= 1.2 ~ 1L,
      creatinine < 1.2  ~ 0L
    )
  ) %>%
  mutate(
    across(starts_with("sofa_"), ~ coalesce(.x, 0L)),
    sofa_total = sofa_cv_97 + sofa_coag + sofa_liver + sofa_resp + sofa_cns + sofa_renal
  )

message("Daily SOFA computed: ", nrow(sofa_daily), " patient-days, ",
        n_distinct(sofa_daily$hospitalization_id), " patients")

# Per-encounter SOFA = day 1 SOFA (for cross-sectional analysis in script 04)
sofa_scores <- sofa_daily %>%
  filter(sofa_day == 0) %>%
  select(-sofa_day)

message("Day-1 SOFA available for ", nrow(sofa_scores), " hospitalizations")

# =============================================================================
# 3d. Derive analysis variables
# =============================================================================



# Build per-timepoint dataset from the waterfall output (IMV rows only)
imv_timepoints <- resp_waterfall %>%
  filter(tolower(device_category) == "imv") %>%
  select(hospitalization_id, recorded_dttm, device_category, mode_category,
         fio2_set, tidal_volume_set, peep_set, plateau_pressure_obs,
         any_of("minute_vent_obs"))

# Join with PBW/PFVC data
analysis_data <- imv_timepoints %>%
  inner_join(
    pbw_pfvc_data %>% select(hospitalization_id, age_at_admission,
                              sex_category, race_category, sex_numeric,
                              race_numeric, height_cm, pbw, pfvc, deceased),
    by = "hospitalization_id"
  ) %>%
  mutate(
    pbwpfvc = pbw / pfvc,
    vtpbw = tidal_volume_set / pbw,
    vtpfvc = tidal_volume_set / pfvc * 0.1,
    dp = if_else(
      !is.na(plateau_pressure_obs) & !is.na(peep_set),
      plateau_pressure_obs - peep_set,
      NA_real_
    ),
    crs = if_else(!is.na(dp) & dp > 0, tidal_volume_set / dp, NA_real_),
    # Elastance in cmH2O/L. crs is mL/cmH2O, so 1/crs is cmH2O/mL; the x1000
    # converts to cmH2O/L so downstream tables don't round elastance to 0.0.
    # ers_pbw / ers_pfvc (elastance normalized to PBW / PFVC) inherit this scale.
    ers = if_else(!is.na(crs) & crs > 0, 1000 / crs, NA_real_),
    ers_pbw = ers * pbw,
    ers_pfvc = ers * pfvc
  )

message("Analysis timepoints: ", nrow(analysis_data), " rows, ",
        n_distinct(analysis_data$hospitalization_id), " hospitalizations")

# --- NE equivalents: compute from cohort_meds and rolling-join to IMV timepoints ---
# All catecholamine doses are first standardized to mcg/kg/min; vasopressin is
# handled separately because it is dosed in units/min (units/hr), not mcg-based.
catecholamines <- c("norepinephrine", "epinephrine", "dopamine",
                    "phenylephrine", "dobutamine")

ne_equiv_cat <- cohort_meds %>%
  filter(med_category %in% catecholamines, !is.na(med_dose)) %>%
  standardize_pressor_dose(
    weights = cohort_weights,
    out_col = "dose_mcg_kg_min",
    label = "NE-equiv catecholamines"
  ) %>%
  filter(!is.na(dose_mcg_kg_min)) %>%
  mutate(
    ne_equiv = case_when(
      med_category == "norepinephrine" ~ dose_mcg_kg_min,
      med_category == "epinephrine"    ~ dose_mcg_kg_min,
      med_category == "dopamine"       ~ dose_mcg_kg_min / 100,
      med_category == "phenylephrine"  ~ dose_mcg_kg_min / 10,
      med_category == "dobutamine"     ~ 0
    )
  ) %>%
  select(hospitalization_id, admin_dttm, ne_equiv)

# Vasopressin: units/min -> NE-equiv uses the raw rate (standard practice),
# so we pass through without weight-based conversion. Units/hr entries are
# normalized to units/min.
ne_equiv_vaso <- cohort_meds %>%
  filter(med_category == "vasopressin", !is.na(med_dose)) %>%
  mutate(
    med_dose = as.numeric(med_dose),
    unit_clean = str_to_lower(str_replace_all(med_dose_unit, "\\s+", "")),
    units_per_min = case_when(
      unit_clean %in% c("units/min", "u/min")   ~ med_dose,
      unit_clean %in% c("units/hr", "u/hr",
                        "units/h",  "u/h")      ~ med_dose / 60,
      TRUE                                       ~ NA_real_
    ),
    ne_equiv = units_per_min * 2.5
  ) %>%
  filter(!is.na(ne_equiv)) %>%
  select(hospitalization_id, admin_dttm, ne_equiv)

ne_equiv <- bind_rows(ne_equiv_cat, ne_equiv_vaso) %>%
  filter(ne_equiv > 0) %>%
  group_by(hospitalization_id, admin_dttm) %>%
  summarise(ne_equiv_total = sum(ne_equiv, na.rm = TRUE), .groups = "drop")

ne_dt <- as.data.table(ne_equiv)
ne_dt[, join_dttm := as.numeric(admin_dttm)]
setkey(ne_dt, hospitalization_id, join_dttm)

analysis_ne_dt <- as.data.table(analysis_data)
analysis_ne_dt[, join_dttm := as.numeric(recorded_dttm)]
setkey(analysis_ne_dt, hospitalization_id, join_dttm)

analysis_ne_joined <- ne_dt[, .(hospitalization_id, join_dttm, ne_equiv_total)][
  analysis_ne_dt, roll = 4 * 3600, on = .(hospitalization_id, join_dttm)
]
analysis_ne_joined[, join_dttm := NULL]

analysis_data <- as_tibble(analysis_ne_joined) %>%
  mutate(ne_equiv_total = coalesce(ne_equiv_total, 0))

message("NE equivalents joined: ", sum(analysis_data$ne_equiv_total > 0),
        " timepoints with vasopressors")

# --- Ventilatory ratio: computed at PaCO2 measurement time ---
# VR = (minute_vent_obs * PaCO2) / (PBW * 100 * 37.5)
# Join minute_vent_obs from waterfall to PaCO2 timestamps, compute VR there,
# then forward-join to IMV timepoints.
pco2_data <- cohort_labs %>%
  filter(lab_category == "pco2_arterial", !is.na(lab_value_numeric)) %>%
  select(hospitalization_id, lab_result_dttm, lab_value_numeric) %>%
  rename(paco2 = lab_value_numeric, pco2_dttm = lab_result_dttm) %>%
  mutate(pco2_dttm_num = as.numeric(pco2_dttm))

if (nrow(pco2_data) > 0) {
  # Get minute_vent_obs from waterfall (with timestamps)
  ve_data <- resp_waterfall %>%
    filter(!is.na(minute_vent_obs)) %>%
    select(hospitalization_id, recorded_dttm, minute_vent_obs) %>%
    mutate(ve_dttm_num = as.numeric(recorded_dttm))

  ve_dt <- as.data.table(ve_data)
  setkey(ve_dt, hospitalization_id, ve_dttm_num)

  pco2_dt_vr <- as.data.table(pco2_data)
  setkey(pco2_dt_vr, hospitalization_id, pco2_dttm_num)

  # Join nearest prior minute_vent_obs to each PaCO2 measurement (within 4h)
  vr_at_pco2 <- ve_dt[, .(hospitalization_id, ve_dttm_num, minute_vent_obs)][
    pco2_dt_vr, roll = 4 * 3600, on = .(hospitalization_id, ve_dttm_num = pco2_dttm_num)
  ]

  # Join PBW for each patient
  vr_at_pco2 <- as_tibble(vr_at_pco2) %>%
    left_join(pbw_pfvc_data %>% select(hospitalization_id, pbw), by = "hospitalization_id") %>%
    mutate(
      vent_ratio = if_else(
        !is.na(minute_vent_obs) & !is.na(paco2) & !is.na(pbw) & pbw > 0,
        (minute_vent_obs * paco2) / (pbw * 0.1 * 37.5),
        NA_real_
      ),
      vr_dttm = ve_dttm_num  # timestamp = when PaCO2 was measured
    ) %>%
    rename(vr_minute_vent = minute_vent_obs, vr_paco2 = paco2) %>%
    filter(!is.na(vent_ratio)) %>%
    select(hospitalization_id, vr_dttm, vent_ratio, vr_minute_vent, vr_paco2)

  # Join VR to nearest IMV timepoint within 1h (no long-range forward fill)
  vr_for_join <- as.data.table(vr_at_pco2)
  setkey(vr_for_join, hospitalization_id, vr_dttm)

  analysis_vr_dt <- as.data.table(analysis_data)
  analysis_vr_dt[, recorded_dttm_num := as.numeric(recorded_dttm)]
  setkey(analysis_vr_dt, hospitalization_id, recorded_dttm_num)

  analysis_vr <- vr_for_join[analysis_vr_dt,
                               roll = 1 * 3600,
                               on = .(hospitalization_id, vr_dttm = recorded_dttm_num)]
  analysis_vr[, vr_dttm := NULL]

  analysis_data <- as_tibble(analysis_vr)

  message("Ventilatory ratio available at ", sum(!is.na(analysis_data$vent_ratio)),
          " / ", nrow(analysis_data), " timepoints")
} else {
  analysis_data$vent_ratio <- NA_real_
  message("No PaCO2 data available — ventilatory ratio set to NA")
}

# --- Join pre-computed SF ratio to IMV timepoints (most recent prior within 4h) ---
sf_for_join <- as.data.table(sf_joined %>% select(hospitalization_id, sf_dttm, sf_ratio))
setkey(sf_for_join, hospitalization_id, sf_dttm)

analysis_dt <- as.data.table(analysis_data)
analysis_dt[, recorded_dttm_num := as.numeric(recorded_dttm)]
setkey(analysis_dt, hospitalization_id, recorded_dttm_num)

analysis_with_sf <- sf_for_join[analysis_dt,
                                 roll = 4 * 3600,
                                 on = .(hospitalization_id, sf_dttm = recorded_dttm_num)]
analysis_with_sf[, sf_dttm := NULL]

# --- Join pre-computed PF ratio to IMV timepoints (most recent prior within 4h) ---
pf_for_join <- as.data.table(pf_joined %>% select(hospitalization_id, pf_dttm, pf_ratio))
setkey(pf_for_join, hospitalization_id, pf_dttm)

# Need numeric key for the join
analysis_with_sf[, recorded_dttm_num := as.numeric(recorded_dttm)]
setkey(analysis_with_sf, hospitalization_id, recorded_dttm_num)

analysis_with_sf_pf <- pf_for_join[analysis_with_sf,
                                     roll = 4 * 3600,
                                     on = .(hospitalization_id, pf_dttm = recorded_dttm_num)]
analysis_with_sf_pf[, pf_dttm := NULL]

analysis_with_sf <- as_tibble(analysis_with_sf_pf)

message("SF ratio available at ", sum(!is.na(analysis_with_sf$sf_ratio)),
        " / ", nrow(analysis_with_sf), " timepoints")
message("PF ratio available at ", sum(!is.na(analysis_with_sf$pf_ratio)),
        " / ", nrow(analysis_with_sf), " timepoints")

# --- Join creatinine to IMV timepoints (most recent within 4h) ---
creat_labs <- cohort_labs %>%
  filter(lab_category == "creatinine", !is.na(lab_value_numeric)) %>%
  select(hospitalization_id, creat_dttm = lab_result_dttm,
         creatinine = lab_value_numeric)

creat_for_join <- as.data.table(creat_labs)
setkey(creat_for_join, hospitalization_id, creat_dttm)

analysis_creat_dt <- as.data.table(analysis_with_sf)
analysis_creat_dt[, recorded_dttm_num := as.numeric(recorded_dttm)]
setkey(analysis_creat_dt, hospitalization_id, recorded_dttm_num)

analysis_with_creat <- creat_for_join[analysis_creat_dt,
                                       roll = 4 * 3600,
                                       on = .(hospitalization_id, creat_dttm = recorded_dttm_num)]
analysis_with_creat[, creat_dttm := NULL]

# Compute baseline creatinine (first value per patient) and delta
analysis_with_creat[, creatinine_baseline := creatinine[which(!is.na(creatinine))[1]],
                     by = hospitalization_id]
analysis_with_creat[, delta_creatinine := creatinine - creatinine_baseline]

analysis_with_sf <- as_tibble(analysis_with_creat)

message("Creatinine available at ", sum(!is.na(analysis_with_sf$creatinine)),
        " / ", nrow(analysis_with_sf), " timepoints")
message("Delta creatinine available at ", sum(!is.na(analysis_with_sf$delta_creatinine)),
        " / ", nrow(analysis_with_sf), " timepoints")

# --- Join platelet count to IMV timepoints (most recent within 4h) ---
platelet_labs <- cohort_labs %>%
  filter(lab_category == "platelet_count", !is.na(lab_value_numeric)) %>%
  select(hospitalization_id, platelet_dttm = lab_result_dttm,
         platelet_count = lab_value_numeric)

platelet_for_join <- as.data.table(platelet_labs)
setkey(platelet_for_join, hospitalization_id, platelet_dttm)

analysis_platelet_dt <- as.data.table(analysis_with_sf)
analysis_platelet_dt[, recorded_dttm_num := as.numeric(recorded_dttm)]
setkey(analysis_platelet_dt, hospitalization_id, recorded_dttm_num)

analysis_with_platelet <- platelet_for_join[analysis_platelet_dt,
                                             roll = 4 * 3600,
                                             on = .(hospitalization_id,
                                                    platelet_dttm = recorded_dttm_num)]
analysis_with_platelet[, platelet_dttm := NULL]

analysis_with_platelet[, platelet_baseline := platelet_count[which(!is.na(platelet_count))[1]],
                        by = hospitalization_id]
analysis_with_platelet[, delta_platelet := platelet_count - platelet_baseline]

analysis_with_sf <- as_tibble(analysis_with_platelet)

message("Platelet count available at ", sum(!is.na(analysis_with_sf$platelet_count)),
        " / ", nrow(analysis_with_sf), " timepoints")
message("Delta platelet available at ", sum(!is.na(analysis_with_sf$delta_platelet)),
        " / ", nrow(analysis_with_sf), " timepoints")

# Join daily SOFA scores by hospitalization_id + day
analysis_with_sf <- analysis_with_sf %>%
  left_join(
    cohort_demographics %>% select(hospitalization_id, admission_dttm),
    by = "hospitalization_id"
  ) %>%
  mutate(
    sofa_day = as.integer(floor(as.numeric(
      difftime(recorded_dttm, admission_dttm, units = "days")
    )))
  ) %>%
  left_join(
    sofa_daily %>% select(hospitalization_id, sofa_day, sofa_total,
                           sofa_cv_97, sofa_coag, sofa_liver,
                           sofa_resp, sofa_cns, sofa_renal),
    by = c("hospitalization_id", "sofa_day")
  )

# Join demographics for survival time
analysis_with_sf <- analysis_with_sf %>%
  left_join(
    cohort_demographics %>%
      select(hospitalization_id, discharge_dttm, death_dttm),
    by = "hospitalization_id"
  ) %>%
  mutate(
    # Time from admission to death (any cause), using the patient-level death_dttm,
    # which captures out-of-hospital (post-discharge) deaths as well as in-hospital
    # deaths. discharge_dttm is deliberately NOT used as the survival endpoint:
    # censoring survivors at hospital discharge discards known post-discharge vital
    # status and is what caused the KM curves to censor most patients.
    death_day = as.numeric(difftime(death_dttm, admission_dttm, units = "days")),
    # All-cause mortality within the 60-day horizon. Patients with no recorded death
    # (or a death after day 60) are alive at the horizon and censored at day 60.
    # This assumes complete vital-status ascertainment to 60 days from the death
    # registry linkage (out-of-hospital deaths captured => no competing risk).
    mortality_event_60 = if_else(
      !is.na(death_day) & death_day >= 0 & death_day <= 60, 1L, 0L
    ),
    surv_time = if_else(mortality_event_60 == 1L, death_day, 60)
  )

# =============================================================================
# 3e. Cross-sectional cohort selection
# =============================================================================
# Cross-sectional design: for each patient, take the FIRST IMV timepoint at
# which all core analysis variables are observed ("first timepoint with all
# available data"). A patient is included only if, AT that index timepoint,
# they were receiving lung-protective tidal volumes (VT/PBW 6-8 mL/kg) AND
# were hypoxemic (SF ratio < 315). Longitudinal / repeated-measures analyses
# are deferred to a future project.

# SF 315 is the Rice-equivalent of PF 300 (mild ARDS / hypoxemia threshold).
SF_HYPOXEMIA_THRESHOLD <- 315

# A timepoint has "all available data" when the core analysis variables used in
# the cross-sectional models are all observed: VT/PBW, VT/PFVC, SF ratio, and a
# daily SOFA score. Driving-pressure-derived measures (dp, crs) are intentionally
# NOT required here, since plateau pressure is frequently unrecorded and would
# otherwise shrink the cohort dramatically.
analysis_with_completeness <- analysis_with_sf %>%
  mutate(
    has_all_data = !is.na(vtpbw) & !is.na(vtpfvc) &
      !is.na(sf_ratio) & !is.na(sofa_total)
  )

# Apply the inclusion predicates BEFORE reducing to one row per patient, matching
# the original analysis (filter qualifying timepoints, then distinct(subject_id)).
# A patient is included if ANY complete-data IMV timepoint is simultaneously
# lung-protective (VT/PBW 6-8) and hypoxemic (SF < threshold); the patient's FIRST
# such qualifying timepoint becomes the index row. The previous behaviour required
# each patient's *first* complete-data timepoint to itself meet the criteria, which
# discarded patients whose later timepoints qualified and shrank the cohort
# relative to the original.
qualifying_timepoints <- analysis_with_completeness %>%
  filter(has_all_data, vtpbw >= 6, vtpbw <= 8, sf_ratio < SF_HYPOXEMIA_THRESHOLD)

message("Patients with >=1 complete-data IMV timepoint: ",
        n_distinct(analysis_with_completeness$hospitalization_id[analysis_with_completeness$has_all_data]))

# Index timepoint = each qualifying patient's FIRST qualifying IMV timepoint
cross_sectional <- qualifying_timepoints %>%
  group_by(hospitalization_id) %>%
  slice_min(recorded_dttm, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # Attach weight and BMI (needed as a covariate in the driving-pressure /
  # elastance demographic-bias models in script 04). Weight is extracted and
  # imputed in script 01 but is not otherwise carried downstream.
  left_join(cohort_weights, by = "hospitalization_id") %>%
  mutate(bmi = if_else(!is.na(weight_kg) & height_cm > 0,
                       weight_kg / (height_cm / 100)^2, NA_real_))

eligible_patients <- cross_sectional$hospitalization_id

message("Included patients (VT/PBW 6-8 AND SF<", SF_HYPOXEMIA_THRESHOLD,
        " at any complete-data timepoint): ", length(eligible_patients))

# All IMV timepoints for the included patients, retained for descriptive
# summaries. Repeated-measures / longitudinal modeling is deferred to a
# future project, so no per-timepoint eligibility flags are computed here.
analysis_all <- analysis_with_sf %>%
  filter(hospitalization_id %in% eligible_patients)

# =============================================================================
# 3f. 28-day ventilator-free days (VFDs)
# =============================================================================

# Count total hours on IMV per patient within first 28 days
imv_hours_28d <- resp_waterfall %>%
  filter(hospitalization_id %in% eligible_patients,
         tolower(device_category) == "imv") %>%
  inner_join(
    cohort_demographics %>% select(hospitalization_id, admission_dttm),
    by = "hospitalization_id"
  ) %>%
  mutate(t_days = as.numeric(difftime(recorded_dttm, admission_dttm, units = "days"))) %>%
  filter(t_days >= 0, t_days <= 28) %>%
  group_by(hospitalization_id) %>%
  summarise(imv_hours = n(), .groups = "drop") %>%
  mutate(imv_days = imv_hours / 24)

# VFD-28: 28 minus ventilator days, deaths within 28 days get 0.
# Death is all-cause (in- or out-of-hospital) within 28 days, from the patient-level
# death_dttm — consistent with the 60-day survival endpoint. The previous in-hospital
# -only condition (deceased == 1) wrongly credited ventilator-free days to patients
# who died out of hospital within 28 days.
vfd_data <- cohort_demographics %>%
  filter(hospitalization_id %in% eligible_patients) %>%
  select(hospitalization_id, deceased, death_dttm, admission_dttm) %>%
  left_join(imv_hours_28d, by = "hospitalization_id") %>%
  mutate(
    imv_days = replace_na(imv_days, 0),
    death_day = as.numeric(difftime(death_dttm, admission_dttm, units = "days")),
    died_within_28 = !is.na(death_day) & death_day >= 0 & death_day <= 28,
    vfd_28 = if_else(died_within_28, 0, pmax(28 - imv_days, 0))
  ) %>%
  select(hospitalization_id, vfd_28)

# Attach VFDs to the cross-sectional cohort (one row per included patient)
cross_sectional <- cross_sectional %>%
  left_join(vfd_data, by = "hospitalization_id")

message("28-day VFDs computed. Median VFD-28: ",
        round(median(vfd_data$vfd_28, na.rm = TRUE), 1))

# =============================================================================
# 3g. VFR terciles
# =============================================================================

cross_sectional <- cross_sectional %>%
  mutate(
    vfr_tercile = ntile(vtpfvc, 3),
    vfr_tercile = factor(vfr_tercile, labels = c("T1 (Low)", "T2 (Mid)", "T3 (High)"))
  )

# =============================================================================
# 3h. Save outputs
# =============================================================================

write_parquet(analysis_all, file.path(output_dir, "analysis_all_timepoints.parquet"))
write_parquet(cross_sectional, file.path(output_dir, "analysis_cross_sectional.parquet"))
write_parquet(sofa_scores, file.path(output_dir, "sofa_scores.parquet"))
write_parquet(sofa_daily, file.path(output_dir, "sofa_daily.parquet"))

# =============================================================================
# 3i. Complete the attrition log (steps 4-7) and write the full CONSORT table
# =============================================================================
# Patient-level, monotonic counts continuing the funnel from script 01.
n_step4 <- n_distinct(pbw_pfvc_data$hospitalization_id)
n_step5 <- n_distinct(
  analysis_with_completeness$hospitalization_id[analysis_with_completeness$has_all_data]
)
n_step6 <- analysis_with_completeness %>%
  filter(has_all_data, vtpbw >= 6, vtpbw <= 8) %>%
  summarise(n = n_distinct(hospitalization_id)) %>% pull(n)
n_step7 <- length(eligible_patients)

partial_path <- file.path(output_dir, "attrition_log_partial.csv")
if (!file.exists(partial_path)) {
  stop("Missing attrition_log_partial.csv from script 01: ", partial_path)
}
attrition <- read_csv(partial_path, show_col_types = FALSE) %>%
  attrition_add(ATTRITION_STEPS[4], n_step4,
                exclusion_reason = "Height outside 150-210 cm or PBW/PFVC missing") %>%
  attrition_add(ATTRITION_STEPS[5], n_step5,
                exclusion_reason = "Incomplete index data (VT/PBW, VT/PFVC, SF, SOFA)") %>%
  attrition_add(ATTRITION_STEPS[6], n_step6,
                exclusion_reason = "Not lung-protective (VT/PBW outside 6-8)") %>%
  attrition_add(ATTRITION_STEPS[7], n_step7,
                exclusion_reason = "Not hypoxemic (SF ratio >= 315)") %>%
  mutate(site = site_name, .before = 1)

write_csv(attrition, file.path(final_dir, paste0("attrition_log_", site_name, ".csv")))
message("Attrition log written (7 steps): ",
        paste(attrition$n_remaining, collapse = " -> "))

# =============================================================================
# 3j. Federated PBW:PFVC distribution exports (site-specific; poolable)
# =============================================================================
# Aggregated summaries of the PBW:PFVC ratio by demographic group only -- no
# row-level data leaves the site. (a) histograms on fixed bins (summable across
# sites into a pooled distribution); (b) quantile summaries for exact per-site
# boxplots. Groups with n < DIST_MIN_CELL are dropped entirely (CLAUDE.md n>=10).
dist_groups <- bind_rows(
  cross_sectional %>% transmute(group_type = "sex",
    group_value = as.character(sex_category), value = pbwpfvc),
  cross_sectional %>% transmute(group_type = "race",
    group_value = as.character(race_category), value = pbwpfvc),
  cross_sectional %>%
    mutate(gv = cut(age_at_admission, AGE_BIN_EDGES, right = FALSE)) %>%
    transmute(group_type = "age_bin", group_value = as.character(gv), value = pbwpfvc),
  cross_sectional %>%
    mutate(gv = cut(height_cm, HEIGHT_BIN_EDGES, right = FALSE)) %>%
    transmute(group_type = "height_bin", group_value = as.character(gv), value = pbwpfvc)
) %>%
  filter(!is.na(value), !is.na(group_value), group_value != "NA")

dist_keep <- dist_groups %>% count(group_type, group_value) %>%
  filter(n >= DIST_MIN_CELL)
dist_groups <- dist_groups %>% semi_join(dist_keep, by = c("group_type", "group_value"))

# (a) histograms with tails clamped into the edge bins so per-group totals == N
.bin_lo <- min(PBWPFVC_BIN_EDGES); .bin_hi <- max(PBWPFVC_BIN_EDGES)
dist_histograms <- dist_groups %>%
  mutate(value = pmin(pmax(value, .bin_lo), .bin_hi - 1e-9),
         bin_i = findInterval(value, PBWPFVC_BIN_EDGES, rightmost.closed = TRUE)) %>%
  count(group_type, group_value, bin_i, name = "count") %>%
  mutate(bin_left  = PBWPFVC_BIN_EDGES[bin_i],
         bin_right = PBWPFVC_BIN_EDGES[bin_i + 1]) %>%
  transmute(site = site_name, group_type, group_value, bin_left, bin_right, count) %>%
  arrange(group_type, group_value, bin_left)

# (b) per-group quantile summaries
dist_quantiles <- dist_groups %>%
  group_by(group_type, group_value) %>%
  summarise(
    n = n(), median = median(value),
    q1 = quantile(value, 0.25), q3 = quantile(value, 0.75),
    p10 = quantile(value, 0.10), p20 = quantile(value, 0.20),
    p30 = quantile(value, 0.30), p40 = quantile(value, 0.40),
    p50 = quantile(value, 0.50), p60 = quantile(value, 0.60),
    p70 = quantile(value, 0.70), p80 = quantile(value, 0.80),
    p90 = quantile(value, 0.90),
    .groups = "drop"
  ) %>%
  mutate(site = site_name, .before = 1)

write_csv(dist_histograms, file.path(final_dir, paste0("dist_histograms_", site_name, ".csv")))
write_csv(dist_quantiles, file.path(final_dir, paste0("dist_quantiles_", site_name, ".csv")))
message("Federated distribution exports written for ",
        n_distinct(dist_groups$group_value), " demographic groups")

message("Script 03 complete.")
message("  Cross-sectional (1 row/included patient): ", nrow(cross_sectional))
message("  All timepoints (included patients): ", nrow(analysis_all))
