# =============================================================================
# Script 01: Cohort Identification
# =============================================================================
rm(list= ls())
library(tidyverse)
library(arrow)
library(here)
library(data.table)
library(lubridate)

# Load config
source("utils/config.R")
source("utils/attrition_log.R")

site_name <- config$site_name
tables_path <- config$tables_path
file_type <- config$file_type

message("Site: ", site_name)
message("Tables path: ", tables_path)
message("File type: ", file_type)

# Expand ~ in path
tables_path <- path.expand(tables_path)

# Fail loudly if the data DIRECTORY itself is unreachable (e.g. the remote drive is not
# mounted). Without this, a missing mount can surface as a confusing per-file "Missing
# table" error -- or, worse, a stale/empty mount passes every check and the pipeline runs
# silently on empty data (see the disconnected-drive incident).
if (!dir.exists(tables_path)) {
  stop("CLIF data path does not exist: '", tables_path,
       "'. Is the remote drive mounted? Check config$tables_path.")
}

# =============================================================================
# Load CLIF tables
# =============================================================================

table_names <- c("patient", "hospitalization", "adt", "respiratory_support",
                 "vitals", "labs", "medication_admin_continuous",
                 "patient_assessments")

# Check that all files are present before loading
for (tbl in table_names) {
  fpath <- file.path(tables_path, paste0("clif_", tbl, ".", file_type))
  if (!file.exists(fpath)) {
    stop("Missing table: ", fpath)
  }
}

# Loader. For PARQUET, open each table as a lazy Arrow dataset so column/row
# predicates push down to the file scan -- only the rows we keep are ever
# materialized, and collect() realizes the query. CSV/FST have no predicate
# pushdown, so they read in full and filter in memory (collect() is then a
# no-op on the local frame). One code path serves all three formats because
# collect() is a no-op on a local data frame.
open_clif <- function(tbl) {
  fpath <- file.path(tables_path, paste0("clif_", tbl, ".", file_type))
  if (file_type == "parquet") return(arrow::open_dataset(fpath))
  if (file_type == "csv")     return(readr::read_csv(fpath, show_col_types = FALSE))
  if (file_type == "fst")     return(fst::read_fst(fpath))
  stop("Unsupported file_type: ", file_type)
}

# Category whitelists for the big event tables. Defined once and used BOTH for the
# load-time predicate pushdown here AND the downstream extraction filters below.
# Pre-filtering at load is safe: each big table is consumed only within its
# whitelist. ph_arterial/ph_venous feed the script-10 [T5b] pH sensitivity.
vitals_categories_needed     <- c("height_cm", "weight_kg", "spo2", "map")
med_categories_needed        <- c("norepinephrine", "epinephrine", "vasopressin",
                                  "dopamine", "phenylephrine", "dobutamine")
lab_categories_needed        <- c("po2_arterial", "pco2_arterial", "creatinine",
                                  "bilirubin_total", "platelet_count",
                                  "ph_arterial", "ph_venous")
assessment_categories_needed <- c("gcs_total")

# Eligibility tables: read in full (used wholesale to derive the cohort).
clif_patient             <- open_clif("patient") %>% collect()
clif_hospitalization     <- open_clif("hospitalization") %>% collect()
clif_adt                 <- open_clif("adt") %>% collect()
clif_respiratory_support <- open_clif("respiratory_support") %>% collect()

# Big event tables: push the category predicate down (parquet), then materialize.
clif_vitals      <- open_clif("vitals") %>%
  filter(vital_category %in% vitals_categories_needed) %>% collect()
clif_labs        <- open_clif("labs") %>%
  filter(lab_category %in% lab_categories_needed) %>% collect()
clif_meds        <- open_clif("medication_admin_continuous") %>%
  filter(med_category %in% med_categories_needed) %>% collect()
clif_assessments <- open_clif("patient_assessments") %>%
  filter(assessment_category %in% assessment_categories_needed) %>% collect()
# Continuous renal replacement therapy (CLIF 2.1 crrt_therapy): the biotrauma
# joint models censor the creatinine trajectory at the first CRRT record, so the
# panel needs the start time. Only the identifiers and the mode category are kept.
# The table is optional in CLIF and absent at some sites (MIMIC among them): then
# an empty table is written, the creatinine trajectory runs uncensored, and the
# joint-model panel summary records crrt_available = FALSE. This is a site
# characteristic reported in the outputs, not a silent fallback.
crrt_path <- file.path(tables_path, paste0("clif_crrt_therapy.", file_type))
crrt_available <- file.exists(crrt_path)
clif_crrt <- if (crrt_available) {
  open_clif("crrt_therapy") %>%
    select(hospitalization_id, recorded_dttm, crrt_mode_category) %>% collect()
} else {
  message("*** No crrt_therapy table at ", tables_path,
          ": RRT censoring of the creatinine trajectory is unavailable at this site. ***")
  # zero rows, with hospitalization_id of the site's own type so later joins agree
  clif_hospitalization %>% slice(0) %>% select(hospitalization_id) %>%
    mutate(recorded_dttm = as.POSIXct(character()), crrt_mode_category = character())
}

message("Loaded: patient=", nrow(clif_patient), " hosp=", nrow(clif_hospitalization),
        " adt=", nrow(clif_adt), " resp=", nrow(clif_respiratory_support))
message("Loaded (category-filtered): vitals=", nrow(clif_vitals), " labs=", nrow(clif_labs),
        " meds=", nrow(clif_meds), " assessments=", nrow(clif_assessments),
        " crrt=", nrow(clif_crrt))

# Fail LOUDLY if any CORE table loaded empty. arrow::open_dataset() does NOT error on an
# unreachable/stale parquet path -- it opens a 0-row dataset and collect() succeeds -- so
# a disconnected drive otherwise runs the whole pipeline on empty data and silently
# reproduces stale outputs. The core eligibility tables are non-empty for any real CLIF
# site, so a 0-row read here is a hard error, not a site characteristic. (labs/meds/
# assessments may be legitimately sparse, so they are reported above but not gated here.)
.core_counts <- c(patient = nrow(clif_patient), hospitalization = nrow(clif_hospitalization),
                  adt = nrow(clif_adt), respiratory_support = nrow(clif_respiratory_support),
                  vitals = nrow(clif_vitals))
if (any(.core_counts == 0L)) {
  stop("Empty CORE CLIF table(s) after load: ",
       paste(names(.core_counts)[.core_counts == 0L], collapse = ", "),
       ". A 0-row read almost always means the data path '", tables_path,
       "' is unreachable or stale (remote drive not mounted?) -- arrow opens an empty ",
       "dataset without erroring. Verify the mount, then rerun.")
}

# =============================================================================
# Cohort filtering
# =============================================================================

# Pull the last recorded hospitalization where the patient received IMV with VC/AC in the ICU

# Identify hospitalizations that received invasive ventilation with a recorded
# set tidal volume.
#
# The original analysis (PBWvsFVC) imposed NO ventilator-mode requirement: a
# patient was eligible if a *set* tidal volume was ever recorded (MIMIC itemid
# 224684). We mirror that here using tidal_volume_set, rather than requiring
# mode_category to map exactly to "assist control-volume control". In CLIF,
# mode_category is recorded only at mode-change events (far sparser than the set
# tidal volume) and depends on ETL-specific string mapping, so the strict mode
# filter dropped many ventilated patients the original analysis retained. The
# downstream VT/PBW 6-8 gate restricts to volume-targeted breaths, matching the
# original's ccperkg 6-8 criterion.
imv_ids <- clif_respiratory_support %>%
  mutate(tidal_volume_set_numeric = suppressWarnings(as.numeric(tidal_volume_set))) %>%
  filter(!is.na(tidal_volume_set_numeric), tidal_volume_set_numeric > 0) %>%
  distinct(hospitalization_id) %>%
  pull(hospitalization_id)

icu_ids <- clif_adt %>%
    filter(tolower(location_category) %in% c("icu")) %>% # in the ICU 
    distinct(hospitalization_id) %>%
    pull(hospitalization_id)

cohort_patient_and_hospitalization_ids <- clif_hospitalization %>%
  filter(age_at_admission >= 18) %>% # Only adults
  filter(hospitalization_id %in% imv_ids) %>% # who received VC ventilation
  filter(hospitalization_id %in% icu_ids) %>% # in the ICU 
  arrange(desc(admission_dttm)) %>% # if multiple hospitalizations, we want the last admission
  distinct(patient_id, .keep_all = T) %>%
  dplyr::select(patient_id, hospitalization_id)

eligible_patients <- cohort_patient_and_hospitalization_ids$patient_id
eligible_hospitalizations <- cohort_patient_and_hospitalization_ids$hospitalization_id

message("Cohort size: ", length(eligible_hospitalizations), " hospitalizations")

# =============================================================================
# Attrition log (steps 1-3): the cohort funnel down to the index hospitalization
# =============================================================================
# Counts are distinct PATIENTS surviving each filter applied cumulatively to the
# SAME hospitalization, so the funnel is monotonic and matches the per-patient
# analytic cohort (steps 4-7 are appended in script 03). Patient-level counts
# (not hospitalization-level) keep the chain consistent with downstream outputs.
funnel <- clif_hospitalization %>% filter(age_at_admission >= 18)
n_adult <- n_distinct(funnel$patient_id)

funnel <- funnel %>% filter(hospitalization_id %in% icu_ids)
n_icu <- n_distinct(funnel$patient_id)

funnel <- funnel %>% filter(hospitalization_id %in% imv_ids)
n_imv <- n_distinct(funnel$patient_id)   # == length(eligible_patients)

attrition <- attrition_init() %>%
  attrition_add(ATTRITION_STEPS[1], n_adult) %>%
  attrition_add(ATTRITION_STEPS[2], n_icu,
                exclusion_reason = "No ICU admission") %>%
  attrition_add(ATTRITION_STEPS[3], n_imv,
                exclusion_reason = "No invasive ventilation with set tidal volume")

message("Attrition (steps 1-3): adults=", n_adult, ", +ICU=", n_icu, ", +IMV=", n_imv)

# =============================================================================
# Run respiratory support waterfall
# =============================================================================

source("utils/process_resp_support_waterfall.R")

resp_support_cohort <- clif_respiratory_support %>%
  filter(hospitalization_id %in% eligible_hospitalizations)

resp_waterfall <- process_resp_support_waterfall(resp_support_cohort)

# =============================================================================
# Extract height
# =============================================================================

#Impute heights for patients with missing height data
all_hospitalizations_for_heights <- 
  clif_hospitalization %>%
  filter(patient_id %in% eligible_patients) %>%
  distinct(patient_id, hospitalization_id)


cohort_heights <- clif_vitals %>%
  filter(vital_category == "height_cm") %>%
  filter(!is.na(vital_value)) %>%
  mutate(height_cm = as.numeric(vital_value)) %>%
  dplyr::select(hospitalization_id, height_cm) %>%
  summarize(height_cm = mean(height_cm), .by = hospitalization_id) %>%
  right_join(all_hospitalizations_for_heights) %>%
  full_join(cohort_patient_and_hospitalization_ids %>% dplyr::rename(eligible_hosp_id = hospitalization_id)) %>%
  arrange(patient_id) %>% 
  mutate(median_height = median(height_cm, na.rm = T), .by = patient_id) %>%
  mutate(height_cm = case_when(
    hospitalization_id == eligible_hosp_id & !is.na(height_cm) ~ height_cm, #If height measured, use height from that admission
    hospitalization_id == eligible_hosp_id & is.na(height_cm) ~ median_height, #If height unavailable, use median of all available heights
    TRUE ~ NA
  )) %>%
  filter(!is.na(height_cm)) %>%
  dplyr::select(hospitalization_id, height_cm)


message("Heights extracted: ", nrow(cohort_heights), " hospitalizations")

# =============================================================================
# Extract weights (for vasopressor dose standardization to mcg/kg/min)
# =============================================================================

# Use mean recorded weight per hospitalization (sanity-bounded to 30-1100 kg).
# Weight is used downstream to convert norepinephrine doses reported in mcg/min
# to mcg/kg/min. Patient-level median is used when no admission weight is
# recorded, matching the height imputation pattern.
cohort_weights <- clif_vitals %>%
  filter(vital_category == "weight_kg", !is.na(vital_value)) %>%
  mutate(weight_kg = as.numeric(vital_value)) %>%
  filter(weight_kg >= 30, weight_kg <= 1100) %>%
  dplyr::select(hospitalization_id, weight_kg) %>%
  summarize(weight_kg = mean(weight_kg, na.rm = TRUE), .by = hospitalization_id) %>%
  right_join(all_hospitalizations_for_heights, by = "hospitalization_id") %>%
  full_join(cohort_patient_and_hospitalization_ids %>%
              dplyr::rename(eligible_hosp_id = hospitalization_id),
            by = "patient_id") %>%
  arrange(patient_id) %>%
  mutate(median_weight = median(weight_kg, na.rm = TRUE), .by = patient_id) %>%
  mutate(weight_kg = case_when(
    hospitalization_id == eligible_hosp_id & !is.na(weight_kg) ~ weight_kg,
    hospitalization_id == eligible_hosp_id & is.na(weight_kg) ~ median_weight,
    TRUE ~ NA_real_
  )) %>%
  filter(hospitalization_id %in% eligible_hospitalizations) %>%
  dplyr::select(hospitalization_id, weight_kg) %>%
  distinct(hospitalization_id, .keep_all = TRUE)

message("Weights extracted: ", sum(!is.na(cohort_weights$weight_kg)),
        " of ", nrow(cohort_weights), " hospitalizations")

# =============================================================================
# Extract SpO2
# =============================================================================

cohort_spo2 <- clif_vitals %>%
  filter(hospitalization_id %in% eligible_hospitalizations,
         vital_category == "spo2") %>%
  mutate(vital_value = as.numeric(vital_value)) %>%
  select(hospitalization_id, recorded_dttm, vital_value) %>%
  rename(spo2_value = vital_value)
  # SpO2 is NO LONGER capped at <= 97 here. The SF-validity bounds (80-97) live where the
  # SF ratio is actually formed -- script 03 (filter >= 80 & <= 97 for the cross-sectional
  # index SF) and script 10 (SpO2 clamped to [80,97] for the daily worst SF). Capping here
  # poisoned the SHARED cohort_vitals intermediate: well-oxygenated patient-days (SpO2
  # always >= 98) lost ALL their SpO2 rows, so script 10's longitudinal panel read them as
  # "no SpO2 charted" and dropped them. The QC outlier threshold (spo2 50-100) now governs
  # the upper bound in cohort_vitals_clean; downstream SF filtering is unchanged.

# =============================================================================
# Extract MAP
# =============================================================================

cohort_map <- clif_vitals %>%
  filter(hospitalization_id %in% eligible_hospitalizations,
         vital_category == "map") %>%
  mutate(vital_value = as.numeric(vital_value)) %>%
  select(hospitalization_id, recorded_dttm, vital_value) %>%
  rename(map_value = vital_value)

# Combine vitals (SpO2 + MAP)
cohort_vitals <- bind_rows(
  cohort_spo2 %>% mutate(vital_category = "spo2") %>% rename(vital_value = spo2_value),
  cohort_map %>% mutate(vital_category = "map") %>% rename(vital_value = map_value)
)

# =============================================================================
# Extract labs (PaO2, creatinine, bilirubin_total, platelets)
# =============================================================================

# lab_categories_needed (incl. ph_arterial/ph_venous for the script-10 [T5b]
# sensitivity) is defined at load above and already pushed down at read time;
# this filter is now a no-op safeguard on the in-memory frame.
cohort_labs <- clif_labs %>%
  filter(hospitalization_id %in% eligible_hospitalizations,
         lab_category %in% lab_categories_needed) %>%
  mutate(lab_value_numeric = as.numeric(lab_value))

message("Labs extracted: ", nrow(cohort_labs), " rows")

# =============================================================================
# Extract vasopressor meds (norepinephrine, vasopressin)
# =============================================================================

cohort_meds <- clif_meds %>%
  filter(hospitalization_id %in% eligible_hospitalizations,
         med_category %in% c("norepinephrine", "epinephrine", "vasopressin",
                              "dopamine", "phenylephrine", "dobutamine"))

# =============================================================================
# Extract GCS assessments
# =============================================================================

cohort_assessments <- clif_assessments %>%
  filter(hospitalization_id %in% eligible_hospitalizations,
         assessment_category %in% c("gcs_total"))

# =============================================================================
# Extract CRRT records
# =============================================================================

cohort_crrt <- clif_crrt %>%
  filter(hospitalization_id %in% eligible_hospitalizations, !is.na(recorded_dttm))
attr(cohort_crrt, "crrt_available") <- crrt_available

message("CRRT extracted: ", nrow(cohort_crrt), " rows, ",
        n_distinct(cohort_crrt$hospitalization_id), " hospitalizations",
        if (!crrt_available) " (table absent at this site)")

# =============================================================================
# Merge demographics
# =============================================================================

cohort_demographics <- clif_hospitalization %>%
  filter(hospitalization_id %in% eligible_hospitalizations) %>%
  distinct(hospitalization_id, .keep_all = T) %>%
  left_join(clif_patient, by = "patient_id") %>%
  mutate(
    # Harmonize CLIF race categories to GLI-compatible groups
    race_category = case_when(
      race_category == "White"                     ~ "WHITE",
      race_category == "Black or African American" ~ "BLACK",
      TRUE                                         ~ "OTHER"
    ),
    # In-hospital mortality 
    deceased = case_when(
      discharge_category == "Expired" ~ 1,
      TRUE ~ 0L
    )
  )

message("Demographics: ", nrow(cohort_demographics), " rows")
message("Mortality rate: ", round(mean(cohort_demographics$deceased) * 100, 1), "%")

# =============================================================================
# Negative-control cohorts: adult ICU patients OUTSIDE the analytic cohort
# =============================================================================
# The analytic cohort is hypoxemic AND ventilated with a set tidal volume, so any
# PBW/PFVC -> mortality association in it could run through tidal-volume dosing or
# through what the ratio's inputs (age, sex, race, height) index on their own. Two
# cohorts in which no lung-protective dosing decision was made test that: if the
# ratio (and height) predict death equally in (a) NON-hypoxemic, NON-ventilated ICU
# adults and (b) NON-hypoxemic VENTILATED adults, the association is not ventilatory.
# Script 03 derives PBW/PFVC for them, script 04 fits the models (section 4j).
#
# One row per patient (last adult ICU admission, as for the analytic cohort).
# Hypoxemia is classified from every SpO2 in the hospitalization: each SpO2 is
# joined to the most recent respiratory-support record within 4 h; FiO2 comes from
# fio2_set (fractions), from room-air device category (0.21), or from nasal-cannula
# flow (0.21 + 0.03 x L/min, capped at 0.60); a SpO2 with no support record within
# 4 h is taken as room air. Hypoxemic = any SF < 315 among SpO2 80-97 with a known
# FiO2, or any SpO2 < 80. Patients with no SpO2 at all cannot be classified and are
# dropped. This is deliberately lighter than the analytic cohort's waterfall (no
# hourly scaffold), because it runs on every adult ICU hospitalization.

nc_hosp <- clif_hospitalization %>%
  filter(age_at_admission >= 18, hospitalization_id %in% icu_ids) %>%
  arrange(desc(admission_dttm)) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  mutate(imv_set_vt = hospitalization_id %in% imv_ids)

# heights: mean in the index hospitalization, else the patient's median elsewhere
nc_height_all <- clif_vitals %>%
  filter(vital_category == "height_cm", !is.na(vital_value)) %>%
  mutate(height_cm = as.numeric(vital_value)) %>%
  summarize(height_cm = mean(height_cm), .by = hospitalization_id) %>%
  inner_join(clif_hospitalization %>% distinct(hospitalization_id, patient_id), by = "hospitalization_id")
nc_heights <- nc_hosp %>% select(patient_id, hospitalization_id) %>%
  left_join(nc_height_all %>% select(hospitalization_id, height_cm), by = "hospitalization_id") %>%
  left_join(nc_height_all %>% summarize(median_height = median(height_cm), .by = patient_id), by = "patient_id") %>%
  mutate(height_cm = coalesce(height_cm, median_height)) %>%
  select(hospitalization_id, height_cm)

# FiO2 timeline from the raw respiratory-support table (all NC hospitalizations)
nc_fio2 <- clif_respiratory_support %>%
  filter(hospitalization_id %in% nc_hosp$hospitalization_id) %>%
  mutate(device_category = tolower(device_category),
         fio2_num = suppressWarnings(as.numeric(fio2_set)),
         fio2_num = if_else(!is.na(fio2_num) & fio2_num > 1, fio2_num / 100, fio2_num),
         lpm_num  = suppressWarnings(as.numeric(lpm_set)),
         fio2_est = case_when(
           !is.na(fio2_num) & fio2_num >= 0.21 & fio2_num <= 1 ~ fio2_num,
           replace_na(device_category == "room air", FALSE)     ~ 0.21,
           replace_na(device_category == "nasal cannula", FALSE) & !is.na(lpm_num) ~ pmin(0.21 + 0.03 * lpm_num, 0.60),
           TRUE ~ NA_real_)) %>%
  filter(!is.na(fio2_est), !is.na(recorded_dttm)) %>%
  transmute(hospitalization_id, fio2_dttm = as.numeric(recorded_dttm), fio2_est)
nc_spo2 <- clif_vitals %>%
  filter(hospitalization_id %in% nc_hosp$hospitalization_id, vital_category == "spo2") %>%
  mutate(spo2 = as.numeric(vital_value)) %>%
  filter(!is.na(spo2), spo2 >= 50, spo2 <= 100, !is.na(recorded_dttm)) %>%
  transmute(hospitalization_id, spo2_dttm = as.numeric(recorded_dttm), spo2)
nc_fio2_dt <- as.data.table(nc_fio2); setkey(nc_fio2_dt, hospitalization_id, fio2_dttm)
nc_spo2_dt <- as.data.table(nc_spo2); setkey(nc_spo2_dt, hospitalization_id, spo2_dttm)
nc_sf <- nc_fio2_dt[nc_spo2_dt, roll = 4 * 3600, on = .(hospitalization_id, fio2_dttm = spo2_dttm)] %>%
  as_tibble() %>%
  mutate(fio2_est = coalesce(fio2_est, 0.21),          # no support record within 4 h = room air
         sf = spo2 / fio2_est,
         hypox_obs = spo2 < 80 | (spo2 <= 97 & sf < 315)) %>%
  summarize(n_spo2 = n(), hypoxemic = any(hypox_obs), min_sf = min(sf), .by = hospitalization_id)

# The ventilated control is valid only if the patient stayed non-hypoxemic for the
# WHOLE ventilated period, and the SpO2-based rule above cannot see hypoxemia masked
# by a high FiO2 (SpO2 99% on FiO2 0.60 yields no computable SF). So, over the IMV
# rows of the waterfall (these patients are in the IMV cohort): every SpO2 must give
# SF >= 315, no PaO2/FiO2 may fall below 300, and FiO2 may never exceed 0.40. IMV
# without any documented FiO2 is unclassifiable (NA) and excluded from the control.
imv_rows <- resp_waterfall %>%
  filter(device_category == "imv", !is.na(recorded_dttm)) %>%
  transmute(hospitalization_id, t = as.numeric(recorded_dttm), fio2_set)
imv_span <- imv_rows %>% summarize(t_start = min(t), t_end = max(t),
                                   fio2_documented = any(!is.na(fio2_set)),
                                   max_fio2 = suppressWarnings(max(fio2_set, na.rm = TRUE)), .by = hospitalization_id)
imv_fio2_dt <- imv_rows %>% filter(!is.na(fio2_set)) %>% rename(fio2_t = t) %>% as.data.table()
setkey(imv_fio2_dt, hospitalization_id, fio2_t)
vent_spo2 <- cohort_spo2 %>%
  transmute(hospitalization_id, spo2_t = as.numeric(recorded_dttm), spo2 = spo2_value) %>%
  filter(!is.na(spo2), spo2 >= 50, spo2 <= 100) %>%
  inner_join(imv_span, by = "hospitalization_id") %>%
  filter(spo2_t >= t_start, spo2_t <= t_end) %>%
  select(hospitalization_id, spo2_t, spo2) %>% as.data.table()
setkey(vent_spo2, hospitalization_id, spo2_t)
vent_spo2 <- imv_fio2_dt[vent_spo2, roll = 4 * 3600, on = .(hospitalization_id, fio2_t = spo2_t)] %>%
  as_tibble() %>%
  summarize(spo2_hypox = any(spo2 < 80 | (spo2 <= 97 & !is.na(fio2_set) & spo2 / fio2_set < 315)),
            .by = hospitalization_id)
vent_pao2 <- clif_labs %>%
  filter(lab_category == "po2_arterial", hospitalization_id %in% imv_span$hospitalization_id) %>%
  transmute(hospitalization_id, pao2_t = as.numeric(lab_result_dttm), pao2 = as.numeric(lab_value)) %>%
  filter(!is.na(pao2), !is.na(pao2_t)) %>%
  inner_join(imv_span, by = "hospitalization_id") %>%
  filter(pao2_t >= t_start, pao2_t <= t_end) %>%
  select(hospitalization_id, pao2_t, pao2) %>% as.data.table()
setkey(vent_pao2, hospitalization_id, pao2_t)
vent_pao2 <- imv_fio2_dt[vent_pao2, roll = 4 * 3600, on = .(hospitalization_id, fio2_t = pao2_t)] %>%
  as_tibble() %>%
  summarize(pf_hypox = any(!is.na(fio2_set) & pao2 / fio2_set < 300), .by = hospitalization_id)
vent_hypox <- imv_span %>%
  left_join(vent_spo2, by = "hospitalization_id") %>%
  left_join(vent_pao2, by = "hospitalization_id") %>%
  transmute(hospitalization_id, vent_fio2_documented = fio2_documented,
            hypoxemic_during_vent = case_when(
              !fio2_documented ~ NA,
              TRUE ~ coalesce(spo2_hypox, FALSE) | coalesce(pf_hypox, FALSE) | max_fio2 > 0.40))

nc_cohort <- nc_hosp %>%
  left_join(clif_patient, by = "patient_id") %>%
  mutate(race_category = case_when(race_category == "White" ~ "WHITE",
                                   race_category == "Black or African American" ~ "BLACK",
                                   TRUE ~ "OTHER"),
         deceased = if_else(discharge_category == "Expired", 1L, 0L, missing = 0L)) %>%
  inner_join(nc_sf, by = "hospitalization_id") %>%        # drops patients with no SpO2
  left_join(nc_heights, by = "hospitalization_id") %>%
  left_join(vent_hypox, by = "hospitalization_id") %>%
  select(hospitalization_id, patient_id, age_at_admission, sex_category, race_category,
         height_cm, admission_dttm, discharge_dttm, death_dttm, deceased,
         imv_set_vt, hypoxemic, hypoxemic_during_vent, vent_fio2_documented, n_spo2, min_sf)
message("Negative-control frame: ", nrow(nc_cohort), " adult ICU patients with SpO2; ",
        sum(!nc_cohort$hypoxemic & !nc_cohort$imv_set_vt), " non-hypoxemic non-ventilated, ",
        sum(!nc_cohort$hypoxemic & nc_cohort$imv_set_vt &
              !coalesce(nc_cohort$hypoxemic_during_vent, TRUE)),
        " ventilated and non-hypoxemic throughout ventilation (",
        sum(nc_cohort$imv_set_vt & is.na(nc_cohort$hypoxemic_during_vent)), " unclassifiable: no FiO2 during IMV)")

# =============================================================================
# Save intermediates
# =============================================================================

output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_parquet(nc_cohort, file.path(output_dir, "nc_cohort.parquet"))

saveRDS(eligible_hospitalizations, file.path(output_dir, "cohort_hospitalization_ids.rds"))
write_parquet(resp_waterfall, file.path(output_dir, "resp_support_waterfall.parquet"))
write_parquet(cohort_demographics, file.path(output_dir, "cohort_demographics.parquet"))
write_parquet(cohort_vitals, file.path(output_dir, "cohort_vitals.parquet"))
write_parquet(cohort_labs, file.path(output_dir, "cohort_labs.parquet"))
write_parquet(cohort_meds, file.path(output_dir, "cohort_meds.parquet"))
write_parquet(cohort_assessments, file.path(output_dir, "cohort_assessments.parquet"))
write_parquet(cohort_crrt, file.path(output_dir, "cohort_crrt.parquet"))
saveRDS(crrt_available, file.path(output_dir, "crrt_available.rds"))
write_parquet(cohort_heights, file.path(output_dir, "cohort_heights.parquet"))
write_parquet(cohort_weights, file.path(output_dir, "cohort_weights.parquet"))

# Partial attrition log (steps 1-3). Script 03 reads this back, appends the
# analytic-filter steps 4-7, and writes the complete log to final/.
write_csv(attrition, file.path(output_dir, "attrition_log_partial.csv"))

message("All intermediates saved to: ", output_dir)
message("Script 01 complete.")

