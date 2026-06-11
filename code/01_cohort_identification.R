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

# =============================================================================
# Load CLIF tables
# =============================================================================

table_names <- c("patient", "hospitalization", "adt", "respiratory_support",
                 "vitals", "labs", "medication_admin_continuous",
                 "patient_assessments")

clif_tables <- list()
# Check if files are available
for (tbl in table_names) {
  fpath <- file.path(tables_path, paste0("clif_", tbl, ".", file_type))
  if (!file.exists(fpath)) {
    stop("Missing table: ", fpath)
  }
}

#If all files are available, load the data: 
if(file_type == "parquet"){
  for (tbl in table_names) {
    fpath <- file.path(tables_path, paste0("clif_", tbl, ".", file_type))
    clif_tables[[tbl]] <- arrow::read_parquet(fpath)
    message("Loaded ", tbl, ": ", nrow(clif_tables[[tbl]]), " rows")
  }
}

if(file_type == "csv"){
  for (tbl in table_names) {
    fpath <- file.path(tables_path, paste0("clif_", tbl, ".", file_type))
    clif_tables[[tbl]] <- readr::read_csv(fpath)
    message("Loaded ", tbl, ": ", nrow(clif_tables[[tbl]]), " rows")
  }
}

if(file_type == "fst"){
  for (tbl in table_names) {
    fpath <- file.path(tables_path, paste0("clif_", tbl, ".", file_type))
    clif_tables[[tbl]] <- fst::read_fst(fpath)
    message("Loaded ", tbl, ": ", nrow(clif_tables[[tbl]]), " rows")
  }
}

#Create R objects for analysis: 
clif_patient <- clif_tables$patient
clif_hospitalization <- clif_tables$hospitalization
clif_adt <- clif_tables$adt
clif_respiratory_support <- clif_tables$respiratory_support
clif_vitals <- clif_tables$vitals
clif_labs <- clif_tables$labs
clif_meds <- clif_tables$medication_admin_continuous
clif_assessments <- clif_tables$patient_assessments
rm(clif_tables)

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
  rename(spo2_value = vital_value) %>% 
  filter(spo2_value <= 97)

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

lab_categories_needed <- c("po2_arterial", "pco2_arterial", "creatinine", "bilirubin_total", "platelet_count")

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
# Save intermediates
# =============================================================================

output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(eligible_hospitalizations, file.path(output_dir, "cohort_hospitalization_ids.rds"))
write_parquet(resp_waterfall, file.path(output_dir, "resp_support_waterfall.parquet"))
write_parquet(cohort_demographics, file.path(output_dir, "cohort_demographics.parquet"))
write_parquet(cohort_vitals, file.path(output_dir, "cohort_vitals.parquet"))
write_parquet(cohort_labs, file.path(output_dir, "cohort_labs.parquet"))
write_parquet(cohort_meds, file.path(output_dir, "cohort_meds.parquet"))
write_parquet(cohort_assessments, file.path(output_dir, "cohort_assessments.parquet"))
write_parquet(cohort_heights, file.path(output_dir, "cohort_heights.parquet"))
write_parquet(cohort_weights, file.path(output_dir, "cohort_weights.parquet"))

# Partial attrition log (steps 1-3). Script 03 reads this back, appends the
# analytic-filter steps 4-7, and writes the complete log to final/.
write_csv(attrition, file.path(output_dir, "attrition_log_partial.csv"))

message("All intermediates saved to: ", output_dir)
message("Script 01 complete.")

