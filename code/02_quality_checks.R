# =============================================================================
# Script 02: Quality Checks
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(data.table)
library(collapse)

source("utils/config.R")

site_name <- config$site_name
output_dir <- here("output", paste0(site_name, "_output"), "intermediate")

# Utility functions from outlier_handler.R (sourced inline to avoid side effects)
replace_outliers_with_na_long <- function(df, df_outlier_thresholds,
                                          category_variable, numeric_variable) {
  df %>%
    left_join(df_outlier_thresholds, by = category_variable) %>%
    mutate(!!sym(numeric_variable) := ifelse(
      !is.na(lower_limit) & !is.na(upper_limit) &
        (get(numeric_variable) < lower_limit | get(numeric_variable) > upper_limit),
      NA,
      get(numeric_variable)
    )) %>%
    select(-lower_limit, -upper_limit)
}

generate_summary_stats <- function(data, category_variable, numeric_variable) {
  data %>%
    group_by({{ category_variable }}) %>%
    summarise(
      N = sum(!is.na({{ numeric_variable }})),
      Min = min({{ numeric_variable }}, na.rm = TRUE),
      Max = max({{ numeric_variable }}, na.rm = TRUE),
      Mean = mean({{ numeric_variable }}, na.rm = TRUE),
      Median = median({{ numeric_variable }}, na.rm = TRUE),
      Q1 = quantile({{ numeric_variable }}, 0.25, na.rm = TRUE),
      Q3 = quantile({{ numeric_variable }}, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange({{ category_variable }})
}

# =============================================================================
# Load intermediate data from Script 01
# =============================================================================

cohort_ids <- readRDS(file.path(output_dir, "cohort_hospitalization_ids.rds"))
resp_waterfall <- read_parquet(file.path(output_dir, "resp_support_waterfall.parquet"))
cohort_vitals <- read_parquet(file.path(output_dir, "cohort_vitals.parquet"))
cohort_labs <- read_parquet(file.path(output_dir, "cohort_labs.parquet"))
cohort_heights <- read_parquet(file.path(output_dir, "cohort_heights.parquet"))

message("Loaded intermediate data")

# =============================================================================
# Labs: apply outlier thresholds
# =============================================================================

labs_thresholds <- read_csv(here("outlier-thresholds", "outlier_thresholds_labs.csv"),
                            show_col_types = FALSE) %>%
  mutate(across(c(lower_limit, upper_limit), as.numeric))

cohort_labs_clean <- replace_outliers_with_na_long(
  cohort_labs, labs_thresholds, "lab_category", "lab_value_numeric"
)

lab_summary <- generate_summary_stats(cohort_labs_clean, lab_category, lab_value_numeric)
message("Lab summary stats:")
print(lab_summary)

# =============================================================================
# Vitals: apply outlier thresholds
# =============================================================================

vitals_thresholds <- read_csv(
  here("outlier-thresholds", "outlier_thresholds_adults_vitals.csv"),
  show_col_types = FALSE
) %>%
  mutate(across(c(lower_limit, upper_limit), as.numeric))

cohort_vitals_clean <- replace_outliers_with_na_long(
  cohort_vitals, vitals_thresholds, "vital_category", "vital_value"
)

vital_summary <- generate_summary_stats(cohort_vitals_clean, vital_category, vital_value)
message("Vital summary stats:")
print(vital_summary)

# =============================================================================
# Heights: apply manual thresholds (120-230 cm)
# =============================================================================

cohort_heights_clean <- cohort_heights %>%
  mutate(height_cm = if_else(
    height_cm < 120 | height_cm > 230, NA_real_, height_cm
  ))

message("Heights before cleaning: ", sum(!is.na(cohort_heights$height_cm)),
        ", after: ", sum(!is.na(cohort_heights_clean$height_cm)))

# =============================================================================
# Respiratory support: apply outlier thresholds to waterfall output
# =============================================================================

resp_thresholds <- read_csv(
  here("outlier-thresholds", "outlier_thresholds_respiratory_support.csv"),
  show_col_types = FALSE
) %>%
  mutate(across(c(lower_limit, upper_limit), as.numeric))

# Apply thresholds to each numeric column in the waterfall output
resp_waterfall_clean <- resp_waterfall

for (i in seq_len(nrow(resp_thresholds))) {
  var_name <- resp_thresholds$variable_name[i]
  lower <- resp_thresholds$lower_limit[i]
  upper <- resp_thresholds$upper_limit[i]

  if (var_name %in% names(resp_waterfall_clean)) {
    n_before <- sum(!is.na(resp_waterfall_clean[[var_name]]))
    resp_waterfall_clean[[var_name]] <- if_else(
      resp_waterfall_clean[[var_name]] < lower | resp_waterfall_clean[[var_name]] > upper,
      NA_real_,
      resp_waterfall_clean[[var_name]]
    )
    n_after <- sum(!is.na(resp_waterfall_clean[[var_name]]))
    message(var_name, ": ", n_before - n_after, " outliers removed")
  }
}

# =============================================================================
# Save cleaned data
# =============================================================================

write_parquet(cohort_labs_clean, file.path(output_dir, "cohort_labs_clean.parquet"))
write_parquet(cohort_vitals_clean, file.path(output_dir, "cohort_vitals_clean.parquet"))
write_parquet(cohort_heights_clean, file.path(output_dir, "cohort_heights_clean.parquet"))
write_parquet(resp_waterfall_clean, file.path(output_dir, "resp_support_waterfall_clean.parquet"))

# Save summary stats
dir.create(file.path(output_dir, "summary_stats"), recursive = TRUE, showWarnings = FALSE)
write_csv(lab_summary, file.path(output_dir, "summary_stats",
                                  paste0("lab_summary_", site_name, ".csv")))
write_csv(vital_summary, file.path(output_dir, "summary_stats",
                                    paste0("vital_summary_", site_name, ".csv")))

message("Script 02 complete. Cleaned data saved.")
