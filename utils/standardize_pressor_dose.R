# =============================================================================
# standardize_pressor_dose
#
# CLIF medication_admin_continuous records vasopressor doses in one of several
# units. Downstream SOFA and NE-equivalent calculations assume mcg/kg/min, so
# this helper normalizes each row against the patient's weight.
#
# Supported units (case-insensitive, whitespace-tolerant):
#   - mcg/kg/min or ug/kg/min       -> as-is
#   - mcg/min or ug/min             -> divide by weight_kg
#   - mg/kg/min                     -> multiply by 1000
#   - mg/min                        -> multiply by 1000 / weight_kg
#   - mcg/kg/hr or ug/kg/hr         -> divide by 60
#   - mcg/hr or ug/hr               -> divide by (60 * weight_kg)
#
# Rows with unrecognized units or missing weight (when weight is required)
# return NA for the standardized dose. The helper prints a summary of any
# rows it could not convert so callers can investigate upstream data issues
# rather than silently dropping them.
# =============================================================================

library(dplyr)
library(stringr)

standardize_pressor_dose <- function(meds, weights,
                                     dose_col = "med_dose",
                                     unit_col = "med_dose_unit",
                                     out_col = "dose_mcg_kg_min",
                                     id_col = "hospitalization_id",
                                     label = "pressor") {
  stopifnot(is.data.frame(meds), is.data.frame(weights))
  stopifnot(all(c(id_col, dose_col, unit_col) %in% names(meds)))
  stopifnot(all(c(id_col, "weight_kg") %in% names(weights)))

  weights_distinct <- weights %>%
    distinct(.data[[id_col]], .keep_all = TRUE) %>%
    select(all_of(id_col), weight_kg)

  meds <- meds %>%
    mutate(
      .dose_numeric = suppressWarnings(as.numeric(.data[[dose_col]])),
      .unit_clean = str_to_lower(str_replace_all(.data[[unit_col]], "\\s+", ""))
    ) %>%
    left_join(weights_distinct, by = id_col)

  meds <- meds %>%
    mutate(
      !!out_col := case_when(
        is.na(.dose_numeric) ~ NA_real_,
        .unit_clean %in% c("mcg/kg/min", "ug/kg/min", "mcg/kgs/min") ~ .dose_numeric,
        .unit_clean %in% c("mg/kg/min") ~ .dose_numeric * 1000,
        .unit_clean %in% c("mcg/kg/hr", "ug/kg/hr", "mcg/kg/h", "ug/kg/h") ~
          .dose_numeric / 60,
        .unit_clean %in% c("mcg/min", "ug/min") & !is.na(weight_kg) & weight_kg > 0 ~
          .dose_numeric / weight_kg,
        .unit_clean %in% c("mg/min") & !is.na(weight_kg) & weight_kg > 0 ~
          .dose_numeric * 1000 / weight_kg,
        .unit_clean %in% c("mcg/hr", "ug/hr", "mcg/h", "ug/h") &
          !is.na(weight_kg) & weight_kg > 0 ~
          .dose_numeric / (weight_kg * 60),
        TRUE ~ NA_real_
      )
    )

  n_total <- nrow(meds)
  n_success <- sum(!is.na(meds[[out_col]]))
  n_drop_weight <- sum(
    !is.na(meds$.dose_numeric) &
      meds$.unit_clean %in% c("mcg/min", "ug/min", "mg/min",
                              "mcg/hr", "ug/hr", "mcg/h", "ug/h") &
      (is.na(meds$weight_kg) | meds$weight_kg <= 0)
  )
  unrecognized_units <- meds %>%
    filter(!is.na(.dose_numeric), is.na(.data[[out_col]])) %>%
    pull(.unit_clean) %>%
    unique() %>%
    setdiff(c(NA_character_, "", "mcg/min", "ug/min", "mg/min",
              "mcg/hr", "ug/hr", "mcg/h", "ug/h"))

  message(sprintf(
    "[%s] standardized %d/%d rows to mcg/kg/min (weight missing: %d)",
    label, n_success, n_total, n_drop_weight
  ))
  if (length(unrecognized_units) > 0) {
    message(sprintf("[%s] unrecognized units: %s", label,
                    paste(unrecognized_units, collapse = ", ")))
  }

  meds %>% select(-.dose_numeric, -.unit_clean, -weight_kg)
}
