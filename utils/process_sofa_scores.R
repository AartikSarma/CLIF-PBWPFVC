#' Compute SOFA scores from a wide CLIF dataset
#'
#' R port of the Python SOFA scoring pipeline. Imputes PaO2 from SpO2 via
#' the Severinghaus equation, aggregates worst (extremal) values per encounter,
#' and computes all 6 SOFA component scores plus the total.

library(tidyverse)

# ------------------------------------------------------------------ #
# Constants                                                          #
# ------------------------------------------------------------------ #

REQUIRED_SOFA_CATEGORIES_BY_TABLE <- list(
  labs = c("creatinine", "platelet_count", "po2_arterial", "bilirubin_total"),
  vitals = c("map", "spo2"),
  patient_assessments = c("gcs_total"),
  medication_admin_continuous = c(
    "norepinephrine", "epinephrine", "dopamine", "dobutamine"
  ),
  respiratory_support = c("device_category", "fio2_set")
)

MAX_ITEMS <- c(
  "norepinephrine_mcg_kg_min", "epinephrine_mcg_kg_min",
  "dopamine_mcg_kg_min", "dobutamine_mcg_kg_min",
  "fio2_set", "creatinine", "bilirubin_total"
)

MIN_ITEMS <- c(
  "map", "spo2", "po2_arterial", "pao2_imputed",
  "platelet_count", "gcs_total"
)

DEVICE_RANK_MAPPING <- tibble(
  device_category = c(
    "IMV", "NIPPV", "CPAP", "High Flow NC", "Face Mask",
    "Trach Collar", "Nasal Cannula", "Other", "Room Air"
  ),
  device_rank = 1:9
)

# ------------------------------------------------------------------ #
# Internal helpers                                                   #
# ------------------------------------------------------------------ #

#' Impute PaO2 from SpO2 using Severinghaus equation.
#' Only applied when SpO2 < 97%.
.impute_pao2_from_spo2 <- function(wide_df) {
  wide_df |>
    mutate(
      .s = spo2 / 100,
      .a = 11700.0 / ((1 / .s) - 1),
      .b = sqrt(50^3 + .a^2),
      pao2_imputed = if_else(
        spo2 < 97,
        (.b + .a)^(1 / 3) - (.b - .a)^(1 / 3),
        NA_real_
      )
    ) |>
    select(-".s", -".a", -".b")
}

#' Aggregate worst (extremal) values per ID for SOFA scoring.
#' MAX for worse-when-higher variables, MIN for worse-when-lower.
.agg_extremal_values_by_id <- function(wide_df, extremal_type, id_name) {
  if (extremal_type == "latest") {
    stop("'latest' is a future feature and currently unavailable")
  }
  if (extremal_type != "worst") {
    stop(sprintf("Invalid extremal type: %s", extremal_type))
  }

  # Join device ranks
  df <- wide_df |>
    left_join(DEVICE_RANK_MAPPING, by = "device_category")

  # Identify which MAX/MIN columns actually exist
  max_cols <- intersect(MAX_ITEMS, names(df))
  min_cols <- intersect(MIN_ITEMS, names(df))

  df |>
    group_by(.data[[id_name]]) |>
    summarise(
      across(all_of(max_cols), \(x) max(x, na.rm = TRUE)),
      across(all_of(min_cols), \(x) min(x, na.rm = TRUE)),
      device_rank = min(device_rank, na.rm = TRUE),
      .groups = "drop"
    ) |>
    # Replace Inf/-Inf from all-NA groups back to NA
    mutate(across(where(is.numeric), \(x) if_else(is.infinite(x), NA_real_, x)))
}

#' Compute the 6 SOFA component scores + total from extremal values.
.compute_sofa_from_extremal_values <- function(extremal_df, id_name) {
  extremal_df |>
    left_join(DEVICE_RANK_MAPPING, by = "device_rank") |>
    mutate(
      p_f = po2_arterial / fio2_set,
      p_f_imputed = pao2_imputed / fio2_set,

      # -- Cardiovascular --
      sofa_cv_97 = case_when(
        dopamine_mcg_kg_min > 15 | epinephrine_mcg_kg_min > 0.1 | norepinephrine_mcg_kg_min > 0.1 ~ 4L,
        dopamine_mcg_kg_min > 5  | epinephrine_mcg_kg_min <= 0.1 | norepinephrine_mcg_kg_min <= 0.1 ~ 3L,
        dopamine_mcg_kg_min <= 5 | dobutamine_mcg_kg_min > 0 ~ 2L,
        map < 70  ~ 1L,
        map >= 70 ~ 0L
      ),

      # -- Coagulation --
      sofa_coag = case_when(
        platelet_count < 20  ~ 4L,
        platelet_count < 50  ~ 3L,
        platelet_count < 100 ~ 2L,
        platelet_count < 150 ~ 1L,
        platelet_count >= 150 ~ 0L
      ),

      # -- Liver --
      sofa_liver = case_when(
        bilirubin_total >= 12 ~ 4L,
        bilirubin_total >= 6  ~ 3L,
        bilirubin_total >= 2  ~ 2L,
        bilirubin_total >= 1.2 ~ 1L,
        bilirubin_total < 1.2 ~ 0L
      ),

      # -- Respiratory --
      sofa_resp = case_when(
        p_f < 100 & device_category %in% c("IMV", "NIPPV", "CPAP") ~ 4L,
        p_f >= 100 & p_f < 200 & device_category %in% c("IMV", "NIPPV", "CPAP") ~ 3L,
        p_f >= 200 & p_f < 300 ~ 2L,
        p_f >= 300 & p_f < 400 ~ 1L,
        p_f >= 400 ~ 0L
      ),

      # -- CNS --
      sofa_cns = case_when(
        gcs_total < 6   ~ 4L,
        gcs_total <= 9  ~ 3L,
        gcs_total <= 12 ~ 2L,
        gcs_total <= 14 ~ 1L,
        gcs_total == 15 ~ 0L
      ),

      # -- Renal --
      sofa_renal = case_when(
        creatinine >= 5   ~ 4L,
        creatinine >= 3.5 ~ 3L,
        creatinine >= 2   ~ 2L,
        creatinine >= 1.2 ~ 1L,
        creatinine < 1.2  ~ 0L
      ),

      sofa_total = sofa_cv_97 + sofa_coag + sofa_liver + sofa_resp + sofa_renal + sofa_cns
    )
}

#' Fill missing SOFA sub-scores with 0 (absence of data = no organ failure).
.fill_na_scores <- function(sofa_df) {
  subscore_cols <- c("sofa_cv_97", "sofa_coag", "sofa_renal",
                     "sofa_liver", "sofa_resp", "sofa_cns")

  # Recalculate total ignoring NAs, then zero-fill subscores
  sofa_df |>
    mutate(
      sofa_total = rowSums(pick(all_of(subscore_cols)), na.rm = TRUE),
      across(all_of(subscore_cols), \(x) replace_na(x, 0L))
    )
}

# ------------------------------------------------------------------ #
# Main entry point                                                   #
# ------------------------------------------------------------------ #

#' Compute SOFA scores from a wide dataset.
#'
#' @param wide_df Wide dataset containing all required SOFA variables.
#'   Medication columns should be pre-converted to standard units
#'   (e.g., \code{norepinephrine_mcg_kg_min}).
#' @param cohort_df Optional tibble with columns
#'   \code{[id_name, start_time, end_time]} to filter observations by
#'   time window.
#' @param extremal_type \code{"worst"} (default) or \code{"latest"} (future).
#' @param id_name Grouping column
#'   (\code{"hospitalization_id"}, \code{"encounter_block"}, etc.).
#' @param fill_na_scores_with_zero Logical. If TRUE fill missing component
#'   scores with 0.
#' @param remove_outliers Logical. If TRUE clamp po2, fio2, spo2 to
#'   physiological ranges.
#'
#' @return A tibble with SOFA component scores and total per ID.
compute_sofa <- function(wide_df,
                         cohort_df = NULL,
                         extremal_type = "worst",
                         id_name = "encounter_block",
                         fill_na_scores_with_zero = TRUE,
                         remove_outliers = TRUE) {

  # -- Validate --
  stopifnot(extremal_type %in% c("worst", "latest"))
  if (!id_name %in% names(wide_df)) {
    stop(sprintf("id_name '%s' not found in wide_df columns", id_name))
  }

  # -- Cohort time filtering --
  if (!is.null(cohort_df)) {
    required_cols <- c(id_name, "start_time", "end_time")
    missing_cols <- setdiff(required_cols, names(cohort_df))
    if (length(missing_cols) > 0) {
      stop(sprintf(
        "cohort_df must contain columns: %s. Missing: %s",
        paste(required_cols, collapse = ", "),
        paste(missing_cols, collapse = ", ")
      ))
    }
    wide_df <- wide_df |>
      inner_join(cohort_df, by = id_name) |>
      filter(event_time >= start_time, event_time <= end_time) |>
      select(-start_time, -end_time)
  }

  # -- Outlier removal --
  if (remove_outliers) {
    message("Removing outliers from wide dataset")
    wide_df <- wide_df |>
      mutate(
        po2_arterial = if_else(
          between(po2_arterial, 0, 700), po2_arterial, NA_real_
        ),
        fio2_set = if_else(
          between(fio2_set, 0.21, 1), fio2_set, NA_real_
        ),
        spo2 = if_else(
          between(spo2, 50, 100), spo2, NA_real_
        )
      )
  }

  # -- Pipeline --
  sofa_scores <- wide_df |>
    .impute_pao2_from_spo2() |>
    .agg_extremal_values_by_id(extremal_type, id_name) |>
    .compute_sofa_from_extremal_values(id_name)

  if (fill_na_scores_with_zero) {
    sofa_scores <- .fill_na_scores(sofa_scores)
  }

  sofa_scores
}
