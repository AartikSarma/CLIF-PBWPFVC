# =============================================================================
# Attrition logging for the PBW-vs-PFVC cohort funnel
# =============================================================================
# A small, dependency-light accumulator so scripts 01 and 03 can record the
# cohort size at each inclusion step and export a structured CONSORT table.
# Because the seven steps are split across two R sessions (01 builds the cohort,
# 03 applies the analytic filters), script 01 writes a partial log and script 03
# reads it back and appends the remaining steps.
#
# The step ordering and labels are hard-coded constants: pooling across sites
# (script 05) sums n by step_order, which is only exact if every site emits the
# identical, ordered set of steps.

library(tidyverse)

# Canonical inclusion steps, in order. Steps 1-3 are logged in script 01,
# steps 4-7 in script 03.
ATTRITION_STEPS <- c(
  "Adults (age >= 18)",
  "ICU admission",
  "Invasive ventilation with set tidal volume",
  "Height 150-210 cm",
  "Complete index data (VT/PBW, VT/PFVC, SF ratio, SOFA)",
  "Lung-protective VT/PBW 6-8 mL/kg",
  "Hypoxemic (SF ratio < 315)"
)

# Empty log with the correct column schema/types.
attrition_init <- function() {
  tibble(
    step_order       = integer(),
    step_label       = character(),
    n_remaining      = integer(),
    n_excluded       = integer(),
    exclusion_reason = character()
  )
}

# Append one step. step_order is the running row count; n_excluded is derived
# from the previous row's n_remaining (NA for the first step).
attrition_add <- function(log_tbl, step_label, n_remaining,
                          exclusion_reason = NA_character_) {
  prev_n <- if (nrow(log_tbl) == 0) NA_integer_ else dplyr::last(log_tbl$n_remaining)
  n_excluded <- if (is.na(prev_n)) NA_integer_ else as.integer(prev_n - n_remaining)
  bind_rows(
    log_tbl,
    tibble(
      step_order       = nrow(log_tbl) + 1L,
      step_label       = step_label,
      n_remaining      = as.integer(n_remaining),
      n_excluded       = n_excluded,
      exclusion_reason = exclusion_reason
    )
  )
}
