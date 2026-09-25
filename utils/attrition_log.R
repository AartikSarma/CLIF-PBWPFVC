# =============================================================================
# Attrition logging for the PBW-vs-PFVC cohort funnel
# =============================================================================
# A small, dependency-light accumulator so scripts 01 and 03 can record the
# cohort size at each inclusion step and export a structured CONSORT table.
# Because the seven steps are split across two R sessions (01 builds the cohort,
# 03 applies the analytic filters), script 01 writes a partial log and script 03
# reads it back and appends the remaining steps.
#
# The step ordering and labels are hard-coded constants, so every site emits the
# identical, ordered set of steps and the logs can be compared or summed across
# sites step by step. utils/site_anonymization.R ranks sites by the last step's n.

library(tidyverse)

# Canonical inclusion steps of the ventilated (imv) cohort, in order. Steps 1-3 are
# logged in script 01, steps 4-7 in script 03. These labels are pooled across sites:
# do not reword them.
ATTRITION_STEPS <- c(
  "Adults (age >= 18)",
  "ICU admission",
  "Invasive ventilation with set tidal volume",
  "Height 150-210 cm",
  "Complete index data (VT/PBW, VT/PFVC, SF ratio)",
  "Lung-protective VT/PBW 6-8 mL/kg",
  "Hypoxemic (SF ratio < 315)"
)

# The nosupport and niv cohorts follow the same seven steps, but steps 3, 5, 6 and 7 mean
# something else there: entry is by first respiratory support, the index needs only an
# SF ratio, no tidal-volume band applies (step 6 excludes no one), and the no-support
# control's last step is its ICU-admission index and 24-hour escalation landmark.
ATTRITION_STEPS_CONTROL <- list(
  nosupport = c(
    ATTRITION_STEPS[1:2],
    "Room air or nasal cannula before any advanced support",
    ATTRITION_STEPS[4],
    "Complete index data (SF ratio)",
    "No tidal-volume band (no set tidal volume)",
    "Indexed at ICU admission, not escalated within 24 h"
  ),
  niv = c(
    ATTRITION_STEPS[1:2],
    "High-flow nasal cannula or non-invasive ventilation as first advanced support",
    ATTRITION_STEPS[4],
    "Complete index data (SF ratio)",
    "No tidal-volume band (no set tidal volume)",
    ATTRITION_STEPS[7]
  )
)
attrition_steps_for <- function(cohort) {
  if (cohort == "imv") ATTRITION_STEPS else ATTRITION_STEPS_CONTROL[[cohort]]
}

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
