# Load necessary libraries
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}

library(jsonlite)

# Load the site configuration from config/config.json, then apply any overrides
# passed through the environment. The pipeline runner (code/00_run_pipeline.R)
# translates its --site_name / --site_path / --file_type arguments into these
# variables so that every numbered script, each a separate Rscript subprocess,
# sees the same override without config.json being edited:
#   PBWPFVC_SITE_NAME    -> config$site_name
#   PBWPFVC_TABLES_PATH  -> config$tables_path
#   PBWPFVC_FILE_TYPE    -> config$file_type
# The same variables work when a single script is run by hand. config.json must
# still exist (it carries any field not overridden).
load_config <- function() {
  json_path <- "config/config.json"
  if (file.exists(json_path)) {
    config <- fromJSON(json_path)
    message("Loaded configuration from config.json")
  } else {
    stop("Configuration file not found. Please create config.json",
         "based on the config_template.")
  }
  overrides <- c(site_name = "PBWPFVC_SITE_NAME", tables_path = "PBWPFVC_TABLES_PATH",
                 file_type = "PBWPFVC_FILE_TYPE")
  for (field in names(overrides)) {
    value <- Sys.getenv(overrides[[field]], unset = "")
    if (nzchar(value)) {
      config[[field]] <- value
      message("  config$", field, " overridden from ", overrides[[field]], ": ", value)
    }
  }
  if (!config$file_type %in% c("parquet", "csv", "fst"))
    stop("config$file_type must be parquet, csv or fst; got '", config$file_type, "'")
  # Cohort (PBWPFVC_COHORT), the strain gradient of the biotrauma suite:
  #   "imv"        the analytic cohort: invasive ventilation with a set tidal volume
  #   "niv"        the middle arm: first advanced support is high-flow nasal cannula
  #                or non-invasive ventilation (spontaneous volumes, titrated
  #                pressures; a weaker, less controlled strain), no invasive
  #                ventilation before it; intubation later is a competing event
  #   "nosupport"  the negative control: room air or nasal cannula only, no
  #                advanced support before the index nor in the 24 h after it, so
  #                strain per lung size cannot act; escalation to any support
  #                later is a competing event
  config$cohort <- Sys.getenv("PBWPFVC_COHORT", "imv")
  if (!config$cohort %in% c("imv", "niv", "nosupport"))
    stop("PBWPFVC_COHORT must be imv, niv or nosupport; got '", config$cohort, "'")
  # Where a cohort's files live. One site has ONE output folder, output/{site}_output/:
  #   intermediate/                      patient-level, never shared
  #   intermediate/controls/{cohort}/    the same for a control cohort
  #   final/                             aggregates only: the folder a site returns,
  #                                      sorted by manuscript block (final_dir_for() below):
  #     final/cross_sectional/           figures 1-3 (scripts 03-05)
  #     final/injury/                    figure 4 (scripts 21-28)
  #     final/supplement/                supplementary analyses (code/supplement/)
  #     final/controls/                  the control cohorts' aggregates, all in one folder
  # A control's FILE NAMES carry {site}_{cohort} (e.g. jm_estimates_pfvc_7d_MIMIC_nosupport.csv),
  # so config$site_name is that tag and config$base_site is the site itself. Setting
  # PBWPFVC_COHORT is enough; a PBWPFVC_SITE_NAME that already ends in _{cohort} (the
  # older convention) is read the same way. The pooling scripts list final/ without
  # recursing, so final/controls/ never leaks into a cross-site pool of the main cohort.
  config$base_site <- config$site_name
  if (config$cohort != "imv") {
    config$base_site <- sub(paste0("_", config$cohort, "$"), "", config$site_name)
    config$site_name <- paste0(config$base_site, "_", config$cohort)
    message("  cohort: ", config$cohort, " (PBWPFVC_COHORT); files tagged ", config$site_name)
  }
  site_root <- file.path(getwd(), "output", paste0(config$base_site, "_output"))
  config$output_dir <- if (config$cohort == "imv") file.path(site_root, "intermediate") else
    file.path(site_root, "intermediate", "controls", config$cohort)
  config$final_root <- file.path(site_root, "final")
  config$final_dir  <- if (config$cohort == "imv") config$final_root else file.path(config$final_root, "controls")
  return(config)
}
# The folder a script writes its aggregates to: final/<block>/ for the ventilated cohort.
# A control cohort keeps everything in final/controls/, whatever the block, because its
# file names already say which cohort they are and the comparison reads one folder.
FINAL_BLOCKS <- c("cross_sectional", "injury", "supplement")
final_dir_for <- function(block) {
  if (!block %in% FINAL_BLOCKS) stop("final_dir_for(): unknown block '", block, "'; blocks are ", paste(FINAL_BLOCKS, collapse = ", "))
  block_dir <- if (config$cohort == "imv") file.path(config$final_root, block) else config$final_dir
  dir.create(block_dir, recursive = TRUE, showWarnings = FALSE)
  block_dir
}
# device categories (CLIF mCIDE, lower case): the middle arm's, the control's, and
# everything that counts as advanced support (escalation)
NIV_DEVICES       <- c("high flow nc", "nippv", "cpap")
NOSUPPORT_DEVICES <- c("room air", "nasal cannula")
SUPPORT_DEVICES   <- c("imv", NIV_DEVICES)
# FiO2 on room air and nasal cannula, for the no-support control only (the
# analytic cohort's SF uses documented FiO2): 0.21 on room air, 0.21 + 0.03 per
# L/min on a cannula capped at 0.60, the rule script 01 uses for its mortality
# controls. Documented fio2_set is kept where present.
estimate_fio2_nosupport <- function(df) {
  if (config$cohort != "nosupport") return(df)
  dev <- tolower(df$device_category)
  lpm <- if ("lpm_set" %in% names(df)) suppressWarnings(as.numeric(df$lpm_set)) else NA_real_
  est <- dplyr::case_when(!is.na(df$fio2_set) ~ as.numeric(df$fio2_set),
                          dev == "room air" ~ 0.21,
                          dev == "nasal cannula" & !is.na(lpm) ~ pmin(0.21 + 0.03 * lpm, 0.60),
                          TRUE ~ NA_real_)
  df$fio2_set <- est
  df
}

# Small cells. Every aggregate a site returns reports no group of 1-9 patients or
# events (project rule: minimum cell size 10). mask_small_counts() blanks such values
# in any count column it knows, and is applied wherever a block writes a table of counts.
SMALL_CELL_MIN <- 10L
COUNT_COLUMNS <- c("n_obs", "n_patients", "n_deaths", "n_extubations", "n_rrt", "n_competing",
                   "patient_days", "patients_any", "patients_with_baseline", "patients_day0_baseline",
                   "patients_ge2_obs", "deaths_ge2", "extubations_ge2", "plateau_subset_ge2",
                   "rrt_before_index", "rrt_within_horizon", "esrd_censored_day0",
                   "creatinine_days_removed_rrt", "nonpositive_set_missing", "lag_missing_rows",
                   "movement_last_n", "patients_pre_ge1", "patients_pre_ge2")
# A count column is one named in COUNT_COLUMNS or one of them with a suffix
# ("n_patients_ventilated" in the difference-in-differences table), so a table that
# widens counts by arm cannot slip past the mask.
mask_small_counts <- function(df) {
  is_count <- vapply(names(df), function(nm) any(nm == COUNT_COLUMNS | startsWith(nm, paste0(COUNT_COLUMNS, "_"))), logical(1))
  for (col in names(df)[is_count]) {
    v <- suppressWarnings(as.numeric(df[[col]]))
    hide <- !is.na(v) & v > 0 & v < SMALL_CELL_MIN
    if (any(hide)) df[[col]][hide] <- NA
  }
  df
}

# Load the configuration
config <- load_config()
