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
  # Cohort: "imv" (the analytic cohort: invasive ventilation with a set tidal
  # volume) or "niv" (the never-intubated negative control for the biotrauma
  # suite: first advanced support is high-flow nasal cannula or non-invasive
  # ventilation, no invasive ventilation before it; intubation later is a
  # competing event). PBWPFVC_COHORT; run the control under its own site name
  # (PBWPFVC_SITE_NAME={site}_niv) so its outputs sit beside, not over, the site's.
  config$cohort <- Sys.getenv("PBWPFVC_COHORT", "imv")
  if (!config$cohort %in% c("imv", "niv"))
    stop("PBWPFVC_COHORT must be imv or niv; got '", config$cohort, "'")
  if (config$cohort == "niv") message("  cohort: niv (never-intubated control; PBWPFVC_COHORT)")
  return(config)
}
# device categories of the never-intubated control (CLIF mCIDE, lower case)
NIV_DEVICES <- c("high flow nc", "nippv", "cpap")

# Load the configuration
config <- load_config()
