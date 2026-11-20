# Packages come from the project library that uvr installs (uvr.lock).
library(jsonlite)

# The data zones of SOP.md (the CLIF cloud-environment template, project-template/ at
# commit 9813e85). utils/clif_paths.R is that template's helper file, vendored
# unchanged; this file builds the project's configuration on top of it.
source(file.path("utils", "clif_paths.R"))   # load_config(), phi_columns(), write_manifest(), ...

# Load the site configuration from config/config.json (keys: config_template.json),
# then apply any overrides passed through the environment. The pipeline runner
# (code/00_run_pipeline.R) translates its --site_name / --site_path / --file_type
# arguments into these variables so that every numbered script, each a separate
# Rscript subprocess, sees the same override without config.json being edited:
#   PBWPFVC_SITE_NAME    -> site_name
#   PBWPFVC_TABLES_PATH  -> data_directory
#   PBWPFVC_FILE_TYPE    -> filetype
# The same variables work when a single script is run by hand. config.json must
# still exist (it carries any field not overridden).
#
# Where files go (SOP.md, the four standards):
#   data_directory           the core CLIF tables, read in place and never written
#   intermediate_directory   row-level working files, outside the project folder
#                            (default ../clif_workspace/<repo>/intermediate)
#   results_directory        aggregates from a run (default ../clif_workspace/<repo>/results)
#   output/final_no_phi/     reviewed aggregates only, copied there by a person with
#                            release_result(); no script writes there
# The scripts keep their internal names: config$tables_path is data_directory,
# config$file_type is filetype, config$output_dir the cohort's intermediate folder,
# config$final_root the site's results folder.
load_project_config <- function() {
  json_path <- file.path("config", "config.json")
  if (!file.exists(json_path))
    stop("Configuration file not found: copy config/config_template.json to config/config.json")
  raw <- fromJSON(json_path)
  old <- intersect(c("tables_path", "file_type"), names(raw))
  if (length(old) && is.null(raw$data_directory))
    stop("config/config.json uses the old keys (", paste(old, collapse = ", "), "). Rename tables_path to ",
         "data_directory and file_type to filetype, and add the other keys of config/config_template.json.")
  overrides <- c(site_name = "PBWPFVC_SITE_NAME", data_directory = "PBWPFVC_TABLES_PATH",
                 filetype = "PBWPFVC_FILE_TYPE")
  for (field in names(overrides)) {
    value <- Sys.getenv(overrides[[field]], unset = "")
    if (nzchar(value)) {
      raw[[field]] <- value
      message("  ", field, " overridden from ", overrides[[field]], ": ", value)
    }
  }
  # the template's loader checks the paths (outside the project, outside the data) and
  # creates the workspace folders; it reads a file, so the overridden copy is written to one
  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp), add = TRUE)
  write_json(raw, tmp, auto_unbox = TRUE)
  cfg <- load_config(tmp)
  message("Loaded configuration from config.json")
  config <- cfg
  config$tables_path <- cfg$data_directory
  config$file_type <- cfg$filetype
  if (!config$file_type %in% c("parquet", "csv", "fst"))
    stop("config$file_type must be parquet, csv or fst; got '", config$file_type, "'")
  # Cohort (PBWPFVC_COHORT; environment variable only):
  #   "imv"        default, the paper's cohort: invasive ventilation with a set tidal volume
  #   "nosupport"  the negative control: room air or nasal cannula at ICU admission,
  #                no advanced support before the index, so strain per lung size
  #                cannot act; escalation to any advanced support later ends
  #                follow-up as a competing event or a censoring point
  #   "niv"        built on request only, not part of the paper: first advanced
  #                support is high-flow nasal cannula or non-invasive ventilation,
  #                with no invasive ventilation before it; intubation later is a
  #                competing event
  config$cohort <- Sys.getenv("PBWPFVC_COHORT", "imv")
  if (!config$cohort %in% c("imv", "niv", "nosupport"))
    stop("PBWPFVC_COHORT must be imv, niv or nosupport; got '", config$cohort, "'")
  # Where a cohort's files live (both folders outside the project; see above):
  #   <intermediate_directory>/{site}/                    patient-level, never shared
  #   <intermediate_directory>/{site}/controls/{cohort}/  the same for a control cohort
  #   <results_directory>/{site}/                         aggregates only, sorted by
  #                                      manuscript block (final_dir_for() below):
  #     cross_sectional/                 figures 1-3 (scripts 03-05)
  #     injury/                          figure 4 (scripts 21-28)
  #     supplement/                      supplementary analyses (code/supplement/)
  #     controls/                        the control cohorts' aggregates, all in one folder
  # A control's FILE NAMES carry {site}_{cohort} (e.g. jm_estimates_pfvc_7d_MIMIC_nosupport.csv),
  # so config$site_name is that tag and config$base_site is the site itself. Setting
  # PBWPFVC_COHORT is enough; a PBWPFVC_SITE_NAME that already ends in _{cohort} is
  # read the same way (the suffix is not doubled). The pooling scripts list final/ without
  # recursing, so final/controls/ never leaks into a cross-site pool of the main cohort.
  config$base_site <- config$site_name
  if (config$cohort != "imv") {
    config$base_site <- sub(paste0("_", config$cohort, "$"), "", config$site_name)
    config$site_name <- paste0(config$base_site, "_", config$cohort)
    message("  cohort: ", config$cohort, " (PBWPFVC_COHORT); files tagged ", config$site_name)
  }
  # The VT/PBW gate (PBWPFVC_VTPBW_GATE; environment variable only):
  #   "1"  default, the paper's cohort: a ventilated patient enters only with a
  #        complete hypoxemic timepoint at VT/PBW 6-8 mL/kg, and is indexed there
  #   "0"  the ungated ventilated cohort, for the gate's sensitivity analyses: the same
  #        rule without the VT/PBW band (code/supplement/xsec_vtpbw_gate_collider.R
  #        shows why the band matters). Script 03 then writes
  #        analysis_cross_sectional_ungated.parquet and nothing downstream of it; the
  #        figure-4 scripts build and fit their panel from that table, and every
  #        aggregate goes to final/ungated/<block>/, apart from the paper's tables.
  config$vtpbw_gate <- Sys.getenv("PBWPFVC_VTPBW_GATE", "1")
  if (!config$vtpbw_gate %in% c("0", "1"))
    stop("PBWPFVC_VTPBW_GATE must be 0 or 1; got '", config$vtpbw_gate, "'")
  if (config$vtpbw_gate == "0" && config$cohort != "imv")
    stop("PBWPFVC_VTPBW_GATE=0 applies to the ventilated cohort only; unset PBWPFVC_COHORT")
  config$cs_suffix <- if (config$vtpbw_gate == "0") "_ungated" else ""
  # one folder per site in each zone: <intermediate_directory>/<site>/ (a control cohort
  # in controls/<cohort>/ below it) and <results_directory>/<site>/
  config$intermediate_site <- file.path(config$intermediate_directory, config$base_site)
  config$output_dir <- if (config$cohort == "imv") config$intermediate_site else
    file.path(config$intermediate_site, "controls", config$cohort)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  config$final_root <- file.path(config$results_directory, config$base_site)
  config$final_dir  <- if (config$cohort == "imv") config$final_root else file.path(config$final_root, "controls")
  return(config)
}
# The folder a script writes its aggregates to: final/<block>/ for the ventilated cohort
# (final/ungated/<block>/ with PBWPFVC_VTPBW_GATE=0).
# A control cohort keeps everything in final/controls/, whatever the block, because its
# file names already say which cohort they are and the comparison reads one folder.
FINAL_BLOCKS <- c("cross_sectional", "injury", "supplement")
final_dir_for <- function(block) {
  if (!block %in% FINAL_BLOCKS) stop("final_dir_for(): unknown block '", block, "'; blocks are ", paste(FINAL_BLOCKS, collapse = ", "))
  block_dir <- if (config$cohort != "imv") config$final_dir else
    if (config$vtpbw_gate == "0") file.path(config$final_root, "ungated", block) else file.path(config$final_root, block)
  dir.create(block_dir, recursive = TRUE, showWarnings = FALSE)
  block_dir
}
# device categories (CLIF mCIDE, lower case): the niv cohort's, the control's, and
# everything that counts as advanced support (escalation)
# The age range of every cohort: adults, and under 90. PFVC, and so the exposure, is a
# function of age, and an age above 89 is withheld under HIPAA Safe Harbor (SOP.md
# top-codes it at 90), so a patient whose PFVC cannot be computed from their true age
# is not entered.
AGE_MIN <- 18   # inclusive
AGE_MAX <- 90   # exclusive: age_at_admission < 90
in_age_range <- function(age) !is.na(age) & age >= AGE_MIN & age < AGE_MAX

NIV_DEVICES       <- c("high flow nc", "nippv", "cpap")
NOSUPPORT_DEVICES <- c("room air", "nasal cannula")
SUPPORT_DEVICES   <- c("imv", NIV_DEVICES)
# An SpO2 or PaO2 is paired with the most recent FiO2 recorded up to this many hours
# before it (scripts 01 and 03). The window is short so that the ratio reflects the
# support in force early in ventilation, not a setting charted hours earlier.
FIO2_LOOKBACK_H <- 4
# A death within this many hours of the last IMV record is a death on the ventilator,
# not a death after liberation (03's ventilator-free days, 10's event clock).
DEATH_ON_VENT_TOL_H <- 1
# FiO2 on room air and nasal cannula, for the no-support control only (the
# analytic cohort's SF uses documented FiO2): 0.21 on room air, 0.21 + 0.03 per
# L/min on a cannula capped at 0.60, the rule script 01 uses for its negative-control
# cohorts. Documented fio2_set is kept where present.
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

# Small cells: no script masks or suppresses a count. Every table in final/
# carries raw counts; masking, whether censoring small cells or a deterministic scheme
# agreed with CLIF, is a separate step applied to final/ before a site shares it.

# Load the configuration
config <- load_project_config()
source(file.path("utils", "intermediate_io.R"))   # read_parquet()/write_parquet() for the intermediate zone
