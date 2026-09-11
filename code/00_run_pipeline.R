# =============================================================================
# Script 00: Run the full pipeline
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Single entry point for running this project at a CLIF consortium site.
#
#   1. Restores the project environment from renv.lock.
#   2. Runs scripts 01-05 in order, each as a clean R subprocess.
#
# Cross-cohort pooling is a separate, centrally-run step (code/pooled_estimates.R)
# and is intentionally not invoked here.
#
# Usage (from the project root, or anywhere — the script locates the repo):
#   Rscript code/00_run_pipeline.R
#   Rscript code/00_run_pipeline.R --site_name my_site --site_path /data/clif [--file_type parquet]
#
# The optional arguments override the matching fields of config/config.json for
# this run only (config.json is not edited): --site_name sets config$site_name,
# which names the output folder output/<site_name>_output/; --site_path sets
# config$tables_path, the folder holding the clif_*.<file_type> tables;
# --file_type sets config$file_type (parquet, csv or fst). Both `--flag value`
# and `--flag=value` forms are accepted. The overrides reach every numbered
# script through environment variables read by utils/config.R
# (PBWPFVC_SITE_NAME, PBWPFVC_TABLES_PATH, PBWPFVC_FILE_TYPE), so a single script
# can be re-run by hand with the same override, e.g.
#   PBWPFVC_SITE_NAME=my_site Rscript code/04_analysis.R
#
# Each numbered script is standalone: it reads its inputs from disk and writes
# its outputs back to disk, so they are run as separate subprocesses rather than
# sourced. This isolates package namespaces and means script 01's
# `rm(list = ls())` cannot wipe state belonging to this runner. If any script
# exits with an error, the pipeline stops and reports which one failed.
# =============================================================================

# --- Locate the repository root ----------------------------------------------
# Find this file's path from the Rscript invocation so the pipeline can be
# launched from any working directory.
get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    return(normalizePath(sub("^--file=", "", file_arg)))
  }
  # Fallback for interactive source()
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
  stop("Unable to determine the path to 00_run_pipeline.R. ",
       "Run it with: Rscript code/00_run_pipeline.R")
}

script_path <- get_script_path()
repo_root   <- normalizePath(file.path(dirname(script_path), ".."))
setwd(repo_root)

# --- Command-line overrides of config.json -----------------------------------
parse_pipeline_args <- function(args) {
  known <- c(site_name = "PBWPFVC_SITE_NAME", site_path = "PBWPFVC_TABLES_PATH",
             file_type = "PBWPFVC_FILE_TYPE")
  usage <- paste0("Usage: Rscript code/00_run_pipeline.R [--site_name NAME] [--site_path DIR] ",
                  "[--file_type parquet|csv|fst]")
  out <- list(); i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (grepl("^--[a-z_]+=", a)) {                       # --flag=value
      key <- sub("^--([a-z_]+)=.*$", "\\1", a); val <- sub("^--[a-z_]+=", "", a); i <- i + 1L
    } else if (grepl("^--[a-z_]+$", a)) {                # --flag value
      key <- sub("^--", "", a)
      if (i == length(args) || grepl("^--", args[[i + 1L]]))
        stop("Missing value for --", key, "\n", usage)
      val <- args[[i + 1L]]; i <- i + 2L
    } else stop("Unrecognized argument: ", a, "\n", usage)
    if (!key %in% names(known)) stop("Unknown option --", key, "\n", usage)
    if (!nzchar(val)) stop("Empty value for --", key, "\n", usage)
    out[[known[[key]]]] <- val
  }
  out
}
cli_overrides <- parse_pipeline_args(commandArgs(trailingOnly = TRUE))
if (length(cli_overrides)) {
  if (!is.null(cli_overrides$PBWPFVC_TABLES_PATH)) {
    p <- path.expand(cli_overrides$PBWPFVC_TABLES_PATH)
    if (!dir.exists(p)) stop("--site_path does not exist: ", p)
    cli_overrides$PBWPFVC_TABLES_PATH <- normalizePath(p)
  }
  if (!is.null(cli_overrides$PBWPFVC_FILE_TYPE) &&
      !cli_overrides$PBWPFVC_FILE_TYPE %in% c("parquet", "csv", "fst"))
    stop("--file_type must be parquet, csv or fst")
  # Child Rscript processes inherit the environment, and utils/config.R applies these.
  do.call(Sys.setenv, cli_overrides)
}

message("=============================================================")
message("PBW vs PFVC pipeline runner")
message("Repository root: ", repo_root)
if (length(cli_overrides)) {
  for (k in names(cli_overrides)) message("Override ", k, " = ", cli_overrides[[k]])
} else message("Config: config/config.json (no command-line overrides)")
message("=============================================================")

# --- 1. Restore the project environment --------------------------------------
message("\n[00] Restoring renv environment from renv.lock ...")
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}
renv::restore(prompt = FALSE)
message("[00] renv environment restored.\n")

# --- 2. Run the numbered pipeline scripts in order ---------------------------
pipeline_scripts <- c(
  "01_cohort_identification.R",
  "02_quality_checks.R",
  "03_variable_derivation.R",
  "04_analysis.R",                 # outcome analyses (replication, survival, bias)
  "05_normalization_analysis.R"    # PBW vs PFVC normalization discordance + prognostics
)
# NOTE: cross-cohort pooling (code/pooled_estimates.R) is NOT part of the per-site
# pipeline. It is run centrally by the study coordinator after every site returns
# its `final/` outputs, and is kept local (not in the repository).

rscript_bin <- file.path(R.home("bin"), "Rscript")

for (script_name in pipeline_scripts) {
  script_file <- file.path("code", script_name)
  if (!file.exists(script_file)) {
    stop("Expected pipeline script not found: ", script_file)
  }

  message("=============================================================")
  message("[00] Running ", script_name, " ...")
  message("=============================================================")

  status <- system2(rscript_bin, args = shQuote(script_file))

  if (!identical(status, 0L)) {
    stop("Pipeline halted: ", script_name, " exited with status ", status,
         ". Fix the error above before re-running.")
  }
  message("[00] ", script_name, " completed successfully.\n")
}

message("=============================================================")
message("[00] Pipeline complete. All scripts ran successfully.")
message("Outputs are under output/<site_name>_output/.")
message("=============================================================")
