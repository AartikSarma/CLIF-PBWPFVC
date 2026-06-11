# =============================================================================
# Script 00: Run the full pipeline
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Single entry point for running this project at a CLIF consortium site.
#
#   1. Restores the project environment from renv.lock.
#   2. Runs scripts 01-05 in order, each as a clean R subprocess.
#
# Usage (from the project root, or anywhere — the script locates the repo):
#   Rscript code/00_run_pipeline.R
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

message("=============================================================")
message("PBW vs PFVC pipeline runner")
message("Repository root: ", repo_root)
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
  "04_analysis.R",
  "04b_elastance_age_interaction.R",
  "05_cross_cohort_forest.R"
)

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
