# =============================================================================
# Script 00: Run the full pipeline
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Single entry point for running this project at a CLIF consortium site.
#
#   1. Restores the project environment from renv.lock.
#   2. Runs the requested stages in order, each script as a clean subprocess.
#
# Everything a stage writes for sharing goes to output/<site>_output/final/, which
# holds aggregates only and is the folder a site returns. Re-running a stage updates
# that folder in place. Stages, by the manuscript figure they feed:
#
#   prep             01-03   cohort, quality checks, derived variables
#   cross_sectional  04-05   figures 1-3: PBW bias by demographics, respiratory
#                            mechanics, mortality
#   injury           20-28   figure 4: organ-injury markers over time (29_run_biotrauma.sh)
#   controls         21-27   figure 4's control cohorts, built together into
#                            final/controls/ (29_run_controls.sh build + anchors; the
#                            matched fits run too when SEV_MIN is set in the environment)
#   causal           30-38   figure 5: target trial emulation and the preference instrument
#
# The default is "prep,cross_sectional", which is what this runner has always done.
# The later stages take hours; ask for them by name, or --stages all.
#
# Cross-cohort pooling (code/pooling/) is a separate, centrally-run step and is
# intentionally not invoked here.
#
# Usage (from the project root, or anywhere — the script locates the repo):
#   Rscript code/00_run_pipeline.R
#   Rscript code/00_run_pipeline.R --site_name my_site --site_path /data/clif [--file_type parquet]
#   Rscript code/00_run_pipeline.R --stages injury,controls
#   Rscript code/00_run_pipeline.R --stages all
#   Rscript code/00_run_pipeline.R --analysis_only          (the same as --stages cross_sectional)
#
# --analysis_only skips the data-preparation scripts (01 cohort, 02 QC, 03 variables)
# and runs only the analyses (04, 05) against the script-03 outputs already on disk
# for the site -- for re-running the analyses after a code change without rebuilding
# the cohort. It fails loudly if those outputs are missing.
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
             file_type = "PBWPFVC_FILE_TYPE", stages = "stages")
  flags <- c("analysis_only")
  usage <- paste0("Usage: Rscript code/00_run_pipeline.R [--site_name NAME] [--site_path DIR] ",
                  "[--file_type parquet|csv|fst] [--stages prep,cross_sectional,injury,controls,causal|all] ",
                  "[--analysis_only]")
  out <- list(); i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% paste0("--", flags)) {                    # boolean flag
      out[[sub("^--", "", a)]] <- TRUE; i <- i + 1L; next
    }
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
cli_args <- parse_pipeline_args(commandArgs(trailingOnly = TRUE))
analysis_only <- isTRUE(cli_args$analysis_only); cli_args$analysis_only <- NULL
ALL_STAGES <- c("prep", "cross_sectional", "injury", "controls", "causal")
stages <- if (!is.null(cli_args$stages)) trimws(strsplit(cli_args$stages, ",")[[1]]) else
  if (analysis_only) "cross_sectional" else c("prep", "cross_sectional")
if (identical(stages, "all")) stages <- ALL_STAGES
if (!all(stages %in% ALL_STAGES))
  stop("Unknown stage: ", paste(setdiff(stages, ALL_STAGES), collapse = ", "), "; stages are ", paste(ALL_STAGES, collapse = ", "))
stages <- ALL_STAGES[ALL_STAGES %in% stages]            # always in pipeline order
analysis_only <- !("prep" %in% stages)
cli_args$stages <- NULL
cli_overrides <- cli_args
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
message("Stages: ", paste(stages, collapse = ", "))
if (analysis_only) message("No prep stage: running on the script-03 outputs already on disk")
message("=============================================================")

# --- 1. Restore the project environment --------------------------------------
message("\n[00] Restoring renv environment from renv.lock ...")
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}
renv::restore(prompt = FALSE)
message("[00] renv environment restored.\n")

# --- 2. Run the numbered pipeline scripts in order ---------------------------
# Each step is a command line: an R script, or one of the shell runners.
r_step  <- function(script) list(label = script, bin = "Rscript", args = shQuote(file.path("code", script)))
sh_step <- function(script, ...) list(label = paste(script, ...), bin = "bash", args = c(shQuote(file.path("code", script)), ...))
stage_steps <- list(
  prep = list(r_step("01_cohort_identification.R"), r_step("02_quality_checks.R"), r_step("03_variable_derivation.R")),
  cross_sectional = list(r_step("04_analysis.R"),                # bias, mechanics, mortality (figures 1-3)
                         r_step("05_normalization_analysis.R")), # PBW vs PFVC normalization of the injury metrics
  injury = list(sh_step("29_run_biotrauma.sh")),
  # the matched fits need a severity floor chosen from the anchors stage's output, so
  # they run only when SEV_MIN is already in the environment
  controls = c(list(sh_step("29_run_controls.sh", "build"), sh_step("29_run_controls.sh", "anchors")),
               if (nzchar(Sys.getenv("SEV_MIN", ""))) list(sh_step("29_run_controls.sh", "fits"))),
  causal = list(r_step("32_tte_run_all.R"), r_step("38_iv_preference.R"))
)
pipeline_steps <- unlist(stage_steps[stages], recursive = FALSE)
if (analysis_only) {
  # the analyses read script 03's outputs; refuse to start if they are not on disk
  site_for_dir <- if (!is.null(cli_overrides$PBWPFVC_SITE_NAME)) cli_overrides$PBWPFVC_SITE_NAME
                  else jsonlite::fromJSON("config/config.json")$site_name
  stage_dir <- file.path("output", paste0(site_for_dir, "_output"), "intermediate")
  needed <- file.path(stage_dir, c("analysis_cross_sectional.parquet", "analysis_negative_control.parquet"))
  missing <- needed[!file.exists(needed)]
  if (length(missing))
    stop("--analysis_only needs the script-03 outputs, which are missing:\n  ",
         paste(missing, collapse = "\n  "), "\nRun the full pipeline first.")
}
# NOTE: cross-cohort pooling (code/pooling/pooled_estimates.R) is NOT part of the per-site
# pipeline. It is run centrally by the study coordinator after every site returns
# its `final/` outputs, and is kept local (not in the repository).

rscript_bin <- file.path(R.home("bin"), "Rscript")

for (step in pipeline_steps) {
  script_name <- step$label
  message("=============================================================")
  message("[00] Running ", script_name, " ...")
  message("=============================================================")

  status <- system2(if (step$bin == "Rscript") rscript_bin else step$bin, args = step$args)

  if (!identical(status, 0L)) {
    stop("Pipeline halted: ", script_name, " exited with status ", status,
         ". Fix the error above before re-running.")
  }
  message("[00] ", script_name, " completed successfully.\n")
}
if ("controls" %in% stages && !nzchar(Sys.getenv("SEV_MIN", "")))
  message("[00] controls: cohorts built and anchor distributions written. Choose the severity floor from\n",
          "     final/injury/jm_severity_anchor_*, then: SEV_MIN=\"platelets=2,bilirubin=1\" bash code/29_run_controls.sh fits")

message("=============================================================")
message("[00] Pipeline complete. All scripts ran successfully.")
message("Shareable aggregates: output/<site_name>_output/final/, sorted into cross_sectional/, injury/, causal/, supplement/ and controls/.")
message("=============================================================")
