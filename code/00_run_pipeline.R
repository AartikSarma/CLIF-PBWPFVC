# =============================================================================
# Script 00: Run the full pipeline
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Single entry point for running this project at a CLIF consortium site.
#
#   1. Installs the packages pinned in uvr.lock (uvr sync).
#   2. Runs the requested stages in order, each script as a clean subprocess.
#
# Everything a stage writes for sharing goes to output/<site>_output/final/, which
# holds aggregates only and is the folder a site returns. Re-running a stage updates
# that folder in place. Stages, by the manuscript figure they feed:
#
#   prep             01-03   cohort, quality checks, derived variables; then 03 again
#                            for the ungated cohort (below)
#   cross_sectional  04-05   figures 1-3: PBW bias by demographics, respiratory
#                            mechanics, mortality; then three supplement scripts, all
#                            written to final/supplement/: xsec_crs_channels.R (does
#                            measured compliance scale like PFVC or PBW; figure 3C),
#                            xsec_dp_vtpfvc_additive.R (driving pressure and VT/PFVC
#                            as additive mortality risks) and xsec_mortality_prediction.R
#                            (which dose or mechanics measure predicts death best);
#                            then, gated pass only, the gate check
#                            (xsec_vtpbw_gate_collider.R) and what precedes a tidal-volume
#                            reduction (xsec_vt_reduction_timing.R)
#   injury           20-29   figure 4: organ-injury markers over 7 days in the ventilated
#                            cohort and the no-support controls (all, and hypoxemic with
#                            SF < 315), each control standardised to the ventilated
#                            cohort's severity rather than matched (29_run_figure4.R:
#                            builds the control cohort, runs the ventilated-vs-control
#                            mortality contrasts, fits, draws)
#
# Two cohorts, both reported. The paper's cohort enters a ventilated patient at a
# hypoxemic timepoint with VT/PBW 6-8 mL/kg. VT/PBW follows both predicted lung size
# (a round-number volume is more mL/kg in a small patient) and illness (clinicians turn
# the volume down in sicker patients), so selecting on it can link the two
# (supplement/xsec_vtpbw_gate_collider.R). Every cross-sectional script therefore runs
# twice: on the gated cohort (final/<block>/) and on the ungated one, the same rule
# without the band (PBWPFVC_VTPBW_GATE=0; final/ungated/<block>/, same file names).
# The runner sets the variable itself; a site never sets it by hand. The ungated pass of
# 03 repeats its derivation (sections 3a-3d) and writes only its cross-sectional table.
#
# The default is "prep,cross_sectional". The injury stage takes hours; ask for it by
# name, or --stages all.
#
# Cohort switch: PBWPFVC_COHORT=imv (default; the paper's cohort) | nosupport
# (negative control) | niv (built on request only). Environment variable only; it
# selects the cohort that scripts 01-03 build, e.g.
#   PBWPFVC_COHORT=nosupport uvr run code/01_cohort_identification.R
#
# Cross-cohort pooling (code/pooling/) is a separate, centrally-run step and is
# intentionally not invoked here.
#
# Usage (from anywhere inside the repository, where uvr finds uvr.toml):
#   uvr run code/00_run_pipeline.R
#   uvr run code/00_run_pipeline.R -- --site_name my_site --site_path /data/clif [--file_type parquet]
#   uvr run code/00_run_pipeline.R -- --stages injury
#   uvr run code/00_run_pipeline.R -- --stages all
#   uvr run code/00_run_pipeline.R -- --analysis_only          (the same as --stages cross_sectional)
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
#   PBWPFVC_SITE_NAME=my_site uvr run code/04_analysis.R
#
# Each numbered script is standalone: it reads its inputs from disk and writes
# its outputs back to disk, so they are run as separate subprocesses rather than
# sourced. This isolates package namespaces and means script 01's
# `rm(list = ls())` cannot wipe state belonging to this runner. If any script
# exits with an error, the pipeline stops and reports which one failed.
# =============================================================================

# --- Locate the repository root ----------------------------------------------
# Find this file's path from the R invocation so the pipeline can be
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
       "Run it with: uvr run code/00_run_pipeline.R")
}

script_path <- get_script_path()
repo_root   <- normalizePath(file.path(dirname(script_path), ".."))
setwd(repo_root)

# --- Command-line overrides of config.json -----------------------------------
parse_pipeline_args <- function(args) {
  known <- c(site_name = "PBWPFVC_SITE_NAME", site_path = "PBWPFVC_TABLES_PATH",
             file_type = "PBWPFVC_FILE_TYPE", stages = "stages")
  flags <- c("analysis_only")
  usage <- paste0("Usage: uvr run code/00_run_pipeline.R -- [--site_name NAME] [--site_path DIR] ",
                  "[--file_type parquet|csv|fst] [--stages prep,cross_sectional,injury|all] ",
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
ALL_STAGES <- c("prep", "cross_sectional", "injury")
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

# --- 1. Sync the project environment -----------------------------------------
# uvr (https://github.com/nbafrank/uvr) installs the packages pinned in uvr.lock
# into .uvr/library/. `uvr run` puts that library first on this process's library
# path through R_LIBS_USER, which every script started below inherits, so the
# runner must be launched with `uvr run`, not bare Rscript.
message("\n[00] Syncing packages from uvr.lock ...")
if (!nzchar(Sys.which("uvr")))
  stop("uvr is not installed; see https://github.com/nbafrank/uvr#installation")
if (!identical(system2("uvr", "sync"), 0L))
  stop("uvr sync failed; fix the error above before re-running.")
uvr_library <- normalizePath(file.path(repo_root, ".uvr", "library"))
if (!uvr_library %in% normalizePath(.libPaths()))
  stop("The project library ", uvr_library, " is not on the library path. ",
       "Launch the runner with: uvr run code/00_run_pipeline.R")
message("[00] Packages synced.\n")

# --- 2. Run the numbered pipeline scripts in order ---------------------------
# Each step is an R script in code/ and the cohort it runs on: gate "1" is the paper's
# cohort, "0" the ungated one (PBWPFVC_VTPBW_GATE, set for that step only).
cross_sectional_scripts <- c(
  "04_analysis.R",                             # bias, mechanics, mortality (figures 1-3)
  "05_normalization_analysis.R",               # PBW vs PFVC normalization of the injury metrics
  "supplement/xsec_crs_channels.R",            # compliance against PFVC and PBW (figure 3C)
  "supplement/xsec_dp_vtpfvc_additive.R",      # driving pressure and VT/PFVC, additive risks
  "supplement/xsec_mortality_prediction.R")    # which measure predicts death best
gate_check_scripts <- c(                       # read both cohorts themselves; gated pass only
  "supplement/xsec_vtpbw_gate_collider.R",     # does the 6-8 gate link predicted size to illness?
  "supplement/xsec_vt_reduction_timing.R")     # what precedes a tidal-volume reduction
step <- function(scripts, gate) data.frame(script = scripts, gate = gate, stringsAsFactors = FALSE)
stage_steps <- list(
  prep = rbind(step(c("01_cohort_identification.R", "02_quality_checks.R", "03_variable_derivation.R"), "1"),
               step("03_variable_derivation.R", "0")),
  cross_sectional = rbind(step(c(cross_sectional_scripts, gate_check_scripts), "1"),
                          step(cross_sectional_scripts, "0")),
  injury = step("29_run_figure4.R", "1")      # sets the gate for each of its own steps
)
pipeline_steps <- do.call(rbind, stage_steps[stages])
if (analysis_only) {
  # the analyses read script 03's outputs; refuse to start if they are not on disk
  site_for_dir <- if (!is.null(cli_overrides$PBWPFVC_SITE_NAME)) cli_overrides$PBWPFVC_SITE_NAME
                  else jsonlite::fromJSON("config/config.json")$site_name
  stage_dir <- file.path("output", paste0(site_for_dir, "_output"), "intermediate")
  needed <- file.path(stage_dir, c("analysis_cross_sectional.parquet", "analysis_cross_sectional_ungated.parquet",
                                   "analysis_ungated_index.parquet", "analysis_negative_control.parquet"))
  missing <- needed[!file.exists(needed)]
  if (length(missing))
    stop("--analysis_only needs the script-03 outputs, which are missing:\n  ",
         paste(missing, collapse = "\n  "), "\nRun the full pipeline first.")
}
# NOTE: cross-cohort pooling (code/pooling/) is NOT part of the per-site pipeline. It is
# run centrally by the study coordinator after every site returns its `final/` outputs.

# --- Clear the cross-sectional block before it is rewritten --------------------
# Scripts 03-05 rewrite final/cross_sectional/ whole, so a file an older version of
# the code wrote there (and this one no longer does) would otherwise travel with the
# site's return. Cleared: everything when prep runs; everything but script 03's own
# tables when only the analyses run. final/injury/ is never cleared: its tables are
# merged on write and cache the joint-model fits. The control cohorts' folder
# (final/controls/) is not touched. The ungated cohort's block
# (final/ungated/cross_sectional/) is cleared whole: 03 writes nothing there.
if (Sys.getenv("PBWPFVC_COHORT", "imv") == "imv" && any(c("prep", "cross_sectional") %in% stages)) {
  site_for_dir <- if (!is.null(cli_overrides$PBWPFVC_SITE_NAME)) cli_overrides$PBWPFVC_SITE_NAME
                  else jsonlite::fromJSON("config/config.json")$site_name
  final_root <- file.path("output", paste0(site_for_dir, "_output"), "final")
  block_dir <- file.path(final_root, "cross_sectional")
  old_files <- list.files(block_dir, full.names = TRUE)
  if (!("prep" %in% stages)) old_files <- old_files[!grepl("^(attrition_log_|dist_)", basename(old_files))]
  old_files <- c(old_files, list.files(file.path(final_root, "ungated", "cross_sectional"), full.names = TRUE))
  if (length(old_files)) {
    message("[00] Clearing ", length(old_files), " file(s) from ", block_dir, " and its ungated twin before the rerun")
    unlink(old_files)
  }
}

rscript_bin <- file.path(R.home("bin"), "Rscript")

# PBWPFVC_VTPBW_GATE is set for each child and removed afterwards, so it never leaks
# into a later step (utils/config.R stops if it meets a control cohort)
for (i in seq_len(nrow(pipeline_steps))) {
  script_name <- pipeline_steps$script[i]
  gate <- pipeline_steps$gate[i]
  step_label <- paste0(script_name, if (gate == "0") " (ungated cohort)" else "")
  message("=============================================================")
  message("[00] Running ", step_label, " ...")
  message("=============================================================")

  if (gate == "0") Sys.setenv(PBWPFVC_VTPBW_GATE = "0") else Sys.unsetenv("PBWPFVC_VTPBW_GATE")
  status <- system2(rscript_bin, shQuote(file.path("code", script_name)))
  Sys.unsetenv("PBWPFVC_VTPBW_GATE")

  if (!identical(status, 0L)) {
    stop("Pipeline halted: ", step_label, " exited with status ", status,
         ". Fix the error above before re-running.")
  }
  message("[00] ", step_label, " completed successfully.\n")
}

message("=============================================================")
message("[00] Pipeline complete. All scripts ran successfully.")
message("Shareable aggregates: output/<site_name>_output/final/, sorted into cross_sectional/, injury/, supplement/ and controls/, with the ungated cohort's twins in ungated/.")
message("=============================================================")
