# =============================================================================
# Find which step of the causal scripts crashes R (the macOS segfault)
# =============================================================================
#
# Repeats the reading steps of causal_triangulation.R several times, printing a line
# before each, so the last line printed before a crash names the step. Prints only
# step names, row counts and thread settings: no patient data. Writes nothing.
#
# Usage: uvr run code/tools/crash_probe.R
#        PBWPFVC_PROBE_REPS=10 uvr run code/tools/crash_probe.R
# Run it without the thread variables to reproduce the crash, then with one of them
# set at a time (for example OMP_NUM_THREADS=1) to see which one stops it.
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(here) })
source("utils/config.R")
source(here::here("code", "20_biotrauma_grid.R"))
source(here::here("code", "supplement", "causal_common.R"))

reps <- as.integer(Sys.getenv("PBWPFVC_PROBE_REPS", "5"))
step <- function(what) { message(format(Sys.time(), "%H:%M:%S"), "  ", what); flush(stderr()) }
step(paste0("arrow CPU threads ", arrow::cpu_count(), ", I/O threads ", arrow::io_thread_count(),
            "; data.table loaded: ", "data.table" %in% loadedNamespaces()))
for (r in seq_len(reps)) {
  step(paste0("pass ", r, ": the patient tables"))
  patients <- load_causal_patients(config)
  for (f in both_cohorts(config, "cohort_labs_clean.parquet")) {
    step(paste0("pass ", r, ": read_parquet ", basename(dirname(f)), "/", basename(f)))
    labs <- arrow::read_parquet(f)
    step(paste0("pass ", r, ":   ", nrow(labs), " rows read"))
    rm(labs); invisible(gc())
  }
  step(paste0("pass ", r, ": the platelet panel"))
  panel <- load_causal_outcome(config, patients, "platelets")
  step(paste0("pass ", r, ":   ", nrow(panel), " rows"))
  step(paste0("pass ", r, ": the designs"))
  fits <- fit_all_designs(causal_designs(TRUE), patients, panel)
  step(paste0("pass ", r, ":   ", sum(!is.na(fits$estimate)), " designs fitted"))
}
step(paste0("no crash in ", reps, " passes; data.table loaded: ", "data.table" %in% loadedNamespaces()))
