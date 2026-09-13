# =============================================================================
# Script 10 (engine loader): a config-aware cache/guard in front of 10_tte_common.R
# =============================================================================
# Iteration speed-up AND a correctness guard. The 11.* leaves build on the shared engine
# in 10_tte_common.R (parquet reads, the SpO2->FiO2 rolling join, the IPCW propensity
# fits, the primary `des`), slow to recompute on every standalone run. This loader
# snapshots the built engine to a gitignored .rds keyed by SITE + cohort window, restores
# it (sub-second) on later runs, and -- crucially -- re-resolves the site on EVERY source
# so a config switch is never missed.
#
# Source THIS instead of 10_tte_common.R, and source it UNCONDITIONALLY (do NOT wrap it in
# `if (!exists("build_design"))` -- that guard would pin a previous site's in-memory engine
# across a config.json switch in a persistent session). The loader self-guards:
#   * in-session no-op  -- engine already loaded for THIS exact config -> returns instantly
#   * cross-process hit -- valid .rds for this site/window -> loads it (sub-second)
#   * miss / changed    -- sources 10_tte_common.R and writes the snapshot
# A cache hit/rebuild is REBUILT (never silently reused) whenever the engine source,
# utils/config.R, ANY input parquet's mtime, the site, the PBWPFVC_TTE_VTPBW window, or the
# R version changes. A hard assertion stops the run if a restored engine's site ever fails
# to match config. PBWPFVC_TTE_NOCACHE=1 bypasses the cache entirely. Every decision prints
# the resolved SITE and cache PATH, so a cross-site mix-up is impossible to miss.
# =============================================================================

# The engine assumes single-threaded BLAS and a fixed library set; a cache hit skips
# 10_tte_common.R, so both must be established here too (run on every source -- idempotent).
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(arrow)
  library(here); library(splines); library(survival); library(parallel)
})
source("utils/config.R")   # re-read config.json EVERY time, so a site switch is seen

.tte_mt      <- function(f) if (file.exists(f)) format(file.info(f)$mtime, "%Y-%m-%d %H:%M:%OS6") else NA_character_
.tte_site    <- config$site_name
.tte_out     <- here("output", paste0(.tte_site, "_output"), "intermediate")
.tte_vtpfvc  <- Sys.getenv("PBWPFVC_TTE_VTPBW", "")
# Sensitivity knobs that change the BUILT engine (base/design) MUST key the cache AND the
# snapshot filename -- else a pfvc_age25 / tighter-ceiling run silently reloads the pfvc /
# default-ceiling engine and writes the wrong result (bug fixed 2026-06: NORM/ceiling were absent).
.tte_norm    <- Sys.getenv("PBWPFVC_TTE_NORM",  "pfvc")
.tte_clow    <- Sys.getenv("PBWPFVC_TTE_CLOW",  "11")
.tte_chigh   <- Sys.getenv("PBWPFVC_TTE_CHIGH", "16")
.tte_suffix  <- paste0(
  if (nzchar(.tte_vtpfvc)) paste0("_vtpbw_", gsub("[^0-9]+", "_", .tte_vtpfvc)) else "",
  if (!identical(.tte_norm, "pfvc")) paste0("_", .tte_norm) else "",
  if (!(identical(.tte_clow, "11") && identical(.tte_chigh, "16"))) paste0("_c", .tte_clow, "_", .tte_chigh) else "")
.tte_cache   <- file.path(.tte_out, paste0("tte_engine_cache_", .tte_site, .tte_suffix, ".rds"))
.tte_common  <- here("code", "10_tte_common.R")
.tte_panel   <- here("code", "10_panel_common.R")   # sourced by common; keys the cache too
.tte_nocache <- nzchar(Sys.getenv("PBWPFVC_TTE_NOCACHE"))

# Every input the engine reads; ANY change in their mtimes invalidates the cache.
.tte_inputs <- file.path(.tte_out, c(
  "analysis_cross_sectional.parquet", "analysis_all_eligible_timepoints.parquet",
  "resp_support_waterfall_clean.parquet", "cohort_vitals_clean.parquet",
  "cohort_meds.parquet", "cohort_labs_clean.parquet", "cohort_demographics.parquet",
  "ne_equiv_admin.parquet"))
.tte_key <- list(
  site = .tte_site, vtpfvc = .tte_vtpfvc,
  norm = .tte_norm, clow = .tte_clow, chigh = .tte_chigh,
  common_mtime = .tte_mt(.tte_common), panel_mtime = .tte_mt(.tte_panel),
  config_mtime = .tte_mt(here("utils", "config.R")),
  inputs = setNames(vapply(.tte_inputs, .tte_mt, character(1)), basename(.tte_inputs)),
  r_version = R.version.string)

# (0) In-session short-circuit: is the engine ALREADY in memory for THIS exact config?
if (exists(".tte_loaded_key", envir = .GlobalEnv) && exists("build_design", envir = .GlobalEnv) &&
    identical(get(".tte_loaded_key", envir = .GlobalEnv), .tte_key)) {
  message("tte engine: already loaded this session for site ", .tte_site, " -- reusing (no-op).")
} else {
  # Large raw scratch no leaf references (verified) -- excluded to keep the snapshot lean.
  .tte_blocklist <- c("wf", "gas", "fio2_dt", "spo2_dt", "sf_daily", "vit", "med", "panel_full")
  .tte_reason <- if (.tte_nocache) "PBWPFVC_TTE_NOCACHE set" else
    if (!file.exists(.tte_cache)) "no cache yet" else {
      .tte_snap <- tryCatch(readRDS(.tte_cache), error = function(e) NULL)
      if (is.null(.tte_snap) || is.null(.tte_snap$meta)) "cache unreadable/corrupt" else
      if (!identical(.tte_snap$meta, .tte_key)) "inputs changed (engine/config/data/site/window/R)" else NULL
    }
  .tte_t0 <- Sys.time()
  if (is.null(.tte_reason)) {
    list2env(.tte_snap$objs, envir = .GlobalEnv)
    message(sprintf("tte engine: LOADED FROM CACHE in %.1fs -- site %s (built %s)\n  %s",
            as.numeric(difftime(Sys.time(), .tte_t0, units = "secs")), .tte_site, .tte_snap$built_at, .tte_cache))
  } else {
    message("tte engine: REBUILDING site ", .tte_site, " (", .tte_reason, ")\n  -> ", .tte_cache)
    source(.tte_common)
    if (!.tte_nocache) {
      .tte_names <- setdiff(ls(.GlobalEnv, all.names = FALSE), .tte_blocklist)
      saveRDS(list(meta = .tte_key, built_at = as.character(Sys.time()),
                   objs = mget(.tte_names, envir = .GlobalEnv)), .tte_cache, compress = FALSE)
      message(sprintf("tte engine: built site %s in %.0fs, cached %d objects",
              .tte_site, as.numeric(difftime(Sys.time(), .tte_t0, units = "secs")), length(.tte_names)))
    }
  }
  # Hard guard against cross-site leakage: the restored/built engine MUST be this site.
  if (!identical(get0("site_name", envir = .GlobalEnv), .tte_site))
    stop("tte engine: restored site_name (", get0("site_name", envir = .GlobalEnv),
         ") != config site (", .tte_site, "). Refusing to proceed -- delete the cache .rds and rebuild.")
  assign(".tte_loaded_key", .tte_key, envir = .GlobalEnv)
}
