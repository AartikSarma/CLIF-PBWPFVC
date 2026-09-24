# =============================================================================
# 29_run_figure4: every analysis behind manuscript figure 4, then the figure
# =============================================================================
# Figure 4 asks whether a smaller predicted lung, at the same VT/PBW, goes with
# organ-injury markers that worsen faster over the first 7 days of ventilation. This
# script produces it for one site, from the CLIF tables to the PDF:
#
#   markers   platelets, bilirubin, creatinine (dialysis, continuous or intermittent,
#             as a third competing cause; ESRD censored at day 0), vasopressor dose on
#             pressor days, oxygen saturation index, and the SF ratio (panel C, the
#             positive control for mechanics: a larger VT/PFVC recruits lung, so SF
#             should be better early in the smaller predicted lung, then reverse)
#   arms      ventilated, all patients           (panels A, B and C)
#             no respiratory support             (panel B, the negative control): every
#                                                patient, the divergence read at the
#                                                ventilated cohort's mean severity
#             no support, hypoxemic on the index day (index-day SF <= 315, the ventilated
#                                                cohort's own gate), read at the ventilated severity:
#                                                the arms then differ in ventilation, not
#                                                hypoxemia (27 writes its DiD separately;
#                                                HYPOXEMIC_CONTROL_MARKERS, the control's markers by default)
#             optional: ventilated by baseline SF class (SF_BANDS; off by default)
#   The control is standardised, not matched (2026-09-21): severity cannot confound a
#   PFVC fixed by height, age, sex and race, but it could MODIFY the divergence, so the
#   control's divergence varies with its severity anchor and is read at the ventilated
#   mean (20_biotrauma_grid.R, PBWPFVC_JM_SEV_CENTER). No control patient is discarded,
#   and the severity x divergence term tests whether sicker controls diverge faster.
#   each fit adjusted and unadjusted for age, sex and race
#
# Steps, each logged to output/{site}_output/logs/figure4_{stamp}/:
#   1 build     scripts 01-03 for the ventilated cohort and the no-support cohort,
#               each only if its derived tables are missing (FORCE_BUILD=1 rebuilds)
#   2 panels    the 7-day panel of both cohorts
#   3 anchors   the severity-anchor distributions of both cohorts, and the ventilated
#               mean anchor per marker (final/injury/jm_severity_anchor_mean_*)
#   4 centres   each control marker's centre = that ventilated mean
#   5 fits      22_biotrauma_fit.R and 23_biotrauma_report.R for every arm
#   6 figure    supplement/xsec_pfvc_age_control.R (death against each control, the
#               checks figure's mortality rows), 27_control_comparison.R (the
#               difference-in-differences), then 24_biotrauma_figures.R: figure 4 and
#               the checks figure (biotrauma_fig_checks_*: every outcome against each
#               control)
#   7 channels  the supplement's channel breakdown: the ventilated divergence read through
#               each GLI piece of log PFVC (height, age, sex, race; 22's channels form, one
#               fit per marker, the pieces standing in for the demographics), and its
#               figure (biotrauma_fig_channels_*). CHANNEL_MARKERS, platelets by default;
#               empty skips it
# Fits already on disk with the same chain settings are reused, so a rerun after a
# failure costs only what failed. A failed step is reported and the run continues.
#
# The oxygen saturation index needs positive-pressure ventilation (mean airway
# pressure), so it has no control arm. A marker or arm with too few patients or
# deaths is skipped by the fit, with the reason in its manifest.
#
# Usage (from anywhere inside the repository, where uvr finds uvr.toml). Launch it with `uvr run`:
# each step is a child Rscript of the same R, and it finds the project library only
# through the library path `uvr run` sets.
#   caffeinate -i nohup uvr run code/29_run_figure4.R > figure4.out 2>&1 &
#   uvr run code/29_run_figure4.R -- --dry-run
# Site default: 28 figure-4 fits (5 markers and creatinine in the ventilated cohort; 3 markers
# and creatinine in each control, the whole one and the hypoxemic one; each adjusted and
# unadjusted) plus 1 channel fit, at
# 2000 / 500 iterations, four at a time. At
# MIMIC the divergence terms the figure rests on converged at 2000 iterations; the
# hazard blocks did not converge at any length tried; the survival submodel was simplified
# on 2026-09-24 for that reason (22_biotrauma_fit.R, HAZARD_SPEC), so every fit refits once.
# Knobs (environment): ITER BURNIN CHAINS THIN (2000 / 500 / 3 / 5), PAR (fits at a
#   time, 4; an earlier estimate put a 7-day fit at a 7,000-patient site at 15-25 GB,
#   so four at once can need 60-100 GB: lower PAR on a smaller machine), MARKERS, CONTROL_MARKERS, CREATININE (1; 0 skips it),
#   SF_BANDS (baseline SF classes of the ventilated cohort, off by default; the lead
#   site runs SF_BANDS="235,315 115,235 0,115"), HYPOXEMIC_CONTROL_MARKERS (the hypoxemic
#   control arm, CONTROL_MARKERS by default, with creatinine as in the other arms; empty
#   skips it),
#   CHANNEL_MARKERS (step 7, platelets by default; empty skips it),
#   FORCE_BUILD, FORCE_PANEL. The VT/PFVC companion (form vtpfvc) was dropped on
#   2026-09-24: at a fixed VT/PBW, VT/PFVC moves only with PBW/PFVC, whose variance
#   after age, sex and race is a few percent (Claim 5a), so the companion re-reads the
#   pfvc form on a scale the data cannot identify. 22-24 still accept the form.
#   Output: final/injury/biotrauma_fig_main_pfvc_7d_{site}.pdf.
# =============================================================================

# --- Locate the repository root ------------------------------------------------
get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1) stop("Run this script with: uvr run code/29_run_figure4.R")
  normalizePath(sub("^--file=", "", file_arg))
}
setwd(normalizePath(file.path(dirname(get_script_path()), "..")))

# --- Knobs -----------------------------------------------------------------------
# knob(): an unset OR empty variable takes the default. knob_or_skip(): only an unset
# variable does; set empty, it switches its step off.
knob <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}
knob_or_skip <- function(name, default) Sys.getenv(name, unset = default)

MARKERS                   <- knob("MARKERS", "osi,sf,pressor_dose,platelets,bilirubin")   # creatinine runs on its own, with RRT as a third cause
CONTROL_MARKERS           <- knob("CONTROL_MARKERS", "pressor_dose,platelets,bilirubin")
HYPOXEMIC_CONTROL_MARKERS <- knob_or_skip("HYPOXEMIC_CONTROL_MARKERS", CONTROL_MARKERS)   # the hypoxemic control arm
CHANNEL_MARKERS           <- knob_or_skip("CHANNEL_MARKERS", "platelets")             # step 7, the channel breakdown (supplement)
CREATININE  <- knob("CREATININE", "1") == "1"
SF_BANDS    <- strsplit(trimws(knob("SF_BANDS", "")), "[[:space:]]+")[[1]]   # e.g. "235,315 115,235 0,115"; off by default
ITER   <- knob("ITER", "2000"); BURNIN <- knob("BURNIN", "500"); CHAINS <- knob("CHAINS", "3")
THIN   <- knob("THIN", "5");    PAR    <- knob("PAR", "4")
FORCE_BUILD <- knob("FORCE_BUILD", "0") == "1"
FORCE_PANEL <- knob("FORCE_PANEL", "0") == "1"
DRY <- identical(commandArgs(trailingOnly = TRUE)[1], "--dry-run")

# the site's name, from the environment or config.json, without a control cohort's suffix
BASE_SITE <- knob("PBWPFVC_SITE_NAME", jsonlite::fromJSON("config/config.json")$site_name)
BASE_SITE <- sub("_nosupport$", "", sub("_niv$", "", BASE_SITE))
if (!nzchar(BASE_SITE)) stop("could not read site_name from config/config.json")
ROOT    <- file.path("output", paste0(BASE_SITE, "_output"))
LOG_DIR <- file.path(ROOT, "logs", paste0("figure4_", format(Sys.time(), "%Y%m%d_%H%M%S")))
if (!DRY) dir.create(LOG_DIR, recursive = TRUE)

# settings every step shares; each step's own settings are added around its run
Sys.unsetenv("PBWPFVC_COHORT")
Sys.setenv(PBWPFVC_SITE_NAME = BASE_SITE,
           PBWPFVC_JM_GRID = "daily", PBWPFVC_JM_HORIZON = "7", PBWPFVC_JM_MODIFIER = "pfvc",
           PBWPFVC_JM_MODELS = "main",
           PBWPFVC_JM_ITER = ITER, PBWPFVC_JM_BURNIN = BURNIN, PBWPFVC_JM_CHAINS = CHAINS,
           PBWPFVC_JM_THIN = THIN, PBWPFVC_JM_PAR = PAR)
message("site ", BASE_SITE, "; markers ", MARKERS, if (CREATININE) ",creatinine", "; controls ", CONTROL_MARKERS,
        "; SF bands '", if (length(SF_BANDS)) paste(SF_BANDS, collapse = " ") else "none", "'; chains ",
        ITER, "/", BURNIN, " x ", CHAINS)

if (!DRY) {
  message("logs -> ", LOG_DIR)
  # the steps are child Rscripts, which see the project library only if this process does
  uvr_library <- normalizePath(".uvr/library", mustWork = FALSE)
  if (!identical(system2("uvr", "sync"), 0L)) stop("uvr sync failed")
  if (!uvr_library %in% normalizePath(.libPaths()))
    stop("The project library ", uvr_library, " is not on the library path. ",
         "Launch the script with: uvr run code/29_run_figure4.R")
}

# --- Running a step --------------------------------------------------------------
rscript_bin <- file.path(R.home("bin"), "Rscript")
timestamp   <- function() format(Sys.time(), "%H:%M:%S")
FAILED <- character(0)

# One R script as a child process, its output in LOG_DIR/<step_name>.log. The cohort
# (imv | nosupport) and the step's own settings are set in the environment for the
# child and restored afterwards. A failure is recorded and the run continues.
run_step <- function(step_name, script, cohort = "imv", step_env = character(0)) {
  step_env <- c(if (cohort != "imv") c(PBWPFVC_COHORT = cohort), step_env)
  if (DRY) {
    message("[dry] ", step_name, ": ", paste(c(if (length(step_env)) paste0(names(step_env), "=", step_env), script), collapse = " "))
    return(invisible(TRUE))
  }
  message("[", timestamp(), "] ", step_name)
  previous <- if (length(step_env)) Sys.getenv(names(step_env), unset = NA, names = TRUE) else character(0)
  if (length(step_env)) do.call(Sys.setenv, as.list(step_env))
  log_file <- file.path(LOG_DIR, paste0(step_name, ".log"))
  status <- system2(rscript_bin, shQuote(script), stdout = log_file, stderr = log_file)
  for (name in names(previous))
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(previous[name]))
  if (identical(status, 0L)) {
    message("    ok")
  } else {
    FAILED <<- c(FAILED, step_name)
    message("    FAILED, see ", log_file)
    message(paste0("    | ", tail(readLines(log_file, warn = FALSE), 5), collapse = "\n"))
  }
  invisible(identical(status, 0L))
}

# one arm: the markers, then creatinine with dialysis as a third competing cause
fit_arm <- function(arm, cohort, markers, arm_env = character(0), with_creatinine = CREATININE) {
  # a cohort whose panel failed to build is not fitted: its old panel is out of date
  if (paste0("panel_", cohort) %in% FAILED) {
    message("[", timestamp(), "] ", arm, ": skipped, the ", cohort, " panel failed to build")
    FAILED <<- c(FAILED, paste0(arm, "_skipped"))
    return(invisible())
  }
  marker_env <- c(arm_env, PBWPFVC_JM_MARKERS = markers)
  run_step(paste0(arm, "_fit"),    "code/22_biotrauma_fit.R",    cohort, marker_env)
  run_step(paste0(arm, "_report"), "code/23_biotrauma_report.R", cohort, marker_env)
  if (with_creatinine) {
    creatinine_env <- c(arm_env, PBWPFVC_JM_MARKERS = "creatinine", PBWPFVC_JM_RRT_EVENT = "1")
    run_step(paste0(arm, "_creatinine_fit"),    "code/22_biotrauma_fit.R",    cohort, creatinine_env)
    run_step(paste0(arm, "_creatinine_report"), "code/23_biotrauma_report.R", cohort, creatinine_env)
  }
}

# ---- 1 build
build_cohort <- function(cohort, derived) {   # cohort, folder holding its derived tables
  # a cohort built before dialysis and ESRD entered the RRT definition (2026-09-21) lacks
  # rrt_sources_available.rds, and its panel cannot be built: rebuild it. A control built
  # before it was indexed at ICU admission (same day) lacks cohort_icu_stays.parquet.
  icu_ok <- cohort != "nosupport" || file.exists(file.path(derived, "cohort_icu_stays.parquet"))
  if (!FORCE_BUILD && icu_ok && file.exists(file.path(derived, "analysis_cross_sectional.parquet")) &&
      file.exists(file.path(derived, "rrt_sources_available.rds"))) {
    message("[", timestamp(), "] ", cohort, ": cohort already built, scripts 01-03 skipped (FORCE_BUILD=1 rebuilds)")
    return(invisible())
  }
  for (script in c("01_cohort_identification", "02_quality_checks", "03_variable_derivation"))
    run_step(paste0("build_", cohort, "_", script), file.path("code", paste0(script, ".R")), cohort)
}
build_cohort("imv",       file.path(ROOT, "intermediate"))
build_cohort("nosupport", file.path(ROOT, "intermediate", "controls", "nosupport"))

# ---- 2 panels, rebuilt only when something they are built from has changed: a fit made
#      on an older panel is refitted (22_biotrauma_fit.R compares the times), so an
#      unconditional rebuild would refit everything on every rerun. FORCE_PANEL=1 rebuilds.
build_panel <- function(cohort, derived) {    # cohort, folder holding its derived tables
  panel <- file.path(derived, "jm_surv_7d.parquet")
  if (!FORCE_PANEL && !DRY && file.exists(panel)) {
    inputs <- c("code/21_biotrauma_panel.R", "code/10_panel_common.R", "code/20_biotrauma_grid.R", "utils/config.R",
                file.path(derived, c("analysis_cross_sectional.parquet", "cohort_dialysis.parquet", "cohort_esrd.parquet")))
    stale <- inputs[!file.exists(inputs) | file.mtime(inputs) > file.mtime(panel)]
    if (!length(stale)) {
      message("[", timestamp(), "] panel_", cohort, ": up to date, kept (FORCE_PANEL=1 rebuilds)")
      return(invisible())
    }
    message("[", timestamp(), "] panel_", cohort, ": ", basename(stale[1]), " is newer than the panel (or missing); rebuilding")
  }
  run_step(paste0("panel_", cohort), "code/21_biotrauma_panel.R", cohort)
}
build_panel("imv",       file.path(ROOT, "intermediate"))
build_panel("nosupport", file.path(ROOT, "intermediate", "controls", "nosupport"))

# ---- 3 anchors
ANCHOR_MARKERS <- paste0("creatinine,", CONTROL_MARKERS)
anchor_env <- c(PBWPFVC_JM_ANCHOR_ONLY = "1", PBWPFVC_JM_MARKERS = ANCHOR_MARKERS)
run_step("anchors_ventilated", "code/22_biotrauma_fit.R", "imv",       anchor_env)
run_step("anchors_nosupport",  "code/22_biotrauma_fit.R", "nosupport", anchor_env)

# ---- 4 centres: the ventilated cohort's mean anchor per control marker, where each
#      control fit reads its divergence (PBWPFVC_JM_SEV_CENTER, 20_biotrauma_grid.R)
CENTER_FILE <- file.path(ROOT, "final", "injury", paste0("jm_severity_anchor_mean_7d_", BASE_SITE, ".csv"))
anchor_markers <- strsplit(ANCHOR_MARKERS, ",")[[1]]
if (DRY) {
  message("[dry] centres: ventilated mean anchor per marker (", ANCHOR_MARKERS, ") from ", CENTER_FILE)
  SEV_CENTER <- "creatinine=M,platelets=M,..."
} else {
  # the file is written by anchors_ventilated; one "marker=mean" pair per control marker
  anchor_means <- if (file.exists(CENTER_FILE)) read.csv(CENTER_FILE) else data.frame(marker = character(0), anchor_mean = numeric(0))
  anchor_means <- anchor_means[anchor_means$marker %in% anchor_markers, ]
  SEV_CENTER <- if (nrow(anchor_means)) paste0(anchor_means$marker, "=", sprintf("%.4f", anchor_means$anchor_mean), collapse = ",") else ""
  if (!setequal(anchor_means$marker, anchor_markers)) {
    message("severity centres missing for some of ", ANCHOR_MARKERS, " (found '", SEV_CENTER, "' in ",
            CENTER_FILE, "); the control arm is skipped")
    FAILED <- c(FAILED, "centres"); SEV_CENTER <- ""
  } else {
    message("[", timestamp(), "] severity centres (ventilated mean anchor): ", SEV_CENTER)
  }
}

# ---- 5 fits, arm by arm
fit_arm("ventilated", "imv", MARKERS)
for (band in SF_BANDS)
  fit_arm(paste0("ventilated_sf", sub(",", "to", band)), "imv", MARKERS, c(PBWPFVC_JM_SF_BAND = band))
if (nzchar(SEV_CENTER)) {
  fit_arm("nosupport", "nosupport", CONTROL_MARKERS, c(PBWPFVC_JM_SEV_CENTER = SEV_CENTER))
  # the hypoxemic control: the same control, index-day SF <= 315, with the same markers
  if (nzchar(HYPOXEMIC_CONTROL_MARKERS))
    fit_arm("nosupport_hypoxemic", "nosupport", HYPOXEMIC_CONTROL_MARKERS,
            c(PBWPFVC_JM_SEV_CENTER = SEV_CENTER, PBWPFVC_JM_SF_BAND = "0,315"))
}

# ---- 6 comparison table and the figure
FIG_MARKERS <- paste0("platelets,bilirubin", if (CREATININE) ",creatinine", ",pressor_dose,osi,sf")
# death against each control: the PFVC association with death in the ventilated cohort
# and in the no-support control, everyone and hypoxemic at the index (the checks
# figure's mortality rows); it reads both cohorts' tables and the control's 7-day panel
run_step("mortality_contrast", "code/supplement/xsec_pfvc_age_control.R")
run_step("comparison", "code/27_control_comparison.R")
run_step("figure", "code/24_biotrauma_figures.R", "imv",
         c(PBWPFVC_JM_WITH_RRT = "1", PBWPFVC_FIG_MARKERS = FIG_MARKERS))

# ---- 7 the channel breakdown (supplement), after figure 4 so it cannot delay it: the
#      ventilated divergence per log unit of each GLI piece, pooled by pooled_biotrauma.R
#      from jm_level_contrast_channels_* (exposures ch_height, ch_age, ch_sex, ch_race)
if (nzchar(CHANNEL_MARKERS) && !"panel_imv" %in% FAILED) {
  channel_env <- c(PBWPFVC_JM_MODIFIER = "channels", PBWPFVC_JM_MARKERS = CHANNEL_MARKERS)
  run_step("channels_fit",    "code/22_biotrauma_fit.R",     "imv", channel_env)
  run_step("channels_report", "code/23_biotrauma_report.R",  "imv", channel_env)
  run_step("channels_figure", "code/24_biotrauma_figures.R", "imv", c(PBWPFVC_JM_MODIFIER = "channels"))
}

if (!DRY) {
  message("\nfigure  -> ", file.path(ROOT, "final", "injury", paste0("biotrauma_fig_main_pfvc_7d_", BASE_SITE, ".pdf")))
  message("tables  -> ", file.path(ROOT, "final", "injury"), "/ and ", file.path(ROOT, "final", "controls"), "/")
  if (length(FAILED)) {
    message("FAILED steps (", length(FAILED), "): ", paste(FAILED, collapse = " "))
    quit(status = 1)
  }
  message("all steps ok")
}
