# =============================================================================
# 29_run_figure4: every analysis behind manuscript figure 4, then the figure
# =============================================================================
# Figure 4 asks whether a smaller predicted lung, at the same VT/PBW, goes with
# organ-injury markers that worsen faster over the first 7 days of ventilation. This
# script produces it for one site, from the CLIF tables to the PDF:
#
#   markers   platelets, bilirubin, creatinine (dialysis, continuous or intermittent,
#             as a third competing cause; ESRD censored at day 0), vasopressors,
#             oxygen saturation index, and the SF ratio (a row of the figure; the
#             positive control for mechanics: a larger VT/PFVC recruits lung, so SF
#             should be better early in the smaller predicted lung, then reverse)
#   vasopressors are a two-part (hurdle) outcome, reported as a pair: on/off
#             (any_pressor, every patient-day, a logistic mixed model) and the dose
#             on the days a pressor runs (pressor_dose, conditional on being on a
#             pressor that day). Being on a pressor is itself an outcome, so the
#             dose part conditions on it and is read only in the full ventilated
#             cohort; every comparison with a control (the ICU-admission arm, both
#             controls, the difference-in-differences) uses the on/off part, which
#             conditions on nothing
#   arms      ventilated, all patients           (panels A, B and C)
#             ventilated, on IMV at ICU admission (the ventilated side of every comparison
#                                                with a control: status is assigned at ICU
#                                                admission in both arms; the control markers
#                                                and creatinine)
#             no respiratory support             (panel B, the negative control): every
#                                                patient, the divergence read at the
#                                                ventilated cohort's mean severity; follow-up
#                                                ends at escalation to invasive ventilation,
#                                                NIPPV or another advanced support (the
#                                                control's competing event), and no patient
#                                                is in both arms (03 drops from the control
#                                                every patient of the ICU-admission arm)
#             no support, hypoxemic at the index (index SF < 315, the ventilated
#                                                cohort's own gate), read at the ventilated severity:
#                                                the arms then differ in ventilation, not
#                                                hypoxemia (27 writes its DiD separately;
#                                                HYPOXEMIC_CONTROL_MARKERS, the control's markers by default)
#             ventilated, without the lags       (sensitivity: the full ventilated cohort
#                                                without the previous-day SF and pressor
#                                                terms; NOLAG_MARKERS and creatinine)
#             optional: ventilated by baseline SF class (SF_BANDS; off by default)
#   The control is standardised, not matched: severity cannot confound a
#   PFVC fixed by height, age, sex and race, but it could MODIFY the divergence, so the
#   control's divergence varies with its severity anchor and is read at the ventilated
#   mean (20_biotrauma_grid.R, PBWPFVC_JM_SEV_CENTER). No control patient is discarded
#   for severity, and the severity x divergence term tests whether sicker controls
#   diverge faster.
#   Every fit runs on one clock: days from the index (the first qualifying ventilator
#   row; ICU admission in the control), each patient entering the survival submodel at
#   their first trajectory day; the baseline marker is its value in the first 24 h after the index; the dose is
#   VT/PBW at the index and the previous day's VT/PBW minus the index (22_biotrauma_fit.R).
#   each fit adjusted and unadjusted for age, sex and race
#
# Steps, each logged to output/{site}_output/logs/figure4_{stamp}/:
#   1 build     scripts 01-03 for the ventilated cohort and the no-support cohort,
#               each only if its derived tables are missing (FORCE_BUILD=1 rebuilds)
#   2 panels    the 7-day panel of both cohorts
#   3 anchors   the severity-anchor distributions of both cohorts, and the ventilated
#               mean anchor per marker (final/injury/jm_severity_anchor_mean_*)
#   4 centres   each control marker's centre = that ventilated mean
#   5 fits      22_biotrauma_fit.R and 23_biotrauma_report.R for every arm: ventilated
#               (all, then on IMV at ICU admission), any SF classes, the two controls,
#               and the sensitivity without the lags
#   6 figure    supplement/xsec_pfvc_age_control.R (60-day death before escalation
#               against each control, the checks figure's mortality rows),
#               27_control_comparison.R (the difference-in-differences, ventilated at ICU
#               admission against each control), then 24_biotrauma_figures.R: figure 4
#               and the checks figure (biotrauma_fig_checks_*: every outcome against each
#               control, and the rates with and without the lags)
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
# Site default: 46 figure-4 fits (6 markers and creatinine in the ventilated cohort; 3 markers
# and creatinine in the ventilated arm at ICU admission and in each control, the whole one
# and the hypoxemic one; 3 markers and creatinine without the lags; each adjusted and
# unadjusted) plus 1 channel fit, at 5000 / 1000 iterations, four at a time; a run takes
# about 2.5 times as long as one at 2,000 iterations. The survival submodel is deliberately
# small (22_biotrauma_fit.R, HAZARD_SPEC).
# Convergence: the figure's estimates are gated on the lung-size terms (the size level and
# divergence, R-hat <= 1.1, 20_biotrauma_grid.R); the hazard-link convergence is reported
# (hazard_rhat in each manifest) and read with the longitudinal-only comparison
# (jm_lme_check_*, and the checks figure's last panel).
# Knobs (environment): ITER BURNIN CHAINS THIN (5000 / 1000 / 3 / 5), PAR (fits at a
#   time, 4; a 7-day fit at a 7,000-patient site needs about 15-25 GB, so four at once
#   can need 60-100 GB: lower PAR on a smaller machine), MARKERS, CONTROL_MARKERS, CREATININE (1; 0 skips it),
#   SF_BANDS (baseline SF classes of the ventilated cohort, off by default; the lead
#   site runs SF_BANDS="235,315 115,235 0,115"), HYPOXEMIC_CONTROL_MARKERS (the hypoxemic
#   control arm, CONTROL_MARKERS by default, with creatinine as in the other arms; empty
#   skips it), NOLAG_MARKERS (the sensitivity without the lags, "platelets,bilirubin,osi"
#   by default, with creatinine as in the other arms; empty skips it),
#   CHANNEL_MARKERS (step 7, platelets by default; empty skips it),
#   FORCE_BUILD, FORCE_PANEL.
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

MARKERS                   <- knob("MARKERS", "osi,sf,any_pressor,pressor_dose,platelets,bilirubin")   # creatinine runs on its own, with RRT as a third cause
CONTROL_MARKERS           <- knob("CONTROL_MARKERS", "any_pressor,platelets,bilirubin")   # the on/off part, never the dose
HYPOXEMIC_CONTROL_MARKERS <- knob_or_skip("HYPOXEMIC_CONTROL_MARKERS", CONTROL_MARKERS)   # the hypoxemic control arm
NOLAG_MARKERS             <- knob_or_skip("NOLAG_MARKERS", "platelets,bilirubin,osi")     # the sensitivity without the lags
CHANNEL_MARKERS           <- knob_or_skip("CHANNEL_MARKERS", "platelets")             # step 7, the channel breakdown (supplement)
CREATININE  <- knob("CREATININE", "1") == "1"
SF_BANDS    <- strsplit(trimws(knob("SF_BANDS", "")), "[[:space:]]+")[[1]]   # e.g. "235,315 115,235 0,115"; off by default
ITER   <- knob("ITER", "5000"); BURNIN <- knob("BURNIN", "1000"); CHAINS <- knob("CHAINS", "3")
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
  # an empty marker list means creatinine only (an unset one would make 22 fit every marker)
  if (nzchar(markers)) {
    marker_env <- c(arm_env, PBWPFVC_JM_MARKERS = markers)
    run_step(paste0(arm, "_fit"),    "code/22_biotrauma_fit.R",    cohort, marker_env)
    run_step(paste0(arm, "_report"), "code/23_biotrauma_report.R", cohort, marker_env)
  }
  if (with_creatinine) {
    creatinine_env <- c(arm_env, PBWPFVC_JM_MARKERS = "creatinine", PBWPFVC_JM_RRT_EVENT = "1")
    run_step(paste0(arm, "_creatinine_fit"),    "code/22_biotrauma_fit.R",    cohort, creatinine_env)
    run_step(paste0(arm, "_creatinine_report"), "code/23_biotrauma_report.R", cohort, creatinine_env)
  }
}

# ---- 1 build
build_cohort <- function(cohort, derived) {   # cohort, folder holding its derived tables
  # a cohort lacking rrt_sources_available.rds (or, for the control,
  # cohort_icu_stays.parquet) was built by an older 01-03, and its panel cannot be
  # built: it is rebuilt.
  # The control is built after the ventilated cohort and without its ICU-admission arm
  # (03), so a ventilated table newer than the control's rebuilds the control.
  icu_ok <- cohort != "nosupport" || file.exists(file.path(derived, "cohort_icu_stays.parquet"))
  own_table <- file.path(derived, "analysis_cross_sectional.parquet")
  ventilated_table <- file.path(ROOT, "intermediate", "analysis_cross_sectional.parquet")
  control_current <- cohort != "nosupport" || !file.exists(own_table) || !file.exists(ventilated_table) ||
    file.mtime(ventilated_table) <= file.mtime(own_table)
  if (!FORCE_BUILD && icu_ok && control_current && file.exists(own_table) &&
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
# creatinine is always anchored, even with CREATININE=0, so a missing creatinine centre
# skips the whole control arm (step 4); any_pressor's anchor drops the cardiovascular
# SOFA (20_biotrauma_grid.R, ANCHOR_DROP)
ANCHOR_MARKERS <- paste0("creatinine,", CONTROL_MARKERS)
anchor_env <- c(PBWPFVC_JM_ANCHOR_ONLY = "1", PBWPFVC_JM_MARKERS = ANCHOR_MARKERS)
# the controls are read at the mean anchor of the ventilated arm they are compared with,
# the patients on IMV at ICU admission (the difference-in-differences' ventilated side)
ventilated_anchor_env <- c(anchor_env, PBWPFVC_JM_ICU_DAY0 = "1")
run_step("anchors_ventilated", "code/22_biotrauma_fit.R", "imv",       ventilated_anchor_env)
run_step("anchors_nosupport",  "code/22_biotrauma_fit.R", "nosupport", anchor_env)

# ---- 4 centres: the ventilated cohort's mean anchor per control marker, where each
#      control fit reads its divergence (PBWPFVC_JM_SEV_CENTER, 20_biotrauma_grid.R)
CENTER_FILE <- file.path(ROOT, "final", "injury", paste0("jm_severity_anchor_mean_day0_7d_", BASE_SITE, ".csv"))
anchor_markers <- strsplit(ANCHOR_MARKERS, ",")[[1]]
if (DRY) {
  message("[dry] centres: ventilated mean anchor per marker (", ANCHOR_MARKERS, ") from ", CENTER_FILE)
  SEV_CENTER <- "creatinine=M,platelets=M,..."
  centred_markers <- anchor_markers
} else {
  # the file is written by anchors_ventilated; one "marker=mean" pair per control marker.
  # A marker without a centre (too few ICU-day-0 patients with its day-0 value) is left
  # out of the control arms and listed as a failed step; the other markers still run.
  anchor_means <- if (file.exists(CENTER_FILE)) read.csv(CENTER_FILE) else data.frame(marker = character(0), anchor_mean = numeric(0))
  anchor_means <- anchor_means[anchor_means$marker %in% anchor_markers, ]
  centred_markers <- anchor_means$marker
  SEV_CENTER <- if (nrow(anchor_means)) paste0(anchor_means$marker, "=", sprintf("%.4f", anchor_means$anchor_mean), collapse = ",") else ""
  uncentred <- setdiff(anchor_markers, centred_markers)
  if (length(uncentred)) {
    message("severity centres missing for ", paste(uncentred, collapse = ", "), " in ", CENTER_FILE,
            "; those markers are left out of the control arms")
    FAILED <- c(FAILED, paste0("centre_", uncentred))
  }
  if (nzchar(SEV_CENTER)) message("[", timestamp(), "] severity centres (ICU-day-0 ventilated mean anchor): ", SEV_CENTER)
}
# the markers each control arm can be read at: those with a centre
centred_list <- function(markers) paste(intersect(strsplit(markers, ",")[[1]], centred_markers), collapse = ",")
CONTROL_MARKERS_CENTRED   <- centred_list(CONTROL_MARKERS)
HYPOXEMIC_MARKERS_CENTRED <- centred_list(HYPOXEMIC_CONTROL_MARKERS)
CONTROL_CREATININE        <- CREATININE && "creatinine" %in% centred_markers

# ---- 5 fits, arm by arm
fit_arm("ventilated", "imv", MARKERS)
# the ventilated side of the comparisons with a control: patients on IMV at ICU
# admission, fitted for the control's markers (and creatinine)
fit_arm("ventilated_day0", "imv", CONTROL_MARKERS, c(PBWPFVC_JM_ICU_DAY0 = "1"))
for (band in SF_BANDS)
  fit_arm(paste0("ventilated_sf", sub(",", "to", band)), "imv", MARKERS, c(PBWPFVC_JM_SF_BAND = band))
if (nzchar(SEV_CENTER)) {
  if (nzchar(CONTROL_MARKERS_CENTRED) || CONTROL_CREATININE)
    fit_arm("nosupport", "nosupport", CONTROL_MARKERS_CENTRED, c(PBWPFVC_JM_SEV_CENTER = SEV_CENTER),
            with_creatinine = CONTROL_CREATININE)
  # the hypoxemic control: the same control, index SF < 315, with the same markers
  if (nzchar(HYPOXEMIC_CONTROL_MARKERS) && (nzchar(HYPOXEMIC_MARKERS_CENTRED) || CONTROL_CREATININE))
    fit_arm("nosupport_hypoxemic", "nosupport", HYPOXEMIC_MARKERS_CENTRED,
            c(PBWPFVC_JM_SEV_CENTER = SEV_CENTER, PBWPFVC_JM_SF_BAND = "0,315"),
            with_creatinine = CONTROL_CREATININE)
}

# the sensitivity without the previous-day SF and pressor terms, full ventilated cohort
if (nzchar(NOLAG_MARKERS))
  fit_arm("ventilated_nolag", "imv", NOLAG_MARKERS, c(PBWPFVC_JM_NO_LAGS = "1"))

# ---- 6 comparison table and the figure
# figure rows are fixed; changing MARKERS does not change them (a marker with no fit
# on disk is left out of the figure)
FIG_MARKERS <- paste0("platelets,bilirubin", if (CREATININE) ",creatinine", ",any_pressor,pressor_dose,osi,sf")
# death against each control: the PFVC association with 60-day death before escalation
# in the ventilated arm at ICU admission and in the no-support control, everyone and
# hypoxemic at the index (the checks figure's mortality rows); it reads both cohorts'
# tables and the control's 7-day panel
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
