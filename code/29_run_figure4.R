# =============================================================================
# 29_run_figure4: every analysis behind manuscript figure 4, then the figure
# =============================================================================
# Figure 4 asks whether a smaller predicted lung, at the same VT/PBW, goes with
# organ-injury markers that worsen faster over the first 7 days of ventilation. This
# script produces it for one site, from the CLIF tables to the PDF:
#
#   markers   platelets, creatinine (dialysis, continuous or intermittent,
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
#   The control is standardised, not matched: severity cannot confound a
#   PFVC fixed by height, age, sex and race, but it could MODIFY the divergence, so the
#   control's divergence varies with its severity anchor and is read at the ventilated
#   mean (20_biotrauma_grid.R, PBWPFVC_JM_SEV_CENTER). No control patient is discarded
#   for severity, and the severity x divergence term tests whether sicker controls
#   diverge faster.
#   Clocks (PBWPFVC_JM_CLOCK, 20_biotrauma_grid.R). Each fit runs on one clock, days
#   from its time zero t0, and the comparisons with a control share one:
#     index  t0 = the index, the first qualifying ventilator row: ventilated all
#            (panel A), the arm without the lags, the channels
#     icu    t0 = ICU admission, the first ICU in_dttm of the stay: ventilated on IMV
#            at ICU admission, its anchors, both no-support arms and their anchors,
#            so both sides of every difference-in-differences start at ICU admission
#   On either clock each patient enters the survival submodel at their first
#   trajectory day, follow-up ends 168 h after t0, and the baseline marker is its value
#   in the first 24 h after t0; the dose is VT/PBW at the index (the first 24 h after
#   index_dttm, on both clocks) and the previous day's VT/PBW minus it (22_biotrauma_fit.R).
#   each fit adjusted and unadjusted for age, sex and race
#
# Steps, each logged to output/{site}_output/logs/figure4_{stamp}/:
#   1 build     scripts 01-03 for the ventilated cohort and the no-support cohort,
#               each only if its derived tables are missing (FORCE_BUILD=1 rebuilds);
#               for the ventilated cohort 03 runs twice, the second time for the
#               ungated cohort (PBWPFVC_VTPBW_GATE=0, below)
#   2 panels    three 7-day panels: the ventilated cohort on the index clock, the
#               ventilated patients on IMV at ICU admission on the icu clock, and the
#               no-support cohort on the icu clock
#   2b contrast the ventilated-vs-control mortality contrasts, which need both cohorts and
#               the control's panel but no fit: supplement/xsec_strain_invariance.R,
#               xsec_mortality_channel_equality.R and xsec_pfvc_age_control.R (60-day
#               death with controls censored at any escalation against each control, also the checks figure's
#               mortality rows). Each runs on the gated and the ungated ventilated cohort
#               (final/supplement/, final/ungated/supplement/), and the last two each
#               with and without VT/PBW in the ventilated arm's models (novtpbw/
#               subfolders; xsec_strain_invariance.R fits both dose versions itself)
#   3 anchors   the severity-anchor distributions of both cohorts, and the ventilated
#               mean anchor per marker (final/injury/jm_severity_anchor_mean_*), both
#               on the icu clock
#   3b overlap  supplement/xsec_intubation_overlap.R: a propensity score for intubation
#               (ventilated at ICU admission against each control, covariates from the
#               24 h before ICU admission) and its overlap, balance and effective sample
#               size, the diagnostic for weighting the controls instead of standardising
#               them; it needs both cohorts built and nothing from the fits
#   4 centres   each control marker's centre = that ventilated mean
#   5 fits      22_biotrauma_fit.R and 23_biotrauma_report.R for every arm: ventilated
#               (all, then on IMV at ICU admission), the two controls,
#               and the sensitivity without the lags
#   5b measure  26_measurement_check.R on both sides of the difference-in-differences:
#               why bilirubin is not a figure-4 marker (its missingness is informative:
#               who enters its analysis, by severity and death, against the routine
#               labs), whether being measured depends on predicted lung size, and
#               whether the platelet divergence survives on bilirubin's sampling
#               frame (panels only, seconds; MEASUREMENT_CHECK=0 skips it)
#   6 figure    27_control_comparison.R (the difference-in-differences, ventilated at ICU
#               admission against each control), then 24_biotrauma_figures.R: figure 4
#               and the checks figure (biotrauma_fig_checks_*: every outcome against each
#               control, and the rates with and without the lags)
#   7 channels  the supplement's channel breakdown: the ventilated divergence read through
#               each GLI piece of log PFVC (height, age, sex, race; 22's channels form, one
#               fit per marker, the pieces standing in for the demographics), and its
#               figure (biotrauma_fig_channels_*). CHANNEL_MARKERS, platelets by default;
#               empty skips it
#   8 ungated   the same divergences and difference-in-differences on the ungated
#               ventilated cohort, from the longitudinal submodels alone
#               (PBWPFVC_JM_LME_ONLY=1; the joint model and its longitudinal submodel
#               agree on the gated cohort, jm_lme_check_*): the ungated panels, the
#               ungated ventilated mean anchors, LME fits of the ventilated cohort (all,
#               and on IMV at ICU admission) and of both controls read at those anchors,
#               then 27. Ventilated tables in final/ungated/injury/; the controls' LME
#               tables (tagged lmeonly_) in final/controls/, read at the ungated anchors.
#               Minutes, not hours. UNGATED_MARKERS (CONTROL_MARKERS by default, with
#               creatinine as in the other arms); empty skips it
#   8b ungated forms  step 8's difference-in-differences again for each of UNGATED_FORMS
#               (log PBW/PFVC, log PFVC at age 25, log PBW/PFVC at age 25 by default;
#               empty skips it): does the age-free part of the strain error move the
#               markers in proportion to the whole? Tables tagged by form
#               (jm_control_did_{form}_lmeonly_7d_*, final/ungated/injury/)
#   8c reference volumes  step 8's difference-in-differences for platelets and creatinine,
#               on age <= 80 in each race stratum (REFVOL_STRATA, White first), for the
#               strain error against PFVC, PFVC at age 25, TLC, alveolar volume and FRC
#               and for log PFVC (REFVOL_FORMS); tables jm_refvol_control_did_*
#   8d lambda sweep  the same for the strain error against log V(lambda) = log PFVC at
#               25 + lambda (log PFVC - log PFVC at 25), lambda in LAMBDA_GRID, on the
#               whole ungated cohort; tables jm_control_did_lambda<NNN>_lmeonly_7d_*
#   8e flexible severity  supplement/xsec_dml_contrast.R on the gated and the ungated
#               cohort: the ventilated-minus-control contrast for mortality and the
#               platelet and creatinine slopes with severity entered by double machine
#               learning and with overlap weights (dml_contrast_*, dml_overlap_*). DML (1;
#               0 skips it)
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
# Site default: 36 figure-4 fits (5 markers and creatinine in the ventilated cohort; 2 markers
# and creatinine in the ventilated arm at ICU admission and in each control, the whole one
# and the hypoxemic one; 2 markers and creatinine without the lags; each adjusted and
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
#   MEASUREMENT_CHECK (1; 0 skips step 5b),
#   HYPOXEMIC_CONTROL_MARKERS (the hypoxemic
#   control arm, CONTROL_MARKERS by default, with creatinine as in the other arms; empty
#   skips it), NOLAG_MARKERS (the sensitivity without the lags, "platelets,osi"
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

MARKERS                   <- knob("MARKERS", "osi,sf,any_pressor,pressor_dose,platelets")   # creatinine runs on its own, with RRT as a third cause
CONTROL_MARKERS           <- knob("CONTROL_MARKERS", "any_pressor,platelets")   # the on/off part, never the dose
HYPOXEMIC_CONTROL_MARKERS <- knob_or_skip("HYPOXEMIC_CONTROL_MARKERS", CONTROL_MARKERS)   # the hypoxemic control arm
NOLAG_MARKERS             <- knob_or_skip("NOLAG_MARKERS", "platelets,osi")     # the sensitivity without the lags
CHANNEL_MARKERS           <- knob_or_skip("CHANNEL_MARKERS", "platelets")             # step 7, the channel breakdown (supplement)
UNGATED_MARKERS           <- knob_or_skip("UNGATED_MARKERS", CONTROL_MARKERS)         # step 8, the ungated cohort's LME fits
UNGATED_FORMS             <- knob_or_skip("UNGATED_FORMS", "disc_level,pfvc25,disc25_level")   # step 8b, the other size forms
REFVOL_FORMS              <- knob_or_skip("REFVOL_FORMS", "pfvc,disc_level,disc25_level,disc_tlc,disc_va,disc_frc")   # step 8c
REFVOL_STRATA             <- knob_or_skip("REFVOL_STRATA", "white,notwhite")   # step 8c, the race strata (age <= 80)
REFVOL_MARKERS            <- knob("REFVOL_MARKERS", "platelets")   # steps 8c and 8d, with creatinine as in the other arms
LAMBDA_GRID               <- knob_or_skip("LAMBDA_GRID", "0,0.1,0.2,0.3,0.5,0.75,1")   # step 8d, the lambda sweep
CREATININE  <- knob("CREATININE", "1") == "1"
MEASUREMENT_CHECK <- knob("MEASUREMENT_CHECK", "1") == "1"   # step 5b, 26_measurement_check.R
DML <- knob("DML", "1") == "1"   # step 8e, supplement/xsec_dml_contrast.R
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
           PBWPFVC_JM_REPLACE_TABLES = "1",   # each arm writes its tables whole: no stale markers
           PBWPFVC_JM_ITER = ITER, PBWPFVC_JM_BURNIN = BURNIN, PBWPFVC_JM_CHAINS = CHAINS,
           PBWPFVC_JM_THIN = THIN, PBWPFVC_JM_PAR = PAR)
message("site ", BASE_SITE, "; markers ", MARKERS, if (CREATININE) ",creatinine", "; controls ", CONTROL_MARKERS,
        "; chains ",
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

# the settings of the ICU-admission clock, for every arm of a comparison with a control
ICU_CLOCK <- c(PBWPFVC_JM_CLOCK = "icu")
# the ungated ventilated cohort (utils/config.R): its own panels, final/ungated/ tables
UNGATED <- c(PBWPFVC_VTPBW_GATE = "0")
# the panel step an arm reads: its cohort, its clock and, for the ventilated cohort, its gate
panel_step_for <- function(cohort, arm_env) {
  clock <- if ("PBWPFVC_JM_CLOCK" %in% names(arm_env)) arm_env[["PBWPFVC_JM_CLOCK"]] else "index"
  ungated <- identical(unname(arm_env["PBWPFVC_VTPBW_GATE"]), "0")
  paste0("panel_", cohort, if (clock == "icu") "_icu" else "", if (ungated) "_ungated" else "")
}

# one arm: the markers, then creatinine with dialysis as a third competing cause. 23
# reports joint models; an LME-only arm (PBWPFVC_JM_LME_ONLY=1) has its tables from 22.
fit_arm <- function(arm, cohort, markers, arm_env = character(0), with_creatinine = CREATININE) {
  with_report <- !identical(unname(arm_env["PBWPFVC_JM_LME_ONLY"]), "1")
  # an arm whose panel failed to build is not fitted: its old panel is out of date
  panel_step <- panel_step_for(cohort, arm_env)
  if (panel_step %in% FAILED) {
    message("[", timestamp(), "] ", arm, ": skipped, its panel (", panel_step, ") failed to build")
    FAILED <<- c(FAILED, paste0(arm, "_skipped"))
    return(invisible())
  }
  # an empty marker list means creatinine only (an unset one would make 22 fit every marker)
  if (nzchar(markers)) {
    marker_env <- c(arm_env, PBWPFVC_JM_MARKERS = markers)
    run_step(paste0(arm, "_fit"),    "code/22_biotrauma_fit.R",    cohort, marker_env)
    if (with_report) run_step(paste0(arm, "_report"), "code/23_biotrauma_report.R", cohort, marker_env)
  }
  if (with_creatinine) {
    creatinine_env <- c(arm_env, PBWPFVC_JM_MARKERS = "creatinine", PBWPFVC_JM_RRT_EVENT = "1")
    run_step(paste0(arm, "_creatinine_fit"),    "code/22_biotrauma_fit.R",    cohort, creatinine_env)
    if (with_report) run_step(paste0(arm, "_creatinine_report"), "code/23_biotrauma_report.R", cohort, creatinine_env)
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
  # the ventilated cohort's ungated table comes from a second run of 03, after the gated one
  ungated_table <- file.path(derived, "analysis_cross_sectional_ungated.parquet")
  ungated_current <- cohort != "imv" || (file.exists(ungated_table) && file.exists(own_table) &&
                                           file.mtime(own_table) <= file.mtime(ungated_table))
  built <- !FORCE_BUILD && icu_ok && control_current && file.exists(own_table) &&
    file.exists(file.path(derived, "rrt_sources_available.rds"))
  if (built) {
    message("[", timestamp(), "] ", cohort, ": cohort already built, scripts 01-03 skipped (FORCE_BUILD=1 rebuilds)")
  } else {
    for (script in c("01_cohort_identification", "02_quality_checks", "03_variable_derivation"))
      run_step(paste0("build_", cohort, "_", script), file.path("code", paste0(script, ".R")), cohort)
  }
  if (cohort == "imv" && !(built && ungated_current))
    run_step("build_imv_03_variable_derivation_ungated", "code/03_variable_derivation.R", cohort, UNGATED)
}
build_cohort("imv",       file.path(ROOT, "intermediate"))
build_cohort("nosupport", file.path(ROOT, "intermediate", "controls", "nosupport"))

# ---- 2 panels, rebuilt only when something they are built from has changed: a fit made
#      on an older panel is refitted (22_biotrauma_fit.R compares the times), so an
#      unconditional rebuild would refit everything on every rerun. FORCE_PANEL=1 rebuilds.
#      Each clock has its own panel files (an icu-clock panel's end in "_icu"), and so
#      does the ungated cohort (then "_ungated").
build_panel <- function(cohort, derived, panel_env = character(0)) {   # cohort, folder holding its derived tables, clock, gate
  step_name <- panel_step_for(cohort, panel_env)
  panel <- file.path(derived, paste0("jm_surv_7d", sub(paste0("^panel_", cohort), "", step_name), ".parquet"))   # "_icu", "_ungated"
  cross_sectional_table <- if (grepl("_ungated$", step_name)) "analysis_cross_sectional_ungated.parquet" else "analysis_cross_sectional.parquet"
  if (!FORCE_PANEL && !DRY && file.exists(panel)) {
    inputs <- c("code/21_biotrauma_panel.R", "code/10_panel_common.R", "code/20_biotrauma_grid.R", "utils/config.R",
                file.path(derived, c(cross_sectional_table, "cohort_dialysis.parquet", "cohort_esrd.parquet")))
    stale <- inputs[!file.exists(inputs) | file.mtime(inputs) > file.mtime(panel)]
    if (!length(stale)) {
      message("[", timestamp(), "] ", step_name, ": up to date, kept (FORCE_PANEL=1 rebuilds)")
      return(invisible())
    }
    message("[", timestamp(), "] ", step_name, ": ", basename(stale[1]), " is newer than the panel (or missing); rebuilding")
  }
  run_step(step_name, "code/21_biotrauma_panel.R", cohort, panel_env)
}
build_panel("imv",       file.path(ROOT, "intermediate"))                              # index clock: ventilated all, no lags, channels
build_panel("imv",       file.path(ROOT, "intermediate"), ICU_CLOCK)                   # icu clock: ventilated on IMV at ICU admission
build_panel("nosupport", file.path(ROOT, "intermediate", "controls", "nosupport"), ICU_CLOCK)   # icu clock: the control

# ---- 2b the ventilated-vs-control mortality contrasts: both cohorts and the control's
#      icu-clock panel, no fit. Gated and ungated ventilated cohort; channel equality and
#      the age control also with and without VT/PBW in the ventilated arm (novtpbw/).
#      The age control's gated run with VT/PBW is the checks figure's mortality rows (step 6).
for (gate in c("1", "0")) {
  gate_env <- if (gate == "0") UNGATED else character(0)
  gate_tag <- if (gate == "0") "_ungated" else ""
  run_step(paste0("strain_invariance", gate_tag), "code/supplement/xsec_strain_invariance.R", "imv", gate_env)
  for (with_vtpbw in c("1", "0")) {
    step_env <- c(gate_env, PBWPFVC_DID_VTPBW = with_vtpbw)
    dose_tag <- if (with_vtpbw == "0") "_novtpbw" else ""
    run_step(paste0("channel_equality", gate_tag, dose_tag), "code/supplement/xsec_mortality_channel_equality.R", "imv", step_env)
    run_step(paste0("mortality_contrast", gate_tag, dose_tag), "code/supplement/xsec_pfvc_age_control.R", "imv", step_env)
  }
}

# ---- 3 anchors
# creatinine is always anchored, even with CREATININE=0, so a missing creatinine centre
# skips the whole control arm (step 4); any_pressor's anchor drops the cardiovascular
# SOFA (20_biotrauma_grid.R, ANCHOR_DROP)
ANCHOR_MARKERS <- paste0("creatinine,", CONTROL_MARKERS)
anchor_env <- c(PBWPFVC_JM_ANCHOR_ONLY = "1", PBWPFVC_JM_MARKERS = ANCHOR_MARKERS)
# the controls are read at the mean anchor of the ventilated arm they are compared with,
# the patients on IMV at ICU admission (the difference-in-differences' ventilated side),
# from the icu-clock panels both arms are fitted on
ventilated_anchor_env <- c(anchor_env, ICU_CLOCK, PBWPFVC_JM_ICU_DAY0 = "1")
run_step("anchors_ventilated", "code/22_biotrauma_fit.R", "imv",       ventilated_anchor_env)
run_step("anchors_nosupport",  "code/22_biotrauma_fit.R", "nosupport", c(anchor_env, ICU_CLOCK))

# ---- 3b overlap: could a propensity score for intubation weight the controls instead?
run_step("intubation_overlap", "code/supplement/xsec_intubation_overlap.R")

# ---- 4 centres: the ventilated cohort's mean anchor per control marker, where each
#      control fit reads its divergence (PBWPFVC_JM_SEV_CENTER, 20_biotrauma_grid.R)
anchor_markers <- strsplit(ANCHOR_MARKERS, ",")[[1]]
# the centres in one anchor file: list(sev_center = "marker=mean,...", markers = those centred)
read_centres <- function(center_file, label = "") {
  if (DRY) {
    message("[dry] centres", label, ": ventilated mean anchor per marker (", ANCHOR_MARKERS, ") from ", center_file)
    return(list(sev_center = "creatinine=M,platelets=M,...", markers = anchor_markers))
  }
  # the file is written by an anchors step; one "marker=mean" pair per control marker.
  # A marker without a centre (too few ICU-day-0 patients with its day-0 value) is left
  # out of the control arms and listed as a failed step; the other markers still run.
  anchor_means <- if (file.exists(center_file)) read.csv(center_file) else data.frame(marker = character(0), anchor_mean = numeric(0))
  anchor_means <- anchor_means[anchor_means$marker %in% anchor_markers, ]
  sev_center <- if (nrow(anchor_means)) paste0(anchor_means$marker, "=", sprintf("%.4f", anchor_means$anchor_mean), collapse = ",") else ""
  uncentred <- setdiff(anchor_markers, anchor_means$marker)
  if (length(uncentred)) {
    message("severity centres missing for ", paste(uncentred, collapse = ", "), " in ", center_file,
            "; those markers are left out of the control arms")
    FAILED <<- c(FAILED, paste0("centre", label, "_", uncentred))
  }
  if (nzchar(sev_center)) message("[", timestamp(), "] severity centres", label, " (ICU-day-0 ventilated mean anchor): ", sev_center)
  list(sev_center = sev_center, markers = anchor_means$marker)
}
centres <- read_centres(file.path(ROOT, "final", "injury", paste0("jm_severity_anchor_mean_day0_7d_", BASE_SITE, ".csv")))
SEV_CENTER <- centres$sev_center
# the markers each control arm can be read at: those with a centre
centred_list <- function(markers, centred) paste(intersect(strsplit(markers, ",")[[1]], centred), collapse = ",")
CONTROL_MARKERS_CENTRED   <- centred_list(CONTROL_MARKERS, centres$markers)
HYPOXEMIC_MARKERS_CENTRED <- centred_list(HYPOXEMIC_CONTROL_MARKERS, centres$markers)
CONTROL_CREATININE        <- CREATININE && "creatinine" %in% centres$markers

# ---- 5 fits, arm by arm (index clock unless the arm sets ICU_CLOCK)
fit_arm("ventilated", "imv", MARKERS)
# the ventilated side of the comparisons with a control: patients on IMV at ICU
# admission, fitted for the control's markers (and creatinine), on the icu clock
fit_arm("ventilated_day0", "imv", CONTROL_MARKERS, c(ICU_CLOCK, PBWPFVC_JM_ICU_DAY0 = "1"))
if (nzchar(SEV_CENTER)) {
  if (nzchar(CONTROL_MARKERS_CENTRED) || CONTROL_CREATININE)
    fit_arm("nosupport", "nosupport", CONTROL_MARKERS_CENTRED, c(ICU_CLOCK, PBWPFVC_JM_SEV_CENTER = SEV_CENTER),
            with_creatinine = CONTROL_CREATININE)
  # the hypoxemic control: the same control, index SF < 315, with the same markers
  if (nzchar(HYPOXEMIC_CONTROL_MARKERS) && (nzchar(HYPOXEMIC_MARKERS_CENTRED) || CONTROL_CREATININE))
    fit_arm("nosupport_hypoxemic", "nosupport", HYPOXEMIC_MARKERS_CENTRED,
            c(ICU_CLOCK, PBWPFVC_JM_SEV_CENTER = SEV_CENTER, PBWPFVC_JM_SF_BAND = "0,315"),
            with_creatinine = CONTROL_CREATININE)
}

# the sensitivity without the previous-day SF and pressor terms, full ventilated cohort
if (nzchar(NOLAG_MARKERS))
  fit_arm("ventilated_nolag", "imv", NOLAG_MARKERS, c(PBWPFVC_JM_NO_LAGS = "1"))

# ---- 5b is a marker missing in a way that biases the comparison? (26). Bilirubin is
#      ordered, not routine, and it reaches a smaller and sicker share of the control
#      arm than of the ventilated one. Panels only, so this costs seconds and does not
#      wait on any fit; both sides of the difference-in-differences, on their own clock.
if (MEASUREMENT_CHECK) {
  run_step("measurement_ventilated", "code/26_measurement_check.R", "imv",       ICU_CLOCK)
  run_step("measurement_nosupport",  "code/26_measurement_check.R", "nosupport", ICU_CLOCK)
  run_step("measurement_ventilated_index", "code/26_measurement_check.R", "imv", character(0))
}

# ---- 6 comparison table and the figure
# figure rows are fixed; changing MARKERS does not change them (a marker with no fit
# on disk is left out of the figure)
FIG_MARKERS <- paste0("platelets", if (CREATININE) ",creatinine", ",any_pressor,pressor_dose,osi,sf")
# the checks figure's mortality rows (death against each control) come from step 2b's
# gated mortality_contrast
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

# ---- 8 the ungated ventilated cohort, longitudinal submodels alone. The same arms as
#      the difference-in-differences above: its own panels and anchors, the ventilated
#      fits (all on the index clock; on IMV at ICU admission on the icu clock), both
#      controls read at the ungated ventilated mean anchors, then 27. The control fits are
#      the gated run's controls refitted by maximum likelihood; only their centre differs.
if (nzchar(UNGATED_MARKERS)) {
  LME <- c(PBWPFVC_JM_LME_ONLY = "1")
  build_panel("imv", file.path(ROOT, "intermediate"), UNGATED)
  build_panel("imv", file.path(ROOT, "intermediate"), c(ICU_CLOCK, UNGATED))
  run_step("anchors_ventilated_ungated", "code/22_biotrauma_fit.R", "imv", c(ventilated_anchor_env, UNGATED))
  ungated_centres <- read_centres(file.path(ROOT, "final", "ungated", "injury",
                                            paste0("jm_severity_anchor_mean_day0_7d_", BASE_SITE, ".csv")), "_ungated")
  ungated_creatinine <- CREATININE && "creatinine" %in% ungated_centres$markers
  fit_arm("ventilated_ungated", "imv", UNGATED_MARKERS, c(UNGATED, LME))
  fit_arm("ventilated_day0_ungated", "imv", UNGATED_MARKERS, c(ICU_CLOCK, PBWPFVC_JM_ICU_DAY0 = "1", UNGATED, LME))
  if (nzchar(ungated_centres$sev_center)) {
    centred_ungated <- centred_list(UNGATED_MARKERS, ungated_centres$markers)
    if (nzchar(centred_ungated) || ungated_creatinine) {
      fit_arm("nosupport_lme", "nosupport", centred_ungated,
              c(ICU_CLOCK, PBWPFVC_JM_SEV_CENTER = ungated_centres$sev_center, LME), with_creatinine = ungated_creatinine)
      if (nzchar(HYPOXEMIC_CONTROL_MARKERS))
        fit_arm("nosupport_hypoxemic_lme", "nosupport", centred_list(HYPOXEMIC_CONTROL_MARKERS, strsplit(centred_ungated, ",")[[1]]),
                c(ICU_CLOCK, PBWPFVC_JM_SEV_CENTER = ungated_centres$sev_center, PBWPFVC_JM_SF_BAND = "0,315", LME),
                with_creatinine = ungated_creatinine)
    }
  }
  run_step("comparison_ungated", "code/27_control_comparison.R", "imv", c(UNGATED, LME))

  # ---- 8b the same difference-in-differences for the strain error and the age-free
  #      forms (UNGATED_FORMS): log PBW/PFVC, and log PFVC and log PBW/PFVC with GLI at
  #      age 25. Strain predicts that the age-free part of the strain error moves the
  #      markers in proportion to the whole; these fits are what that comparison reads.
  #      Same panels, anchors and arms as step 8; only the size term changes.
  if (nzchar(ungated_centres$sev_center) && nzchar(UNGATED_FORMS)) {
    for (form in strsplit(UNGATED_FORMS, ",")[[1]]) {
      form_env <- c(UNGATED, LME, PBWPFVC_JM_MODIFIER = form)
      fit_arm(paste0("ventilated_day0_ungated_", form), "imv", UNGATED_MARKERS,
              c(ICU_CLOCK, PBWPFVC_JM_ICU_DAY0 = "1", form_env))
      if (nzchar(centred_ungated) || ungated_creatinine) {
        fit_arm(paste0("nosupport_lme_", form), "nosupport", centred_ungated,
                c(ICU_CLOCK, PBWPFVC_JM_SEV_CENTER = ungated_centres$sev_center, LME, PBWPFVC_JM_MODIFIER = form),
                with_creatinine = ungated_creatinine)
        if (nzchar(HYPOXEMIC_CONTROL_MARKERS))
          fit_arm(paste0("nosupport_hypoxemic_lme_", form), "nosupport",
                  centred_list(HYPOXEMIC_CONTROL_MARKERS, strsplit(centred_ungated, ",")[[1]]),
                  c(ICU_CLOCK, PBWPFVC_JM_SEV_CENTER = ungated_centres$sev_center, PBWPFVC_JM_SF_BAND = "0,315",
                    LME, PBWPFVC_JM_MODIFIER = form),
                  with_creatinine = ungated_creatinine)
      }
      run_step(paste0("comparison_ungated_", form), "code/27_control_comparison.R", "imv", form_env)
    }
  }

  # ---- 8c reference lung volumes: the same difference-in-differences, platelets and
  #      creatinine (REFVOL_MARKERS), for the strain error against each candidate
  #      denominator (REFVOL_FORMS: PFVC, PFVC at age 25, TLC, alveolar volume, FRC, and
  #      log PFVC itself) on one sample per stratum (REFVOL_STRATA): age <= 80, White,
  #      and everyone else, pre-specified. The control is read at the ungated centres.
  #      Tables jm_refvol_control_did_{stratum}80_{form}_lmeonly_7d_*.
  if (nzchar(ungated_centres$sev_center) && nzchar(REFVOL_FORMS) && nzchar(REFVOL_STRATA)) {
    refvol_markers <- centred_list(REFVOL_MARKERS, ungated_centres$markers)
    for (stratum in strsplit(REFVOL_STRATA, ",")[[1]]) for (form in strsplit(REFVOL_FORMS, ",")[[1]]) {
      # the gate applies to the ventilated cohort only; the control has none
      form_env <- c(LME, PBWPFVC_JM_MODIFIER = form, PBWPFVC_JM_REF_STRATUM = stratum)
      step_tag <- paste0(stratum, "_", form)
      fit_arm(paste0("ventilated_day0_refvol_", step_tag), "imv", refvol_markers,
              c(ICU_CLOCK, PBWPFVC_JM_ICU_DAY0 = "1", UNGATED, form_env), with_creatinine = ungated_creatinine)
      fit_arm(paste0("nosupport_lme_refvol_", step_tag), "nosupport", refvol_markers,
              c(ICU_CLOCK, PBWPFVC_JM_SEV_CENTER = ungated_centres$sev_center, form_env), with_creatinine = ungated_creatinine)
      run_step(paste0("comparison_refvol_", step_tag), "code/27_control_comparison.R", "imv", c(UNGATED, form_env))
    }
  }

  # ---- 8d the lambda sweep: the same difference-in-differences, platelets and
  #      creatinine, for the strain error against a denominator keeping the share
  #      lambda of GLI's age decline in FVC (LAMBDA_GRID; 22's lambda<NNN> forms), on the
  #      whole ungated cohort. 27 runs straight after each lambda's fits, because the
  #      scale table it reads holds the last lambda fitted.
  if (nzchar(ungated_centres$sev_center) && nzchar(LAMBDA_GRID)) {
    lambda_markers <- centred_list(REFVOL_MARKERS, ungated_centres$markers)
    for (lambda in as.numeric(strsplit(LAMBDA_GRID, ",")[[1]])) {
      form <- sprintf("lambda%03d", round(100 * lambda))
      form_env <- c(LME, PBWPFVC_JM_MODIFIER = form)
      fit_arm(paste0("ventilated_day0_", form), "imv", lambda_markers,
              c(ICU_CLOCK, PBWPFVC_JM_ICU_DAY0 = "1", UNGATED, form_env), with_creatinine = ungated_creatinine)
      fit_arm(paste0("nosupport_lme_", form), "nosupport", lambda_markers,
              c(ICU_CLOCK, PBWPFVC_JM_SEV_CENTER = ungated_centres$sev_center, form_env), with_creatinine = ungated_creatinine)
      run_step(paste0("comparison_", form), "code/27_control_comparison.R", "imv", c(UNGATED, form_env))
    }
  }
}

# ---- 8e flexible severity: the mortality and marker-slope contrasts with severity
#      entered by double machine learning and by overlap weights. Both cohorts' icu-clock
#      panels are needed; build_panel() keeps those already up to date.
if (DML) {
  build_panel("imv", file.path(ROOT, "intermediate"), c(ICU_CLOCK, UNGATED))
  run_step("dml_contrast", "code/supplement/xsec_dml_contrast.R", "imv", character(0))
  run_step("dml_contrast_ungated", "code/supplement/xsec_dml_contrast.R", "imv", UNGATED)
}

if (!DRY) {
  message("\nfigure  -> ", file.path(ROOT, "final", "injury", paste0("biotrauma_fig_main_pfvc_7d_", BASE_SITE, ".pdf")))
  message("tables  -> ", file.path(ROOT, "final", "injury"), "/, ", file.path(ROOT, "final", "controls"), "/, ",
          file.path(ROOT, "final", "supplement"), "/ and their ungated twins in ", file.path(ROOT, "final", "ungated"), "/")
  if (length(FAILED)) {
    message("FAILED steps (", length(FAILED), "): ", paste(FAILED, collapse = " "))
    quit(status = 1)
  }
  message("all steps ok")
}
