# =============================================================================
# Script 13 (fit): Biotrauma joint models -- one shared-random-effects JM per marker
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# For each organ-injury marker, fits a joint model (JMbayes2) that links
#   longitudinal submodel   log marker on day t  ~  spline(day) + previous-day
#                           VT/PFVC as a WITHIN-patient deviation from the
#                           patient's mean + that mean (between) + baseline
#                           marker + previous-day confounders + non-respiratory
#                           SOFA, BMI and demographics; random intercept and
#                           slope per patient
#   survival submodel       cause-specific stratified Cox, death vs extubation,
#                           with index VT/PBW (dose), log PFVC (size) and
#                           baseline covariates: the paper's primary exposure set
#   association             current value and current slope of the marker on
#                           each cause-specific hazard
#
# Three questions (docs/joint_model_plan_2026-09.md, section 3):
#   Q1  the previous-day strain coefficient in the longitudinal submodel: a
#       conditional, within-patient dose-response (not a policy effect; the 11.*
#       g-methods are the causal version)
#   Q2  the value and slope association parameters, per cause
#   Markers: creatinine, platelets, bilirubin, sf, dp, ne_equiv_peak (log dose per kg,
#   flagged: it carries -2 log(height)), any_pressor (the hurdle's binary part:
#   a logistic mixed model of any vasoactive running, the vasopressor read).
#   Q3  the log PFVC (and VT/PBW) coefficient on the death hazard in the plain
#       cause-specific Cox (no linkage) versus inside the JM (with linkage):
#       "association absorbed by the trajectory", not proportion mediated
#
# Model set per marker:
#   main      full cohort, adjusted (age spline, sex, race) and unadjusted
#   hetero    plateau-measured subset, adds baseline specific elastance
#             (Ers x PFVC) x previous-day strain; adjusted only
# Each marker's model excludes its OWN lag (the SF model drops lagged SF, the
# NE-equivalent model drops the lagged pressor flag). Random-effects structure is
# fixed per marker in advance; a fit that fails is reported as failed, never
# refit with a smaller model.
#
# Inputs:  intermediate/jm_long_{H}d.parquet, jm_surv_{H}d.parquet (13_biotrauma_panel.R)
# Outputs: final/jm_estimates_{H}d_{site}.csv   every coefficient of every fit (poolable)
#          final/jm_absorption_{H}d_{site}.csv  Q3 table
#          final/jm_manifest_{H}d_{site}.csv    fit status, counts, R-hat gate
#          intermediate/jm_fit_{marker}_{model}_{adj}_{H}d.rds   fit bundles
#          (patient-level rows inside, so never in final/)
#
# Environment knobs: PBWPFVC_JM_GRID (6h | daily) with PBWPFVC_JM_HORIZON_H (48) or
#   PBWPFVC_JM_HORIZON (7 days), PBWPFVC_JM_MARKERS (comma list),
#   PBWPFVC_JM_MODELS (main,hetero), PBWPFVC_JM_BASELINE (free | offset; offset
#   fixes the baseline coefficient at 1 = the log percent-change outcome, written
#   with an offset_ prefix), PBWPFVC_JM_ITER / _BURNIN / _CHAINS (3500 / 500 / 3;
#   lower them only for plumbing runs), PBWPFVC_CORES, PBWPFVC_JM_PILOT (0 skips
#   the timing pilot), PBWPFVC_JM_HEARTBEAT (seconds between progress lines; 0 off).
#
# Usage: Rscript code/13_biotrauma_fit.R
# =============================================================================

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(here)
  library(splines)
  library(nlme)
  library(survival)
  library(JMbayes2)
  library(parallel)
})
rm(list = ls())
source("utils/config.R")

site_name  <- config$site_name
output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

source(here("code", "13_biotrauma_grid.R"))   # JM_GRID, STEP, JM_HORIZON, N_PERIODS, h_suffix
N_ITER     <- as.integer(Sys.getenv("PBWPFVC_JM_ITER",   "3500"))
N_BURNIN   <- as.integer(Sys.getenv("PBWPFVC_JM_BURNIN", "500"))
N_CHAINS   <- as.integer(Sys.getenv("PBWPFVC_JM_CHAINS", "3"))
N_CORES    <- suppressWarnings(as.integer(Sys.getenv("PBWPFVC_CORES", unset = NA)))
if (is.na(N_CORES)) N_CORES <- max(1L, detectCores() - 2L)
JM_CORES   <- min(N_CHAINS, N_CORES)
# How the day-0 baseline enters the longitudinal submodel.
#   free   (primary): log y_t ~ ... + b * log y_0, b estimated. Nests the ratio model
#          and lets a noisy single-day baseline regress to the mean instead of
#          putting that noise into the outcome.
#   offset (sensitivity): b fixed at 1, so the outcome is log(y_t / y_0), the
#          log percent change from baseline exactly (the archive's outcome).
BASELINE_FORM <- Sys.getenv("PBWPFVC_JM_BASELINE", "free")
stopifnot(BASELINE_FORM %in% c("free", "offset"))
# Association structure and sampler. ASSOC = "value" (current value on each
# cause-specific hazard) or "value_slope" (adds the current slope). MALA = 1 uses
# JMbayes2's gradient-based (MALA) update for the fixed effects, which mixes
# better on the beta-random-effect ridge that random intercepts create.
ASSOC_FORM <- Sys.getenv("PBWPFVC_JM_ASSOC", if (JM_GRID == "6h") "value" else "value_slope")
stopifnot(ASSOC_FORM %in% c("value", "value_slope"))
USE_MALA   <- identical(Sys.getenv("PBWPFVC_JM_MALA", "0"), "1")
# With the within-between decomposition the patient mean already carries the
# dose level, so the default adds no further cumulative term; "mean" adds the
# running mean VT/PFVC through the previous day, "days" the count of prior days
# above 11% (which grows with time and fights the day spline).
CUM_FORM <- Sys.getenv("PBWPFVC_JM_CUM", "none")
stopifnot(CUM_FORM %in% c("none", "mean", "days"))
CUM_TERM <- switch(CUM_FORM, none = NULL, mean = "mean_prior_vtpfvc", days = "cum_days_above")
# Effect modifier of the dose slope (PRIMARY = "disc"): the within-patient VT/PBW
# change interacts with centred log PBW/PFVC discordance and, in the adjusted
# model, with the age spline (the 11.Z adjudicator: discordance is 99%
# demographics, so a discordance interaction only means something if it survives
# an age interaction). "saturated" interacts the dose change with log PBW and log
# PFVC separately (mirrors 4k). "none" fits the dose change alone.
# "pfvc" (the PFVC-level question, 2026-09-14): is a lower PFVC, at a given
# VT/PBW, associated with a worse marker at the horizon? log PFVC (per SD) enters
# as a between-patient level term and as a divergence over time, beside the
# clinician's dose (patient mean and within-patient change); the read is the
# marker difference per SD of log PFVC at H from the joint posterior, which the
# shared random effects correct for death before H. "disc_level" is its
# companion with log PBW/PFVC. Log PBW and log PFVC are never entered together.
# "channels" (2026-09-15): the pfvc form with log PFVC replaced by its four GLI
# pieces (height, age, sex, race; 13_biotrauma_grid.R), each with a level and a
# divergence term, in BOTH submodels, in place of the demographic covariates.
# One arm only (the pieces are the demographics). The report tests whether the
# four horizon contrasts are equal: if lung size is the operative quantity they
# are, and the form collapses to the pfvc form unadjusted.
MOD_FORM <- Sys.getenv("PBWPFVC_JM_MODIFIER", "disc")
# "vtpfvc" (2026-09-17): the pfvc form told the reader's way round. Log VT/PFVC
# (the patient's mean over the window, per SD) enters beside the clinician's dose,
# so the contrast is "patients at the same VT/PBW with a different VT/PFVC"; at a
# given VT/PBW it equals the log PBW/PFVC contrast (disc_level). The hazard
# carries the index VT/PFVC in place of log PFVC.
stopifnot(MOD_FORM %in% c("disc", "saturated", "none", "pfvc", "disc_level", "channels", "vtpfvc"))
adj_label <- function(adjusted) if (MOD_FORM == "channels") "channels" else if (adjusted) "adjusted" else "unadjusted"
# Hazard interaction VT/PBW x log PFVC (secondary; 0 = the paper's main-effects set)
HAZARD_INT <- identical(Sys.getenv("PBWPFVC_JM_HAZARD_INT", "0"), "1")
# Age in the HAZARD: linear (default) or the 4-df spline. With sex and race also
# in the hazard, log PFVC is nearly a linear combination of a spline in age, so
# the two sit on a posterior ridge the sampler crawls along (MIMIC: log PFVC x
# death R-hat 3.2 with 472 deaths). The longitudinal submodel keeps the spline.
HAZARD_AGE <- Sys.getenv("PBWPFVC_JM_HAZARD_AGE", "linear")
stopifnot(HAZARD_AGE %in% c("linear", "spline"))
# Every finished fit leaves a small result file (jm_result_*.rds) and a slim
# bundle (jm_fit_*.rds: the posterior draws the report needs, not the model
# object). A fit whose result file exists with the same chain settings is
# reused, so a rerun after a crash costs only the fits that had not finished.
#   PBWPFVC_JM_FRESH=1   ignore the cache and refit everything
#   PBWPFVC_JM_RESUME=1  reuse a result file even if its chain settings differ
USE_RESUME <- identical(Sys.getenv("PBWPFVC_JM_RESUME", "0"), "1")
USE_FRESH  <- identical(Sys.getenv("PBWPFVC_JM_FRESH", "0"), "1")
# Memory: each fit runs its chains as separate processes, and a 7,000-patient
# joint model is large, so fits at a time is capped (PBWPFVC_JM_PAR, default 2)
# below the core budget, and the stored draws are thinned (PBWPFVC_JM_THIN).
N_FITS_MAX <- max(1L, as.integer(Sys.getenv("PBWPFVC_JM_PAR", "1")))
N_THIN     <- max(1L, as.integer(Sys.getenv("PBWPFVC_JM_THIN", "5")))
# Terms whose convergence the paper depends on; the manifest reports their R-hat
# beside the all-parameter maximum so a nuisance term cannot hide a converged read.
KEY_TERMS <- c("l_vtpbw_within", "l_vtpbw_within:ldisc_c", "l_vtpbw_within:age10_c",
               "^log_pfvc_sd", "^ldisc_sd", "vent_day:log_pfvc_sd", "vent_day:ldisc_sd", "^ch_", "vent_day:ch_",
               "^log_vtpfvc_sd", "vent_day:log_vtpfvc_sd", "log_vtpfvc_idx:strata\\(strata\\)death",
               "value\\(log_y\\):stratadeath", "log_pfvc:strata\\(strata\\)death",
               "vtpbw_idx:strata\\(strata\\)death")
# Progress reporting (see the MCMC block in fit_one). The pilot costs about
# PILOT_ITER / N_ITER of one chain's time.
USE_PILOT     <- !identical(Sys.getenv("PBWPFVC_JM_PILOT", "1"), "0")
PILOT_ITER    <- 300L   # burn-in 100: shorter pilots fail the adaptive-covariance Cholesky
HEARTBEAT_SEC <- as.integer(Sys.getenv("PBWPFVC_JM_HEARTBEAT", "60"))
# Heartbeat: a detached shell loop that prints elapsed time (and the pilot's
# expected finish) to this process's stderr every HEARTBEAT_SEC while the MCMC
# runs, and exits when the R process ends or stop_heartbeat() kills it.
start_heartbeat <- function(tag, t0, eta_txt) {
  if (HEARTBEAT_SEC <= 0L) return(NULL)
  pidfile <- tempfile("jm_heartbeat_")
  cmd <- sprintf(
    "echo $$ > %s; while kill -0 %d 2>/dev/null; do sleep %d; now=$(date +%%s); el=$(( now - %d )); printf '  [%%s] %s: MCMC %%d:%%02d elapsed%s\\n' \"$(date +%%H:%%M:%%S)\" $(( el / 60 )) $(( el %% 60 )) >&2; done",
    shQuote(pidfile), Sys.getpid(), HEARTBEAT_SEC, as.integer(as.numeric(t0)), tag,
    if (nzchar(eta_txt)) paste0(" (", eta_txt, ")") else "")
  system(paste("bash -c", shQuote(cmd)), wait = FALSE)   # R appends the & itself
  Sys.sleep(0.2)
  pid <- if (file.exists(pidfile)) suppressWarnings(as.integer(readLines(pidfile, n = 1))) else NA_integer_
  list(pid = pid, pidfile = pidfile)
}
stop_heartbeat <- function(hb) {
  if (is.null(hb)) return(invisible(NULL))
  if (is.finite(hb$pid)) suppressWarnings(system2("kill", as.character(hb$pid), stderr = FALSE, stdout = FALSE))
  unlink(hb$pidfile)
  invisible(NULL)
}
RHAT_GATE  <- 1.1
MIN_PATIENTS <- 20L; MIN_DEATHS <- 5L
message("=== 13_biotrauma_fit: horizon ", JM_HORIZON, "d, site ", site_name,
        ", MCMC ", N_ITER, "/", N_BURNIN, " x ", N_CHAINS, " chains on ", JM_CORES, " cores ===")
if (N_ITER < 3000L) message("*** PLUMBING setting: N_ITER < 3000; raise PBWPFVC_JM_ITER for any reported fit ***")

long_all <- read_parquet(file.path(output_dir, paste0("jm_long_", h_suffix, ".parquet")))
surv_all <- read_parquet(file.path(output_dir, paste0("jm_surv_", h_suffix, ".parquet")))
meta     <- readRDS(file.path(output_dir, paste0("jm_meta_", h_suffix, ".rds")))
message("Loaded ", nrow(long_all), " patient-days, ", nrow(surv_all), " patients")

# =============================================================================
# 13e. Marker specification
# =============================================================================
# y      : the daily column;  y0 : its index-day baseline (surv table)
# own_lag: the lagged confounder that IS this marker's own lag, dropped from its model
# random : pdDiag for the sparse plateau-measured mechanics marker, unstructured otherwise.
#          On the 6h grid the labs (creatinine, platelets, bilirubin) carry one or
#          two values in 48 hours; JMbayes2 refuses a random intercept alone for a
#          single outcome, so they get the pdDiag structure (independent intercept
#          and slope variances, the slope variance shrinking toward zero where the
#          data cannot support it), as the sparse plateau marker already does.
# offset : added before the log for markers with true zeros (NE-equivalent dose)
markers <- list(
  creatinine    = list(y = "creatinine",    y0 = "creatinine_0", own_lag = NULL,       random = if (JM_GRID == "6h") "pddiag" else "unstructured", offset = 0,    label = "Creatinine"),
  platelets     = list(y = "platelets",     y0 = "platelet_0",   own_lag = NULL,       random = if (JM_GRID == "6h") "pddiag" else "unstructured", offset = 0,    label = "Platelets"),
  bilirubin     = list(y = "bilirubin",     y0 = "bilirubin_0",  own_lag = NULL,       random = if (JM_GRID == "6h") "pddiag" else "unstructured", offset = 0,    label = "Bilirubin"),
  sf            = list(y = "sf",            y0 = "sf_0",         own_lag = "l_log_sf", random = "unstructured", offset = 0,    label = "SF ratio"),
  dp            = list(y = "dp",            y0 = "dp_0",         own_lag = NULL,       random = "pddiag",       offset = 0,    label = "Driving pressure"),
  ne_equiv_peak = list(y = "ne_equiv_peak", y0 = "ne_equiv_0",   own_lag = "l_pressor",random = "unstructured", offset = 0.01, label = "NE-equivalent dose"),
  # The hurdle's binary part: any vasoactive running in the period (NE-equivalent
  # dose > 0), a logistic mixed model (GLMMadaptive) linked to the hazards through
  # its logit. Unit-invariant, so free of the per-kg height artefact that makes
  # the dose part (ne_equiv_peak) uninterpretable against PFVC. The two together
  # are the hurdle model, fitted as two joint models on the same cohort.
  any_pressor   = list(y = "ne_equiv_peak", y0 = "ne_equiv_0",   own_lag = "l_pressor",random = "pddiag",       offset = 0,    label = "Any vasopressor", binary = TRUE)
)
for (nm in names(markers)) markers[[nm]]$name <- nm   # the output name; any_pressor shares the dose column
want_markers <- Sys.getenv("PBWPFVC_JM_MARKERS", "")
if (nzchar(want_markers)) {
  want <- trimws(strsplit(want_markers, ",")[[1]])
  if (!all(want %in% names(markers))) stop("PBWPFVC_JM_MARKERS names unknown markers: ",
                                          paste(setdiff(want, names(markers)), collapse = ", "))
  markers <- markers[want]
}
want_models <- trimws(strsplit(Sys.getenv("PBWPFVC_JM_MODELS", "main,hetero"), ",")[[1]])
stopifnot(all(want_models %in% c("main", "hetero")))

DEMO_RHS  <- "ns(age10, 4) + sex_category + race_category"
DEMO_RHS_HAZARD <- function() paste(if (HAZARD_AGE == "spline") "ns(age10, 4)" else "age10",
                                    "+ sex_category + race_category")
# Severity covariates of the longitudinal submodel. Non-respiratory SOFA (log SF
# carries the respiratory component). BMI ONLY for the pressure-derived marker:
# BMI is weight over height squared, so it carries height, which is what
# identifies log PFVC once age, sex and race are in; its chest-wall rationale
# applies to driving pressure and elastance, not to labs, oxygenation or
# vasopressors (user, 2026-09-15). The hazard keeps the paper's mortality set.
BASE_RHS  <- "np_sofa"
PRESSURE_MARKERS <- c("dp")
base_rhs_for <- function(y) if (y %in% PRESSURE_MARKERS) paste(BASE_RHS, "+ bmi") else BASE_RHS

# =============================================================================
# 13f. One fit
# =============================================================================
fit_one <- function(mk, model = c("main", "hetero"), adjusted = TRUE) {
  model <- match.arg(model)
  adj_lab <- adj_label(adjusted)
  tag <- paste(mk$name, model, adj_lab, sep = "_")
  stamp <- function(...) message(sprintf("  [%s] %s: %s", format(Sys.time(), "%H:%M:%S"), tag, paste0(...)))
  stamp("start")
  bundle_file <- file.path(output_dir, paste0("jm_fit_", tag, "_", BASELINE_FORM,
                                              if (MOD_FORM != "disc") paste0("_", MOD_FORM) else "", "_", h_suffix, ".rds"))
  rf <- result_file(mk$name, model, adj_lab)
  if (!USE_FRESH && file.exists(rf) && file.exists(bundle_file)) {
    r <- readRDS(rf)
    same <- identical(as.integer(r$n_iter), N_ITER) && identical(as.integer(r$n_burnin), N_BURNIN)
    if (same || USE_RESUME) {
      stamp("cached: ", basename(rf), if (same) "" else " (different chain settings; PBWPFVC_JM_RESUME=1)", "; no MCMC")
      return(r)
    }
    stamp("result on disk has other chain settings (", r$n_iter, "/", r$n_burnin, "); refitting")
  }
  # --- longitudinal rows: day >= 1 (day 0 is the baseline covariate), marker and lag observed
  ld <- long_all %>%
    filter(period >= 1L, !is.na(.data[[mk$y]]), !is.na(l_vtpfvc), !is.na(l_sf), !is.na(l_pressor)) %>%
    mutate(log_y = if (isTRUE(mk$binary)) as.numeric(.data[[mk$y]] > 0) else log(.data[[mk$y]] + mk$offset),
           l_log_sf = log(l_sf)) %>%
    inner_join(surv_all %>% select(hospitalization_id, np_sofa, bmi, age10, sex_category,
                                   race_category, ers_pfvc_0, vtpfvc_pt_mean, vtpbw_pt_mean,
                                   ldisc_c, log_pbw, log_pfvc, log_pfvc_sd, ldisc_sd, log_vtpfvc_sd, log_vtpfvc_idx,
                                   all_of(CHANNELS), all_of(mk$y0)),
               by = "hospitalization_id") %>%
    filter(!is.na(np_sofa), !is.na(vtpbw_pt_mean), !is.na(l_vtpbw_within),
           if (mk$y %in% PRESSURE_MARKERS) !is.na(bmi) else TRUE)
  # centred age for the dose x age interaction: uncentred, the interaction and the
  # dose main effect are collinear (age10 has a large mean relative to its spread)
  age_med <- median(ld %>% distinct(hospitalization_id, age10) %>% pull(age10))
  ld <- ld %>% mutate(age10_c = age10 - age_med)
  if (!is.null(mk$y0)) ld <- ld %>% filter(!is.na(.data[[mk$y0]])) %>%
    mutate(log_y0 = if (isTRUE(mk$binary)) as.numeric(.data[[mk$y0]] > 0) else log(.data[[mk$y0]] + mk$offset))
  # offset form: JMbayes2 rejects offset() terms, so the fixed unit coefficient is
  # applied by hand -- the response becomes log(y_t / y_0), the log percent change.
  if (!is.null(mk$y0) && BASELINE_FORM == "offset") ld <- ld %>% mutate(log_y = log_y - log_y0)
  if (model == "hetero") ld <- ld %>% filter(!is.na(ers_pfvc_0))
  n_per <- ld %>% count(hospitalization_id) %>% filter(n >= 2L)
  ld <- ld %>% filter(hospitalization_id %in% n_per$hospitalization_id)

  # --- survival rows for those patients; hazard exposures = index VT/PBW (dose)
  #     and log PFVC (size), the paper's primary parameterization
  sd_ <- surv_all %>%
    filter(hospitalization_id %in% ld$hospitalization_id) %>%
    filter(!is.na(vtpbw_idx), !is.na(log_pfvc), !is.na(sf_0), !is.na(bmi), !is.na(ch_height),
           if (MOD_FORM == "vtpfvc") is.finite(log_vtpfvc_idx) else TRUE) %>%   # the hazard keeps BMI (paper's set)
    mutate(log_sf_0 = log(sf_0))
  ld <- ld %>% filter(hospitalization_id %in% sd_$hospitalization_id)
  lv <- sort(unique(ld$hospitalization_id))
  ld$id  <- factor(ld$hospitalization_id,  levels = lv)
  sd_$id <- factor(sd_$hospitalization_id, levels = lv)
  sd_ <- sd_ %>% arrange(id)
  n_pts <- length(lv); n_deaths <- sum(sd_$event == 1L); n_extub <- sum(sd_$event == 2L)
  stamp(nrow(ld), " rows, ", n_pts, " patients, ", n_deaths, " deaths, ", n_extub, " extubations")
  # within-patient spread of the dose: a null dose slope on an exposure that
  # barely moves is a power statement, not a finding
  within_sd <- sd(ld$l_vtpbw_within)
  frac_moved <- mean(abs(ld$l_vtpbw_within) > 0.5)
  counts <- tibble(marker = mk$name, model = model, adjustment = adj_lab,
                   n_obs = nrow(ld), n_patients = n_pts, n_deaths = n_deaths, n_extubations = n_extub,
                   dose_within_sd = within_sd, frac_days_dose_moved_gt_0.5 = frac_moved)
  if (n_pts < MIN_PATIENTS || n_deaths < MIN_DEATHS) {
    stamp("skipped: too few patients or deaths")
    return(list(status = "skipped", reason = "too few patients or deaths", counts = counts))
  }

  # --- every modelled column must be finite; name the offender instead of letting
  #     nlme fail with "NA/NaN/Inf in foreign function call"
  num_cols <- intersect(c("log_y", "log_y0", "l_vtpbw_within", "vtpbw_pt_mean", "ldisc_c",
                          "log_pbw", "log_pfvc", "log_pfvc_sd", "ldisc_sd", CHANNELS, CUM_TERM,
                          if (MOD_FORM == "vtpfvc") "log_vtpfvc_sd",
                          "l_log_sf", "l_pressor", "np_sofa", if (mk$y %in% PRESSURE_MARKERS) "bmi",
                          "age10", "ers_pfvc_0"), names(ld))
  if (model != "hetero") num_cols <- setdiff(num_cols, "ers_pfvc_0")
  n_bad <- vapply(num_cols, function(v) sum(!is.finite(ld[[v]])), integer(1))
  if (any(n_bad > 0))
    stop("non-finite values in the longitudinal design: ",
         paste(sprintf("%s (%d rows)", names(n_bad)[n_bad > 0], n_bad[n_bad > 0]), collapse = ", "),
         ". Check the marker's non-positive values and the baseline covariates in 13_biotrauma_panel.R.")
  s_bad <- vapply(c("vtpbw_idx", "log_pfvc", "np_sofa", "log_sf_0", "bmi", "age10", "event_time", CHANNELS,
                    if (MOD_FORM == "vtpfvc") "log_vtpfvc_idx"),
                  function(v) sum(!is.finite(sd_[[v]])), integer(1))
  if (any(s_bad > 0))
    stop("non-finite values in the survival design: ",
         paste(sprintf("%s (%d rows)", names(s_bad)[s_bad > 0], s_bad[s_bad > 0]), collapse = ", "))

  # --- longitudinal submodel
  lag_terms <- setdiff(c("l_log_sf", "l_pressor"), mk$own_lag)
  # PRIMARY: PFVC as an effect modifier of the clinician's dose. l_vtpbw_within =
  # yesterday's VT/PBW minus the patient's mean over the course (the dose change,
  # identified within patient); its slope is modified by centred log PBW/PFVC
  # discordance (and, adjusted, by the age spline). The ratio of the interaction
  # to the main effect is the scaling exponent gamma: 0 = injury per mL/kg PBW
  # does not depend on true lung size (PBW normalizer), 1 = it scales in
  # proportion to PBW/PFVC (PFVC normalizer, the VT/PFVC model as a special case).
  # vtpbw_pt_mean = the between-patient dose level.
  mod_terms <- switch(MOD_FORM,
    # the age modification of the dose slope is LINEAR in age: a 4-df spline
    # interaction is four weakly identified parameters that fail the R-hat gate
    # even on synthetic data; the age main effect keeps its spline
    disc       = c("l_vtpbw_within * ldisc_c", if (adjusted) "l_vtpbw_within:age10_c"),
    saturated  = c("l_vtpbw_within * (log_pbw + log_pfvc)", if (adjusted) "l_vtpbw_within:age10_c"),
    none       = "l_vtpbw_within",
    pfvc       = c("l_vtpbw_within", "log_pfvc_sd", "log_pfvc_sd:vent_day"),
    disc_level = c("l_vtpbw_within", "ldisc_sd", "ldisc_sd:vent_day"),
    vtpfvc     = c("l_vtpbw_within", "log_vtpfvc_sd", "log_vtpfvc_sd:vent_day"),
    channels   = c("l_vtpbw_within", CHANNELS, paste0(CHANNELS, ":vent_day")))
  # time: linear over the 48-hour grid (the plausible shape there); a 3-df
  # natural spline over the 7-day daily grid
  time_term <- if (JM_GRID == "6h") "vent_day" else "ns(vent_day, 3)"
  rhs <- c(time_term, mod_terms, "vtpbw_pt_mean", CUM_TERM,
           if (!is.null(mk$y0) && BASELINE_FORM == "free") "log_y0",
           if (model == "hetero") "ers_pfvc_0 * l_vtpbw_within",
           lag_terms, base_rhs_for(mk$y), if (adjusted && MOD_FORM != "channels") DEMO_RHS)
  lme_formula <- as.formula(paste("log_y ~", paste(rhs, collapse = " + ")))
  random_spec <- switch(mk$random,
    pddiag       = list(id = nlme::pdDiag(~ vent_day)),
    intercept    = ~ 1 | id,
    unstructured = ~ vent_day | id)
  if (mk$random == "intercept" && ASSOC_FORM == "value_slope")
    stop("marker ", mk$y, " has a random intercept only on this grid; the slope association needs a random slope. ",
         "Use PBWPFVC_JM_ASSOC=value.")
  if (isTRUE(mk$binary)) {
    # logistic mixed model, independent intercept and slope variances (the || form)
    # Student-t penalty on the fixed effects, and a wide coefficient ceiling:
    # pressor status is persistent within patient, so the baseline-status term is
    # large; on the synthetic site it is near-deterministic (0.2% of patients off a
    # pressor at the index are on one later) and the coefficient reaches 30, which
    # is a property of the synthetic data, not of the model. A coefficient near the
    # ceiling on real data means separation and is reported as such.
    lme_fit <- GLMMadaptive::mixed_model(fixed = lme_formula, random = ~ vent_day || id, data = ld,
                                         family = binomial(), penalized = TRUE,
                                         control = list(max_coef_value = 100))
    big <- GLMMadaptive::fixef(lme_fit)[abs(GLMMadaptive::fixef(lme_fit)) > 15]
    if (length(big)) stamp("WARNING: logistic coefficients beyond 15 (near separation): ",
                           paste(sprintf("%s %.1f", names(big), big), collapse = ", "))
    stamp("logistic mixed model converged")
  } else {
    lme_fit <- lme(lme_formula, random = random_spec, data = ld,
                   control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200))
    stamp("LME converged")
  }

  # --- survival submodel: cause-specific stratified Cox (death vs extubation)
  surv_cr <- crisk_setup(as.data.frame(sd_), statusVar = "event_factor", censLevel = "censored")
  surv_cr$id <- factor(surv_cr$id, levels = lv)
  # Baseline covariates of the hazard. The SF model drops log_sf_0: it is that
  # marker's own baseline, collinear with value(log_y) on day 1 (the same
  # own-lag rule as the longitudinal submodel).
  # channels form: the size term of the hazard is the four pieces too, in place of log PFVC and the demographics
  size_haz <- if (MOD_FORM == "channels") CHANNELS else if (MOD_FORM == "vtpfvc") "log_vtpfvc_idx" else
              if (HAZARD_INT) "vtpbw_idx * log_pfvc" else "log_pfvc"
  cox_rhs <- paste(c(if (!HAZARD_INT || MOD_FORM == "channels") "vtpbw_idx", size_haz,
                     "np_sofa", if (mk$y != "sf") "log_sf_0", "bmi",
                     if (adjusted && MOD_FORM != "channels") DEMO_RHS_HAZARD()), collapse = " + ")
  cox_formula <- as.formula(paste0("Surv(event_time, status2) ~ (", cox_rhs, "):strata(strata)"))
  cox_cr <- coxph(cox_formula, data = surv_cr, x = TRUE)
  stamp("Cox converged")

  # --- joint model: value + slope association on each cause-specific hazard
  ff <- if (ASSOC_FORM == "value_slope") list(log_y = ~ value(log_y):strata + slope(log_y):strata)
        else list(log_y = ~ value(log_y):strata)
  # JMbayes2 runs the chains in C++ with no per-iteration hook, so progress has
  # to be inferred. A short single-chain pilot times the iterations and gives an
  # expected finish; a heartbeat (a detached shell loop) then prints the elapsed
  # time every HEARTBEAT_SEC while the real run is silent. PBWPFVC_JM_PILOT=0
  # skips the pilot.
  fit_jm <- function(n_iter, n_burnin, n_chains, cores)
    jm(cox_cr, lme_fit, time_var = "vent_day", data_Surv = surv_cr, id_var = "id",
       functional_forms = ff, n_iter = n_iter, n_burnin = n_burnin, n_thin = N_THIN,
       n_chains = n_chains, cores = cores, control = list(MALA = USE_MALA))
  eta_txt <- ""
  if (USE_PILOT) {
    t_pilot <- Sys.time()
    invisible(fit_jm(PILOT_ITER, 100L, 1L, 1L))
    sec_per_iter <- as.numeric(difftime(Sys.time(), t_pilot, units = "secs")) / PILOT_ITER
    # chains run in parallel across JM_CORES; parallel chains slow each other a little
    est_sec <- sec_per_iter * N_ITER * ceiling(N_CHAINS / JM_CORES) * 1.15
    eta_txt <- sprintf("expected %.1f min, finish about %s", est_sec / 60,
                       format(Sys.time() + est_sec, "%H:%M"))
    stamp(sprintf("pilot: %.3f s per iteration on one chain; %s", sec_per_iter, eta_txt))
  }
  t_jm <- Sys.time()
  stamp(sprintf("MCMC running: %d iterations x %d chains on %d cores (silent; heartbeat every %ds)",
                N_ITER, N_CHAINS, JM_CORES, HEARTBEAT_SEC))
  hb <- start_heartbeat(tag, t_jm, eta_txt)
  jm_fit <- tryCatch(fit_jm(N_ITER, N_BURNIN, N_CHAINS, JM_CORES), finally = stop_heartbeat(hb))
  acc_b <- mean(jm_fit$acc_rates$b, na.rm = TRUE)
  stamp(sprintf("JM done (%.1f min); random-effects acceptance %.3f, fixed-effects acceptance %.3f",
                as.numeric(difftime(Sys.time(), t_jm, units = "mins")), acc_b,
                mean(unlist(jm_fit$acc_rates$betas), na.rm = TRUE)))

  # --- estimates: every block of summary(jm) that carries a coefficient table
  s <- summary(jm_fit)
  pull_block <- function(tbl, block) {
    if (!(is.matrix(tbl) || is.data.frame(tbl)) || !"Mean" %in% colnames(tbl)) return(NULL)
    tibble(block = block, term = rownames(tbl),
           estimate = tbl[, "Mean"], sd = tbl[, "StDev"],
           lo = tbl[, "2.5%"], hi = tbl[, "97.5%"],
           rhat = if ("Rhat" %in% colnames(tbl)) tbl[, "Rhat"] else NA_real_)
  }
  est <- bind_rows(
    pull_block(s$Survival, "survival"),
    map_dfr(grep("^Outcome", names(s), value = TRUE), function(nm) pull_block(s[[nm]], "longitudinal"))
  ) %>%
    mutate(block = if_else(block == "survival" & grepl("value\\(|slope\\(", term), "association", block)) %>%
    bind_cols(counts[rep(1L, nrow(.)), ]) %>%
    mutate(baseline_form = BASELINE_FORM, horizon_days = JM_HORIZON, site = site_name)
  # --- scaling exponent gamma = (dose x discordance) / dose, from the joint posterior
  bd <- do.call(rbind, jm_fit$mcmc$betas1)
  if (is.null(colnames(bd))) colnames(bd) <- names(fixef(lme_fit))
  scaling <- if (MOD_FORM == "disc" && all(c("l_vtpbw_within", "l_vtpbw_within:ldisc_c") %in% colnames(bd))) {
    g <- bd[, "l_vtpbw_within:ldisc_c"] / bd[, "l_vtpbw_within"]
    tibble(marker = mk$name, model = model, adjustment = adj_lab,
           dose_slope = mean(bd[, "l_vtpbw_within"]),
           dose_slope_lo = quantile(bd[, "l_vtpbw_within"], 0.025), dose_slope_hi = quantile(bd[, "l_vtpbw_within"], 0.975),
           p_dose_slope_gt0 = mean(bd[, "l_vtpbw_within"] > 0),
           interaction = mean(bd[, "l_vtpbw_within:ldisc_c"]),
           interaction_lo = quantile(bd[, "l_vtpbw_within:ldisc_c"], 0.025), interaction_hi = quantile(bd[, "l_vtpbw_within:ldisc_c"], 0.975),
           gamma_median = median(g), gamma_lo = quantile(g, 0.025), gamma_hi = quantile(g, 0.975),
           p_gamma_gt0 = mean(g > 0), p_gamma_gt_half = mean(g > 0.5), p_gamma_lt1 = mean(g < 1),
           dose_within_sd = within_sd, frac_days_dose_moved_gt_0.5 = frac_moved,
           n_patients = n_pts, horizon_days = JM_HORIZON, site = site_name)
  } else NULL
  max_rhat <- max(est$rhat, na.rm = TRUE)
  gate <- is.finite(max_rhat) && max_rhat <= RHAT_GATE
  worst <- est %>% slice_max(rhat, n = 3, with_ties = FALSE)
  key <- est %>% filter(grepl(paste(KEY_TERMS, collapse = "|"), term))
  key_rhat <- if (nrow(key)) max(key$rhat, na.rm = TRUE) else NA_real_
  key_gate <- is.finite(key_rhat) && key_rhat <= RHAT_GATE
  stamp(sprintf("max R-hat %.3f (%s); key terms %.3f (%s); worst: %s",
                max_rhat, if (gate) "passes" else "FAILS gate",
                key_rhat, if (key_gate) "pass" else "fail",
                paste(sprintf("%s %.2f", worst$term, worst$rhat), collapse = ", ")))

  # --- Q3: the hazard exposures on the DEATH hazard, plain Cox vs inside the JM.
  #     log PFVC is the paper's primary size term (the absorption read); VT/PBW
  #     is the dose and is reported beside it.
  cox_tbl <- summary(cox_cr)$coefficients
  absorption <- map_dfr(c("log_pfvc", "vtpbw_idx", "log_vtpfvc_idx"), function(tm) {
    cox_row <- grep(paste0("^", tm, ":strata\\(strata\\)death$"), rownames(cox_tbl))
    jm_row  <- est %>% filter(block == "survival", grepl(paste0("^", tm, ":"), term), grepl("death", term))
    tibble(marker = mk$name, model = model, adjustment = adj_lab, term = tm,
           cox_log_hr = if (length(cox_row) == 1L) cox_tbl[cox_row, "coef"] else NA_real_,
           cox_se     = if (length(cox_row) == 1L) cox_tbl[cox_row, "se(coef)"] else NA_real_,
           jm_log_hr  = if (nrow(jm_row) == 1L) jm_row$estimate else NA_real_,
           jm_sd      = if (nrow(jm_row) == 1L) jm_row$sd else NA_real_)
  }) %>%
    mutate(absorbed = cox_log_hr - jm_log_hr,
           absorbed_frac = if_else(is.finite(cox_log_hr) & cox_log_hr != 0, absorbed / cox_log_hr, NA_real_),
           n_patients = n_pts, n_deaths = n_deaths, baseline_form = BASELINE_FORM,
           horizon_days = JM_HORIZON, site = site_name)

  # terms with predvars: carries the ns() knots so the report can rebuild the
  # fixed-effects design on a prediction grid without re-deriving the basis
  # slim bundle: the posterior draws the report uses (fixed effects, association),
  # the LME (for the coefficient names), the long data and the terms object; not
  # the joint-model object, whose design matrices and quadrature arrays are what
  # made the full bundle gigabytes
  mf_terms <- terms(model.frame(lme_formula, data = ld))
  saveRDS(list(jm = list(mcmc = jm_fit$mcmc[c("betas1", "alphas")], acc_rates = jm_fit$acc_rates),
               lme = lme_fit, marker = mk, model = model, binary = isTRUE(mk$binary),
               adjusted = adjusted, counts = counts, long_data = ld,
               lme_formula = lme_formula, mf_terms = mf_terms, assoc_form = ASSOC_FORM,
               mod_form = MOD_FORM, baseline_form = BASELINE_FORM, horizon = JM_HORIZON,
               n_iter = N_ITER, n_burnin = N_BURNIN, n_thin = N_THIN),
          bundle_file)
  rm(jm_fit, surv_cr, cox_cr); invisible(gc())
  result <- list(status = if (gate) "converged" else "rhat_fail", reason = NA_character_,
                 counts = counts, estimates = est, absorption = absorption, scaling = scaling,
                 max_rhat = max_rhat, key_rhat = key_rhat,
                 worst_terms = paste(sprintf("%s %.2f", worst$term, worst$rhat), collapse = "; "),
                 acc_b = acc_b, n_iter = N_ITER, n_burnin = N_BURNIN, n_thin = N_THIN)
  # the small result list also goes to disk, so a cluster failure after the fits
  # finished loses nothing (the master collects these files if the cluster dies)
  saveRDS(result, result_file(mk$name, model, adj_lab))
  result
}
# per-fit result file (aggregates only; beside the bundles in intermediate/)
result_file <- function(marker, model, adj_lab)
  file.path(output_dir, paste0("jm_result_", marker, "_", model, "_", adj_lab, "_", BASELINE_FORM,
                               if (MOD_FORM != "disc") paste0("_", MOD_FORM) else "", "_", h_suffix, ".rds"))
RUN_START <- Sys.time()

# =============================================================================
# 13g. Run the model set
# =============================================================================
jobs <- expand_grid(marker = names(markers), model = want_models, adjusted = c(TRUE, FALSE)) %>%
  filter(!(model == "hetero" & !adjusted)) %>%          # heterogeneity: adjusted only
  filter(!(MOD_FORM == "channels" & adjusted))          # channels: one arm, the pieces are the demographics
run_job <- function(marker, model, adjusted) {
  tryCatch(fit_one(markers[[marker]], model, adjusted),
           error = function(e) {
             message(sprintf("  %s/%s/%s FAILED: %s", marker, model, adj_label(adjusted), conditionMessage(e)))
             list(status = "failed", reason = conditionMessage(e),
                  counts = tibble(marker = markers[[marker]]$name, model = model,
                                  adjustment = adj_label(adjusted)))
           })
}
# Fits run in parallel across PSOCK workers, each fit using one core per chain,
# as many fits at once as the core budget allows (N_CORES / N_CHAINS). Worker
# output is forwarded to this console (outfile = ""), so stamps and heartbeats
# from concurrent fits interleave, each prefixed by its fit tag.
N_FITS_PAR <- max(1L, min(nrow(jobs), N_CORES %/% N_CHAINS, N_FITS_MAX))
if (N_FITS_PAR > 1L) {
  message("Running ", nrow(jobs), " fits, ", N_FITS_PAR, " at a time (", N_CHAINS, " chains each)")
  cl <- makeCluster(N_FITS_PAR, type = "PSOCK", outfile = "")
  clusterEvalQ(cl, suppressPackageStartupMessages({
    library(tidyverse); library(splines); library(nlme); library(survival); library(JMbayes2)
  }))
  clusterExport(cl, setdiff(ls(envir = .GlobalEnv), "cl"), envir = .GlobalEnv)
  results <- tryCatch(clusterMap(cl, run_job, jobs$marker, jobs$model, jobs$adjusted, SIMPLIFY = FALSE, USE.NAMES = FALSE),
                      error = function(e) {
                        message("\nCLUSTER FAILED while collecting results (", conditionMessage(e),
                                "); collecting the per-fit result files written during this run instead")
                        NULL
                      })
  try(stopCluster(cl), silent = TRUE)
  if (is.null(results)) results <- pmap(jobs, function(marker, model, adjusted) {
    f <- result_file(markers[[marker]]$name, model, adj_label(adjusted))
    if (file.exists(f)) readRDS(f) else
      list(status = "failed", reason = "cluster failed before this fit's result was written; rerun (finished fits are cached)",
           counts = tibble(marker = markers[[marker]]$name, model = model, adjustment = adj_label(adjusted)))
  })
} else {
  results <- pmap(jobs, run_job)
}

manifest <- map_dfr(results, function(r)
  r$counts %>% mutate(status = r$status, reason = r$reason,
                      max_rhat = if (is.null(r$max_rhat)) NA_real_ else r$max_rhat,
                      key_terms_rhat = if (is.null(r$key_rhat)) NA_real_ else r$key_rhat,
                      worst_terms = if (is.null(r$worst_terms)) NA_character_ else r$worst_terms,
                      acc_random_effects = if (is.null(r$acc_b)) NA_real_ else r$acc_b)) %>%
  mutate(grid = JM_GRID, baseline_form = BASELINE_FORM, assoc_form = ASSOC_FORM, modifier_form = MOD_FORM,
         hazard_age = HAZARD_AGE, mala = USE_MALA, horizon_days = JM_HORIZON,
         n_iter = N_ITER, n_burnin = N_BURNIN, n_chains = N_CHAINS, n_thin = N_THIN, site = site_name)
estimates  <- map_dfr(results, "estimates")
absorption <- map_dfr(results, "absorption")
scaling    <- map_dfr(results, "scaling")

# Output tag: non-default forms (baseline offset, saturated or no modifier) get
# their own files so a sensitivity run never overwrites the primary's rows.
out_tag <- paste0(if (BASELINE_FORM == "offset") "offset_" else "",
                  if (MOD_FORM != "disc") paste0(MOD_FORM, "_") else "",
                  h_suffix, "_", site_name)
# Merge on write: a run restricted to some markers (PBWPFVC_JM_MARKERS) replaces
# only its own marker/model/adjustment rows in each table and keeps the rest, so
# the full set can be assembled from several runs (a lab-only rerun, a longer-
# chain rerun of one marker). PBWPFVC_JM_FRESH=1 discards the existing tables.
if (identical(Sys.getenv("PBWPFVC_JM_FRESH", "0"), "1"))
  for (nm in c("manifest", "estimates", "absorption", "scaling"))
    unlink(file.path(final_dir, paste0("jm_", nm, "_", out_tag, ".csv")))
merge_write <- function(new, name) {
  path <- file.path(final_dir, paste0("jm_", name, "_", out_tag, ".csv"))
  if (file.exists(path) && !identical(Sys.getenv("PBWPFVC_JM_FRESH", "0"), "1") && nrow(new)) {
    old <- read_csv(path, show_col_types = FALSE)
    keys <- new %>% distinct(marker, model, adjustment)
    old  <- old %>% anti_join(keys, by = c("marker", "model", "adjustment"))
    # both as text: the CSV holds text, and a column that is all NA reads back as
    # logical, which bind_rows refuses to combine with the new table's type
    as_text <- function(d) d %>% mutate(across(everything(), as.character))
    new  <- bind_rows(as_text(old), as_text(new))
    message("  ", name, ": kept ", nrow(old), " rows from other markers")
  }
  if (nrow(new)) write_csv(new, path)
}
merge_write(manifest,   "manifest")
merge_write(estimates,  "estimates")
merge_write(absorption, "absorption")
merge_write(scaling,    "scaling")

message("\n========== 13_biotrauma_fit SUMMARY (", h_suffix, ") ==========")
print(as.data.frame(manifest %>% select(any_of(c("marker", "model", "adjustment", "status", "reason",
                                                 "n_patients", "n_deaths", "max_rhat", "key_terms_rhat",
                                                 "worst_terms")))),
      row.names = FALSE)
message("Estimates: ", nrow(estimates), " rows -> ", final_dir)
