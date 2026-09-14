# =============================================================================
# Script 13 (fit): Biotrauma joint models -- one shared-random-effects JM per marker
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# For each organ-injury marker, fits a joint model (JMbayes2) that links
#   longitudinal submodel   log marker on day t  ~  spline(day) + previous-day
#                           VT/PFVC + prior days above 11% + baseline marker +
#                           previous-day confounders + baseline severity and
#                           demographics; random intercept and slope per patient
#   survival submodel       cause-specific stratified Cox, death vs extubation,
#                           with the index-day VT/PFVC and baseline covariates
#   association             current value and current slope of the marker on
#                           each cause-specific hazard
#
# Three questions (docs/joint_model_plan_2026-09.md, section 3):
#   Q1  the previous-day strain coefficient in the longitudinal submodel: a
#       conditional, within-patient dose-response (not a policy effect; the 11.*
#       g-methods are the causal version)
#   Q2  the value and slope association parameters, per cause
#   Q3  the index-day strain coefficient on the death hazard in the plain
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
# Environment knobs: PBWPFVC_JM_HORIZON (7), PBWPFVC_JM_MARKERS (comma list),
#   PBWPFVC_JM_MODELS (main,hetero), PBWPFVC_JM_BASELINE (free | offset; offset
#   fixes the baseline coefficient at 1 = the log percent-change outcome, written
#   with an offset_ prefix), PBWPFVC_JM_ITER / _BURNIN / _CHAINS (3500 / 500 / 3;
#   lower them only for plumbing runs), PBWPFVC_CORES.
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

JM_HORIZON <- as.integer(Sys.getenv("PBWPFVC_JM_HORIZON", "7"))
h_suffix   <- paste0(JM_HORIZON, "d")
N_ITER     <- as.integer(Sys.getenv("PBWPFVC_JM_ITER",   "3500"))
N_BURNIN   <- as.integer(Sys.getenv("PBWPFVC_JM_BURNIN", "500"))
N_CHAINS   <- as.integer(Sys.getenv("PBWPFVC_JM_CHAINS", "3"))
N_CORES    <- suppressWarnings(as.integer(Sys.getenv("PBWPFVC_CORES", unset = NA)))
if (is.na(N_CORES)) N_CORES <- max(1L, detectCores() - 1L)
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
ASSOC_FORM <- Sys.getenv("PBWPFVC_JM_ASSOC", "value_slope")
stopifnot(ASSOC_FORM %in% c("value", "value_slope"))
USE_MALA   <- identical(Sys.getenv("PBWPFVC_JM_MALA", "0"), "1")
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
# random : pdDiag for the sparse plateau-measured mechanics marker, unstructured otherwise
# offset : added before the log for markers with true zeros (NE-equivalent dose)
markers <- list(
  creatinine    = list(y = "creatinine",    y0 = "creatinine_0", own_lag = NULL,       random = "unstructured", offset = 0,    label = "Creatinine"),
  platelets     = list(y = "platelets",     y0 = "platelet_0",   own_lag = NULL,       random = "unstructured", offset = 0,    label = "Platelets"),
  bilirubin     = list(y = "bilirubin",     y0 = "bilirubin_0",  own_lag = NULL,       random = "unstructured", offset = 0,    label = "Bilirubin"),
  sf            = list(y = "sf",            y0 = "sf_0",         own_lag = "l_log_sf", random = "unstructured", offset = 0,    label = "SF ratio"),
  dp            = list(y = "dp",            y0 = "dp_0",         own_lag = NULL,       random = "pddiag",       offset = 0,    label = "Driving pressure"),
  ne_equiv_peak = list(y = "ne_equiv_peak", y0 = "ne_equiv_0",   own_lag = "l_pressor",random = "unstructured", offset = 0.01, label = "NE-equivalent dose")
)
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
BASE_RHS  <- "sofa_total + bmi"

# =============================================================================
# 13f. One fit
# =============================================================================
fit_one <- function(mk, model = c("main", "hetero"), adjusted = TRUE) {
  model <- match.arg(model)
  adj_lab <- if (adjusted) "adjusted" else "unadjusted"
  tag <- paste(mk$y, model, adj_lab, sep = "_")
  stamp <- function(...) message(sprintf("  [%s] %s: %s", format(Sys.time(), "%H:%M:%S"), tag, paste0(...)))
  stamp("start")

  # --- longitudinal rows: day >= 1 (day 0 is the baseline covariate), marker and lag observed
  ld <- long_all %>%
    filter(vent_day >= 1L, !is.na(.data[[mk$y]]), !is.na(l_vtpfvc), !is.na(l_sf), !is.na(l_pressor)) %>%
    mutate(log_y = log(.data[[mk$y]] + mk$offset), l_log_sf = log(l_sf)) %>%
    inner_join(surv_all %>% select(hospitalization_id, sofa_total, bmi, age10, sex_category,
                                   race_category, ers_pfvc_0, all_of(mk$y0)),
               by = "hospitalization_id") %>%
    filter(!is.na(sofa_total), !is.na(bmi))
  if (!is.null(mk$y0)) ld <- ld %>% filter(!is.na(.data[[mk$y0]])) %>%
    mutate(log_y0 = log(.data[[mk$y0]] + mk$offset))
  # offset form: JMbayes2 rejects offset() terms, so the fixed unit coefficient is
  # applied by hand -- the response becomes log(y_t / y_0), the log percent change.
  if (!is.null(mk$y0) && BASELINE_FORM == "offset") ld <- ld %>% mutate(log_y = log_y - log_y0)
  if (model == "hetero") ld <- ld %>% filter(!is.na(ers_pfvc_0))
  n_per <- ld %>% count(hospitalization_id) %>% filter(n >= 2L)
  ld <- ld %>% filter(hospitalization_id %in% n_per$hospitalization_id)

  # --- survival rows for those patients; index-day strain from day 0
  idx_strain <- long_all %>% filter(vent_day == 0L) %>% select(hospitalization_id, vtpfvc_idx = vtpfvc)
  sd_ <- surv_all %>%
    filter(hospitalization_id %in% ld$hospitalization_id) %>%
    inner_join(idx_strain, by = "hospitalization_id") %>%
    filter(!is.na(vtpfvc_idx), !is.na(sf_0)) %>%
    mutate(log_sf_0 = log(sf_0))
  ld <- ld %>% filter(hospitalization_id %in% sd_$hospitalization_id)
  lv <- sort(unique(ld$hospitalization_id))
  ld$id  <- factor(ld$hospitalization_id,  levels = lv)
  sd_$id <- factor(sd_$hospitalization_id, levels = lv)
  sd_ <- sd_ %>% arrange(id)
  n_pts <- length(lv); n_deaths <- sum(sd_$event == 1L); n_extub <- sum(sd_$event == 2L)
  stamp(nrow(ld), " rows, ", n_pts, " patients, ", n_deaths, " deaths, ", n_extub, " extubations")
  counts <- tibble(marker = mk$y, model = model, adjustment = adj_lab,
                   n_obs = nrow(ld), n_patients = n_pts, n_deaths = n_deaths, n_extubations = n_extub)
  if (n_pts < MIN_PATIENTS || n_deaths < MIN_DEATHS) {
    stamp("skipped: too few patients or deaths")
    return(list(status = "skipped", reason = "too few patients or deaths", counts = counts))
  }

  # --- longitudinal submodel
  lag_terms <- setdiff(c("l_log_sf", "l_pressor"), mk$own_lag)
  rhs <- c("ns(vent_day, 3)", "l_vtpfvc", "cum_days_above",
           if (!is.null(mk$y0) && BASELINE_FORM == "free") "log_y0",
           if (model == "hetero") "ers_pfvc_0 * l_vtpfvc",
           lag_terms, BASE_RHS, if (adjusted) DEMO_RHS)
  lme_formula <- as.formula(paste("log_y ~", paste(rhs, collapse = " + ")))
  random_spec <- if (mk$random == "pddiag") list(id = nlme::pdDiag(~ vent_day)) else ~ vent_day | id
  lme_fit <- lme(lme_formula, random = random_spec, data = ld,
                 control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200))
  stamp("LME converged")

  # --- survival submodel: cause-specific stratified Cox (death vs extubation)
  surv_cr <- crisk_setup(as.data.frame(sd_), statusVar = "event_factor", censLevel = "censored")
  surv_cr$id <- factor(surv_cr$id, levels = lv)
  # Baseline covariates of the hazard. The SF model drops log_sf_0: it is that
  # marker's own baseline, collinear with value(log_y) on day 1 (the same
  # own-lag rule as the longitudinal submodel).
  cox_rhs <- paste(c("vtpfvc_idx", "sofa_total", if (mk$y != "sf") "log_sf_0", "bmi",
                     if (adjusted) DEMO_RHS), collapse = " + ")
  cox_formula <- as.formula(paste0("Surv(event_time, status2) ~ (", cox_rhs, "):strata(strata)"))
  cox_cr <- coxph(cox_formula, data = surv_cr, x = TRUE)
  stamp("Cox converged")

  # --- joint model: value + slope association on each cause-specific hazard
  ff <- if (ASSOC_FORM == "value_slope") list(log_y = ~ value(log_y):strata + slope(log_y):strata)
        else list(log_y = ~ value(log_y):strata)
  t_jm <- Sys.time()
  jm_fit <- jm(cox_cr, lme_fit, time_var = "vent_day", data_Surv = surv_cr, id_var = "id",
               functional_forms = ff,
               n_iter = N_ITER, n_burnin = N_BURNIN, n_thin = 1L,
               n_chains = N_CHAINS, cores = JM_CORES,
               control = list(MALA = USE_MALA))
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
  max_rhat <- max(est$rhat, na.rm = TRUE)
  gate <- is.finite(max_rhat) && max_rhat <= RHAT_GATE
  worst <- est %>% slice_max(rhat, n = 3, with_ties = FALSE)
  stamp(sprintf("max R-hat %.3f (%s); worst: %s", max_rhat, if (gate) "passes" else "FAILS gate",
                paste(sprintf("%s %.2f", worst$term, worst$rhat), collapse = ", ")))

  # --- Q3: index-day strain on the DEATH hazard, plain Cox vs inside the JM
  cox_tbl <- summary(cox_cr)$coefficients
  cox_row <- grep("^vtpfvc_idx:strata\\(strata\\)death$|^vtpfvc_idx:strata\\(strata\\)death", rownames(cox_tbl))
  jm_row  <- est %>% filter(block == "survival", grepl("^vtpfvc_idx", term), grepl("death", term))
  absorption <- tibble(
    marker = mk$y, model = model, adjustment = adj_lab,
    cox_log_hr = if (length(cox_row) == 1L) cox_tbl[cox_row, "coef"] else NA_real_,
    cox_se     = if (length(cox_row) == 1L) cox_tbl[cox_row, "se(coef)"] else NA_real_,
    jm_log_hr  = if (nrow(jm_row) == 1L) jm_row$estimate else NA_real_,
    jm_sd      = if (nrow(jm_row) == 1L) jm_row$sd else NA_real_) %>%
    mutate(absorbed = cox_log_hr - jm_log_hr,
           absorbed_frac = if_else(is.finite(cox_log_hr) & cox_log_hr != 0, absorbed / cox_log_hr, NA_real_),
           n_patients = n_pts, n_deaths = n_deaths, baseline_form = BASELINE_FORM,
           horizon_days = JM_HORIZON, site = site_name)

  # terms with predvars: carries the ns() knots so the report can rebuild the
  # fixed-effects design on a prediction grid without re-deriving the basis
  mf_terms <- terms(model.frame(lme_formula, data = ld))
  saveRDS(list(jm = jm_fit, lme = lme_fit, cox = cox_cr, marker = mk, model = model,
               adjusted = adjusted, counts = counts, long_data = ld, surv_cr = surv_cr,
               lme_formula = lme_formula, mf_terms = mf_terms, assoc_form = ASSOC_FORM,
               baseline_form = BASELINE_FORM, horizon = JM_HORIZON),
          file.path(output_dir, paste0("jm_fit_", tag, "_", BASELINE_FORM, "_", h_suffix, ".rds")))
  list(status = if (gate) "converged" else "rhat_fail", reason = NA_character_,
       counts = counts, estimates = est, absorption = absorption, max_rhat = max_rhat,
       acc_b = acc_b)
}

# =============================================================================
# 13g. Run the model set
# =============================================================================
jobs <- expand_grid(marker = names(markers), model = want_models, adjusted = c(TRUE, FALSE)) %>%
  filter(!(model == "hetero" & !adjusted))   # heterogeneity: adjusted only
results <- pmap(jobs, function(marker, model, adjusted) {
  tryCatch(fit_one(markers[[marker]], model, adjusted),
           error = function(e) {
             message(sprintf("  %s/%s/%s FAILED: %s", marker, model,
                             if (adjusted) "adjusted" else "unadjusted", conditionMessage(e)))
             list(status = "failed", reason = conditionMessage(e),
                  counts = tibble(marker = markers[[marker]]$y, model = model,
                                  adjustment = if (adjusted) "adjusted" else "unadjusted"))
           })
})

manifest <- map_dfr(results, function(r)
  r$counts %>% mutate(status = r$status, reason = r$reason,
                      max_rhat = if (is.null(r$max_rhat)) NA_real_ else r$max_rhat,
                      acc_random_effects = if (is.null(r$acc_b)) NA_real_ else r$acc_b)) %>%
  mutate(baseline_form = BASELINE_FORM, assoc_form = ASSOC_FORM, mala = USE_MALA, horizon_days = JM_HORIZON,
         n_iter = N_ITER, n_burnin = N_BURNIN, n_chains = N_CHAINS, site = site_name)
estimates  <- map_dfr(results, "estimates")
absorption <- map_dfr(results, "absorption")

out_tag <- paste0(if (BASELINE_FORM == "offset") "offset_" else "", h_suffix, "_", site_name)
write_csv(manifest,   file.path(final_dir, paste0("jm_manifest_",   out_tag, ".csv")))
write_csv(estimates,  file.path(final_dir, paste0("jm_estimates_",  out_tag, ".csv")))
write_csv(absorption, file.path(final_dir, paste0("jm_absorption_", out_tag, ".csv")))

message("\n========== 13_biotrauma_fit SUMMARY (", h_suffix, ") ==========")
print(as.data.frame(manifest %>% select(any_of(c("marker", "model", "adjustment", "status", "reason",
                                                 "n_patients", "n_deaths", "max_rhat", "acc_random_effects")))),
      row.names = FALSE)
message("Estimates: ", nrow(estimates), " rows -> ", final_dir)
