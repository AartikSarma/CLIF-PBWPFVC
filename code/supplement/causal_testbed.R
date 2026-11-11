# =============================================================================
# Supplement (causal): a planted-truth simulation study of the complementary designs
# =============================================================================
#
# causal_triangulation.R lays fifteen designs side by side, each failing under a
# different assumption. This script reports what each returns, and what the
# pre-specified decision rule concludes, when the truth is known. The protocol
# (aims, data-generating mechanisms, estimand, methods, performance measures) is
# code/supplement/CAUSAL_TESTBED_PROTOCOL.md; this header is the short form.
#
# The generator (causal_generator.R) keeps the site's patients and generates the
# marker, its measurement times and the patients' exits. It is tuned so that the
# analysis reads the same noise, mean slopes and sampling off a simulated site as off
# the real one (causal_testbed_realism_). Stages, each a set of scenarios:
#
#   matched     the generator reduced to the analysis model, under eight truths (no
#               effect; strain; five confounders with no strain effect, each sized to
#               look to the unadjusted model exactly like the strain effect; strain
#               with a confounder). A check that the estimators are right.
#   realistic   the same truths, and more confounders (several inputs at once, the
#               inputs in GLI's own proportions, an age effect that starts at 65, half
#               and twice the size, and the age gradient observed in the controls),
#               from a generator with curved trajectories, autocorrelated heavy-tailed
#               errors, an unmeasured severity, unmodelled variation in lung size, and
#               measurement, liberation and death that depend on the marker.
#   stress      one departure at a time from the realistic base, with and without a
#               strain effect: strongly informative dropout or measurement, heavier
#               tails, stronger autocorrelation, a threshold in dose, an effect that
#               saturates or is delayed, an effect that differs in the patients a
#               within-person design can use, positive pressure before the recorded
#               intubation, a tidal volume set by size and severity (with and without
#               the 6-8 gate), error in charted height, more variation in true lung
#               size, and strain that follows another denominator (GLI at age 25, GLI
#               without race).
#   power       the strain effect at 0.25 to 1.5 times the planted size.
#   timecourse  the segment-by-segment estimate of fit_time_course() and its
#               cumulative effect, beside the straight line, when the effect accrues
#               steadily, saturates, or starts late.
#   illness_death  causal_illness_death.R's analysis on simulated sites: the marker as an
#               injury state, death following the marker, strain directly, or age
#               alone. There is no closed-form truth for a hazard gradient, so the
#               table gives the mean estimate, whether the standard error is right,
#               and how often each transition's gradient is called non-zero.
#   multisite   four sites that differ in effect, confounding and noise, pooled by
#               common-effect and by random-effects inverse variance.
#   crosscheck  the within estimator against a linear mixed model (nlme) on the same
#               simulated sites.
#
# Every design is fitted with the index VT/PBW held and with it dropped (column vtpbw),
# each scored against its own truth, and the decision rule requires both.
# Every figure in the tables carries its Monte Carlo standard error. Only the log-scale
# labs are simulated: the vasopressor and OSI outcomes are not.
#
# Writes to final/supplement/ (final/ungated/supplement/ with PBWPFVC_VTPBW_GATE=0):
#   causal_testbed_scenarios_<site>.csv    every scenario's parameters
#   causal_testbed_operating_<site>.csv    marker x scenario x design: truth, bias,
#                                          empirical and model SE, coverage, rejection
#                                          rate, MDE, verdict; each with its MCSE
#   causal_testbed_decision_<site>.csv     how often the decision rule, and two
#                                          single-analysis comparators, say "supports"
#   causal_testbed_power_<site>.csv        rejection rate by effect size
#   causal_testbed_timecourse_<site>.csv   segments, cumulative effect and straight line
#   causal_testbed_illness_death_<site>.csv  the illness-death gradients under each truth
#   causal_testbed_multisite_<site>.csv    pooled bias, coverage and rejection rate
#   causal_testbed_realism_<site>.csv      observed against simulated summaries
#   causal_testbed_crosscheck_<site>.csv   within estimator against the mixed model
#   causal_testbed_<site>.pdf              the tables as figures
# No patient-level row is written.
#
# Settings (environment):
#   PBWPFVC_CAUSAL_REPS        replicates per scenario (default 1000)
#   PBWPFVC_CAUSAL_EFFECT      the planted strain effect, log units per day per log
#                              unit of VT/PFVC, in the direction of injury (0.15)
#   PBWPFVC_CAUSAL_MARKERS     platelets,creatinine
#   PBWPFVC_CAUSAL_STAGES      matched,realistic,stress,power,timecourse,illness_death,multisite,crosscheck
#   PBWPFVC_CAUSAL_WORKERS     parallel R processes (default: cores - 1, at most 8)
#   PBWPFVC_CAUSAL_RESUME      1 = reuse the stage results a stopped run left behind
# Each marker's stages are saved under causal_testbed_checkpoints/ as they finish (Monte
# Carlo summaries only, no patient-level row). A run that stops late is finished with
# PBWPFVC_CAUSAL_RESUME=1 and the same settings. A saved stage is reused only if the
# replicates, effect, seed and stages match; it is NOT checked against the code, so
# after changing the generator or an estimator run without RESUME. Do not edit this
# file while it runs: Rscript reads it as it goes.
# Runtime: about 0.4 s per simulated site per process; the default is roughly 50,000
# sites per marker (a few hours on four cores). PBWPFVC_CAUSAL_REPS=100 is a fast look.
# Needs scripts 01-03 for the ventilated cohort and for PBWPFVC_COHORT=nosupport.
#
# Usage: uvr run code/supplement/causal_testbed.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(here)
})

rm(list = ls())
source("utils/config.R")
site_name <- config$site_name
final_dir <- final_dir_for("supplement")
source(here::here("code", "20_biotrauma_grid.R"))            # height_fingerprint(), pfvc_channels(), sorted_term()
source(here::here("code", "supplement", "causal_common.R"))
source(here::here("code", "supplement", "causal_generator.R"))

env_list <- function(name, default) strsplit(Sys.getenv(name, default), ",")[[1]]
N_REPS   <- as.integer(Sys.getenv("PBWPFVC_CAUSAL_REPS", "1000"))
EFFECT   <- as.numeric(Sys.getenv("PBWPFVC_CAUSAL_EFFECT", "0.15"))
MARKERS  <- env_list("PBWPFVC_CAUSAL_MARKERS", "platelets,creatinine")
STAGES   <- env_list("PBWPFVC_CAUSAL_STAGES", "matched,realistic,stress,power,timecourse,illness_death,multisite,crosscheck")
WORKERS  <- as.integer(Sys.getenv("PBWPFVC_CAUSAL_WORKERS", min(8L, max(1L, parallel::detectCores() - 1L))))
ARM_SCALE <- 0.3                 # age_sicker: the control's share of the ventilated age gradient
POWER_GRID <- c(0.25, 0.5, 0.75, 1, 1.5)
N_SITES <- 4L                    # multisite
SITE_EFFECT_CV <- 0.3            # between-site SD of the strain effect, as a share of it
N_REPS_MULTISITE  <- max(N_REPS %/% 4L, min(N_REPS, 100L))   # each replicate is N_SITES sites
N_REPS_CROSSCHECK <- min(100L, N_REPS)
N_REALISM <- 20L                 # simulated sites averaged for the realism table
TESTBED_SEED <- 20261111L
ALPHA <- 0.05
RESUME <- Sys.getenv("PBWPFVC_CAUSAL_RESUME", "0") == "1"
CHECKPOINT_DIR <- file.path(final_dir, "causal_testbed_checkpoints")
dir.create(CHECKPOINT_DIR, showWarnings = FALSE, recursive = TRUE)
SCENARIO_STAGES <- c("matched", "realistic", "stress", "power")
z_crit <- qnorm(1 - ALPHA / 2)
# the generator makes a log-scale lab; the vasopressor and OSI outcomes are not simulated
stopifnot(all(MARKERS %in% CAUSAL_MARKERS$marker[CAUSAL_MARKERS$source == "lab"]),
          all(STAGES %in% c("matched", "realistic", "stress", "power", "timecourse", "illness_death", "multisite", "crosscheck")))
TIMECOURSE_DESIGNS <- c("unadjusted", "control_did", "free_age_did", "within_onoff", "off_reversal")

patients <- as.data.frame(load_causal_patients(config))
# every design with the index VT/PBW held and with it dropped
designs_both <- list(held = causal_designs(TRUE), dropped = causal_designs(FALSE))
designs <- designs_both$held

# ---- The scenarios -------------------------------------------------------------------
# One row per scenario: its stage, whether a strain effect is planted, and its
# parameters. beta is signed (negative for platelets).
scenario_table <- function(beta, control_age_gradient) {
  d <- function(index) confounding_equal_to(beta, patients, index)
  truths <- function(base) list(
    null            = scenario(base = base),
    strain          = scenario(beta = beta, base = base),
    age_everywhere  = scenario(d_all = d("age"), base = base),
    age_sicker      = scenario(d_all = d("age"), arm_scale = ARM_SCALE, base = base),
    age_acute       = scenario(d_acute = d("age"), base = base),
    age_ventilation = scenario(d_vent = d("age"), base = base),
    sex_ventilation = scenario(index = "sex", d_vent = d("sex"), base = base),
    strain_plus_age = scenario(beta = beta, d_all = d("age"), base = base))
  realistic_more <- list(
    composite_ventilation = scenario(index = "composite", d_vent = d("composite")),
    gli_ventilation       = scenario(index = "gli", d_vent = d("gli")),
    age65_everywhere      = scenario(index = "age_hinge", d_all = d("age_hinge")),
    age_everywhere_half   = scenario(d_all = 0.5 * d("age")),
    age_everywhere_double = scenario(d_all = 2 * d("age")),
    # the age gradient of the marker's slope observed in the site's own controls
    age_as_in_controls    = scenario(d_all = control_age_gradient))
  # one departure from the realistic base at a time; `null` says whether the scenario
  # is also run without a strain effect (those that can move a null result)
  stress <- list(
    dropout_strong     = list(over = list(info_off = 1, info_death = 1.5), null = TRUE),
    measurement_strong = list(over = list(info_measure = 1), null = TRUE),
    heavy_tails        = list(over = list(error_df = 3), null = TRUE),
    autocorrelation    = list(over = list(ar_share = 0.8, ar_rho = 0.9), null = TRUE),
    vt_by_indication   = list(over = list(vt_model = "indication", u_slope = 1.5), null = TRUE),
    vt_gate_collider   = list(over = list(vt_model = "indication", gate = TRUE, u_slope = 1.5), null = TRUE),
    dose_threshold     = list(over = list(dose_shape = "threshold"), null = FALSE),
    time_saturating    = list(over = list(time_shape = "saturating"), null = FALSE),
    time_delayed       = list(over = list(time_shape = "delayed"), null = FALSE),
    effect_differs_pre = list(over = list(het_pre = 0.3), null = FALSE),
    pressure_before_t0 = list(over = list(pre_exposure = 0.4), null = FALSE),
    height_error       = list(over = list(height_error_cm = 3), null = FALSE),
    lung_size_varies   = list(over = list(sd_lung = 0.25), null = FALSE),
    denominator_age25  = list(over = list(denominator = "age25"), null = FALSE),
    denominator_norace = list(over = list(denominator = "norace"), null = FALSE))
  rows <- list()
  add <- function(stage, name, sc, effect = sc$beta != 0, multiple = NA_real_)
    rows[[length(rows) + 1]] <<- tibble(stage = stage, scenario = name, has_effect = effect, effect_multiple = multiple, params = list(sc))
  if ("matched" %in% STAGES) for (nm in names(truths("matched"))) add("matched", nm, truths("matched")[[nm]])
  if ("realistic" %in% STAGES) {
    for (nm in names(truths("realistic"))) add("realistic", nm, truths("realistic")[[nm]])
    for (nm in names(realistic_more)) add("realistic", nm, realistic_more[[nm]])
  }
  if ("stress" %in% STAGES) {
    for (nm in names(stress)) {
      over <- stress[[nm]]$over
      if (stress[[nm]]$null) add("stress", paste0(nm, " | null"), do.call(scenario, over))
      add("stress", paste0(nm, " | strain"), do.call(scenario, c(list(beta = beta), over)))
    }
    add("stress", "dropout_strong | age_everywhere", scenario(d_all = d("age"), info_off = 1, info_death = 1.5))
  }
  if ("power" %in% STAGES) for (k in POWER_GRID)
    add("power", paste0("strain x ", k), scenario(beta = k * beta), multiple = k)
  bind_rows(rows) %>% mutate(cell = row_number(), .before = 1)
}

# the marker's age gradient in the site's controls, per ventilated SD of age: a
# confounder of an observed size, not one sized to the hypothesis
control_age_gradient <- function(panel) {
  d <- tibble(design = "control_age", arms = list("no support"), phases = list("on"),
              formula = "~ age10 + sofa_z + sf_z", target = "age10", sign = 1, fails_if = "")
  prep <- prepare_design(d, patients, panel)
  if (!prep$ok) return(0)
  fit_design(prep, align_outcome(prep, panel))$estimate * sd(patients$age10[patients$vent == 1])
}

# ---- One replicate --------------------------------------------------------------------
# simulate a site, fit every design, apply the decision rule
# (both = FALSE: VT/PBW held only, no rule)
one_replicate <- function(sc, calib, sign, seed, pats = patients, both = TRUE) {
  sim <- simulate_site(pats, calib, sc, sign, seed)
  held <- fit_all_designs(designs_both$held, sim$patients, sim$panel)
  fits <- data.frame(vtpbw = "held", design = held$design, estimate = held$estimate, se = held$se, truth = sim$truth)
  if (!both) return(list(fits = fits))
  dropped <- fit_all_designs(designs_both$dropped, sim$patients, sim$panel)
  list(fits = rbind(fits, data.frame(vtpbw = "dropped", design = dropped$design, estimate = dropped$estimate,
                                     se = dropped$se, truth = sim$truth_dropped)),
       rule = triangulation_supports_both(held, dropped, sign))
}
# the time course of one simulated site: each segment, the cumulative effect, and the
# straight line over the whole window beside them
one_timecourse_replicate <- function(sc, calib, sign, seed) {
  sim <- simulate_site(patients, calib, sc, sign, seed)
  by_day <- sim$truth_by_day; before <- c(0, head(by_day, -1)); days <- as.numeric(names(by_day))
  do.call(rbind, lapply(TIMECOURSE_DESIGNS, function(dn) {
    design <- designs[designs$design == dn, ]
    tc <- fit_time_course(design, sim$patients, sim$panel)
    k <- match(tc$to, days)
    tc$truth <- ifelse(tc$quantity == "cumulative", by_day[k], (by_day[k] - before[k]) / (tc$to - tc$from))
    line <- fit_all_designs(design, sim$patients, sim$panel)
    rbind(data.frame(design = dn, quantity = tc$quantity, window = tc$window, estimate = tc$estimate, se = tc$se, truth = tc$truth),
          data.frame(design = dn, quantity = "straight line", window = sprintf("days 0 to %g", ON_DAYS),
                     estimate = line$estimate, se = line$se, truth = sim$truth))
  }))
}
# the illness-death analysis of one simulated site. The marker at t0 is the baseline;
# death is followed over the simulated window.
ILLNESS_DEATH_DESIGNS <- c("unadjusted", "control_did", "free_age_did")
one_illness_death_replicate <- function(sc, calib, sign, seed, marker) {
  sim <- simulate_site(patients, calib, sc, sign, seed)
  values <- rbind(data.frame(id = sim$patients$id, t = 0, value = exp(sim$patients$level0)),
                  data.frame(id = sim$panel$id, t = sim$panel$t, value = exp(sim$panel$y)))
  injury <- load_injury_times(NULL, sim$patients, marker, values = values)
  states <- illness_death_states(sim$patients, injury, horizon = ON_DAYS)
  out <- fit_illness_death(designs[designs$design %in% ILLNESS_DEATH_DESIGNS, ], states)
  data.frame(design = out$design, transition = out$transition, estimate = out$estimate, se = out$se, n_events = out$n_events)
}
# A site of a multisite replicate: the site's patients resampled, its own effect,
# confounder size and noise
one_multisite_replicate <- function(sc, calib, sign, seed) {
  set.seed(seed)
  site_fits <- lapply(seq_len(N_SITES), function(k) {
    pats <- patients[sample.int(nrow(patients), replace = TRUE), ]
    pats$id <- sprintf("s%05d", seq_len(nrow(pats)))
    sck <- sc
    sck$beta  <- sc$beta * (1 + SITE_EFFECT_CV * rnorm(1))
    sck$d_all <- sc$d_all * exp(rnorm(1, 0, 0.3))
    cal <- calib; cal$scale$resid <- cal$scale$resid * runif(1, 0.8, 1.25)
    out <- one_replicate(sck, cal, sign, seed * 10L + k, pats, both = FALSE)$fits
    out$site <- k; out
  })
  d <- do.call(rbind, site_fits); d <- d[!is.na(d$estimate), ]
  do.call(rbind, lapply(split(d, d$design), function(x) {
    w <- 1 / x$se^2; fixed <- sum(w * x$estimate) / sum(w); k <- nrow(x)
    q <- sum(w * (x$estimate - fixed)^2)
    tau2 <- if (k > 1) max((q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)), 0) else 0      # DerSimonian-Laird
    wr <- 1 / (x$se^2 + tau2)
    data.frame(design = x$design[1], n_sites = k, fixed = fixed, fixed_se = sqrt(1 / sum(w)),
               random = sum(wr * x$estimate) / sum(wr), random_se = sqrt(1 / sum(wr)))
  }))
}
# the same simulated site by the within estimator and by a linear mixed model
CROSSCHECK_DESIGNS <- c("unadjusted", "control_did", "free_age_did", "within_onoff")
one_crosscheck_replicate <- function(sc, calib, sign, seed) {
  sim <- simulate_site(patients, calib, sc, sign, seed)
  do.call(rbind, lapply(CROSSCHECK_DESIGNS, function(dn) {
    design <- designs[designs$design == dn, ]
    prep <- prepare_design(design, sim$patients, sim$panel)
    if (!prep$ok) return(NULL)
    within <- fit_design(prep, align_outcome(prep, sim$panel))
    d <- cbind(sim$panel[prep$rows, ], sim$patients[match(sim$panel$id[prep$rows], sim$patients$id), setdiff(names(sim$patients), "id")])
    d$on <- as.numeric(d$phase == "on"); d$lv_on <- d$lv0 * d$on
    # levels and slopes both follow the design's covariates; a random level and slope per patient
    form <- as.formula(paste0("y ~ (", sub("^~", "", design$formula), ") * t"))
    mixed <- tryCatch(nlme::lme(form, random = ~ t | id, data = d, method = "REML",
                                control = nlme::lmeControl(opt = "optim", returnObject = TRUE)), error = function(e) NULL)
    est <- NA_real_; se <- NA_real_
    if (!is.null(mixed)) {
      tab <- summary(mixed)$tTable
      hit <- which(sorted_term(rownames(tab)) == sorted_term(paste0(design$target, ":t")))
      if (length(hit) == 1) { est <- tab[hit, "Value"]; se <- tab[hit, "Std.Error"] }
    }
    data.frame(design = dn, within = within$estimate, within_se = within$se, mixed = est, mixed_se = se, truth = sim$truth)
  }))
}

# ---- Parallel workers -------------------------------------------------------------------
run_jobs <- function(jobs, fun, ...) {
  if (is.null(cluster)) lapply(jobs, fun, ...) else parallel::parLapplyLB(cluster, jobs, fun, ..., chunk.size = 10)
}
# replicate r of cell c always has the same seed, whatever the number of workers
seed_for <- function(cell, rep) TESTBED_SEED + 100000L * cell + rep

# Monte Carlo summaries of one design in one scenario
performance <- function(estimate, se, truth) {
  ok <- !is.na(estimate) & !is.na(se); n <- sum(ok)
  if (n < 2) return(tibble(n_reps = n))
  e <- estimate[ok]; s <- se[ok]; tr <- truth[ok]
  covered <- abs(e - tr) <= z_crit * s; rejected <- abs(e / s) > z_crit
  emp <- sd(e); model <- sqrt(mean(s^2))
  tibble(n_reps = n, truth = mean(tr), mean_estimate = mean(e), bias = mean(e - tr), bias_mcse = sd(e - tr) / sqrt(n),
         empirical_se = emp, empirical_se_mcse = emp / sqrt(2 * (n - 1)), model_se = model,
         # the model SE relative to the empirical one: below 1, the intervals are too narrow
         se_ratio = model / emp, coverage = mean(covered), coverage_mcse = sqrt(mean(covered) * (1 - mean(covered)) / n),
         rejection_rate = mean(rejected), rejection_mcse = sqrt(mean(rejected) * (1 - mean(rejected)) / n),
         mde_80 = (z_crit + qnorm(0.80)) * emp)
}

# Separate R processes (not forks), so the same code runs on Windows, macOS and Linux.
cluster <- NULL
if (WORKERS > 1) {
  # one thread per process: a multithreaded BLAS in every worker makes them fight for
  # the same cores and run many times slower than one process alone
  thread_vars <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
  old_threads <- Sys.getenv(thread_vars, unset = NA)
  do.call(Sys.setenv, as.list(setNames(rep("1", length(thread_vars)), thread_vars)))
  cluster <- parallel::makePSOCKcluster(WORKERS)
  for (v in thread_vars) if (is.na(old_threads[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, setNames(list(old_threads[[v]]), v))
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  invisible(parallel::clusterCall(cluster, function(paths) .libPaths(paths), .libPaths()))
  invisible(parallel::clusterEvalQ(cluster, suppressPackageStartupMessages({ library(dplyr); library(tibble) })))
  parallel::clusterExport(cluster, setdiff(ls(globalenv()), c("cluster", "config")))
}

# A stage's result for one marker: read from the checkpoint when resuming and the
# settings match, otherwise computed (value is evaluated only then) and saved.
checkpoint <- function(marker, stage, value) {
  path <- file.path(CHECKPOINT_DIR, paste0(marker, "_", stage, ".rds"))
  settings <- list(reps = N_REPS, effect = EFFECT, seed = TESTBED_SEED, sites = N_SITES, site_cv = SITE_EFFECT_CV,
                   stages = if (stage == "scenarios") sort(intersect(STAGES, SCENARIO_STAGES)))
  if (RESUME && file.exists(path)) {
    saved <- readRDS(path)
    if (identical(saved$settings, settings)) { message("  ", stage, ": reused from ", path); return(saved$value) }
    message("  ", stage, ": the saved result has other settings; running it again")
  }
  saveRDS(list(settings = settings, value = value), path)
  value
}

# ---- Run ----------------------------------------------------------------------------------
operating <- list(); decision <- list(); scenario_log <- list(); realism <- list(); multisite <- list(); crosscheck <- list()
timecourse <- list(); illness_death <- list()
started <- Sys.time()
for (m in which(CAUSAL_MARKERS$marker %in% MARKERS)) {
  marker <- CAUSAL_MARKERS$marker[m]; sign <- CAUSAL_MARKERS$injury_sign[m]
  panel_raw <- as.data.frame(load_causal_panel(config, patients, CAUSAL_MARKERS$lab_category[m], complete = FALSE))
  panel <- complete_phases(panel_raw)
  calib0 <- calibrate_site(patients, panel_raw, panel)
  calib <- list(matched   = tune_generator(patients, calib0, "matched", sign, TESTBED_SEED),
                realistic = tune_generator(patients, calib0, "realistic", sign, TESTBED_SEED))
  cells <- scenario_table(sign * EFFECT, control_age_gradient(panel))
  message("\n", marker, ": generator tuned after ", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1), " min; ",
          nrow(cells), " scenarios x ", N_REPS, " replicates on ", WORKERS, " process(es)")
  if (!is.null(cluster)) parallel::clusterExport(cluster, c("calib", "cells", "sign"), envir = environment())

  # --- realism: what the analysis reads off the real site and off simulated ones ---
  observed <- panel_summaries(panel, patients)
  for (base in names(calib)) {
    sims <- lapply(seq_len(N_REALISM), function(r) {
      sim <- simulate_site(patients, calib[[base]], scenario(base = base), sign, TESTBED_SEED + r)
      rbind(panel_summaries(sim$panel, sim$patients),
            data.frame(arm = c("ventilated", "ventilated"), phase = "on",
                       statistic = c("liberated within the window", "died within the window"), value = unname(sim$ended)))
    })
    expo <- phase_exposure(patients); v <- patients$vent == 1
    observed_all <- rbind(observed, data.frame(arm = "ventilated", phase = "on",
                                               statistic = c("liberated within the window", "died within the window"),
                                               value = c(mean(expo$ended[v] == "off"), mean(expo$ended[v] == "death"))))
    realism[[paste(marker, base)]] <- bind_rows(sims) %>% group_by(arm, phase, statistic) %>%
      summarise(simulated = mean(value, na.rm = TRUE), simulated_sd = sd(value, na.rm = TRUE), .groups = "drop") %>%
      inner_join(observed_all %>% rename(observed = value), by = c("arm", "phase", "statistic")) %>%
      mutate(marker = marker, generator = base, .before = 1)
  }

  if (nrow(cells)) {
  # --- the scenarios ---
  saved <- checkpoint(marker, "scenarios", {
  jobs <- expand.grid(rep = seq_len(N_REPS), cell = cells$cell)
  results <- run_jobs(split(jobs, seq_len(nrow(jobs))), function(job) {
    sc <- cells$params[[job$cell]]
    out <- one_replicate(sc, calib[[sc$base]], sign, seed_for(job$cell, job$rep))
    out$fits$cell <- job$cell; out$rule <- c(cell = job$cell, out$rule); out
  })
  fits <- do.call(rbind, lapply(results, `[[`, "fits"))
  rules <- as.data.frame(do.call(rbind, lapply(results, `[[`, "rule")))
  list(operating = fits %>% group_by(cell, vtpbw, design) %>%
    group_modify(~ performance(.x$estimate, .x$se, .x$truth)) %>% ungroup() %>%
    left_join(fits %>% group_by(cell, vtpbw, design) %>% summarise(n_not_estimable = sum(is.na(estimate)), .groups = "drop"), by = c("cell", "vtpbw", "design")) %>%
    left_join(cells %>% select(cell, stage, scenario, has_effect, effect_multiple), by = "cell") %>% mutate(marker = marker),
  decision = rules %>% pivot_longer(-cell, names_to = "rule", values_to = "supports") %>%
    group_by(cell, rule) %>% summarise(n_reps = n(), supports_rate = mean(supports),
                                       supports_mcse = sqrt(supports_rate * (1 - supports_rate) / n_reps), .groups = "drop") %>%
    left_join(cells %>% select(cell, stage, scenario, has_effect, effect_multiple), by = "cell") %>% mutate(marker = marker),
  scenario_log = cells %>% mutate(marker = marker, .before = 1) %>%
    bind_cols(bind_rows(lapply(cells$params, as_tibble))) %>% select(-params))
  })
  operating[[marker]] <- saved$operating; decision[[marker]] <- saved$decision; scenario_log[[marker]] <- saved$scenario_log
  message("  scenarios done after ", round(as.numeric(difftime(Sys.time(), started, units = "mins"))), " min")
  }

  # --- the time course: is the segment-by-segment estimate right, and what does it cost? ---
  if ("timecourse" %in% STAGES) {
    timecourse[[marker]] <- checkpoint(marker, "timecourse", {
    truths <- list(null = scenario(), linear = scenario(beta = sign * EFFECT),
                   saturating = scenario(beta = sign * EFFECT, time_shape = "saturating"),
                   delayed = scenario(beta = sign * EFFECT, time_shape = "delayed"))
    if (!is.null(cluster)) parallel::clusterExport(cluster, "truths", envir = environment())
    jobs <- expand.grid(rep = seq_len(N_REPS), truth = seq_along(truths))
    courses <- run_jobs(split(jobs, seq_len(nrow(jobs))), function(job) {
      out <- one_timecourse_replicate(truths[[job$truth]], calib$realistic, sign, 7000000L + 10000L * job$truth + job$rep)
      out$scenario <- names(truths)[job$truth]; out
    })
    bind_rows(courses) %>%
      mutate(window = factor(window, levels = unique(window))) %>%
      group_by(scenario, design, quantity, window) %>% group_modify(~ performance(.x$estimate, .x$se, .x$truth)) %>% ungroup() %>%
      mutate(marker = marker, .before = 1)
    })
    message("  timecourse done after ", round(as.numeric(difftime(Sys.time(), started, units = "mins"))), " min")
  }

  # --- injury and death as states ---
  if ("illness_death" %in% STAGES) {
    illness_death[[marker]] <- checkpoint(marker, "illness_death", {
    age <- confounding_equal_to(sign * EFFECT, patients, "age")
    truths <- list(null = scenario(),
                   # death follows the marker (the realistic base), with strain at 1x and 3x the planted effect
                   strain = scenario(beta = sign * EFFECT), strain_x3 = scenario(beta = 3 * sign * EFFECT),
                   # and strain also kills by a route the marker does not see
                   strain_x3_direct_death = scenario(beta = 3 * sign * EFFECT, death_direct = 3),
                   # no strain: age raises death in both cohorts; age steepens the marker in both
                   frailty_death = scenario(death_age = 0.5), age_on_marker = scenario(d_all = 3 * age))
    if (!is.null(cluster)) parallel::clusterExport(cluster, c("truths", "marker"), envir = environment())
    jobs <- expand.grid(rep = seq_len(N_REPS), truth = seq_along(truths))
    states <- run_jobs(split(jobs, seq_len(nrow(jobs))), function(job) {
      out <- one_illness_death_replicate(truths[[job$truth]], calib$realistic, sign, 8000000L + 10000L * job$truth + job$rep, marker)
      out$scenario <- names(truths)[job$truth]; out
    })
    bind_rows(states) %>%
      mutate(scenario = factor(scenario, levels = names(truths)), transition = factor(transition, levels = unique(transition))) %>%
      group_by(scenario, design, transition) %>%
      summarise(n_reps = sum(!is.na(estimate)), n_not_estimable = sum(is.na(estimate)), mean_events = mean(n_events, na.rm = TRUE),
                mean_estimate = mean(estimate, na.rm = TRUE), mean_mcse = sd(estimate, na.rm = TRUE) / sqrt(n_reps),
                empirical_se = sd(estimate, na.rm = TRUE), model_se = sqrt(mean(se^2, na.rm = TRUE)), se_ratio = model_se / empirical_se,
                rejection_rate = mean(abs(estimate / se) > z_crit, na.rm = TRUE),
                rejection_mcse = sqrt(rejection_rate * (1 - rejection_rate) / n_reps),
                positive_rate = mean(estimate / se > z_crit, na.rm = TRUE), .groups = "drop") %>%
      mutate(marker = marker, .before = 1)
    })
    message("  illness_death done after ", round(as.numeric(difftime(Sys.time(), started, units = "mins"))), " min")
  }

  # --- multisite ---
  if ("multisite" %in% STAGES) {
    multisite[[marker]] <- checkpoint(marker, "multisite", {
    truths <- list(null = scenario(), strain = scenario(beta = sign * EFFECT),
                   strain_plus_age = scenario(beta = sign * EFFECT, d_all = confounding_equal_to(sign * EFFECT, patients, "age")),
                   age_sicker = scenario(d_all = confounding_equal_to(sign * EFFECT, patients, "age"), arm_scale = ARM_SCALE))
    if (!is.null(cluster)) parallel::clusterExport(cluster, "truths", envir = environment())
    jobs <- expand.grid(rep = seq_len(N_REPS_MULTISITE), truth = seq_along(truths))
    pooled <- run_jobs(split(jobs, seq_len(nrow(jobs))), function(job) {
      out <- one_multisite_replicate(truths[[job$truth]], calib$realistic, sign, 5000000L + 10000L * job$truth + job$rep)
      out$scenario <- names(truths)[job$truth]; out$truth <- truths[[job$truth]]$beta; out
    })
    bind_rows(pooled) %>%
      pivot_longer(c(fixed, random), names_to = "pooling", values_to = "estimate") %>%
      mutate(se = if_else(pooling == "fixed", fixed_se, random_se),
             pooling = if_else(pooling == "fixed", "common effect", "random effects (DerSimonian-Laird)")) %>%
      group_by(scenario, design, pooling) %>% group_modify(~ performance(.x$estimate, .x$se, .x$truth)) %>% ungroup() %>%
      mutate(marker = marker, n_sites = N_SITES, site_effect_cv = SITE_EFFECT_CV)
    })
    message("  multisite done after ", round(as.numeric(difftime(Sys.time(), started, units = "mins"))), " min")
  }

  # --- the within estimator against a mixed model ---
  if ("crosscheck" %in% STAGES) {
    crosscheck[[marker]] <- checkpoint(marker, "crosscheck", {
    truths <- list(null = scenario(), strain = scenario(beta = sign * EFFECT))
    if (!is.null(cluster)) parallel::clusterExport(cluster, "truths", envir = environment())
    jobs <- expand.grid(rep = seq_len(N_REPS_CROSSCHECK), truth = seq_along(truths))
    pairs <- run_jobs(split(jobs, seq_len(nrow(jobs))), function(job) {
      out <- one_crosscheck_replicate(truths[[job$truth]], calib$realistic, sign, 9000000L + 10000L * job$truth + job$rep)
      if (!is.null(out)) out$scenario <- names(truths)[job$truth]
      out
    })
    bind_rows(pairs) %>% group_by(scenario, design) %>%
      summarise(n_reps = n(), mixed_failed = sum(is.na(mixed)), truth = mean(truth),
                within_mean = mean(within[!is.na(mixed)]), mixed_mean = mean(mixed, na.rm = TRUE),
                mean_difference = mean(mixed - within, na.rm = TRUE),
                difference_mcse = sd(mixed - within, na.rm = TRUE) / sqrt(sum(!is.na(mixed))),
                correlation = cor(mixed, within, use = "complete.obs"),
                within_se = sqrt(mean(within_se^2)), mixed_se = sqrt(mean(mixed_se^2, na.rm = TRUE)), .groups = "drop") %>%
      mutate(marker = marker, .before = 1)
    })
    message("  crosscheck done after ", round(as.numeric(difftime(Sys.time(), started, units = "mins"))), " min")
  }
}

# ---- Tables ---------------------------------------------------------------------------------
tag <- function(d) d %>% mutate(site = site_name, cohort = if (config$vtpbw_gate == "0") "ungated" else "gated", .before = 1)
# the scenario tables and figures, when a scenario stage was run (a run of the other
# stages alone leaves the earlier tables in place)
has_scenarios <- length(operating) > 0
if (has_scenarios) {
operating <- bind_rows(operating) %>%
  mutate(planted_effect = EFFECT,
         # A design's verdict on a scenario.
         #   biased        off by more than a quarter of the planted effect and by more
         #                 than three Monte Carlo standard errors
         #   undercovers   not biased, but its 95% interval covers the truth less than 90%
         #                 of the time (beyond Monte Carlo error)
         #   underpowered  neither, but this site would detect the planted effect less
         #                 than 80% of the time: its null is not evidence against strain
         verdict = case_when(is.na(bias) ~ "not estimable",
                             abs(bias) > 0.25 * EFFECT & abs(bias) > 3 * bias_mcse ~ "biased",
                             coverage + 3 * coverage_mcse < 0.90 ~ "undercovers",
                             mde_80 > EFFECT ~ "underpowered",
                             TRUE ~ "ok")) %>%
  left_join(designs %>% select(design, fails_if), by = "design") %>%
  select(marker, stage, scenario, has_effect, effect_multiple, vtpbw, design, truth, mean_estimate, bias, bias_mcse,
         empirical_se, empirical_se_mcse, model_se, se_ratio, coverage, coverage_mcse, rejection_rate, rejection_mcse,
         mde_80, verdict, n_reps, n_not_estimable, planted_effect, fails_if) %>% tag()
decision <- bind_rows(decision) %>%
  mutate(rule = recode(rule, rule_both = "decision rule", rule_held = "rule, VT/PBW held only", rule_dropped = "rule, VT/PBW dropped only",
                       unadjusted_alone = "unadjusted model alone", control_did_alone = "control contrast alone"),
         # with an effect, "supports" is a correct call; without one, a false positive
         reading = if_else(has_effect, "power", "false positive")) %>%
  select(marker, stage, scenario, has_effect, effect_multiple, rule, reading, supports_rate, supports_mcse, n_reps) %>% tag()
scenario_log <- bind_rows(scenario_log) %>% mutate(planted_effect = EFFECT, n_reps = N_REPS) %>% tag()
realism <- bind_rows(realism) %>% tag()

write_csv(operating, file.path(final_dir, paste0("causal_testbed_operating_", site_name, ".csv")))
write_csv(decision, file.path(final_dir, paste0("causal_testbed_decision_", site_name, ".csv")))
write_csv(scenario_log, file.path(final_dir, paste0("causal_testbed_scenarios_", site_name, ".csv")))
write_csv(realism, file.path(final_dir, paste0("causal_testbed_realism_", site_name, ".csv")))
if ("power" %in% STAGES)
  write_csv(operating %>% filter(stage == "power", vtpbw == "held") %>%
              select(site, cohort, marker, design, effect_multiple, truth, rejection_rate, rejection_mcse, bias, bias_mcse, n_reps) %>%
              bind_rows(decision %>% filter(stage == "power") %>%
                          transmute(site, cohort, marker, design = rule, effect_multiple, rejection_rate = supports_rate,
                                    rejection_mcse = supports_mcse, n_reps)),
            file.path(final_dir, paste0("causal_testbed_power_", site_name, ".csv")))
}
if (length(multisite)) write_csv(bind_rows(multisite) %>% tag(), file.path(final_dir, paste0("causal_testbed_multisite_", site_name, ".csv")))
if (length(timecourse)) write_csv(bind_rows(timecourse) %>% tag(), file.path(final_dir, paste0("causal_testbed_timecourse_", site_name, ".csv")))
if (length(illness_death)) write_csv(bind_rows(illness_death) %>% tag(), file.path(final_dir, paste0("causal_testbed_illness_death_", site_name, ".csv")))
if (length(crosscheck)) write_csv(bind_rows(crosscheck) %>% tag(), file.path(final_dir, paste0("causal_testbed_crosscheck_", site_name, ".csv")))

if (has_scenarios) {
# ---- Figures -----------------------------------------------------------------------------------
verdict_colours <- c(ok = "#cfe8d5", underpowered = "#f3e9c0", undercovers = "#f6d6b0", biased = "#f1c9c5", "not estimable" = "grey90")
pdf(file.path(final_dir, paste0("causal_testbed_", site_name, ".pdf")), width = 14, height = 8.5)
for (st in intersect(c("matched", "realistic", "stress"), STAGES)) {
  d <- operating %>% filter(stage == st, vtpbw == "held") %>%
    mutate(design = factor(design, levels = rev(designs$design)), scenario = factor(scenario, levels = unique(scenario)),
           label = if_else(is.na(bias), "", sprintf("%.2f\n(%+.2f)", rejection_rate, bias / EFFECT)))
  print(ggplot(d, aes(scenario, design, fill = verdict)) + geom_tile(colour = "white") +
          geom_text(aes(label = label), size = 2, lineheight = 0.85) + facet_wrap(~ marker, ncol = 1) +
          scale_fill_manual(values = verdict_colours, name = NULL) +
          labs(x = NULL, y = NULL, title = paste0("Stage: ", st, ", VT/PBW held (", site_name, ", ", N_REPS, " replicates)"),
               subtitle = paste0("Rejection rate at 5%, and (bias as a share of the planted effect, ", EFFECT, " log units/day per log unit)")) +
          theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom", panel.grid = element_blank()))
}
d <- decision %>% filter(stage != "power") %>% mutate(scenario = factor(paste(stage, scenario, sep = ": "), levels = rev(unique(paste(stage, scenario, sep = ": ")))))
print(ggplot(d, aes(supports_rate, scenario, colour = rule, shape = rule)) +
        geom_vline(xintercept = ALPHA, linetype = 2, colour = "grey50") +
        geom_pointrange(aes(xmin = supports_rate - 1.96 * supports_mcse, xmax = supports_rate + 1.96 * supports_mcse), size = 0.25) +
        facet_grid(reading ~ marker, scales = "free_y", space = "free_y") + scale_colour_brewer(palette = "Dark2", name = NULL) + scale_shape(name = NULL) +
        labs(x = "Share of simulated sites called \"supports the hypothesis\" (95% Monte Carlo interval)", y = NULL,
             title = "The decision rule against single analyses",
             subtitle = "The decision rule requires the rule to hold with VT/PBW held and dropped. Top: scenarios with no strain effect (every call is a false positive). Bottom: scenarios with one (power).") +
        theme_bw(base_size = 9) + theme(legend.position = "bottom"))
if ("power" %in% STAGES) {
  d <- bind_rows(operating %>% filter(stage == "power", vtpbw == "held") %>% transmute(marker, line = design, effect_multiple, rate = rejection_rate, mcse = rejection_mcse),
                 decision %>% filter(stage == "power", rule == "decision rule") %>% transmute(marker, line = rule, effect_multiple, rate = supports_rate, mcse = supports_mcse)) %>%
    filter(line %in% c("decision rule", "control_did", "within_onoff", "within_triple", "off_reversal", "free_age_did", "adjusted"))
  print(ggplot(d, aes(effect_multiple * EFFECT, rate, colour = line)) + geom_line() + geom_point(size = 1) +
          geom_errorbar(aes(ymin = rate - 1.96 * mcse, ymax = rate + 1.96 * mcse), width = 0) +
          geom_hline(yintercept = 0.8, linetype = 2, colour = "grey50") + facet_wrap(~ marker) +
          scale_colour_brewer(palette = "Dark2", name = NULL) +
          labs(x = "Planted strain effect (log units per day per log unit of VT/PFVC)", y = "Rejection rate (decision rule: share called \"supports\")",
               title = "Power by effect size, realistic generator") + theme_bw(base_size = 10) + theme(legend.position = "bottom"))
}
d <- realism %>% filter(generator == "realistic", is.finite(observed), is.finite(simulated))
print(ggplot(d, aes(observed, simulated, colour = paste(arm, phase))) + geom_abline(linetype = 2, colour = "grey50") +
        geom_pointrange(aes(ymin = simulated - 2 * simulated_sd, ymax = simulated + 2 * simulated_sd), size = 0.3) +
        facet_wrap(~ paste(marker, statistic, sep = ": "), scales = "free", ncol = 5) + scale_colour_brewer(palette = "Dark2", name = NULL) +
        labs(title = "Realism: what the analysis reads off the real site and off simulated ones (realistic generator, no effect)",
             x = "Observed", y = paste0("Simulated (mean of ", N_REALISM, " sites, +/- 2 SD)")) +
        theme_bw(base_size = 8) + theme(legend.position = "bottom"))
invisible(dev.off())
}

# ---- Console summary -------------------------------------------------------------------------------
if (has_scenarios) {
message("\nDecision rule: share of simulated sites called \"supports\" (MCSE)")
print(as.data.frame(decision %>% filter(stage != "power") %>%
                      transmute(marker, stage, scenario, reading, rule, cell = sprintf("%.3f (%.3f)", supports_rate, supports_mcse)) %>%
                      pivot_wider(names_from = rule, values_from = cell)), row.names = FALSE)
message("\nDesigns that are biased or undercover, by stage:")
print(as.data.frame(operating %>% filter(verdict %in% c("biased", "undercovers"), stage != "power") %>%
                      count(marker, stage, vtpbw, design, verdict) %>% pivot_wider(names_from = stage, values_from = n, values_fill = 0)), row.names = FALSE)
}
message("\ncausal_testbed complete after ", round(as.numeric(difftime(Sys.time(), started, units = "mins"))), " min: ", final_dir)
