# =============================================================================
# Script 06: Target trial emulation — optimal mechanical-metric thresholds
#            and PBW-vs-PFVC normalization
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# FEATURE BRANCH (feature/target-trial-emulation): NOT part of 00_run_pipeline.R.
# Standalone, in the explore_*.R tradition — reads the script-03 cross-sectional
# dataset and refits its own models.
#
# WHY THIS SCRIPT EXISTS -------------------------------------------------------
# A direct VT/PFVC-vs-VT/PBW strategy comparison fails on positivity: VT/PFVC is
# nearly deterministic in demographics (height/age/sex/race via GLI-2012), so
# conditioning on demographics collapses the exposure contrast. The mechanical
# metrics (driving pressure, mechanical power, elastance) carry variation driven
# by lung pathology, not body size, so a threshold-strategy target trial
# emulation is feasible there. This script:
#   (A) IPW threshold emulation  — grid-of-trials + spline + maxstat per metric
#                                   => the OPERATIONAL threshold.
#   (B) IV / natural experiment  — anthropometric mis-sizing (PBW/PFVC | demo) as
#                                   an instrument for mechanical stress, to address
#                                   UNMEASURED severity confounding (2SRI).
#   (C) Normalization-absorption — the correct lung-size scale is the one that
#                                   renders the instrument conditionally
#                                   independent of mortality (causal headline).
#
# LIMITATIONS (see header notes): achieved-not-assigned point exposure; IPW
# controls only MEASURED severity (causal weight carried by the IV arm + bounded
# by E-values); the IV assumes anthropometry affects mortality only via lung
# mechanics (exclusion restriction; reported with first-stage F + Conley bounds).
# =============================================================================

# Pin BLAS to a single thread in THIS process and (by inheritance) in every PSOCK
# worker spawned later. Without this, N workers each run all-core multithreaded
# BLAS, oversubscribing the CPU (N x cores threads) and slowing the run to a
# crawl. Must be set before any BLAS use and before the workers are launched.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")

library(tidyverse)
library(arrow)
library(here)
library(survival)
library(splines)
library(maxstat)
library(boot)
library(EValue)
library(patchwork)
library(parallel)

source("utils/config.R")
site_name <- config$site_name

output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

# Okabe-Ito (discrete); viridis used for continuous fills via scale_*_viridis_*.
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

# L'Ecuyer streams give reproducible bootstraps across the parallel metric jobs.
RNGkind("L'Ecuyer-CMRG")
set.seed(20260615)

# Parallelism: the 9 metrics are independent, so the per-metric battery (grid +
# maxstat/c-index/2SRI bootstraps x 2 instruments) runs across a PSOCK cluster
# (one worker per metric). Override the worker count with PBWPFVC_CORES; set it
# to 1 for a serial run with a clean traceback.
N_CORES <- suppressWarnings(as.integer(Sys.getenv("PBWPFVC_CORES", unset = NA)))
if (is.na(N_CORES)) N_CORES <- max(1L, detectCores() - 1L)

# -----------------------------------------------------------------------------
# Run-size knobs (smaller bootstrap on synthetic, where we only check plumbing).
# -----------------------------------------------------------------------------
HORIZON       <- 60      # days, primary 60-day all-cause mortality
N_GRID        <- 30      # candidate thresholds per metric (10th-90th pctile)
MIN_ARM       <- 50      # skip a threshold if either arm < MIN_ARM
MIN_CELL      <- 10      # never report any group smaller than this (CLIF rule)
is_synthetic  <- identical(site_name, "synthetic_clif")
# Bootstrap reps (override the real-data IV count via PBWPFVC_NBOOT). 500 is
# ample for percentile CIs and roughly halves the dominant cost (the 2SRI loop
# refits an lm + Cox per rep per instrument); maxstat/c-index are trimmed too.
N_BOOT        <- if (is_synthetic) 200L else
  suppressWarnings(as.integer(Sys.getenv("PBWPFVC_NBOOT", unset = "500")))
if (is.na(N_BOOT)) N_BOOT <- 500L
N_BOOT_MAXSTAT<- if (is_synthetic) 200L else 300L
N_BOOT_CSTAT  <- if (is_synthetic) 100L else 200L
N_BOOT_STRATUM<- if (is_synthetic) 100L else 300L  # 2SRI reps within each age tertile

message("Site: ", site_name, " | synthetic = ", is_synthetic,
        " | bootstrap reps = ", N_BOOT)

# =============================================================================
# 6a. Load cross-sectional cohort and derive the two missing DP-normalized metrics
# =============================================================================
cross_sectional <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet")) %>%
  mutate(
    # The "VT/PBW x Ers" and "VT/PFVC x Ers" metrics (= DP normalized by lung size).
    dp_pbw  = dp / pbw,
    dp_pfvc = dp / pfvc,
    # Per-10-unit covariates, matching the script-04 convention.
    age10 = age_at_admission / 10,
    sf10  = sf_ratio / 10
  )

message("Cohort: ", nrow(cross_sectional), " patients")

# =============================================================================
# 6b. Time-zero re-anchoring at the index ventilation timepoint
# =============================================================================
# The published surv_time is anchored at ADMISSION; we re-anchor at the index
# IMV timepoint (recorded_dttm) so follow-up starts when the exposure is
# ascertained (removes the small immortal-time gap).
#
# *** SYNTHETIC-ONLY MORTALITY WORKAROUND *************************************
# The synthetic CLIF mortality fields (death_dttm / death_day and everything
# derived) are KNOWN to be buggy (upstream fix in progress on
# AartikSarma/synthetic_clif). For the SYNTHETIC SITE ONLY we replace survival
# with draws from a long-tailed distribution typical of critically ill
# mechanically ventilated patients (~35% 60-day mortality; decedent time-to-death
# truncated-lognormal, median ~9 d, tail to 60 d), drawn INDEPENDENTLY of the
# exposures — so the machinery is exercised without manufacturing a fake effect.
# For ANY REAL SITE the dataset's own death_dttm is used unchanged.
# REMOVE this block once the synthetic-CLIF mortality fix lands.
# ****************************************************************************
rtrunc_lnorm <- function(n_needed, meanlog, sdlog, lo, hi) {
  acc <- numeric(0)
  while (length(acc) < n_needed) {
    cand <- rlnorm(max(n_needed * 2L, 1000L), meanlog, sdlog)
    cand <- cand[cand > lo & cand <= hi]
    acc  <- c(acc, cand)
  }
  acc[seq_len(n_needed)]
}

if (is_synthetic) {
  message("*** SYNTHETIC SITE: replacing buggy mortality with simulated ",
          "long-tailed survival (see header). ***")
  set.seed(20260615)
  n <- nrow(cross_sectional)
  died60 <- rbinom(n, 1L, 0.35)                       # ~35% 60-day mortality
  tte    <- rep(NA_real_, n)
  n_dec  <- sum(died60 == 1L)
  tte[died60 == 1L] <- rtrunc_lnorm(n_dec, meanlog = log(9), sdlog = 0.95,
                                    lo = 0.04, hi = HORIZON)
  cross_sectional <- cross_sectional %>%
    mutate(
      event      = as.integer(died60),
      time       = if_else(died60 == 1L, tte, as.numeric(HORIZON)),
      # in-hospital mortality proxy for the synthetic secondary endpoint
      deceased   = as.integer(died60)
    )
} else {
  cross_sectional <- cross_sectional %>%
    mutate(
      idx_to_death = as.numeric(difftime(death_dttm, recorded_dttm, units = "days")),
      event = as.integer(!is.na(idx_to_death) & idx_to_death >= 0 &
                           idx_to_death <= HORIZON),
      time  = if_else(event == 1L, idx_to_death, as.numeric(HORIZON))
    )
}

# Guard against zero/negative follow-up (e.g. death recorded same instant as t0).
cross_sectional <- cross_sectional %>%
  mutate(time = pmax(time, 1 / 24))

message("Re-anchored 60-day mortality: ", sum(cross_sectional$event), " / ",
        nrow(cross_sectional), " (",
        round(100 * mean(cross_sectional$event), 1), "%)")

# =============================================================================
# 6c. Metric registry and confounder sets
# =============================================================================
# Three families x three normalizations. norm: "raw" | "pbw" | "pfvc".
metric_registry <- tribble(
  ~family,            ~norm,  ~var,               ~label,
  "Driving pressure", "raw",  "dp",               "DP",
  "Driving pressure", "pbw",  "dp_pbw",           "DP / PBW",
  "Driving pressure", "pfvc", "dp_pfvc",          "DP / PFVC",
  "Mechanical power", "raw",  "mechanical_power", "MP",
  "Mechanical power", "pbw",  "mp_pbw",           "MP / PBW",
  "Mechanical power", "pfvc", "mp_pfvc",          "MP / PFVC",
  "Elastance",        "raw",  "ers",              "Ers",
  "Elastance",        "pbw",  "ers_pbw",          "Ers x PBW",
  "Elastance",        "pfvc", "ers_pfvc",         "Ers x PFVC"
)

# Confounder strings (reuse script-04 conventions: severity + dose + habitus +/-
# demographics). vtpbw (dose) is included so threshold variation reflects
# mechanics, not delivered volume (the DP-adjust-vtpbw convention).
demo_terms   <- "age10 + sex_category + race_category"
sev_terms    <- "sofa_total + sf10 + bmi + vtpbw"
cov_adjusted   <- paste(sev_terms, "+", demo_terms)
cov_unadjusted <- sev_terms
# Covariate set for the mechanic x age interaction test (age enters via the
# interaction, so it is dropped from the main-effect block to avoid duplication).
cov_noage      <- "sofa_total + sf10 + bmi + vtpbw + sex_category + race_category"

# Global age-tertile cutpoints (shared across all metrics so the strata are
# comparable). Stratum-specific optimal thresholds are estimated within these.
AGE_TERTILE_BREAKS <- quantile(cross_sectional$age_at_admission, c(1/3, 2/3),
                               na.rm = TRUE)
AGE_TERTILE_LABELS <- c("T1_young", "T2_mid", "T3_old")

# Instruments for the IV arm (anthropometric mis-sizing). PRIMARY = height: it
# drives delivered VT (via PBW) and hence raw DP/MP/Ers, shares NO algebra with
# the exposures, and on real (UCSF) data is a strong instrument for the raw
# metrics (F: raw MP ~55, raw Ers ~370) where PBW/PFVC is near-zero. SECONDARY =
# pbwpfvc, kept as a sensitivity check (its strong first stage on PFVC-normalized
# metrics is largely a shared-1/PFVC-denominator artifact). INSTRUMENT aliases
# the primary for the absorption test and the complete-case requirement.
INSTRUMENT_PRIMARY   <- "height_cm"
INSTRUMENT_SECONDARY <- "pbwpfvc"
INSTRUMENT           <- INSTRUMENT_PRIMARY

# Complete-case frame for a given metric: needs the metric, all confounders and
# both instruments observed. (DP-derived metrics are missing where plateau
# pressure was not recorded — a legitimate, logged reduction, not a fallback.)
metric_frame <- function(var) {
  needed <- c("time", "event", "deceased", var,
              INSTRUMENT_PRIMARY, INSTRUMENT_SECONDARY,
              "age10", "sex_category", "race_category",
              "sofa_total", "sf10", "bmi", "vtpbw",
              "pbw", "pfvc", "age_at_admission", "height_cm")
  df <- cross_sectional %>%
    select(any_of(needed)) %>%
    filter(if_all(all_of(c("time", "event", var,
                           INSTRUMENT_PRIMARY, INSTRUMENT_SECONDARY,
                           "sofa_total", "sf10", "bmi", "vtpbw",
                           "age10", "sex_category", "race_category")),
                  ~ !is.na(.)))
  df
}

# =============================================================================
# 6d. IPW threshold emulation: stabilized weights + weighted Cox across a grid
# =============================================================================
# Treatment A = 1{metric <= tau} ("protective target achieved"); the reported
# contrast is the IPW-weighted hazard ratio for EXCEEDING tau (harm), with an
# unweighted covariate-adjusted Cox as corroboration.

stabilized_weights <- function(df, treat, ps_formula) {
  ps_fit <- glm(ps_formula, data = df, family = binomial)
  ps     <- as.numeric(predict(ps_fit, type = "response"))
  p_marg <- mean(treat)
  sw <- ifelse(treat == 1, p_marg / ps, (1 - p_marg) / (1 - ps))
  # Truncate at 1st/99th percentile (reported) to tame extreme weights.
  qs <- quantile(sw, c(0.01, 0.99), na.rm = TRUE)
  list(sw = pmin(pmax(sw, qs[1]), qs[2]),
       trunc_frac = mean(sw < qs[1] | sw > qs[2]),
       ps = ps)
}

fit_threshold <- function(df, var, tau, cov_str) {
  treat <- as.integer(df[[var]] <= tau)         # 1 = protective (<= tau)
  n_below <- sum(treat == 1); n_above <- sum(treat == 0)
  if (n_below < MIN_ARM || n_above < MIN_ARM) return(NULL)

  df2 <- df %>% mutate(treat = treat, above = 1L - treat)
  ps_formula <- as.formula(paste("treat ~", cov_str))
  wts <- tryCatch(stabilized_weights(df2, treat, ps_formula), error = function(e) NULL)
  if (is.null(wts)) return(NULL)
  df2$sw <- wts$sw

  ipw_cox <- tryCatch(
    coxph(Surv(time, event) ~ above, data = df2, weights = sw, robust = TRUE),
    error = function(e) NULL)
  adj_cox <- tryCatch(
    coxph(as.formula(paste("Surv(time, event) ~ above +", cov_str)), data = df2),
    error = function(e) NULL)
  if (is.null(ipw_cox) || is.null(adj_cox)) return(NULL)

  s_ipw <- summary(ipw_cox)$conf.int
  s_adj <- summary(adj_cox)$conf.int["above", ]
  tibble(
    tau = tau,
    # ECDF percentile of this threshold within the metric — a COMMON x-axis so
    # the (incomparable native-unit) normalizations can be overlaid fairly.
    tau_pctile = mean(df[[var]] <= tau, na.rm = TRUE),
    n_below = n_below, n_above = n_above,
    ipw_hr = s_ipw["above", "exp(coef)"],
    ipw_lo = s_ipw["above", "lower .95"],
    ipw_hi = s_ipw["above", "upper .95"],
    adj_hr = s_adj["exp(coef)"],
    adj_lo = s_adj["lower .95"],
    adj_hi = s_adj["upper .95"],
    weight_trunc_frac = wts$trunc_frac,
    ess = sum(df2$sw)^2 / sum(df2$sw^2)         # effective sample size
  )
}

run_threshold_grid <- function(df, var, cov_str, adjustment) {
  grid <- quantile(df[[var]], probs = seq(0.10, 0.90, length.out = N_GRID),
                   na.rm = TRUE)
  grid <- unique(round(grid, 6))
  out <- map_dfr(grid, ~ {
    r <- fit_threshold(df, var, .x, cov_str)
    if (!is.null(r)) r
  })
  if (nrow(out) == 0) return(out)
  out %>% mutate(adjustment = adjustment, .before = 1)
}

# Optimal threshold from the grid: lowest tau whose IPW HR for exceedance is
# significantly > 1 (harm emerges); fall back to the max-HR tau if none clears.
optimal_from_grid <- function(grid_df) {
  if (nrow(grid_df) == 0) return(tibble(opt_tau = NA_real_, opt_hr = NA_real_,
                                        opt_rule = "none"))
  harmful <- grid_df %>% filter(ipw_lo > 1) %>% arrange(tau)
  if (nrow(harmful) > 0) {
    return(tibble(opt_tau = harmful$tau[1], opt_hr = harmful$ipw_hr[1],
                  opt_rule = "first tau with IPW HR CI > 1"))
  }
  top <- grid_df %>% arrange(desc(ipw_hr)) %>% slice(1)
  tibble(opt_tau = top$tau, opt_hr = top$ipw_hr,
         opt_rule = "max IPW HR (no threshold reached significance)")
}

# =============================================================================
# 6e. maxstat data-driven cutpoint + bootstrap stability
# =============================================================================
maxstat_cutpoint <- function(df, var) {
  fml <- as.formula(paste("Surv(time, event) ~", var))
  ms  <- tryCatch(maxstat.test(fml, data = df, smethod = "LogRank",
                               pmethod = "Lau92"),
                  error = function(e) NULL)
  if (is.null(ms)) return(c(cut = NA_real_, p = NA_real_))
  c(cut = unname(ms$estimate), p = unname(ms$p.value))
}

maxstat_bootstrap <- function(df, var, B) {
  cuts <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(nrow(df), replace = TRUE)
    cuts[b] <- maxstat_cutpoint(df[idx, , drop = FALSE], var)["cut"]
  }
  cuts[is.finite(cuts)]
}

# =============================================================================
# 6f. Optimism-corrected discrimination (Harrell C via bootstrap)
# =============================================================================
cstat_apparent <- function(df, var) {
  fit <- coxph(as.formula(paste("Surv(time, event) ~", var)), data = df)
  as.numeric(summary(fit)$concordance["C"])
}

cstat_optimism_corrected <- function(df, var, B) {
  app <- cstat_apparent(df, var)
  opt <- numeric(0)
  for (b in seq_len(B)) {
    idx <- sample.int(nrow(df), replace = TRUE)
    boot_df <- df[idx, , drop = FALSE]
    fit_b <- tryCatch(
      coxph(as.formula(paste("Surv(time, event) ~", var)), data = boot_df),
      error = function(e) NULL)
    if (is.null(fit_b)) next
    c_boot <- as.numeric(summary(fit_b)$concordance["C"])
    # evaluate the bootstrap model on the original data
    lp_orig <- as.numeric(predict(fit_b, newdata = df, type = "lp"))
    c_orig  <- survival::concordance(Surv(df$time, df$event) ~ lp_orig,
                                     reverse = TRUE)$concordance
    opt <- c(opt, c_boot - c_orig)
  }
  list(c_apparent = app, optimism = mean(opt),
       c_corrected = app - mean(opt))
}

# =============================================================================
# 6g. IV / natural-experiment arm (2SRI + LPM + weak-IV-robust + sensitivity)
# =============================================================================
# PRIMARY instrument Z = height (anthropometric mis-sizing). SECONDARY = pbwpfvc.
# Height drives delivered VT (via PBW) and hence raw DP/MP/Ers, shares no algebra
# with the exposures, and on UCSF is strong for the raw metrics (the causal
# targets). The IV covariate set is severity + DOSE (vtpbw) + demographics — i.e.
# the same adjustment set as every other model here. vtpbw IS included: VT =
# (VT/PBW) x PBW and PBW = f(height), so the instrument enters the exposure
# through PBW, NOT through VT/PBW. VT/PBW is the clinician's dose setting (a
# separate input), so conditioning on it does not block height's pathway (height
# still moves VT via PBW with the dose fixed); it isolates the mechanics/size
# channel from the dose, per the project's DP-adjust-vtpbw convention. Age sits
# on BOTH sides — it shifts predicted volumes (PBW & PFVC) AND directly modifies
# elastic recoil (lung-dominated in hypoxemic failure) — so it is conditioned on
# in every IV/absorption model, keeping the instrument's residual net of age.
#
# Why pbwpfvc is only the SECONDARY/sensitivity instrument: the PFVC-normalized
# metrics (dp_pfvc, mp_pfvc) share the 1/pfvc denominator with pbw/pfvc, which
# mechanically inflates its first-stage F there, and it is near-zero for the raw
# metrics. Height has neither problem. Weak-IV-robust (Anderson-Rubin) intervals
# are reported for both so weak first stages surface as unbounded sets.
#
# All metrics are z-scored before the IV so effects are PER 1 SD of the metric,
# i.e. comparable across the (wildly different-scaled) exposures.
iv_cov <- cov_adjusted   # severity + dose (vtpbw) + demographics

# Just-identified 2SLS (one endogenous regressor d, one instrument z, exogenous
# X incl. intercept) with HC1 heteroskedasticity-robust SEs, base R (no IV pkg).
iv2sls_robust <- function(y, d, z, X) {
  W  <- cbind(d, X)                 # endogenous + exogenous (X has intercept)
  Zf <- cbind(z, X)                 # instrument + exogenous
  What <- Zf %*% (solve(crossprod(Zf)) %*% crossprod(Zf, W))  # stage-1 fitted
  bread <- solve(crossprod(What))
  beta  <- bread %*% crossprod(What, y)
  resid <- as.numeric(y - W %*% beta)
  n <- length(y); k <- ncol(W)
  meat <- crossprod(What * resid)                 # HC0 meat
  vcov <- (n / (n - k)) * (bread %*% meat %*% bread)  # HC1
  list(beta_d = beta[1, 1], se_d = sqrt(vcov[1, 1]))
}

iv_first_stage_F <- function(df, var, instr) {
  full <- lm(as.formula(paste(var, "~", instr, "+", iv_cov)), data = df)
  red  <- lm(as.formula(paste(var, "~", iv_cov)), data = df)
  an   <- anova(red, full)
  list(F = an$F[2], p = an$`Pr(>F)`[2], resid = resid(full))
}

iv_2sri <- function(df, var, B, instr) {
  fs <- iv_first_stage_F(df, var, instr)
  df2 <- df %>% mutate(.ctrl = fs$resid)
  stage2 <- coxph(as.formula(paste("Surv(time, event) ~", var, "+ .ctrl +", iv_cov)),
                  data = df2)
  beta <- coef(stage2)[[var]]
  bvec <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(nrow(df), replace = TRUE)
    bdf <- df[idx, , drop = FALSE]
    fb  <- tryCatch(lm(as.formula(paste(var, "~", instr, "+", iv_cov)), data = bdf),
                    error = function(e) NULL)
    if (is.null(fb)) { bvec[b] <- NA; next }
    bdf$.ctrl <- resid(fb)
    s2 <- tryCatch(
      coxph(as.formula(paste("Surv(time, event) ~", var, "+ .ctrl +", iv_cov)), data = bdf),
      error = function(e) NULL)
    bvec[b] <- if (is.null(s2)) NA else coef(s2)[[var]]
  }
  bvec <- bvec[is.finite(bvec)]
  tibble(first_stage_F = fs$F, weak_instrument = fs$F < 10,
         sri_hr_sd = exp(beta),
         sri_lo = exp(quantile(bvec, 0.025)),
         sri_hi = exp(quantile(bvec, 0.975)))
}

# Anderson-Rubin weak-IV-robust 95% confidence set for the LPM coefficient
# (risk difference per 1 SD). Inverts the AR test on a grid; an accepted region
# touching the grid edge is reported unbounded (the honest signature of a weak
# instrument). Uses FWL residualization on the exogenous controls.
iv_ar_ci <- function(df, var, instr, k_se = 40, n_grid = 4001) {
  X <- model.matrix(as.formula(paste("~", iv_cov)), data = df)
  y <- df$deceased; d <- df[[var]]; z <- df[[instr]]
  ry <- lm.fit(X, y)$residuals
  rd <- lm.fit(X, d)$residuals
  rz <- lm.fit(X, z)$residuals
  szz <- sum(rz^2)
  if (szz <= 0 || sum(rz * rd) == 0)
    return(tibble(ar_lo = NA_real_, ar_hi = NA_real_, ar_unbounded = NA))
  b_iv <- sum(rz * ry) / sum(rz * rd)
  se   <- iv2sls_robust(y, d, z, X)$se_d
  grid <- seq(b_iv - k_se * se, b_iv + k_se * se, length.out = n_grid)
  n <- length(y)
  accept <- vapply(grid, function(b0) {
    r <- ry - b0 * rd
    slope <- sum(rz * r) / szz
    e <- r - slope * rz
    var_slope <- (n / (n - 1)) * sum(rz^2 * e^2) / szz^2   # HC1 robust
    abs(slope / sqrt(var_slope)) < 1.96
  }, logical(1))
  if (!any(accept)) return(tibble(ar_lo = NA_real_, ar_hi = NA_real_, ar_unbounded = TRUE))
  acc <- which(accept)
  tibble(ar_lo = grid[min(acc)], ar_hi = grid[max(acc)],
         ar_unbounded = (min(acc) == 1 || max(acc) == n_grid))
}

# LPM 2SLS point estimate (risk difference per 1 SD) + AR CI + Conley
# "plausibly exogenous" union bound (assumed direct Z->Y effect up to the
# reduced-form coefficient = most adversarial).
iv_lpm <- function(df, var, instr) {
  X <- model.matrix(as.formula(paste("~", iv_cov)), data = df)
  y <- df$deceased; d <- df[[var]]; z <- df[[instr]]
  fit <- tryCatch(iv2sls_robust(y, d, z, X), error = function(e) NULL)
  if (is.null(fit)) return(tibble(lpm_beta_sd = NA, lpm_lo = NA, lpm_hi = NA,
                                  ar_lo = NA, ar_hi = NA, ar_unbounded = NA,
                                  conley_lo = NA, conley_hi = NA))
  est <- fit$beta_d; se <- fit$se_d
  ar <- iv_ar_ci(df, var, instr)
  rf <- lm(as.formula(paste("deceased ~", instr, "+", iv_cov)), data = df)
  gamma_max <- abs(coef(rf)[[instr]])
  lo <- Inf; hi <- -Inf
  for (delta in seq(0, gamma_max, length.out = 25)) {
    f2 <- tryCatch(iv2sls_robust(y - delta * z, d, z, X), error = function(e) NULL)
    if (is.null(f2)) next
    lo <- min(lo, f2$beta_d - 1.96 * f2$se_d)
    hi <- max(hi, f2$beta_d + 1.96 * f2$se_d)
  }
  bind_cols(tibble(lpm_beta_sd = est, lpm_lo = est - 1.96 * se,
                   lpm_hi = est + 1.96 * se), ar,
            tibble(conley_lo = lo, conley_hi = hi))
}

# Run the full IV battery for one metric under one instrument (metric z-scored).
iv_battery <- function(df, var, instr, B) {
  df_std <- df %>% mutate(.metric_std = as.numeric(scale(.data[[var]])))
  bind_cols(iv_2sri(df_std, ".metric_std", B, instr),
            iv_lpm(df_std, ".metric_std", instr))
}

# =============================================================================
# 6h. Normalization-absorption test (causal headline)
# =============================================================================
# Does the anthropometric instrument stay predictive of mortality once we
# condition on the (normalized) metric? The normalization that ABSORBS the
# instrument (drives its conditional HR to ~1) is the better proxy for true lung
# size. Called with the primary instrument (height), which shares no denominator
# with any metric, so the test is clean across all normalizations.
absorption_test <- function(df, var, instr) {
  fml <- as.formula(paste("Surv(time, event) ~", var, "+", instr, "+", cov_adjusted))
  fit <- coxph(fml, data = df)
  s   <- summary(fit)$conf.int[instr, ]
  pv  <- summary(fit)$coefficients[instr, "Pr(>|z|)"]
  tibble(z_hr = s["exp(coef)"], z_lo = s["lower .95"], z_hi = s["upper .95"],
         z_p = pv)
}

# =============================================================================
# 6h2. Age-stratified thresholds + mechanic x age interaction
# =============================================================================
# Does the optimal mechanical-stress threshold differ by age? The interaction
# test (Cox LRT of metric x age, adjusted) is REPORTED but does NOT gate the
# stratified thresholds (always produced). Within each age tertile we report the
# maxstat optimal cutpoint (prognostic) and the height-IV / AR-robust effect
# (severity-robust causal support). Per the project's age-bias framing, the RAW
# metrics are the cleaner test (PFVC already absorbs age via its reference eqs).

# Mechanic x age interaction on mortality (metric z-scored; manual LRT so it is
# robust to anova column-naming). Returns the LRT p and the per-decade shift in
# the (1-SD) metric HR.
mechanic_age_interaction <- function(df, var) {
  d <- df %>% mutate(.mc = as.numeric(scale(.data[[var]])))
  full <- coxph(as.formula(paste("Surv(time, event) ~ .mc * age10 +", cov_noage)), data = d)
  red  <- coxph(as.formula(paste("Surv(time, event) ~ .mc + age10 +", cov_noage)), data = d)
  lrt  <- 2 * as.numeric(logLik(full) - logLik(red))
  ddf  <- length(coef(full)) - length(coef(red))
  tibble(age_int_lrt_p = pchisq(lrt, ddf, lower.tail = FALSE),
         age_int_hr_per_decade = exp(unname(coef(full)[".mc:age10"])))
}

# Optimal threshold + severity-robust effect within one age stratum. Returns an
# all-NA row (preserving columns) when the stratum is too small to estimate.
stratum_threshold <- function(sdf, var) {
  na_row <- tibble(n = nrow(sdf), maxstat_cut = NA_real_, maxstat_boot_lo = NA_real_,
                   maxstat_boot_hi = NA_real_, iv_first_stage_F = NA_real_,
                   iv_weak = NA, sri_hr_sd = NA_real_, ar_lo = NA_real_,
                   ar_hi = NA_real_, ar_unbounded = NA)
  if (nrow(sdf) < 2 * MIN_ARM) return(na_row)
  ms  <- maxstat_cutpoint(sdf, var)
  msb <- maxstat_bootstrap(sdf, var, N_BOOT_MAXSTAT)
  ivb <- iv_battery(sdf, var, INSTRUMENT_PRIMARY, N_BOOT_STRATUM)
  tibble(n = nrow(sdf),
         maxstat_cut = unname(ms["cut"]),
         maxstat_boot_lo = quantile(msb, 0.025, na.rm = TRUE),
         maxstat_boot_hi = quantile(msb, 0.975, na.rm = TRUE),
         iv_first_stage_F = ivb$first_stage_F, iv_weak = ivb$weak_instrument,
         sri_hr_sd = ivb$sri_hr_sd, ar_lo = ivb$ar_lo, ar_hi = ivb$ar_hi,
         ar_unbounded = ivb$ar_unbounded)
}

# Demographic patterning of threshold-crossing: when a metric's (median-)threshold
# label "harmful" is regressed on demographics, how much of WHO is flagged tracks
# age/sex/race? (PFVC normalization injects demographic structure via PFVC itself.)
demographic_patterning <- function(df, var) {
  tau <- median(df[[var]], na.rm = TRUE)
  d2 <- df %>% mutate(harmful = as.integer(.data[[var]] > tau))
  if (length(unique(d2$harmful)) < 2) return(tibble(demo_lrt_p = NA, pseudo_r2 = NA))
  full <- glm(as.formula(paste("harmful ~", demo_terms)), data = d2, family = binomial)
  null <- glm(harmful ~ 1, data = d2, family = binomial)
  lrt  <- anova(null, full, test = "LRT")
  tibble(demo_lrt_p = lrt$`Pr(>Chi)`[2],
         pseudo_r2  = 1 - as.numeric(logLik(full)) / as.numeric(logLik(null)))
}

# =============================================================================
# 6i. Run everything over the metric registry
# =============================================================================
# All work for one metric, returned with an explicit status so failures are
# never silent: "ok" (results), "skipped" (too few complete cases), or "error"
# (the worker hit an R error — surfaced, not swallowed). Independent across
# metrics -> dispatched over a PSOCK cluster below.
analyze_metric <- function(i) {
  m <- metric_registry[i, ]
  var <- m$var; lab <- m$label
  tryCatch({
    df <- metric_frame(var)
    if (nrow(df) < 2 * MIN_ARM)
      return(list(status = "skipped", metric = lab, n = nrow(df)))
    meta <- tibble(family = m$family, norm = m$norm, metric = lab)

    # --- (A) IPW threshold grid, adjusted + unadjusted ----------------------
    g_adj   <- run_threshold_grid(df, var, cov_adjusted,   "adjusted")
    g_unadj <- run_threshold_grid(df, var, cov_unadjusted, "unadjusted")
    grid_m  <- bind_rows(g_adj, g_unadj)
    grid_out <- if (nrow(grid_m) > 0)
      grid_m %>% mutate(family = m$family, norm = m$norm, metric = lab, .before = 1) else NULL

    # --- optimal threshold + maxstat ----------------------------------------
    opt <- optimal_from_grid(g_adj)
    ms  <- maxstat_cutpoint(df, var)
    ms_boot <- maxstat_bootstrap(df, var, N_BOOT_MAXSTAT)
    optimal_out <- bind_cols(meta, tibble(
      grid_opt_tau = opt$opt_tau, grid_opt_hr = opt$opt_hr, grid_rule = opt$opt_rule,
      maxstat_cut = unname(ms["cut"]), maxstat_p = unname(ms["p"]),
      # percentile of the maxstat cutpoint within the metric, for the common-axis plot
      maxstat_cut_pctile = mean(df[[var]] <= ms["cut"], na.rm = TRUE),
      maxstat_boot_median = median(ms_boot, na.rm = TRUE),
      maxstat_boot_lo = quantile(ms_boot, 0.025, na.rm = TRUE),
      maxstat_boot_hi = quantile(ms_boot, 0.975, na.rm = TRUE)))

    # --- discrimination -----------------------------------------------------
    cs <- cstat_optimism_corrected(df, var, N_BOOT_CSTAT)
    discrim_out <- bind_cols(meta, tibble(
      c_apparent = cs$c_apparent, optimism = cs$optimism, c_corrected = cs$c_corrected))

    # --- (B) IV: primary (height) + secondary (pbwpfvc) instruments ---------
    iv_out  <- bind_cols(meta, instrument = INSTRUMENT_PRIMARY,
                         iv_battery(df, var, INSTRUMENT_PRIMARY,   N_BOOT))
    iv2_out <- bind_cols(meta, instrument = INSTRUMENT_SECONDARY,
                         iv_battery(df, var, INSTRUMENT_SECONDARY, N_BOOT))

    # --- (C) absorption test + demographic patterning -----------------------
    abs_out  <- bind_cols(meta, absorption_test(df, var, INSTRUMENT))
    demo_out <- bind_cols(meta, demographic_patterning(df, var))

    # --- (D) age-stratified thresholds + mechanic x age interaction ---------
    df_age <- df %>%
      mutate(age_tertile = cut(age_at_admission,
                               breaks = c(-Inf, AGE_TERTILE_BREAKS, Inf),
                               labels = AGE_TERTILE_LABELS))
    by_age_out <- map_dfr(AGE_TERTILE_LABELS, function(lv) {
      sdf <- df_age %>% filter(age_tertile == lv)
      bind_cols(tibble(age_tertile = lv), stratum_threshold(sdf, var))
    }) %>% mutate(family = m$family, norm = m$norm, metric = lab, .before = 1)
    age_int_out <- bind_cols(meta, mechanic_age_interaction(df, var))

    list(status = "ok", metric = lab,
         grid = grid_out, optimal = optimal_out, discrim = discrim_out,
         iv = iv_out, iv_secondary = iv2_out, absorption = abs_out, demo = demo_out,
         by_age = by_age_out, age_int = age_int_out)
  },
  error = function(e) list(status = "error", metric = lab,
                           message = conditionMessage(e)))
}

n_cores_used <- min(N_CORES, nrow(metric_registry))
message("Running ", nrow(metric_registry), " metrics across ", n_cores_used,
        " worker(s) [", if (n_cores_used > 1) "PSOCK" else "serial", "]")

# PSOCK (separate R processes) avoids the fork + multithreaded-BLAS instability
# that makes mclapply crash nondeterministically on macOS, and works on Windows.
# A worker that dies raises a loud error here instead of returning a silent NULL.
if (n_cores_used > 1) {
  cl <- makeCluster(n_cores_used, type = "PSOCK")
  clusterEvalQ(cl, {
    library(tidyverse); library(survival); library(splines); library(maxstat)
  })
  clusterExport(cl, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
  clusterSetRNGStream(cl, 20260615)   # reproducible parallel bootstrap streams
  # finally{} guarantees the workers are torn down on normal completion, error,
  # OR Ctrl-C — so an interrupted run doesn't leave orphaned idle R processes.
  results <- tryCatch(
    parLapply(cl, seq_len(nrow(metric_registry)), analyze_metric),
    finally = stopCluster(cl))
} else {
  results <- lapply(seq_len(nrow(metric_registry)), analyze_metric)
}

# Each result carries an explicit status: "error" (surfaced loudly), "skipped"
# (too few complete cases — logged), or "ok".
status_of <- function(x) if (is.list(x)) x$status else NA_character_
errs <- Filter(function(x) identical(status_of(x), "error"), results)
if (length(errs) > 0)
  stop("Metric analysis failed:\n",
       paste0("  ", vapply(errs, `[[`, "", "metric"), ": ",
              vapply(errs, `[[`, "", "message"), collapse = "\n"))

skipped <- Filter(function(x) identical(status_of(x), "skipped"), results)
if (length(skipped) > 0)
  message("Skipped (< ", 2 * MIN_ARM, " complete cases): ",
          paste(vapply(skipped, `[[`, "", "metric"), collapse = ", "))

results <- Filter(function(x) identical(status_of(x), "ok"), results)
if (length(results) == 0) {
  # No metric had enough complete cases: report per-field completeness so the
  # offending column (e.g. bmi / pressures) is obvious.
  req_cov <- c("time", "event", "pbwpfvc", "sofa_total", "sf10", "bmi", "vtpbw",
               "age10", "sex_category", "race_category")
  cov_n   <- vapply(req_cov, function(v) sum(!is.na(cross_sectional[[v]])), integer(1))
  met_n   <- vapply(metric_registry$var,
                    function(v) sum(!is.na(cross_sectional[[v]])), integer(1))
  stop("No metric produced output: every metric had < ", 2 * MIN_ARM,
       " complete cases (cohort n = ", nrow(cross_sectional), ").\n",
       "Non-missing per required covariate:\n",
       paste0("  ", names(cov_n), ": ", cov_n, collapse = "\n"), "\n",
       "Non-missing per metric:\n",
       paste0("  ", names(met_n), " (", metric_registry$label, "): ", met_n,
              collapse = "\n"),
       "\nThe near-zero field above is the bottleneck (commonly bmi/weight or ",
       "plateau/peak pressure). Fix that upstream rather than dropping it here.")
}

grid_tbl       <- bind_rows(map(results, "grid"))
optimal_tbl    <- bind_rows(map(results, "optimal"))
discrim_tbl    <- bind_rows(map(results, "discrim"))
iv_tbl         <- bind_rows(map(results, "iv"))
iv_secondary_tbl <- bind_rows(map(results, "iv_secondary"))
absorption_tbl <- bind_rows(map(results, "absorption"))
demo_tbl       <- bind_rows(map(results, "demo"))
by_age_tbl     <- bind_rows(map(results, "by_age"))
age_int_tbl    <- bind_rows(map(results, "age_int"))

# =============================================================================
# 6j. Positivity contrast: mechanical metric overlap vs VT/PFVC negative control
# =============================================================================
# For one representative threshold per metric (the grid optimum), summarise the
# propensity overlap by demographic stratum, and contrast against a VT/PFVC
# threshold (expected near-deterministic separation = positivity failure).
positivity_summary <- function(df, var, tau, cov_str, label) {
  treat <- as.integer(df[[var]] <= tau)
  if (length(unique(treat)) < 2) {
    return(tibble(exposure = label, tau = tau, ps_min = NA, ps_max = NA,
                  frac_extreme = NA, ess = NA))
  }
  ps_fit <- glm(as.formula(paste("treat ~", cov_str)),
                data = df %>% mutate(treat = treat), family = binomial)
  ps <- as.numeric(predict(ps_fit, type = "response"))
  tibble(exposure = label, tau = tau,
         ps_min = min(ps), ps_max = max(ps),
         frac_extreme = mean(ps < 0.05 | ps > 0.95),
         ess = {
           p_marg <- mean(treat)
           sw <- ifelse(treat == 1, p_marg / ps, (1 - p_marg) / (1 - ps))
           sum(sw)^2 / sum(sw^2)
         })
}

pos_rows <- list()
for (i in seq_len(nrow(metric_registry))) {
  m <- metric_registry[i, ]; var <- m$var
  df <- metric_frame(var)
  if (nrow(df) < 2 * MIN_ARM) next
  tau <- optimal_tbl$grid_opt_tau[optimal_tbl$metric == m$label]
  if (length(tau) == 0 || is.na(tau)) tau <- median(df[[var]], na.rm = TRUE)
  pos_rows[[var]] <- positivity_summary(df, var, tau, cov_adjusted, m$label)
}
# Negative control: a VT/PFVC threshold on the full cohort.
vtpfvc_df <- cross_sectional %>%
  filter(if_all(all_of(c("vtpfvc", "age10", "sex_category", "race_category",
                         "sofa_total", "sf10", "bmi", "vtpbw")), ~ !is.na(.)))
pos_rows[["vtpfvc_negctrl"]] <- positivity_summary(
  vtpfvc_df, "vtpfvc", median(vtpfvc_df$vtpfvc, na.rm = TRUE),
  cov_adjusted, "VT/PFVC (neg. control)")
positivity_tbl <- bind_rows(pos_rows)

# =============================================================================
# 6k. Normalization head-to-head (assembled by joins from the per-metric battery)
# =============================================================================
# Raw vs PBW vs PFVC within each family: discrimination, instrument absorption,
# positivity quality, and the demographic patterning of threshold-crossing (the
# demographic_patterning is computed in the parallel battery above).
normcomp_tbl <- discrim_tbl %>%
  select(family, norm, metric, c_corrected) %>%
  left_join(absorption_tbl %>% select(metric, z_hr, z_p),              by = "metric") %>%
  left_join(positivity_tbl %>% select(metric = exposure, ess, frac_extreme), by = "metric") %>%
  left_join(demo_tbl %>% select(metric, demo_lrt_p, pseudo_r2),        by = "metric")

# =============================================================================
# 6l. Write tables
# =============================================================================
sfx <- function(stub) file.path(final_dir, paste0(stub, "_", site_name, ".csv"))
write_csv(grid_tbl,       sfx("tte_threshold_grid"))
write_csv(optimal_tbl,    sfx("tte_optimal_thresholds"))
write_csv(discrim_tbl,    sfx("tte_discrimination"))
write_csv(iv_tbl,           sfx("tte_iv_estimates"))            # primary: height
write_csv(iv_secondary_tbl, sfx("tte_iv_estimates_secondary"))  # sensitivity: pbwpfvc
write_csv(absorption_tbl, sfx("tte_normalization_absorption"))
write_csv(positivity_tbl, sfx("tte_positivity"))
write_csv(normcomp_tbl,   sfx("tte_normalization_comparison"))
write_csv(by_age_tbl,     sfx("tte_thresholds_by_age"))         # age-tertile thresholds + IV
write_csv(age_int_tbl,    sfx("tte_age_interaction"))           # mechanic x age LRT (reported, not gating)
message("Wrote 10 result tables to ", final_dir)

# =============================================================================
# 6m. Figures (Okabe-Ito discrete, viridis continuous)
# =============================================================================
pdf_path <- function(stub) file.path(final_dir, paste0(stub, "_", site_name, ".pdf"))
fam_levels <- c("Driving pressure", "Mechanical power", "Elastance")
norm_cols  <- c(raw = okabe[1], pbw = okabe[2], pfvc = okabe[3])

# --- dose-response: IPW HR vs metric PERCENTILE, so the (incomparable native-
# unit) normalizations overlay fairly on one 10-90% axis. maxstat cutpoints are
# rugged on the same axis. Native-unit thresholds live in tte_optimal_thresholds.
if (nrow(grid_tbl) > 0) {
  dr <- grid_tbl %>% filter(adjustment == "adjusted") %>%
    mutate(family = factor(family, levels = fam_levels), tau_pct = 100 * tau_pctile)
  ms_rug <- optimal_tbl %>%
    mutate(family = factor(family, levels = fam_levels), tau_pct = 100 * maxstat_cut_pctile)
  p_dr <- ggplot(dr, aes(tau_pct, ipw_hr, colour = norm, fill = norm)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
    geom_ribbon(aes(ymin = ipw_lo, ymax = ipw_hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_rug(data = ms_rug, aes(x = tau_pct, colour = norm),
             sides = "b", linewidth = 0.9, inherit.aes = FALSE) +
    facet_wrap(~ family) +
    scale_colour_manual(values = norm_cols, name = "Normalization") +
    scale_fill_manual(values = norm_cols, name = "Normalization") +
    scale_y_log10() +
    labs(x = "Threshold percentile within cohort (%)",
         y = "IPW hazard ratio for exceeding threshold (60-day mortality)",
         title = "Threshold dose-response by metric family and normalization",
         subtitle = paste0(site_name,
                           if (is_synthetic) " (SYNTHETIC - simulated survival)" else "",
                           " - x = within-cohort percentile (common axis); ticks = maxstat optimal cutpoint")) +
    theme_minimal(base_size = 11)
  ggsave(pdf_path("tte_doseresponse"), p_dr, width = 11, height = 4.5)
}

# --- IV panel: dual-instrument strength + AR-robust effect + absorption --------
if (nrow(iv_tbl) > 0) {
  met_lvls <- rev(metric_registry$label)
  facit <- function(x) factor(x, levels = met_lvls)

  # (1) Instrument strength: primary (height) vs secondary (pbwpfvc). Height is
  # strong for the RAW metrics (the causal targets); pbwpfvc is weak there and
  # only strong on PFVC-normalized metrics via the shared 1/PFVC denominator.
  iv_strength <- bind_rows(iv_tbl, iv_secondary_tbl) %>%
    transmute(metric = facit(metric), instrument, first_stage_F)
  p_f <- ggplot(iv_strength, aes(first_stage_F, metric, fill = instrument)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_vline(xintercept = 10, linetype = 2, colour = "grey40") +
    scale_fill_manual(values = c(height_cm = okabe[4], pbwpfvc = okabe[6]),
                      name = "Instrument", breaks = c("height_cm", "pbwpfvc")) +
    labs(y = NULL, x = "First-stage F (F<10 = weak)",
         title = "Instrument strength: height (primary) vs PBW/PFVC (secondary)",
         subtitle = "Height is strong for the raw metrics; PBW/PFVC's strength on PFVC-normalized is a shared-denominator artifact") +
    theme_minimal(base_size = 10)

  # (2) Weak-IV-robust effect (Anderson-Rubin), primary instrument. Risk
  # difference per 1 SD; unbounded AR sets are clipped to the panel (= weak).
  ar_plot <- iv_tbl %>% mutate(metric = facit(metric))
  fin <- c(ar_plot$ar_lo[!ar_plot$ar_unbounded %in% TRUE],
           ar_plot$ar_hi[!ar_plot$ar_unbounded %in% TRUE], ar_plot$lpm_beta_sd)
  fin <- fin[is.finite(fin)]
  xr <- if (length(fin)) range(fin) else c(-0.1, 0.1)
  xr <- xr + c(-1, 1) * diff(xr) * 0.1
  p_ar <- ggplot(ar_plot, aes(y = metric, colour = norm)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_segment(aes(x = pmax(ar_lo, xr[1]), xend = pmin(ar_hi, xr[2]),
                     y = metric, yend = metric), na.rm = TRUE) +
    geom_point(aes(x = lpm_beta_sd), na.rm = TRUE) +
    scale_colour_manual(values = norm_cols, guide = "none") +
    coord_cartesian(xlim = xr) +
    labs(x = "LPM risk difference per 1 SD (Anderson-Rubin 95% CI)", y = NULL,
         title = "IV effect, weak-instrument-robust (instrument = height)",
         subtitle = "CIs reaching a panel edge are unbounded = not identified") +
    theme_minimal(base_size = 10)

  # (3) Normalization-absorption test (primary instrument).
  p_abs <- ggplot(absorption_tbl %>% mutate(metric = facit(metric)),
                  aes(z_hr, metric, colour = norm)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = z_lo, xmax = z_hi)) +
    scale_colour_manual(values = norm_cols, guide = "none") +
    labs(x = "Instrument (height) HR | metric (~1 = absorbed)", y = NULL,
         title = "Normalization-absorption test",
         subtitle = "Height shares no denominator with any metric, so this is clean across normalizations") +
    theme_minimal(base_size = 10)

  ggsave(pdf_path("tte_iv"), p_f / p_ar / p_abs, width = 8.5, height = 11)
}

# --- positivity contrast ------------------------------------------------------
if (nrow(positivity_tbl) > 0) {
  pos_plot <- positivity_tbl %>%
    mutate(is_negctrl = grepl("neg", exposure),
           exposure = fct_reorder(exposure, frac_extreme))
  p_pos <- ggplot(pos_plot, aes(frac_extreme, exposure, fill = is_negctrl)) +
    geom_col() +
    scale_fill_manual(values = c(`FALSE` = okabe[3], `TRUE` = okabe[5]),
                      labels = c("Mechanical metric", "VT/PFVC neg. control"),
                      name = NULL) +
    labs(x = "Fraction with extreme propensity (<0.05 or >0.95)", y = NULL,
         title = "Positivity: mechanical metrics vs VT/PFVC negative control",
         subtitle = "Higher = worse overlap (positivity failure)") +
    theme_minimal(base_size = 11)
  ggsave(pdf_path("tte_positivity"), p_pos, width = 8, height = 5)
}

# --- age-stratified optimal thresholds (maxstat) ------------------------------
if (nrow(by_age_tbl) > 0 && any(!is.na(by_age_tbl$maxstat_cut))) {
  ba <- by_age_tbl %>%
    filter(!is.na(maxstat_cut)) %>%
    mutate(metric = factor(metric, levels = metric_registry$label),
           age_tertile = factor(age_tertile, levels = AGE_TERTILE_LABELS))
  # annotate each metric facet with the mechanic x age interaction p (reported,
  # not gating).
  p_lab <- age_int_tbl %>%
    mutate(metric = factor(metric, levels = metric_registry$label),
           lab = paste0("age-int p=", signif(age_int_lrt_p, 2)))
  p_ba <- ggplot(ba, aes(age_tertile, maxstat_cut, group = 1)) +
    geom_line(colour = okabe[2]) +
    geom_point(colour = okabe[2]) +
    geom_errorbar(aes(ymin = maxstat_boot_lo, ymax = maxstat_boot_hi),
                  width = 0.15, colour = okabe[2]) +
    geom_text(data = p_lab, aes(x = 2, y = Inf, label = lab), inherit.aes = FALSE,
              vjust = 1.4, size = 3, colour = "grey30") +
    facet_wrap(~ metric, scales = "free_y") +
    labs(x = "Age tertile", y = "maxstat optimal cutpoint (metric units)",
         title = "Age-stratified optimal thresholds (maxstat) per metric",
         subtitle = paste0(site_name,
           " - bootstrap 95% CI; prognostic (see tte_thresholds_by_age IV columns for causal support)")) +
    theme_minimal(base_size = 10)
  ggsave(pdf_path("tte_thresholds_by_age"), p_ba, width = 11, height = 8)
}

message("Wrote figures to ", final_dir)
message("Script 06 complete.")
