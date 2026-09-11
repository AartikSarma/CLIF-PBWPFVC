# =============================================================================
# Script 09: Height-instrumented PFVC-anchoring policy (modified treatment policy)
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# FEATURE BRANCH (feature/target-trial-emulation): NOT part of 00_run_pipeline.R.
# Standalone; reads the script-03 cross-sectional dataset.
#
# WHY THIS SCRIPT EXISTS ------------------------------------------------------
# A demographics-adjusted VT/PFVC contrast is UNIDENTIFIABLE: VT/PFVC is
# near-deterministic in demographics (script 07 positivity: c~0.996, ESS~0),
# and a static "everyone to <11%" target is positivity-violating in the
# bias-prone cells (extrapolation audit). The escape is two moves the height
# instrument and a modified treatment policy (MTP) make possible:
#
#   (A) INSTRUMENT, not adjust. Height supplies exogenous variation in mechanical
#       stress; identification rests on its exclusion restriction, NOT on
#       demographic overlap -- so the positivity wall does not apply. Height
#       drives PBW -> delivered VT -> raw DP/MP (and is strong there: UCSF
#       first-stage F ~ 74-657), sharing no algebra with the per-kg dose setting.
#
#   (B) SHIFT, not set. Instead of the static target, estimate a FEASIBLE dose
#       SHIFT toward PFVC-anchoring, clamped to each demographic cell's observed
#       VT/PFVC support [p5,p95]. The MTP positivity condition is only that the
#       shifted value lie in observed support -- far weaker than "both arms in
#       every cell", and exactly what the extrapolation audit checks.
#
# ESTIMAND. Under a (locally) linear structural model Y = beta*D + controls, the
# effect of an additive shift D -> D - delta is beta*E[delta]. We report the
# excess 60-day mortality of PBW dosing vs the PFVC-anchoring policy as
#   excess = beta_IV * mean(D_actual - D_policy),   (>0 = PBW dosing kills more)
# overall and by subgroup, with beta_IV from a height-instrumented LPM (2SLS,
# HC1). GENERALIZATION from the height-driven slice of mis-sizing to the full
# policy assumes the stress->mortality effect is invariant to the SOURCE of
# mis-sizing -- i.e. VT/PFVC is the structural dose -- which script 08
# establishes (homogeneous dose-response across strata and sites).
#
# LIMITATIONS. beta_IV is a local (complier) effect; exclusion (height affects
# mortality only via mechanics | covariates) is bounded, not tested, here (see
# script 06 Conley/E-value). DP/MP missing without plateau pressure (logged
# reduction). Expanded VT/PBW range needs the height IV for confounding by
# indication AND a wider upstream export than this protocol-restricted cohort.
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(splines)

source("utils/config.R")
site_name <- config$site_name

output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

okabe <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7",
           "#F0E442", "#999999")
RNGkind("L'Ecuyer-CMRG")
set.seed(20260616)

HORIZON       <- 60
MIN_CELL      <- 10
VTPFVC_TARGET <- 11        # static PFVC-anchoring target (VT/PFVC <11%, ARMA ~p75)
is_synthetic  <- identical(site_name, "synthetic_clif")
N_BOOT        <- if (is_synthetic) 200L else
  suppressWarnings(as.integer(Sys.getenv("PBWPFVC_NBOOT", unset = "500")))
if (is.na(N_BOOT)) N_BOOT <- 500L
message("Site: ", site_name, " | synthetic = ", is_synthetic, " | boot = ", N_BOOT)

# =============================================================================
# 9a. Load + re-anchor 60-day mortality (identical to scripts 06/07/08)
# =============================================================================
cross_sectional <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"))
rtrunc_lnorm <- function(n_needed, meanlog, sdlog, lo, hi) {
  acc <- numeric(0)
  while (length(acc) < n_needed) {
    cand <- rlnorm(max(n_needed * 2L, 1000L), meanlog, sdlog)
    cand <- cand[cand > lo & cand <= hi]; acc <- c(acc, cand)
  }
  acc[seq_len(n_needed)]
}
if (is_synthetic) {
  message("*** SYNTHETIC SITE: simulated long-tailed survival (synthetic CLIF mortality is unreliable). ***")
  set.seed(20260615)
  n <- nrow(cross_sectional); died60 <- rbinom(n, 1L, 0.35)
  tte <- rep(NA_real_, n); n_dec <- sum(died60 == 1L)
  tte[died60 == 1L] <- rtrunc_lnorm(n_dec, log(9), 0.95, 0.04, HORIZON)
  cross_sectional <- cross_sectional %>%
    mutate(event = as.integer(died60),
           time = if_else(died60 == 1L, tte, as.numeric(HORIZON)), deceased = as.integer(died60))
} else {
  cross_sectional <- cross_sectional %>%
    mutate(idx_to_death = as.numeric(difftime(death_dttm, recorded_dttm, units = "days")),
           event = as.integer(!is.na(idx_to_death) & idx_to_death >= 0 & idx_to_death <= HORIZON),
           time  = if_else(event == 1L, idx_to_death, as.numeric(HORIZON)), deceased = event)
}
cross_sectional <- cross_sectional %>% mutate(time = pmax(time, 1 / 24))
message("60-day mortality: ", sum(cross_sectional$event), " / ", nrow(cross_sectional),
        " (", round(100 * mean(cross_sectional$event), 1), "%)")

# =============================================================================
# 9b. Analytic frame, subgroup factors, and per-cell VT/PFVC support
# =============================================================================
age_breaks <- quantile(cross_sectional$age_at_admission, c(1/3, 2/3), na.rm = TRUE)
analytic <- cross_sectional %>%
  filter(!is.na(vtpbw), !is.na(vtpfvc), !is.na(pbwpfvc), !is.na(pbw), !is.na(pfvc),
         !is.na(height_cm), !is.na(sex_category), !is.na(race_category),
         !is.na(age_at_admission), !is.na(sofa_total), !is.na(sf_ratio), !is.na(bmi),
         !is.na(time), !is.na(event), vtpfvc > 0) %>%
  group_by(sex_category) %>% mutate(height_z = as.numeric(scale(height_cm))) %>% ungroup() %>%
  mutate(age_grp = cut(age_at_admission, c(-Inf, age_breaks, Inf),
                       labels = c("Young", "Middle", "Old")),
         height_grp = cut(height_z, c(-Inf, quantile(height_z, c(1/3, 2/3), na.rm = TRUE), Inf),
                          labels = c("Short", "Middle", "Tall")),
         sex_grp = factor(sex_category), race_grp = factor(race_category))

# Per demographic-cell observed VT/PFVC support, for the feasible (MTP) target.
cells <- analytic %>% group_by(age_grp, sex_grp, race_grp, height_grp) %>%
  summarise(vt_p05 = quantile(vtpfvc, 0.05), vt_p95 = quantile(vtpfvc, 0.95), .groups = "drop")
analytic <- analytic %>% left_join(cells, by = c("age_grp", "sex_grp", "race_grp", "height_grp"))

# MP/Crs (from script 03: MP * ers / 1000) = power referenced to MEASURED (not
# predicted) lung size. Positivity-safe (Crs has non-demographic variation, unlike
# PFVC) and ~ DP-like (MP/Crs ~ DP^2), so it may rescue the paradoxical raw-MP
# instrument. Caveat: more compliance-laden than DP, so the mediator/collider
# concern is correspondingly stronger.

vtpbw_rng <- range(analytic$vtpbw)
message("Analytic n = ", nrow(analytic), " | VT/PBW range ", round(vtpbw_rng[1], 2), "-",
        round(vtpbw_rng[2], 2))

subgroup_vars <- tribble(
  ~svar,        ~slabel,                     ~ref_level,
  "sex_grp",    "Sex",                       "Male",
  "race_grp",   "Race",                      "WHITE",
  "age_grp",    "Age tertile",               "Young",
  "height_grp", "Height tertile (w/in sex)", "Tall"
)
# Candidate exposures. Shiftable policy metrics (DP, MP, MP/Crs) move with VT
# under re-dosing; Ers is diagnostic-only (elastance is a lung property, unchanged
# by re-dosing VT). MP/Crs and MP share the VT^2 mechanics map.
candidates <- tribble(
  ~var,               ~label,    ~vt_exp, ~shiftable,
  "dp",               "DP",      1,       TRUE,
  "mechanical_power", "MP",      2,       TRUE,
  "mp_crs",           "MP/Crs",  2,       TRUE,
  "ers",              "Ers",     NA_real_, FALSE
)

# Instrument + control set: severity + dose + spline age + sex + race (vtpbw kept
# so height enters via PBW, not the dose lever).
INSTRUMENT <- "height_cm"
iv_cov <- "sofa_total + sf_ratio + bmi + vtpbw + ns(age_at_admission, df = 3) + sex_category + race_category"

# =============================================================================
# 9c. Just-identified 2SLS (HC1 SEs) + Anderson-Rubin weak-IV-robust CIs (base R)
# =============================================================================
iv2sls_robust <- function(y, d, z, X) {
  W  <- cbind(d, X); Zf <- cbind(z, X)
  What <- Zf %*% (solve(crossprod(Zf)) %*% crossprod(Zf, W))
  bread <- solve(crossprod(What)); beta <- bread %*% crossprod(What, y)
  resid <- as.numeric(y - W %*% beta); n <- length(y); k <- ncol(W)
  meat <- crossprod(What * resid); vcov <- (n / (n - k)) * (bread %*% meat %*% bread)
  list(beta_d = beta[1, 1], se_d = sqrt(vcov[1, 1]))
}
first_stage_F <- function(df, var) {
  full <- lm(as.formula(paste(var, "~", INSTRUMENT, "+", iv_cov)), data = df)
  red  <- lm(as.formula(paste(var, "~", iv_cov)), data = df)
  anova(red, full)$F[2]
}

# Anderson-Rubin confidence set for the structural coefficient beta (LPM scale,
# same scale as iv2sls_robust's beta_d). WEAK-INSTRUMENT-ROBUST: correct coverage
# regardless of first-stage strength, where the Wald CI is unreliable when F is
# small. Just-identified (1 instrument), so AR is the efficient weak-robust test.
# After partialling the covariates X out of y, d (endogenous) and z (instrument),
# the test AR(b0) <= F_crit reduces to a quadratic A*b^2 + B*b + C <= 0 in b0; the
# set is a bounded interval (A>0, disc>0) OR unbounded / whole-line / disconnected
# (A<=0) -- the latter honestly signals the data cannot pin beta down.
anderson_rubin <- function(y, d, z, X, alpha = 0.05) {
  rX <- function(v) as.numeric(v - X %*% solve(crossprod(X), crossprod(X, v)))  # partial out X
  yt <- rX(y); dt <- rX(d); zt <- rX(z)
  Szz <- sum(zt * zt); Szy <- sum(zt * yt); Szd <- sum(zt * dt)
  Syy <- sum(yt * yt); Syd <- sum(yt * dt); Sdd <- sum(dt * dt)
  m  <- length(y) - ncol(X) - 1                       # df2 of the AR F(1, m)
  cc <- qf(1 - alpha, 1, m)
  A  <- (m + cc) * Szd^2 - cc * Szz * Sdd
  B  <- -2 * (m + cc) * Szy * Szd + 2 * cc * Szz * Syd
  C  <- (m + cc) * Szy^2 - cc * Szz * Syy
  disc <- B^2 - 4 * A * C
  if (is.finite(A) && A > 0 && disc > 0) {
    r <- sort(c((-B - sqrt(disc)) / (2 * A), (-B + sqrt(disc)) / (2 * A)))
    list(ar_lo = r[1], ar_hi = r[2], ar_unbounded = FALSE)
  } else {
    list(ar_lo = NA_real_, ar_hi = NA_real_, ar_unbounded = TRUE)
  }
}

# =============================================================================
# 9d. Policy estimation: excess mortality of PBW vs PFVC-anchoring (per metric,
#     per policy), height-instrumented, by subgroup + net population
# =============================================================================
# Two policies: "static" (target = 11% for everyone) and "feasible" (target
# clamped to the cell's observed [p5,p95] support -> a valid MTP).
shift_frame <- function(df, var, vt_exp, policy) {
  df %>% mutate(
    tau_i = if (policy == "feasible") pmin(pmax(VTPFVC_TARGET, vt_p05), vt_p95) else VTPFVC_TARGET,
    r     = tau_i / vtpfvc,                     # VT_policy / VT_actual
    d_cf  = .data[[var]] * r^vt_exp,            # mechanics map (DP ~ VT, MP ~ VT^2)
    shift = .data[[var]] - d_cf)                # >0 = de-escalation under policy
}
by_subgroup_shift <- function(d) {
  sub <- map_dfr(seq_len(nrow(subgroup_vars)), function(i) {
    d %>% group_by(level = .data[[subgroup_vars$svar[i]]]) %>%
      summarise(n = n(), mean_shift = mean(shift), .groups = "drop") %>%
      filter(n >= MIN_CELL) %>%
      mutate(subgroup = subgroup_vars$slabel[i], level = as.character(level))
  })
  ov <- d %>% summarise(n = n(), mean_shift = mean(shift)) %>%
    mutate(subgroup = "Overall", level = "All")
  bind_rows(sub, ov)
}

estimate_policy <- function(var, lbl, vt_exp, policy) {
  d <- analytic %>% filter(!is.na(.data[[var]])) %>% shift_frame(var, vt_exp, policy)
  X <- model.matrix(as.formula(paste("~", iv_cov)), data = d)
  iv  <- iv2sls_robust(d$deceased, d[[var]], d[[INSTRUMENT]], X)
  Fst <- first_stage_F(d, var)
  ols <- coef(lm(as.formula(paste("deceased ~", var, "+", iv_cov)), data = d))[[var]]
  point <- by_subgroup_shift(d) %>% mutate(excess_pp = 100 * iv$beta_d * mean_shift)

  acc <- vector("list", N_BOOT)
  for (b in seq_len(N_BOOT)) {
    idx <- sample.int(nrow(d), replace = TRUE); bd <- d[idx, , drop = FALSE]
    Xb  <- model.matrix(as.formula(paste("~", iv_cov)), data = bd)
    ivb <- tryCatch(iv2sls_robust(bd$deceased, bd[[var]], bd[[INSTRUMENT]], Xb),
                    error = function(e) NULL)
    if (is.null(ivb)) next
    acc[[b]] <- by_subgroup_shift(bd) %>%
      transmute(subgroup, level, excess_b = 100 * ivb$beta_d * mean_shift)
  }
  ci <- bind_rows(acc) %>% group_by(subgroup, level) %>%
    summarise(excess_lo = quantile(excess_b, 0.025, na.rm = TRUE),
              excess_hi = quantile(excess_b, 0.975, na.rm = TRUE), .groups = "drop")
  # MTP support: fraction whose STATIC target lies in their cell's observed range
  frac_support <- mean(d$vt_p05 <= VTPFVC_TARGET & VTPFVC_TARGET <= d$vt_p95)
  point %>% left_join(ci, by = c("subgroup", "level")) %>%
    mutate(metric = lbl, policy = policy, first_stage_F = Fst, weak_instrument = Fst < 10,
           beta_iv_per_unit = iv$beta_d, beta_ols_per_unit = ols,
           frac_static_in_support = frac_support, .before = 1)
}

# --- Exposure-IV validity diagnostic (gate the policy on this) -----------------
# For each candidate exposure: first-stage F, first-stage SIGN (height -> exposure),
# the IV vs OLS coefficient, and a validity flag. "valid" = strong (F>=10) AND
# IV/OLS sign-concordant. A strong-F exposure whose IV flips sign vs OLS (e.g. raw
# MP: height raises VT but also raises compliance, opposing channels) is an
# INVALID instrument despite a large F, and must not enter the policy.
iv_diagnostic <- function(var, lbl) {
  d  <- analytic %>% filter(!is.na(.data[[var]]))
  X  <- model.matrix(as.formula(paste("~", iv_cov)), data = d)
  iv <- iv2sls_robust(d$deceased, d[[var]], d[[INSTRUMENT]], X)
  ar <- tryCatch(anderson_rubin(d$deceased, d[[var]], d[[INSTRUMENT]], X),
                 error = function(e) list(ar_lo = NA_real_, ar_hi = NA_real_, ar_unbounded = NA))
  fs <- lm(as.formula(paste(var, "~", INSTRUMENT, "+", iv_cov)), data = d)
  Fst <- first_stage_F(d, var)
  ols <- coef(lm(as.formula(paste("deceased ~", var, "+", iv_cov)), data = d))[[var]]
  tibble(metric = lbl, n = nrow(d), first_stage_F = Fst,
         first_stage_sign = unname(sign(coef(fs)[[INSTRUMENT]])),
         beta_iv = iv$beta_d, beta_ols = unname(ols),
         wald_lo = iv$beta_d - 1.959964 * iv$se_d,      # Wald CI (unreliable if F small)
         wald_hi = iv$beta_d + 1.959964 * iv$se_d,
         ar_lo = ar$ar_lo, ar_hi = ar$ar_hi,            # Anderson-Rubin (weak-IV-robust)
         ar_unbounded = ar$ar_unbounded,
         sign_concordant = sign(iv$beta_d) == sign(ols),
         valid = Fst >= 10 & (sign(iv$beta_d) == sign(ols)))
}
diag_tbl <- pmap_dfr(candidates %>% select(var, label),
                     function(var, label) iv_diagnostic(var, label))
write_csv(diag_tbl, file.path(final_dir, paste0("ivpolicy_exposure_diagnostic_", site_name, ".csv")))

# AR vs Wald: they should coincide for strong-F (valid) metrics, and AR stays valid
# (often wide / unbounded) where the first stage is weak and Wald is unreliable.
message("\nIV diagnostic -- Anderson-Rubin (weak-IV-robust) vs Wald CIs:")
diag_tbl %>% pwalk(function(metric, first_stage_F, beta_iv, wald_lo, wald_hi,
                            ar_lo, ar_hi, ar_unbounded, valid, ...)
  message(sprintf("  %-7s F=%6.1f  beta_IV=%+.3f  Wald[%+.3f, %+.3f]  AR%s  [%s]",
                  metric, first_stage_F, beta_iv, wald_lo, wald_hi,
                  if (isTRUE(ar_unbounded)) "=unbounded" else
                    sprintf("[%+.3f, %+.3f]", ar_lo, ar_hi),
                  if (valid) "valid" else "weak/invalid")))
message("Exposure-IV diagnostic (valid = F>=10 AND IV/OLS sign-concordant):")
walk(seq_len(nrow(diag_tbl)), ~ message(sprintf(
  "  %-7s F=%6.1f  fs_sign=%+d  beta_IV=%+.4f  beta_OLS=%+.4f  valid=%s",
  diag_tbl$metric[.x], diag_tbl$first_stage_F[.x], diag_tbl$first_stage_sign[.x],
  diag_tbl$beta_iv[.x], diag_tbl$beta_ols[.x], diag_tbl$valid[.x])))

# --- Gate: estimate the policy ONLY for shiftable exposures that pass validity --
policy_metrics <- candidates %>% filter(shiftable) %>%
  left_join(diag_tbl %>% select(label = metric, valid), by = "label") %>%
  filter(valid)
if (nrow(policy_metrics) == 0) {
  message("No shiftable exposure passed the validity gate (F>=10 & sign-concordant); ",
          "policy NOT estimated. See ivpolicy_exposure_diagnostic_", site_name, ".csv.")
  policy_tbl <- tibble()
} else {
  message("Policy estimated for valid exposure(s): ",
          paste(policy_metrics$label, collapse = ", "))
  policy_tbl <- pmap_dfr(policy_metrics %>% select(var, label, vt_exp),
    function(var, label, vt_exp)
      map_dfr(c("static", "feasible"), function(pl) estimate_policy(var, label, vt_exp, pl)))
  write_csv(policy_tbl, file.path(final_dir, paste0("ivpolicy_estimates_", site_name, ".csv")))
  ov <- policy_tbl %>% filter(subgroup == "Overall")
  message("Net PFVC-anchoring effect (excess pp of PBW dosing, >0 = PBW worse):")
  walk(seq_len(nrow(ov)), ~ message(sprintf(
    "  %s %-8s F=%.0f  excess=%.2f [%.2f, %.2f]", ov$metric[.x], ov$policy[.x],
    ov$first_stage_F[.x], ov$excess_pp[.x], ov$excess_lo[.x], ov$excess_hi[.x])))
}

# =============================================================================
# 9e. Expanded VT/PBW range IV sensitivity (only if the cohort actually has range)
# =============================================================================
# Off-protocol dosing is confounded by indication, which the height IV -- not OLS
# -- is built to handle. We run it ONLY if the export carries VT/PBW beyond the
# protocol band; otherwise we log that it needs a wider upstream export rather
# than fabricate range.
frac_offband <- mean(analytic$vtpbw < 6 | analytic$vtpbw > 8)
if (frac_offband > 0.10 && nrow(policy_metrics) > 0) {
  message("Expanded-range IV sensitivity: ", round(100 * frac_offband, 1),
          "% off-band; refitting beta on full VT/PBW range.")
  expand_tbl <- pmap_dfr(policy_metrics %>% select(var, label, vt_exp),
                         function(var, label, vt_exp) {
    d <- analytic %>% filter(!is.na(.data[[var]])) %>% shift_frame(var, vt_exp, "feasible")
    X <- model.matrix(as.formula(paste("~", iv_cov)), data = d)
    iv <- iv2sls_robust(d$deceased, d[[var]], d[[INSTRUMENT]], X)
    by_subgroup_shift(d) %>% filter(subgroup == "Overall") %>%
      mutate(metric = label, first_stage_F = first_stage_F(d, var),
             beta_iv_per_unit = iv$beta_d, excess_pp = 100 * iv$beta_d * mean_shift)
  })
  write_csv(expand_tbl, file.path(final_dir, paste0("ivpolicy_expanded_range_", site_name, ".csv")))
} else {
  message("Expanded-range IV sensitivity SKIPPED: only ", round(100 * frac_offband, 1),
          "% of patients off the 6-8 band (cohort is protocol-restricted). ",
          "Requires a wider upstream export from script 03.")
}

# =============================================================================
# 9f. Figure: excess mortality by subgroup, static vs feasible policy, per metric
# =============================================================================
pdf_path <- function(stub) file.path(final_dir, paste0(stub, "_", site_name, ".pdf"))
sub_levels <- c("Male","Female","WHITE","BLACK","OTHER","Young","Middle","Old","Tall","Short","All")
if (nrow(policy_tbl) > 0) {
  pt <- policy_tbl %>%
    mutate(metric = factor(metric, levels = intersect(candidates$label, unique(policy_tbl$metric))),
           subgroup = factor(subgroup, levels = c(subgroup_vars$slabel, "Overall")),
           level = factor(level, levels = sub_levels),
           policy = factor(policy, levels = c("static", "feasible")))
  weaklab <- pt %>% distinct(metric, policy, first_stage_F, weak_instrument) %>%
    filter(weak_instrument) %>% mutate(lab = paste0("weak IV (F=", round(first_stage_F), ")"))
  p <- ggplot(pt, aes(level, excess_pp, fill = policy)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = excess_lo, ymax = excess_hi),
                  position = position_dodge(width = 0.8), width = 0.25) +
    facet_grid(metric ~ subgroup, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c(static = okabe[1], feasible = okabe[3]),
                      name = "PFVC-anchoring policy",
                      labels = c(static = "static (target 11%, may extrapolate)",
                                 feasible = "feasible (MTP, clamped to support)")) +
    labs(x = NULL, y = "Excess 60-day mortality, PBW vs PFVC-anchoring (pct points)",
         title = "Height-instrumented PFVC-anchoring policy effect (>0 = PBW dosing harms)",
         subtitle = paste0(site_name, if (is_synthetic) " (SYNTHETIC)" else "",
           " - IV identifies off exogenous height variation (no demographic positivity needed); ",
           "feasible policy stays in observed support")) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")
  ggsave(pdf_path("ivpolicy_excess_mortality"), p, width = 12, height = 6.5)
}

# =============================================================================
# 9g. Smooth age-varying policy effect + age-gradient flattening
# =============================================================================
# Manipulation-respecting integration of age. We do NOT estimate "the effect of
# age" (not modifiable). We estimate the effect of the MODIFIABLE PFVC-anchoring
# policy as a smooth function of age (effect modification), and how much of the
# observed age-mortality gradient is REMOVED by that modifiable policy (the
# flattening) -- the manipulation-safe substitute for mediation, since the
# manipulated quantity is the dose, not age or the height instrument.
#
# beta(age) = beta0 + beta1*age comes from a 2-endogenous 2SLS: exposure and
# exposure x age, instrumented by height and height x age (just-identified).
iv2sls_2 <- function(y, D, Z, X) {               # D,Z: n x k endogenous/instruments
  W <- cbind(D, X); Zf <- cbind(Z, X)
  What <- Zf %*% solve(crossprod(Zf), crossprod(Zf, W))
  beta <- solve(crossprod(What), crossprod(What, y))
  beta[seq_len(ncol(D)), 1]                       # endogenous coefs c(beta0, beta1)
}
AGE_CENTER <- 60; AGE_WIN <- 7.5                  # decade-centered; +/- yr smoothing window

age_curve_one <- function(var, label, vt_exp) {
  d <- analytic %>% filter(!is.na(.data[[var]])) %>% shift_frame(var, vt_exp, "feasible") %>%
    mutate(age_c = (age_at_admission - AGE_CENTER) / 10)
  fit_beta <- function(dd) {
    Xd <- model.matrix(as.formula(paste("~", iv_cov)), data = dd)
    D  <- cbind(dd[[var]], dd[[var]] * dd$age_c)
    # exposure x age depends on height x age^2 in the first stage, so the
    # interaction needs height x age^2 as a 3rd instrument (over-identified);
    # height x age alone under-identifies beta1.
    Z  <- cbind(dd[[INSTRUMENT]], dd[[INSTRUMENT]] * dd$age_c, dd[[INSTRUMENT]] * dd$age_c^2)
    iv2sls_2(dd$deceased, D, Z, Xd)
  }
  b <- fit_beta(d)
  ages <- seq(quantile(d$age_at_admission, 0.05), quantile(d$age_at_admission, 0.95), length.out = 40)
  local_mean <- function(dd, col, a) mean(dd[[col]][abs(dd$age_at_admission - a) <= AGE_WIN])
  point <- tibble(metric = label, age = ages) %>% rowwise() %>%
    mutate(beta_age = b[1] + b[2] * ((age - AGE_CENTER) / 10),
           mean_shift = local_mean(d, "shift", age),
           obs_mort   = local_mean(d, "event", age),
           excess_pp  = 100 * beta_age * mean_shift,
           cf_mort    = obs_mort - excess_pp / 100) %>% ungroup()
  boot <- matrix(NA_real_, N_BOOT, length(ages))
  for (k in seq_len(N_BOOT)) {
    bd <- d[sample.int(nrow(d), replace = TRUE), , drop = FALSE]
    bb <- tryCatch(fit_beta(bd), error = function(e) NULL); if (is.null(bb)) next
    ms <- vapply(ages, function(a) local_mean(bd, "shift", a), numeric(1))
    boot[k, ] <- 100 * (bb[1] + bb[2] * ((ages - AGE_CENTER) / 10)) * ms
  }
  point %>% mutate(excess_lo = apply(boot, 2, quantile, 0.025, na.rm = TRUE),
                   excess_hi = apply(boot, 2, quantile, 0.975, na.rm = TRUE))
}

if (nrow(policy_metrics) > 0) {
  age_curve_tbl <- pmap_dfr(policy_metrics %>% select(var, label, vt_exp),
    function(var, label, vt_exp) age_curve_one(var, label, vt_exp))
  write_csv(age_curve_tbl, file.path(final_dir, paste0("ivpolicy_age_curve_", site_name, ".csv")))

  ac <- age_curve_tbl %>% mutate(metric = factor(metric, levels = unique(metric)))
  # (A) excess pp vs age, with bootstrap ribbon and the young/old crossover
  pA <- ggplot(ac, aes(age, excess_pp)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_ribbon(aes(ymin = excess_lo, ymax = excess_hi, fill = metric), alpha = 0.2) +
    geom_line(aes(colour = metric), linewidth = 1) +
    scale_colour_manual(values = okabe[c(3, 1, 2)], name = NULL, aesthetics = c("colour", "fill")) +
    labs(x = "Age (years)", y = "Excess 60-day mortality, PBW vs PFVC-anchoring (pp)",
         title = "Age-varying PFVC-anchoring policy effect (modifiable dose; age as modifier)",
         subtitle = paste0(site_name, " - benefit rises with age; crossover where young are relatively under-dosed")) +
    theme_minimal(base_size = 10)
  # (B) age-mortality gradient: observed (PBW) vs counterfactual (PFVC-anchored)
  grad <- ac %>% select(metric, age, Observed = obs_mort, `PFVC-anchored` = cf_mort) %>%
    pivot_longer(c(Observed, `PFVC-anchored`), names_to = "dosing", values_to = "mort")
  pB <- ggplot(grad, aes(age, 100 * mort, colour = dosing, linetype = dosing)) +
    geom_line(linewidth = 1) + facet_wrap(~ metric) +
    scale_colour_manual(values = c(Observed = okabe[1], `PFVC-anchored` = okabe[3]), name = NULL) +
    scale_linetype_manual(values = c(Observed = 1, `PFVC-anchored` = 2), name = NULL) +
    labs(x = "Age (years)", y = "60-day mortality (%)",
         title = "Age-mortality gradient flattening under the modifiable PFVC-anchoring policy",
         subtitle = "Gap = age-mortality removable by changing the dose (no claim about the effect of age)") +
    theme_minimal(base_size = 10)
  ggsave(pdf_path("ivpolicy_age_curve"), pA / pB, width = 9, height = 8)
}

n_tbl <- 1L + (nrow(policy_tbl) > 0) + (nrow(policy_metrics) > 0) +
  (frac_offband > 0.10 && nrow(policy_metrics) > 0)               # diagnostic always
n_fig <- (nrow(policy_tbl) > 0) + (nrow(policy_metrics) > 0)
message("Wrote ", n_tbl, " table(s) + ", n_fig, " figure(s) to ", final_dir)
message("Script 09 complete.")
