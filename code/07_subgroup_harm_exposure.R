# =============================================================================
# Script 07: Does PBW-based dosing concentrate HARMFUL mechanical-stress
#            exposure in the bias-prone subgroups?
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# FEATURE BRANCH (feature/target-trial-emulation): NOT part of 00_run_pipeline.R.
# Standalone, in the explore_*.R / script-06 tradition — reads the script-03
# cross-sectional dataset and refits its own models.
#
# THE QUESTION ----------------------------------------------------------------
# Standard of care doses tidal volume per PREDICTED BODY WEIGHT (PBW). But PBW
# over-estimates lung size relative to predicted FVC (PFVC), and script 05 shows
# the over-sizing PBW/PFVC is *patterned*: larger in older, female, Black/Other,
# and (via PBW) shorter patients. So a constant, guideline-correct VT/PBW
# (6-8 mL/kg) delivers a LARGER fraction of true lung volume — and higher
# absolute mechanical stress (DP, MP, VT/PFVC) — in exactly those subgroups.
# Script 09 (height-IV) separately shows absolute stress causally raises mortality.
# This script JOINS the two halves and asks the subgroup question directly:
#
#   (A) HARMFUL-EXPOSURE BURDEN — among patients dosed by the book at
#       VT/PBW 6-8 mL/kg, what fraction cross the harmful stress threshold
#       (DP, MP, Ers, VT/PFVC), stratified by sex, race, age tertile, and
#       within-sex height tertile? The bias predicts: higher in female,
#       non-white, older, shorter patients.
#
#   (B) COUNTERFACTUAL RE-DOSING — for the same guideline-dosed patients,
#       contrast their actual PBW-anchored dose against a PFVC-anchored dose
#       (everyone re-dosed to the cohort-median VT/PFVC, holding the average
#       dose constant), map the VT change through lung mechanics to DP/MP, and
#       report the EXCESS predicted 60-day mortality of PBW- vs PFVC-dosing by
#       subgroup. This is the on-the-nose target-trial estimand: "PBW dosing
#       kills more, and the excess lands on subgroups X."
#
# KEY ASSUMPTIONS / LIMITATIONS (B) -------------------------------------------
#   * DP scales linearly with VT at fixed elastance (DP = Ers x VT); MP scales
#     ~ VT^2 at fixed RR/PEEP (MP ~ VT x DP). Partial-equilibrium: Ers, RR, PEEP
#     held fixed under re-dosing. Elastance is a lung property, unchanged by VT,
#     so Ers is NOT a re-dosing target (reported in A only).
#   * The mortality model is assumed causal in the stress metric conditional on
#     severity + demographics (same residual-confounding framing as the rest of
#     the project; E-values elsewhere bound the unmeasured channel).
#   * target VT/PFVC = cohort median is one PFVC-dosing rule (redistribution at
#     constant mean dose); a sensitivity target is reported.
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(survival)
library(splines)
library(EValue)
library(patchwork)

source("utils/config.R")
site_name <- config$site_name

output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

# Okabe-Ito (discrete); viridis for any continuous fill.
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7",
           "#F0E442", "#999999")

set.seed(20260616)

HORIZON      <- 60     # days, primary 60-day all-cause mortality
VT_PBW_LO    <- 6      # guideline-correct lung-protective dose window (mL/kg)
VT_PBW_HI    <- 8
MIN_CELL     <- 10     # CLIF rule: never report any group smaller than this
# Single interpretable, externally-anchored harm cut for VT/PFVC: tidal volume
# delivering >11% of predicted FVC. ~75th percentile of the low-VT (6 mL/kg PBW)
# arm of ARDSNet ARMA on individual-patient data. Fixed (NOT site-specific
# maxstat) so the harm cut and the protective re-dosing target mean the same
# thing at every site and pool coherently.
VTPFVC_TARGET <- 11
is_synthetic <- identical(site_name, "synthetic_clif")
N_BOOT_B     <- if (is_synthetic) 200L else 1000L  # re-dosing excess-mortality CIs

message("Site: ", site_name, " | synthetic = ", is_synthetic)

# =============================================================================
# 7a. Load cohort + re-anchor 60-day mortality at the index ventilation time
#     (synthetic-only workaround for the known-buggy synthetic-CLIF mortality
#     fields; real sites use death_dttm unchanged).
# =============================================================================
cross_sectional <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"))
message("Cohort: ", nrow(cross_sectional), " patients")

# *** SYNTHETIC-ONLY MORTALITY WORKAROUND (synthetic-CLIF mortality is buggy). REMOVE once
# the synthetic-CLIF mortality fix lands. Real sites use death_dttm unchanged.
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
          "long-tailed survival (synthetic CLIF mortality is unreliable). ***")
  set.seed(20260615)
  n <- nrow(cross_sectional)
  died60 <- rbinom(n, 1L, 0.35)
  tte    <- rep(NA_real_, n)
  n_dec  <- sum(died60 == 1L)
  tte[died60 == 1L] <- rtrunc_lnorm(n_dec, meanlog = log(9), sdlog = 0.95,
                                    lo = 0.04, hi = HORIZON)
  cross_sectional <- cross_sectional %>%
    mutate(event = as.integer(died60),
           time  = if_else(died60 == 1L, tte, as.numeric(HORIZON)),
           deceased = as.integer(died60))
} else {
  cross_sectional <- cross_sectional %>%
    mutate(idx_to_death = as.numeric(difftime(death_dttm, recorded_dttm, units = "days")),
           event = as.integer(!is.na(idx_to_death) & idx_to_death >= 0 &
                                idx_to_death <= HORIZON),
           time  = if_else(event == 1L, idx_to_death, as.numeric(HORIZON)),
           deceased = event)
}
cross_sectional <- cross_sectional %>% mutate(time = pmax(time, 1 / 24))

message("Re-anchored 60-day mortality: ", sum(cross_sectional$event), " / ",
        nrow(cross_sectional), " (",
        round(100 * mean(cross_sectional$event), 1), "%)")

# =============================================================================
# 7b. Metric registry + harmful thresholds (raw absolute stress)
# =============================================================================
# Harmful exposure is defined on the dose-as-fraction-of-true-lung-size (VT/PFVC)
# and on the stress the lung experiences (raw DP/MP/Ers) — NOT on a PBW/PFVC-
# normalized metric, because the question is precisely whether the raw physiologic
# insult is patterned by demographics.
#
# The burden-by-subgroup analysis reports the headline size-relative dose VT/PFVC
# at the FIXED clinical cut (VTPFVC_TARGET) — interpretable and cross-site-poolable.
# The size-relative-vs-extensive mechanic comparators (DP/MP/Ers) and their
# discrimination have moved to the physiology/normalization analysis (scripts 05/08);
# this removes the dependency on script 06's (retired) maxstat cutpoints. The
# counterfactual re-dosing + E-value models below still use raw DP/MP as
# absolute-stress exposures — that is independent of any harm threshold.
metric_specs <- tribble(
  ~var,     ~label,    ~unit,  ~scale_type,
  "vtpfvc", "VT/PFVC", "mL/L", "size-relative"
) %>% mutate(threshold = VTPFVC_TARGET)

message("Harmful threshold (VT/PFVC = fixed clinical cut):")
walk(seq_len(nrow(metric_specs)), ~ message(
  "  ", metric_specs$label[.x], " > ", round(metric_specs$threshold[.x], 3),
  " ", metric_specs$unit[.x]))

# =============================================================================
# 7c. Guideline-dosed analytic frame + subgroup factors
# =============================================================================
# Whole-cohort age tertiles (comparable strata); within-sex height tertiles so
# "shorter patients" is isolated from the sex effect (height is sexually
# dimorphic, so pooled height tertiles would just re-encode sex).
age_breaks <- quantile(cross_sectional$age_at_admission, c(1/3, 2/3), na.rm = TRUE)

analytic <- cross_sectional %>%
  filter(!is.na(vtpbw), vtpbw >= VT_PBW_LO, vtpbw <= VT_PBW_HI,
         !is.na(sex_category), !is.na(race_category),
         !is.na(age_at_admission), !is.na(height_cm),
         !is.na(vtpfvc), !is.na(pbw), !is.na(pfvc),
         !is.na(sofa_total), !is.na(sf_ratio), !is.na(bmi),
         !is.na(time), !is.na(event)) %>%
  group_by(sex_category) %>%
  mutate(height_z = as.numeric(scale(height_cm))) %>%
  ungroup() %>%
  mutate(
    age_grp = cut(age_at_admission, breaks = c(-Inf, age_breaks, Inf),
                  labels = c("Young", "Middle", "Old")),
    height_grp = cut(height_z,
                     breaks = c(-Inf, quantile(height_z, c(1/3, 2/3), na.rm = TRUE), Inf),
                     labels = c("Short", "Middle", "Tall")),
    sex_grp  = factor(sex_category),
    race_grp = factor(race_category)
  )

message("Guideline-dosed (VT/PBW ", VT_PBW_LO, "-", VT_PBW_HI, "): ",
        nrow(analytic), " of ", nrow(cross_sectional))

# Subgroup variables with their bias-prone (reference-contrast) ordering. The
# first level is the LOW-bias reference; later levels are the predicted-harmed.
subgroup_vars <- tribble(
  ~svar,        ~slabel,                 ~ref_level,
  "sex_grp",    "Sex",                   "Male",
  "race_grp",   "Race",                  "WHITE",
  "age_grp",    "Age tertile",           "Young",
  "height_grp", "Height tertile (w/in sex)", "Tall"
)

# =============================================================================
# 7d. (A) Harmful-exposure burden by subgroup
# =============================================================================
wilson_ci <- function(k, n) {                       # 95% Wilson interval
  if (n == 0) return(c(NA_real_, NA_real_))
  z <- 1.959964; p <- k / n
  ctr <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(ctr - hw, ctr + hw)
}

# Crude fraction-over-threshold per subgroup level, with cells < MIN_CELL
# suppressed (CLIF rule). Returns one row per (metric x subgroup-var x level).
burden_one <- function(df, var, lbl, thr, svar, slabel) {
  df %>%
    filter(!is.na(.data[[var]])) %>%             # DP/MP missing w/o plateau pressure
    mutate(over = as.integer(.data[[var]] > thr)) %>%
    group_by(level = .data[[svar]]) %>%
    summarise(n = n(), n_over = sum(over), .groups = "drop") %>%
    filter(n >= MIN_CELL) %>%
    rowwise() %>%
    mutate(frac_over = n_over / n,
           ci_lo = wilson_ci(n_over, n)[1],
           ci_hi = wilson_ci(n_over, n)[2]) %>%
    ungroup() %>%
    mutate(metric = lbl, threshold = thr, subgroup = slabel,
           level = as.character(level), .before = 1)
}

# Severity-adjusted disparity: logistic over_threshold ~ subgroup + severity,
# OR per non-reference level vs the bias-low reference. Shows the gap survives
# (or not) adjustment for measured severity (SOFA, S/F, BMI).
burden_adjusted <- function(df, var, lbl, thr, svar, slabel, ref_level) {
  d <- df %>%
    filter(!is.na(.data[[var]])) %>%
    mutate(over = as.integer(.data[[var]] > thr),
           grp = relevel(factor(.data[[svar]]), ref = ref_level))
  if (length(unique(d$over)) < 2 || nlevels(droplevels(d$grp)) < 2)
    return(tibble())
  fit <- tryCatch(
    glm(over ~ grp + sofa_total + sf_ratio + bmi, data = d, family = binomial),
    error = function(e) NULL)
  if (is.null(fit)) return(tibble())
  s <- summary(fit)$coefficients
  ci <- suppressMessages(confint.default(fit))
  grp_rows <- grep("^grp", rownames(s), value = TRUE)
  tibble(metric = lbl, subgroup = slabel,
         level = sub("^grp", "", grp_rows),
         adj_or = exp(s[grp_rows, "Estimate"]),
         adj_lo = exp(ci[grp_rows, 1]),
         adj_hi = exp(ci[grp_rows, 2]),
         adj_p  = s[grp_rows, "Pr(>|z|)"])
}

burden_tbl <- pmap_dfr(metric_specs, function(var, label, unit, scale_type, threshold) {
  pmap_dfr(subgroup_vars, function(svar, slabel, ref_level) {
    crude <- burden_one(analytic, var, label, threshold, svar, slabel)
    adj   <- burden_adjusted(analytic, var, label, threshold, svar, slabel, ref_level)
    crude %>% left_join(adj, by = c("metric", "subgroup", "level"))
  })
})

write_csv(burden_tbl,
          file.path(final_dir, paste0("harm_burden_by_subgroup_", site_name, ".csv")))

# =============================================================================
# 7e. (B) Counterfactual re-dosing: PBW-anchored (actual) vs PFVC-anchored
# =============================================================================
# Re-dose every patient to a COMMON fraction of true lung size (a target VT/PFVC),
# map the resulting VT change through lung mechanics to DP and MP, and push both
# the actual and counterfactual stress through a 60-day mortality model. The
# excess (actual - counterfactual) predicted death is the mortality attributable
# to PBW- vs PFVC-anchored dosing.
#
# THREE TARGETS (sensitivity to the dosing rule):
#   * "median" = cohort-median VT/PFVC. Holds MEAN dose fixed, so this is a pure
#                REDISTRIBUTION (net population delta ~ 0 by construction); it
#                isolates the EQUITY question (who currently bears the harm).
#   * "p25"    = 25th-percentile VT/PFVC. A distributional PROTECTIVE target that
#                lowers most patients' dose -> breaks the zero-sum and shows NET
#                benefit (site-specific value, defined identically per site).
#   * "arma11" = the FIXED clinical cut, VT/PFVC < 11% of predicted FVC (ARMA LTVV
#                ~p75). The headline protective rule: one number, identical across
#                sites, so the protective arm pools coherently and the harm cut (A)
#                and re-dosing target (B) are the same line.
# The "Overall / All" row in every block is the NET POPULATION delta.
B_targets <- tibble(
  target_name   = c("median", "p25", "arma11"),
  target_vtpfvc = c(median(analytic$vtpfvc, na.rm = TRUE),
                    unname(quantile(analytic$vtpfvc, 0.25, na.rm = TRUE)),
                    VTPFVC_TARGET)
)
message("Re-dosing targets (VT/PFVC mL/L): ",
        paste(sprintf("%s=%.2f", B_targets$target_name, B_targets$target_vtpfvc),
              collapse = ", "))

# DP = Ers x VT (linear in VT); MP ~ VT x DP (~ VT^2) at fixed RR/PEEP/Ers.
B_metrics <- tribble(
  ~var,               ~label,  ~vt_exp,
  "dp",               "DP",    1,
  "mechanical_power", "MP",    2
)
# Mortality model: 60-day death ~ absolute stress metric + severity + demo.
# Demographics ARE included (the metric is the exposure; demographics confound the
# metric->death link, they are not mediators of the dose change).
# Age enters as a natural spline (df=3), not linearly: age's NON-mechanical paths
# to mortality (frailty, reserve, withdrawal of care) are nonlinear, and nonlinear
# age is the only enumerable confounding channel here -- so we enumerate it fully
# rather than leave curvature as residual. ns() knots are stored in the fit, so
# predict() on the counterfactual frame uses the same basis.
mort_cov <- "sofa_total + sf_ratio + bmi + ns(age_at_admission, df = 3) + sex_category + race_category"

# Subgroup means PLUS an "Overall / All" net-population row. Expects .delta and
# .vt_drop already attached to d.
excess_by_sub <- function(d) {
  sub <- map_dfr(seq_len(nrow(subgroup_vars)), function(i) {
    sv <- subgroup_vars$svar[i]; sl <- subgroup_vars$slabel[i]
    d %>% group_by(level = .data[[sv]]) %>%
      summarise(n = n(), mean_delta = mean(.delta),
                mean_vt_drop_ml = mean(.vt_drop), .groups = "drop") %>%
      filter(n >= MIN_CELL) %>%
      mutate(subgroup = sl, level = as.character(level))
  })
  overall <- d %>% summarise(n = n(), mean_delta = mean(.delta),
                             mean_vt_drop_ml = mean(.vt_drop)) %>%
    mutate(subgroup = "Overall", level = "All")
  bind_rows(sub, overall)
}

# Excess predicted mortality for one metric under one re-dosing target; percentile
# bootstrap CI (refit each rep so model uncertainty is included).
excess_one <- function(df, var, vt_exp, lbl, target, target_name) {
  df <- df %>%
    filter(!is.na(.data[[var]])) %>%                  # metric-complete
    mutate(.vt_ratio = target / vtpfvc,               # VT_cf / VT_actual
           .cf = .data[[var]] * .vt_ratio^vt_exp,      # mechanics map to DP/MP
           .vt_drop = (1 - .vt_ratio) * vtpfvc * pfvc) # mL removed (>0 = de-escalated)
  fit <- glm(as.formula(paste("event ~", var, "+", mort_cov)),
             data = df, family = binomial)
  base <- df %>% mutate(.p_actual = predict(fit, newdata = df, type = "response"))
  cf_df <- df; cf_df[[var]] <- df$.cf
  base$.p_cf  <- predict(fit, newdata = cf_df, type = "response")
  base$.delta <- base$.p_actual - base$.p_cf           # >0 => PBW dosing kills more
  point <- excess_by_sub(base)

  acc <- vector("list", N_BOOT_B)
  for (b in seq_len(N_BOOT_B)) {
    idx <- sample.int(nrow(df), replace = TRUE)
    bd  <- df[idx, , drop = FALSE]
    fb  <- tryCatch(glm(as.formula(paste("event ~", var, "+", mort_cov)),
                        data = bd, family = binomial), error = function(e) NULL)
    if (is.null(fb)) next
    bd$.p_actual <- predict(fb, newdata = bd, type = "response")
    cf_b <- bd; cf_b[[var]] <- bd$.cf
    bd$.p_cf  <- predict(fb, newdata = cf_b, type = "response")
    bd$.delta <- bd$.p_actual - bd$.p_cf
    acc[[b]] <- excess_by_sub(bd) %>% select(subgroup, level, mean_delta)
  }
  ci <- bind_rows(acc) %>% group_by(subgroup, level) %>%
    summarise(delta_lo = quantile(mean_delta, 0.025, na.rm = TRUE),
              delta_hi = quantile(mean_delta, 0.975, na.rm = TRUE),
              .groups = "drop")
  point %>% left_join(ci, by = c("subgroup", "level")) %>%
    mutate(metric = lbl, target_name = target_name, target_vtpfvc = target,
           .before = 1)
}

excess_tbl <- pmap_dfr(B_metrics, function(var, label, vt_exp)
  pmap_dfr(B_targets, function(target_name, target_vtpfvc)
    excess_one(analytic, var, vt_exp, label, target_vtpfvc, target_name)))

write_csv(excess_tbl,
          file.path(final_dir, paste0("harm_redose_excess_mortality_", site_name, ".csv")))

# --- E-value on the per-SD stress -> mortality association (DP, MP) ------------
# Quantifies how strong an UNMEASURED confounder would have to be -- associated
# with both the stress metric and mortality, above the measured set -- to explain
# away the adjusted association. The dominant such residual here is code-status /
# withdrawal of life support (age-correlated, not in CLIF for this cohort); the
# E-value bounds it without measuring it. Non-rare OR transform (~35% mortality).
evalue_one <- function(df, var, lbl) {
  d <- df %>% filter(!is.na(.data[[var]])) %>%
    mutate(.z = as.numeric(scale(.data[[var]])))
  fit <- glm(as.formula(paste("event ~ .z +", mort_cov)), data = d, family = binomial)
  s   <- summary(fit)$coefficients[".z", ]
  ci  <- suppressMessages(confint.default(fit))[".z", ]
  or  <- exp(unname(s["Estimate"])); lo <- exp(unname(ci[1])); hi <- exp(unname(ci[2]))
  ev  <- EValue::evalues.OR(est = or, lo = lo, hi = hi, rare = FALSE)
  ci_col <- if (or > 1) "lower" else "upper"      # E-value for the CI bound nearer null
  tibble(metric = lbl, or_per_sd = or, or_lo = lo, or_hi = hi,
         evalue_point = ev["E-values", "point"],
         evalue_ci    = ev["E-values", ci_col])
}
evalue_tbl <- pmap_dfr(B_metrics, function(var, label, vt_exp) evalue_one(analytic, var, label))
write_csv(evalue_tbl, file.path(final_dir, paste0("harm_evalue_", site_name, ".csv")))
message("E-value (per-SD stress->mortality): ",
        paste(sprintf("%s OR=%.2f point=%.2f CI=%.2f", evalue_tbl$metric,
                      evalue_tbl$or_per_sd, evalue_tbl$evalue_point, evalue_tbl$evalue_ci),
              collapse = " | "))

# =============================================================================
# 7f. Positivity diagnostics for the single VT/PFVC < 11% threshold
# =============================================================================
# The fixed threshold makes positivity AUDITABLE, not absent. VT/PFVC is near-
# deterministic in demographics (PFVC = f(age,sex,race,height)) within the 6-8
# VT/PBW window, so conditioning on demographics COLLAPSES the exposure contrast
# (propensity -> 0/1). We document that deliberately, and quantify how much of the
# counterfactual (B) rests on extrapolation beyond observed support.
treat_below <- as.integer(analytic$vtpfvc < VTPFVC_TARGET)   # 1 = protective achieved

auc_fn <- function(score, y) {                                # Mann-Whitney AUC, no deps
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(score)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# (1) Propensity overlap: the positivity-respecting set (severity + the clinical
# dose lever) vs the demographics-augmented set that collapses the contrast.
pos_specs <- tribble(
  ~spec,               ~cov_str,
  "severity_dose",     "sofa_total + sf_ratio + bmi + vtpbw",
  "plus_demographics", "sofa_total + sf_ratio + bmi + vtpbw + age_at_admission + sex_category + race_category + height_cm"
)
ps_long <- list(); pos_overlap <- list()
for (i in seq_len(nrow(pos_specs))) {
  fit <- glm(as.formula(paste("y ~", pos_specs$cov_str[i])),
             data = analytic %>% mutate(y = treat_below), family = binomial)
  ps     <- as.numeric(predict(fit, type = "response"))
  p_marg <- mean(treat_below)
  sw     <- ifelse(treat_below == 1, p_marg / ps, (1 - p_marg) / (1 - ps))
  ess    <- sum(sw)^2 / sum(sw^2)
  pos_overlap[[i]] <- tibble(
    spec = pos_specs$spec[i], n = length(ps), n_protective = sum(treat_below),
    ps_min = min(ps), ps_max = max(ps),
    frac_ps_under05 = mean(ps < 0.05), frac_ps_over95 = mean(ps > 0.95),
    frac_extreme = mean(ps < 0.05 | ps > 0.95),
    ess = ess, ess_frac = ess / length(ps), auc = auc_fn(ps, treat_below))
  ps_long[[i]] <- tibble(spec = pos_specs$spec[i], ps = ps,
                         arm = ifelse(treat_below == 1, "VT/PFVC <11%", "VT/PFVC >=11%"))
}
pos_overlap_tbl <- bind_rows(pos_overlap)
ps_long         <- bind_rows(ps_long)
write_csv(pos_overlap_tbl, file.path(final_dir, paste0("pos_overlap_", site_name, ".csv")))

# (2) Per-stratum support: fraction above/below 11% in each demographic cell.
# Cells ~all on one side carry no within-cell contrast (identified only by the
# model, not by data) -> flagged near_deterministic.
pos_support_tbl <- pmap_dfr(subgroup_vars, function(svar, slabel, ref_level) {
  analytic %>%
    mutate(below = as.integer(vtpfvc < VTPFVC_TARGET)) %>%
    group_by(level = .data[[svar]]) %>%
    summarise(n = n(), n_below = sum(below), .groups = "drop") %>%
    filter(n >= MIN_CELL) %>%
    mutate(subgroup = slabel, level = as.character(level),
           frac_below = n_below / n, frac_above = 1 - frac_below,
           near_deterministic = pmin(frac_below, frac_above) < 0.05)
})
write_csv(pos_support_tbl, file.path(final_dir, paste0("pos_stratum_support_", site_name, ".csv")))

# (3) Counterfactual-move / extrapolation audit (for B). Within demographic cells
# (age x sex x race x within-sex height tertile), does the 11% target fall inside
# the observed VT/PFVC support [p5,p95]? Patients whose target sits outside are
# moved by MODEL EXTRAPOLATION, not by data; small cells are unestimable.
cells <- analytic %>%
  group_by(age_grp, sex_grp, race_grp, height_grp) %>%
  summarise(cell_n = n(),
            vt_p05 = quantile(vtpfvc, 0.05, na.rm = TRUE),
            vt_p95 = quantile(vtpfvc, 0.95, na.rm = TRUE), .groups = "drop")
analytic_sup <- analytic %>%
  left_join(cells, by = c("age_grp", "sex_grp", "race_grp", "height_grp")) %>%
  mutate(small_cell  = cell_n < MIN_CELL,
         supported   = !small_cell & VTPFVC_TARGET >= vt_p05 & VTPFVC_TARGET <= vt_p95,
         vt_drop_ml  = (vtpfvc - VTPFVC_TARGET) * pfvc)   # same convention as 7e

extrap_summ <- function(d) {
  d %>% summarise(n = n(),
                  frac_supported    = mean(supported),
                  frac_extrapolated = mean(!supported & !small_cell),
                  frac_small_cell   = mean(small_cell),
                  frac_deescalated  = mean(vtpfvc > VTPFVC_TARGET),
                  mean_vt_drop_ml   = mean(vt_drop_ml), .groups = "drop")
}
pos_extrap_overall <- extrap_summ(analytic_sup) %>%
  mutate(subgroup = "Overall", level = "All", .before = 1)
pos_extrap_sub <- pmap_dfr(subgroup_vars, function(svar, slabel, ref_level) {
  analytic_sup %>% group_by(level = .data[[svar]]) %>% extrap_summ() %>%
    filter(n >= MIN_CELL) %>%
    mutate(subgroup = slabel, level = as.character(level), .before = 1)
})
pos_extrap_tbl <- bind_rows(pos_extrap_overall, pos_extrap_sub)
write_csv(pos_extrap_tbl, file.path(final_dir, paste0("pos_extrapolation_audit_", site_name, ".csv")))

message("Positivity: +demographics frac_extreme = ",
        round(pos_overlap_tbl$frac_extreme[pos_overlap_tbl$spec == "plus_demographics"], 3),
        " (vs severity-only ",
        round(pos_overlap_tbl$frac_extreme[pos_overlap_tbl$spec == "severity_dose"], 3),
        "); ESS frac = ",
        round(pos_overlap_tbl$ess_frac[pos_overlap_tbl$spec == "plus_demographics"], 3),
        "; B extrapolation frac = ", round(pos_extrap_overall$frac_extrapolated, 3))

# =============================================================================
# 7g. Figures
# =============================================================================
pdf_path <- function(stub) file.path(final_dir, paste0(stub, "_", site_name, ".pdf"))
sub_levels <- c("Male", "Female", "WHITE", "BLACK", "OTHER",
                "Young", "Middle", "Old", "Tall", "Short")  # for stable ordering
order_lvl <- function(x) factor(x, levels = intersect(sub_levels, unique(x)))

# --- (A) burden: fraction over harmful threshold, metric x subgroup grid -------
if (nrow(burden_tbl) > 0) {
  bt <- burden_tbl %>%
    mutate(metric = factor(metric, levels = metric_specs$label),
           subgroup = factor(subgroup, levels = subgroup_vars$slabel),
           level = order_lvl(level))
  p_burden <- ggplot(bt, aes(level, frac_over, fill = subgroup)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
    facet_grid(metric ~ subgroup, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = okabe[c(1, 2, 3, 4)], guide = "none") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = NULL,
         y = "Fraction over harmful threshold (95% Wilson CI)",
         title = paste0("Harmful mechanical-stress exposure among guideline-dosed (VT/PBW ",
                        VT_PBW_LO, "-", VT_PBW_HI, ") patients"),
         subtitle = paste0(site_name,
           if (is_synthetic) " (SYNTHETIC - simulated survival)" else "",
           " - size-relative VT/PFVC harm is higher in Female, non-White, Older, ",
           "and Shorter patients (the demographic bias of a fixed VT/PBW dose)")) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(pdf_path("harm_burden_by_subgroup"), p_burden, width = 11, height = 8)
}

# --- (B) counterfactual excess mortality, 3 re-dosing targets + net population --
if (nrow(excess_tbl) > 0) {
  sub_order <- c(subgroup_vars$slabel, "Overall")
  et <- excess_tbl %>%
    mutate(metric = factor(metric, levels = B_metrics$label),
           subgroup = factor(subgroup, levels = sub_order),
           level = factor(level, levels = c(sub_levels, "All")),
           target_name = factor(target_name, levels = B_targets$target_name))
  p_excess <- ggplot(et, aes(level, 100 * mean_delta, fill = target_name)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = 100 * delta_lo, ymax = 100 * delta_hi),
                  position = position_dodge(width = 0.8), width = 0.25) +
    facet_grid(metric ~ subgroup, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c(median = okabe[5], p25 = okabe[3], arma11 = okabe[4]),
                      name = "PFVC re-dosing target",
                      labels = c(median = "median (redistribution)",
                                 p25 = "25th pctile (protective)",
                                 arma11 = "VT/PFVC <11% (ARMA LTVV, fixed)")) +
    labs(x = NULL,
         y = "Excess 60-day mortality, PBW minus PFVC dosing (pct points)",
         title = "Counterfactual re-dosing: mortality attributable to PBW (vs PFVC) dosing",
         subtitle = paste0(site_name,
           if (is_synthetic) " (SYNTHETIC - simulated survival)" else "",
           " - >0 = PBW dosing harms (PFVC better). 'Overall/All' = net population delta: ",
           "median target nets ~0 (redistribution), protective targets net > 0 (PFVC saves lives)")) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "top")
  ggsave(pdf_path("harm_redose_excess_mortality"), p_excess, width = 12, height = 6.5)
}

# --- (7f-1) propensity overlap: severity+dose vs +demographics (collapse) ------
if (nrow(ps_long) > 0) {
  pl <- ps_long %>%
    mutate(spec = factor(spec, levels = pos_specs$spec,
                         labels = c("severity + dose\n(positivity-respecting)",
                                    "+ demographics\n(collapses contrast)")))
  p_ps <- ggplot(pl, aes(ps, fill = arm, colour = arm)) +
    geom_density(alpha = 0.4) +
    facet_wrap(~ spec) +
    geom_vline(xintercept = c(0.05, 0.95), linetype = 3, colour = "grey50") +
    scale_fill_manual(values = c("VT/PFVC <11%" = okabe[3], "VT/PFVC >=11%" = okabe[1]),
                      name = NULL) +
    scale_colour_manual(values = c("VT/PFVC <11%" = okabe[3], "VT/PFVC >=11%" = okabe[1]),
                        name = NULL) +
    labs(x = "Propensity to achieve VT/PFVC <11%", y = "Density",
         title = "Positivity: propensity overlap for the VT/PFVC <11% strategy",
         subtitle = paste0(site_name,
           " - adding demographics piles mass at 0/1 (positivity failure); this is WHY ",
           "demographics are not balanced on the VT/PFVC contrast")) +
    theme_minimal(base_size = 10)
  ggsave(pdf_path("pos_overlap"), p_ps, width = 10, height = 4.5)
}

# --- (7f-2) per-stratum support: fraction over the 11% harm cut ----------------
if (nrow(pos_support_tbl) > 0) {
  st <- pos_support_tbl %>%
    mutate(subgroup = factor(subgroup, levels = subgroup_vars$slabel),
           level = order_lvl(level))
  p_sup <- ggplot(st, aes(level, frac_above)) +
    geom_col(aes(fill = near_deterministic), width = 0.7) +
    geom_hline(yintercept = c(0.05, 0.95), linetype = 3, colour = "grey50") +
    facet_grid(~ subgroup, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c(`FALSE` = okabe[2], `TRUE` = okabe[5]),
                      name = "near-deterministic\n(<5% either side)") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Fraction >= VT/PFVC 11% (over harm cut)",
         title = "Per-stratum support for the VT/PFVC <11% contrast",
         subtitle = paste0(site_name,
           " - cells near 0% or 100% carry no within-cell contrast (identified by model only)")) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(pdf_path("pos_stratum_support"), p_sup, width = 11, height = 4)
}

message("Wrote 6 tables + 4 figures to ", final_dir)
message("Script 07 complete.")
