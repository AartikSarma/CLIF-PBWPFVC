# =============================================================================
# Script 08: Dose heterogeneity under PBW dosing + VT/PFVC as the structural dose
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# FEATURE BRANCH (feature/target-trial-emulation): NOT part of 00_run_pipeline.R.
# Standalone; reads the script-03 cross-sectional dataset.
#
# THE CLAIM (prong-2 PREMISE) -------------------------------------------------
# Patients given the "same treatment" (a guideline-protective VT/PBW of 6-8
# mL/kg) receive WILDLY DIFFERENT physiologic doses (VT as a fraction of
# predicted FVC), because PBW/PFVC mis-sizing varies by body size/demographics.
# This script shows, in observational CLIF data:
#   (1) DECOMPOSITION   - at fixed VT/PBW, the variance in the delivered dose
#                         (VT/PFVC) is dominated by mis-sizing (PBW/PFVC), NOT by
#                         the clinician's protocol setting.
#   (2) WITHIN-BAND SPREAD - the VT/PFVC distribution at a fixed VT/PBW is wide
#                         and demographically shifted (Female/non-White/Older/
#                         Shorter dosed higher as a fraction of lung size).
#   (3) STRUCTURAL-DOSE INVARIANCE - the VT/PFVC -> mortality dose-response is
#                         HOMOGENEOUS across demographic strata (and, run at both
#                         sites, across cohorts). An invariant effect across
#                         environments is the signature of the STRUCTURAL dose.
#
# SCOPE NOTE: CLIF is fixed at VT/PBW 6-8 (no randomized dose contrast), so the
# *consequence* for trials -- heterogeneity of the randomized VT/PBW treatment
# effect -- is NOT testable here; it belongs to the separate Bayesian RCT
# analysis. This script establishes only the PREMISE. The invariance shown here
# is also the identifying assumption the height-IV policy (script 09) relies on
# to generalize from height-driven mis-sizing to full PFVC-anchoring.
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(splines)
library(patchwork)

source("utils/config.R")
site_name <- config$site_name

output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

okabe <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7",
           "#F0E442", "#999999")
set.seed(20260616)

HORIZON       <- 60
MIN_CELL      <- 10
VTPFVC_TARGET <- 11        # VT/PFVC <11% of predicted FVC (ARMA LTVV ~p75); ref line
is_synthetic  <- identical(site_name, "synthetic_clif")
message("Site: ", site_name, " | synthetic = ", is_synthetic)

# =============================================================================
# 8a. Load + re-anchor 60-day mortality (identical handling to scripts 06/07,
#     incl. the synthetic-only mortality workaround).
# =============================================================================
cross_sectional <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"))
message("Cohort: ", nrow(cross_sectional), " patients")

rtrunc_lnorm <- function(n_needed, meanlog, sdlog, lo, hi) {
  acc <- numeric(0)
  while (length(acc) < n_needed) {
    cand <- rlnorm(max(n_needed * 2L, 1000L), meanlog, sdlog)
    cand <- cand[cand > lo & cand <= hi]; acc <- c(acc, cand)
  }
  acc[seq_len(n_needed)]
}
if (is_synthetic) {
  message("*** SYNTHETIC SITE: simulated long-tailed survival (see script 06 header). ***")
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
# 8b. Analytic frame + subgroup factors (same construction as script 07)
# =============================================================================
age_breaks <- quantile(cross_sectional$age_at_admission, c(1/3, 2/3), na.rm = TRUE)
analytic <- cross_sectional %>%
  filter(!is.na(vtpbw), !is.na(vtpfvc), !is.na(pbwpfvc), !is.na(pbw), !is.na(pfvc),
         !is.na(sex_category), !is.na(race_category), !is.na(age_at_admission),
         !is.na(height_cm), !is.na(sofa_total), !is.na(sf_ratio), !is.na(bmi),
         !is.na(time), !is.na(event), vtpbw > 0, vtpfvc > 0, pbwpfvc > 0) %>%
  group_by(sex_category) %>% mutate(height_z = as.numeric(scale(height_cm))) %>% ungroup() %>%
  mutate(age_grp = cut(age_at_admission, c(-Inf, age_breaks, Inf),
                       labels = c("Young", "Middle", "Old")),
         height_grp = cut(height_z, c(-Inf, quantile(height_z, c(1/3, 2/3), na.rm = TRUE), Inf),
                          labels = c("Short", "Middle", "Tall")),
         sex_grp = factor(sex_category), race_grp = factor(race_category))

vtpbw_rng <- range(analytic$vtpbw)
message("Analytic n = ", nrow(analytic), " | observed VT/PBW range: ",
        round(vtpbw_rng[1], 2), "-", round(vtpbw_rng[2], 2),
        if (diff(vtpbw_rng) < 2.5) "  (narrow: cohort is protocol-restricted)" else "")

subgroup_vars <- tribble(
  ~svar,        ~slabel,                     ~ref_level,
  "sex_grp",    "Sex",                       "Male",
  "race_grp",   "Race",                      "WHITE",
  "age_grp",    "Age tertile",               "Young",
  "height_grp", "Height tertile (w/in sex)", "Tall"
)

# =============================================================================
# 8c. (1) Variance decomposition of the delivered dose
# =============================================================================
# Identity (script 03): VT/PFVC = (VT/PBW) x (PBW/PFVC) x 0.1, so in logs
#   Var(log VT/PFVC) = Var(log VT/PBW) + Var(log PBW/PFVC) + 2 Cov(.,.).
# The PBW/PFVC (mis-sizing) share is the fraction of delivered-dose variance the
# clinician does NOT control at a fixed protocol.
decomp_one <- function(df, lbl) {
  lv <- log(df$vtpfvc); lp <- log(df$vtpbw); lm <- log(df$pbwpfvc)
  vtot <- var(lp + lm)                       # = Var(log VT/PFVC) up to the 0.1 constant
  tibble(cohort = lbl, n = nrow(df),
         var_total = vtot,
         var_protocol  = var(lp), var_missizing = var(lm), cov2 = 2 * cov(lp, lm),
         share_protocol  = var(lp) / vtot,
         share_missizing = var(lm) / vtot,
         share_cov       = 2 * cov(lp, lm) / vtot,
         identity_max_abs_dev = max(abs(lv - (lp + lm - log(10)))))   # sanity check
}
decomp_tbl <- bind_rows(
  decomp_one(analytic, "guideline_6_8"),
  decomp_one(cross_sectional %>% filter(vtpbw > 0, vtpfvc > 0, pbwpfvc > 0,
                                        !is.na(vtpbw), !is.na(vtpfvc), !is.na(pbwpfvc)),
             "all_with_dose"))
write_csv(decomp_tbl, file.path(final_dir, paste0("dose_variance_decomposition_", site_name, ".csv")))
message("Dose variance from mis-sizing (guideline 6-8): ",
        round(100 * decomp_tbl$share_missizing[decomp_tbl$cohort == "guideline_6_8"], 1), "%")

# =============================================================================
# 8d. (2) Within-band VT/PFVC spread by demographic subgroup
# =============================================================================
spread_tbl <- pmap_dfr(subgroup_vars, function(svar, slabel, ref_level) {
  analytic %>% group_by(level = .data[[svar]]) %>%
    summarise(n = n(),
              vtpfvc_median = median(vtpfvc), vtpfvc_p5 = quantile(vtpfvc, 0.05),
              vtpfvc_p95 = quantile(vtpfvc, 0.95),
              frac_over_target = mean(vtpfvc > VTPFVC_TARGET), .groups = "drop") %>%
    filter(n >= MIN_CELL) %>%
    mutate(spread_ratio_p95_p5 = vtpfvc_p95 / vtpfvc_p5,
           subgroup = slabel, level = as.character(level))
})
write_csv(spread_tbl, file.path(final_dir, paste0("dose_within_band_spread_", site_name, ".csv")))

# =============================================================================
# 8e. (3) VT/PFVC dose-response + homogeneity across strata (structural dose)
# =============================================================================
# Per-stratum slope of 60-day mortality on VT/PFVC (global-SD units, severity-
# adjusted) + a pooled interaction LRT. A HOMOGENEOUS slope across strata (null
# interaction) is the invariance signature that VT/PFVC is the structural dose.
# NB: this tests SLOPE HOMOGENEITY, which is estimable within strata; it is NOT
# the demographics-adjusted causal magnitude (unidentifiable by positivity --
# that effect is script 09's height-IV job).
sev <- "sofa_total + sf_ratio + bmi"
analytic <- analytic %>% mutate(z_vtpfvc = as.numeric(scale(vtpfvc)))

stratum_slopes <- pmap_dfr(subgroup_vars, function(svar, slabel, ref_level) {
  d <- analytic %>% mutate(stratum = .data[[svar]])
  # pooled interaction LRT (does the VT/PFVC slope differ across this axis?)
  full <- glm(as.formula(paste("event ~ z_vtpfvc * stratum +", sev)), data = d, family = binomial)
  red  <- glm(as.formula(paste("event ~ z_vtpfvc + stratum +", sev)), data = d, family = binomial)
  lrt_p <- anova(red, full, test = "LRT")$`Pr(>Chi)`[2]
  # per-stratum slope (per global SD of VT/PFVC)
  d %>% group_by(level = stratum) %>% group_modify(~ {
    if (nrow(.x) < 2 * MIN_CELL || length(unique(.x$event)) < 2)
      return(tibble(n = nrow(.x), slope = NA_real_, slope_lo = NA_real_, slope_hi = NA_real_))
    f <- glm(as.formula(paste("event ~ z_vtpfvc +", sev)), data = .x, family = binomial)
    s <- summary(f)$coefficients["z_vtpfvc", ]
    ci <- suppressMessages(confint.default(f))["z_vtpfvc", ]
    tibble(n = nrow(.x), slope = unname(s["Estimate"]),
           slope_lo = unname(ci[1]), slope_hi = unname(ci[2]))
  }) %>% ungroup() %>%
    filter(n >= MIN_CELL) %>%
    mutate(subgroup = slabel, level = as.character(level), interaction_lrt_p = lrt_p)
})
write_csv(stratum_slopes, file.path(final_dir, paste0("dose_response_homogeneity_", site_name, ".csv")))
message("VT/PFVC dose-response homogeneity (interaction LRT p, null = homogeneous):")
walk(unique(stratum_slopes$subgroup), ~ message("  ", .x, ": p = ",
     signif(stratum_slopes$interaction_lrt_p[stratum_slopes$subgroup == .x][1], 3)))

# =============================================================================
# 8f. Recoil-absorption: age x dose interaction across the normalization ladder
# =============================================================================
# Age's elastic-recoil channel is EXPECTED to modify the dose-response. VT/PFVC
# normalizes by PREDICTED lung size and cannot see the patient's actual
# age-changed recoil, so its dose-response should carry an age interaction. DP
# (= VT/Crs) and MP/Crs embed MEASURED compliance, so they already absorb recoil
# and should be more age-invariant. PREDICTION: the per-SD age-slope difference
# (old - young) shrinks along VT/PFVC -> DP -> MP/Crs. A residual interaction in
# MP/Crs would be susceptibility/baseline-risk, NOT recoil.
recoil_dat <- analytic %>% mutate(mp_crs = mechanical_power * ers / 1000)
ladder  <- tribble(~var, ~label, "vtpfvc", "VT/PFVC", "dp", "DP", "mp_crs", "MP/Crs")
sev_age <- "sofa_total + sf_ratio + bmi"

# per-1-SD logit slope of mortality at a given age (includes the interaction)
slope_at_age <- function(fit, a) {
  base <- tibble(sofa_total = mean(recoil_dat$sofa_total, na.rm = TRUE),
                 sf_ratio = mean(recoil_dat$sf_ratio, na.rm = TRUE),
                 bmi = mean(recoil_dat$bmi, na.rm = TRUE), age_at_admission = a)
  predict(fit, base %>% mutate(.z = 1)) - predict(fit, base %>% mutate(.z = 0))
}
recoil_fit <- function(var) {
  d <- recoil_dat %>% filter(!is.na(.data[[var]])) %>% mutate(.z = as.numeric(scale(.data[[var]])))
  list(d = d,
       full = glm(as.formula(paste("event ~ .z * ns(age_at_admission, df = 3) +", sev_age)),
                  data = d, family = binomial),
       red  = glm(as.formula(paste("event ~ .z + ns(age_at_admission, df = 3) +", sev_age)),
                  data = d, family = binomial))
}
age_q <- quantile(recoil_dat$age_at_admission, c(0.05, 0.10, 0.90, 0.95), na.rm = TRUE)
recoil_tbl <- pmap_dfr(ladder, function(var, label) {
  f <- recoil_fit(var)
  tibble(metric = label, n = nrow(f$d),
         interaction_lrt_p = anova(f$red, f$full, test = "LRT")$`Pr(>Chi)`[2],
         slope_young = slope_at_age(f$full, age_q[2]),
         slope_old   = slope_at_age(f$full, age_q[3]),
         slope_diff  = slope_at_age(f$full, age_q[3]) - slope_at_age(f$full, age_q[2]))
}) %>% mutate(metric = factor(metric, levels = ladder$label))
write_csv(recoil_tbl, file.path(final_dir, paste0("dose_recoil_absorption_", site_name, ".csv")))
message("Recoil-absorption (|slope_diff| should shrink VT/PFVC -> MP/Crs):")
walk(seq_len(nrow(recoil_tbl)), ~ message(sprintf(
  "  %-8s LRTp=%.3g  slope young=%.3f old=%.3f  diff=%+.3f",
  as.character(recoil_tbl$metric[.x]), recoil_tbl$interaction_lrt_p[.x],
  recoil_tbl$slope_young[.x], recoil_tbl$slope_old[.x], recoil_tbl$slope_diff[.x])))

# age-slope curve across the ladder (for the figure)
age_grid <- seq(age_q[1], age_q[4], length.out = 40)
recoil_curve <- pmap_dfr(ladder, function(var, label) {
  f <- recoil_fit(var)
  tibble(metric = label, age = age_grid,
         slope = vapply(age_grid, function(a) slope_at_age(f$full, a), numeric(1)))
}) %>% mutate(metric = factor(metric, levels = ladder$label))

# =============================================================================
# 8g. Elastance decomposition: where does the age-modification live?
# =============================================================================
# The composite metrics (DP = Ers x VT; MP/Crs = MP x Ers ~ Ers^2 x VT^2) are
# PRODUCTS of two physically independent components -- volume (clinician-set) and
# elastance (patient physiology). Testing "does age x Ers mediate age x MP/Crs"
# is CIRCULAR (MP/Crs contains Ers, and Crs = 1/Ers in our data). Instead we
# enter the two components SEPARATELY, each x ns(age), and ask which one carries
# the age-modification. Signed predictions:
#   * Volume is the (age-invariant) Dreyfuss mediator -> age-interaction VANISHES
#     under decomposition (both components age-flat); the composite age-steepening
#     was a proxy artifact of carrying elastance.
#   * Elastance/pressure has its own age-modified effect -> age-interaction LOADS
#     on log(Ers).
# MEASUREMENT CAVEAT: Ers = dP_rs/VT is RESPIRATORY-SYSTEM elastance, NOT specific
# elastance (no V0, no transpulmonary pressure); it conflates tissue stiffness,
# lung volume, and chest wall. So a log(Ers) age-interaction localizes the
# modification to the elastance component but does NOT isolate specific elastance.
# Volume component is the SIZE-NORMALIZED strain proxy log(VT/PFVC), NOT raw
# log(VT): raw tidal volume is body-size confounded (bigger patients -> bigger VT
# AND lower mortality), which contaminates the volume slope (negative) and its
# age-interaction. VT/PFVC removes that via predicted-size normalization and
# matches the §8f recoil exposure.
decomp_dat <- analytic %>%
  mutate(log_vol = log(vtpfvc), log_Ers = log(ers)) %>%
  filter(is.finite(log_vol), is.finite(log_Ers), ers > 0, vtpfvc > 0)
sev_demo <- "sofa_total + sf_ratio + bmi + sex_category + race_category"
f_full  <- as.formula(paste("event ~ log_vol*ns(age_at_admission,3) + log_Ers*ns(age_at_admission,3) +", sev_demo))
f_noVT  <- as.formula(paste("event ~ log_vol + log_Ers*ns(age_at_admission,3) +", sev_demo))
f_noErs <- as.formula(paste("event ~ log_vol*ns(age_at_admission,3) + log_Ers +", sev_demo))
decomp_fit <- glm(f_full, data = decomp_dat, family = binomial)

modal <- function(x) names(sort(table(x), decreasing = TRUE))[1]
base_row <- function(d, a) tibble(
  log_vol = mean(d$log_vol), log_Ers = mean(d$log_Ers),
  sofa_total = mean(d$sofa_total), sf_ratio = mean(d$sf_ratio), bmi = mean(d$bmi),
  sex_category = modal(d$sex_category), race_category = modal(d$race_category),
  age_at_admission = a)
slope_comp <- function(fit, d, comp, a) {           # per +1 log-unit logit slope at age a
  b0 <- base_row(d, a); b1 <- b0; b1[[comp]] <- b0[[comp]] + 1
  as.numeric(predict(fit, b1) - predict(fit, b0))
}
ages_yo <- quantile(decomp_dat$age_at_admission, c(0.10, 0.90))
decomp_tbl <- tibble(
  component = c("log(VT/PFVC) (strain)", "log_Ers (elastance)"),
  age_interaction_lrt_p = c(
    anova(glm(f_noVT, data = decomp_dat, family = binomial), decomp_fit, test = "LRT")$`Pr(>Chi)`[2],
    anova(glm(f_noErs, data = decomp_dat, family = binomial), decomp_fit, test = "LRT")$`Pr(>Chi)`[2]),
  slope_young = c(slope_comp(decomp_fit, decomp_dat, "log_vol", ages_yo[1]),
                  slope_comp(decomp_fit, decomp_dat, "log_Ers", ages_yo[1])),
  slope_old   = c(slope_comp(decomp_fit, decomp_dat, "log_vol", ages_yo[2]),
                  slope_comp(decomp_fit, decomp_dat, "log_Ers", ages_yo[2]))) %>%
  mutate(slope_diff = slope_old - slope_young)
# brief in-cohort CONTEXT (not a proof; healthy-lung recoil-age is established):
# net direction of measured Ers with age, severity-adjusted.
ers_age_fit <- lm(log_Ers ~ ns(age_at_admission, 3) + sofa_total + sf_ratio + bmi, data = decomp_dat)
ers_dir <- predict(ers_age_fit, base_row(decomp_dat, ages_yo[2])) -
           predict(ers_age_fit, base_row(decomp_dat, ages_yo[1]))
decomp_tbl <- decomp_tbl %>% mutate(ers_log_change_old_minus_young = unname(ers_dir))
write_csv(decomp_tbl, file.path(final_dir, paste0("dose_elastance_decomposition_", site_name, ".csv")))
message("Elastance decomposition (which component carries the age-modification):")
walk(seq_len(nrow(decomp_tbl)), ~ message(sprintf(
  "  %-20s LRTp=%.3g  slope young=%.3f old=%.3f  diff=%+.3f",
  decomp_tbl$component[.x], decomp_tbl$age_interaction_lrt_p[.x],
  decomp_tbl$slope_young[.x], decomp_tbl$slope_old[.x], decomp_tbl$slope_diff[.x])))
message(sprintf("  (context) measured log(Ers) old-young, severity-adj = %+.3f", unname(ers_dir)))

# per-component age-slope curve (for the figure)
decomp_curve <- bind_rows(
  tibble(component = "log(VT/PFVC) (strain)", age = age_grid,
         slope = vapply(age_grid, function(a) slope_comp(decomp_fit, decomp_dat, "log_vol", a), numeric(1))),
  tibble(component = "log_Ers (elastance)", age = age_grid,
         slope = vapply(age_grid, function(a) slope_comp(decomp_fit, decomp_dat, "log_Ers", a), numeric(1))))

# -----------------------------------------------------------------------------
# ABSOLUTE-SCALE strain effect by age (where the "younger survive biotrauma"
# claim lives). The multiplicative slope (above) can be age-flat while the
# ABSOLUTE risk difference rises with age, because RD = beta_strain(age) x p(1-p)
# scales with baseline risk. We compute the average marginal RD (g-computation):
# per-patient predicted-risk difference at observed strain vs strain + 1 SD,
# averaged in a moving age window. This integrates BOTH the multiplicative
# interaction AND the baseline-risk channel, so it directly tests the claim on
# the absolute scale rather than inferring it from the OR.
# Use a STRAIN-MARGINAL model (VT/PFVC + severity, ns(age) interaction; the §8f
# parameterization), NOT the decomposition: the biotrauma claim is about the
# TOTAL effect of strain, so we must NOT condition on elastance (which is
# entangled with strain via DP and flips the slope's sign by over-adjustment).
WIN       <- 7.5                                     # +/- yr smoothing window
N_BOOT_RD <- if (is_synthetic) 100L else 300L
sd_vol    <- sd(decomp_dat$log_vol)
f_rd   <- as.formula("event ~ log_vol*ns(age_at_admission,3) + sofa_total + sf_ratio + bmi")
rd_fit <- glm(f_rd, data = decomp_dat, family = binomial)
rd_dat <- decomp_dat %>% mutate(.p0 = predict(rd_fit, ., type = "response"))
rd_dat$.p1 <- predict(rd_fit, decomp_dat %>% mutate(log_vol = log_vol + sd_vol),
                      type = "response")
rd_dat <- rd_dat %>% mutate(.rd = .p1 - .p0)         # absolute RD per +1 SD strain
win_mean <- function(d, col, a) mean(d[[col]][abs(d$age_at_admission - a) <= WIN])
rd_curve <- tibble(age = age_grid) %>% rowwise() %>%
  mutate(rd_pp = 100 * win_mean(rd_dat, ".rd", age),
         baseline_mort_pct = 100 * win_mean(rd_dat, "event", age)) %>% ungroup()
boot_rd <- matrix(NA_real_, N_BOOT_RD, length(age_grid))
for (b in seq_len(N_BOOT_RD)) {
  bd <- decomp_dat[sample.int(nrow(decomp_dat), replace = TRUE), , drop = FALSE]
  fb <- tryCatch(glm(f_rd, data = bd, family = binomial), error = function(e) NULL)
  if (is.null(fb)) next
  bd$.rdb <- predict(fb, bd %>% mutate(log_vol = log_vol + sd_vol), type = "response") -
             predict(fb, bd, type = "response")
  boot_rd[b, ] <- vapply(age_grid, function(a) 100 * win_mean(bd, ".rdb", a), numeric(1))
}
rd_curve <- rd_curve %>%
  mutate(rd_lo = apply(boot_rd, 2, quantile, 0.025, na.rm = TRUE),
         rd_hi = apply(boot_rd, 2, quantile, 0.975, na.rm = TRUE))
write_csv(rd_curve, file.path(final_dir, paste0("dose_strain_rd_by_age_", site_name, ".csv")))
message(sprintf("Absolute strain RD per +1 SD (pp): young(%.0f)=%.2f  old(%.0f)=%.2f",
                age_grid[1], rd_curve$rd_pp[1],
                age_grid[length(age_grid)], rd_curve$rd_pp[length(age_grid)]))

# =============================================================================
# 8h. Figures
# =============================================================================
pdf_path <- function(stub) file.path(final_dir, paste0(stub, "_", site_name, ".pdf"))
sub_levels <- c("Male","Female","WHITE","BLACK","OTHER","Young","Middle","Old","Tall","Short")
order_lvl <- function(x) factor(x, levels = intersect(sub_levels, unique(x)))

# (8f-1) within-band VT/PFVC spread: violin per subgroup level at fixed VT/PBW 6-8
long_sg <- map_dfr(seq_len(nrow(subgroup_vars)), function(i) {
  analytic %>% transmute(subgroup = subgroup_vars$slabel[i],
                         level = as.character(.data[[subgroup_vars$svar[i]]]), vtpfvc)
})
if (nrow(long_sg) > 0) {
  ls <- long_sg %>% mutate(subgroup = factor(subgroup, levels = subgroup_vars$slabel),
                           level = order_lvl(level))
  p_spread <- ggplot(ls, aes(level, vtpfvc, fill = subgroup)) +
    geom_violin(scale = "width", alpha = 0.7, colour = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.6) +
    geom_hline(yintercept = VTPFVC_TARGET, linetype = 2, colour = "grey40") +
    facet_grid(~ subgroup, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = okabe[c(1, 2, 3, 4)], guide = "none") +
    labs(x = NULL, y = "Delivered dose VT/PFVC (% predicted FVC)",
         title = "Same treatment, different dose: VT/PFVC spread at fixed VT/PBW 6-8",
         subtitle = paste0(site_name,
           if (is_synthetic) " (SYNTHETIC)" else "",
           " - dashed = 11% harm cut; a single protocol delivers a wide, demographically-shifted dose")) +
    theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(pdf_path("dose_within_band_spread"), p_spread, width = 11, height = 4.5)
}

# (8f-2) dose-response homogeneity: per-stratum VT/PFVC slope (forest)
if (nrow(stratum_slopes) > 0 && any(!is.na(stratum_slopes$slope))) {
  ss <- stratum_slopes %>% filter(!is.na(slope)) %>%
    mutate(subgroup = factor(subgroup, levels = subgroup_vars$slabel), level = order_lvl(level))
  lab <- ss %>% distinct(subgroup, interaction_lrt_p) %>%
    mutate(lab = paste0("interaction p=", signif(interaction_lrt_p, 2)))
  p_homog <- ggplot(ss, aes(slope, level, colour = subgroup)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = slope_lo, xmax = slope_hi)) +
    geom_text(data = lab, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE,
              hjust = 1.05, vjust = 1.4, size = 3, colour = "grey30") +
    facet_grid(subgroup ~ ., scales = "free_y", space = "free_y") +
    scale_colour_manual(values = okabe[c(1, 2, 3, 4)], guide = "none") +
    labs(x = "VT/PFVC -> 60-day mortality log-OR per SD (severity-adjusted)", y = NULL,
         title = "VT/PFVC dose-response is homogeneous across strata (structural dose)",
         subtitle = paste0(site_name,
           " - overlapping per-stratum slopes + null interaction = invariant effect (run both sites)")) +
    theme_minimal(base_size = 10)
  ggsave(pdf_path("dose_response_homogeneity"), p_homog, width = 8.5, height = 9)
}

# (8g-3) recoil-absorption: per-SD age-slope curve across the normalization ladder
if (nrow(recoil_curve) > 0) {
  p_recoil <- ggplot(recoil_curve, aes(age, slope, colour = metric)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = c("VT/PFVC" = okabe[5], "DP" = okabe[1], "MP/Crs" = okabe[3]),
                        name = "Dose metric (predicted -> measured size)") +
    labs(x = "Age (years)", y = "Per-SD mortality log-OR slope",
         title = "Recoil-absorption: age-modification of the dose-response across the ladder",
         subtitle = paste0(site_name, if (is_synthetic) " (SYNTHETIC)" else "",
           " - flatter = more age-invariant; measured-compliance metrics (DP, MP/Crs) should absorb the recoil channel")) +
    theme_minimal(base_size = 10)
  ggsave(pdf_path("dose_recoil_absorption"), p_recoil, width = 8, height = 4.5)
}

# (8h-4) elastance decomposition: which component's slope is age-modified
if (nrow(decomp_curve) > 0) {
  p_decomp <- ggplot(decomp_curve, aes(age, slope, colour = component)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = c("log(VT/PFVC) (strain)" = okabe[3], "log_Ers (elastance)" = okabe[6]),
                        name = NULL) +
    labs(x = "Age (years)", y = "Per +1 log-unit mortality logit slope",
         title = "Elastance decomposition: age-modification by component (volume vs elastance)",
         subtitle = paste0(site_name, if (is_synthetic) " (SYNTHETIC)" else "",
           " - flat volume slope = Dreyfuss mediator; an age-sloped elastance = pressure has its own age effect")) +
    theme_minimal(base_size = 10)
  ggsave(pdf_path("dose_elastance_decomposition"), p_decomp, width = 8, height = 4.5)
}

# (8h-5) absolute-scale strain effect by age (the biotrauma-survival claim)
if (nrow(rd_curve) > 0) {
  p_rd <- ggplot(rd_curve, aes(age, rd_pp)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    geom_ribbon(aes(ymin = rd_lo, ymax = rd_hi), alpha = 0.2, fill = okabe[3]) +
    geom_line(colour = okabe[3], linewidth = 1) +
    labs(x = "Age (years)", y = "Absolute mortality RD per +1 SD strain (pct points)",
         title = "Absolute-scale strain effect by age (biotrauma survival)",
         subtitle = paste0(site_name, if (is_synthetic) " (SYNTHETIC)" else "",
           " - RD = beta_strain(age) x p(1-p); rises with age via baseline risk even if the OR is age-flat")) +
    theme_minimal(base_size = 10)
  ggsave(pdf_path("dose_strain_rd_by_age"), p_rd, width = 8, height = 4.5)
}

message("Wrote 6 tables + 5 figures to ", final_dir)
message("Script 08 complete.")
