# =============================================================================
# Script 28 (height fingerprint): does the marker follow the ratio's
# sex-reversed height curve?
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# The claim tested (manuscript argument, Claim 5c.2). At a fixed VT/PBW, strain is
# proportional to PBW/PFVC. Devine PBW is a line in height with an intercept; GLI
# FVC is a power law in height. Their ratio is therefore an inverted U in height,
# and the U differs by sex: it peaks near 166 cm in men and near 183 cm in women
# (range 0.05 and 0.07 log units from 150 to 200 cm). If strain injures, the
# marker's course follows that curve, sex by sex. A direct effect of height, through
# body size or socioeconomic position, has no reason to reverse between the sexes or
# to peak where Devine's intercept stops dominating.
#
# The fingerprint F is the ratio's height piece at each patient's own sex
# (height_fingerprint() in 20_biotrauma_grid.R). GLI's height coefficient does not
# depend on age or race, so F is exactly the part of log PBW/PFVC that height and
# sex determine. The model gives the height function every sex shares, and sex
# itself, the same terms F gets, level and rate; what is left of F is the
# sex-specific departure from the shared curve, and its coefficient carries the
# contrast:
#
#   log marker ~ ns(day, 3) + F + F:day
#                + ns(height, k) + ns(height, k):day + sex + sex:day   (shared)
#                + VT/PBW at the index + (previous day's VT/PBW - index VT/PBW)
#                + baseline marker (the index-day value) + lagged SF and pressor
#                + non-respiratory SOFA
#                [+ ns(age, 4) + race]                                 (adjusted)
#   random intercept and slope on day per patient (unstructured)
#
# Leaving the shared height function out of the rate would let any common effect
# of height on the marker's course load onto F:day, so it is in both.
#
# What is reported, per marker, shared-smooth df k (3 primary; 4, 5) and adjustment:
#   F:day            divergence per day per log unit of PBW/PFVC moved by height
#                    within sex. The primary read.
#   ldisc:day        the predicted value F:day is read against: the same per log unit
#                    of log PBW/PFVC, from figure 4's model with log PBW/PFVC in
#                    place of log PFVC and sex held fixed (no shared height smooth).
#                    Adjusted, age and race are held fixed too, so the ratio moves
#                    only through height within sex: the fingerprint's own variation
#                    without the shared smooth taken out. Unadjusted, the ratio moves
#                    through height, age and race within sex. Sex is in both pairs
#                    because the fingerprint model always holds it. Under strain the
#                    two agree in sign and size: the ratio should act the same
#                    whichever input moved it. These rows carry model "whole ratio"
#                    (the name the pooling reads).
#   identifying SD   the SD of F left after ns(height, k) and sex (and, adjusted,
#                    age and race). The formula's whole range is 0.05 to 0.07; the
#                    lever left after the shared curve is a small fraction of the
#                    benchmark's, and it is reported at k = 3, 4 and 5 so a
#                    reader can see that the limit is the formula, not the smooth.
#                    For the ldisc rows the column holds the SD of log PBW/PFVC
#                    left after the same terms the benchmark holds fixed (sex; and,
#                    adjusted, age and race). The table therefore carries the minimum detectable
#                    effect (80% power, two-sided 5%) beside the predicted value, so
#                    a null can be read as uninformative rather than as a refutation.
#                    The comparison is the ventilated arm's read; in a control the
#                    predicted value is the control's own ratio rate.
#   ladder           maximum-likelihood fits: shared only; shared + F (1 df at
#                    level, 1 at rate, and 2 more for F x anchor in a
#                    severity-standardised control); sex-specific free height curves
#                    (ns(height, k) x sex, level and rate). Shared vs free asks whether there is
#                    any sex-specific height shape; free vs fingerprint (AIC, not
#                    nested) asks whether the formula's shape is enough to describe
#                    it. A direct height effect can make a sex-specific shape, but
#                    not this one.
#   curves           the rate by height and sex from the free model and from the
#                    fingerprint model, beside F itself: the plot a reader checks.
#   sev_center       every table carries the severity centre the fit was read at (NA
#                    when none), so the ventilated run can refuse a control table that
#                    was not read at the ventilated severity.
#
# "Unadjusted" drops age and race only. Height and sex stay in every model: they
# define the contrast.
#
# The negative control. Run with PBWPFVC_COHORT=nosupport (no ventilator, so no
# dose terms), the same model reads the fingerprint where no PBW-scaled volume is
# applied; strain predicts null there. With PBWPFVC_JM_SEV_CENTER set, the control's
# rate is read at the ventilated cohort's severity, as in figure 4 (the anchor, its
# rate, and its modification of F level and rate are added). Run the control first:
# the ventilated run then writes the difference in differences if the control's
# table is on disk, and stops if that table was fitted without a severity centre.
#
# Death and extubation before day 7 are not modelled: this is a linear mixed model
# alone (nlme), so the estimate is conditional on remaining under observation. Days
# count from the index (the first qualifying ventilator row; ICU admission in the
# control), and the panel holds no measurement after a patient's death, extubation or
# escalation (21_biotrauma_panel.R).
#
# Inputs: the 7-day daily panel of 21_biotrauma_panel.R (jm_long_7d, jm_surv_7d).
# Outputs, in final/injury/ (a control cohort's in final/controls/):
#   fingerprint_{marker}_{site}.csv          coefficients, identifying SD, MDE
#   fingerprint_ladder_{marker}_{site}.csv   likelihood-ratio tests and AIC
#   fingerprint_curves_{marker}_{site}.csv   rate by height and sex, and F
#   fingerprint_did_{marker}_{site}.csv      ventilated minus control F:day
#   fingerprint_{marker}_{site}.pdf
#
# Usage (from the repo root; the panel first if it is missing):
#   export PBWPFVC_JM_GRID=daily PBWPFVC_JM_HORIZON=7
#   uvr run code/21_biotrauma_panel.R
#   PBWPFVC_COHORT=nosupport uvr run code/21_biotrauma_panel.R
#   PBWPFVC_COHORT=nosupport PBWPFVC_JM_SEV_CENTER=platelets=<ventilated mean anchor> \
#     uvr run code/28_height_fingerprint.R
#   uvr run code/28_height_fingerprint.R
# PBWPFVC_INJ_MARKER picks the marker (platelets by default; creatinine, bilirubin).
# The ventilated mean anchor is in final/injury/jm_severity_anchor_mean_7d_{site}.csv.
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(arrow); library(here); library(splines); library(nlme); library(patchwork)
})
rm(list = ls())
options(width = 200)
source("utils/config.R")
site_name  <- config$site_name
output_dir <- config$output_dir
final_dir  <- final_dir_for("injury")
source(here("code", "20_biotrauma_grid.R"))
if (JM_GRID != "daily" || JM_HORIZON != 7L)
  stop("the fingerprint reads figure 4's panel: set PBWPFVC_JM_GRID=daily PBWPFVC_JM_HORIZON=7")

MARKER <- Sys.getenv("PBWPFVC_INJ_MARKER", "platelets")
stopifnot(MARKER %in% c("platelets", "creatinine", "bilirubin"))
y_col  <- MARKER
y0_col <- c(platelets = "platelet_0", creatinine = "creatinine_0", bilirubin = "bilirubin_0")[[MARKER]]
HAS_DOSE <- config$cohort == "imv"
SHARED_DF <- c(3L, 4L, 5L)            # df of the shared height smooth; 4 and 5 are sensitivity analyses
PRIMARY_DF <- 3L                      # the primary df of the shared height smooth
# minimum counts for a fit: the CLIF minimum-count standard
MIN_PATIENTS <- 100L
MIN_PATIENTS_PER_SEX <- 50L
# the minimum detectable effect is MDE_Z standard errors: 1.96 + 0.84, a two-sided
# 5% test with 80% power
MDE_Z <- qnorm(0.975) + qnorm(0.80)
DEMO <- c("ns(age10, 4)", "race_category")
okabe <- c(Male = "#0072B2", Female = "#D55E00")
message("=== 28_height_fingerprint: ", MARKER, ", site ", site_name, " (cohort ", config$cohort, ") ===")

# =============================================================================
# 1. The panel, with the fingerprint and the whole ratio per patient
# =============================================================================
panel_path <- file.path(output_dir, paste0("jm_long_", h_suffix, ".parquet"))
if (!file.exists(panel_path))
  stop("no ", h_suffix, " panel for this cohort: run  PBWPFVC_JM_GRID=daily PBWPFVC_JM_HORIZON=7 uvr run code/21_biotrauma_panel.R")
long <- read_parquet(panel_path)
surv <- read_parquet(file.path(output_dir, paste0("jm_surv_", h_suffix, ".parquet")))

# severity standardisation of the control (20_biotrauma_grid.R): the marker's own
# anchor, centred at the ventilated cohort's mean
sev_center <- sev_center_for(MARKER)
# each exclusion in turn, logged with the number of patients it drops
drop_step <- function(dat, reason, ...) {
  kept <- dat %>% filter(...)
  message(sprintf("  %-58s dropped %5d, %d left", reason, n_distinct(dat$hospitalization_id) - n_distinct(kept$hospitalization_id),
                  n_distinct(kept$hospitalization_id)))
  kept
}
message(sprintf("  %d patients in the panel", nrow(surv)))
pt <- surv %>%
  drop_step("no baseline marker", !is.na(.data[[y0_col]])) %>%
  drop_step("no non-respiratory SOFA", !is.na(np_sofa)) %>%
  drop_step("no height", !is.na(height_cm)) %>%
  drop_step("no age", !is.na(age10)) %>%
  drop_step("no race", !is.na(race_category)) %>%
  drop_step("no index VT/PBW (ventilated cohort only)", if (HAS_DOSE) !is.na(vtpbw_idx) else TRUE) %>%
  mutate(fingerprint = height_fingerprint(height_cm, sex_category),
         ldisc = log(pbw / pfvc_gli),
         sex_female = as.numeric(sex_category == "Female"),
         log_y0 = log(.data[[y0_col]]),
         sev_anchor_c = if (is.na(sev_center)) NA_real_ else
           rowSums(as.matrix(pick(all_of(anchor_components(MARKER))))) - sev_center)
if (!is.na(sev_center)) pt <- pt %>% drop_step("no severity anchor", !is.na(sev_anchor_c))

d_joined <- long %>%
  filter(period >= 1L, !is.na(.data[[y_col]]), .data[[y_col]] > 0, !is.na(l_sf), !is.na(l_pressor),
         if (HAS_DOSE) !is.na(l_vtpbw_within) else TRUE) %>%
  select(hospitalization_id, vent_day, y = all_of(y_col), l_sf, l_pressor, any_of("l_vtpbw_within")) %>%
  inner_join(pt %>% select(hospitalization_id, fingerprint, ldisc, height_cm, sex_category, sex_female,
                           age10, race_category, np_sofa, log_y0, sev_anchor_c, any_of("vtpbw_idx")),
             by = "hospitalization_id") %>%
  mutate(log_y = log(y), l_log_sf = log(l_sf))
message(sprintf("  %-58s dropped %5d, %d left", "no usable post-baseline day (marker, lagged SF, pressor, dose)",
                n_distinct(pt$hospitalization_id) - n_distinct(d_joined$hospitalization_id),
                n_distinct(d_joined$hospitalization_id)))
d_all <- d_joined %>%
  group_by(hospitalization_id) %>% drop_step("fewer than two post-baseline days", n() >= 2L) %>% ungroup() %>%
  mutate(id = factor(hospitalization_id))
pt_used <- pt %>% filter(hospitalization_id %in% d_all$hospitalization_id)
n_patients <- nrow(pt_used)
n_by_sex <- count(pt_used, sex_category)
message(sprintf("  %d patients with a baseline and two or more post-baseline days, %d patient-days; %s",
                n_patients, nrow(d_all), paste(n_by_sex$sex_category, n_by_sex$n, collapse = ", ")))
if (n_patients < MIN_PATIENTS || any(n_by_sex$n < MIN_PATIENTS_PER_SEX))
  stop("too few patients for the fingerprint (need ", MIN_PATIENTS, ", and ", MIN_PATIENTS_PER_SEX, " of each sex)")

# =============================================================================
# 2. The model
# =============================================================================
# The shared height smooth is built once, on the patients, as numeric basis
# columns, so the curves in section 4 evaluate it at the same knots
height_basis <- map(set_names(SHARED_DF), ~ ns(pt_used$height_cm, df = .x))
add_basis <- function(dat, k) {
  B <- predict(height_basis[[as.character(k)]], dat$height_cm)
  colnames(B) <- paste0("hs", k, "_", seq_len(k))
  bind_cols(dat, as_tibble(B))
}
hs_names <- function(k) paste0("hs", k, "_", seq_len(k))
ctrl <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)

# size_terms: exposures given a level and a rate; height: "shared", "free" or "none"
rhs_of <- function(size_terms, k, height = "shared", adjusted = TRUE) {
  hs <- hs_names(k)
  height_terms <- switch(height,
    # the benchmark: log PBW/PFVC as figure 4 enters log PFVC, with sex held fixed in
    # both adjustments, as the fingerprint model holds it
    none   = "sex_female",
    shared = c(hs, paste0(hs, ":vent_day"), "sex_female", "sex_female:vent_day"),
    free   = c(hs, paste0(hs, ":vent_day"), "sex_female", "sex_female:vent_day",
               paste0(hs, ":sex_female"), paste0(hs, ":sex_female:vent_day")))
  sev_terms <- if (!is.na(sev_center)) c("sev_anchor_c", "sev_anchor_c:vent_day",
                                         if (length(size_terms)) c(paste0(size_terms, ":sev_anchor_c"),
                                                                   paste0(size_terms, ":vent_day:sev_anchor_c")))
  paste(c("ns(vent_day, 3)", size_terms, if (length(size_terms)) paste0(size_terms, ":vent_day"),
          height_terms, sev_terms, if (HAS_DOSE) c("l_vtpbw_within", "vtpbw_idx"),
          "log_y0", "l_log_sf", "l_pressor", "np_sofa", if (adjusted) DEMO), collapse = " + ")
}
fit_lme <- function(rhs, dat, method = "REML")
  lme(as.formula(paste("log_y ~", rhs)), random = ~ vent_day | id, data = dat, control = ctrl, method = method)
# a fixed-effect row by its components, whatever order R wrote them in
term_by <- function(b, comps) names(b)[sapply(strsplit(names(b), ":"), setequal, comps)]
coef_row <- function(fit, comps) {
  b <- fixef(fit); nm <- term_by(b, comps); se <- sqrt(diag(as.matrix(vcov(fit))))[nm]
  tibble(term = paste(comps, collapse = ":"), estimate = unname(b[nm]), se = unname(se))
}

# the identifying variation: what is left of F after the shared height smooth and
# sex (and, adjusted, age and race), on the patients
identifying_sd <- function(k, adjusted) {
  dd <- add_basis(pt_used, k)
  f <- as.formula(paste("fingerprint ~", paste(c(hs_names(k), "sex_female", if (adjusted) DEMO), collapse = " + ")))
  sd(resid(lm(f, data = dd)))
}
# the benchmark's variation: what is left of log PBW/PFVC after the terms its model
# holds fixed (sex; and, adjusted, age and race)
benchmark_sd <- function(adjusted) {
  f <- as.formula(paste("ldisc ~", paste(c("sex_female", if (adjusted) DEMO), collapse = " + ")))
  sd(resid(lm(f, data = pt_used)))
}

results <- list(); ladder <- list(); fits_keep <- list()
for (k in SHARED_DF) for (adjusted in c(TRUE, FALSE)) {
  adj_lab <- if (adjusted) "adjusted" else "unadjusted"
  dk <- add_basis(d_all, k)
  fp <- fit_lme(rhs_of("fingerprint", k, "shared", adjusted), dk)
  id_sd <- identifying_sd(k, adjusted)
  rows <- bind_rows(coef_row(fp, c("fingerprint", "vent_day")) %>% mutate(model = "fingerprint", quantity = "rate"),
                    coef_row(fp, "fingerprint") %>% mutate(model = "fingerprint", quantity = "level"))
  if (k == PRIMARY_DF) {
    # the predicted value: figure 4's model with log PBW/PFVC in place of log PFVC and
    # sex held fixed; height is on the ratio's causal path, so no height smooth here
    rt <- fit_lme(rhs_of("ldisc", k, "none", adjusted), dk)
    rows <- bind_rows(rows,
                      coef_row(rt, c("ldisc", "vent_day")) %>% mutate(model = "whole ratio", quantity = "rate"),
                      coef_row(rt, "ldisc") %>% mutate(model = "whole ratio", quantity = "level"))
    fits_keep[[adj_lab]] <- fp
  }
  results[[length(results) + 1]] <- rows %>%
    mutate(marker = MARKER, shared_df = k, adjustment = adj_lab,
           lo = estimate - 1.96 * se, hi = estimate + 1.96 * se, p = 2 * pnorm(-abs(estimate / se)),
           # the fingerprint's residual SD; for the whole-ratio rows, the SD of log PBW/PFVC
           # left after the terms the benchmark holds fixed
           identifying_sd = if_else(model == "fingerprint", id_sd, benchmark_sd(adjusted)),
           mde_80 = MDE_Z * se,
           n_patients = n_patients, n_rows = nrow(dk))
  # the ladder, by maximum likelihood so the likelihoods compare
  ml <- list(shared      = fit_lme(rhs_of(NULL, k, "shared", adjusted), dk, "ML"),
             fingerprint = fit_lme(rhs_of("fingerprint", k, "shared", adjusted), dk, "ML"),
             free        = fit_lme(rhs_of(NULL, k, "free", adjusted), dk, "ML"))
  if (k == PRIMARY_DF) fits_keep[[paste0("free_", adj_lab)]] <- fit_lme(rhs_of(NULL, k, "free", adjusted), dk)
  lr <- function(small, big, test) {
    ll_s <- logLik(ml[[small]]); ll_b <- logLik(ml[[big]])
    stat <- as.numeric(2 * (ll_b - ll_s)); df <- attr(ll_b, "df") - attr(ll_s, "df")
    tibble(test = test, small = small, big = big, lr = stat, df = df, p = pchisq(stat, df, lower.tail = FALSE),
           d_aic = AIC(ml[[big]]) - AIC(ml[[small]]))
  }
  ladder[[length(ladder) + 1]] <- bind_rows(
    lr("shared", "fingerprint", "the fingerprint adds to the shared height curve (level and rate; + F x anchor in a standardised control)"),
    lr("shared", "free", "any sex-specific height shape (level and rate)"),
    tibble(test = "free vs fingerprint: not nested, AIC only (< 0 favours the free curves)",
           small = "fingerprint", big = "free", lr = NA_real_, df = NA_real_, p = NA_real_,
           d_aic = AIC(ml$free) - AIC(ml$fingerprint))) %>%
    mutate(marker = MARKER, shared_df = k, adjustment = adj_lab, n_patients = n_patients)
}
res <- bind_rows(results) %>%
  # the predicted value beside every fingerprint row: the benchmark's rate at the
  # same adjustment, and whether it clears the fingerprint's minimum detectable effect
  left_join(bind_rows(results) %>% filter(model == "whole ratio", quantity == "rate") %>%
              select(adjustment, predicted_rate = estimate), by = "adjustment") %>%
  mutate(predicted_rate = if_else(model == "fingerprint" & quantity == "rate", predicted_rate, NA_real_),
         detectable = if_else(model == "fingerprint" & quantity == "rate", abs(predicted_rate) >= mde_80, NA),
         cohort = config$cohort, sev_center = sev_center, site = site_name) %>%
  select(marker, cohort, model, quantity, term, shared_df, adjustment, estimate, se, lo, hi, p,
         predicted_rate, mde_80, detectable, identifying_sd, n_patients, n_rows, sev_center, site)
ladder <- bind_rows(ladder) %>% mutate(cohort = config$cohort, sev_center = sev_center, site = site_name)

# =============================================================================
# 3. Difference in differences (the ventilated run, when the control's table exists)
# =============================================================================
did <- NULL
if (HAS_DOSE) {
  ctrl_file <- file.path(config$final_root, "controls",
                         paste0("fingerprint_", MARKER, "_", config$base_site, "_nosupport.csv"))
  if (file.exists(ctrl_file)) {
    ctl <- read_csv(ctrl_file, show_col_types = FALSE)
    # the control must be read at the ventilated severity, as in figure 4
    if (!"sev_center" %in% names(ctl) || all(is.na(ctl$sev_center)))
      stop("the control's fingerprint table ", ctrl_file, " was fitted without a severity centre: rerun it with ",
           "PBWPFVC_COHORT=nosupport PBWPFVC_JM_SEV_CENTER=", MARKER, "=<ventilated mean anchor>")
    did <- res %>% filter(model == "fingerprint", quantity == "rate") %>%
      select(shared_df, adjustment, estimate_ventilated = estimate, se_ventilated = se, n_patients_ventilated = n_patients) %>%
      inner_join(ctl %>% filter(model == "fingerprint", quantity == "rate") %>%
                   select(shared_df, adjustment, estimate_control = estimate, se_control = se, n_patients_control = n_patients,
                          control_sev_center = sev_center),
                 by = c("shared_df", "adjustment")) %>%
      mutate(did_estimate = estimate_ventilated - estimate_control,
             did_se = sqrt(se_ventilated^2 + se_control^2),
             did_lo = did_estimate - 1.96 * did_se, did_hi = did_estimate + 1.96 * did_se,
             did_p = 2 * pnorm(-abs(did_estimate / did_se)),
             marker = MARKER, site = site_name, .before = 1)
  } else message("  no control table at ", ctrl_file, ": run the no-support cohort first for the difference in differences")
}

# =============================================================================
# 4. Curves: the rate by height and sex, from the free model and the fingerprint
# =============================================================================
# The rate of change attributable to height and sex (log marker per day) on a grid
# spanning the central 95% of each sex's heights, both curves centred at the
# men's median height so they share an origin. Adjusted, primary df.
male_median_height <- median(pt_used$height_cm[pt_used$sex_category == "Male"])
curve_grid <- map_dfr(c("Male", "Female"), function(sx) {
  hh <- pt_used$height_cm[pt_used$sex_category == sx]
  tibble(height_cm = seq(quantile(hh, 0.025), quantile(hh, 0.975), length.out = 60), sex_category = sx)
}) %>% bind_rows(tibble(height_cm = male_median_height, sex_category = "Male", reference = TRUE)) %>%
  mutate(reference = coalesce(reference, FALSE))
rate_curve <- function(fit, model_label, k = PRIMARY_DF) {
  b <- fixef(fit); V <- as.matrix(vcov(fit))
  B <- predict(height_basis[[as.character(k)]], curve_grid$height_cm)
  female <- as.numeric(curve_grid$sex_category == "Female")
  # one row of weights on the fixed effects per grid point: the regressor of every
  # rate term that involves height or sex
  W <- matrix(0, nrow(curve_grid), length(b), dimnames = list(NULL, names(b)))
  for (j in seq_len(k)) {
    W[, term_by(b, c(hs_names(k)[j], "vent_day"))] <- B[, j]
    free_term <- term_by(b, c(hs_names(k)[j], "sex_female", "vent_day"))
    if (length(free_term)) W[, free_term] <- B[, j] * female
  }
  W[, term_by(b, c("sex_female", "vent_day"))] <- female
  fp_term <- term_by(b, c("fingerprint", "vent_day"))
  if (length(fp_term)) W[, fp_term] <- height_fingerprint(curve_grid$height_cm, curve_grid$sex_category)
  # relative to men at their median height, with the uncertainty of the difference
  D <- sweep(W, 2, W[curve_grid$reference, ])
  curve_grid %>%
    mutate(rate = as.numeric(D %*% b), se = sqrt(rowSums((D %*% V) * D)),
           lo = rate - 1.96 * se, hi = rate + 1.96 * se, model = model_label) %>%
    filter(!reference) %>% select(-reference)
}
curves <- bind_rows(rate_curve(fits_keep[["free_adjusted"]], "free sex-specific height curves"),
                    rate_curve(fits_keep[["adjusted"]], "shared curve + fingerprint")) %>%
  mutate(fingerprint = height_fingerprint(height_cm, sex_category), marker = MARKER,
         cohort = config$cohort, sev_center = sev_center, site = site_name)

# =============================================================================
# 5. Print, write, draw
# =============================================================================
message("\nThe fingerprint (F:day) against the benchmark (ldisc:day, sex held fixed), log ", MARKER,
        " per day per log unit of PBW/PFVC.")
message("Strain predicts the two agree; detectable = |predicted| >= the 80%-power minimum detectable effect.")
print(as.data.frame(res %>% filter(quantity == "rate") %>%
                      select(model, shared_df, adjustment, estimate, lo, hi, p, predicted_rate, mde_80, detectable, identifying_sd) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
message("\nLadder (ML):")
print(as.data.frame(ladder %>% select(shared_df, adjustment, test, lr, df, p, d_aic) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
if (!is.null(did)) {
  message("\nDifference in differences, ventilated minus no-support F:day:")
  print(as.data.frame(did %>% select(shared_df, adjustment, did_estimate, did_lo, did_hi, did_p) %>%
                        mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
}
out <- function(x, stem) write_csv(x, file.path(final_dir, paste0(stem, "_", MARKER, "_", site_name, ".csv")))
out(res, "fingerprint"); out(ladder, "fingerprint_ladder"); out(curves, "fingerprint_curves")
if (!is.null(did)) out(did, "fingerprint_did")

formula_panel <- ggplot(curves %>% distinct(sex_category, height_cm, fingerprint),
                        aes(height_cm, fingerprint, colour = sex_category)) +
  geom_line(linewidth = 0.9) + scale_colour_manual(values = okabe, name = NULL) +
  labs(title = "A. The ratio's height piece (formula)", x = "Height, cm",
       y = "log PBW/PFVC, relative to a 170 cm man") + theme_minimal(base_size = 10)
rate_panel <- ggplot(curves, aes(height_cm, rate, colour = sex_category, fill = sex_category)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) + facet_wrap(~ model) +
  scale_colour_manual(values = okabe, name = NULL) + scale_fill_manual(values = okabe, guide = "none") +
  labs(title = paste0("B. Rate of change in log ", MARKER, " by height and sex"),
       subtitle = paste0("adjusted, shared smooth ns(height, ", PRIMARY_DF, "); relative to men at their median height.\n",
                         "Strain predicts B mirrors A, inverted for a marker that falls with injury"),
       x = "Height, cm", y = "log marker per day") + theme_minimal(base_size = 10)
ggsave(file.path(final_dir, paste0("fingerprint_", MARKER, "_", site_name, ".pdf")),
       formula_panel + rate_panel + plot_layout(widths = c(1, 2), guides = "collect"), width = 13, height = 4.8)
message("28_height_fingerprint complete -> ", final_dir)
