# =============================================================================
# Supplement (cross-sectional): does linear age leave curvature for PFVC to carry?
# =============================================================================
# Scripts 04 and 05 adjust for age as a straight line (age10). GLI predicts FVC from
# age non-linearly, so whatever curvature a straight line misses is still in log
# PFVC, and a model with both would credit it to predicted lung size. After
# demographic adjustment PBW/PFVC keeps 6.1% of its variance with linear age and
# 2.4% with a spline (negative_control_identifying_variation_*), so the concern is
# not hypothetical. The figure-4 longitudinal models already use ns(age10, 4); the
# cross-sectional ones do not.
#
# This refits the models whose exposure contains predicted lung size, once with
# 04/05's own linear age and once with ns(age10, 4), changing NOTHING else -- the
# same data, outcomes, exposures, covariates and BMI rule -- and compares:
#
#   coefficients   every exposure term under both age forms; the change is given in
#                  the spline fit's standard errors. Both fits use the same patients,
#                  so this is a measure of how far the adjustment moves the estimate,
#                  not a test of whether it did.
#   fit            AIC under each form; the spline's gain over the line says whether
#                  the curvature is there at all.
#   ranking        within each outcome, whether the exposure specifications keep their
#                  AIC order -- in particular whether VT/PBW + log PFVC still beats
#                  VT/PFVC, and MP/PFVC still beats MP/PBW.
#
# Covered:
#   04  in-hospital death, and elastance, compliance and static DP, for every
#       exposure specification 04 fits (4b-4d)
#   05  in-hospital death for the normalization head-to-head (the mechanic alone,
#       PBW-locked, PFVC-locked, and mechanic plus size, for Ers and for MP)
# Adjusted models only: the unadjusted ones have no age term to change.
#
# Writes to final/supplement/:
#   age_form_coefficients_{site}.csv  every exposure term, both age forms
#   age_form_fit_{site}.csv           AIC per model under both forms, and the ranks
#   age_form_contrasts_{site}.csv     the AIC gaps the manuscript's claims rest on
#                                     (VT/PBW + log PFVC vs VT/PFVC, MP/PFVC vs MP/PBW,
#                                     ...), under each form, and how much survives
#
# Usage: Rscript code/supplement/xsec_age_form_check.R
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(splines) })
rm(list = ls())
source("utils/config.R")
site_name  <- config$site_name
output_dir <- config$output_dir
final_dir  <- final_dir_for("supplement")

AGE_FORMS <- c(linear = "age10", spline = "ns(age10, 4)")

# --- data, as 04 prepares it (04_analysis.R, "Load data")
cross_sectional <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet")) %>%
  mutate(sex_category  = factor(sex_category,  levels = c("Male", "Female")),
         race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
         age10 = age_at_admission / 10, sf10 = sf_ratio / 10, log_pfvc = log(pfvc))
message("=== xsec_age_form_check: ", site_name, ", ", nrow(cross_sectional), " patients")

# --- 04: the exposure specifications and the covariate rule (04_analysis.R, 4b)
EXPOSURES_04 <- c(
  "VT/PFVC"               = "vtpfvc",
  "VT/PBW"                = "vtpbw",
  "VT/PFVC + VT/PBW"      = "vtpfvc + vtpbw",
  "VT/PBW + PFVC"         = "vtpbw + pfvc",
  "VT/PBW + log PFVC"     = "vtpbw + log_pfvc",
  "VT/PBW + PBW/PFVC"     = "vtpbw + pbwpfvc",
  "VT/PBW + VT excess (mL)" = "vtpbw + vt_excess_ml")
OUTCOMES_04 <- c("Mortality" = "deceased", "Elastance" = "ers", "Compliance" = "crs", "Static DP" = "dp")
DP_DERIVED <- c("dp", "ers", "crs", "ers_pbw", "ers_pfvc", "mechanical_power", "mp_pbw", "mp_pfvc",
                "mp_crs", "mp_elastic", "mp_elastic_pbw", "mp_elastic_pfvc")
covariates_04 <- function(outcome, exposure, age_term) {
  vars <- trimws(unlist(strsplit(paste(outcome, exposure, sep = " + "), "\\+")))
  paste("race_category +", age_term, "+ sex_category + sofa_total + sf10",
        if (any(vars %in% DP_DERIVED)) "+ bmi" else "")
}

# --- 05: the normalization head-to-head, on 05's own frames (05, PART 2)
base_05 <- cross_sectional %>%
  filter(!is.na(dp), dp > 0, !is.na(pfvc), pfvc > 0, !is.na(pbw), pbw > 0,
         !is.na(vtpbw), !is.na(bmi), !is.na(sofa_total), !is.na(sf_ratio), !is.na(deceased)) %>%
  mutate(pbwpfvc = pbw / pfvc)
DATA_05 <- list(
  "Elastance"        = base_05 %>% filter(!is.na(ers), ers > 0) %>% mutate(ers_pbw = ers * pbw, ers_pfvc = ers * pfvc),
  "Mechanical power" = base_05 %>% filter(!is.na(mechanical_power), mechanical_power > 0,
                                          !is.na(mp_pbw), !is.na(mp_pfvc)))
SPECS_05 <- tribble(
  ~family,            ~spec,                    ~exposure,
  "Elastance",        "Mechanic only",          "log(ers)",
  "Elastance",        "PBW-locked (Goligher)",  "log(ers_pbw)",
  "Elastance",        "PFVC-locked",            "log(ers_pfvc)",
  "Elastance",        "Separate (Ers + PBW)",   "log(ers) + log(pbw)",
  "Elastance",        "Separate (Ers + PFVC)",  "log(ers) + log(pfvc)",
  "Mechanical power", "Mechanic only",          "log(mechanical_power)",
  "Mechanical power", "PBW-locked (Gattinoni)", "log(mp_pbw)",
  "Mechanical power", "PFVC-locked",            "log(mp_pfvc)",
  "Mechanical power", "Separate (MP + PBW)",    "log(mechanical_power) + log(pbw)",
  "Mechanical power", "Separate (MP + PFVC)",   "log(mechanical_power) + log(pfvc)")
covariates_05 <- function(age_term) paste("vtpbw + sofa_total + sf10 + bmi +", age_term, "+ sex_category + race_category")

# --- one fit, reduced to what is compared: its exposure terms and its AIC. The
#     exposure's own terms are the ones compared; the covariates are not reported.
exposure_terms <- function(model, exposure) {
  wanted <- trimws(strsplit(exposure, "\\+")[[1]])
  co <- summary(model)$coefficients
  co <- co[rownames(co) %in% wanted, , drop = FALSE]
  tibble(term = rownames(co), estimate = co[, 1], se = co[, 2])
}
fit_one <- function(source, outcome, spec, exposure, formula, data, family) {
  model <- if (family == "binomial") glm(as.formula(formula), data = data, family = binomial)
           else lm(as.formula(formula), data = data)
  list(terms = exposure_terms(model, exposure) %>% mutate(source = source, outcome = outcome, spec = spec, .before = 1),
       fit = tibble(source = source, outcome = outcome, spec = spec, aic = AIC(model), n = stats::nobs(model)))
}

fits <- list()
for (age_form in names(AGE_FORMS)) {
  age_term <- AGE_FORMS[[age_form]]
  for (outcome_label in names(OUTCOMES_04)) {
    y <- OUTCOMES_04[[outcome_label]]
    for (spec in names(EXPOSURES_04)) {
      x <- EXPOSURES_04[[spec]]
      f <- paste(y, "~", x, "+", covariates_04(y, x, age_term))
      fits[[length(fits) + 1L]] <- c(fit_one("04", outcome_label, spec, x, f, cross_sectional,
                                             if (y == "deceased") "binomial" else "gaussian"),
                                     list(age_form = age_form))
    }
  }
  for (i in seq_len(nrow(SPECS_05))) {
    s <- SPECS_05[i, ]
    f <- paste("deceased ~", s$exposure, "+", covariates_05(age_term))
    fits[[length(fits) + 1L]] <- c(fit_one("05", paste("Mortality,", s$family), s$spec, s$exposure, f,
                                           DATA_05[[s$family]], "binomial"),
                                   list(age_form = age_form))
  }
}
terms <- map_dfr(fits, ~ .x$terms %>% mutate(age_form = .x$age_form))
fit   <- map_dfr(fits, ~ .x$fit   %>% mutate(age_form = .x$age_form))

# Both forms fit the same rows or the comparison is not like with like: a spline and
# a line drop no different patients, so a mismatch here is a bug, not a caveat.
n_check <- fit %>% select(source, outcome, spec, age_form, n) %>% pivot_wider(names_from = age_form, values_from = n)
stopifnot("the two age forms fitted different patients" = all(n_check$linear == n_check$spline))

# --- coefficients: linear beside spline, the move in the spline fit's SEs
coefficients <- terms %>%
  pivot_wider(id_cols = c(source, outcome, spec, term), names_from = age_form,
              values_from = c(estimate, se)) %>%
  mutate(lo_linear = estimate_linear - 1.96 * se_linear, hi_linear = estimate_linear + 1.96 * se_linear,
         lo_spline = estimate_spline - 1.96 * se_spline, hi_spline = estimate_spline + 1.96 * se_spline,
         change = estimate_spline - estimate_linear,
         change_in_spline_se = change / se_spline,
         pct_change = 100 * change / abs(estimate_linear),
         site = site_name)

# --- fit and ranking: the spline's AIC gain per model, and each outcome's order of
#     exposure specifications under each form
ranking <- fit %>%
  group_by(source, outcome, age_form) %>%
  mutate(rank = rank(aic, ties.method = "min"), delta_aic_within_outcome = aic - min(aic)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(source, outcome, spec, n), names_from = age_form,
              values_from = c(aic, rank, delta_aic_within_outcome)) %>%
  mutate(spline_minus_linear_aic = aic_spline - aic_linear,
         rank_changed = rank_linear != rank_spline, site = site_name)

# --- the comparisons the manuscript's claims rest on, as AIC gaps under each age form.
#     A rank change between two specifications a fraction of an AIC point apart means
#     nothing; a claimed winner's gap shrinking toward zero, or flipping, does.
CONTRASTS <- tribble(
  ~source, ~outcome,                      ~better,                  ~worse,
  "04",    "Mortality",                   "VT/PBW + log PFVC",      "VT/PFVC",
  "04",    "Mortality",                   "VT/PBW + log PFVC",      "VT/PBW",
  "04",    "Elastance",                   "VT/PBW + log PFVC",      "VT/PFVC",
  "04",    "Compliance",                  "VT/PBW + log PFVC",      "VT/PFVC",
  "05",    "Mortality, Mechanical power", "PFVC-locked",            "PBW-locked (Gattinoni)",
  "05",    "Mortality, Mechanical power", "Separate (MP + PFVC)",   "Separate (MP + PBW)",
  "05",    "Mortality, Elastance",        "PFVC-locked",            "PBW-locked (Goligher)",
  "05",    "Mortality, Elastance",        "Separate (Ers + PFVC)",  "Separate (Ers + PBW)")
aic_of <- function(src, out, sp, form) ranking %>% filter(source == src, outcome == out, spec == sp) %>%
  pull(paste0("aic_", form))
contrasts <- CONTRASTS %>% rowwise() %>%
  mutate(gap_linear = aic_of(source, outcome, better, "linear") - aic_of(source, outcome, worse, "linear"),
         gap_spline = aic_of(source, outcome, better, "spline") - aic_of(source, outcome, worse, "spline")) %>%
  ungroup() %>%
  # a share of a gap that was never there is noise: below 2 AIC points under linear
  # age there is no claimed winner to retain, and the share is left blank
  mutate(still_better = gap_spline < 0,
         gap_retained_pct = if_else(abs(gap_linear) >= 2, round(100 * gap_spline / gap_linear), NA_real_),
         site = site_name)

write_csv(coefficients, file.path(final_dir, paste0("age_form_coefficients_", site_name, ".csv")))
write_csv(ranking,      file.path(final_dir, paste0("age_form_fit_", site_name, ".csv")))
write_csv(contrasts,    file.path(final_dir, paste0("age_form_contrasts_", site_name, ".csv")))

# --- what a reader needs first: whether the claims survive, then the size terms that
#     moved most, then whether the curvature is there at all
options(width = 200)
message("\nThe claimed comparisons: AIC of the better specification minus the worse, under each age form")
message("(negative = still better; gap_retained_pct = how much of the linear-age gap survives the spline)")
print(as.data.frame(contrasts %>% transmute(source, outcome, better, worse,
                                            gap_linear = round(gap_linear, 1), gap_spline = round(gap_spline, 1),
                                            still_better, gap_retained_pct)), row.names = FALSE)
size_terms <- coefficients %>% filter(str_detect(term, "pfvc|pbw"))
message("\nPredicted-lung-size terms, linear age -> spline age (largest moves first, in the spline fit's SEs):")
print(as.data.frame(size_terms %>% arrange(desc(abs(change_in_spline_se))) %>% head(12) %>%
  transmute(source, outcome, spec, term, linear = signif(estimate_linear, 3),
            spline = signif(estimate_spline, 3), moved_se = round(change_in_spline_se, 2))), row.names = FALSE)
message("\nSpline age's AIC gain over linear age, by outcome (negative = the curvature is real):")
print(as.data.frame(ranking %>% group_by(source, outcome) %>%
  summarise(median_gain = round(median(spline_minus_linear_aic), 1), .groups = "drop")), row.names = FALSE)
message("\n=== xsec_age_form_check complete -> ", final_dir)
