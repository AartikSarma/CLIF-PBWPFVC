# =============================================================================
# Exploratory: does PBW/PFVC survive FLEXIBLE demographic adjustment?
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Exploratory, standalone analysis (NOT part of the 00 pipeline runner).
#
# IMPORTANT (causal framing): this is an IDENTIFIABILITY / collinearity diagnostic,
# NOT a sequence of better-specified models. Two things to keep straight:
#   * HEIGHT is NOT a confounder (no path to mortality except via lung size), so
#     adjusting for it is OVER-ADJUSTMENT that removes the lung-size mechanism. The
#     "+ height" and spline-height steps are diagnostic, not preferred models.
#   * AGE and RACE are the genuine confounders (they have non-lung paths to
#     mortality), but they are ALSO determinants of lung size -- i.e. they sit on
#     BOTH the confounding path and the causal (size) path. So flexibly adjusting for
#     them simultaneously removes confounding AND removes the age/race-driven size
#     mechanism. The collapse below therefore OVERSTATES confounding; you cannot
#     separate the two by adjustment, because the confounders ARE the size mechanism.
#
# The height-sensitivity check showed "VT/PBW + PFVC"'s advantage was largely height,
# while "VT/PBW + PBW/PFVC" (which cancels height) survived. PBW/PFVC is itself a
# deterministic function of age/height/sex/race, entered only linearly/categorically.
#
# This walks a covariate-flexibility ladder and re-checks the evidence ratios:
#   1. Linear (base)             -- race + age10 + sex + SOFA + SF
#   2. + spline age (NO height)  -- race + ns(age,4) + sex + SOFA + SF   <-- KEY rung
#   3. + spline age & height     -- ...+ ns(height,4)   (over-adjusted contrast only)
# Rung 2 is the causally-correct specification: it flexibly adjusts the genuine
# confounder (age) without over-adjusting for height (a non-confounder on the causal
# path -- height -> lung size -> strain -> mortality, with no other path). If the
# PFVC / PBW-PFVC signal SURVIVES rung 2, it is robust to nonlinear age confounding
# (the real concern). Rung 3 is included only to show that adding height collapses it
# via over-adjustment, NOT to argue the signal is spurious.
#
# Also includes a quick elastance-anomaly check: refit the lone surviving >1000 cell
# (VT/PBW + PFVC for Elastance) on the log scale to test whether it was a heavy-tail
# / scale artifact of fitting elastance (= 1000/compliance) linearly.
#
# Inputs : output/<site>/intermediate/analysis_cross_sectional.parquet (script 03)
# Outputs: final/spline_sensitivity_<site>.csv     (AIC / evidence ratios, all levels)
#          final/spline_sensitivity_<site>.pdf      (heatmaps across the covariate ladder)
#          final/spline_sensitivity_vif_<site>.csv  (VIF of the size terms per level)
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(survival)
library(splines)

source("utils/config.R")
site_name <- config$site_name

output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

cross_sectional <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"))

mech <- cross_sectional %>%
  mutate(
    sex_category  = factor(sex_category,  levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
    age10    = age_at_admission / 10,
    height10 = height_cm / 10,
    sf10     = sf_ratio / 10
  )

specs <- c(
  "VT/PFVC"           = "vtpfvc",
  "VT/PBW"            = "vtpbw",
  "VT/PFVC + VT/PBW"  = "vtpfvc + vtpbw",
  "VT/PBW + PFVC"     = "vtpbw + pfvc",
  "VT/PBW + PBW/PFVC" = "vtpbw + pbwpfvc"
)
outcomes <- list(
  list(key = "Mortality",        var = "deceased",         type = "logistic", bmi = FALSE),
  list(key = "28-day VFDs",      var = NA,                 type = "finegray", bmi = FALSE),
  list(key = "Compliance",       var = "crs",              type = "linear",   bmi = TRUE),
  list(key = "Elastance",        var = "ers",              type = "linear",   bmi = TRUE),
  list(key = "Static DP",        var = "dp",               type = "linear",   bmi = TRUE),
  list(key = "Mechanical power", var = "mechanical_power", type = "linear",   bmi = TRUE)
)

# The KEY rung is "+ spline age (no height)": it flexibly adjusts the genuine
# confounder (age, which has non-lung paths to mortality) WITHOUT over-adjusting for
# height (a non-confounder on the causal path). "+ spline age & height" is shown only
# as the over-adjusted contrast -- the combined panel cannot tell whether a collapse
# is due to proper age adjustment or improper height over-adjustment.
covar_levels <- c("Linear (base)", "+ spline age (no height)", "+ spline age & height")
build_covars <- function(outcome, level) {
  base <- switch(level,
    "Linear (base)"            = "race_category + age10 + sex_category + sofa_total + sf10",
    "+ spline age (no height)" = paste("race_category + ns(age_at_admission, 4) +",
                                       "sex_category + sofa_total + sf10"),
    "+ spline age & height"    = paste("race_category + ns(age_at_admission, 4) + sex_category +",
                                       "sofa_total + sf10 + ns(height_cm, 4)"))
  if (outcome$bmi) base <- paste(base, "+ bmi")
  base
}

aic_one <- function(outcome, spec_rhs, level) {
  rhs <- paste(spec_rhs, "+", build_covars(outcome, level))
  if (outcome$type == "logistic") {
    AIC(glm(as.formula(paste("deceased ~", rhs)), data = mech, family = binomial))
  } else if (outcome$type == "linear") {
    AIC(lm(as.formula(paste(outcome$var, "~", rhs)), data = mech))
  } else { # competing-risks Fine-Gray for 28-day VFDs
    needed <- all.vars(as.formula(paste("~", rhs)))
    df <- mech %>%
      mutate(vfd_status_f = factor(vfd_status, levels = c(0, 1, 2),
                                   labels = c("censored", "extubation", "death"))) %>%
      select(vfd_time, vfd_status_f, all_of(needed)) %>%
      filter(!is.na(vfd_time), vfd_time > 0, !is.na(vfd_status_f)) %>%
      drop_na(all_of(needed))
    fg <- finegray(Surv(vfd_time, vfd_status_f) ~ ., data = df, etype = "extubation")
    AIC(coxph(as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ ", rhs)),
              weights = fgwt, data = fg))
  }
}

# =============================================================================
# Evidence ratios across the covariate-flexibility ladder
# =============================================================================
aic_df <- expand_grid(level = covar_levels, outcome_i = seq_along(outcomes),
                      spec = names(specs)) %>%
  mutate(
    outcome = map_chr(outcome_i, ~ outcomes[[.x]]$key),
    aic = pmap_dbl(list(outcome_i, spec, level), function(oi, sp, lv) {
      aic_one(outcomes[[oi]], specs[[sp]], lv)
    })
  ) %>%
  group_by(level, outcome) %>%
  mutate(
    evidence_ratio = exp(-0.5 * (aic - aic[spec == "VT/PBW"])),
    er_trunc = pmin(pmax(evidence_ratio, 0.001), 1000),
    er_label = case_when(evidence_ratio > 1000 ~ ">1000",
                         evidence_ratio < 0.001 ~ "<0.001",
                         TRUE ~ formatC(evidence_ratio, format = "g", digits = 2))
  ) %>%
  ungroup() %>%
  mutate(
    level   = factor(level, levels = covar_levels),
    outcome = factor(outcome, levels = c("Mortality", "28-day VFDs", "Compliance",
                                         "Elastance", "Static DP", "Mechanical power")),
    spec    = factor(spec, levels = c("VT/PBW", "VT/PFVC", "VT/PFVC + VT/PBW",
                                      "VT/PBW + PBW/PFVC", "VT/PBW + PFVC"))
  ) %>%
  select(-outcome_i)

write_csv(aic_df, file.path(final_dir, paste0("spline_sensitivity_", site_name, ".csv")))

message("Evidence ratio (vs VT/PBW) for VT/PBW + PBW/PFVC across the ladder:")
aic_df %>%
  filter(spec == "VT/PBW + PBW/PFVC", outcome %in% c("Mortality", "28-day VFDs")) %>%
  arrange(outcome, level) %>%
  pwalk(function(outcome, level, evidence_ratio, ...) {
    message("  [", outcome, "] ", level, ": ", signif(evidence_ratio, 3))
  })

heat <- ggplot(aic_df, aes(outcome, spec, fill = log10(er_trunc))) +
  geom_tile(color = "grey80", linewidth = 0.4) +
  geom_text(aes(label = er_label), size = 2.8) +
  facet_wrap(~ level) +
  scale_fill_gradient2(
    name = "Evidence ratio\n(vs VT/PBW)",
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(log10(0.001), log10(1000)),
    breaks = -3:3, labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000")) +
  labs(
    title = "Flexible age adjustment without height over-adjustment",
    subtitle = paste0(site_name,
      " — middle panel (spline age, NO height) is the causally-correct adjustment: if ",
      "the PFVC/ratio signal survives there, it is robust to nonlinear age confounding"),
    x = "Outcome", y = "Exposure specification") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(final_dir, paste0("spline_sensitivity_", site_name, ".pdf")),
       heat, width = 16, height = 6)

# =============================================================================
# Elastance anomaly: raw vs log scale (linear + height covariate set)
# =============================================================================
el_covars  <- "race_category + age10 + sex_category + sofa_total + sf10 + height10 + bmi"
el_data    <- mech %>% filter(is.finite(ers), ers > 0)
er_for <- function(response) {
  aic_full <- AIC(lm(as.formula(paste(response, "~ vtpbw + pfvc +", el_covars)), data = el_data))
  aic_ref  <- AIC(lm(as.formula(paste(response, "~ vtpbw +", el_covars)), data = el_data))
  exp(-0.5 * (aic_full - aic_ref))
}
message("Elastance anomaly (VT/PBW + PFVC vs VT/PBW, +height covariates):")
message("  raw elastance:  evidence ratio = ", signif(er_for("ers"), 3))
message("  log elastance:  evidence ratio = ", signif(er_for("log(ers)"), 3))

# =============================================================================
# Variance inflation across the covariate-flexibility ladder
# =============================================================================
# The evidence-ratio collapse and variance inflation are the SAME phenomenon: as
# the covariates capture a size term's functional form, that term becomes redundant
# (collinear) and adds no independent information -- its evidence ratio falls AND its
# VIF rises. VIF for a single continuous predictor = 1 / (1 - R^2) from the auxiliary
# regression of that predictor on all the other predictors (clinical, no-BMI set;
# the auxiliary regression does not involve the outcome).
vif_term <- function(target, others) {
  r2 <- summary(lm(as.formula(paste(target, "~", others)), data = mech))$r.squared
  1 / (1 - r2)
}
vif_tbl <- map_dfr(covar_levels, function(lv) {
  covs <- build_covars(list(bmi = FALSE), lv)
  tibble(site = site_name, level = lv,
         vif_PFVC    = vif_term("pfvc",    paste("vtpbw +", covs)),
         vif_PBWPFVC = vif_term("pbwpfvc", paste("vtpbw +", covs)),
         vif_VTPFVC  = vif_term("vtpfvc",  paste("vtpbw +", covs)))
})
write_csv(vif_tbl, file.path(final_dir, paste0("spline_sensitivity_vif_", site_name, ".csv")))
message("Variance inflation (VIF) of the size terms across the ladder:")
pwalk(vif_tbl, function(level, vif_PFVC, vif_PBWPFVC, vif_VTPFVC, ...) {
  message("  ", level, ": PFVC=", round(vif_PFVC, 1),
          ", PBW/PFVC=", round(vif_PBWPFVC, 1), ", VT/PFVC=", round(vif_VTPFVC, 1))
})

message("Spline-sensitivity outputs written.")
message("Exploratory spline-sensitivity analysis complete.")
