# =============================================================================
# Exploratory: age functional-form sensitivity (linear vs spline age)
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Exploratory, standalone analysis (NOT part of the 00 pipeline runner).
#
# Sensitivity analysis for the modeling choice: does entering AGE linearly (as in
# the main models) vs flexibly (a spline) change the observed exposure signal? Age
# is a genuine confounder (it has non-lung paths to mortality), so the relevant
# robustness check is whether a more flexible age adjustment absorbs the signal.
# Two specifications, BOTH WITHOUT height:
#   1. Linear age -- race + age10 + sex + SOFA + SF
#   2. Spline age -- race + ns(age, 4) + sex + SOFA + SF
# Height is deliberately NOT adjusted for in either: in the DAG height has no path to
# mortality except through lung size (height -> lung size -> strain -> mortality), so
# it is a determinant of the exposure on the causal path, and conditioning on it is
# over-adjustment that removes the mechanism -- not a confounder to control.
#
# Read: a size-carrying exposure (PFVC, VT/PFVC) that survives the spline-age
# specification is robust to nonlinear age confounding. A term that collapses under
# spline age (e.g. the PBW/PFVC ratio, which cancels height and so carries mostly
# age) was largely age. The VIF table makes the same point from the collinearity
# side (a term's signal vanishes as it becomes collinear with the flexible age basis).
#
# Inputs : output/<site>/intermediate/analysis_cross_sectional.parquet (script 03)
# Outputs: final/spline_sensitivity_<site>.csv     (AIC / evidence ratios, both age forms)
#          final/spline_sensitivity_<site>.pdf      (heatmaps: linear vs spline age)
#          final/spline_sensitivity_vif_<site>.csv  (VIF of the size terms per age form)
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

# Two age specifications, both without height (height is over-adjustment; see header).
covar_levels <- c("Linear age", "Spline age")
build_covars <- function(outcome, level) {
  base <- switch(level,
    "Linear age" = "race_category + age10 + sex_category + sofa_total + sf10",
    "Spline age" = paste("race_category + ns(age_at_admission, 4) + sex_category +",
                         "sofa_total + sf10"))
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

message("Evidence ratio (vs VT/PBW) for VT/PBW + PBW/PFVC, linear vs spline age:")
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
    title = "Age functional-form sensitivity: linear vs spline age (no height)",
    subtitle = paste0(site_name,
      " — an exposure that survives the spline-age panel is robust to nonlinear age ",
      "confounding; a term that collapses (e.g. the PBW/PFVC ratio) was largely age"),
    x = "Outcome", y = "Exposure specification") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(final_dir, paste0("spline_sensitivity_", site_name, ".pdf")),
       heat, width = 11, height = 6)

# =============================================================================
# Variance inflation across the two age specifications
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
message("Variance inflation (VIF) of the size terms, linear vs spline age:")
pwalk(vif_tbl, function(level, vif_PFVC, vif_PBWPFVC, vif_VTPFVC, ...) {
  message("  ", level, ": PFVC=", round(vif_PFVC, 1),
          ", PBW/PFVC=", round(vif_PBWPFVC, 1), ", VT/PFVC=", round(vif_VTPFVC, 1))
})

message("Spline-sensitivity outputs written.")
message("Exploratory spline-sensitivity analysis complete.")
