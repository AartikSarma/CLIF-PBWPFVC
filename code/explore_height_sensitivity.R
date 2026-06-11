# =============================================================================
# Exploratory: height as a COLLINEARITY / OVER-ADJUSTMENT diagnostic (not a model)
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Exploratory, standalone analysis (NOT part of the 00 pipeline runner).
#
# IMPORTANT (causal framing): height is NOT a confounder of strain -> mortality. In
# the DAG, height affects mortality only through lung size (height -> lung size ->
# strain -> mortality); there is no plausible path from height to mortality that
# bypasses lung size. So height is a determinant of the exposure ON the causal path,
# and conditioning on it is OVER-ADJUSTMENT (it removes the lung-size variation the
# analysis is about), not confounder control. Height was therefore deliberately
# EXCLUDED from the script-04 models, correctly. This script is a DIAGNOSTIC: it
# shows what over-adjusting for height does to the evidence-ratio structure (it
# collapses the PFVC terms because they are largely height), which is a collinearity
# probe -- NOT a "better-specified" or preferred model. Do not read the +height
# panel as the right adjustment.
#
# Strain IS height-dependent: Devine PBW is ~linear in height while true lung volume
# scales ~height^2.5-3, so PBW/V0 falls with height and VT-to-PBW dosing delivers
# higher strain to shorter patients -- the allometric mismatch at the core of the
# prior bias work. So height's predictive signal here is the STRAIN MECHANISM
# (height -> PBW mis-sizes the lung -> wrong strain -> mortality), and conditioning
# on height removes that mechanism, not confounding.
#
# Refits the AIC comparison (Mortality, 28-day VFD competing-risks, Compliance,
# Elastance, Static DP, Mechanical power x the five exposure specs) with vs without
# height (same rows; height present for every subject). Also runs a height-mediation
# test (below): does height predict mortality CONDITIONAL on the strain estimate
# (VT/PFVC) and severity? Because strain is height-dependent, conditioning on VT/PFVC
# blocks the height -> strain path, so a RESIDUAL height effect would mean VT/PFVC
# under-captures the height allometry of true lung size (PFVC is an imperfect V0
# surrogate), or absolute size / a non-lung path -- NOT that height is a confounder.
#
# Inputs : output/<site>/intermediate/analysis_cross_sectional.parquet (script 03)
# Outputs: final/height_sensitivity_<site>.csv     (AIC / evidence ratios, both covsets)
#          final/height_sensitivity_<site>.pdf      (side-by-side heatmaps)
#          final/height_mediation_<site>.csv         (height | strain mortality test)
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

# Exposure specifications (same as script 04) and the heatmap outcomes (the
# non-normalized ones). BMI is added for the driving-pressure-derived outcomes.
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

base_covars <- "race_category + age10 + sex_category + sofa_total + sf10"

aic_one <- function(outcome, spec_rhs, with_height) {
  covars <- base_covars
  if (outcome$bmi)   covars <- paste(covars, "+ bmi")
  if (with_height)   covars <- paste(covars, "+ height10")
  rhs <- paste(spec_rhs, "+", covars)
  if (outcome$type == "logistic") {
    AIC(glm(as.formula(paste("deceased ~", rhs)), data = mech, family = binomial))
  } else if (outcome$type == "linear") {
    AIC(lm(as.formula(paste(outcome$var, "~", rhs)), data = mech))
  } else { # competing-risks Fine-Gray for 28-day VFDs
    rhs_vars <- unique(trimws(unlist(strsplit(rhs, "\\+"))))
    df <- mech %>%
      mutate(vfd_status_f = factor(vfd_status, levels = c(0, 1, 2),
                                   labels = c("censored", "extubation", "death"))) %>%
      select(vfd_time, vfd_status_f, all_of(rhs_vars)) %>%
      filter(!is.na(vfd_time), vfd_time > 0, !is.na(vfd_status_f))
    fg <- finegray(Surv(vfd_time, vfd_status_f) ~ ., data = df, etype = "extubation")
    AIC(coxph(as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ ", rhs)),
              weights = fgwt, data = fg))
  }
}

# =============================================================================
# Fit every outcome x exposure-spec x covariate-set and compute evidence ratios
# =============================================================================
aic_df <- expand_grid(
  covset = c("Without height", "With height"),
  outcome_i = seq_along(outcomes),
  spec = names(specs)
) %>%
  mutate(
    outcome = map_chr(outcome_i, ~ outcomes[[.x]]$key),
    aic = pmap_dbl(list(outcome_i, spec, covset), function(oi, sp, cs) {
      aic_one(outcomes[[oi]], specs[[sp]], cs == "With height")
    })
  ) %>%
  group_by(covset, outcome) %>%
  mutate(
    delta_aic = aic - aic[spec == "VT/PBW"],
    evidence_ratio = exp(-0.5 * delta_aic),
    er_trunc = pmin(pmax(evidence_ratio, 0.001), 1000),
    er_label = case_when(evidence_ratio > 1000 ~ ">1000",
                         evidence_ratio < 0.001 ~ "<0.001",
                         TRUE ~ formatC(evidence_ratio, format = "g", digits = 2))
  ) %>%
  ungroup() %>%
  mutate(
    covset  = factor(covset, levels = c("Without height", "With height")),
    outcome = factor(outcome, levels = c("Mortality", "28-day VFDs", "Compliance",
                                         "Elastance", "Static DP", "Mechanical power")),
    spec    = factor(spec, levels = c("VT/PBW", "VT/PFVC", "VT/PFVC + VT/PBW",
                                      "VT/PBW + PBW/PFVC", "VT/PBW + PFVC"))
  ) %>%
  select(-outcome_i)

write_csv(aic_df, file.path(final_dir, paste0("height_sensitivity_", site_name, ".csv")))

# Focused contrast: the two specs of interest, evidence ratio with vs without height.
key <- aic_df %>%
  filter(spec %in% c("VT/PBW + PFVC", "VT/PBW + PBW/PFVC")) %>%
  select(outcome, spec, covset, evidence_ratio) %>%
  pivot_wider(names_from = covset, values_from = evidence_ratio)
message("Evidence ratio (vs VT/PBW) for the two key specs, without vs with height:")
pwalk(key, function(outcome, spec, `Without height`, `With height`) {
  message("  [", outcome, "] ", spec, ": ", signif(`Without height`, 3),
          "  ->  ", signif(`With height`, 3))
})

# =============================================================================
# Side-by-side heatmaps (without height | with height)
# =============================================================================
heat <- ggplot(aic_df, aes(outcome, spec, fill = log10(er_trunc))) +
  geom_tile(color = "grey80", linewidth = 0.4) +
  geom_text(aes(label = er_label), size = 3) +
  facet_wrap(~ covset) +
  scale_fill_gradient2(
    name = "Evidence ratio\n(vs VT/PBW)",
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(log10(0.001), log10(1000)),
    breaks = -3:3, labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000")
  ) +
  labs(
    title = "Over-adjusting for height (a non-confounder): collinearity diagnostic",
    subtitle = paste0(site_name,
      " — height is NOT a confounder (no path to mortality except via lung size), so ",
      "the +height panel is over-adjustment that removes the PFVC mechanism, not a ",
      "better model"),
    x = "Outcome", y = "Exposure specification"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(final_dir, paste0("height_sensitivity_", site_name, ".pdf")),
       heat, width = 13, height = 6)

# =============================================================================
# Height-mediation test: does height predict mortality CONDITIONAL on strain?
# =============================================================================
# Strain is height-dependent (PBW's allometric mis-scaling), so conditioning on the
# strain estimate (VT/PFVC, flexible ns) blocks the height -> strain path. Under the
# DAG, height should then add little: a residual height effect means VT/PFVC
# under-captures the true height allometry of lung size (PFVC imperfectly estimates
# V0), or absolute size / a non-lung path -- it is NOT confounder control. Tested in
# two adjustment sets: minimal (strain + severity + sex) and + age + race.
mort <- cross_sectional %>%
  filter(!is.na(deceased), is.finite(vtpfvc), is.finite(sofa_total), is.finite(sf_ratio)) %>%
  mutate(
    sex_category  = factor(sex_category,  levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
    age10 = age_at_admission / 10, height10 = height_cm / 10, sf10 = sf_ratio / 10
  )

height_test <- function(label, base_rhs) {
  m0 <- glm(as.formula(paste("deceased ~", base_rhs)), data = mort, family = binomial)
  m1 <- glm(as.formula(paste("deceased ~", base_rhs, "+ height10")), data = mort, family = binomial)
  lrt <- anova(m0, m1, test = "LRT")
  ht  <- broom::tidy(m1, conf.int = TRUE, exponentiate = TRUE) %>% filter(term == "height10")
  tibble(site = site_name, adjustment = label, n_obs = stats::nobs(m1),
         height_OR_per10cm = ht$estimate, conf_low = ht$conf.low, conf_high = ht$conf.high,
         p_height = ht$p.value, lrt_p = lrt$`Pr(>Chi)`[2],
         er_add_height = exp(-0.5 * (AIC(m1) - AIC(m0))))
}
height_mediation <- bind_rows(
  height_test("strain + severity + sex",
              "ns(vtpfvc, 4) + sofa_total + sf10 + sex_category"),
  height_test("+ age + race",
              "ns(vtpfvc, 4) + sofa_total + sf10 + sex_category + age10 + race_category")
)
write_csv(height_mediation, file.path(final_dir, paste0("height_mediation_", site_name, ".csv")))

message("Height conditional on strain (VT/PFVC) + severity -- does height still predict mortality?")
pwalk(height_mediation, function(adjustment, height_OR_per10cm, p_height, er_add_height, ...) {
  message("  [", adjustment, "] height OR/10cm = ", round(height_OR_per10cm, 2),
          " (p = ", signif(p_height, 2), "), evidence ratio for adding height = ",
          signif(er_add_height, 3))
})

message("Height diagnostics written.")
message("Exploratory height diagnostic complete.")
