# =============================================================================
# Exploratory: does the VT/PBW + PFVC fit advantage survive height adjustment?
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Exploratory, standalone analysis (NOT part of the 00 pipeline runner).
#
# In the script-04 AIC heatmap, "VT/PBW + PFVC" fits dramatically better than every
# other exposure specification (evidence ratio >1000) for nearly every outcome.
# Hypothesis: this is largely an OMITTED-COVARIATE effect -- PFVC = GLI(age, height,
# sex, race) injects HEIGHT (a strong body/lung-size axis absent from the covariate
# set), whereas the PBW/PFVC ratio cancels height and keeps only the relative
# discrepancy. If so, adding height to the covariates should collapse the
# "VT/PBW + PFVC" vs "VT/PBW + PBW/PFVC" gap.
#
# This refits the AIC comparison (Mortality, 28-day VFD competing-risks, Compliance,
# Elastance, Static DP, Mechanical power x the five exposure specs) TWICE: with the
# standard covariates, and with height added. Height is present for every cohort
# subject (cohort height range 150-210 cm), so the two fits are on the SAME rows
# within each outcome -- the AIC change is purely the height covariate.
#
# Inputs : output/<site>/intermediate/analysis_cross_sectional.parquet (script 03)
# Outputs: final/height_sensitivity_<site>.csv     (AIC / evidence ratios, both covsets)
#          final/height_sensitivity_<site>.pdf      (side-by-side heatmaps)
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(survival)

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
    title = "Height-adjustment sensitivity of the evidence-ratio heatmap",
    subtitle = paste0(site_name,
      " — if adding height collapses the VT/PBW + PFVC advantage, the >1000 was ",
      "largely omitted size adjustment"),
    x = "Outcome", y = "Exposure specification"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(final_dir, paste0("height_sensitivity_", site_name, ".pdf")),
       heat, width = 13, height = 6)

message("Height-sensitivity outputs written.")
message("Exploratory height-sensitivity analysis complete.")
