# =============================================================================
# Script 04c: Mechanical power normalized to PBW vs PFVC
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Just as tidal volume and elastance can be scaled by predicted body weight (PBW)
# or predicted FVC (PFVC), so can mechanical power (MP). This script compares MP
# normalized to PBW (J/min/kg) against MP normalized to PFVC (J/min/L) as a
# predictor of in-hospital mortality. MP itself is derived per index timepoint in
# script 03 using a mode-aware simplified power equation (volume control uses the
# Gattinoni driving-pressure form; pressure-targeted modes use the rectangular
# Becher form), so MP is defined only for controlled modes with recorded peak /
# plateau pressures.
#
# Inputs : output/<site>/intermediate/analysis_cross_sectional.parquet (script 03)
# Outputs: final/regression_mp_normalization_<site>.html  (model table)
#          final/mp_normalization_<site>.csv              (poolable summary)
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(gtsummary)
library(gt)
library(broom)

source("utils/config.R")
site_name <- config$site_name

output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

cross_sectional <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"))
message("Cross-sectional: ", nrow(cross_sectional), " observations")

# =============================================================================
# Prepare variables (mirrors script 04)
# =============================================================================
# age10 / sf10 are per-10-unit; mp_pbw / mp_pfvc are z-scaled (per SD) so the two
# normalizations are on a common scale and their odds ratios are comparable
# despite the different raw units (J/min/kg vs J/min/L). Standardizing is a linear
# rescaling, so model fit / AIC are unchanged.
zscore <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

model_data <- cross_sectional %>%
  mutate(
    sex_category  = factor(sex_category,  levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
    age10     = age_at_admission / 10,
    sf10      = sf_ratio / 10,
    mp_pbw_z  = zscore(mp_pbw),
    mp_pfvc_z = zscore(mp_pfvc)
  )

n_mp <- sum(!is.na(model_data$mechanical_power))
message("Encounters with a computable mechanical power: ", n_mp,
        " of ", nrow(model_data))
if (n_mp == 0) {
  stop("No encounters have a computable mechanical power. Check that the ",
       "respiratory_support table provides peak_inspiratory_pressure_obs, ",
       "plateau_pressure_obs, resp_rate_set, and a controlled mode_category.")
}

# VT/PBW is retained as a severity-of-illness adjustment (as in the script 04
# elastance models); BMI is added because MP is a driving-pressure derivative.
adjustment <- "vtpbw + race_category + age10 + sex_category + sofa_total + sf10 + bmi"

mp_specs <- c(mp_pbw_z = "MP / PBW (per SD)", mp_pfvc_z = "MP / PFVC (per SD)")

# =============================================================================
# Fit the two mortality models
# =============================================================================
has_mortality_variation <- length(unique(na.omit(model_data$deceased))) > 1
if (!has_mortality_variation) {
  stop("Cannot fit mechanical-power mortality models: outcome `deceased` has no ",
       "variation (all = ", unique(na.omit(model_data$deceased)), ").")
}

mp_models <- list()
for (v in names(mp_specs)) {
  fstr <- paste("deceased ~", v, "+", adjustment)
  mp_models[[v]] <- glm(as.formula(fstr), data = model_data, family = binomial)
}

# =============================================================================
# AIC / evidence-ratio comparison (referenced to the PBW-scaled model)
# =============================================================================
# MP/PBW and MP/PFVC are fit on the same rows (both require mechanical_power), so
# their AICs are directly comparable. The evidence ratio is referenced to MP/PBW
# (the PBW-scaled model), matching the PBW-reference convention used elsewhere:
# ER = exp(-1/2 * (AIC - AIC_MP/PBW)); ER > 1 favors MP/PFVC over MP/PBW.
aic_pbw <- AIC(mp_models[["mp_pbw_z"]])

mp_summary <- imap_dfr(mp_models, function(m, v) {
  or_row <- broom::tidy(m, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == v)
  tibble(
    site           = site_name,
    exposure       = mp_specs[[v]],
    n_obs          = stats::nobs(m),
    or_per_sd      = or_row$estimate,
    conf_low       = or_row$conf.low,
    conf_high      = or_row$conf.high,
    p_value        = or_row$p.value,
    aic            = AIC(m),
    delta_aic      = AIC(m) - aic_pbw,
    evidence_ratio = exp(-0.5 * (AIC(m) - aic_pbw))
  )
})

write_csv(mp_summary, file.path(final_dir, paste0("mp_normalization_", site_name, ".csv")))

message("Mechanical power normalization (mortality, OR per SD; ",
        "evidence ratio vs MP/PBW):")
pwalk(mp_summary, function(exposure, or_per_sd, evidence_ratio, ...) {
  message("  ", exposure, ": OR = ", round(or_per_sd, 3),
          ", evidence ratio = ", signif(evidence_ratio, 3))
})

# =============================================================================
# Regression table (odds ratios)
# =============================================================================
mp_tables <- imap(mp_models, ~ {
  tbl_regression(
    .x, exponentiate = TRUE,
    label = c(
      setNames(list(mp_specs[[.y]]), .y),
      list(age10 = "Age (per 10 yr)", sf10 = "SF ratio (per 10)", vtpbw = "VT/PBW")
    )
  ) %>% bold_p()
})
tbl_merge(mp_tables, tab_spanner = unname(mp_specs)) %>%
  as_gt() %>%
  gt::gtsave(file.path(final_dir,
                       paste0("regression_mp_normalization_", site_name, ".html")))

message("Mechanical-power normalization table written")
message("Script 04c complete.")
