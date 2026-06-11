# =============================================================================
# Exploratory: stress-strain mortality surface (predicted risk of death)
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Exploratory, standalone analysis (NOT part of the 00 pipeline runner).
#
# In the Gattinoni stress-strain framework, lung injury is driven by the joint
# action of STRAIN (VT / aerated volume) and STRESS (transpulmonary driving
# pressure). This reproduces the original paper's view: model in-hospital mortality
# on the strain and stress estimates together, then plot the predicted risk of
# death over the strain-stress plane (marginal effects), strain on x, stress on y.
#
#   * Strain estimate = VT/PFVC. VT/PBW is deliberately NOT used: it is confounded
#     by severity of illness (clinicians lower set VT in sicker patients), whereas
#     VT/PFVC is the strain estimate and PBW/PFVC (below) is the purely demographic
#     size-surrogate discrepancy -- the two exposures used in the paper.
#   * Stress estimate = driving pressure (DP = Pplat - PEEP), the airway-pressure
#     surrogate for transpulmonary stress (= E_spec x strain). Requires a measured
#     plateau.
#   * PBW/PFVC is included in the model (the size-surrogate discrepancy, which is
#     demographic and not severity-confounded) and held at its typical value for
#     the surface.
#
# Strain and stress both scale with VT, so they are correlated; the surface is only
# trustworthy where data exist (observed points are overlaid), and off-support
# corners (e.g. high strain + low stress) are extrapolation. Adjusted for the
# standard covariates + BMI; covariates held at typical values for the grid.
#
# Inputs : output/<site>/intermediate/analysis_cross_sectional.parquet (script 03)
# Outputs: final/stress_strain_mortality_<site>.pdf    (predicted-risk surface)
#          final/stress_strain_models_<site>.csv         (strain / stress / PBW:PFVC ORs)
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(broom)
library(marginaleffects)

source("utils/config.R")
site_name <- config$site_name

output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

cross_sectional <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"))

MIN_N <- 30  # surface fitting needs a reasonable sample

# Stress (DP) requires a measured plateau, so restrict to the mechanics subset.
mech <- cross_sectional %>%
  filter(is.finite(dp), dp > 0, is.finite(vtpfvc), is.finite(pbwpfvc),
         !is.na(deceased), is.finite(sofa_total), is.finite(sf_ratio),
         is.finite(bmi)) %>%
  mutate(
    sex_category  = factor(sex_category,  levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
    age10 = age_at_admission / 10,
    sf10  = sf_ratio / 10
  )

message("Mechanics subset (finite DP + covariates): ", nrow(mech), " subjects, ",
        sum(mech$deceased), " deaths.")
if (nrow(mech) < MIN_N || length(unique(mech$deceased)) < 2) {
  stop("Too few subjects (", nrow(mech), ") or no mortality variation for the ",
       "stress-strain surface.")
}

message("cor(VT/PFVC, DP) = ", round(cor(mech$vtpfvc, mech$dp), 2),
        " (strain and stress share VT)")

# =============================================================================
# Strain + stress mortality model (PBW/PFVC and covariates adjusted)
# =============================================================================
# DP is a driving-pressure derivative, so adjust for BMI as elsewhere.
adjustment <- "pbwpfvc + age10 + sex_category + race_category + sofa_total + sf10 + bmi"
m <- glm(as.formula(paste("deceased ~ vtpfvc + dp +", adjustment)),
         data = mech, family = binomial)

coefs <- broom::tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term %in% c("vtpfvc", "dp", "pbwpfvc")) %>%
  transmute(
    site = site_name,
    exposure = recode(term, vtpfvc = "strain (VT/PFVC)", dp = "stress (DP)",
                      pbwpfvc = "PBW/PFVC"),
    odds_ratio = estimate, conf_low = conf.low, conf_high = conf.high,
    p_value = p.value, n_obs = stats::nobs(m), n_deaths = sum(mech$deceased)
  )
write_csv(coefs, file.path(final_dir, paste0("stress_strain_models_", site_name, ".csv")))

message("Odds ratios:")
pwalk(coefs, function(exposure, odds_ratio, p_value, ...) {
  message("  ", exposure, ": OR = ", round(odds_ratio, 2), " (p = ", signif(p_value, 2), ")")
})

# =============================================================================
# Predicted-risk surface over strain (x) and stress (y)
# =============================================================================
# Grid over the 5th-95th percentile of each axis; PBW/PFVC and covariates held at
# their means / modes.
strain_seq <- seq(quantile(mech$vtpfvc, 0.05), quantile(mech$vtpfvc, 0.95), length.out = 40)
dp_seq     <- seq(quantile(mech$dp, 0.05),     quantile(mech$dp, 0.95),     length.out = 40)
grid  <- datagrid(model = m, vtpfvc = strain_seq, dp = dp_seq)
preds <- predictions(m, newdata = grid) %>% as_tibble()

surface_plot <- ggplot(preds, aes(vtpfvc, dp, fill = estimate)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = estimate), color = "white", linewidth = 0.3, alpha = 0.6) +
  # observed subjects (where the surface is actually supported)
  geom_point(data = mech, aes(x = vtpfvc, y = dp), inherit.aes = FALSE,
             alpha = 0.10, size = 0.5, color = "white") +
  scale_fill_viridis_c(name = "Predicted\nmortality", limits = c(0, NA)) +
  labs(
    title = "Predicted mortality over the stress-strain plane",
    subtitle = paste0(site_name,
      " — strain (VT/PFVC) on x, stress (driving pressure) on y; PBW/PFVC and ",
      "covariates adjusted; white points = observed subjects (else extrapolation)"),
    x = "Strain estimate: VT/PFVC",
    y = "Stress estimate: driving pressure (cmH2O)"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(final_dir, paste0("stress_strain_mortality_", site_name, ".pdf")),
       surface_plot, width = 8, height = 6)

message("Stress-strain mortality surface written.")
message("Exploratory stress-strain mortality analysis complete.")
