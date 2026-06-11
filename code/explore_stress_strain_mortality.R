# =============================================================================
# Exploratory: stress-strain mortality surface (predicted risk of death)
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Exploratory, standalone analysis (NOT part of the 00 pipeline runner).
#
# In the Gattinoni stress-strain framework, lung injury is driven by the joint
# action of STRAIN (VT / aerated volume) and STRESS (transpulmonary driving
# pressure). This reproduces the original paper's view: model in-hospital mortality
# on a strain estimate and a stress estimate together, then plot the predicted risk
# of death over the strain-stress plane (marginal effects), with strain on the
# x-axis and stress on the y-axis.
#
#   * Strain estimate = VT / predicted size. We make the surface for BOTH size
#     surrogates: VT/PBW and VT/PFVC. (Strain is purely geometric, so the strain
#     surrogate carries only the size-surrogate error -- the clean axis to compare
#     PBW vs PFVC.)
#   * Stress estimate = driving pressure (DP = Pplat - PEEP), the airway-pressure
#     surrogate for transpulmonary stress (= E_spec x strain). Shared across both
#     panels. Requires a measured plateau.
#
# Strain and stress both scale with VT, so they are correlated; the surface is only
# trustworthy where data exist (observed points are overlaid), and the off-support
# corners (e.g. high strain + low stress) are extrapolation. Adjusted for the
# standard covariates + BMI; covariates held at their typical values for the grid.
#
# Inputs : output/<site>/intermediate/analysis_cross_sectional.parquet (script 03)
# Outputs: final/stress_strain_mortality_<site>.pdf    (risk surfaces, VT/PBW & VT/PFVC)
#          final/stress_strain_models_<site>.csv         (strain & stress odds ratios)
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(patchwork)
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
  filter(is.finite(dp), dp > 0, is.finite(vtpbw), is.finite(vtpfvc),
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

# DP is a driving-pressure derivative, so adjust for BMI as elsewhere.
adjustment <- "age10 + sex_category + race_category + sofa_total + sf10 + bmi"

# Report the strain-stress correlation (they share VT) so the collinearity is explicit.
message("cor(VT/PBW, DP) = ", round(cor(mech$vtpbw, mech$dp), 2),
        "; cor(VT/PFVC, DP) = ", round(cor(mech$vtpfvc, mech$dp), 2))

# =============================================================================
# Fit a strain + stress mortality model and build the predicted-risk surface
# =============================================================================
build_surface <- function(strain_var, strain_label) {
  f <- as.formula(paste("deceased ~", strain_var, "+ dp +", adjustment))
  m <- glm(f, data = mech, family = binomial)

  # Grid over the 5th-95th percentile of each axis; covariates at their means/modes.
  s_seq  <- seq(quantile(mech[[strain_var]], 0.05), quantile(mech[[strain_var]], 0.95),
                length.out = 40)
  dp_seq <- seq(quantile(mech$dp, 0.05), quantile(mech$dp, 0.95), length.out = 40)
  grid <- do.call(datagrid, c(list(model = m, dp = dp_seq),
                              setNames(list(s_seq), strain_var)))
  preds <- predictions(m, newdata = grid) %>% as_tibble()
  preds$strain <- preds[[strain_var]]

  coefs <- broom::tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% c(strain_var, "dp")) %>%
    transmute(site = site_name, strain_axis = strain_label,
              term = if_else(term == strain_var, "strain (VT/size)", "stress (DP)"),
              odds_ratio = estimate, conf_low = conf.low, conf_high = conf.high,
              p_value = p.value, n_obs = stats::nobs(m), n_deaths = sum(mech$deceased))

  list(preds = preds, coefs = coefs, strain_var = strain_var, strain_label = strain_label)
}

surfaces <- list(
  build_surface("vtpbw",  "Strain estimate: VT/PBW (mL/kg)"),
  build_surface("vtpfvc", "Strain estimate: VT/PFVC")
)

# =============================================================================
# Outputs
# =============================================================================
models_tbl <- bind_rows(map(surfaces, "coefs"))
write_csv(models_tbl, file.path(final_dir, paste0("stress_strain_models_", site_name, ".csv")))

message("Strain / stress odds ratios:")
pwalk(models_tbl, function(strain_axis, term, odds_ratio, p_value, ...) {
  message("  [", strain_axis, "] ", term, ": OR = ", round(odds_ratio, 2),
          " (p = ", signif(p_value, 2), ")")
})

make_surface_plot <- function(s) {
  ggplot(s$preds, aes(strain, dp, fill = estimate)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = estimate), color = "white", linewidth = 0.3, alpha = 0.6) +
    # observed subjects (where the surface is actually supported)
    geom_point(data = mech, aes(x = .data[[s$strain_var]], y = dp),
               inherit.aes = FALSE, alpha = 0.10, size = 0.5, color = "white") +
    scale_fill_viridis_c(name = "Predicted\nmortality", limits = c(0, NA)) +
    labs(x = s$strain_label, y = "Stress estimate: driving pressure (cmH2O)") +
    theme_minimal(base_size = 11)
}

surface_fig <- (make_surface_plot(surfaces[[1]]) | make_surface_plot(surfaces[[2]])) +
  plot_annotation(
    title = "Predicted mortality over the stress-strain plane",
    subtitle = paste0(site_name,
      " — strain (VT/size) on x, stress (driving pressure) on y; white points are ",
      "observed subjects (surface is extrapolation away from them)"))

ggsave(file.path(final_dir, paste0("stress_strain_mortality_", site_name, ".pdf")),
       surface_fig, width = 13, height = 5.5)

message("Stress-strain mortality surfaces written.")
message("Exploratory stress-strain mortality analysis complete.")
