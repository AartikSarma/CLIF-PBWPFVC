# =============================================================================
# Exploratory: stress-strain mortality surfaces (predicted risk of death)
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Exploratory, standalone analysis (NOT part of the 00 pipeline runner).
#
# In the Gattinoni stress-strain framework, lung injury is driven jointly by STRAIN
# (VT / aerated volume) and STRESS (transpulmonary driving pressure). This
# reproduces the original paper's view: model in-hospital mortality on the strain
# estimate and a y-axis mechanics estimate together, then plot predicted risk of
# death over the plane (marginal effects), strain on x, the mechanics estimate on y.
#
#   * Strain estimate (x) = VT/PFVC. VT/PBW is NOT used as the strain axis: it is
#     confounded by severity of illness (clinicians lower set VT in sicker
#     patients).
#   * y-axis mechanics estimate, made for THREE quantities:
#       - driving pressure (DP = Pplat - PEEP), the stress estimate;
#       - Ers x PBW and Ers x PFVC, the two elastance normalizations (specific-
#         elastance / injury-extent estimates).
#   * Adjusted for VT/PBW -- the observed bedside dosing variable (VT/PFVC is NOT
#     observable at the bedside, so VT/PBW is the right confounder for the actual
#     ventilation decision) -- plus the standard covariates + BMI.
#
# NOTE on collinearity: VT/PFVC ~= VT/PBW x PBW/PFVC, so VT/PFVC, VT/PBW, and
# PBW/PFVC are algebraically linked; including all three is redundant. We keep
# VT/PFVC (x) and adjust for VT/PBW, which IMPLICITLY captures PBW/PFVC (sweeping
# VT/PFVC at fixed VT/PBW is sweeping the size-surrogate discrepancy). Strain and
# the y-axis quantities also share VT, so the surface is trustworthy only where data
# exist (observed points overlaid); off-support corners are extrapolation.
#
# Inputs : output/<site>/intermediate/analysis_cross_sectional.parquet (script 03)
# Outputs: final/stress_strain_mortality_<site>.pdf    (three predicted-risk surfaces)
#          final/stress_strain_models_<site>.csv         (exposure odds ratios)
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

# DP and the elastance normalizations require a measured plateau (mechanics subset).
mech <- cross_sectional %>%
  filter(is.finite(dp), dp > 0, is.finite(vtpfvc), is.finite(vtpbw),
         is.finite(ers_pbw), is.finite(ers_pfvc),
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
       "stress-strain surfaces.")
}

message("cor(VT/PFVC, VT/PBW) = ", round(cor(mech$vtpfvc, mech$vtpbw), 2),
        " (strain estimate vs the observed-dosing adjustment)")

# Adjust for observed dosing (VT/PBW) + covariates; BMI because DP / elastance are
# driving-pressure derivatives.
adjustment <- "vtpbw + age10 + sex_category + race_category + sofa_total + sf10 + bmi"

# =============================================================================
# Fit a strain + (y-axis mechanics) mortality model and build the risk surface
# =============================================================================
build_surface <- function(y_var, y_label) {
  m <- glm(as.formula(paste("deceased ~ vtpfvc +", y_var, "+", adjustment)),
           data = mech, family = binomial)

  strain_seq <- seq(quantile(mech$vtpfvc, 0.05), quantile(mech$vtpfvc, 0.95),
                    length.out = 40)
  y_seq <- seq(quantile(mech[[y_var]], 0.05), quantile(mech[[y_var]], 0.95),
               length.out = 40)
  grid <- do.call(datagrid, c(list(model = m, vtpfvc = strain_seq),
                              setNames(list(y_seq), y_var)))
  preds <- predictions(m, newdata = grid) %>% as_tibble()
  preds$yval <- preds[[y_var]]

  coefs <- broom::tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% c("vtpfvc", y_var, "vtpbw")) %>%
    transmute(site = site_name, y_axis = y_label,
              exposure = case_when(term == "vtpfvc" ~ "strain (VT/PFVC)",
                                   term == "vtpbw"  ~ "VT/PBW (adjustment)",
                                   TRUE             ~ y_label),
              odds_ratio = estimate, conf_low = conf.low, conf_high = conf.high,
              p_value = p.value, n_obs = stats::nobs(m), n_deaths = sum(mech$deceased))

  list(preds = preds, coefs = coefs, y_var = y_var, y_label = y_label)
}

surfaces <- list(
  build_surface("dp",       "Stress: driving pressure (cmH2O)"),
  build_surface("ers_pbw",  "Elastance x PBW"),
  build_surface("ers_pfvc", "Elastance x PFVC")
)

# =============================================================================
# Outputs
# =============================================================================
models_tbl <- bind_rows(map(surfaces, "coefs"))
write_csv(models_tbl, file.path(final_dir, paste0("stress_strain_models_", site_name, ".csv")))

message("Odds ratios (each model: strain VT/PFVC + y-axis quantity + VT/PBW + covars):")
pwalk(models_tbl, function(y_axis, exposure, odds_ratio, p_value, ...) {
  message("  [y=", y_axis, "] ", exposure, ": OR = ", round(odds_ratio, 2),
          " (p = ", signif(p_value, 2), ")")
})

make_surface_plot <- function(s) {
  ggplot(s$preds, aes(vtpfvc, yval, fill = estimate)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = estimate), color = "white", linewidth = 0.3, alpha = 0.6) +
    geom_point(data = mech, aes(x = vtpfvc, y = .data[[s$y_var]]), inherit.aes = FALSE,
               alpha = 0.08, size = 0.4, color = "white") +
    scale_fill_viridis_c(name = "Predicted\nmortality", limits = c(0, NA)) +
    labs(x = "Strain estimate: VT/PFVC", y = s$y_label) +
    theme_minimal(base_size = 10)
}

surface_fig <- (make_surface_plot(surfaces[[1]]) |
                make_surface_plot(surfaces[[2]]) |
                make_surface_plot(surfaces[[3]])) +
  plot_annotation(
    title = "Predicted mortality: strain (VT/PFVC) vs stress / normalized elastance",
    subtitle = paste0(site_name,
      " — strain (VT/PFVC) on x; y is driving pressure, then Ers x PBW, Ers x PFVC. ",
      "All models adjusted for VT/PBW (observed dosing) + covariates. White points = ",
      "observed subjects (else extrapolation)."))

ggsave(file.path(final_dir, paste0("stress_strain_mortality_", site_name, ".pdf")),
       surface_fig, width = 15, height = 5)

message("Stress-strain mortality surfaces written.")
message("Exploratory stress-strain mortality analysis complete.")
