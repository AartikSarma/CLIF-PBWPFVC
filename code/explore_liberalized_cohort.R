# =============================================================================
# Exploratory: does liberalizing the VT/PBW band rescue the strain signal?
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Exploratory, standalone analysis (NOT part of the 00 pipeline runner).
#
# The main cohort gates VT/PBW to 6-8, so VT is ~proportional to PBW and strain
# (VT/PFVC = VT/PBW x PBW/PFVC) is nearly a deterministic function of demographics
# -- which is why flexible demographic adjustment absorbed it (collinearity / VIF
# explosion). Liberalizing the VT/PBW band introduces dose variation that is NOT
# purely demographic, which should (a) relieve the collinearity (lower VIF) and
# (b) let the strain -> outcome signal survive flexible adjustment. The cost is
# confounding by indication (VT/PBW also varies with illness severity).
#
# This re-derives the index cohort under the restricted (6-8) and a liberalized
# (4-12) VT/PBW band from the pre-gate per-timepoint dataset, then for each:
#   - VIF of the size terms across the covariate-flexibility ladder, and
#   - the mortality evidence ratio for VT/PFVC (strain) vs VT/PBW,
# and fits a flexible strain dose-response (ns(VT/PFVC)) on mortality in the
# liberalized cohort. Mortality only (VFD competing-risks is not re-derived here).
#
# NOTE: index = the patient's FIRST qualifying timepoint (a simplification of the
# pipeline's two-tier rule; applied identically to both bands so the comparison is
# apples-to-apples). Severity is adjusted for but residual indication confounding
# remains -- the liberalized estimate is the strain dose-response, NOT a clean
# "PBW-mis-sizing" effect.
#
# NOTE (causal): the covariate ladder below includes height and spline-height as a
# COLLINEARITY diagnostic. Height is NOT a confounder (no path to mortality except
# via lung size), so adjusting for it is over-adjustment, not a preferred model; the
# VIF/ER "with height/spline" columns show what over-adjustment does, not the right
# specification. The dose-response model adjusts for age/height flexibly for the same
# diagnostic reason -- it is conservative (likely biased toward the null).
#
# Inputs : output/<site>/intermediate/analysis_all_eligible_timepoints.parquet (script 03)
#          output/<site>/intermediate/cohort_weights.parquet                    (script 01)
# Outputs: final/liberalized_cohort_<site>.csv     (VIF / evidence ratios, both cohorts)
#          final/liberalized_cohort_<site>.pdf      (VIF, ER, and strain dose-response)
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(splines)
library(patchwork)
library(marginaleffects)

source("utils/config.R")
site_name <- config$site_name

output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

RESTRICTED_BAND <- c(6, 8)
LIBERAL_BAND    <- c(4, 12)
SF_THRESHOLD    <- 315

pre <- read_parquet(file.path(output_dir, "analysis_all_eligible_timepoints.parquet"))
weights <- read_parquet(file.path(output_dir, "cohort_weights.parquet")) %>%
  select(hospitalization_id, weight_kg)

make_cohort <- function(band, label) {
  pre %>%
    filter(has_all_data, vtpbw >= band[1], vtpbw <= band[2], sf_ratio < SF_THRESHOLD) %>%
    group_by(hospitalization_id) %>%
    slice_min(recorded_dttm, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    left_join(weights, by = "hospitalization_id") %>%
    mutate(
      cohort = label,
      bmi = if_else(!is.na(weight_kg) & height_cm > 0, weight_kg / (height_cm / 100)^2, NA_real_),
      sex_category  = factor(sex_category,  levels = c("Male", "Female")),
      race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
      age10 = age_at_admission / 10, height10 = height_cm / 10, sf10 = sf_ratio / 10
    )
}
cohort_list <- list(
  "Restricted (6-8)"   = make_cohort(RESTRICTED_BAND, "Restricted (6-8)"),
  "Liberalized (4-12)" = make_cohort(LIBERAL_BAND,    "Liberalized (4-12)")
)

walk2(cohort_list, names(cohort_list), function(d, nm) {
  message(nm, ": ", nrow(d), " patients, ", sum(d$deceased), " deaths; VT/PBW range ",
          round(min(d$vtpbw), 1), "-", round(max(d$vtpbw), 1),
          " (SD ", round(sd(d$vtpbw), 2), ")")
})

# =============================================================================
# VIF and mortality evidence ratio across the covariate-flexibility ladder
# =============================================================================
covar_levels <- c("Linear (base)", "+ height", "+ spline age & height")
build_covars <- function(level) switch(level,
  "Linear (base)"         = "race_category + age10 + sex_category + sofa_total + sf10",
  "+ height"              = "race_category + age10 + sex_category + sofa_total + sf10 + height10",
  "+ spline age & height" = paste("race_category + ns(age_at_admission, 4) + sex_category +",
                                  "sofa_total + sf10 + ns(height_cm, 4)"))

vif_term <- function(d, target, others)
  1 / (1 - summary(lm(as.formula(paste(target, "~", others)), data = d))$r.squared)

mort_er <- function(d, spec_rhs, covars) {
  full <- AIC(glm(as.formula(paste("deceased ~", spec_rhs, "+", covars)), data = d, family = binomial))
  ref  <- AIC(glm(as.formula(paste("deceased ~ vtpbw +", covars)), data = d, family = binomial))
  exp(-0.5 * (full - ref))
}

results <- map_dfr(names(cohort_list), function(cn) {
  d <- cohort_list[[cn]]
  map_dfr(covar_levels, function(lv) {
    covs <- build_covars(lv)
    tibble(
      cohort = cn, level = lv,
      vif_VTPFVC  = vif_term(d, "vtpfvc", paste("vtpbw +", covs)),
      vif_PBWPFVC = vif_term(d, "pbwpfvc", paste("vtpbw +", covs)),
      er_VTPFVC        = mort_er(d, "vtpfvc", covs),
      er_VTPBW_PBWPFVC = mort_er(d, "vtpbw + pbwpfvc", covs)
    )
  })
}) %>%
  mutate(level = factor(level, levels = covar_levels),
         cohort = factor(cohort, levels = names(cohort_list)))

write_csv(results, file.path(final_dir, paste0("liberalized_cohort_", site_name, ".csv")))

message("VIF(VT/PFVC) and mortality ER(VT/PFVC vs VT/PBW) across the ladder:")
pwalk(results, function(cohort, level, vif_VTPFVC, er_VTPFVC, ...) {
  message("  [", cohort, "] ", level, ": VIF=", round(vif_VTPFVC, 1),
          ", ER=", signif(er_VTPFVC, 3))
})

# =============================================================================
# Flexible strain dose-response on mortality (liberalized cohort)
# =============================================================================
lib <- cohort_list[["Liberalized (4-12)"]]
dr <- glm(deceased ~ ns(vtpfvc, 4) + ns(age_at_admission, 4) + sex_category +
            race_category + sofa_total + sf10 + ns(height_cm, 4),
          data = lib, family = binomial)
vtpfvc_grid <- seq(quantile(lib$vtpfvc, 0.05), quantile(lib$vtpfvc, 0.95), length.out = 50)
dr_preds <- predictions(dr, newdata = datagrid(model = dr, vtpfvc = vtpfvc_grid)) %>%
  as_tibble()

# =============================================================================
# Figures
# =============================================================================
okabe <- c("Restricted (6-8)" = "#0072B2", "Liberalized (4-12)" = "#D55E00")

p_vif <- ggplot(results, aes(level, vif_VTPFVC, color = cohort, group = cohort)) +
  geom_line() + geom_point(size = 2) +
  scale_color_manual(values = okabe, name = NULL) +
  scale_y_log10() +
  labs(x = NULL, y = "VIF of VT/PFVC (log scale)",
       title = "Collinearity of strain with demographics") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_er <- ggplot(results, aes(level, er_VTPFVC, color = cohort, group = cohort)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_line() + geom_point(size = 2) +
  scale_color_manual(values = okabe, name = NULL) +
  scale_y_log10() +
  labs(x = NULL, y = "Evidence ratio: VT/PFVC vs VT/PBW (log)",
       title = "Does strain survive flexible adjustment?") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_dr <- ggplot(dr_preds, aes(vtpfvc, estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#D55E00", alpha = 0.15) +
  geom_line(color = "#D55E00") +
  geom_rug(data = lib, aes(x = vtpfvc), inherit.aes = FALSE, alpha = 0.1, sides = "b") +
  labs(x = "Strain estimate: VT/PFVC", y = "Predicted mortality",
       title = "Flexible strain dose-response (liberalized cohort)") +
  theme_minimal(base_size = 10)

fig <- (p_vif | p_er) / p_dr +
  plot_annotation(
    title = "Liberalizing VT/PBW: does it rescue the strain signal from collinearity?",
    subtitle = paste0(site_name,
      " — if liberalizing lowers VT/PFVC's VIF and keeps its evidence ratio above 1 ",
      "under spline adjustment, the strain dose-response is identifiable (confounded ",
      "by severity, not demographics)"))

ggsave(file.path(final_dir, paste0("liberalized_cohort_", site_name, ".pdf")),
       fig, width = 11, height = 9)

message("Liberalized-cohort outputs written.")
message("Exploratory liberalized-cohort analysis complete.")
