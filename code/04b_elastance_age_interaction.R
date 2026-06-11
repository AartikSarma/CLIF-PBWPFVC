# =============================================================================
# Script 04b: Normalized elastance x age interaction
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Motivation: normalized elastance (Ers x PBW, Ers x PFVC) is used as a surrogate
# for severity of lung injury, but the elastic recoil of lung tissue changes with
# age independently of injury -- so age confounds elastance as a severity marker.
# This script tests whether the elastance-mortality association differs by age by
# adding an elastance x age interaction to the mortality models from script 04,
# and uses marginaleffects to express the interaction on the interpretable
# probability scale (the elastance slope at representative ages).
#
# Inputs : output/<site>/intermediate/analysis_cross_sectional.parquet (script 03)
# Outputs: final/regression_ers_age_interaction_<site>.html  (model table)
#          final/ers_age_interaction_<site>.csv              (poolable summary)
#          final/ers_age_interaction_slopes_<site>.pdf       (slope-by-age plot)
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(gtsummary)
library(gt)
library(broom)
library(marginaleffects)

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
# age10 / sf10 are per-10-unit; ers_pbw / ers_pfvc are z-scaled (per SD) so the
# elastance main effect and its age interaction are on a common, comparable scale
# across the two normalizations.
zscore <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

model_data <- cross_sectional %>%
  mutate(
    sex_category  = factor(sex_category,  levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
    age10      = age_at_admission / 10,
    sf10       = sf_ratio / 10,
    ers_pbw_z  = zscore(ers_pbw),
    ers_pfvc_z = zscore(ers_pfvc)
  )

# Severity-of-illness and demographic adjustments shared by every model. Age
# enters through the interaction term, so it is NOT repeated here. BMI is included
# because the elastance-normalized exposures are driving-pressure derivatives
# (matching the original paper and the script 04 ers models). VT/PBW is retained
# as a severity adjustment, as in the script 04 ers mortality models.
adjustment <- "vtpbw + race_category + sex_category + sofa_total + sf10 + bmi"

ers_specs <- c(ers_pbw_z = "Ers x PBW (per SD)", ers_pfvc_z = "Ers x PFVC (per SD)")

# Representative ages for the interpretable slope/prediction summaries.
AGE_POINTS <- c(40, 60, 80)

# =============================================================================
# Fit interaction models and their nested no-interaction counterparts
# =============================================================================
has_mortality_variation <- length(unique(na.omit(model_data$deceased))) > 1
if (!has_mortality_variation) {
  stop("Cannot fit elastance x age interaction: outcome `deceased` has no ",
       "variation (all = ", unique(na.omit(model_data$deceased)), ").")
}

interaction_models    <- list()
no_interaction_models <- list()
for (z in names(ers_specs)) {
  f_int  <- as.formula(paste0("deceased ~ ", z, " * age10 + ", adjustment))
  f_main <- as.formula(paste0("deceased ~ ", z, " + age10 + ", adjustment))
  interaction_models[[z]]    <- glm(f_int,  data = model_data, family = binomial)
  no_interaction_models[[z]] <- glm(f_main, data = model_data, family = binomial)
}

# =============================================================================
# Likelihood-ratio test: does the age interaction improve fit?
# =============================================================================
# The no-interaction model is nested in the interaction model (1 extra df), so a
# LRT directly tests the interaction. AIC is reported alongside for comparability.
lrt_summary <- imap_dfr(interaction_models, function(m_int, z) {
  m_main <- no_interaction_models[[z]]
  lrt <- anova(m_main, m_int, test = "LRT")
  tibble(
    site            = site_name,
    exposure        = ers_specs[[z]],
    term            = paste0(z, ":age10"),
    n_obs           = stats::nobs(m_int),
    aic_no_interact = AIC(m_main),
    aic_interact    = AIC(m_int),
    delta_aic       = AIC(m_int) - AIC(m_main),
    lrt_chisq       = lrt$Deviance[2],
    lrt_df          = lrt$Df[2],
    lrt_p           = lrt$`Pr(>Chi)`[2]
  )
})

# Interaction odds ratio (per SD elastance, per 10 yr of age) from each model.
interaction_or <- imap_dfr(interaction_models, function(m, z) {
  broom::tidy(m, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == paste0(z, ":age10")) %>%
    transmute(site = site_name, exposure = ers_specs[[z]], term,
              or = estimate, conf_low = conf.low, conf_high = conf.high,
              p_value = p.value)
})

interaction_out <- lrt_summary %>%
  left_join(interaction_or %>% select(exposure, or, conf_low, conf_high, p_value),
            by = "exposure")

write_csv(interaction_out,
          file.path(final_dir, paste0("ers_age_interaction_", site_name, ".csv")))

message("Elastance x age interaction (LRT vs no-interaction model):")
pwalk(interaction_out, function(exposure, or, lrt_p, ...) {
  message("  ", exposure, ": interaction OR = ", round(or, 3),
          ", LRT p = ", signif(lrt_p, 3))
})

# =============================================================================
# Regression table (interaction models, odds ratios)
# =============================================================================
interaction_tables <- imap(interaction_models, ~ {
  tbl_regression(
    .x, exponentiate = TRUE,
    label = list(
      age10 ~ "Age (per 10 yr)",
      sf10  ~ "SF ratio (per 10)",
      vtpbw ~ "VT/PBW"
    )
  ) %>% bold_p()
})
tbl_merge(interaction_tables, tab_spanner = unname(ers_specs)) %>%
  as_gt() %>%
  gt::gtsave(file.path(final_dir,
                       paste0("regression_ers_age_interaction_", site_name, ".html")))
message("Interaction regression table written")

# =============================================================================
# Interpretable interaction: elastance slope on the probability scale, by age
# =============================================================================
# marginaleffects::slopes gives the change in predicted mortality probability per
# 1 SD of elastance, evaluated at each representative age with all other
# covariates held at their observed values (average slope). A non-flat slope-vs-
# age profile is the interaction expressed on the clinically meaningful scale.
slopes_by_age <- imap_dfr(interaction_models, function(m, z) {
  slopes(
    m,
    variables = z,
    newdata = datagrid(age10 = AGE_POINTS / 10, grid_type = "counterfactual"),
    by = "age10"
  ) %>%
    as_tibble() %>%
    transmute(
      exposure  = ers_specs[[z]],
      age_years = age10 * 10,
      slope     = estimate,
      conf_low  = conf.low,
      conf_high = conf.high,
      p_value   = p.value
    )
})

write_csv(slopes_by_age,
          file.path(final_dir, paste0("ers_age_interaction_slopes_", site_name, ".csv")))

# Okabe-Ito palette for the discrete elastance metric (per project convention).
okabe_ito <- c("#E69F00", "#0072B2")
slopes_plot <- ggplot(slopes_by_age,
                      aes(x = age_years, y = slope, color = exposure)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line() +
  geom_pointrange(aes(ymin = conf_low, ymax = conf_high)) +
  scale_color_manual(values = okabe_ito, name = "Normalized elastance") +
  labs(
    title = "Elastance effect on mortality by age",
    subtitle = paste0(site_name,
      " — change in predicted mortality probability per 1 SD elastance (95% CI)"),
    x = "Age (years)",
    y = "Marginal effect on mortality probability (per SD)"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(final_dir, paste0("ers_age_interaction_slopes_", site_name, ".pdf")),
       slopes_plot, width = 8, height = 5)

message("Slope-by-age plot written")
message("Script 04b complete.")
