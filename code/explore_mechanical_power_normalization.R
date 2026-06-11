# =============================================================================
# Exploratory: mechanical power normalized to PBW vs PFVC
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Exploratory, standalone analysis (NOT part of the 00 pipeline runner). Reads the
# script 03 cross-sectional dataset and refits its own models.
#
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
library(survival)
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

# Representative ages for the interpretable interaction slopes (see below).
AGE_POINTS <- c(40, 60, 80)

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

# =============================================================================
# Mechanical power x age interaction
# =============================================================================
# As with elastance (explore_elastance_age_interaction.R), lung and chest-wall mechanics change with age
# independently of injury, so age can modify the mechanical-power / mortality
# association. Re-fit each model with an MP x age interaction, test it by
# likelihood-ratio test against the no-interaction model, and use marginaleffects
# to express the result on the probability scale (the MP slope at each age).
# Age enters through the interaction, so it is dropped from the adjustment set.
adjustment_no_age <- "vtpbw + race_category + sex_category + sofa_total + sf10 + bmi"

interaction_models    <- list()
no_interaction_models <- list()
for (v in names(mp_specs)) {
  f_int  <- as.formula(paste0("deceased ~ ", v, " * age10 + ", adjustment_no_age))
  f_main <- as.formula(paste0("deceased ~ ", v, " + age10 + ", adjustment_no_age))
  interaction_models[[v]]    <- glm(f_int,  data = model_data, family = binomial)
  no_interaction_models[[v]] <- glm(f_main, data = model_data, family = binomial)
}

# Likelihood-ratio test (the no-interaction model is nested, 1 extra df) plus the
# interaction odds ratio (per SD MP, per 10 yr of age).
mp_age_interaction <- imap_dfr(interaction_models, function(m_int, v) {
  m_main <- no_interaction_models[[v]]
  lrt <- anova(m_main, m_int, test = "LRT")
  or_row <- broom::tidy(m_int, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == paste0(v, ":age10"))
  tibble(
    site            = site_name,
    exposure        = mp_specs[[v]],
    term            = paste0(v, ":age10"),
    n_obs           = stats::nobs(m_int),
    interaction_or  = or_row$estimate,
    conf_low        = or_row$conf.low,
    conf_high       = or_row$conf.high,
    p_value         = or_row$p.value,
    aic_no_interact = AIC(m_main),
    aic_interact    = AIC(m_int),
    delta_aic       = AIC(m_int) - AIC(m_main),
    lrt_chisq       = lrt$Deviance[2],
    lrt_df          = lrt$Df[2],
    lrt_p           = lrt$`Pr(>Chi)`[2]
  )
})

write_csv(mp_age_interaction,
          file.path(final_dir, paste0("mp_age_interaction_", site_name, ".csv")))

message("Mechanical power x age interaction (LRT vs no-interaction model):")
pwalk(mp_age_interaction, function(exposure, interaction_or, lrt_p, ...) {
  message("  ", exposure, ": interaction OR = ", round(interaction_or, 3),
          ", LRT p = ", signif(lrt_p, 3))
})

# Regression table for the interaction models.
interaction_tables <- imap(interaction_models, ~ {
  tbl_regression(
    .x, exponentiate = TRUE,
    label = c(
      setNames(list(mp_specs[[.y]]), .y),
      list(age10 = "Age (per 10 yr)", sf10 = "SF ratio (per 10)", vtpbw = "VT/PBW")
    )
  ) %>% bold_p()
})
tbl_merge(interaction_tables, tab_spanner = unname(mp_specs)) %>%
  as_gt() %>%
  gt::gtsave(file.path(final_dir,
                       paste0("regression_mp_age_interaction_", site_name, ".html")))

# Interpretable interaction: MP slope on the probability scale, by age. A non-flat
# slope-vs-age profile is the interaction expressed in clinical terms.
slopes_by_age <- imap_dfr(interaction_models, function(m, v) {
  slopes(
    m,
    variables = v,
    newdata = datagrid(age10 = AGE_POINTS / 10, grid_type = "counterfactual"),
    by = "age10"
  ) %>%
    as_tibble() %>%
    transmute(
      exposure  = mp_specs[[v]],
      age_years = age10 * 10,
      slope     = estimate,
      conf_low  = conf.low,
      conf_high = conf.high,
      p_value   = p.value
    )
})

write_csv(slopes_by_age,
          file.path(final_dir, paste0("mp_age_interaction_slopes_", site_name, ".csv")))

okabe_ito <- c("#E69F00", "#0072B2")
slopes_plot <- ggplot(slopes_by_age,
                      aes(x = age_years, y = slope, color = exposure)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line() +
  geom_pointrange(aes(ymin = conf_low, ymax = conf_high)) +
  scale_color_manual(values = okabe_ito, name = "Normalized mechanical power") +
  labs(
    title = "Mechanical power effect on mortality by age",
    subtitle = paste0(site_name,
      " — change in predicted mortality probability per 1 SD MP (95% CI)"),
    x = "Age (years)",
    y = "Marginal effect on mortality probability (per SD)"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(final_dir, paste0("mp_age_interaction_slopes_", site_name, ".pdf")),
       slopes_plot, width = 8, height = 5)

message("Mechanical power x age interaction outputs written")

# =============================================================================
# Second outcome: 28-day VFDs (competing risks)
# =============================================================================
# Repeat the MP/PBW-vs-MP/PFVC comparison and the MP x age interaction for the
# competing-risks VFD outcome (extubation = event, death = competing risk; Yehya &
# Harhay 2019), using Fine-Gray subdistribution models. Effects are subdistribution
# hazard ratios for liberation (SHR > 1 = faster liberation).
vfd_status_levels <- c(0, 1, 2)
vfd_status_labels <- c("censored", "extubation", "death")

fit_vfd_finegray <- function(rhs) {
  rhs_vars <- unique(trimws(unlist(strsplit(rhs, "\\+|\\*"))))
  df <- model_data %>%
    mutate(vfd_status_f = factor(vfd_status, levels = vfd_status_levels,
                                 labels = vfd_status_labels)) %>%
    select(vfd_time, vfd_status_f, all_of(rhs_vars)) %>%
    filter(!is.na(vfd_time), vfd_time > 0, !is.na(vfd_status_f))
  fg <- finegray(Surv(vfd_time, vfd_status_f) ~ ., data = df, etype = "extubation")
  coxph(as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ ", rhs)),
        weights = fgwt, data = fg)
}

# --- Main comparison: MP/PBW vs MP/PFVC (SHR per SD), referenced to MP/PBW ------
mp_vfd_models <- map(names(mp_specs), ~ fit_vfd_finegray(paste(.x, "+", adjustment)))
names(mp_vfd_models) <- names(mp_specs)
aic_pbw_vfd <- AIC(mp_vfd_models[["mp_pbw_z"]])

mp_vfd_summary <- imap_dfr(mp_vfd_models, function(m, v) {
  shr_row <- broom::tidy(m, conf.int = TRUE, exponentiate = TRUE) %>% filter(term == v)
  tibble(
    site = site_name, exposure = mp_specs[[v]], n_obs = m$n,
    shr_per_sd = shr_row$estimate, conf_low = shr_row$conf.low,
    conf_high = shr_row$conf.high, p_value = shr_row$p.value,
    aic = AIC(m), delta_aic = AIC(m) - aic_pbw_vfd,
    evidence_ratio = exp(-0.5 * (AIC(m) - aic_pbw_vfd))
  )
})
write_csv(mp_vfd_summary,
          file.path(final_dir, paste0("mp_normalization_vfd_", site_name, ".csv")))

mp_vfd_tables <- imap(mp_vfd_models, ~ {
  tbl_regression(.x, exponentiate = TRUE,
    label = c(setNames(list(mp_specs[[.y]]), .y),
              list(age10 = "Age (per 10 yr)", sf10 = "SF ratio (per 10)",
                   vtpbw = "VT/PBW"))) %>% bold_p()
})
tbl_merge(mp_vfd_tables, tab_spanner = unname(mp_specs)) %>%
  as_gt() %>%
  gt::gtsave(file.path(final_dir, paste0("regression_mp_normalization_vfd_", site_name, ".html")))

message("Mechanical power normalization on 28-day VFDs (SHR per SD; ",
        "evidence ratio vs MP/PBW):")
pwalk(mp_vfd_summary, function(exposure, shr_per_sd, evidence_ratio, ...) {
  message("  ", exposure, ": SHR = ", round(shr_per_sd, 3),
          ", evidence ratio = ", signif(evidence_ratio, 3))
})

# --- MP x age interaction on the VFD outcome (SHR + LRT) ------------------------
mp_vfd_int_models  <- list()
mp_vfd_main_models <- list()
for (v in names(mp_specs)) {
  mp_vfd_int_models[[v]]  <- fit_vfd_finegray(paste0(v, " * age10 + ", adjustment_no_age))
  mp_vfd_main_models[[v]] <- fit_vfd_finegray(paste0(v, " + age10 + ", adjustment_no_age))
}

mp_vfd_interaction <- imap_dfr(mp_vfd_int_models, function(m_int, v) {
  lrt <- anova(mp_vfd_main_models[[v]], m_int)
  shr <- broom::tidy(m_int, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == paste0(v, ":age10"))
  tibble(
    site = site_name, exposure = mp_specs[[v]], term = paste0(v, ":age10"),
    n_obs = m_int$n,
    interaction_shr = shr$estimate, conf_low = shr$conf.low, conf_high = shr$conf.high,
    p_value = shr$p.value,
    lrt_chisq = lrt$Chisq[2], lrt_df = lrt$Df[2], lrt_p = lrt$`Pr(>|Chi|)`[2]
  )
})
write_csv(mp_vfd_interaction,
          file.path(final_dir, paste0("mp_age_interaction_vfd_", site_name, ".csv")))

message("Mechanical power x age interaction on 28-day VFDs (SHR; LRT):")
pwalk(mp_vfd_interaction, function(exposure, interaction_shr, lrt_p, ...) {
  message("  ", exposure, ": interaction SHR = ", round(interaction_shr, 3),
          ", LRT p = ", signif(lrt_p, 3))
})

message("VFD competing-risks outputs written")
message("Exploratory mechanical-power normalization complete.")
