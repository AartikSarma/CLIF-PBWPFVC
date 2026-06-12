# =============================================================================
# Script 04: Analysis
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(gtsummary)
library(gt)
library(survival)
library(survminer)
library(splines)
library(patchwork)
library(broom)
library(algorithmDiagnostics)
library(EValue)

source("utils/config.R")
source("utils/consort_diagram.R")
site_name <- config$site_name

output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# Load data
# =============================================================================

cross_sectional <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"))
analysis_all <- read_parquet(file.path(output_dir, "analysis_all_timepoints.parquet"))

message("Cross-sectional: ", nrow(cross_sectional), " observations")
message("All timepoints: ", nrow(analysis_all), " observations")

# Stratify Table 1 by PBW/PFVC tercile (the analysis now centers on the
# PBW:PFVC discrepancy rather than delivered VT/PFVC).
cross_sectional <- cross_sectional %>%
  mutate(
    pbwpfvc_tercile = ntile(pbwpfvc, 3),
    pbwpfvc_tercile = factor(pbwpfvc_tercile, labels = c("T1 (Low)", "T2 (Mid)", "T3 (High)"))
  )

# Add age_group for downstream use
cross_sectional <- cross_sectional %>%
  mutate(age_group = cut(age_at_admission,
                          breaks = c(18, 40, 60, 80, Inf),
                          labels = c("18-39", "40-59", "60-79", "80+"),
                          right = FALSE))

# Reference categories: Male and White. Set in script 03, but re-applied here so
# the reference is explicit and robust to the parquet round-trip; all regressions
# below report effects relative to male / white patients.
cross_sectional <- cross_sectional %>%
  mutate(
    sex_category  = factor(sex_category,  levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
    # Per-10-unit covariates so the adjusted age and SF-ratio coefficients are
    # reported per 10 years / 10 SF units (Table 1 keeps the raw scales).
    age10 = age_at_admission / 10,
    sf10  = sf_ratio / 10
  )

# =============================================================================
# 4a. Table 1
# =============================================================================

table1_data <- cross_sectional %>%
  select(pbwpfvc_tercile, age_at_admission, sex_category, race_category,
         height_cm, pbw, pfvc, sofa_total, sf_ratio, deceased,
         vtpbw, vtpfvc, pbwpfvc, crs, ers, vfd_28) %>%
  mutate(deceased = factor(deceased, levels = c(0, 1), labels = c("Alive", "Deceased")))

table1 <- table1_data %>%
  tbl_summary(
    by = pbwpfvc_tercile,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      age_at_admission ~ "Age (years)",
      sex_category ~ "Sex",
      race_category ~ "Race",
      height_cm ~ "Height (cm)",
      pbw ~ "PBW (kg)",
      pfvc ~ "PFVC (L)",
      sofa_total ~ "SOFA Total",
      sf_ratio ~ "SF Ratio",
      deceased ~ "Mortality",
      vtpbw ~ "VT/PBW (mL/kg)",
      vtpfvc ~ "VT/PFVC (%)",
      pbwpfvc ~ "PBW/PFVC",
      crs ~ "Compliance (mL/cmH2O)",
      ers ~ "Elastance (cmH2O/L)",
      vfd_28 ~ "28-day VFDs"
    )
  ) %>%
  add_p() %>%
  add_overall() %>%
  # The stratifying columns are terciles of PBW/PFVC; without a spanning header
  # the bare "T1 (Low) / T2 (Mid) / T3 (High)" labels don't say tercile of what.
  modify_spanning_header(
    all_stat_cols(stat_0 = FALSE) ~ "**Predicted body weight / predicted FVC (PBW/PFVC), tercile**"
  ) %>%
  modify_header(
    label  ~ "**Characteristic**",
    stat_0 ~ "**Overall**, N = {N}"
  )

message("Table 1 generated")

table1 %>%
  as_gt() %>%
  gt::gtsave(file.path(final_dir, paste0("table1_", site_name, ".pdf")))

table1 %>%
  as_gt() %>%
  gt::gtsave(file.path(final_dir, paste0("table1_", site_name, ".html")))

# =============================================================================
# 4b. Model specifications
# =============================================================================

# Common covariates for all models
covariates <- "race_category + age10 + sex_category + sofa_total + sf10"

# Readable labels for the per-10-unit covariates in the rendered regression
# tables (age10 / sf10 appear in every model's covariate set).
covar_labels <- list(age10 ~ "Age (per 10 yr)", sf10 ~ "SF ratio (per 10)")

# Per the original paper, any model whose outcome OR exposure is derived from
# driving pressure (static DP, elastance, compliance, and the elastance-normalized
# Ers x PBW / Ers x PFVC) is additionally adjusted for BMI. Models without a
# driving-pressure component use the standard covariate set.
DP_DERIVED <- c("dp", "ers", "crs", "ers_pbw", "ers_pfvc",
                "mechanical_power", "mp_pbw", "mp_pfvc")
uses_dp <- function(...) {
  vars <- trimws(unlist(strsplit(paste(c(...), collapse = " + "), "\\+")))
  any(vars %in% DP_DERIVED)
}
model_covariates <- function(...) {
  if (uses_dp(...)) paste(covariates, "+ bmi") else covariates
}

# Define all exposure specifications
exposure_specs <- list(
  vtpfvc        = "vtpfvc",
  vtpbw         = "vtpbw",
  vtpfvc_vtpbw  = "vtpfvc + vtpbw",
  vtpbw_pfvc    = "vtpbw + pfvc",
  vtpbw_pbwpfvc = "vtpbw + pbwpfvc"
)

exposure_labels <- c(
  vtpfvc        = "VT/PFVC",
  vtpbw         = "VT/PBW",
  vtpfvc_vtpbw  = "VT/PFVC + VT/PBW",
  vtpbw_pfvc    = "VT/PBW + PFVC",
  vtpbw_pbwpfvc = "VT/PBW + PBW/PFVC"
)

# =============================================================================
# 4c. Logistic regression — mortality
# =============================================================================

has_mortality_variation <- length(unique(na.omit(cross_sectional$deceased))) > 1

if (has_mortality_variation) {
  mortality_models <- map(exposure_specs, ~ {
    formula_str <- paste("deceased ~", .x, "+", model_covariates(.x))
    glm(as.formula(formula_str), data = cross_sectional, family = binomial)
  })

  mortality_tables <- map(mortality_models, ~ {
    tbl_regression(.x, exponentiate = TRUE, label = covar_labels) %>% bold_p()
  })

  message("Logistic regression (mortality) — AIC:")
  iwalk(mortality_models, ~ message("  ", exposure_labels[.y], ": ", round(AIC(.x), 1)))

  tbl_merge(mortality_tables, tab_spanner = exposure_labels) %>%
    as_gt() %>%
    gt::gtsave(file.path(final_dir, paste0("regression_mortality_", site_name, ".html")))
} else {
  mortality_models <- NULL
  message("Skipping mortality regression: no variation in outcome (all deceased = ",
          unique(na.omit(cross_sectional$deceased)), ")")
}

# =============================================================================
# 4d. Linear regression — continuous outcomes (elastance, compliance, VFD-28, DP)
# =============================================================================

# NOTE: 28-day VFDs are NOT modeled here as a continuous outcome. Per Yehya &
# Harhay (AJRCCM 2019), VFDs are analyzed as a competing-risks outcome (extubation
# vs death) in section 4d2 below; mortality (section 4c) is the other component.
continuous_outcomes <- list(
  ers      = list(var = "ers",              label = "Elastance"),
  crs      = list(var = "crs",              label = "Compliance"),
  dp       = list(var = "dp",               label = "Static DP"),
  ers_pbw  = list(var = "ers_pbw",          label = "Ers x PBW"),
  ers_pfvc = list(var = "ers_pfvc",         label = "Ers x PFVC"),
  mp       = list(var = "mechanical_power", label = "Mechanical power"),
  mp_pbw   = list(var = "mp_pbw",           label = "MP / PBW"),
  mp_pfvc  = list(var = "mp_pfvc",          label = "MP / PFVC")
)

continuous_models <- list()
continuous_tables <- list()

for (outcome_name in names(continuous_outcomes)) {
  outcome_var <- continuous_outcomes[[outcome_name]]$var
  outcome_label <- continuous_outcomes[[outcome_name]]$label

  covars_used <- model_covariates(outcome_var)
  models_for_outcome <- map(exposure_specs, ~ {
    formula_str <- paste(outcome_var, "~", .x, "+", covars_used)
    lm(as.formula(formula_str), data = cross_sectional)
  })

  tables_for_outcome <- map(models_for_outcome, ~ {
    tbl_regression(.x, label = covar_labels) %>% bold_p()
  })

  continuous_models[[outcome_name]] <- models_for_outcome
  continuous_tables[[outcome_name]] <- tables_for_outcome

  message(outcome_label, " models — AIC:")
  iwalk(models_for_outcome, ~ message("  ", exposure_labels[.y], ": ", round(AIC(.x), 1)))

  tbl_merge(tables_for_outcome, tab_spanner = exposure_labels) %>%
    as_gt() %>%
    gt::gtsave(file.path(final_dir, paste0("regression_", outcome_name, "_", site_name, ".html")))
}

# =============================================================================
# 4d2. 28-day VFDs — competing-risks (Fine-Gray) regression
# =============================================================================
# Per Yehya & Harhay (AJRCCM 2019), VFDs are analyzed as a competing-risks
# outcome: event of interest = extubation (vfd_status 1), competing risk = death
# within 28 days (vfd_status 2), censored at day 28 if still ventilated
# (vfd_status 0). Each exposure specification is fit with a Fine-Gray
# subdistribution hazard model (survival::finegray weights + coxph), so effects
# are subdistribution hazard ratios for liberation (SHR > 1 = faster liberation).
# The mortality component is modeled in section 4c.
vfd_cr_covariates <- covariates  # VFD + these exposures are never DP-derived (no BMI)

fit_vfd_finegray <- function(exposure_spec) {
  model_rhs <- paste(exposure_spec, "+", vfd_cr_covariates)
  rhs_vars  <- unique(trimws(unlist(strsplit(model_rhs, "\\+"))))
  df <- cross_sectional %>%
    mutate(vfd_status_f = factor(vfd_status, levels = c(0, 1, 2),
                                 labels = c("censored", "extubation", "death"))) %>%
    select(vfd_time, vfd_status_f, all_of(rhs_vars)) %>%
    filter(!is.na(vfd_time), vfd_time > 0, !is.na(vfd_status_f))
  fg <- survival::finegray(survival::Surv(vfd_time, vfd_status_f) ~ ., data = df,
                           etype = "extubation")
  survival::coxph(
    as.formula(paste0("survival::Surv(fgstart, fgstop, fgstatus) ~ ", model_rhs)),
    weights = fgwt, data = fg
  )
}

vfd_cr_models <- map(exposure_specs, fit_vfd_finegray)

vfd_cr_tables <- map(vfd_cr_models, ~ {
  tbl_regression(.x, exponentiate = TRUE, label = covar_labels) %>% bold_p()
})

message("28-day VFDs (Fine-Gray, extubation SHR) — AIC:")
iwalk(vfd_cr_models, ~ message("  ", exposure_labels[.y], ": ", round(AIC(.x), 1)))

tbl_merge(vfd_cr_tables, tab_spanner = exposure_labels) %>%
  as_gt() %>%
  gt::gtsave(file.path(final_dir, paste0("regression_vfd28_", site_name, ".html")))

# =============================================================================
# 4e. AIC comparison across all models and outcomes
# =============================================================================

aic_results <- list()

# Mortality
if (!is.null(mortality_models)) {
  aic_results[["Mortality"]] <- tibble(
    exposure = exposure_labels,
    AIC = map_dbl(mortality_models, AIC),
    is_reference = exposure == "VT/PBW"
  )
}

# Continuous outcomes. Exclude the normalized-mechanics outcomes (Ers x PBW/PFVC,
# MP/PBW/PFVC) from the AIC comparison/heatmap: the outcome shares a predicted-size
# variable with the size-containing exposures, so their evidence ratios are inflated
# by the shared denominator rather than a dosing-outcome relationship. Their
# standalone regression tables and long-format results are still produced.
AIC_EXCLUDE <- c("ers_pbw", "ers_pfvc", "mp_pbw", "mp_pfvc")
for (outcome_name in setdiff(names(continuous_outcomes), AIC_EXCLUDE)) {
  aic_results[[continuous_outcomes[[outcome_name]]$label]] <- tibble(
    exposure = exposure_labels,
    AIC = map_dbl(continuous_models[[outcome_name]], AIC),
    is_reference = exposure == "VT/PBW"
  )
}

# 28-day VFDs (competing-risks Fine-Gray models). AICs are comparable within the
# outcome and referenced to the VT/PBW model, as elsewhere.
aic_results[["28-day VFDs"]] <- tibble(
  exposure = exposure_labels,
  AIC = map_dbl(vfd_cr_models, AIC),
  is_reference = exposure == "VT/PBW"
)

# Evidence ratios are all referenced to the VT/PBW-alone model WITHIN each
# outcome: ER = exp(-0.5 * (AIC_model - AIC_VT/PBW)). ER > 1 means more support
# than VT/PBW alone, ER = 1 for VT/PBW itself. Truncated to [0.001, 1000].
ER_FLOOR <- 0.001
ER_CEIL  <- 1000
aic_all <- bind_rows(aic_results, .id = "outcome") %>%
  group_by(outcome) %>%
  mutate(
    aic_ref = AIC[is_reference][1],
    delta_AIC = AIC - aic_ref,
    evidence_ratio = exp(-0.5 * delta_AIC),
    evidence_ratio_trunc = pmin(pmax(evidence_ratio, ER_FLOOR), ER_CEIL),
    er_label = case_when(
      evidence_ratio > ER_CEIL  ~ ">1000",
      evidence_ratio < ER_FLOOR ~ "<0.001",
      TRUE                      ~ formatC(evidence_ratio, format = "g", digits = 2)
    )
  ) %>%
  ungroup() %>%
  # Column order: clinical outcomes first, then compliance/elastance, then static
  # driving pressure and mechanical power.
  mutate(outcome = factor(outcome, levels = intersect(
    c("Mortality", "28-day VFDs", "Compliance", "Elastance",
      "Static DP", "Mechanical power"),
    unique(outcome))))

# Row order: exposure specs by the strongest evidence ratio they reach in any
# outcome (best at the top). Raw AIC is not comparable across outcomes, so the
# AIC-based evidence ratio (truncated) is used; ties are broken by the summed
# log10 ratio.
exposure_order <- aic_all %>%
  group_by(exposure) %>%
  summarise(max_er = max(evidence_ratio_trunc, na.rm = TRUE),
            sum_lr = sum(log10(evidence_ratio_trunc), na.rm = TRUE),
            .groups = "drop") %>%
  arrange(max_er, sum_lr) %>%
  pull(exposure)
aic_all <- aic_all %>% mutate(exposure = factor(exposure, levels = exposure_order))

message("AIC comparison (evidence ratios vs VT/PBW-alone within each outcome):")
print(aic_all)

write_csv(aic_all, file.path(final_dir, paste0("aic_comparison_all_", site_name, ".csv")))

# Evidence ratio heatmap: divergent log10 colour scale, white = 1 (no difference
# from VT/PBW), blue = less support, red = more support; truncated to [0.001, 1000].
er_heatmap <- ggplot(aic_all,
                     aes(x = outcome, y = exposure,
                         fill = log10(evidence_ratio_trunc))) +
  geom_tile(color = "grey80", linewidth = 0.5) +
  geom_text(aes(label = er_label), size = 3.5) +
  scale_fill_gradient2(
    name = "Evidence ratio\n(vs VT/PBW)",
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(log10(ER_FLOOR), log10(ER_CEIL)),
    breaks = -3:3, labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000")
  ) +
  labs(
    title = "Evidence ratios across models and outcomes",
    subtitle = "Each cell vs the VT/PBW-alone model within that outcome; truncated to [0.001, 1000]",
    x = "Outcome",
    y = "Exposure specification"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(final_dir, paste0("evidence_ratio_heatmap_all_", site_name, ".pdf")),
       er_heatmap, width = 10, height = 8)

# =============================================================================
# 4f. Survival analysis
# =============================================================================

# Event = all-cause death within 60 days, in- or out-of-hospital, derived in
# script 03 (mortality_event_60) from the patient-level death_dttm. Using this
# instead of the in-hospital-only `deceased` flag stops survivors from being
# censored at hospital discharge and counts post-discharge deaths as events.
surv_data <- cross_sectional %>%
  filter(!is.na(surv_time), surv_time > 0, !is.na(mortality_event_60)) %>%
  mutate(event = as.integer(mortality_event_60))

n_deaths <- sum(surv_data$event, na.rm = TRUE)
message("Survival analysis: ", nrow(surv_data), " patients, ", n_deaths,
        " deaths within 60 days (in- and out-of-hospital), ",
        nrow(surv_data) - n_deaths, " censored at day 60")

if (n_deaths > 0 && length(unique(surv_data$event)) > 1) {
  cox_model <- coxph(
    Surv(surv_time, event) ~ pbwpfvc + vtpbw + age10 +
      sex_category + race_category + sf10 + sofa_total,
    data = surv_data
  )

  # Companion model with PFVC as the scaling exposure (mirrors the VT/PBW + PFVC
  # mortality model), so the survival analysis carries both headline exposures.
  cox_model_pfvc <- coxph(
    Surv(surv_time, event) ~ pfvc + vtpbw + age10 +
      sex_category + race_category + sf10 + sofa_total,
    data = surv_data
  )

  message("Cox model (PBW/PFVC):")
  print(summary(cox_model))
  message("Cox model (PFVC):")
  print(summary(cox_model_pfvc))

  # Fit on the labelled pbwpfvc_tercile factor (not ntile()) so the strata carry
  # informative names instead of "ntile(pbwpfvc, 3)=1".
  km_fit <- survfit(Surv(surv_time, event) ~ pbwpfvc_tercile, data = surv_data)

  km_plot <- ggsurvplot(
    km_fit,
    data = surv_data,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    palette = c("#CC5555", "#266cae", "#d47e0e"),
    xlab = "Days from Admission",
    ylab = "Survival Probability",
    title = "Kaplan-Meier Survival by PBW/PFVC Tercile",
    legend.title = "PBW/PFVC tercile (lower = PBW closer to PFVC)",
    legend.labs = c("Lowest tercile (T1)", "Middle tercile (T2)", "Highest tercile (T3)"),
    ggtheme = theme_minimal(),
    # Clean risk table: drop the background grid behind the at-risk counts.
    tables.theme = theme_cleantable()
  )

  pdf(file.path(final_dir, paste0("km_curves_", site_name, ".pdf")),
      width = 10, height = 8)
  print(km_plot)
  dev.off()

  message("KM curves saved")

  sink(file.path(final_dir, paste0("cox_model_summary_", site_name, ".txt")))
  cat("=== Cox model: PBW/PFVC + VT/PBW ===\n")
  print(summary(cox_model))
  cat("\n=== Cox model: VT/PBW + PFVC ===\n")
  print(summary(cox_model_pfvc))
  sink()
} else {
  message("Skipping survival analysis: no mortality events in data")
}

# =============================================================================
# 4f2. Demographic-bias models (outcome ~ demographics only)
# =============================================================================
# How each dosing / driving-pressure / elastance metric varies by demographics
# (the algorithmic-bias view). Predictors are scaled per 10 units (age, height,
# SF ratio). Outcomes for the dosing and elastance-normalized metrics are
# z-scored WITHIN SITE (so standardized betas are comparable across cohorts in SD
# units, and no row-level data is needed to pool). Static DP is on the raw cmH2O
# scale and additionally adjusts for BMI; Mortality is logistic (OR). Reference
# categories are Male and White.

demo_data <- cross_sectional %>%
  mutate(age10 = age_at_admission / 10,
         height10 = height_cm / 10,
         sf10 = sf_ratio / 10)

demo_covars     <- "age10 + sex_category + race_category + height10 + sf10 + sofa_total"
demo_covars_bmi <- paste(demo_covars, "+ bmi")

# Each entry: fitted model + metadata for the unified long table.
demo_models <- list()

# z-scored linear outcomes. The elastance-normalized outcomes (Ers x PBW / PFVC)
# are driving-pressure derivatives, so they are adjusted for BMI (matching the
# original paper); VT/PBW and VT/PFVC are not.
demo_z_outcomes <- c(vtpbw = "VT/PBW", vtpfvc = "VT/PFVC (%)",
                     ers_pbw = "Ers x PBW", ers_pfvc = "Ers x PFVC")
for (v in names(demo_z_outcomes)) {
  cov_v <- if (v %in% c("ers_pbw", "ers_pfvc")) demo_covars_bmi else demo_covars
  fstr <- paste0("scale(", v, ") ~ ", cov_v)
  demo_models[[demo_z_outcomes[[v]]]] <- list(
    model = lm(as.formula(fstr), data = demo_data),
    type = "Beta", family = "linear", formula = fstr
  )
}

# Mortality logistic (no BMI)
if (has_mortality_variation) {
  fstr <- paste("deceased ~", demo_covars)
  demo_models[["Mortality"]] <- list(
    model = glm(as.formula(fstr), data = demo_data, family = binomial),
    type = "OR", family = "logistic", formula = fstr
  )
}

# Static DP raw linear + BMI (guard on minimum driving-pressure N)
if (sum(!is.na(demo_data$dp)) >= 30) {
  fstr <- paste("dp ~", demo_covars_bmi)
  demo_models[["Static DP"]] <- list(
    model = lm(as.formula(fstr), data = demo_data),
    type = "Beta", family = "linear", formula = fstr
  )
} else {
  message("Skipping Static DP demographic-bias model: < 30 driving-pressure observations")
}

# --- Combined table: rows = covariates, columns = outcomes -------------------
demo_star <- function(p) dplyr::case_when(
  p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ ""
)
demo_term_labels <- c(
  age10              = "Age (per 10 yr)",
  sex_categoryFemale = "Female vs male",
  race_categoryOTHER = "Other vs white",
  race_categoryBLACK = "Black vs white",
  height10           = "Height (per 10 cm)",
  sf10               = "SF ratio (per 10)",
  sofa_total         = "SOFA",
  bmi                = "BMI"
)
demo_outcome_order <- c("VT/PBW", "VT/PFVC (%)", "Mortality",
                        "Static DP", "Ers x PBW", "Ers x PFVC")

demo_cells <- imap_dfr(demo_models, function(m, outcome) {
  broom::tidy(m$model, conf.int = TRUE, exponentiate = m$type == "OR") %>%
    filter(term != "(Intercept)") %>%
    transmute(
      outcome = outcome,
      term,
      cell = sprintf("%.2f [%.2f, %.2f]%s",
                     estimate, conf.low, conf.high, demo_star(p.value))
    )
})

demo_n_row <- tibble(
  term_label = "N",
  outcome = names(demo_models),
  cell = map_chr(demo_models, ~ format(stats::nobs(.x$model), big.mark = ","))
)

demo_table <- demo_cells %>%
  mutate(term_label = recode(term, !!!demo_term_labels)) %>%
  select(term_label, outcome, cell) %>%
  bind_rows(demo_n_row) %>%
  mutate(
    term_label = factor(term_label, levels = c(unname(demo_term_labels), "N")),
    outcome = factor(outcome, levels = demo_outcome_order)
  ) %>%
  arrange(term_label, outcome) %>%
  pivot_wider(names_from = outcome, values_from = cell) %>%
  arrange(term_label)

demo_gt <- demo_table %>%
  mutate(term_label = as.character(term_label)) %>%
  gt::gt(rowname_col = "term_label") %>%
  gt::tab_header(
    title = "Demographic variation in dosing, driving pressure, and elastance",
    subtitle = paste0(site_name,
      " — z-scored outcomes (SD units) except Mortality (OR) and Static DP (cmH2O); ",
      "predictors per 10 units; reference = male, white")
  ) %>%
  gt::sub_missing(missing_text = "")

gt::gtsave(demo_gt, file.path(final_dir, paste0("table_demographic_bias_", site_name, ".html")))
gt::gtsave(demo_gt, file.path(final_dir, paste0("table_demographic_bias_", site_name, ".pdf")))
message("Demographic-bias table written (", length(demo_models), " outcome models)")

# =============================================================================
# 4f3. Predicted FVC vs predicted body weight
# =============================================================================
# Does PBW capture predicted lung size? PFVC regressed on PBW + demographics on
# the BROAD cohort (all eligible patients with height/age/sex/race/PFVC, not just
# the ventilated cross-sectional cohort). Significant age/sex/race coefficients
# indicate PBW alone does not capture predicted lung size.
broad_pfvc <- read_parquet(file.path(output_dir, "analysis_broad_pfvc.parquet")) %>%
  mutate(
    sex_category  = factor(sex_category,  levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
    age10 = age_at_admission / 10
  )

pfvc_vs_pbw_model <- lm(pfvc ~ pbw + age10 + sex_category + race_category, data = broad_pfvc)
pfvc_vs_pbw_formula <- "pfvc ~ pbw + age10 + sex_category + race_category"

pfvc_vs_pbw_gt <- tbl_regression(
  pfvc_vs_pbw_model,
  label = list(
    pbw ~ "PBW (kg)",
    age10 ~ "Age (per 10 yr)",
    sex_category ~ "Sex",
    race_category ~ "Race"
  )
) %>%
  bold_p() %>%
  modify_caption(paste0(site_name,
                        " — predicted FVC vs. predicted body weight (all subjects, N = ",
                        nrow(broad_pfvc), ")"))

pfvc_vs_pbw_gt %>%
  as_gt() %>%
  gt::gtsave(file.path(final_dir, paste0("table_pfvc_vs_pbw_", site_name, ".html")))
pfvc_vs_pbw_gt %>%
  as_gt() %>%
  gt::gtsave(file.path(final_dir, paste0("table_pfvc_vs_pbw_", site_name, ".pdf")))
message("PFVC-vs-PBW table written (N = ", nrow(broad_pfvc), ")")

# =============================================================================
# 4g. Unified long-format regression results table
# =============================================================================
# One row per (model, term). Every estimate is reported on its natural scale —
# OR for logistic regression, HR for the Cox model, Beta for linear regression —
# alongside its 95% CI, p-value, estimate type, variable name, and the full model
# specification. This is the per-cohort artifact that script 05 stacks across
# cohorts to build cross-cohort forest plots, so the column schema must stay
# stable across sites.

extract_model_results <- function(model, estimate_type, analysis,
                                  model_spec, model_family, formula_str) {
  # OR / HR are reported on the exponentiated (ratio) scale; Beta is the raw
  # linear coefficient. broom::tidy applies the matching transform to estimate
  # and CI together so they stay internally consistent.
  exponentiate <- estimate_type %in% c("OR", "HR")
  broom::tidy(model, conf.int = TRUE, exponentiate = exponentiate) %>%
    transmute(
      term,
      estimate,
      conf_low      = conf.low,
      conf_high     = conf.high,
      std_error     = std.error,
      statistic,
      p_value       = p.value,
      estimate_type = estimate_type,
      analysis      = analysis,
      model_spec    = model_spec,
      model_family  = model_family,
      formula       = formula_str,
      n_obs         = stats::nobs(model)
    )
}

results_long <- list()

# Mortality (logistic regression -> odds ratios)
if (!is.null(mortality_models)) {
  results_long <- c(results_long, imap(mortality_models, ~ {
    formula_str <- paste("deceased ~", exposure_specs[[.y]], "+",
                         model_covariates(exposure_specs[[.y]]))
    extract_model_results(.x, "OR", "Mortality",
                          exposure_labels[[.y]], "logistic", formula_str)
  }))
}

# Continuous outcomes (linear regression -> beta coefficients)
for (outcome_name in names(continuous_outcomes)) {
  outcome_var   <- continuous_outcomes[[outcome_name]]$var
  outcome_label <- continuous_outcomes[[outcome_name]]$label
  results_long <- c(results_long, imap(continuous_models[[outcome_name]], ~ {
    formula_str <- paste(outcome_var, "~", exposure_specs[[.y]], "+",
                         model_covariates(outcome_var))
    extract_model_results(.x, "Beta", outcome_label,
                          exposure_labels[[.y]], "linear", formula_str)
  }))
}

# 28-day VFDs (competing-risks Fine-Gray -> subdistribution hazard ratios for
# extubation; estimate_type "HR" so the cross-cohort forest plots it on the log
# scale, like the other ratio outcomes).
results_long <- c(results_long, imap(vfd_cr_models, ~ {
  formula_str <- paste("finegray(Surv(vfd_time, vfd_status) [extubation vs death]) ~",
                       exposure_specs[[.y]], "+", vfd_cr_covariates)
  extract_model_results(.x, "HR", "28-day VFDs",
                        exposure_labels[[.y]], "finegray", formula_str)
}))

# 60-day death hazard (Cox proportional hazards -> hazard ratios). Two exposure
# specs, mirroring the mortality models; HR > 1 = higher death hazard (worse), the
# same direction as the mortality OR.
if (exists("cox_model")) {
  cox_covars  <- "vtpbw + age10 + sex_category + race_category + sf10 + sofa_total"
  results_long <- c(results_long, list(
    # Same vtpbw + pbwpfvc exposure spec as the mortality model — label it with
    # the shared convention so it collapses into one column cross-cohort.
    extract_model_results(cox_model, "HR", "Survival", "VT/PBW + PBW/PFVC", "cox",
                          paste("Surv(surv_time, event) ~ pbwpfvc +", cox_covars)),
    extract_model_results(cox_model_pfvc, "HR", "Survival", "VT/PBW + PFVC", "cox",
                          paste("Surv(surv_time, event) ~ pfvc +", cox_covars))
  ))
}

# Demographic-bias family (one analysis name per outcome, single "Demographics" spec)
results_long <- c(results_long, imap(demo_models, ~ {
  extract_model_results(.x$model, if (.x$type == "OR") "OR" else "Beta",
                        paste0("Demo bias: ", .y), "Demographics",
                        .x$family, .x$formula)
}))

# Predicted FVC vs PBW (broad cohort)
results_long <- c(results_long, list(
  extract_model_results(pfvc_vs_pbw_model, "Beta", "PFVC vs PBW",
                        "PFVC ~ PBW", "linear", pfvc_vs_pbw_formula)
))

# Drop intercepts (not a reportable effect) and stamp the cohort name so the
# table is self-describing once aggregated across sites in script 05.
regression_results_long <- bind_rows(results_long) %>%
  filter(term != "(Intercept)") %>%
  mutate(site = site_name, .before = 1)

write_csv(regression_results_long,
          file.path(final_dir, paste0("regression_results_long_", site_name, ".csv")))
write_parquet(regression_results_long,
              file.path(final_dir, paste0("regression_results_long_", site_name, ".parquet")))

message("Unified regression results table: ", nrow(regression_results_long),
        " rows across ", n_distinct(regression_results_long$analysis), " analyses; saved to ",
        file.path(final_dir, paste0("regression_results_long_", site_name, ".csv")))

# =============================================================================
# 4g2. Residual confounding: E-values + age functional-form (spline) sensitivity
# =============================================================================
# Because PBW and PFVC are DETERMINISTIC functions of {height, age, sex, race},
# those four variables are the complete parent set of every PBW/PFVC-derived
# exposure, and the parents of treatment are a sufficient backdoor adjustment set.
# The confounder space is therefore closed (not open-ended), and residual
# confounding can only enter through the gap between that full parent set and what
# the models actually adjust for. Height is NOT in that gap: in the DAG it reaches
# mortality only through predicted lung size (height -> lung size -> strain ->
# mortality), so it lies on the causal pathway / is the identifying variation in
# the exposure, not a backdoor — conditioning on it is over-adjustment. The single
# enumerable residual confounder is therefore the FUNCTIONAL FORM of age (a genuine
# confounder, entered linearly in the main models). Anything beyond that would have
# to be a truly unmeasured common cause acting outside this deterministic structure
# — exactly what the E-value bounds.
#
# This section reports, for every ratio-scale exposure estimate (mortality OR,
# 60-day survival HR, 28-day liberation Fine-Gray SHR):
#   1. the per-1-SD estimate under LINEAR age (the main-model specification),
#   2. the per-1-SD estimate under SPLINE age (race + ns(age, 4) + sex + SOFA + SF):
#      the direct bound on nonlinear-age confounding — an estimate unchanged under
#      flexible age is robust to it; one that collapses was carrying age, and
#   3. the E-value (point and CI-limit) for the linear-age estimate: the minimum
#      association, on the risk-ratio scale, an unmeasured confounder would need
#      with BOTH the exposure and the outcome to explain the estimate (or, for the
#      CI E-value, to move the CI to include the null).
#
# Estimates are rescaled to a 1-SD increase in the exposure (E-values on the raw
# per-unit scale would sit artificially close to 1). Mortality and liberation are
# common outcomes, so each OR/HR is converted to an approximate risk ratio before
# the E-value (rare = FALSE; VanderWeele & Ding, Ann Intern Med 2017); the
# Fine-Gray subdistribution HR uses the same common-outcome HR conversion.

EXPOSURE_TERMS <- c("vtpfvc", "vtpbw", "pfvc", "pbwpfvc")
exposure_term_labels_ev <- c(vtpfvc = "VT/PFVC", vtpbw = "VT/PBW",
                             pfvc = "PFVC", pbwpfvc = "PBW/PFVC")

# Spline-age covariate set (mirrors the linear-age set, age10 -> ns(age, 4)).
# Reused for the mortality and Fine-Gray refits below; the Cox refit is spelled
# out separately because it carries its own exposure terms.
spline_covars <- "race_category + ns(age_at_admission, 4) + sex_category + sofa_total + sf10"

# Per-SD-rescaled ratio estimate + 95% CI from a fitted model. beta and its SE
# are the log-OR / log-HR (link scale); multiplying by the exposure SD gives the
# per-SD log-ratio, exponentiated back to the ratio scale.
persd_ratio <- function(model, term, exposure_sd) {
  beta <- coef(model)[[term]]
  se   <- sqrt(diag(vcov(model)))[[term]]
  list(
    estimate  = exp(beta * exposure_sd),
    conf_low  = exp((beta - 1.96 * se) * exposure_sd),
    conf_high = exp((beta + 1.96 * se) * exposure_sd)
  )
}

# Per-SD estimate + CI for every exposure term in a model (no E-value). The SD is
# the marginal, patient-level SD of the exposure in sd_data, so the linear- and
# spline-age fits are compared on an identical contrast.
persd_terms <- function(model, sd_data, terms = EXPOSURE_TERMS) {
  present <- intersect(terms, names(coef(model)))
  map_dfr(present, function(tm) {
    exposure_sd <- sd(sd_data[[tm]], na.rm = TRUE)
    pr <- persd_ratio(model, tm, exposure_sd)
    tibble(term = tm, exposure_sd = exposure_sd,
           estimate = pr$estimate, conf_low = pr$conf_low, conf_high = pr$conf_high)
  })
}

# Point and CI E-values for one per-SD ratio estimate. type is "OR" or "HR";
# both use rare = FALSE (common outcomes). The CI E-value is the non-NA bound
# the package returns (the limit nearest the null; 1 if the CI crosses it).
evalue_for <- function(est, lo, hi, type) {
  ev_mat <- switch(type,
    OR = EValue::evalues.OR(est, lo, hi, rare = FALSE),
    HR = EValue::evalues.HR(est, lo, hi, rare = FALSE)
  )
  ev    <- ev_mat["E-values", ]
  ci_ev <- ev[c("lower", "upper")]
  ci_ev <- ci_ev[!is.na(ci_ev)]
  list(evalue_point = unname(ev[["point"]]),
       evalue_ci    = if (length(ci_ev)) unname(ci_ev[1]) else NA_real_)
}

# One model's exposure rows: per-SD estimate under linear age (model_lin) and
# under spline age (model_spl) side by side, plus the E-value for the linear-age
# estimate. Both fits use the same sd_data so the per-SD contrast is identical.
residual_conf_rows <- function(model_lin, model_spl, model_spec, analysis, type, sd_data) {
  lin <- persd_terms(model_lin, sd_data)
  spl <- persd_terms(model_spl, sd_data)
  if (nrow(lin) == 0) return(tibble())
  ev  <- pmap_dfr(list(lin$estimate, lin$conf_low, lin$conf_high),
                  function(e, l, h) as_tibble(evalue_for(e, l, h, type)))
  spl_i <- match(lin$term, spl$term)
  tibble(
    analysis      = analysis,
    model_spec    = model_spec,
    term          = lin$term,
    term_label    = unname(exposure_term_labels_ev[lin$term]),
    estimate_type = type,
    exposure_sd   = lin$exposure_sd,
    est_linage    = lin$estimate,
    lo_linage     = lin$conf_low,
    hi_linage     = lin$conf_high,
    est_splineage = spl$estimate[spl_i],
    lo_splineage  = spl$conf_low[spl_i],
    hi_splineage  = spl$conf_high[spl_i],
    evalue_point  = ev$evalue_point,
    evalue_ci     = ev$evalue_ci
  )
}

# --- Spline-age refits of each ratio model (age10 -> ns(age_at_admission, 4)) ---
# The exposure specs are dosing ratios, never DP-derived, so no BMI is added.

mortality_models_spline <- if (!is.null(mortality_models)) {
  map(exposure_specs, ~ glm(
    as.formula(paste("deceased ~", .x, "+", spline_covars)),
    data = cross_sectional, family = binomial))
} else NULL

cox_model_spline <- if (exists("cox_model")) {
  coxph(Surv(surv_time, event) ~ pbwpfvc + vtpbw + ns(age_at_admission, 4) +
          sex_category + race_category + sf10 + sofa_total, data = surv_data)
} else NULL
cox_model_pfvc_spline <- if (exists("cox_model_pfvc")) {
  coxph(Surv(surv_time, event) ~ pfvc + vtpbw + ns(age_at_admission, 4) +
          sex_category + race_category + sf10 + sofa_total, data = surv_data)
} else NULL

# Fine-Gray refit with spline age (mirrors fit_vfd_finegray from section 4d2).
fit_vfd_finegray_spline <- function(exposure_spec) {
  model_rhs <- paste(exposure_spec, "+", spline_covars)
  rhs_vars  <- all.vars(as.formula(paste("~", model_rhs)))
  df <- cross_sectional %>%
    mutate(vfd_status_f = factor(vfd_status, levels = c(0, 1, 2),
                                 labels = c("censored", "extubation", "death"))) %>%
    select(vfd_time, vfd_status_f, all_of(rhs_vars)) %>%
    filter(!is.na(vfd_time), vfd_time > 0, !is.na(vfd_status_f)) %>%
    drop_na(all_of(rhs_vars))
  fg <- survival::finegray(survival::Surv(vfd_time, vfd_status_f) ~ ., data = df,
                           etype = "extubation")
  survival::coxph(
    as.formula(paste0("survival::Surv(fgstart, fgstop, fgstatus) ~ ", model_rhs)),
    weights = fgwt, data = fg
  )
}
vfd_cr_models_spline <- map(exposure_specs, fit_vfd_finegray_spline)

# --- Assemble residual-confounding rows across all ratio models ----------------
residual_list <- list()

# In-hospital mortality (logistic, OR) — one model per exposure specification.
if (!is.null(mortality_models)) {
  residual_list <- c(residual_list, imap(mortality_models,
    ~ residual_conf_rows(.x, mortality_models_spline[[.y]], exposure_labels[[.y]],
                         "Mortality (in-hospital)", "OR", cross_sectional)))
}

# 60-day death hazard (Cox, HR) — PBW/PFVC and PFVC exposure specs.
if (!is.null(cox_model_spline)) {
  residual_list <- c(residual_list, list(
    residual_conf_rows(cox_model, cox_model_spline, "VT/PBW + PBW/PFVC",
                       "Survival (60-day)", "HR", surv_data),
    residual_conf_rows(cox_model_pfvc, cox_model_pfvc_spline, "VT/PBW + PFVC",
                       "Survival (60-day)", "HR", surv_data)))
}

# 28-day ventilator liberation (Fine-Gray subdistribution HR). The finegray()
# expansion duplicates rows per subject, so the per-SD contrast uses the
# patient-level SD from the cross-sectional cohort (rows with a valid VFD time),
# not the expanded risk set.
vfd_sd_data <- cross_sectional %>% filter(!is.na(vfd_time), vfd_time > 0)
residual_list <- c(residual_list, imap(vfd_cr_models,
  ~ residual_conf_rows(.x, vfd_cr_models_spline[[.y]], exposure_labels[[.y]],
                       "28-day VFDs (liberation)", "HR", vfd_sd_data)))

residual_confounding <- bind_rows(residual_list) %>%
  mutate(analysis = factor(analysis, levels = c(
    "Survival (60-day)", "Mortality (in-hospital)", "28-day VFDs (liberation)"))) %>%
  arrange(analysis) %>%               # stable: preserves model & term order within analysis
  mutate(analysis = as.character(analysis)) %>%
  mutate(site = site_name, .before = 1)

# Filename kept as `evalues_<site>` for continuity; the table now also carries the
# linear- vs spline-age estimates alongside the E-values.
write_csv(residual_confounding, file.path(final_dir, paste0("evalues_", site_name, ".csv")))

# Rendered table: grouped by analysis, one row per (model spec, exposure).
evalue_gt <- residual_confounding %>%
  transmute(
    analysis,
    Model                 = model_spec,
    Exposure              = term_label,
    Type                  = estimate_type,
    `Linear-age estimate` = sprintf("%.2f (%.2f, %.2f)", est_linage, lo_linage, hi_linage),
    `Spline-age estimate` = sprintf("%.2f (%.2f, %.2f)",
                                    est_splineage, lo_splineage, hi_splineage),
    `E-value (estimate)`  = sprintf("%.2f", evalue_point),
    `E-value (95% CI)`    = ifelse(is.na(evalue_ci), "—", sprintf("%.2f", evalue_ci))
  ) %>%
  gt::gt(groupname_col = "analysis") %>%
  gt::tab_header(
    title = "Residual confounding: per-SD estimates, age functional-form sensitivity, and E-values",
    subtitle = paste0(
      site_name,
      " — per-1-SD exposure contrasts. Linear- vs spline-age estimates bound the ",
      "only enumerable residual confounder (nonlinear age); the exposures are ",
      "deterministic in {height, age, sex, race} and height lies on the causal ",
      "pathway, not a backdoor. The E-value is the minimum risk-ratio association ",
      "an unmeasured confounder would need with both the exposure and the outcome ",
      "to explain the estimate (CI E-value: to move the CI to include the null).")
  ) %>%
  gt::cols_align("left", columns = c(Model, Exposure)) %>%
  gt::sub_missing(missing_text = "—")

gt::gtsave(evalue_gt, file.path(final_dir, paste0("table_evalues_", site_name, ".html")))
gt::gtsave(evalue_gt, file.path(final_dir, paste0("table_evalues_", site_name, ".pdf")))
message("Residual-confounding table written (", nrow(residual_confounding),
        " ratio-scale exposure estimates across ",
        n_distinct(residual_confounding$analysis),
        " analyses; linear/spline age + E-values)")

# =============================================================================
# 4h. Conditional bias diagnostic plots (algorithmDiagnostics)
# =============================================================================

# Prepare data with all needed variables
diag_data <- cross_sectional %>%
  filter(!is.na(age_group))

# --- Plot 1: PBW as prediction of PFVC, stratified by demographics ---
bias_pbw_pfvc <- conditional_bias_plot(
  data = diag_data,
  dependent_vars = "pfvc",
  independent_vars = "pbw",
  grouping_vars = c("race_category", "sex_category", "age_group"),
  dep_var_labels = c("PFVC (L)"), min_n = 20
)

ggsave(file.path(final_dir, paste0("bias_pbw_vs_pfvc_", site_name, ".pdf")),
       bias_pbw_pfvc, width = 14, height = 12)

# --- Plot 2: Mortality by VT/PFVC and VT/PBW, stratified by demographics ---
bias_mortality <- conditional_bias_plot(
  data = diag_data,
  dependent_vars = "deceased",
  independent_vars = c("vtpfvc", "vtpbw", "pbwpfvc"),
  grouping_vars = c("race_category", "sex_category", "age_group"),
  dep_var_labels = "Mortality",
  x_labels = c("Percentile of VT/PFVC", "Percentile of VT/PBW","Percentile of PBW/PFVC")
)

ggsave(file.path(final_dir, paste0("bias_mortality_", site_name, ".pdf")),
       bias_mortality, width = 14, height = 18)

# --- Plot 3: Elastance by exposures, stratified by demographics ---
bias_ers <- conditional_bias_plot(
    data = diag_data,
    dependent_vars = "ers",
    independent_vars = c("vtpfvc", "vtpbw"),
    grouping_vars = c("race_category", "sex_category", "age_group"),
    dep_var_labels = "Elastance (cmH2O/L)",
    x_labels = c("Percentile of VT/PFVC", "Percentile of VT/PBW")
  )

  ggsave(file.path(final_dir, paste0("bias_elastance_", site_name, ".pdf")),
         bias_ers, width = 14, height = 18)
  message("Elastance bias plots saved")

# --- Plot 4: Compliance by exposures, stratified by demographics ---
bias_crs <- conditional_bias_plot(
    data = diag_data,
    dependent_vars = "crs",
    independent_vars = c("pbwpfvc"),
    grouping_vars = c("race_category", "sex_category", "age_group"),
    dep_var_labels = "Compliance (mL/cmH2O)",
    x_labels = c("Percentile of PBW/PFVC")
  )

  ggsave(file.path(final_dir, paste0("bias_compliance_", site_name, ".pdf")),
         bias_crs, width = 14, height = 18)
  message("Compliance bias plots saved")

# --- Plot 5: VFD-28 by exposures, stratified by demographics ---
bias_vfd <- conditional_bias_plot(
  data = diag_data,
  dependent_vars = "vfd_28",
  independent_vars = c("pbwpfvc"),
  grouping_vars = c("race_category", "sex_category", "age_group"),
  dep_var_labels = "28-day VFDs", min_n = 20
)

ggsave(file.path(final_dir, paste0("bias_vfd28_", site_name, ".pdf")),
       bias_vfd, width = 14, height = 18)

# --- Plot 6: PBW/PFVC ratio by demographics ---
bias_pbwpfvc <- conditional_bias_plot(
  data = diag_data,
  dependent_vars = "pbwpfvc",
  independent_vars = c("pfvc", "pbw"),
  grouping_vars = c("race_category", "sex_category", "age_group"),
  dep_var_labels = "PBW/PFVC Ratio",
  x_labels = c("Percentile of PFVC", "Percentile of PBW")
)

ggsave(file.path(final_dir, paste0("bias_pbwpfvc_ratio_", site_name, ".pdf")),
       bias_pbwpfvc, width = 14, height = 18)

message("All bias diagnostic plots saved")

# NOTE: the federated per-percentile conditional-bias export that drives the POOLED
# cross-cohort bias plots lives in a SEPARATE, deferred script — code/
# cbias_federated_export.R — not this pipeline. After stratification some
# (stratum x percentile) cells fall below n = 10, so those aggregates must be run
# through the consortium's deterministic additive-masking pipeline before they can
# leave a site, which is held until all sites confirm participation. The pooled
# plots are likewise deferred in code/cbias_pooled_plots.R.

# =============================================================================
# 4i. Inclusion CONSORT diagram + PBW:PFVC-by-demographics figure (site QC)
# =============================================================================
# CONSORT flow from the 7-step attrition log written in script 03.
attrition <- read_csv(file.path(final_dir, paste0("attrition_log_", site_name, ".csv")),
                      show_col_types = FALSE)
consort_fig <- render_consort(attrition, title = paste0("Cohort inclusion - ", site_name))
ggsave(file.path(final_dir, paste0("consort_diagram_", site_name, ".pdf")),
       consort_fig, width = 9, height = 11)
message("CONSORT diagram saved")

# PBW:PFVC by demographics — site QC figure from local data. The federated
# histogram/quantile exports written in script 03 drive the POOLED figure in
# script 05; this local version uses raw points (never leaves the site).
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
dist_base <- cross_sectional %>% filter(!is.na(pbwpfvc))

p_sex <- ggplot(dist_base, aes(sex_category, pbwpfvc, fill = sex_category)) +
  geom_violin(alpha = 0.5, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = okabe, guide = "none") +
  labs(x = "Sex", y = "PBW:PFVC ratio (kg/L)") + theme_minimal(base_size = 11)

p_race <- ggplot(dist_base, aes(race_category, pbwpfvc, fill = race_category)) +
  geom_violin(alpha = 0.5, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = okabe, guide = "none") +
  labs(x = "Race", y = "PBW:PFVC ratio (kg/L)") + theme_minimal(base_size = 11)

p_age <- ggplot(dist_base, aes(age_at_admission, pbwpfvc)) +
  geom_point(alpha = 0.15, colour = "#D55E00") +
  geom_smooth(method = "loess", colour = "#0072B2", se = TRUE) +
  labs(x = "Age (years)", y = "PBW:PFVC ratio (kg/L)") + theme_minimal(base_size = 11)

p_height <- ggplot(dist_base, aes(height_cm, pbwpfvc)) +
  geom_point(alpha = 0.15, colour = "#D55E00") +
  geom_smooth(method = "loess", colour = "#0072B2", se = TRUE) +
  labs(x = "Height (cm)", y = "PBW:PFVC ratio (kg/L)") + theme_minimal(base_size = 11)

dist_fig <- (p_sex | p_race) / (p_age | p_height) +
  plot_annotation(title = paste0("PBW:PFVC by demographics - ", site_name),
                  theme = theme(plot.title = element_text(face = "bold")))
ggsave(file.path(final_dir, paste0("distribution_pbwpfvc_", site_name, ".pdf")),
       dist_fig, width = 11, height = 9)
message("PBW:PFVC distribution figure saved")

message("All outputs saved to: ", final_dir)
message("Script 04 complete.")
