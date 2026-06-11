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
library(patchwork)
library(broom)
library(algorithmDiagnostics)

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
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER"))
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
covariates <- "race_category + age_at_admission + sex_category + sofa_total + sf_ratio"

# Per the original paper, any model whose outcome OR exposure is derived from
# driving pressure (static DP, elastance, compliance, and the elastance-normalized
# Ers x PBW / Ers x PFVC) is additionally adjusted for BMI. Models without a
# driving-pressure component use the standard covariate set.
DP_DERIVED <- c("dp", "ers", "crs", "ers_pbw", "ers_pfvc")
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
    tbl_regression(.x, exponentiate = TRUE) %>% bold_p()
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

# --- Mortality vs elastance normalized to PBW / PFVC -------------------------
# Two additional logistic models with elastance-normalized exposures (Ers x PBW =
# ers * pbw, Ers x PFVC = ers * pfvc), adjusted for the standard covariates.
ers_mortality_specs <- c(ers_pbw = "Ers x PBW", ers_pfvc = "Ers x PFVC")
ers_mortality_models <- list()
if (has_mortality_variation) {
  for (v in names(ers_mortality_specs)) {
    fstr <- paste("deceased ~", v, "+", model_covariates(v))
    ers_mortality_models[[v]] <- glm(as.formula(fstr), data = cross_sectional,
                                      family = binomial)
  }

  message("Logistic regression (elastance-normalized mortality) — AIC:")
  iwalk(ers_mortality_models,
        ~ message("  ", ers_mortality_specs[.y], ": ", round(AIC(.x), 1)))

  # Standalone merged table, mirroring the primary mortality regression output.
  ers_mortality_tables <- map(ers_mortality_models,
                              ~ tbl_regression(.x, exponentiate = TRUE) %>% bold_p())
  tbl_merge(ers_mortality_tables,
            tab_spanner = unname(ers_mortality_specs[names(ers_mortality_models)])) %>%
    as_gt() %>%
    gt::gtsave(file.path(final_dir, paste0("regression_ers_mortality_", site_name, ".html")))
}

# =============================================================================
# 4d. Linear regression — continuous outcomes (elastance, compliance, VFD-28, DP)
# =============================================================================

continuous_outcomes <- list(
  ers    = list(var = "ers",    label = "Elastance"),
  crs    = list(var = "crs",    label = "Compliance"),
  vfd_28 = list(var = "vfd_28", label = "28-day VFDs"),
  dp     = list(var = "dp",     label = "Static DP")
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
    tbl_regression(.x) %>% bold_p()
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

# Continuous outcomes
for (outcome_name in names(continuous_outcomes)) {
  aic_results[[continuous_outcomes[[outcome_name]]$label]] <- tibble(
    exposure = exposure_labels,
    AIC = map_dbl(continuous_models[[outcome_name]], AIC),
    is_reference = exposure == "VT/PBW"
  )
}

# Mortality, elastance-normalized exposures. Ers x PBW and Ers x PFVC are fit on
# the same support (both require ers), so their AICs are mutually comparable. The
# evidence ratio is referenced to the PBW-scaled model (Ers x PBW) — the same
# PBW-reference convention used for the VT/PBW columns in the other outcomes.
if (length(ers_mortality_models) > 0) {
  aic_results[["Mortality (Ers-normalized)"]] <- tibble(
    exposure = unname(ers_mortality_specs[names(ers_mortality_models)]),
    AIC = map_dbl(ers_mortality_models, AIC),
    is_reference = exposure == "Ers x PBW"
  )
}

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
  ungroup()

message("AIC comparison (evidence ratios vs the PBW-scaled reference within each outcome):")
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
    name = "Evidence ratio\n(vs PBW reference)",
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(log10(ER_FLOOR), log10(ER_CEIL)),
    breaks = -3:3, labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000")
  ) +
  labs(
    title = "Evidence ratios across models and outcomes",
    subtitle = "Each cell vs the PBW-scaled reference within that outcome (VT/PBW; Ers x PBW for the elastance-normalized mortality models); truncated to [0.001, 1000]",
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
    Surv(surv_time, event) ~ pbwpfvc + vtpbw + age_at_admission +
      sex_category + race_category + sf_ratio + sofa_total,
    data = surv_data
  )

  message("Cox model:")
  print(summary(cox_model))

  km_fit <- survfit(Surv(surv_time, event) ~ ntile(pbwpfvc,3), data = surv_data)

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
    legend.title = "PBW/PFVC Tercile",
    ggtheme = theme_minimal()
  )

  pdf(file.path(final_dir, paste0("km_curves_", site_name, ".pdf")),
      width = 10, height = 8)
  print(km_plot)
  dev.off()

  message("KM curves saved")

  sink(file.path(final_dir, paste0("cox_model_summary_", site_name, ".txt")))
  print(summary(cox_model))
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
# 4f3. H1 table: predicted FVC vs predicted body weight
# =============================================================================
# PFVC regressed on PBW + demographics on the BROAD cohort (all eligible patients
# with height/age/sex/race/PFVC, not just the ventilated cross-sectional cohort).
broad_pfvc <- read_parquet(file.path(output_dir, "analysis_broad_pfvc.parquet")) %>%
  mutate(
    sex_category  = factor(sex_category,  levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
    age10 = age_at_admission / 10
  )

h1_model <- lm(pfvc ~ pbw + age10 + sex_category + race_category, data = broad_pfvc)
h1_formula <- "pfvc ~ pbw + age10 + sex_category + race_category"

h1_gt <- tbl_regression(
  h1_model,
  label = list(
    pbw ~ "PBW (kg)",
    age10 ~ "Age (per 10 yr)",
    sex_category ~ "Sex",
    race_category ~ "Race"
  )
) %>%
  bold_p() %>%
  modify_caption(paste0("H1. ", site_name,
                        " — predicted FVC vs. predicted body weight (all subjects, N = ",
                        nrow(broad_pfvc), ")"))

h1_gt %>%
  as_gt() %>%
  gt::gtsave(file.path(final_dir, paste0("table_h1_pfvc_vs_pbw_", site_name, ".html")))
h1_gt %>%
  as_gt() %>%
  gt::gtsave(file.path(final_dir, paste0("table_h1_pfvc_vs_pbw_", site_name, ".pdf")))
message("H1 PFVC-vs-PBW table written (N = ", nrow(broad_pfvc), ")")

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

# Survival (Cox proportional hazards -> hazard ratios)
if (exists("cox_model")) {
  cox_formula <- paste(
    "Surv(surv_time, event) ~ pbwpfvc + vtpbw + age_at_admission +",
    "sex_category + race_category + sf_ratio + sofa_total"
  )
  results_long <- c(results_long, list(
    extract_model_results(cox_model, "HR", "Survival",
                          "PBW/PFVC + VT/PBW", "cox", cox_formula)
  ))
}

# Mortality vs elastance normalized to PBW / PFVC (logistic -> odds ratios)
if (length(ers_mortality_models) > 0) {
  results_long <- c(results_long, imap(ers_mortality_models, ~ {
    fstr <- paste("deceased ~", .y, "+", model_covariates(.y))
    extract_model_results(.x, "OR", "Mortality",
                          ers_mortality_specs[[.y]], "logistic", fstr)
  }))
}

# Demographic-bias family (one analysis name per outcome, single "Demographics" spec)
results_long <- c(results_long, imap(demo_models, ~ {
  extract_model_results(.x$model, if (.x$type == "OR") "OR" else "Beta",
                        paste0("Demo bias: ", .y), "Demographics",
                        .x$family, .x$formula)
}))

# H1: PFVC vs PBW (broad cohort)
results_long <- c(results_long, list(
  extract_model_results(h1_model, "Beta", "H1: PFVC ~ PBW",
                        "PFVC ~ PBW", "linear", h1_formula)
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
