# =============================================================================
# Script 05: Cross-cohort aggregation and forest plots
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Aggregates the per-cohort long-format regression tables written by script 04
# (one `regression_results_long_<site>.csv` per cohort, under each site's
# output/<site>_output/final/ directory) and produces, for each analysis, a
# forest plot comparing effect estimates across cohorts.
#
# Unlike scripts 01-04, this script is intentionally site-agnostic: it discovers
# every cohort's results table on disk and stacks them. Run it after script 04
# has been run for each cohort whose results you want to compare (the tables can
# be collected from multiple sites into the local output/ tree).

library(tidyverse)
library(arrow)
library(here)
library(patchwork)

source("utils/consort_diagram.R")

# =============================================================================
# Discover and load every cohort's results table
# =============================================================================

result_files <- Sys.glob(
  here("output", "*_output", "final", "regression_results_long_*.csv")
)

if (length(result_files) == 0) {
  stop("No per-cohort results tables found. Run script 04 first. Expected files ",
       "matching output/*_output/final/regression_results_long_*.csv")
}

message("Found ", length(result_files), " cohort result table(s):")
walk(result_files, ~ message("  ", .x))

all_results <- map_dfr(result_files, read_csv, show_col_types = FALSE)

# Defend against a stale table written before the `site` column existed.
if (!"site" %in% names(all_results)) {
  stop("Results tables are missing the `site` column. Re-run script 04 to ",
       "regenerate them with the current schema.")
}

# The Cox survival model is the same vtpbw + pbwpfvc spec as the mortality model,
# but older script-04 versions labelled it with the operands reversed. Canonicalise
# so it collapses into a single "VT/PBW + PBW/PFVC" column across cohorts.
all_results <- all_results %>%
  mutate(model_spec = recode(model_spec, "PBW/PFVC + VT/PBW" = "VT/PBW + PBW/PFVC"))

message("Aggregated ", nrow(all_results), " rows across ",
        n_distinct(all_results$site), " cohort(s): ",
        paste(sort(unique(all_results$site)), collapse = ", "))

# =============================================================================
# Output locations
# =============================================================================

cross_dir <- here("output", "cross_cohort")
forest_dir <- file.path(cross_dir, "forest_plots")
dir.create(forest_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(all_results, file.path(cross_dir, "regression_results_all_cohorts.csv"))
write_parquet(all_results, file.path(cross_dir, "regression_results_all_cohorts.parquet"))

# =============================================================================
# Forest plot setup
# =============================================================================

# Exposures of interest. The forest grid shows one row per exposure and one column
# per model specification; covariate (adjustment) coefficients are not plotted but
# remain available in the aggregated CSV.
exposure_term_labels <- c(
  vtpfvc     = "VT/PFVC",
  vtpbw      = "VT/PBW",
  pfvc       = "PFVC",
  pbwpfvc    = "PBW/PFVC"
)

# Row order (exposures) and column order (models, simplest to most complex).
term_order  <- c("VT/PFVC", "VT/PBW", "PFVC", "PBW/PFVC")
model_order <- c("VT/PFVC", "VT/PBW", "VT/PFVC + VT/PBW", "VT/PBW + PFVC",
                 "VT/PBW + PBW/PFVC")

# Reader-facing outcome labels (the regression tables store the terse mortality /
# survival names; spell out the definition wherever an outcome is shown). The two
# mortality endpoints differ: the logistic model is in-hospital mortality, the Cox
# model is 60-day all-cause survival.
outcome_display <- c(
  "Mortality" = "Mortality (in-hospital)",
  "Survival"  = "Survival (60-day)"
)
relabel_outcome <- function(x) unname(coalesce(outcome_display[as.character(x)], as.character(x)))

forest_data <- all_results %>%
  filter(term %in% names(exposure_term_labels)) %>%
  mutate(
    term_label = recode(term, !!!exposure_term_labels),
    term_label = factor(term_label, levels = term_order),
    model_spec = factor(model_spec,
                        levels = unique(c(model_order, sort(unique(model_spec)))))
  )

if (nrow(forest_data) == 0) {
  stop("No exposure-term rows to plot. Check that script 04 produced models for ",
       "the expected exposures: ", paste(names(exposure_term_labels), collapse = ", "))
}

# Okabe-Ito palette for discrete cohorts (per project convention). Extend by
# interpolation only if there are more cohorts than palette colors.
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#000000")
sites <- sort(unique(forest_data$site))
site_colors <- if (length(sites) > length(okabe_ito)) {
  setNames(grDevices::colorRampPalette(okabe_ito)(length(sites)), sites)
} else {
  setNames(okabe_ito[seq_along(sites)], sites)
}

# =============================================================================
# Random-effects meta-analysis of every regression coefficient across cohorts
# =============================================================================
# Pool the per-site coefficients for each (analysis, model_spec, term) with a
# random-effects (REML) model. Ratio estimates (OR/HR) are pooled on the log
# scale and back-transformed; linear Beta coefficients are pooled on their natural
# scale. The SE is the link-scale std_error from each site (Wald), with the CI
# width as a fallback. Cross-site heterogeneity (I^2, tau^2, Q p-value) is reported
# alongside each pooled estimate. The pooled estimates are written to a table and
# overlaid as diamonds on the per-cohort forest plots below.
z975 <- qnorm(0.975)
pool_estimates <- function(df) {
  is_ratio <- df$estimate_type[1] %in% c("OR", "HR")
  yi    <- if (is_ratio) log(df$estimate) else df$estimate
  ci_se <- (if (is_ratio) log(df$conf_high) - log(df$conf_low)
            else df$conf_high - df$conf_low) / (2 * z975)
  sei   <- ifelse(is.finite(df$std_error) & df$std_error > 0, df$std_error, ci_se)
  keep  <- is.finite(yi) & is.finite(sei) & sei > 0
  yi <- yi[keep]; sei <- sei[keep]; k <- length(yi)
  if (k == 0) return(tibble())
  bt   <- if (is_ratio) exp else identity
  base <- tibble(estimate_type = df$estimate_type[1])
  if (k == 1) {
    return(bind_cols(base, tibble(
      k = 1L, pooled = bt(yi), pooled_lo = bt(yi - z975 * sei),
      pooled_hi = bt(yi + z975 * sei), I2 = NA_real_, tau2 = NA_real_, het_p = NA_real_)))
  }
  fit <- metafor::rma(yi = yi, sei = sei, method = "REML")
  bind_cols(base, tibble(
    k = k, pooled = bt(as.numeric(fit$beta)), pooled_lo = bt(fit$ci.lb),
    pooled_hi = bt(fit$ci.ub), I2 = fit$I2, tau2 = fit$tau2, het_p = fit$QEp))
}

pooled_coefs <- all_results %>%
  filter(term != "(Intercept)") %>%
  group_by(analysis, model_spec, term) %>%
  group_modify(~ pool_estimates(.x)) %>%
  ungroup()

write_csv(pooled_coefs, file.path(cross_dir, "pooled_coefficients_all_cohorts.csv"))
write_parquet(pooled_coefs, file.path(cross_dir, "pooled_coefficients_all_cohorts.parquet"))

# Headline table: pooled exposure coefficients, grouped by analysis.
pooled_exposure_gt <- pooled_coefs %>%
  filter(term %in% names(exposure_term_labels)) %>%
  mutate(
    term_label = recode(term, !!!exposure_term_labels),
    term_label = factor(term_label, levels = term_order),
    model_spec = factor(model_spec, levels = unique(c(model_order, sort(unique(model_spec))))),
    analysis   = factor(analysis, levels = sort(unique(analysis)))
  ) %>%
  arrange(analysis, model_spec, term_label) %>%
  transmute(
    analysis = relabel_outcome(analysis),
    Model = model_spec, Exposure = term_label, Type = estimate_type, Sites = k,
    `Pooled estimate (95% CI)` = sprintf("%.2f (%.2f, %.2f)", pooled, pooled_lo, pooled_hi),
    `I-squared (%)` = ifelse(is.na(I2),   "—", sprintf("%.0f", I2)),
    `tau-squared`   = ifelse(is.na(tau2), "—", sprintf("%.3f", tau2)),
    `Het. p`        = ifelse(is.na(het_p), "—", formatC(het_p, format = "g", digits = 2))
  ) %>%
  gt::gt(groupname_col = "analysis") %>%
  gt::tab_header(
    title = "Pooled exposure coefficients (random-effects meta-analysis)",
    subtitle = paste0(
      "Per-site coefficients pooled across ", length(sites), " cohort(s) (REML); ",
      "OR/HR pooled on the log scale, Beta on the natural scale. High I-squared / ",
      "low heterogeneity p flags cross-site effect heterogeneity.")
  ) %>%
  gt::cols_align("left", columns = c(Model, Exposure)) %>%
  gt::sub_missing(missing_text = "—")

gt::gtsave(pooled_exposure_gt, file.path(cross_dir, "table_pooled_coefficients.html"))
gt::gtsave(pooled_exposure_gt, file.path(cross_dir, "table_pooled_coefficients.pdf"))
message("Pooled coefficient table written (", nrow(pooled_coefs),
        " pooled estimates across ", n_distinct(pooled_coefs$analysis), " analyses)")

# Pooled exposure coefficients prepared for overlay on the forest plots — same
# term_label / model_spec factor levels as forest_data so the facets line up.
pooled_forest_data <- pooled_coefs %>%
  filter(term %in% names(exposure_term_labels)) %>%
  mutate(
    term_label = recode(term, !!!exposure_term_labels),
    term_label = factor(term_label, levels = term_order),
    model_spec = factor(model_spec, levels = levels(forest_data$model_spec))
  )

# Build a forest plot for one analysis (a single estimate_type by construction).
# Per-site estimates are coloured points; the random-effects pooled estimate is
# overlaid as a black diamond on a "Pooled (RE)" row at the foot of each panel.
make_forest_plot <- function(analysis_data, analysis_name, pooled_data = NULL) {
  estimate_type <- analysis_data$estimate_type[1]
  is_ratio <- estimate_type %in% c("OR", "HR")
  null_value <- if (is_ratio) 1 else 0

  combined <- analysis_data %>%
    transmute(term_label, model_spec, group = site,
              estimate, conf_low, conf_high, kind = "Site")
  if (!is.null(pooled_data) && nrow(pooled_data) > 0) {
    combined <- bind_rows(combined, pooled_data %>%
      transmute(term_label, model_spec, group = "Pooled (RE)",
                estimate = pooled, conf_low = pooled_lo, conf_high = pooled_hi,
                kind = "Pooled"))
  }
  # "Pooled (RE)" first → rendered at the foot of each panel, below the cohorts.
  group_levels <- c("Pooled (RE)", setdiff(sort(unique(combined$group)), "Pooled (RE)"))
  combined <- combined %>% mutate(group = factor(group, levels = group_levels))
  pal <- c(site_colors, `Pooled (RE)` = "#000000")

  plot <- ggplot(combined, aes(x = estimate, y = group, color = group)) +
    geom_vline(xintercept = null_value, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.25) +
    geom_point(aes(size = kind, shape = kind)) +
    # Grid: one row per exposure, one column per model specification. Empty cells
    # mark exposures a given model does not include. scales = "free_x" lets each
    # model column take its own estimate range; drop = FALSE keeps the full grid.
    facet_grid(term_label ~ model_spec, scales = "free_x", drop = FALSE) +
    scale_color_manual(values = pal, guide = "none") +
    scale_shape_manual(values = c(Site = 16, Pooled = 18), guide = "none") +
    scale_size_manual(values = c(Site = 2.3, Pooled = 3.4), guide = "none") +
    labs(
      title = paste0("Cross-cohort results: ", relabel_outcome(analysis_name)),
      subtitle = paste0(estimate_type, " (95% CI); black diamond = random-effects ",
                        "pooled estimate; dashed line = null (", null_value, ")"),
      x = paste0(estimate_type, " (95% CI)"),
      y = "Cohort"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text.x = element_text(face = "bold"),
      strip.text.y = element_text(face = "bold", angle = 0),
      panel.spacing = unit(0.6, "lines"),
      panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.5)
    )

  # Ratio estimates (OR/HR) read correctly on a log axis.
  if (is_ratio) plot <- plot + scale_x_log10()
  plot
}

# =============================================================================
# One forest plot per analysis, plus a combined multi-page PDF
# =============================================================================

analyses <- sort(unique(forest_data$analysis))

slugify <- function(x) tolower(gsub("[^A-Za-z0-9]+", "_", x))

for (analysis_name in analyses) {
  analysis_data <- forest_data %>% filter(analysis == analysis_name)
  pooled_data   <- pooled_forest_data %>% filter(analysis == analysis_name)
  plot <- make_forest_plot(analysis_data, analysis_name, pooled_data)

  # Grid dimensions: columns = model specifications, rows = exposure terms. The
  # +1 row leaves space for the pooled diamond beneath the cohorts.
  n_col   <- n_distinct(analysis_data$model_spec)
  n_row   <- length(term_order)
  n_sites <- n_distinct(analysis_data$site) + 1

  ggsave(
    file.path(forest_dir, paste0("forest_", slugify(analysis_name), ".pdf")),
    plot,
    width  = 2.6 * n_col + 2,
    height = (0.35 * n_sites + 0.9) * n_row + 1.5,
    limitsize = FALSE
  )
  message("Forest plot saved: ", analysis_name, " (", n_row, " exposures x ",
          n_col, " models, ", n_sites - 1, " cohorts + pooled)")
}

pdf(file.path(forest_dir, "forest_all_analyses.pdf"), width = 12, height = 9)
for (analysis_name in analyses) {
  print(make_forest_plot(forest_data %>% filter(analysis == analysis_name),
                         analysis_name,
                         pooled_forest_data %>% filter(analysis == analysis_name)))
}
invisible(dev.off())

message("Combined forest plots saved to: ",
        file.path(forest_dir, "forest_all_analyses.pdf"))

# =============================================================================
# Key-exposure summary: outcomes in columns, PBW/PFVC and PFVC in rows
# =============================================================================
# A compact "money figure": for the two headline scaling exposures — PBW/PFVC and
# PFVC — show the effect across every outcome at a glance. Each exposure is taken
# from its canonical adjusted model (PBW/PFVC from the VT/PBW + PBW/PFVC model;
# PFVC from the VT/PBW + PFVC model), so each cell is that exposure's coefficient
# adjusted for delivered dose. Rows = exposure, columns = outcome, with per-site
# points and the random-effects pooled diamond. Ratio outcomes (mortality OR,
# survival HR, VFD SHR) and linear mechanics outcomes (Beta) use different null
# values and scales, so they are drawn as two stacked blocks. The normalized
# Ers x / MP / outcomes are omitted: they share a predicted-size term with these
# exposures, so the relationship would be partly mechanical rather than clinical.
key_exposures <- tibble(
  term       = c("pbwpfvc", "pfvc"),
  term_label = c("PBW/PFVC", "PFVC"),
  src_model  = c("VT/PBW + PBW/PFVC", "VT/PBW + PFVC")
)
ratio_outcomes  <- c("Mortality", "Survival", "28-day VFDs")
linear_outcomes <- c("Compliance", "Elastance", "Static DP", "Mechanical power")

key_site <- all_results %>%
  inner_join(key_exposures, by = c("term", "model_spec" = "src_model")) %>%
  transmute(analysis, term_label, estimate_type, group = site,
            estimate, conf_low, conf_high, kind = "Site")
key_pooled <- pooled_coefs %>%
  inner_join(key_exposures, by = c("term", "model_spec" = "src_model")) %>%
  transmute(analysis, term_label, estimate_type, group = "Pooled (RE)",
            estimate = pooled, conf_low = pooled_lo, conf_high = pooled_hi, kind = "Pooled")

make_key_block <- function(outcomes, x_label, is_ratio) {
  null_value <- if (is_ratio) 1 else 0
  site_b   <- key_site   %>% filter(analysis %in% outcomes)
  pooled_b <- key_pooled %>% filter(analysis %in% outcomes)
  if (nrow(site_b) + nrow(pooled_b) == 0) return(NULL)

  # "Pooled (RE)" is the first y level → numeric position 1 (foot of each panel).
  group_levels <- c("Pooled (RE)", sort(setdiff(unique(site_b$group), "Pooled (RE)")))
  set_factors <- function(d) d %>% mutate(
    group      = factor(group, levels = group_levels),
    analysis   = factor(analysis, levels = outcomes, labels = relabel_outcome(outcomes)),
    term_label = factor(term_label, levels = c("PBW/PFVC", "PFVC")))
  site_b   <- set_factors(site_b)
  pooled_b <- set_factors(pooled_b)

  # Pooled estimate as the classic forest-plot diamond: a horizontal diamond
  # spanning the 95% CI, centred on the pooled estimate, on the pooled row (y = 1).
  # One polygon per facet cell (analysis x term_label).
  yh <- 0.34
  cell_id <- function(d) interaction(d$analysis, d$term_label, drop = TRUE)
  diamonds <- bind_rows(
    transmute(pooled_b, analysis, term_label, cell = cell_id(pooled_b), x = conf_low,  y = 1),
    transmute(pooled_b, analysis, term_label, cell = cell_id(pooled_b), x = estimate,  y = 1 + yh),
    transmute(pooled_b, analysis, term_label, cell = cell_id(pooled_b), x = conf_high, y = 1),
    transmute(pooled_b, analysis, term_label, cell = cell_id(pooled_b), x = estimate,  y = 1 - yh)
  )

  pal <- c(site_colors, `Pooled (RE)` = "#000000")
  p <- ggplot(site_b, aes(x = estimate, y = group, color = group)) +
    geom_vline(xintercept = null_value, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.2) +
    geom_point(size = 2.3) +
    geom_polygon(data = diamonds, aes(x = x, y = y, group = cell),
                 inherit.aes = FALSE, fill = "black", color = "black") +
    facet_grid(term_label ~ analysis, scales = "free_x") +
    scale_color_manual(values = pal, guide = "none") +
    scale_y_discrete(limits = group_levels) +
    labs(x = x_label, y = "Cohort") +
    theme_minimal(base_size = 11) +
    theme(strip.text.x = element_text(face = "bold"),
          strip.text.y = element_text(face = "bold", angle = 0),
          panel.spacing = unit(0.6, "lines"),
          panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.5))
  if (is_ratio) p <- p + scale_x_log10()
  p
}

key_blocks <- list(
  make_key_block(ratio_outcomes,  "OR / HR (95% CI), log scale", TRUE),
  make_key_block(linear_outcomes, "Beta (95% CI)", FALSE)
)
key_blocks <- key_blocks[!map_lgl(key_blocks, is.null)]
if (length(key_blocks) > 0) {
  key_fig <- wrap_plots(key_blocks, ncol = 1) +
    plot_annotation(
      title = "Key scaling exposures across outcomes: PBW/PFVC and PFVC",
      subtitle = paste0("Each exposure from its dose-adjusted model; per-site points + ",
                        "black random-effects pooled diamond; dashed line = null"),
      theme = theme(plot.title = element_text(face = "bold")))
  ggsave(file.path(forest_dir, "forest_key_exposures.pdf"), key_fig,
         width = 11, height = 7, limitsize = FALSE)
  message("Key-exposure summary figure saved: forest_key_exposures.pdf")
}

# =============================================================================
# Covariate forests for the demographic-bias and PFVC-vs-PBW analyses
# =============================================================================
# These analyses regress an outcome on demographics only, so the comparison of
# interest is the covariate coefficients (not exposures). Plot all non-intercept
# terms across cohorts, one facet per term.
covariate_analyses <- all_results %>%
  filter(str_detect(analysis, "^Demo bias:") | analysis == "PFVC vs PBW")

if (nrow(covariate_analyses) > 0) {
  cov_term_labels <- c(
    pbw                = "PBW",
    age10              = "Age (per 10 yr)",
    sex_categoryFemale = "Female vs male",
    race_categoryOTHER = "Other vs white",
    race_categoryBLACK = "Black vs white",
    height10           = "Height (per 10 cm)",
    sf10               = "SF ratio (per 10)",
    sofa_total         = "SOFA",
    bmi                = "BMI"
  )
  cov_data <- covariate_analyses %>%
    mutate(term_label = recode(term, !!!cov_term_labels))

  # Pooled covariate coefficients for overlay (same term labels as cov_data).
  pooled_cov_data <- pooled_coefs %>%
    filter(term %in% names(cov_term_labels)) %>%
    mutate(term_label = recode(term, !!!cov_term_labels))

  make_covariate_forest <- function(df, analysis_name, pooled_df = NULL) {
    est_type <- df$estimate_type[1]
    is_ratio <- est_type %in% c("OR", "HR")
    null_value <- if (is_ratio) 1 else 0
    combined <- df %>%
      transmute(term_label, group = site, estimate, conf_low, conf_high, kind = "Site")
    if (!is.null(pooled_df) && nrow(pooled_df) > 0) {
      combined <- bind_rows(combined, pooled_df %>%
        transmute(term_label, group = "Pooled (RE)",
                  estimate = pooled, conf_low = pooled_lo, conf_high = pooled_hi,
                  kind = "Pooled"))
    }
    group_levels <- c("Pooled (RE)", setdiff(sort(unique(combined$group)), "Pooled (RE)"))
    combined <- combined %>% mutate(group = factor(group, levels = group_levels))
    pal <- c(site_colors, `Pooled (RE)` = "#000000")
    p <- ggplot(combined, aes(x = estimate, y = group, color = group)) +
      geom_vline(xintercept = null_value, linetype = "dashed", color = "grey50") +
      geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.25) +
      geom_point(aes(size = kind, shape = kind)) +
      facet_wrap(~ term_label, scales = "free_x") +
      scale_color_manual(values = pal, guide = "none") +
      scale_shape_manual(values = c(Site = 16, Pooled = 18), guide = "none") +
      scale_size_manual(values = c(Site = 2.3, Pooled = 3.4), guide = "none") +
      labs(title = paste0("Cross-cohort: ", analysis_name),
           subtitle = paste0(est_type, " (95% CI); black diamond = random-effects ",
                             "pooled estimate; dashed line = null (", null_value, ")"),
           x = paste0(est_type, " (95% CI)"), y = "Cohort") +
      theme_minimal(base_size = 11) +
      theme(strip.text = element_text(face = "bold"),
            panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.5))
    if (is_ratio) p <- p + scale_x_log10()
    p
  }

  for (analysis_name in sort(unique(cov_data$analysis))) {
    df <- cov_data %>% filter(analysis == analysis_name)
    pooled_df <- pooled_cov_data %>% filter(analysis == analysis_name)
    n_terms <- n_distinct(df$term_label)
    n_sites <- n_distinct(df$site) + 1
    p <- make_covariate_forest(df, analysis_name, pooled_df)
    ggsave(
      file.path(forest_dir, paste0("covforest_", slugify(analysis_name), ".pdf")),
      p,
      width  = 3 * min(3, n_terms) + 1,
      height = (0.4 * n_sites + 1) * ceiling(n_terms / 3) + 1.5,
      limitsize = FALSE
    )
  }
  message("Covariate forests saved for ", n_distinct(cov_data$analysis), " analyses")
}

# =============================================================================
# Pooled robustness: random-effects meta-analysis + meta-analytic E-values
# =============================================================================
# For each ratio-scale exposure estimate, pool the per-SD estimates across sites
# with a random-effects (REML) meta-analysis, then take the E-value of the POOLED
# estimate. Inputs are the per-site evalues_<site>.csv tables written by script 04
# (section 4g2), which carry the per-1-SD linear- and spline-age estimates.
#
# Cross-site heterogeneity (I^2, tau^2, Q p-value) is the load-bearing statistic:
# a confounder that VARIES across sites — e.g. local ventilation practice for the
# VT-containing exposures, or demographic/SES composition for the deterministic
# ones — would induce heterogeneity, so consistency across heterogeneous sites
# argues against confounding by anything site-specific. The pooled CI is tighter
# than any single site, so the pooled CI E-value is the fairer consortium-level
# robustness number. Replication CANNOT defend against a confounder shared by
# construction at every site (e.g. the linear-age specification); the spline-age
# estimate, also pooled here, is the complementary check for that channel.

evalue_files <- Sys.glob(here("output", "*_output", "final", "evalues_*.csv"))
if (length(evalue_files) == 0) {
  message("No per-site evalues_*.csv found; skipping pooled meta-analytic E-values. ",
          "Run script 04 (section 4g2) for each cohort first.")
} else {
  per_site_ev <- map_dfr(evalue_files, read_csv, show_col_types = FALSE) %>%
    mutate(model_spec = recode(model_spec, "PBW/PFVC + VT/PBW" = "VT/PBW + PBW/PFVC"))

  required_cols <- c("site", "analysis", "model_spec", "term", "term_label",
                     "estimate_type", "est_linage", "lo_linage", "hi_linage",
                     "est_splineage", "lo_splineage", "hi_splineage")
  missing_cols <- setdiff(required_cols, names(per_site_ev))
  if (length(missing_cols) > 0) {
    stop("evalues_*.csv is missing columns: ", paste(missing_cols, collapse = ", "),
         ". Re-run script 04 to regenerate them with the per-SD linear/spline schema.")
  }

  # Point and CI E-values for one ratio estimate (same common-outcome conversions
  # as script 04: rare = FALSE). CI E-value is the bound nearest the null (1 if the
  # CI crosses it).
  evalue_for <- function(est, lo, hi, type) {
    ev_mat <- switch(type,
      OR = EValue::evalues.OR(est, lo, hi, rare = FALSE),
      HR = EValue::evalues.HR(est, lo, hi, rare = FALSE))
    ev    <- ev_mat["E-values", ]
    ci_ev <- ev[c("lower", "upper")]
    ci_ev <- ci_ev[!is.na(ci_ev)]
    list(point = unname(ev[["point"]]),
         ci    = if (length(ci_ev)) unname(ci_ev[1]) else NA_real_)
  }

  # Inverse-variance random-effects (REML) pool of one per-SD ratio column across
  # sites. The log-scale SE is recovered from each site's 95% CI. With a single
  # site the "pool" is that site's estimate and heterogeneity is undefined.
  z975 <- qnorm(0.975)
  pool_ratio <- function(est, lo, hi) {
    yi   <- log(est)
    sei  <- (log(hi) - log(lo)) / (2 * z975)
    keep <- is.finite(yi) & is.finite(sei) & sei > 0
    yi   <- yi[keep]; sei <- sei[keep]; k <- length(yi)
    if (k == 0) return(NULL)
    if (k == 1) {
      return(tibble(k = 1L, est = exp(yi), lo = exp(yi - z975 * sei),
                    hi = exp(yi + z975 * sei),
                    I2 = NA_real_, tau2 = NA_real_, QEp = NA_real_))
    }
    fit <- metafor::rma(yi = yi, sei = sei, method = "REML")
    tibble(k = k, est = exp(as.numeric(fit$beta)), lo = exp(fit$ci.lb),
           hi = exp(fit$ci.ub), I2 = fit$I2, tau2 = fit$tau2, QEp = fit$QEp)
  }

  pooled_ev <- per_site_ev %>%
    group_by(analysis, model_spec, term, term_label, estimate_type) %>%
    group_modify(function(df, key) {
      lin <- pool_ratio(df$est_linage, df$lo_linage, df$hi_linage)
      if (is.null(lin)) return(tibble())
      spl <- pool_ratio(df$est_splineage, df$lo_splineage, df$hi_splineage)
      ev  <- evalue_for(lin$est, lin$lo, lin$hi, key$estimate_type)
      tibble(
        k_sites             = lin$k,
        pooled_lin_est      = lin$est, pooled_lin_lo = lin$lo, pooled_lin_hi = lin$hi,
        I2                  = lin$I2,  tau2 = lin$tau2, het_p = lin$QEp,
        pooled_spline_est   = if (is.null(spl)) NA_real_ else spl$est,
        pooled_spline_lo    = if (is.null(spl)) NA_real_ else spl$lo,
        pooled_spline_hi    = if (is.null(spl)) NA_real_ else spl$hi,
        pooled_evalue_point = ev$point,
        pooled_evalue_ci    = ev$ci
      )
    }) %>%
    ungroup()

  write_csv(pooled_ev, file.path(cross_dir, "pooled_evalues.csv"))

  # Ordered, formatted table grouped by analysis.
  pooled_tbl <- pooled_ev %>%
    mutate(
      analysis   = factor(analysis, levels = c(
        "Survival (60-day)", "Mortality (in-hospital)", "28-day VFDs (liberation)")),
      term_label = factor(term_label, levels = term_order),
      model_spec = factor(model_spec, levels = unique(c(model_order, model_spec)))
    ) %>%
    arrange(analysis, model_spec, term_label) %>%
    transmute(
      analysis = as.character(analysis),
      Model = model_spec, Exposure = term_label, Type = estimate_type,
      Sites = k_sites,
      `Pooled (linear age)` = sprintf("%.2f (%.2f, %.2f)",
                                      pooled_lin_est, pooled_lin_lo, pooled_lin_hi),
      `Pooled (spline age)` = ifelse(is.na(pooled_spline_est), "—",
                                     sprintf("%.2f (%.2f, %.2f)",
                                             pooled_spline_est, pooled_spline_lo, pooled_spline_hi)),
      `I-squared (%)` = ifelse(is.na(I2),   "—", sprintf("%.0f", I2)),
      `tau-squared`   = ifelse(is.na(tau2), "—", sprintf("%.3f", tau2)),
      `Het. p`        = ifelse(is.na(het_p), "—", formatC(het_p, format = "g", digits = 2)),
      `Pooled E-value (est)` = sprintf("%.2f", pooled_evalue_point),
      `Pooled E-value (CI)`  = ifelse(is.na(pooled_evalue_ci), "—",
                                      sprintf("%.2f", pooled_evalue_ci))
    )

  pooled_gt <- pooled_tbl %>%
    gt::gt(groupname_col = "analysis") %>%
    gt::tab_header(
      title = "Pooled robustness: random-effects meta-analysis and meta-analytic E-values",
      subtitle = paste0(
        "Per-1-SD exposure contrasts pooled across ", n_distinct(per_site_ev$site),
        " cohort(s) (REML). Low I-squared / tau-squared across heterogeneous sites ",
        "argues against site-specific confounding; the pooled E-value is the minimum ",
        "risk-ratio association an unmeasured confounder would need with both the ",
        "exposure and the outcome to explain the pooled estimate.")
    ) %>%
    gt::cols_align("left", columns = c(Model, Exposure)) %>%
    gt::sub_missing(missing_text = "—")

  gt::gtsave(pooled_gt, file.path(cross_dir, "table_pooled_evalues.html"))
  gt::gtsave(pooled_gt, file.path(cross_dir, "table_pooled_evalues.pdf"))

  # Forest of pooled per-SD estimates (linear age) with per-site points overlaid
  # and the pooled E-value annotated. One panel per analysis; ratio scale (log x).
  forest_df <- pooled_ev %>%
    mutate(
      row_label = paste0(model_spec, " · ", term_label),
      ev_label  = sprintf("E %.2f / %.2f",
                          pooled_evalue_point,
                          ifelse(is.na(pooled_evalue_ci), 1, pooled_evalue_ci))
    )
  site_pts <- per_site_ev %>%
    left_join(distinct(forest_df, analysis, model_spec, term, row_label),
              by = c("analysis", "model_spec", "term")) %>%
    filter(!is.na(row_label), is.finite(est_linage))

  pooled_forest <- ggplot(forest_df, aes(x = pooled_lin_est, y = row_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_point(data = site_pts, aes(x = est_linage, y = row_label),
               inherit.aes = FALSE, color = "grey65", size = 1.4,
               position = position_nudge(y = 0.18)) +
    geom_errorbarh(aes(xmin = pooled_lin_lo, xmax = pooled_lin_hi), height = 0.22) +
    geom_point(size = 2.6, shape = 18) +
    geom_text(aes(label = ev_label), vjust = -0.9, size = 2.6) +
    facet_wrap(~ analysis, scales = "free", ncol = 1) +
    scale_x_log10() +
    labs(
      title = "Pooled per-SD estimates with meta-analytic E-values",
      subtitle = paste0("Diamonds = REML pooled estimate (95% CI); grey points = per-site; ",
                        "labels = pooled E-value (estimate / CI). Dashed line = null (1)."),
      x = "Pooled per-SD estimate (OR/HR, 95% CI)", y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(face = "bold"))

  n_rows_forest <- n_distinct(forest_df$row_label)
  ggsave(file.path(forest_dir, "forest_pooled_evalues.pdf"), pooled_forest,
         width = 9, height = 0.4 * n_rows_forest + 2.5, limitsize = FALSE)

  message("Pooled meta-analytic E-values written for ", nrow(pooled_ev),
          " ratio-scale exposure estimates across ", n_distinct(per_site_ev$site),
          " cohort(s); see ", file.path(cross_dir, "table_pooled_evalues.html"))
}

# =============================================================================
# Pooled CONSORT diagram (sum attrition across cohorts)
# =============================================================================
attrition_files <- Sys.glob(here("output", "*_output", "final", "attrition_log_*.csv"))
if (length(attrition_files) > 0) {
  pooled_attrition <- map_dfr(attrition_files, read_csv, show_col_types = FALSE) %>%
    group_by(step_order, step_label) %>%
    summarise(n_remaining = sum(n_remaining),
              # The first attrition step (starting population) has no exclusion
              # reason. read_csv types an all-blank column as logical, so coerce to
              # character and reduce by hand — this always returns a size-1 string
              # (the reason, or NA) regardless of how each file parsed the column.
              exclusion_reason = {
                reasons <- as.character(exclusion_reason)
                reasons <- reasons[!is.na(reasons)]
                if (length(reasons) == 0L) NA_character_ else reasons[[1L]]
              },
              .groups = "drop") %>%
    arrange(step_order) %>%
    mutate(n_excluded = if_else(step_order == min(step_order), NA_integer_,
                                as.integer(lag(n_remaining) - n_remaining)))
  consort_pooled <- render_consort(pooled_attrition,
                                   title = "Cohort inclusion - pooled across cohorts")
  ggsave(file.path(cross_dir, "consort_diagram_pooled.pdf"),
         consort_pooled, width = 9, height = 11)
  message("Pooled CONSORT diagram saved")
}

# =============================================================================
# Pooled PBW:PFVC-by-demographics figure (sum histograms across cohorts)
# =============================================================================
# Reconstructs violins and binned medians from the FEDERATED histogram exports
# only -- no row-level data. Per-(group, bin) counts are summed across sites for
# a true pooled distribution.
hist_files <- Sys.glob(here("output", "*_output", "final", "dist_histograms_*.csv"))
if (length(hist_files) > 0) {
  pooled_hist <- map_dfr(hist_files, read_csv, show_col_types = FALSE) %>%
    group_by(group_type, group_value, bin_left, bin_right) %>%
    summarise(count = sum(count), .groups = "drop")

  # Quantile of a single group from its binned counts (empirical CDF at centres).
  hist_quantile <- function(h, p) {
    h <- h[order(h$bin_left), ]
    stats::approx(cumsum(h$count) / sum(h$count), (h$bin_left + h$bin_right) / 2,
                  xout = p, ties = "ordered", rule = 2)$y
  }
  # Mirrored violin polygon from binned counts.
  violin_poly <- function(h, halfwidth = 0.42) {
    h <- h[order(h$bin_left), ]
    mid  <- (h$bin_left + h$bin_right) / 2
    dens <- h$count / sum(h$count)
    half <- halfwidth * dens / max(dens)
    tibble(x_off = c(-half, rev(half)), y = c(mid, rev(mid)))
  }
  bin_mid <- function(label) {
    m  <- str_match(label, "\\[([-0-9.]+),\\s*([-0-9.eInf+]+)\\)")
    lo <- as.numeric(m[, 2]); hi <- as.numeric(m[, 3])
    ifelse(is.infinite(hi), lo + 5, (lo + hi) / 2)
  }

  pooled_box <- pooled_hist %>%
    group_by(group_type, group_value) %>%
    summarise(q25 = hist_quantile(pick(everything()), 0.25),
              q50 = hist_quantile(pick(everything()), 0.50),
              q75 = hist_quantile(pick(everything()), 0.75),
              .groups = "drop")

  cat_panel <- function(group, xlab) {
    h <- pooled_hist %>% filter(group_type == group)
    if (nrow(h) == 0) return(NULL)
    lv <- sort(unique(h$group_value)); idx <- setNames(seq_along(lv), lv)
    polys <- h %>% group_by(group_value) %>% group_modify(~ violin_poly(.x)) %>%
      ungroup() %>% mutate(x = idx[group_value] + x_off)
    box <- pooled_box %>% filter(group_type == group) %>% mutate(xc = idx[group_value])
    ggplot() +
      geom_polygon(data = polys, aes(x, y, group = group_value, fill = group_value),
                   alpha = 0.5, colour = NA) +
      geom_segment(data = box, aes(x = xc - 0.12, xend = xc + 0.12, y = q50, yend = q50),
                   linewidth = 0.7) +
      geom_linerange(data = box, aes(x = xc, ymin = q25, ymax = q75), linewidth = 0.7) +
      scale_x_continuous(breaks = idx, labels = names(idx)) +
      scale_fill_manual(values = okabe_ito, guide = "none") +
      labs(x = xlab, y = "PBW:PFVC ratio (kg/L)") + theme_minimal(base_size = 11)
  }
  cont_panel <- function(group, xlab) {
    h <- pooled_hist %>% filter(group_type == group)
    if (nrow(h) == 0) return(NULL)
    box <- pooled_box %>% filter(group_type == group) %>%
      mutate(xmid = bin_mid(group_value)) %>% filter(!is.na(xmid)) %>% arrange(xmid)
    ggplot(box, aes(xmid, q50)) +
      geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2, fill = "#D55E00") +
      geom_line(colour = "#0072B2", linewidth = 0.7) +
      geom_point(colour = "#0072B2", size = 1.8) +
      labs(x = xlab, y = "PBW:PFVC ratio (kg/L)") + theme_minimal(base_size = 11)
  }

  panels <- list(cat_panel("sex", "Sex"), cat_panel("race", "Race"),
                 cont_panel("age_bin", "Age (years)"),
                 cont_panel("height_bin", "Height (cm)"))
  panels <- panels[!map_lgl(panels, is.null)]
  if (length(panels) > 0) {
    pooled_dist_fig <- wrap_plots(panels, ncol = 2) +
      plot_annotation(title = "PBW:PFVC by demographics - pooled across cohorts",
                      theme = theme(plot.title = element_text(face = "bold")))
    ggsave(file.path(cross_dir, "distribution_pbwpfvc_pooled.pdf"),
           pooled_dist_fig, width = 11, height = 9)
    message("Pooled PBW:PFVC distribution figure saved")
  }
}

message("Script 05 complete. Aggregated results and figures in: ", cross_dir)
