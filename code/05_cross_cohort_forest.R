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
  pbwpfvc    = "PBW/PFVC",
  ers_pbw_z  = "Ers x PBW",   # z-scaled in the mortality models (per-SD OR)
  ers_pfvc_z = "Ers x PFVC"
)

# Row order (exposures) and column order (models, simplest to most complex).
term_order  <- c("VT/PFVC", "VT/PBW", "PFVC", "PBW/PFVC", "Ers x PBW", "Ers x PFVC")
model_order <- c("VT/PFVC", "VT/PBW", "VT/PFVC + VT/PBW", "VT/PBW + PFVC",
                 "VT/PBW + PBW/PFVC", "PBW/PFVC + VT/PBW", "Ers x PBW", "Ers x PFVC")

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

# Build a forest plot for one analysis (a single estimate_type by construction).
make_forest_plot <- function(analysis_data, analysis_name) {
  estimate_type <- analysis_data$estimate_type[1]
  is_ratio <- estimate_type %in% c("OR", "HR")
  null_value <- if (is_ratio) 1 else 0

  plot <- ggplot(analysis_data,
                 aes(x = estimate, y = site, color = site)) +
    geom_vline(xintercept = null_value, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.25) +
    geom_point(size = 2.5) +
    # Grid: one row per exposure, one column per model specification. Empty cells
    # mark exposures a given model does not include. scales = "free_x" lets each
    # model column take its own estimate range; drop = FALSE keeps the full grid.
    facet_grid(term_label ~ model_spec, scales = "free_x", drop = FALSE) +
    scale_color_manual(values = site_colors, guide = "none") +
    labs(
      title = paste0("Cross-cohort results: ", analysis_name),
      subtitle = paste0(estimate_type, " (95% CI); dashed line = null effect (",
                        null_value, ")"),
      x = paste0(estimate_type, " (95% CI)"),
      y = "Cohort"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text.x = element_text(face = "bold"),
      strip.text.y = element_text(face = "bold", angle = 0),
      panel.spacing = unit(0.6, "lines")
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
  plot <- make_forest_plot(analysis_data, analysis_name)

  # Grid dimensions: columns = model specifications, rows = exposure terms.
  n_col   <- n_distinct(analysis_data$model_spec)
  n_row   <- length(term_order)
  n_sites <- n_distinct(analysis_data$site)

  ggsave(
    file.path(forest_dir, paste0("forest_", slugify(analysis_name), ".pdf")),
    plot,
    width  = 2.6 * n_col + 2,
    height = (0.35 * n_sites + 0.9) * n_row + 1.5,
    limitsize = FALSE
  )
  message("Forest plot saved: ", analysis_name, " (", n_row, " exposures x ",
          n_col, " models, ", n_sites, " cohorts)")
}

pdf(file.path(forest_dir, "forest_all_analyses.pdf"), width = 12, height = 9)
for (analysis_name in analyses) {
  print(make_forest_plot(forest_data %>% filter(analysis == analysis_name),
                         analysis_name))
}
invisible(dev.off())

message("Combined forest plots saved to: ",
        file.path(forest_dir, "forest_all_analyses.pdf"))

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

  make_covariate_forest <- function(df, analysis_name) {
    est_type <- df$estimate_type[1]
    is_ratio <- est_type %in% c("OR", "HR")
    null_value <- if (is_ratio) 1 else 0
    p <- ggplot(df, aes(x = estimate, y = site, color = site)) +
      geom_vline(xintercept = null_value, linetype = "dashed", color = "grey50") +
      geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.25) +
      geom_point(size = 2.5) +
      facet_wrap(~ term_label, scales = "free_x") +
      scale_color_manual(values = site_colors, guide = "none") +
      labs(title = paste0("Cross-cohort: ", analysis_name),
           subtitle = paste0(est_type, " (95% CI); dashed line = null (", null_value, ")"),
           x = paste0(est_type, " (95% CI)"), y = "Cohort") +
      theme_minimal(base_size = 11) +
      theme(strip.text = element_text(face = "bold"))
    if (is_ratio) p <- p + scale_x_log10()
    p
  }

  for (analysis_name in sort(unique(cov_data$analysis))) {
    df <- cov_data %>% filter(analysis == analysis_name)
    n_terms <- n_distinct(df$term_label)
    n_sites <- n_distinct(df$site)
    p <- make_covariate_forest(df, analysis_name)
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
# Pooled CONSORT diagram (sum attrition across cohorts)
# =============================================================================
attrition_files <- Sys.glob(here("output", "*_output", "final", "attrition_log_*.csv"))
if (length(attrition_files) > 0) {
  pooled_attrition <- map_dfr(attrition_files, read_csv, show_col_types = FALSE) %>%
    group_by(step_order, step_label) %>%
    summarise(n_remaining = sum(n_remaining),
              exclusion_reason = first(na.omit(exclusion_reason)),
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
