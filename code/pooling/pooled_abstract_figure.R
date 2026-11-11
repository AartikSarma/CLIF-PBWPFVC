# =============================================================================
# pooled_abstract_figure.R: the two-part figure for the abstract
# =============================================================================
# Run centrally by the coordinator after every site has returned its final/ folder,
# like pooled_estimates.R, and from the same results root (one subfolder per site).
#
#   Panel A  the PBW:PFVC ratio by sex, race, age and height, pooled across cohorts
#            by summing each site's binned counts (dist_histograms_<site>.csv, script
#            04). No row-level data: violins and medians are rebuilt from the bins.
#   Panel B  the PBW:PFVC coefficient for in-hospital mortality (logistic, OR) and
#            respiratory-system compliance (linear, mL/cmH2O), per cohort and pooled
#            by random-effects (REML) meta-analysis. Both come from the
#            "VT/PBW + PBW/PFVC" model in regression_results_long_<site>.csv (script
#            04), demographic-UNADJUSTED: vtpbw + pbwpfvc + SOFA + SF/10, with BMI
#            added for compliance. Set FIGURE_ADJUSTMENT to "adjusted" to draw the
#            model that also carries age, sex and race.
#
# The script prints the pooled estimate and 95% CI for both outcomes at both
# adjustment levels, so the abstract can quote them.
#
# Usage: PBWPFVC_RESULTS_ROOT=/path/to/results uvr run code/pooling/pooled_abstract_figure.R
#   (default root: results/; output: <root>/All sites/abstract_figure.{pdf,png})
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
source("utils/site_anonymization.R")

FIGURE_ADJUSTMENT <- "unadjusted"
SOURCE_MODEL      <- "VT/PBW + PBW/PFVC"
# The figure is drawn at its printed size (7.5 in wide, a letter page inside its
# margins), so every size below is the size on the page, in points.
FIGURE_WIDTH_IN       <- 7.5
FIGURE_HEIGHT_IN      <- 5.4
ANNOTATION_PT         <- 8    # in-plot text: group differences, R-squared
AXIS_LABEL_PT         <- 9   # tick labels
AXIS_TITLE_PT         <- 9
PLOT_TITLE_PT         <- 10   # outcome titles over the forest plots
PANEL_TAG_PT          <- 16   # the A and B tags

results_root <- Sys.getenv("PBWPFVC_RESULTS_ROOT", unset = "results")
if (!dir.exists(results_root)) stop("Results root not found: ", normalizePath(results_root, mustWork = FALSE))
output_dir <- file.path(results_root, "All sites")
dir.create(output_dir, showWarnings = FALSE)

site_dirs <- list.dirs(results_root, recursive = FALSE, full.names = TRUE)
site_dirs <- site_dirs[basename(site_dirs) != "All sites" & !startsWith(basename(site_dirs), ".")]
site_aliases <- build_site_aliases(site_dirs)
print_site_alias_key(site_aliases)

# A site's tables sit flat in its folder, in final/, or in a cross_sectional/ block
# folder. Exactly one copy per site is expected; a second is a stale return.
read_site_tables <- function(pattern) {
  map_dfr(site_dirs, function(site_dir) {
    files <- unlist(map(c("", "final", "cross_sectional", "final/cross_sectional"),
                        ~ Sys.glob(file.path(site_dir, .x, pattern))))
    if (length(files) == 0) stop("No ", pattern, " in ", site_dir)
    if (length(files) > 1) stop(site_dir, " holds ", length(files), " copies of ", pattern,
                                ":\n  ", paste(files, collapse = "\n  "))
    read_csv(files, show_col_types = FALSE)
  })
}

OKABE_ITO <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#000000")
theme_abstract <- theme_minimal(base_size = AXIS_LABEL_PT) +
  theme(axis.title = element_text(size = AXIS_TITLE_PT),
        axis.text  = element_text(size = AXIS_LABEL_PT, colour = "grey20"),
        plot.title = element_text(size = PLOT_TITLE_PT, face = "bold"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(2, 4, 2, 2),
        aspect.ratio = 1)   # every plot square
ANNOTATION_TEXT_SIZE <- ANNOTATION_PT / .pt   # geom_text sizes are in mm
RATIO_LABEL <- "PBW:PFVC (kg/L)"
# One whole-number y scale shared by all four panel-A plots, wide enough for the
# violin tails, so the groups can be compared across plots. The headroom above 30
# holds the annotations.
ratio_y_scale <- scale_y_continuous(limits = c(12, 33), breaks = seq(12, 30, by = 3))

# =============================================================================
# Panel A: the ratio by demographics, from summed histogram bins
# =============================================================================
# A site bin of 1-9 patients arrives blanked (NA count) and is left out.
site_histograms <- read_site_tables("dist_histograms_*.csv") %>% filter(!is.na(count), count > 0)
pooled_histograms <- site_histograms %>%
  group_by(group_type, group_value, bin_left, bin_right) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  filter(count > 0) %>%
  group_by(group_type, group_value) %>% filter(n() >= 2) %>% ungroup()

# Quantile of one group from its binned counts (empirical CDF at the bin centres).
binned_quantile <- function(bins, probability) {
  bins <- bins[order(bins$bin_left), ]
  stats::approx(cumsum(bins$count) / sum(bins$count), (bins$bin_left + bins$bin_right) / 2,
                xout = probability, ties = "ordered", rule = 2)$y
}
# Mirrored violin outline from binned counts.
violin_outline <- function(bins, half_width = 0.42) {
  bins <- bins[order(bins$bin_left), ]
  bin_centre <- (bins$bin_left + bins$bin_right) / 2
  half <- half_width * bins$count / max(bins$count)
  tibble(x_offset = c(-half, rev(half)), y = c(bin_centre, rev(bin_centre)))
}
# Midpoint of an "[lo, hi)" age or height bin label; the open top bin gets lo + 5.
bin_label_midpoint <- function(label) {
  bounds <- str_match(label, "\\[([-0-9.]+),\\s*([-0-9.eInf+]+)\\)")
  lower <- as.numeric(bounds[, 2]); upper <- as.numeric(bounds[, 3])
  ifelse(is.infinite(upper), lower + 5, (lower + upper) / 2)
}

group_quartiles <- pooled_histograms %>%
  group_by(group_type, group_value) %>%
  summarise(q25 = binned_quantile(pick(everything()), 0.25),
            q50 = binned_quantile(pick(everything()), 0.50),
            q75 = binned_quantile(pick(everything()), 0.75), .groups = "drop")

format_p <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))

# Group differences for the categorical plots. Each site's binned counts give the
# group mean and variance of the ratio (bin centres); the difference in means from
# the reference group (male, White) and its SE are taken within each site, then
# pooled across sites by random-effects (REML) meta-analysis, as in panel B.
group_differences <- function(group, reference_level) {
  moments <- site_histograms %>%
    filter(group_type == group) %>%
    mutate(bin_centre = (bin_left + bin_right) / 2) %>%
    group_by(site, group_value) %>%
    summarise(patients = sum(count),
              mean_ratio = sum(count * bin_centre) / patients,
              variance_ratio = sum(count * (bin_centre - mean_ratio)^2) / (patients - 1),
              .groups = "drop")
  reference <- moments %>% filter(group_value == reference_level) %>%
    select(site, ref_patients = patients, ref_mean = mean_ratio, ref_variance = variance_ratio)
  moments %>%
    filter(group_value != reference_level) %>%
    inner_join(reference, by = "site") %>%
    mutate(difference = mean_ratio - ref_mean,
           difference_se = sqrt(variance_ratio / patients + ref_variance / ref_patients)) %>%
    group_by(group_value) %>%
    group_modify(function(rows, key) {
      fit <- metafor::rma(yi = rows$difference, sei = rows$difference_se, method = "REML")
      tibble(sites = fit$k, difference = as.numeric(fit$beta),
             difference_low = fit$ci.lb, difference_high = fit$ci.ub, p_value = fit$pval)
    }) %>% ungroup() %>%
    mutate(group_type = group, reference = reference_level)
}
categorical_tests <- bind_rows(group_differences("sex", "Male"),
                               group_differences("race", "WHITE"))

# Fit of the ratio on age (quadratic) and height (linear) across all pooled patients.
# Each patient is represented by the centre of their ratio bin and of their age or
# height bin, so the R-squared is that of the binned data (binning age and height
# loses within-bin spread, which lowers R-squared somewhat). The F-test uses the
# pooled patient count.
binned_fit <- function(group, polynomial_degree) {
  bins <- pooled_histograms %>% filter(group_type == group) %>%
    mutate(predictor = bin_label_midpoint(group_value),
           ratio = (bin_left + bin_right) / 2)
  fit <- lm(ratio ~ poly(predictor, polynomial_degree, raw = TRUE), data = bins, weights = count)
  ss_total    <- sum(bins$count * (bins$ratio - weighted.mean(bins$ratio, bins$count))^2)
  ss_residual <- sum(bins$count * residuals(fit)^2)
  r_squared <- 1 - ss_residual / ss_total
  patients  <- sum(bins$count)
  f_statistic <- (r_squared / polynomial_degree) / ((1 - r_squared) / (patients - polynomial_degree - 1))
  tibble(group_type = group, model = if (polynomial_degree == 2) "Quadratic" else "Linear",
         patients = patients, r_squared = r_squared,
         p_value = pf(f_statistic, polynomial_degree, patients - polynomial_degree - 1, lower.tail = FALSE))
}
continuous_fits <- bind_rows(binned_fit("age_bin", 2), binned_fit("height_bin", 1))

cat("\nPBW:PFVC differences in mean ratio (kg/L) vs reference, random effects across sites\n")
categorical_tests %>%
  transmute(line = sprintf("  %-6s vs %-6s %+.2f (95%% CI %.2f to %.2f), %s; %d sites",
                           str_to_title(group_value), str_to_title(reference), difference,
                           difference_low, difference_high, format_p(p_value), sites)) %>%
  pull(line) %>% cat(sep = "\n")
cat("PBW:PFVC fits on binned pooled data\n")
continuous_fits %>%
  transmute(line = sprintf("  %-10s %-9s R2 = %.3f, %s; %s patients", group_type, model,
                           r_squared, format_p(p_value), format(patients, big.mark = ","))) %>%
  pull(line) %>% cat(sep = "\n")
cat("\n")

categorical_panel <- function(group, axis_title) {
  bins <- pooled_histograms %>% filter(group_type == group)
  levels_present <- sort(unique(bins$group_value))
  position <- setNames(seq_along(levels_present), levels_present)
  outlines <- bins %>% group_by(group_value) %>% group_modify(~ violin_outline(.x)) %>%
    ungroup() %>% mutate(x = position[group_value] + x_offset)
  quartiles <- group_quartiles %>% filter(group_type == group) %>%
    mutate(x_centre = position[group_value])
  # Above each non-reference group: its pooled mean difference from the reference
  # group (kg/L) and p-value; the reference group is marked as such.
  labels <- categorical_tests %>% filter(group_type == group) %>%
    transmute(group_value, label = sprintf("%+.1f\n%s", difference, format_p(p_value)))
  reference_level <- unique(categorical_tests$reference[categorical_tests$group_type == group])
  labels <- bind_rows(labels, tibble(group_value = reference_level, label = "Reference")) %>%
    mutate(x_centre = position[group_value])
  ggplot() +
    geom_polygon(data = outlines, aes(x, y, group = group_value, fill = group_value),
                 alpha = 0.6, colour = NA) +
    geom_linerange(data = quartiles, aes(x = x_centre, ymin = q25, ymax = q75), linewidth = 0.45) +
    geom_segment(data = quartiles, aes(x = x_centre - 0.14, xend = x_centre + 0.14,
                                       y = q50, yend = q50), linewidth = 0.45) +
    geom_text(data = labels, aes(x = x_centre, y = 31.4, label = label),
              size = ANNOTATION_TEXT_SIZE, lineheight = 0.9) +
    scale_x_continuous(breaks = position, labels = str_to_title(names(position))) +
    scale_fill_manual(values = OKABE_ITO, guide = "none") +
    ratio_y_scale +
    labs(x = axis_title, y = RATIO_LABEL) + theme_abstract +
    theme(panel.grid.major.x = element_blank())
}
continuous_panel <- function(group, axis_title) {
  fit <- continuous_fits %>% filter(group_type == group)
  quartiles <- group_quartiles %>% filter(group_type == group) %>%
    mutate(bin_midpoint = bin_label_midpoint(group_value)) %>%
    filter(!is.na(bin_midpoint)) %>% arrange(bin_midpoint)
  ggplot(quartiles, aes(bin_midpoint, q50)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#D55E00", alpha = 0.2) +
    geom_line(colour = "#0072B2", linewidth = 0.5) +
    geom_point(colour = "#0072B2", size = 1.2) +
    annotate("text", x = min(quartiles$bin_midpoint), y = 31.4, hjust = 0, lineheight = 0.9, size = ANNOTATION_TEXT_SIZE,
             label = with(fit, sprintf("%s fit: R\u00b2 = %.2f\n%s", model, r_squared, format_p(p_value)))) +
    ratio_y_scale +
    labs(x = axis_title, y = RATIO_LABEL) + theme_abstract
}

# The four plots share one y scale, so only the left column carries the y title.
no_y_title <- labs(y = NULL)
panel_a <- (categorical_panel("sex", "Sex") | categorical_panel("race", "Race") + no_y_title) /
  (continuous_panel("age_bin", "Age (years)") | continuous_panel("height_bin", "Height (cm)") + no_y_title)

# =============================================================================
# Panel B: per-cohort and pooled PBW:PFVC coefficients
# =============================================================================
regression_results <- read_site_tables("regression_results_long_*.csv") %>%
  anonymize_site_column(site_aliases$aliases) %>%
  filter(term == "pbwpfvc", model_spec == SOURCE_MODEL,
         analysis %in% c("Mortality", "Compliance"))

# Random-effects (REML) pooling, as in pooled_estimates.R: odds ratios on the log
# scale, back-transformed; compliance betas on their natural scale. The SE is each
# site's Wald std_error.
pool_random_effects <- function(site_rows) {
  is_ratio <- site_rows$estimate_type[1] == "OR"
  effect   <- if (is_ratio) log(site_rows$estimate) else site_rows$estimate
  fit <- metafor::rma(yi = effect, sei = site_rows$std_error, method = "REML")
  back_transform <- if (is_ratio) exp else identity
  tibble(cohorts = fit$k, pooled = back_transform(as.numeric(fit$beta)),
         pooled_low = back_transform(fit$ci.lb), pooled_high = back_transform(fit$ci.ub),
         i_squared = fit$I2, heterogeneity_p = fit$QEp)
}
pooled_estimates <- regression_results %>%
  group_by(analysis, adjustment) %>%
  group_modify(~ pool_random_effects(.x)) %>% ungroup()

cat("\nPooled PBW:PFVC coefficient (per 1 kg/L), ", SOURCE_MODEL, " model, random effects (REML)\n", sep = "")
pooled_estimates %>%
  arrange(desc(analysis), desc(adjustment)) %>%
  transmute(line = sprintf("  %-10s %-10s %s %.3f (95%% CI %.3f to %.3f); %d cohorts, I2 = %.0f%%, heterogeneity p = %.2g%s",
                           analysis, adjustment,
                           ifelse(analysis == "Mortality", "OR  ", "beta"),
                           pooled, pooled_low, pooled_high, cohorts, i_squared, heterogeneity_p,
                           ifelse(adjustment == FIGURE_ADJUSTMENT, "   <- in figure", ""))) %>%
  pull(line) %>% cat(sep = "\n")
cat("  (compliance beta in mL/cmH2O per 1 kg/L)\n\n")

site_order  <- site_aliases$table$site_label
site_colors <- setNames(grDevices::colorRampPalette(OKABE_ITO)(length(site_order)), site_order)
row_labels  <- c("Pooled", site_order)   # row 0 = pooled, at the foot

forest_panel <- function(outcome, panel_title, axis_title, is_ratio) {
  site_rows <- regression_results %>%
    filter(analysis == outcome, adjustment == FIGURE_ADJUSTMENT) %>%
    mutate(row = match(site, site_order))
  pooled_row <- pooled_estimates %>% filter(analysis == outcome, adjustment == FIGURE_ADJUSTMENT)
  # Classic forest diamond spanning the pooled 95% CI on row 0. Every y here is
  # numeric so the diamond and the cohort rows share one continuous axis.
  diamond_half_height <- 0.38
  diamond <- tibble(x = c(pooled_row$pooled_low, pooled_row$pooled, pooled_row$pooled_high, pooled_row$pooled),
                    y = c(0, diamond_half_height, 0, -diamond_half_height))
  null_value <- if (is_ratio) 1 else 0
  panel <- ggplot(site_rows, aes(y = row)) +
    geom_vline(xintercept = null_value, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high, colour = site),
                  orientation = "y", width = 0.3, linewidth = 0.4) +
    geom_point(aes(x = estimate, colour = site), size = 1.3) +
    geom_polygon(data = diamond, aes(x, y), fill = "black") +
    scale_colour_manual(values = site_colors, guide = "none") +
    scale_y_continuous(breaks = seq_along(row_labels) - 1, labels = row_labels,
                       expand = expansion(add = 0.6)) +
    labs(title = panel_title, x = axis_title, y = NULL) +
    theme_abstract +
    theme(panel.grid.major.y = element_blank(),
          panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.5),
          plot.title.position = "panel")
  panel <- panel + if (is_ratio) scale_x_log10() else scale_x_continuous(breaks = seq(-3, 3, by = 1))
  panel
}

panel_b <- forest_panel("Mortality", "In-hospital mortality",
                        "Odds ratio\n per 1 kg/L PBW:PFVC (95% CI)", TRUE) /
  forest_panel("Compliance", "Compliance",
               "Compliance difference, mL/cmH2O,\nper 1 kg/L PBW:PFVC (95% CI)", FALSE)

# Each half is wrapped as one element so the tags label the halves, not the plots.
abstract_figure <- (wrap_elements(full = panel_a) | wrap_elements(full = panel_b)) +
  plot_layout(widths = c(1.5, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = PANEL_TAG_PT, face = "bold"),
        plot.tag.location = "margin")

ggsave(file.path(output_dir, "abstract_figure.pdf"), abstract_figure, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN)
ggsave(file.path(output_dir, "abstract_figure.png"), abstract_figure, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = 300, bg = "white")
message("Abstract figure saved: ", file.path(output_dir, "abstract_figure.{pdf,png}"))
