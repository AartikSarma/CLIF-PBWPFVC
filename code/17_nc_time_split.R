# =============================================================================
# Script 17: is the height (lung-size) mortality gradient ventilation-specific once
# hazards may vary over time? Time-split Cox models in the negative-control cohorts
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Standalone; NOT in 00_run_pipeline.R. Reads the analytic cohort (analysis_cross_sectional)
# and script 03's negative-control cohorts (analysis_negative_control). Writes aggregates only.
#
# THE QUESTION. Script 04 (4j) fits one 60-day Cox HR per cohort. At MIMIC, height's HR was
# 0.86 in the ventilated analytic cohort and 0.91 in the never-ventilated cohort, which reads
# as "most of height's effect does not need a ventilator". That reading assumes proportional
# hazards. The cohorts die on different schedules (at MIMIC 86% of the ventilated cohort's
# 60-day deaths occur in hospital, against 32% of the never-ventilated cohort's), so one
# 60-day HR averages over different time windows in each. If height acts EARLY through the
# ventilator and LATE through something shared (frailty, body size), the single HRs would
# look similar while the mechanism differs. This script lets height's HR vary by window.
#
# WHAT IT FITS, per cohort x exposure x age form (linear, ns(age, 4)), adjusted for sex and race:
#   1. The single 60-day Cox model, with the Schoenfeld-residual test (cox.zph, KM transform)
#      for the exposure and globally: does the exposure's HR change over follow-up?
#   2. A time-split Cox model (survSplit at CUTS, default days 7 and 28): one HR for the
#      exposure in each window (0-7, 7-28, 28-60 days), covariates held proportional; plus
#      the likelihood-ratio test of window-specific against one constant HR.
#   3. Between cohorts, in each window: the ventilated-minus-control difference in log HR.
#      And the difference-in-differences: (early vent - early control) - (late vent - late
#      control), with the within-cohort covariance of the window HRs taken from each model's
#      variance matrix. Cohorts are disjoint patients, so between-cohort terms are independent.
#
# HOW TO READ IT. The framework (ventilation-specific harm of a small lung) predicts:
#   early windows: ventilated HR < 1, not shared by the never-ventilated cohort (a contrast
#                  ratio below 1); late windows: the two cohorts' HRs agree; a Schoenfeld test
#                  flagging the exposure as non-proportional in the ventilated cohort; and a
#                  difference-in-differences below 0 (per SD of height; above 0 for PBW/PFVC,
#                  whose harm runs the other way). If the HR is constant and similar across
#                  windows in both cohorts, height acts mostly outside the ventilator.
#
# TIME ZERO. The never-ventilated cohort is defined over the whole stay, so its only natural
# origin is admission. The analytic cohort's is its index (the first qualifying ventilated,
# hypoxemic timepoint). DEFAULT (--origin index): each cohort from its own entry, which avoids
# immortal time between admission and the analytic index. --origin admission measures the
# analytic cohort from admission too, reproducing script 04's convention.
#
# DESIGN CAVEAT, stated because it bears on the early window. The never-ventilated cohort is
# selected on the future: its patients were never ventilated and never hypoxemic during the
# stay, so patients who deteriorated are excluded. Its early hazard is therefore low by
# construction, and an early height effect there is conditional on not deteriorating.
#
# Synthetic CLIF mortality is unreliable, so on the synthetic site death times are SIMULATED
# for every cohort: plumbing only, never a result.
#
# Outputs (aggregates only; cells with fewer than MIN_CELL deaths are suppressed):
#   final/nc_timesplit_hr_{site}.csv        HR per cohort x exposure x age form x window
#   final/nc_timesplit_contrast_{site}.csv  ventilated-minus-control per window, and the DiD
#   final/nc_ph_test_{site}.csv             Schoenfeld tests of the single 60-day models
#   final/nc_timesplit_counts_{site}.csv    deaths per cohort x window
#   final/nc_timesplit_{site}.pdf
#   A non-default origin or cut set adds a suffix before the site name.
# Usage: Rscript code/17_nc_time_split.R [--site_name NAME] [--output_root DIR]
#          [--origin index|admission] [--cuts 7,28]
# =============================================================================
rm(list = ls())

# --- Command-line arguments (parsed the way code/00_run_pipeline.R parses its own) ------
parse_script_args <- function(args) {
  valued <- c("site_name", "output_root", "origin", "cuts")
  usage <- paste("Usage: Rscript code/17_nc_time_split.R [--site_name NAME] [--output_root DIR]",
                 "[--origin index|admission] [--cuts 7,28]")
  parsed <- list(); i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("--help", "-h")) { message(usage); quit(save = "no", status = 0) }
    if (grepl("^--[a-z_]+=", a)) {                       # --flag=value
      key <- sub("^--([a-z_]+)=.*$", "\\1", a); val <- sub("^--[a-z_]+=", "", a); i <- i + 1L
    } else if (grepl("^--[a-z_]+$", a)) {                # --flag value
      key <- sub("^--", "", a)
      if (!key %in% valued) stop("Unknown option --", key, "\n", usage)
      if (i == length(args) || grepl("^--", args[[i + 1L]]))
        stop("Missing value for --", key, "\n", usage)
      val <- args[[i + 1L]]; i <- i + 2L
    } else stop("Unrecognized argument: ", a, "\n", usage)
    if (!key %in% valued) stop("Unknown option --", key, "\n", usage)
    if (!nzchar(val)) stop("Empty value for --", key, "\n", usage)
    parsed[[key]] <- val
  }
  parsed
}
cli_args <- parse_script_args(commandArgs(trailingOnly = TRUE))

HORIZON_DAYS <- 60
ORIGIN <- if (is.null(cli_args$origin)) "index" else cli_args$origin
if (!ORIGIN %in% c("index", "admission")) stop("--origin must be index or admission; got '", ORIGIN, "'")
CUTS <- if (is.null(cli_args$cuts)) c(7, 28) else suppressWarnings(as.numeric(strsplit(cli_args$cuts, ",")[[1]]))
if (any(!is.finite(CUTS)) || any(CUTS <= 0) || any(CUTS >= HORIZON_DAYS) || is.unsorted(CUTS, strictly = TRUE))
  stop("--cuts must be increasing days strictly between 0 and ", HORIZON_DAYS, " (e.g. 7,28)")
run_suffix <- paste0(if (ORIGIN != "index") "_admission" else "",
                     if (!identical(CUTS, c(7, 28))) paste0("_cuts", paste(CUTS, collapse = "_")) else "")

# Run from the repository root whatever the caller's working directory. --output_root is
# resolved against the CALLER's directory first, before the setwd.
if (!is.null(cli_args$output_root))
  cli_args$output_root <- normalizePath(path.expand(cli_args$output_root), mustWork = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(script_file) == 1) setwd(normalizePath(file.path(dirname(script_file), "..")))
if (!is.null(cli_args$site_name)) Sys.setenv(PBWPFVC_SITE_NAME = cli_args$site_name)

suppressPackageStartupMessages({
  library(tidyverse); library(arrow); library(here); library(splines); library(survival)
})
source("utils/config.R")
site_name    <- config$site_name
is_synthetic <- identical(site_name, "synthetic_clif")
output_root  <- if (is.null(cli_args$output_root)) here("output") else cli_args$output_root
output_dir   <- file.path(output_root, paste0(site_name, "_output"), "intermediate")
final_dir    <- file.path(output_root, paste0(site_name, "_output"), "final")
analytic_path <- file.path(output_dir, "analysis_cross_sectional.parquet")
control_path  <- file.path(output_dir, "analysis_negative_control.parquet")
for (p in c(analytic_path, control_path))
  if (!file.exists(p)) stop("Input not found: ", p,
                            "\nRun scripts 01-03 for site '", site_name, "' first, or check --site_name / --output_root.")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- function(stub, ext = "csv") file.path(final_dir, paste0(stub, run_suffix, "_", site_name, ".", ext))
message(sprintf("Site: %s | origin %s | windows cut at day %s | horizon %d d",
                site_name, ORIGIN, paste(CUTS, collapse = ", "), HORIZON_DAYS))

MIN_CELL <- 10
COHORT_LEVELS <- c("Hypoxemic, ventilated (analytic)",
                   "Ventilated, non-hypoxemic (dosed, uninjured lung)",
                   "Not ventilated (no tidal volume)")
ANALYTIC <- COHORT_LEVELS[1]
EXPOSURES <- c(height_cm = "Height", pfvc = "PFVC", pbwpfvc = "PBW/PFVC")
AGE_FORMS <- c(linear = "age10", spline = "ns(age_at_admission, 4)")
WINDOW_BREAKS <- c(0, CUTS, HORIZON_DAYS)
WINDOW_LABELS <- sprintf("%g-%g d", head(WINDOW_BREAKS, -1), WINDOW_BREAKS[-1])
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
COHORT_COLOURS <- setNames(okabe[c(6, 1, 5)], COHORT_LEVELS)

# =============================================================================
# 1. Cohorts and survival from each cohort's time zero
# =============================================================================
analytic <- read_parquet(analytic_path)
controls <- read_parquet(control_path)
required_analytic <- c("hospitalization_id", "recorded_dttm", "death_dttm", "age_at_admission", "sex_category",
                       "race_category", "height_cm", "pfvc", "pbwpfvc",
                       if (ORIGIN == "admission") "admission_dttm")
required_controls <- c("hospitalization_id", "nc_cohort", "admission_dttm", "death_dttm", "age_at_admission",
                       "sex_category", "race_category", "height_cm", "pfvc", "pbwpfvc")
missing_analytic <- setdiff(required_analytic, names(analytic))
missing_controls <- setdiff(required_controls, names(controls))
if (length(missing_analytic)) stop("analysis_cross_sectional lacks: ", paste(missing_analytic, collapse = ", "))
if (length(missing_controls)) stop("analysis_negative_control lacks: ", paste(missing_controls, collapse = ", "))

cohort_frame <- bind_rows(
  analytic %>% transmute(hospitalization_id, cohort = ANALYTIC,
                         origin_dttm = if (ORIGIN == "index") recorded_dttm else admission_dttm,
                         death_dttm, age_at_admission, sex_category = as.character(sex_category),
                         race_category = as.character(race_category), height_cm, pfvc, pbwpfvc),
  controls %>% transmute(hospitalization_id, cohort = nc_cohort, origin_dttm = admission_dttm,
                         death_dttm, age_at_admission, sex_category = as.character(sex_category),
                         race_category = as.character(race_category), height_cm, pfvc, pbwpfvc)) %>%
  mutate(days_to_death = as.numeric(difftime(death_dttm, origin_dttm, units = "days")))

if (is_synthetic) {
  message("*** SYNTHETIC SITE: simulated death times for every cohort (plumbing only). ***")
  set.seed(20260918)
  death_probability <- if_else(cohort_frame$cohort == ANALYTIC, 0.35, 0.12)
  died <- rbinom(nrow(cohort_frame), 1L, death_probability) == 1L
  cohort_frame$days_to_death <- if_else(died, pmin(rlnorm(nrow(cohort_frame), log(9), 1.1), 365), NA_real_)
}

n_death_before_origin <- sum(cohort_frame$days_to_death < 0, na.rm = TRUE)
if (n_death_before_origin > 0)
  message("Dropping ", n_death_before_origin, " patients whose recorded death precedes their time zero (data error)")
cohort_frame <- cohort_frame %>%
  filter(is.na(days_to_death) | days_to_death >= 0,
         is.finite(age_at_admission), is.finite(height_cm), is.finite(pfvc), is.finite(pbwpfvc)) %>%
  mutate(event = as.integer(!is.na(days_to_death) & days_to_death <= HORIZON_DAYS),
         # a death at the time-zero timestamp gets one hour, so it enters the first window
         time = if_else(event == 1L, pmax(days_to_death, 1 / 24), HORIZON_DAYS),
         cohort = factor(cohort, COHORT_LEVELS),
         age10 = age_at_admission / 10,
         sex_category = factor(sex_category, c("Male", "Female")),
         race_category = factor(race_category, c("WHITE", "BLACK", "OTHER"))) %>%
  filter(!is.na(sex_category), !is.na(race_category), !is.na(cohort))
# exposures on one scale: per SD of the ANALYTIC cohort (script 04's convention)
analytic_sd <- cohort_frame %>% filter(cohort == ANALYTIC) %>%
  summarise(across(all_of(names(EXPOSURES)), ~ sd(.x)))

deaths_by_window <- cohort_frame %>% filter(event == 1L) %>%
  mutate(window = factor(cut(time, WINDOW_BREAKS, labels = WINDOW_LABELS, include.lowest = TRUE), WINDOW_LABELS)) %>%
  count(cohort, window, name = "deaths", .drop = FALSE) %>%
  mutate(window = paste("deaths", window)) %>%
  pivot_wider(names_from = window, values_from = deaths, values_fill = 0L)
counts <- cohort_frame %>% group_by(cohort) %>%
  summarise(n_patients = n(), deaths_60d = sum(event), .groups = "drop") %>%
  left_join(deaths_by_window, by = "cohort") %>%
  mutate(across(where(is.numeric), ~ if_else(.x < MIN_CELL, NA_integer_, as.integer(.x))), site = site_name)
write_csv(counts, out_file("nc_timesplit_counts"))
message("Deaths by cohort and window (NA = fewer than ", MIN_CELL, "):")
print(as.data.frame(counts %>% select(-site)), row.names = FALSE)

# =============================================================================
# 2. Models: single 60-day Cox + Schoenfeld test; time-split Cox by window
# =============================================================================
fit_cohort <- function(cohort_label, exposure, age_form) {
  d <- cohort_frame %>% filter(cohort == cohort_label) %>%
    mutate(z = .data[[exposure]] / analytic_sd[[exposure]])
  if (sum(d$event) < MIN_CELL * length(WINDOW_LABELS)) return(NULL)
  adjusters <- paste(AGE_FORMS[[age_form]], "+ sex_category + race_category")
  if (n_distinct(d$race_category) < 2) adjusters <- sub(" \\+ race_category", "", adjusters)

  single <- coxph(as.formula(paste("Surv(time, event) ~ z +", adjusters)), data = d)
  zph <- cox.zph(single, transform = "km")
  ph_row <- tibble(cohort = cohort_label, exposure = EXPOSURES[[exposure]], age_form = age_form,
                   hr_single = exp(coef(single)[["z"]]),
                   zph_chisq_exposure = zph$table["z", "chisq"], zph_p_exposure = zph$table["z", "p"],
                   zph_p_global = zph$table["GLOBAL", "p"], n = nrow(d), deaths = sum(d$event))

  split_data <- survSplit(Surv(time, event) ~ ., data = d, cut = CUTS, episode = "window_index") %>%
    mutate(window = factor(WINDOW_LABELS[window_index], WINDOW_LABELS))
  deaths_by_window <- split_data %>% group_by(window) %>% summarise(deaths = sum(event), .groups = "drop")
  by_window <- coxph(as.formula(paste("Surv(tstart, time, event) ~ z:window +", adjusters)), data = split_data)
  constant  <- coxph(as.formula(paste("Surv(tstart, time, event) ~ z +", adjusters)), data = split_data)
  lrt <- anova(constant, by_window)
  terms <- paste0("z:window", WINDOW_LABELS)
  log_hr <- coef(by_window)[terms]
  covariance <- vcov(by_window)[terms, terms]
  hr_rows <- tibble(cohort = cohort_label, exposure = EXPOSURES[[exposure]], age_form = age_form,
                    window = factor(WINDOW_LABELS, WINDOW_LABELS), log_hr = unname(log_hr),
                    se = sqrt(diag(covariance))) %>%
    left_join(deaths_by_window, by = "window") %>%
    mutate(hr = exp(log_hr), lo = exp(log_hr - 1.96 * se), hi = exp(log_hr + 1.96 * se),
           lrt_p_window_heterogeneity = lrt[["Pr(>|Chi|)"]][2])
  list(ph = ph_row, hr = hr_rows, covariance = covariance)
}

grid <- expand_grid(cohort_label = COHORT_LEVELS, exposure = names(EXPOSURES), age_form = names(AGE_FORMS))
fits <- pmap(grid, fit_cohort)
names(fits) <- pmap_chr(grid, paste, sep = "|")
fitted <- compact(fits)
skipped <- setdiff(names(fits), names(fitted))
if (length(skipped)) message("Skipped (fewer than ", MIN_CELL, " deaths per window on average): ",
                             paste(unique(sub("\\|.*", "", skipped)), collapse = "; "))

ph_table <- map_dfr(fitted, "ph") %>% mutate(site = site_name)
hr_table <- map_dfr(fitted, "hr") %>%
  # a window with fewer than MIN_CELL deaths is reported without its estimate
  mutate(across(c(log_hr, se, hr, lo, hi), ~ if_else(deaths < MIN_CELL, NA_real_, .x)),
         deaths = if_else(deaths < MIN_CELL, NA_integer_, as.integer(deaths)),
         origin = ORIGIN, site = site_name)
write_csv(ph_table, out_file("nc_ph_test"))
write_csv(hr_table, out_file("nc_timesplit_hr"))

# =============================================================================
# 3. Between cohorts: per-window contrasts and the early-minus-late difference-in-differences
# =============================================================================
contrast_rows <- function(control_label, exposure, age_form) {
  key_vent <- paste(ANALYTIC, exposure, age_form, sep = "|")
  key_ctrl <- paste(control_label, exposure, age_form, sep = "|")
  if (!all(c(key_vent, key_ctrl) %in% names(fitted))) return(NULL)
  vent <- fitted[[key_vent]]; ctrl <- fitted[[key_ctrl]]
  if (any(vent$hr$deaths < MIN_CELL) || any(ctrl$hr$deaths < MIN_CELL)) return(NULL)
  per_window <- tibble(comparison = paste(control_label, "vs ventilated (analytic)"),
                       exposure = EXPOSURES[[exposure]], age_form = age_form,
                       contrast = paste("window", WINDOW_LABELS),
                       log_ratio = vent$hr$log_hr - ctrl$hr$log_hr,
                       se = sqrt(vent$hr$se^2 + ctrl$hr$se^2))
  # early (first window) minus late (last window), ventilated minus control
  first_last <- c(1, rep(0, length(WINDOW_LABELS) - 2), -1)
  did <- sum(first_last * vent$hr$log_hr) - sum(first_last * ctrl$hr$log_hr)
  did_se <- sqrt(drop(t(first_last) %*% vent$covariance %*% first_last) +
                 drop(t(first_last) %*% ctrl$covariance %*% first_last))
  bind_rows(per_window,
            tibble(comparison = per_window$comparison[1], exposure = EXPOSURES[[exposure]], age_form = age_form,
                   contrast = sprintf("DiD: (%s minus %s), ventilated minus control",
                                      WINDOW_LABELS[1], WINDOW_LABELS[length(WINDOW_LABELS)]),
                   log_ratio = did, se = did_se))
}
contrast_table <- expand_grid(control_label = COHORT_LEVELS[-1], exposure = names(EXPOSURES),
                              age_form = names(AGE_FORMS)) %>%
  pmap_dfr(contrast_rows)
if (nrow(contrast_table)) {
  contrast_table <- contrast_table %>%
    mutate(ratio = exp(log_ratio), lo = exp(log_ratio - 1.96 * se), hi = exp(log_ratio + 1.96 * se),
           p = 2 * pnorm(-abs(log_ratio / se)), origin = ORIGIN, site = site_name)
}
write_csv(contrast_table, out_file("nc_timesplit_contrast"))

message("\nHeight, per analytic-cohort SD, by window (hazard ratios):")
print(as.data.frame(hr_table %>% filter(exposure == "Height") %>%
  transmute(cohort = substr(as.character(cohort), 1, 32), age_form, window, deaths,
            hr = round(hr, 3), lo = round(lo, 3), hi = round(hi, 3),
            p_heterogeneity = signif(lrt_p_window_heterogeneity, 2))), row.names = FALSE)
message("\nSchoenfeld tests of the single 60-day models (height):")
print(as.data.frame(ph_table %>% filter(exposure == "Height") %>%
  transmute(cohort = substr(cohort, 1, 32), age_form, hr_single = round(hr_single, 3),
            zph_p_exposure = signif(zph_p_exposure, 2), zph_p_global = signif(zph_p_global, 2))), row.names = FALSE)
if (nrow(contrast_table)) {
  message("\nVentilated vs never ventilated (ratio of HRs; < 1 = height more protective when ventilated):")
  print(as.data.frame(contrast_table %>% filter(exposure == "Height", grepl("^Not ventilated", comparison)) %>%
    transmute(age_form, contrast, ratio = round(ratio, 3), lo = round(lo, 3), hi = round(hi, 3),
              p = signif(p, 2))), row.names = FALSE)
}

# =============================================================================
# 4. Figure: HR by window, per cohort, for each exposure and age form
# =============================================================================
plot_data <- hr_table %>% filter(!is.na(hr))
if (nrow(plot_data)) {
  figure <- ggplot(plot_data, aes(window, hr, colour = cohort)) +
    geom_hline(yintercept = 1, colour = "grey70") +
    geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.5), size = 0.3) +
    facet_grid(exposure ~ age_form, scales = "free_y",
               labeller = labeller(age_form = c(linear = "linear age", spline = "ns(age, 4)"))) +
    scale_y_log10() +
    scale_colour_manual(values = COHORT_COLOURS, name = NULL) +
    labs(title = sprintf("Mortality HR per analytic-cohort SD, by time window (%s%s)", site_name,
                         if (is_synthetic) ", SYNTHETIC: simulated deaths" else ""),
         subtitle = paste0("Time zero: ", if (ORIGIN == "index") "analytic index / control admission" else "admission (all cohorts)",
                           ". A ventilation-specific effect shows as early divergence and late agreement."),
         x = "days from time zero", y = "hazard ratio (log scale)") +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom", legend.direction = "vertical")
  ggsave(out_file("nc_timesplit", "pdf"), figure, width = 10, height = 9)
}
message("17_nc_time_split complete -> ", final_dir)
