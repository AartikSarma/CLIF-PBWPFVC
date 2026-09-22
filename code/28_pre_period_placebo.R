# =============================================================================
# Script 28 (placebo pre-period): does the lung-size divergence exist before intubation?
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# The difference-in-differences behind figure 4 (27_control_comparison.R) reads the
# ventilated divergence by predicted lung size against the no-support control. Its
# pre-trend check: if the divergence is the ventilator's, patients with smaller
# predicted lungs should NOT already be diverging in the days before intubation.
#
#   log marker ~ day + log_pfvc_sd + log_pfvc_sd:day   [+ ns(age, 4) + sex + race]
#   random intercept and slope per patient (unstructured), maximum likelihood
#
# on the pre-intubation panel (21_biotrauma_panel.R, jm_pre_7d: days -PRE_DAYS to -1
# before the first IMV record), patients with two or more pre-intubation days.
# log_pfvc_sd:day is the pre-period divergence, on the scale of the post-intubation
# rate (log marker per day per SD of log PFVC; the SD is the whole cohort's, as in the
# fits). The model is deliberately minimal: no day-0 baseline and no index-day SOFA,
# because both are the end of the pre-trajectory, so adjusting for them adjusts for
# the outcome; no ventilator terms, because none exist before intubation. Nothing
# dies or is extubated before intubation, so no joint model is needed.
#
# The post-intubation rate beside it is the longitudinal model alone from the shape
# check (22_biotrauma_fit.R, PBWPFVC_JM_SHAPE_ONLY=1, jm_shape_*), the same kind of
# model, when that table exists. It is estimated on the whole cohort, while the
# pre-period sample is only the patients intubated after a day or more in hospital:
# the within-patient version (the pre and post rows of the same patients stacked,
# with a divergence x post-intubation term) is the follow-on.
#
# Output: final/injury/jm_pre_placebo_{panel}_{site}.csv (counts of 1-9 blanked)
# Usage:  PBWPFVC_JM_GRID=daily PBWPFVC_JM_HORIZON=7 PBWPFVC_JM_MARKERS=platelets Rscript code/28_pre_period_placebo.R
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(here); library(splines); library(nlme) })
rm(list = ls())
source("utils/config.R")
site_name  <- config$site_name
output_dir <- config$output_dir
final_dir  <- final_dir_for("injury")
if (config$cohort != "imv") stop("the pre-intubation placebo is for the ventilated cohort: unset PBWPFVC_COHORT")
source(here("code", "20_biotrauma_grid.R"))   # JM_GRID, h_suffix
if (JM_GRID != "daily") stop("the pre-intubation panel is built on the daily grid: set PBWPFVC_JM_GRID=daily")

MARKERS <- trimws(strsplit(Sys.getenv("PBWPFVC_JM_MARKERS", "platelets"), ",")[[1]])
PRE_COLUMNS <- c(creatinine = "creatinine", platelets = "platelets", bilirubin = "bilirubin")
if (!all(MARKERS %in% names(PRE_COLUMNS)))
  stop("the pre-intubation panel holds creatinine, platelets and bilirubin; got ", paste(MARKERS, collapse = ", "))
MIN_PATIENTS <- 50L
DEMO_RHS <- "ns(age10, 4) + sex_category + race_category"

pre_path <- file.path(output_dir, paste0("jm_pre_", h_suffix, ".parquet"))
if (!file.exists(pre_path)) stop("no pre-intubation panel: rebuild it with the current 21_biotrauma_panel.R (daily grid)")
pre  <- read_parquet(pre_path)
surv <- read_parquet(file.path(output_dir, paste0("jm_surv_", h_suffix, ".parquet")))
PRE_DAYS <- -min(pre$vent_day, 0)

# the post-intubation rate from the shape check, where it has been run (the same model class)
# (creatinine from its dialysis-as-third-cause twin, as figure 4 reads it)
read_shape <- function(tag) {
  path <- file.path(final_dir, paste0("jm_shape_", tag, "pfvc_", h_suffix, "_", site_name, ".csv"))
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>% filter(shape == "linear", quantity == "rate", segment == "all") %>%
    transmute(marker, adjustment, post_estimate = estimate, post_lo = lo, post_hi = hi, post_n_patients = n_patients)
}
shape_main <- read_shape(""); shape_rrt <- read_shape("rrtcause_")
post_rates <- bind_rows(if (!is.null(shape_main)) shape_main %>% filter(!(marker == "creatinine" & !is.null(shape_rrt))),
                        if (!is.null(shape_rrt)) shape_rrt %>% filter(marker == "creatinine"))
if (!nrow(post_rates)) post_rates <- NULL

rows <- map_dfr(MARKERS, function(m) {
  d <- pre %>% filter(!is.na(.data[[PRE_COLUMNS[[m]]]])) %>%
    inner_join(surv %>% select(hospitalization_id, log_pfvc_sd, age10, sex_category, race_category), by = "hospitalization_id") %>%
    filter(!is.na(log_pfvc_sd)) %>%
    group_by(hospitalization_id) %>% filter(n() >= 2L) %>% ungroup() %>%
    mutate(log_y = log(.data[[PRE_COLUMNS[[m]]]]), id = factor(hospitalization_id))
  n_pts <- n_distinct(d$id)
  message(m, ": ", nrow(d), " pre-intubation patient-days, ", n_pts, " patients with 2+ days (",
          round(100 * n_pts / nrow(surv), 1), "% of the ventilated cohort)")
  if (n_pts < MIN_PATIENTS) {
    message("  too few patients with a pre-intubation trajectory for ", m, " (", n_pts, " < ", MIN_PATIENTS, "): no fit")
    return(tibble(marker = m, adjustment = c("adjusted", "unadjusted"), status = "skipped",
                  reason = paste0("fewer than ", MIN_PATIENTS, " patients with 2+ pre-intubation days"),
                  n_patients = n_pts, n_obs = nrow(d)))
  }
  imap_dfr(c(adjusted = TRUE, unadjusted = FALSE), function(adj, adj_lab) {
    rhs <- paste(c("vent_day", "log_pfvc_sd", "log_pfvc_sd:vent_day", if (adj) DEMO_RHS), collapse = " + ")
    f <- lme(as.formula(paste("log_y ~", rhs)), random = ~ vent_day | id, data = d, method = "ML",
             control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200))
    b <- fixef(f); se <- sqrt(diag(as.matrix(vcov(f))))
    term <- names(b)[vapply(strsplit(names(b), ":"), setequal, logical(1), c("log_pfvc_sd", "vent_day"))]   # R orders the pair by appearance
    stopifnot(length(term) == 1L)
    tibble(marker = m, adjustment = adj_lab, status = "fitted", reason = NA_character_,
           pre_estimate = b[[term]], pre_se = se[[term]],
           pre_lo = b[[term]] - 1.96 * se[[term]], pre_hi = b[[term]] + 1.96 * se[[term]],
           pre_p = 2 * pnorm(-abs(b[[term]] / se[[term]])),
           pre_level_day_minus1 = b[["log_pfvc_sd"]] - b[[term]],   # the size gap on the last day before intubation
           n_patients = n_pts, n_obs = nrow(d))
  })
}) %>%
  mutate(pre_days = PRE_DAYS, cohort_patients = nrow(surv),
         unit = "log marker per day per SD of log PFVC", panel = h_suffix, site = site_name)
if (!is.null(post_rates)) rows <- rows %>% left_join(post_rates, by = c("marker", "adjustment"))

message("\nPre-intubation divergence (placebo) beside the post-intubation rate",
        if (is.null(post_rates)) " (no shape-check table: run 22 with PBWPFVC_JM_SHAPE_ONLY=1 for the post rate)" else "")
print(as.data.frame(rows %>% select(any_of(c("marker", "adjustment", "status", "pre_estimate", "pre_lo", "pre_hi", "pre_p",
                                             "post_estimate", "post_lo", "post_hi", "n_patients"))) %>%
                      mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
out_path <- file.path(final_dir, paste0("jm_pre_placebo_", h_suffix, "_", site_name, ".csv"))
write_csv(mask_small_counts(rows), out_path)
message("28_pre_period_placebo complete -> ", out_path)
