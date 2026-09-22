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
# The post-intubation rate beside it is figure 4's own estimate, the joint model's
# log_pfvc_sd:vent_day (jm_estimates_pfvc_*, creatinine from its dialysis-as-third-cause
# twin), when those tables exist; at MIMIC it matched the longitudinal model alone to
# the fourth decimal for platelets, so the two are comparable. It is estimated on the whole
# cohort, while the pre-period sample is only the patients intubated after a day or more
# in hospital, so the two numbers come from different samples and different models.
#
# The STACKED within-patient version removes both differences: the pre rows (days -7 to -1
# before the first IMV record) and the post rows (days 1 to 7 of the figure-4 panel) of
# the same patients, two or more of each, in one mixed model,
#   log y ~ post + day_pre + day_post
#           + log_pfvc_sd + log_pfvc_sd:post + log_pfvc_sd:day_pre + log_pfvc_sd:day_post
#           [+ ns(age, 4) + sex + race]
#   random intercept, post-intubation jump and both slopes per patient (pdDiag), ML
# where log_pfvc_sd:day_pre and log_pfvc_sd:day_post are the divergence before and after
# intubation and their difference is the change at intubation, within patient: the test.
# The model is symmetric, so it carries no ventilator-period covariates (lags, dose,
# baseline, non-respiratory SOFA) and does not model death or extubation, unlike figure
# 4's joint model; its post rate is read beside the joint model's, not in place of it.
# The two clocks differ by the first-IMV-to-index gap (a median 4 h at MIMIC).
# Creatinine excludes patients on renal replacement before the index (ESRD included):
# their pre-intubation creatinine is set by dialysis, and the post panel already has none.
#
# Output: final/injury/jm_pre_placebo_{panel}_{site}.csv  the pre-period placebo
#         final/injury/jm_pre_stacked_{panel}_{site}.csv  the stacked within-patient model
#         (counts of 1-9 blanked)
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

# the post-intubation rate: figure 4's own estimate, the joint model's log_pfvc_sd:vent_day
# (creatinine from its dialysis-as-third-cause twin, as figure 4 reads it)
read_post <- function(tag) {
  path <- file.path(final_dir, paste0("jm_estimates_", tag, "pfvc_", h_suffix, "_", site_name, ".csv"))
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>%
    filter(block == "longitudinal", model == "main",
           vapply(strsplit(term, ":"), setequal, logical(1), c("log_pfvc_sd", "vent_day"))) %>%
    transmute(marker, adjustment, post_estimate = estimate, post_lo = lo, post_hi = hi, post_rhat = rhat,
              post_n_patients = n_patients)
}
post_main <- read_post(""); post_rrt <- read_post("rrtcause_")
post_rates <- bind_rows(if (!is.null(post_main)) post_main %>% filter(!(marker == "creatinine" & !is.null(post_rrt))),
                        if (!is.null(post_rrt)) post_rrt %>% filter(marker == "creatinine"))
if (!nrow(post_rates)) post_rates <- NULL

# creatinine: no patient on renal replacement before the index (ESRD included)
rrt_before <- surv$hospitalization_id[surv$rrt_before_index]
eligible_pre <- function(m) if (m == "creatinine") pre %>% filter(!hospitalization_id %in% rrt_before) else pre

rows <- map_dfr(MARKERS, function(m) {
  d <- eligible_pre(m) %>% filter(!is.na(.data[[PRE_COLUMNS[[m]]]])) %>%
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
        if (is.null(post_rates)) " (no joint-model estimates yet: run the figure 4 fits for the post rate)" else "")
print(as.data.frame(rows %>% select(any_of(c("marker", "adjustment", "status", "pre_estimate", "pre_lo", "pre_hi", "pre_p",
                                             "post_estimate", "post_lo", "post_hi", "n_patients"))) %>%
                      mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
out_path <- file.path(final_dir, paste0("jm_pre_placebo_", h_suffix, "_", site_name, ".csv"))
write_csv(mask_small_counts(rows), out_path)

# =============================================================================
# The stacked within-patient model (see the header)
# =============================================================================
long <- read_parquet(file.path(output_dir, paste0("jm_long_", h_suffix, ".parquet")),
                     col_select = c("hospitalization_id", "vent_day", unname(PRE_COLUMNS)))
coef_named <- function(b, comps) {
  hit <- names(b)[vapply(strsplit(names(b), ":"), setequal, logical(1), comps)]
  stopifnot(length(hit) == 1L); hit
}
stacked <- map_dfr(MARKERS, function(m) {
  col <- PRE_COLUMNS[[m]]
  pre_m  <- eligible_pre(m) %>% filter(!is.na(.data[[col]])) %>%
    transmute(hospitalization_id, day = vent_day, y = .data[[col]], post = 0)
  post_m <- long %>% filter(vent_day >= 1, vent_day <= JM_HORIZON, !is.na(.data[[col]]), .data[[col]] > 0) %>%
    transmute(hospitalization_id, day = vent_day, y = .data[[col]], post = 1)
  d <- bind_rows(pre_m, post_m) %>%
    group_by(hospitalization_id) %>% filter(sum(post == 0) >= 2L, sum(post == 1) >= 2L) %>% ungroup() %>%
    inner_join(surv %>% select(hospitalization_id, log_pfvc_sd, age10, sex_category, race_category), by = "hospitalization_id") %>%
    filter(!is.na(log_pfvc_sd)) %>%
    mutate(log_y = log(y), day_pre = day * (1 - post), day_post = day * post, id = factor(hospitalization_id))
  n_pts <- n_distinct(d$id)
  message(m, " (stacked): ", n_pts, " patients with 2+ days both before and after intubation, ", nrow(d), " patient-days")
  if (n_pts < MIN_PATIENTS)
    return(tibble(marker = m, adjustment = c("adjusted", "unadjusted"), status = "skipped",
                  reason = paste0("fewer than ", MIN_PATIENTS, " patients with 2+ days before and after intubation"),
                  n_patients = n_pts, n_obs = nrow(d)))
  imap_dfr(c(adjusted = TRUE, unadjusted = FALSE), function(adj, adj_lab) {
    rhs <- paste(c("post", "day_pre", "day_post", "log_pfvc_sd", "log_pfvc_sd:post",
                   "log_pfvc_sd:day_pre", "log_pfvc_sd:day_post", if (adj) DEMO_RHS), collapse = " + ")
    f <- lme(as.formula(paste("log_y ~", rhs)), random = list(id = pdDiag(~ post + day_pre + day_post)), data = d,
             method = "ML", control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200))
    b <- fixef(f); V <- as.matrix(vcov(f))
    pre_t <- coef_named(b, c("log_pfvc_sd", "day_pre")); post_t <- coef_named(b, c("log_pfvc_sd", "day_post"))
    w <- setNames(numeric(length(b)), names(b))
    contrast <- function(weights) { w[names(weights)] <- weights; c(est = sum(w * b), se = sqrt(as.numeric(t(w) %*% V %*% w))) }
    parts <- list(before = contrast(setNames(1, pre_t)), after = contrast(setNames(1, post_t)),
                  change = contrast(setNames(c(-1, 1), c(pre_t, post_t))))
    imap_dfr(parts, ~ tibble(quantity = .y, estimate = .x[["est"]], se = .x[["se"]])) %>%
      mutate(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se, p = 2 * pnorm(-abs(estimate / se)),
             marker = m, adjustment = adj_lab, status = "fitted", reason = NA_character_,
             n_patients = n_pts, n_obs = nrow(d), .before = 1)
  })
}) %>%
  mutate(unit = "log marker per day per SD of log PFVC; change = after - before, within patient",
         panel = h_suffix, site = site_name)
message("\nStacked within-patient model: the divergence before and after intubation in the same patients")
print(as.data.frame(stacked %>% select(any_of(c("marker", "adjustment", "quantity", "estimate", "lo", "hi", "p", "n_patients"))) %>%
                      mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
stacked_path <- file.path(final_dir, paste0("jm_pre_stacked_", h_suffix, "_", site_name, ".csv"))
write_csv(mask_small_counts(stacked), stacked_path)
message("28_pre_period_placebo complete -> ", out_path, ", ", stacked_path)
