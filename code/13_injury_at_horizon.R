# =============================================================================
# Script 13 (injury at horizon): is a lower PFVC, at a given VT/PBW, associated
# with a worse injury marker at 48 hours?
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# A fixed-horizon, between-patient contrast, the paper's primary exposure set
# (VT/PBW as the clinician's dose, log PFVC as the size term, with PBW/PFVC as
# the companion) applied to an organ-injury marker instead of mortality:
#
#   log marker at H  ~  log marker at baseline + log PFVC (per SD)
#                        + mean VT/PBW over [0, H) + non-respiratory SOFA
#                        + log baseline SF [+ BMI, driving pressure only]
#                        [+ ns(age, 4) + sex + race]
# BMI carries height (weight over height squared) and enters only for the
# pressure-derived marker.
#
# fitted adjusted (with the demographics, which are PFVC's parents) and
# unadjusted, exactly as every exposure model in script 04 is. The companion
# model replaces log PFVC by log PBW/PFVC. Log PBW and log PFVC are never
# entered together: after age, sex and race, PBW is height, and the pair is not
# identifiable.
#
# What the primary conditions on, stated because it matters: the patient is
# alive at H, not on renal replacement before H, and has a baseline draw in the
# first 12 hours and a draw in the outcome window (the last value in (H - 24, H]
# for H >= 48, (H - 12, H] for H = 24). Death before H is not ignorable if PFVC
# protects against it, so a composite-rank sensitivity ranks death before H
# worst, RRT before H next, then the marker value, and fits the same models on
# the rank (0 to 1).
#
# Horizons: 48 h PRIMARY; 24 and 72 h sensitivities. One marker at a time
# (PBWPFVC_INJ_MARKER = creatinine | ne_equiv | platelets | bilirubin | sf | dp;
# creatinine by default). NE-equivalents are a two-part outcome (any pressor at
# H, and the log dose given a pressor) because most patients are at zero. The
# dose per kg carries -log(weight) = -log(BMI) - 2 log(height), and height is
# what identifies log PFVC once age, sex and race are in, so the per-kg dose has
# a mechanical negative association with PFVC; the absolute dose (mcg/min)
# carries the opposite sign. Both are reported and neither is height-neutral;
# the binary part is the read for this marker.
#
# Inputs: the shared daily panel's sources (10_panel_common.R: base, wf, labs,
# NE-equivalent administrations, CRRT). No dependence on the joint-model tables.
# Outputs: final/injury_at_horizon_{marker}_{site}.csv  one row per horizon x
#          exposure x adjustment x outcome type, with counts of who was excluded
#          final/injury_at_horizon_{marker}_{site}.pdf  the forest
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(arrow)
  library(here)
  library(splines)
})
rm(list = ls())
source("utils/config.R")

site_name  <- config$site_name
output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

HORIZON      <- 28L
MAX_VENT_DAY <- 27L
is_synthetic <- identical(site_name, "synthetic_clif")
PANEL_NORM   <- "pfvc"
source(here("code", "10_panel_common.R"))

MARKER   <- Sys.getenv("PBWPFVC_INJ_MARKER", "creatinine")
stopifnot(MARKER %in% c("creatinine", "ne_equiv", "platelets", "bilirubin", "sf", "dp"))
HORIZONS <- as.numeric(strsplit(Sys.getenv("PBWPFVC_INJ_HORIZONS_H", "48,24,72"), ",")[[1]])
BASE_WINDOW_H <- 12
okabe <- c("#0072B2", "#E69F00", "#009E73", "#D55E00")
message("=== 13_injury_at_horizon: marker ", MARKER, ", horizons ", paste(HORIZONS, collapse = "/"), " h, site ", site_name, " ===")

# =============================================================================
# Per-measurement marker series in hours since the index
# =============================================================================
b0 <- base %>% select(hospitalization_id, t0, pbw)
hrs <- function(dttm, t0) as.numeric(difftime(dttm, t0, units = "hours"))
series <- switch(MARKER,
  creatinine = , platelets = , bilirubin = {
    cat_name <- c(creatinine = "creatinine", platelets = "platelet_count", bilirubin = "bilirubin_total")[[MARKER]]
    read_parquet(file.path(output_dir, "cohort_labs_clean.parquet")) %>%
      filter(lab_category == cat_name, !is.na(lab_value_numeric), lab_value_numeric > 0) %>%
      inner_join(b0, by = "hospitalization_id") %>%
      transmute(hospitalization_id, h = hrs(lab_result_dttm, t0), value = lab_value_numeric)
  },
  ne_equiv = {
    # the dose in force at each administration record; hours with no record are zero
    read_parquet(file.path(output_dir, "ne_equiv_admin.parquet")) %>%
      inner_join(b0, by = "hospitalization_id") %>%
      transmute(hospitalization_id, h = hrs(admin_dttm, t0), value = ne_equiv_total)
  },
  sf = {
    t0n <- b0 %>% transmute(hospitalization_id, t0n = as.numeric(t0))
    fio2_dt[spo2_dt, roll = 4 * 3600, on = .(hospitalization_id, t)] %>% as_tibble() %>%
      filter(!is.na(fio2_set)) %>%
      mutate(fio2_frac = if_else(fio2_set > 1.5, fio2_set / 100, fio2_set), value = spo2_clamped / fio2_frac) %>%
      filter(is.finite(value)) %>% inner_join(t0n, by = "hospitalization_id") %>%
      transmute(hospitalization_id, h = (t - t0n) / 3600, value)
  },
  dp = wf %>% filter(!is.na(plateau_pressure_obs), !is.na(peep_set), plateau_pressure_obs - peep_set > 0) %>%
    transmute(hospitalization_id, h = hrs(recorded_dttm, t0), value = plateau_pressure_obs - peep_set)
)
worse_is <- if (MARKER %in% c("sf")) "lower" else "higher"   # direction of injury on the marker scale

# baseline: the first value in the first BASE_WINDOW_H hours (NE: the peak dose in that window, zero if none)
baseline <- if (MARKER == "ne_equiv") {
  series %>% filter(h >= 0, h < BASE_WINDOW_H) %>% group_by(hospitalization_id) %>%
    summarise(y0 = max(value), .groups = "drop")
} else {
  series %>% filter(h >= 0, h < BASE_WINDOW_H) %>% group_by(hospitalization_id) %>%
    slice_min(h, n = 1, with_ties = FALSE) %>% ungroup() %>% select(hospitalization_id, y0 = value)
}
if (MARKER == "ne_equiv") baseline <- b0 %>% select(hospitalization_id) %>%
  left_join(baseline, by = "hospitalization_id") %>% mutate(y0 = coalesce(y0, 0))

# index-day worst SF (the severity covariate of the hazard models), from the shared daily panel
sf0 <- sf_daily %>% filter(vent_day == 0L) %>% select(hospitalization_id, sf_0 = sf)

# RRT start (hours); NA where the site has no CRRT table or the patient none
crrt_available <- readRDS(file.path(output_dir, "crrt_available.rds"))
rrt <- read_parquet(file.path(output_dir, "cohort_crrt.parquet")) %>%
  inner_join(b0, by = "hospitalization_id") %>%
  group_by(hospitalization_id) %>% summarise(rrt_h = min(hrs(recorded_dttm, t0)), .groups = "drop")

# =============================================================================
# One horizon
# =============================================================================
fit_horizon <- function(H) {
  win <- if (H >= 48) 24 else 12
  # outcome: the last value in (H - win, H]; NE: the peak dose in that window, zero if none
  yH <- if (MARKER == "ne_equiv") {
    # every patient has a dose in the window: zero when no pressor was given
    # (without this the "any pressor" model saw only patients on a pressor and
    # separated completely)
    b0 %>% select(hospitalization_id) %>%
      left_join(series %>% filter(h > H - win, h <= H) %>% group_by(hospitalization_id) %>%
                  summarise(yH = max(value), .groups = "drop"), by = "hospitalization_id") %>%
      mutate(yH = coalesce(yH, 0))
  } else {
    series %>% filter(h > H - win, h <= H) %>% group_by(hospitalization_id) %>%
      slice_max(h, n = 1, with_ties = FALSE) %>% ungroup() %>% select(hospitalization_id, yH = value)
  }
  # dose over [0, H): mean VT/PBW of the set tidal volumes in the window
  dose <- wf %>% mutate(h = hrs(recorded_dttm, t0)) %>% filter(h >= 0, h < H) %>%
    inner_join(b0 %>% select(hospitalization_id, pbw), by = "hospitalization_id") %>%
    group_by(hospitalization_id) %>% summarise(vtpbw_H = mean(tidal_volume_set / pbw), .groups = "drop")
  d <- base %>%
    left_join(baseline, by = "hospitalization_id") %>%
    left_join(yH, by = "hospitalization_id") %>%
    left_join(dose, by = "hospitalization_id") %>%
    left_join(rrt, by = "hospitalization_id") %>%
    left_join(sf0, by = "hospitalization_id") %>%
    mutate(
      dead_before_H = !is.na(death_time_days) & death_time_days * 24 <= H,
      rrt_before_H  = !is.na(rrt_h) & rrt_h <= H,
      log_pfvc = log(pfvc_gli), ldisc = log(pbw / pfvc_gli),
      log_pfvc_sd = as.numeric(scale(log_pfvc)), ldisc_sd = as.numeric(scale(ldisc)),
      log_sf_0 = log(sf_0)
    )
  # the survival-side SF baseline and np_sofa come from base via the shared panel;
  # sf_0 is the index-day worst SF (day-0 value)
  if (MARKER == "ne_equiv") d <- b0 %>% select(hospitalization_id) %>% left_join(d, by = "hospitalization_id")
  n_all <- nrow(d)
  cc <- d %>% filter(!dead_before_H, !rrt_before_H | MARKER != "creatinine",
                     !is.na(yH), !is.na(y0), !is.na(vtpbw_H), !is.na(np_sofa), !is.na(sf_0),
                     if (MARKER == "dp") !is.na(bmi) else TRUE,
                     !is.na(age10), !is.na(sex_category), !is.na(race_category))
  counts <- tibble(horizon_h = H, n_cohort = n_all, n_dead_before_H = sum(d$dead_before_H),
                   n_rrt_before_H = sum(d$rrt_before_H),
                   n_no_baseline = sum(!d$dead_before_H & is.na(d$y0)),
                   n_no_outcome_value = sum(!d$dead_before_H & !is.na(d$y0) & is.na(d$yH)),
                   n_analysed = nrow(cc))
  message(sprintf("  H = %g h: %d in cohort, %d dead before H, %d on RRT before H, %d without baseline, %d without a value in the window, %d analysed",
                  H, n_all, counts$n_dead_before_H, counts$n_rrt_before_H, counts$n_no_baseline,
                  counts$n_no_outcome_value, nrow(cc)))
  if (nrow(cc) < 50) return(list(counts = counts, rows = NULL))

  base_rhs <- paste("vtpbw_H + np_sofa + log_sf_0", if (MARKER == "dp") "+ bmi" else "")
  demo_rhs <- "ns(age10, 4) + sex_category + race_category"
  # expo: the coefficient reported; extra: further right-hand-side terms (the baseline)
  one <- function(dat, lhs, expo, adjusted, outcome_type, family = "gaussian", extra = "log_y0") {
    f <- as.formula(paste(lhs, "~", expo, "+", extra, "+", base_rhs, if (adjusted) paste("+", demo_rhs) else ""))
    fit <- if (family == "gaussian") lm(f, data = dat) else glm(f, data = dat, family = binomial)
    co <- summary(fit)$coefficients[expo, ]
    ci <- if (family == "gaussian") confint(fit)[expo, ] else suppressMessages(confint.default(fit)[expo, ])
    tibble(horizon_h = H, marker = MARKER, outcome_type = outcome_type, exposure = expo,
           adjustment = if (adjusted) "adjusted" else "unadjusted",
           estimate = unname(co[1]), se = unname(co[2]), lo = unname(ci[1]), hi = unname(ci[2]),
           p = unname(co[4]), n = nrow(dat), worse_is = worse_is,
           note = if (family == "gaussian") "log marker at H, baseline as covariate" else "logistic: any pressor at H")
  }
  rows <- list()
  if (MARKER == "ne_equiv") {
    # actual weight from BMI and height (kg), for the absolute dose in mcg/min
    cc <- cc %>% mutate(any_H = as.integer(yH > 0), log_y0 = log(y0 + 0.01),
                        weight_kg = bmi * (height_cm / 100)^2)
    on <- cc %>% filter(yH > 0) %>% mutate(log_yH = log(yH), log_yH_abs = log(yH * weight_kg),
                                           log_y0_abs = log(y0 * weight_kg + 0.01))
    for (expo in c("log_pfvc_sd", "ldisc_sd")) for (adj in c(TRUE, FALSE)) {
      rows[[length(rows) + 1]] <- one(cc, "any_H", expo, adj, "any pressor at H", "binomial")
      if (nrow(on) >= 50) {
        rows[[length(rows) + 1]] <- one(on, "log_yH", expo, adj, "log dose per kg given any") %>%
          mutate(note = "mcg/kg/min; carries -2 log(height): mechanically NEGATIVE in PFVC")
        rows[[length(rows) + 1]] <- one(on, "log_yH_abs", expo, adj, "log absolute dose given any", extra = "log_y0_abs") %>%
          mutate(note = "mcg/min; heavier patients need more drug: mechanically POSITIVE in PFVC")
      }
    }
  } else {
    cc <- cc %>% mutate(log_yH = log(yH), log_y0 = log(y0))
    for (expo in c("log_pfvc_sd", "ldisc_sd")) for (adj in c(TRUE, FALSE))
      rows[[length(rows) + 1]] <- one(cc, "log_yH", expo, adj, "log marker at H")
  }
  # composite-rank sensitivity: death before H worst, RRT before H next, then the marker
  # (worse direction first), on everyone with a baseline; rank scaled to (0, 1)
  comp <- d %>% filter(!is.na(y0), !is.na(vtpbw_H), !is.na(np_sofa), !is.na(sf_0),
                       if (MARKER == "dp") !is.na(bmi) else TRUE,
                       !is.na(age10), !is.na(sex_category), !is.na(race_category)) %>%
    mutate(score = case_when(dead_before_H ~ Inf,
                             rrt_before_H & MARKER == "creatinine" ~ 1e9,
                             is.na(yH) ~ NA_real_,
                             worse_is == "higher" ~ yH, TRUE ~ -yH)) %>%
    filter(!is.na(score)) %>%
    mutate(rank01 = (rank(score) - 0.5) / n(), log_y0 = if (MARKER == "ne_equiv") log(y0 + 0.01) else log(y0))
  if (nrow(comp) >= 50)
    for (expo in c("log_pfvc_sd", "ldisc_sd")) for (adj in c(TRUE, FALSE))
      rows[[length(rows) + 1]] <- one(comp, "rank01", expo, adj, "composite rank (death, RRT, marker)") %>%
        mutate(note = sprintf("rank 0-1; %d deaths and %d RRT before H ranked worst",
                                              sum(comp$dead_before_H), sum(comp$rrt_before_H & MARKER == "creatinine")))
  list(counts = counts, rows = bind_rows(rows))
}

res <- map(HORIZONS, fit_horizon)
counts  <- map_dfr(res, "counts") %>% mutate(marker = MARKER, site = site_name)
results <- map_dfr(res, "rows") %>% mutate(site = site_name)
if (nrow(results) == 0) stop("no horizon had enough patients to fit")

write_csv(results, file.path(final_dir, paste0("injury_at_horizon_", MARKER, "_", site_name, ".csv")))
write_csv(counts,  file.path(final_dir, paste0("injury_at_horizon_counts_", MARKER, "_", site_name, ".csv")))
message("\n--- counts"); print(as.data.frame(counts), row.names = FALSE)
message("\n--- log PFVC per SD (direction of a LOWER PFVC is the negative of this)")
print(as.data.frame(results %>% filter(exposure == "log_pfvc_sd") %>%
                      select(horizon_h, outcome_type, adjustment, estimate, lo, hi, p, n) %>%
                      mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)

# forest: log PFVC per SD by horizon, outcome type and adjustment
fp <- results %>% filter(exposure == "log_pfvc_sd") %>%
  mutate(adjustment = factor(adjustment, c("adjusted", "unadjusted")),
         horizon = factor(paste0(horizon_h, " h"), paste0(sort(unique(horizon_h)), " h")))
p <- ggplot(fp, aes(estimate, horizon, colour = adjustment)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5)) +
  facet_wrap(~ outcome_type, scales = "free_x", ncol = 1) +
  scale_colour_manual(values = okabe[1:2]) +
  labs(title = sprintf("%s at the horizon per SD of log PFVC, at a given VT/PBW (%s)", MARKER, site_name),
       subtitle = sprintf("worse injury is %s on this marker; a protective PFVC is %s",
                          worse_is, if (worse_is == "higher") "negative" else "positive"),
       x = "Change per SD of log PFVC (log units; rank units for the composite; log-odds for any pressor)", y = NULL) +
  theme_minimal(base_size = 11)
ggsave(file.path(final_dir, paste0("injury_at_horizon_", MARKER, "_", site_name, ".pdf")), p, width = 9, height = 7)
message("13_injury_at_horizon complete -> ", final_dir)
