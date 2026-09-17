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
#          final/injury_channels_{marker}_{site}.csv/.pdf  the channel decomposition
#          (13_biotrauma_grid.R): the size effect identified through each GLI
#          input separately, with a test that the four agree
#          Supports for the channel read (2026-09-15), on the same fits:
#          final/injury_nested_{marker}_{site}.csv        nested ladder: demographics
#              only, size only, both, the pieces free, pieces + free demographics;
#              likelihood-ratio tests and AIC (does PFVC add to age/sex/race?)
#          final/injury_dose_channels_{marker}_{site}.csv  dose x piece: does each
#              piece's effect scale with the delivered VT/PBW (a lung-size
#              mechanism does, a direct age effect does not)
#          final/injury_sf_channels_{marker}_{site}.csv    severity x piece: the same
#              with baseline SF (the baby-lung gradient; skipped for sf)
#          final/injury_negctrl_{marker}_{site}.csv        negative control: the pieces
#              on the BASELINE marker, before ventilation can act, on everyone with
#              a baseline; each piece's direct effect, read against the at-H result
#          final/injury_supports_{marker}_{site}.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(arrow)
  library(here)
  library(splines)
  library(patchwork)
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
source(here("code", "13_biotrauma_grid.R"))   # pfvc_channels(), channels_equal_p()

MARKER   <- Sys.getenv("PBWPFVC_INJ_MARKER", "creatinine")
stopifnot(MARKER %in% c("creatinine", "ne_equiv", "platelets", "bilirubin", "sf", "dp"))
# the never-intubated control has no ventilator dose: the dose term and the
# dose x piece interaction are dropped, intubation before H excludes like death
HAS_DOSE <- config$cohort == "imv"
if (!HAS_DOSE && MARKER == "dp") stop("driving pressure does not exist outside the ventilated cohort")
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
    left_join(base %>% select(hospitalization_id, pfvc_gli), by = "hospitalization_id") %>%
    group_by(hospitalization_id) %>% summarise(vtpbw_H = mean(tidal_volume_set / pbw),
                                               vtpfvc_H = mean(tidal_volume_set / pfvc_gli), .groups = "drop")
  d <- base %>%
    left_join(baseline, by = "hospitalization_id") %>%
    left_join(yH, by = "hospitalization_id") %>%
    left_join(dose, by = "hospitalization_id") %>%
    left_join(rrt, by = "hospitalization_id") %>%
    left_join(sf0, by = "hospitalization_id") %>%
    mutate(
      dead_before_H = !is.na(death_time_days) & death_time_days * 24 <= H,
      esc_before_H  = !is.na(escalation_time_days) & escalation_time_days * 24 <= H,   # control: intubated before H
      rrt_before_H  = !is.na(rrt_h) & rrt_h <= H,
      log_pfvc = log(pfvc_gli), ldisc = log(pbw / pfvc_gli),
      log_pfvc_sd = as.numeric(scale(log_pfvc)), ldisc_sd = as.numeric(scale(ldisc)),
      # VT/PFVC over the window in percent of predicted FVC, centred, per point: at a
      # given VT/PBW, the discordance contrast scaled by the dose, under its clinical name
      vtpfvc_c = if (HAS_DOSE) vtpfvc_H * 0.1 - median(vtpfvc_H * 0.1, na.rm = TRUE) else NA_real_,
      log_sf_0 = log(sf_0)
    )
  # the survival-side SF baseline and np_sofa come from base via the shared panel;
  # sf_0 is the index-day worst SF (day-0 value)
  if (MARKER == "ne_equiv") d <- b0 %>% select(hospitalization_id) %>% left_join(d, by = "hospitalization_id")
  n_all <- nrow(d)
  cc <- d %>% filter(!dead_before_H, !esc_before_H, !rrt_before_H | MARKER != "creatinine",
                     !is.na(yH), !is.na(y0), if (HAS_DOSE) !is.na(vtpbw_H) else TRUE, !is.na(np_sofa), !is.na(sf_0),
                     if (MARKER == "dp") !is.na(bmi) else TRUE,
                     !is.na(age10), !is.na(sex_category), !is.na(race_category))
  counts <- tibble(horizon_h = H, n_cohort = n_all, n_dead_before_H = sum(d$dead_before_H),
                   n_escalated_before_H = sum(d$esc_before_H),
                   n_rrt_before_H = sum(d$rrt_before_H),
                   n_no_baseline = sum(!d$dead_before_H & is.na(d$y0)),
                   n_no_outcome_value = sum(!d$dead_before_H & !is.na(d$y0) & is.na(d$yH)),
                   n_analysed = nrow(cc))
  message(sprintf("  H = %g h: %d in cohort, %d dead before H, %d on RRT before H, %d without baseline, %d without a value in the window, %d analysed",
                  H, n_all, counts$n_dead_before_H, counts$n_rrt_before_H, counts$n_no_baseline,
                  counts$n_no_outcome_value, nrow(cc)))
  if (nrow(cc) < 50) return(list(counts = counts, rows = NULL))

  base_terms <- c(if (HAS_DOSE) "vtpbw_H", "np_sofa", "log_sf_0", if (MARKER == "dp") "bmi")
  base_rhs <- paste(base_terms, collapse = " + ")
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
  # ---- channel decomposition: the exposure split into its height, age, sex and
  #      race pieces, each with its own coefficient (per log unit of the exposure),
  #      beside the constrained one-beta model on their sum; Wald test of equality.
  #      The pieces replace the demographic covariates, so there is one form only.
  chan_rows <- list()
  channels <- function(dat, lhs, expo_base, outcome_type, family = "gaussian", extra = "log_y0") {
    ch <- pfvc_channels(dat, expo_base)
    dd <- bind_cols(dat, ch)
    f_free <- as.formula(paste(lhs, "~", paste(CHANNELS, collapse = " + "), "+", extra, "+", base_rhs))
    f_one  <- as.formula(paste(lhs, "~ ch_sum +", extra, "+", base_rhs))
    fit  <- function(f) if (family == "gaussian") lm(f, data = dd) else glm(f, data = dd, family = binomial)
    free <- fit(f_free); one_fit <- fit(f_one)
    b <- coef(free)[CHANNELS]; V <- vcov(free)[CHANNELS, CHANNELS]
    co1 <- summary(one_fit)$coefficients["ch_sum", ]
    sd_expo <- sd(dd$ch_sum + dd$ch_remainder)
    bind_rows(
      tibble(channel = c("height", "age", "sex", "race"), estimate = unname(b), se = sqrt(diag(V))),
      tibble(channel = "all (one beta)", estimate = unname(co1[1]), se = unname(co1[2]))) %>%
      mutate(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se, per_sd = estimate * sd_expo,
             # the same coefficient read per 100 mL of PFVC at the median PFVC (d log PFVC = dPFVC / PFVC)
             per_100ml_at_median = if (expo_base == "log_pfvc") estimate * 0.1 / median(dd$pfvc_gli) else NA_real_,
             horizon_h = H, marker = MARKER, outcome_type = outcome_type,
             exposure = if (expo_base == "log_pfvc") "log PFVC" else "log PBW/PFVC (VT/PFVC at a given VT/PBW)",
             p_equal = channels_equal_p(b, V), n = nrow(dd),
             channel_sd = c(sapply(dd[CHANNELS], sd), sd_expo), remainder_sd = sd(dd$ch_remainder),
             worse_is = worse_is)
  }
  if (MARKER == "ne_equiv") {
    # actual weight from BMI and height (kg), for the absolute dose in mcg/min
    cc <- cc %>% mutate(any_H = as.integer(yH > 0), log_y0 = log(y0 + 0.01),
                        weight_kg = bmi * (height_cm / 100)^2)
    on <- cc %>% filter(yH > 0) %>% mutate(log_yH = log(yH), log_yH_abs = log(yH * weight_kg),
                                           log_y0_abs = log(y0 * weight_kg + 0.01))
    for (eb in c("log_pfvc", "ldisc"))
      chan_rows[[length(chan_rows) + 1]] <- channels(cc, "any_H", eb, "any pressor at H", "binomial")
    for (expo in c("log_pfvc_sd", "ldisc_sd", if (HAS_DOSE) "vtpfvc_c")) for (adj in c(TRUE, FALSE)) {
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
    for (expo in c("log_pfvc_sd", "ldisc_sd", if (HAS_DOSE) "vtpfvc_c")) for (adj in c(TRUE, FALSE))
      rows[[length(rows) + 1]] <- one(cc, "log_yH", expo, adj, "log marker at H")
    for (eb in c("log_pfvc", "ldisc"))
      chan_rows[[length(chan_rows) + 1]] <- channels(cc, "log_yH", eb, "log marker at H")
  }
  # ---- supports for the channel read: nested ladder, dose x piece, severity x
  #      piece (all on cc, the complete cases at H) and the baseline negative
  #      control (on d, everyone with a baseline, computed once)
  lhs   <- if (MARKER == "ne_equiv") "any_H" else "log_yH"
  fam   <- if (MARKER == "ne_equiv") "binomial" else "gaussian"
  otype <- if (MARKER == "ne_equiv") "any pressor at H" else "log marker at H"
  ccs <- bind_cols(cc, pfvc_channels(cc, "log_pfvc")) %>%
    mutate(vtpbw_H_c = if (HAS_DOSE) vtpbw_H - median(vtpbw_H) else NA_real_, log_sf_0_c = log_sf_0 - median(log_sf_0))
  fitf <- function(rhs, dat = ccs, y = lhs, family = fam) {
    f <- as.formula(paste(y, "~", rhs))
    if (family == "gaussian") lm(f, data = dat) else glm(f, data = dat, family = binomial)
  }
  # fit metric: R^2 for the linear models, the rank (Wilcoxon) AUC for the logistic
  fit_metric <- function(m) if (inherits(m, "glm")) {
    y <- m$y; r <- rank(fitted(m)); n1 <- sum(y == 1); n0 <- sum(y == 0)
    (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  } else summary(m)$r.squared
  coef_rows <- function(m, terms, model) {
    b <- coef(m)[terms]; V <- vcov(m)[terms, terms, drop = FALSE]
    tibble(model = model, term = terms, estimate = unname(b), se = sqrt(diag(V))) %>%
      mutate(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se, p = 2 * pnorm(-abs(estimate / se)))
  }
  # a term of a fitted model by its components, whatever order R put them in
  term_by <- function(m, comps) { nm <- names(coef(m)); nm[sapply(strsplit(nm, ":"), setequal, comps)] }

  # nested ladder (every LR pair exactly nested; ch_sum is the size term so that
  # one_beta sits inside pieces_free; pieces_demo omits ch_sex/ch_race, which
  # the sex and race factors already span)
  tail_rhs <- paste("log_y0 +", base_rhs)
  fits <- list(
    demo_only     = fitf(paste(demo_rhs, "+", tail_rhs)),
    one_beta      = fitf(paste("ch_sum +", tail_rhs)),
    one_beta_demo = fitf(paste("ch_sum +", demo_rhs, "+", tail_rhs)),
    pieces_free   = fitf(paste(paste(CHANNELS, collapse = " + "), "+", tail_rhs)),
    pieces_demo   = fitf(paste("ch_height + ch_age +", demo_rhs, "+", tail_rhs)),
    log_pfvc_demo = fitf(paste("log_pfvc +", demo_rhs, "+", tail_rhs)))
  lr_row <- function(small, big, test, nested = TRUE) {
    ll <- function(nm) logLik(fits[[nm]])
    lr <- if (nested) as.numeric(2 * (ll(big) - ll(small))) else NA_real_
    df <- if (nested) attr(ll(big), "df") - attr(ll(small), "df") else NA_real_
    tibble(test = test, small = small, big = big, lr = lr, df = df,
           p = if (nested) pchisq(lr, df, lower.tail = FALSE) else NA_real_,
           d_aic = AIC(fits[[big]]) - AIC(fits[[small]]), d_fit = fit_metric(fits[[big]]) - fit_metric(fits[[small]]),
           aic_small = AIC(fits[[small]]), aic_big = AIC(fits[[big]]))
  }
  nested <- bind_rows(
    lr_row("demo_only", "one_beta_demo", "a_vs_c: size adds to demographics"),
    lr_row("one_beta", "one_beta_demo", "b_vs_c: demographics add to size"),
    lr_row("one_beta", "pieces_free", "b_vs_d: the pieces disagree"),
    lr_row("pieces_free", "pieces_demo", "d_vs_top: free demographic shape beyond the pieces"),
    lr_row("one_beta_demo", "pieces_demo", "c_vs_top: pieces beyond one beta + demographics"),
    lr_row("one_beta_demo", "pieces_free", "d_vs_c: not nested, AIC only", nested = FALSE),
    imap_dfr(fits, ~ tibble(test = "model", small = NA_character_, big = .y, lr = NA_real_, df = NA_real_,
                            p = NA_real_, d_aic = NA_real_, d_fit = fit_metric(.x), aic_small = NA_real_, aic_big = AIC(.x)))) %>%
    mutate(horizon_h = H, marker = MARKER, outcome_type = otype, fit_metric = if (fam == "gaussian") "r2" else "auc",
           n = nrow(ccs), remainder_sd = sd(ccs$ch_remainder))

  # dose x piece and severity x piece: the pieces' effects as a function of the
  # delivered VT/PBW over [0, H) and of the baseline SF (each centred at its median)
  interactions <- function(modifier, source_term, label) {
    base_here <- paste(setdiff(base_terms, source_term), collapse = " + ")
    free <- fitf(paste0("(", paste(CHANNELS, collapse = " + "), ") * ", modifier, " + log_y0 + ", base_here))
    one  <- fitf(paste0("log_pfvc * ", modifier, " + ", demo_rhs, " + log_y0 + ", base_here))
    it   <- sapply(CHANNELS, function(ch) term_by(free, c(ch, modifier)))
    b <- coef(free)[it]; V <- vcov(free)[it, it]
    hm <- b[1] - b[2]; hm_se <- sqrt(V[1, 1] + V[2, 2] - 2 * V[1, 2])
    bind_rows(
      coef_rows(free, unname(it), "pieces_free") %>% mutate(p_int_equal = channels_equal_p(b, V)),
      tibble(model = "pieces_free", term = "height_minus_age", estimate = unname(hm), se = hm_se,
             lo = hm - 1.96 * hm_se, hi = hm + 1.96 * hm_se, p = 2 * pnorm(-abs(hm / hm_se)), p_int_equal = NA_real_),
      coef_rows(one, term_by(one, c("log_pfvc", modifier)), "single_index_adjusted") %>% mutate(p_int_equal = NA_real_)) %>%
      mutate(horizon_h = H, marker = MARKER, outcome_type = otype, modifier = label,
             modifier_median = median(ccs[[source_term]]), modifier_sd = sd(ccs[[source_term]]), n = nrow(ccs))
  }
  dose_ch <- if (HAS_DOSE) interactions("vtpbw_H_c", "vtpbw_H", "mean VT/PBW over [0, H), mL/kg, centred") else
    tibble(horizon_h = H, marker = MARKER, note = "skipped: no ventilator dose outside the ventilated cohort")
  sf_ch   <- if (MARKER == "sf") tibble(horizon_h = H, marker = MARKER, note = "skipped: SF is this marker's own baseline") else
    interactions("log_sf_0_c", "log_sf_0", "log baseline SF, centred")

  # baseline negative control, once: the pieces on the baseline marker, on everyone
  # with a baseline (not the survivors at H); no dose term (nothing delivered yet)
  negctrl <- NULL
  if (H == HORIZONS[1]) {
    nc <- d %>% filter(!is.na(y0), !is.na(np_sofa), !is.na(sf_0), !is.na(age10), !is.na(sex_category),
                       !is.na(race_category), !is.na(height_cm), if (MARKER == "dp") !is.na(bmi) else TRUE)
    nc <- bind_cols(nc, pfvc_channels(nc, "log_pfvc"))
    nc_rhs <- paste(c("np_sofa", if (MARKER != "sf") "log_sf_0", if (MARKER == "dp") "bmi"), collapse = " + ")
    if (MARKER == "ne_equiv") {
      nc$any_0 <- as.integer(nc$y0 > 0); nc_lhs <- "any_0"; nc_fam <- "binomial"; nc_out <- "any pressor at baseline"
    } else { nc$log_y0 <- log(nc$y0); nc_lhs <- "log_y0"; nc_fam <- "gaussian"; nc_out <- "log baseline marker" }
    nfit <- list(
      pieces_free           = fitf(paste(paste(CHANNELS, collapse = " + "), "+", nc_rhs), nc, nc_lhs, nc_fam),
      single_index_adjusted = fitf(paste("log_pfvc +", nc_rhs, "+", demo_rhs), nc, nc_lhs, nc_fam),
      pieces_only           = fitf(paste(CHANNELS, collapse = " + "), nc, nc_lhs, nc_fam))
    negctrl <- imap_dfr(nfit, function(m, nm) {
      terms <- if (nm == "single_index_adjusted") "log_pfvc" else CHANNELS
      coef_rows(m, terms, nm) %>%
        mutate(p_equal = if (nm == "single_index_adjusted") NA_real_ else
                 channels_equal_p(coef(m)[CHANNELS], vcov(m)[CHANNELS, CHANNELS]),
               separation_flag = any(abs(coef(m)) > 15, na.rm = TRUE))
    }) %>% mutate(marker = MARKER, outcome = nc_out, n = nrow(nc),
                  n_events = if (nc_fam == "binomial") sum(nc[[nc_lhs]]) else NA_integer_)
  }

  # composite-rank sensitivity: death before H worst, RRT before H next, then the marker
  # (worse direction first), on everyone with a baseline; rank scaled to (0, 1)
  comp <- d %>% filter(!is.na(y0), if (HAS_DOSE) !is.na(vtpbw_H) else TRUE, !esc_before_H, !is.na(np_sofa), !is.na(sf_0),
                       if (MARKER == "dp") !is.na(bmi) else TRUE,
                       !is.na(age10), !is.na(sex_category), !is.na(race_category)) %>%
    mutate(score = case_when(dead_before_H ~ Inf,
                             rrt_before_H & MARKER == "creatinine" ~ 1e9,
                             is.na(yH) ~ NA_real_,
                             worse_is == "higher" ~ yH, TRUE ~ -yH)) %>%
    filter(!is.na(score)) %>%
    mutate(rank01 = (rank(score) - 0.5) / n(), log_y0 = if (MARKER == "ne_equiv") log(y0 + 0.01) else log(y0))
  if (nrow(comp) >= 50)
    for (expo in c("log_pfvc_sd", "ldisc_sd", if (HAS_DOSE) "vtpfvc_c")) for (adj in c(TRUE, FALSE))
      rows[[length(rows) + 1]] <- one(comp, "rank01", expo, adj, "composite rank (death, RRT, marker)") %>%
        mutate(note = sprintf("rank 0-1; %d deaths and %d RRT before H ranked worst",
                                              sum(comp$dead_before_H), sum(comp$rrt_before_H & MARKER == "creatinine")))
  list(counts = counts, rows = bind_rows(rows), chan = bind_rows(chan_rows),
       nested = nested, dose_ch = dose_ch, sf_ch = sf_ch, negctrl = negctrl)
}

res <- map(HORIZONS, fit_horizon)
counts  <- map_dfr(res, "counts") %>% mutate(marker = MARKER, site = site_name)
results <- map_dfr(res, "rows") %>% mutate(site = site_name)
chan    <- map_dfr(res, "chan") %>% mutate(site = site_name)
nested  <- map_dfr(res, "nested")  %>% mutate(site = site_name)
dose_ch <- map_dfr(res, "dose_ch") %>% mutate(site = site_name)
sf_ch   <- map_dfr(res, "sf_ch")   %>% mutate(site = site_name)
negctrl <- map_dfr(res, "negctrl") %>% mutate(site = site_name)
if (nrow(results) == 0) stop("no horizon had enough patients to fit")

write_csv(results, file.path(final_dir, paste0("injury_at_horizon_", MARKER, "_", site_name, ".csv")))
write_csv(counts,  file.path(final_dir, paste0("injury_at_horizon_counts_", MARKER, "_", site_name, ".csv")))
write_csv(chan,    file.path(final_dir, paste0("injury_channels_", MARKER, "_", site_name, ".csv")))
write_csv(nested,  file.path(final_dir, paste0("injury_nested_", MARKER, "_", site_name, ".csv")))
write_csv(dose_ch, file.path(final_dir, paste0("injury_dose_channels_", MARKER, "_", site_name, ".csv")))
write_csv(sf_ch,   file.path(final_dir, paste0("injury_sf_channels_", MARKER, "_", site_name, ".csv")))
write_csv(negctrl, file.path(final_dir, paste0("injury_negctrl_", MARKER, "_", site_name, ".csv")))
message("\n--- nested ladder at ", HORIZONS[1], " h (LR tests; d_aic < 0 favours the bigger model)")
print(as.data.frame(nested %>% filter(horizon_h == HORIZONS[1], test != "model") %>%
                      select(test, lr, df, p, d_aic, d_fit) %>% mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
if (HAS_DOSE) {
  message("\n--- dose x piece at ", HORIZONS[1], " h (per log unit of the piece per mL/kg)")
  print(as.data.frame(dose_ch %>% filter(horizon_h == HORIZONS[1]) %>% select(model, term, estimate, lo, hi, p, p_int_equal) %>%
                        mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
}
message("\n--- baseline negative control (", unique(negctrl$outcome), "): the pieces' direct effects")
print(as.data.frame(negctrl %>% select(model, term, estimate, lo, hi, p, p_equal, separation_flag) %>%
                      mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)

# supports figure: ladder (dAIC), dose x piece, severity x piece, negative control beside the at-H pieces
hz <- function(x) factor(paste0(x, " h"), paste0(sort(unique(x)), " h"))
p_lad <- ggplot(nested %>% filter(test != "model") %>% mutate(horizon = hz(horizon_h)),
                aes(d_aic, test, fill = d_aic < 0)) +
  geom_col() + facet_wrap(~ horizon, nrow = 1) +
  scale_fill_manual(values = okabe[c(4, 3)], guide = "none") +
  labs(title = "Nested ladder: change in AIC (negative favours the bigger model)", x = "dAIC (big - small)", y = NULL)
int_plot <- function(dd, title) ggplot(dd %>% mutate(horizon = hz(horizon_h), term = sub(":.*$", "", term)),
                                       aes(estimate, term, colour = model)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5)) +
  facet_wrap(~ horizon, nrow = 1, scales = "free_x") + scale_colour_manual(values = okabe[1:2], name = NULL) +
  labs(title = title, x = "interaction per log unit of the piece per unit of the modifier", y = NULL)
p_dose <- if ("estimate" %in% names(dose_ch)) int_plot(dose_ch, "Dose x piece: does the piece's effect scale with the delivered VT/PBW?") else plot_spacer()
p_sf   <- if ("estimate" %in% names(sf_ch)) int_plot(sf_ch, "Severity x piece: does it scale with baseline SF (baby lung)?") else plot_spacer()
nc_plot <- bind_rows(
  negctrl %>% filter(model == "pieces_free") %>% transmute(term, estimate, lo, hi, when = "baseline (negative control)"),
  chan %>% filter(exposure == "log PFVC", horizon_h == HORIZONS[1], channel != "all (one beta)") %>%
    transmute(term = paste0("ch_", channel), estimate, lo, hi, when = paste0("at ", HORIZONS[1], " h")))
p_nc <- ggplot(nc_plot, aes(estimate, term, colour = when)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = okabe[c(4, 1)], name = NULL) +
  labs(title = "Negative control: the pieces on the baseline marker vs at the horizon", x = "per log unit of the piece", y = NULL)
p_sup <- (p_lad / p_dose / p_sf / p_nc) + plot_layout(heights = c(1.2, 1, 1, 1)) +
  plot_annotation(title = sprintf("%s: supports for the channel read (%s)", MARKER, site_name)) & theme_minimal(base_size = 10)
ggsave(file.path(final_dir, paste0("injury_supports_", MARKER, "_", site_name, ".pdf")), p_sup, width = 11, height = 14)
message("\n--- channel decomposition: the exposure effect per log unit, identified through each GLI input",
        " (equal coefficients = lung size is the operative quantity; p_equal tests that)")
print(as.data.frame(chan %>% filter(horizon_h == HORIZONS[1]) %>%
                      select(exposure, channel, estimate, lo, hi, per_sd, per_100ml_at_median, p_equal, n) %>%
                      mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
message("remainder of the decomposition (SD, log units): ", signif(max(chan$remainder_sd), 2))
pc <- chan %>%
  mutate(channel = factor(channel, c("height", "age", "sex", "race", "all (one beta)")),
         horizon = factor(paste0(horizon_h, " h"), paste0(sort(unique(horizon_h)), " h")),
         strip = sprintf("%s\np(equal) = %.2g", exposure, p_equal))
p_ch <- ggplot(pc, aes(estimate, channel, colour = channel == "all (one beta)")) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_vline(data = pc %>% filter(channel == "all (one beta)"), aes(xintercept = estimate), colour = okabe[2]) +
  geom_pointrange(aes(xmin = lo, xmax = hi)) +
  facet_wrap(~ strip + horizon, scales = "free_x", ncol = n_distinct(pc$horizon), dir = "h") +
  scale_colour_manual(values = okabe[c(1, 2)], guide = "none") +
  labs(title = sprintf("%s at the horizon: the size effect identified through each input to PFVC (%s)", MARKER, site_name),
       subtitle = "one coefficient per GLI piece (log units of the exposure); the line is the one-beta model on their sum",
       x = "Change per log unit of the exposure (log marker; log-odds for any pressor)", y = NULL) +
  theme_minimal(base_size = 11)
ggsave(file.path(final_dir, paste0("injury_channels_", MARKER, "_", site_name, ".pdf")), p_ch, width = 11, height = 8)
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
