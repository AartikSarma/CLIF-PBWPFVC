# =============================================================================
# Script 13 (7-day feasibility): can the window be extended, and what shape?
# =============================================================================
# Before spending hours on a 7-day joint model, four cheap questions, all from
# the daily-grid 7-day panel and a handful of mixed models (minutes, no MCMC):
#
#   1. Is there a cohort left?  Patients still ventilated on each day, with the
#      deaths and extubations that removed the rest, and the marker rows they
#      still contribute.  A 7-day window on a cohort that has mostly extubated
#      by day 4 estimates a trajectory in the sick tail, not in the cohort.
#   2. Is the exposure still the exposure?  VT/PBW and VT/PFVC by day, and the
#      share of patient-days that leave the 5-8 mL/kg band.  As patients move to
#      pressure support the tidal volume stops being a dose the clinician sets,
#      and "VT/PFVC at a given VT/PBW" stops meaning what it means on day 1.
#   3. Does the dropout differ by exposure?  Cumulative extubation and death by
#      baseline VT/PFVC tercile.  The SF crossover this script exists to look
#      for (an early oxygenation advantage that erodes) is exactly what
#      differential extubation manufactures: if the patients with a large lung
#      for the dose leave first, the ones left behind in that tercile are its
#      sickest, and their marker converges toward the others for no physiologic
#      reason.  This is the diagnostic that decides whether a 7-day estimate
#      needs the joint model's competing-risk correction (it will).
#   4. What functional form does the marker difference need?  The contrast is
#      fitted three ways over the same rows -- linear in time (the current
#      level + divergence x t), a natural spline in time (3 df), and a free
#      contrast per day -- and compared by AIC.  The free-per-day fit is the
#      assumption-light picture of the shape; if a line fits it, the joint model
#      can stay linear and cheap.  Each is fitted twice: on everyone at risk
#      (the estimand, contaminated by dropout) and on the balanced subset still
#      ventilated at day 7 (dropout-free, but a selected and sicker cohort).
#      Neither alone is the answer; the gap between them is the dropout signal.
#
# Needs the 7-day daily-grid panel:
#   PBWPFVC_JM_GRID=daily PBWPFVC_JM_HORIZON=7 Rscript code/13_biotrauma_panel.R
# Usage:
#   PBWPFVC_JM_GRID=daily PBWPFVC_JM_HORIZON=7 PBWPFVC_INJ_MARKER=sf \
#     Rscript code/13_seven_day_feasibility.R
#
# Writes to final/ (aggregates only, every reported cell >= MIN_CELL patients):
#   sevenday_riskset_{marker}_{site}.csv   day, at risk, deaths, extubations,
#       marker rows, VT/PBW and VT/PFVC quartiles, share outside the band
#   sevenday_dropout_{marker}_{site}.csv   cumulative extubation/death by day and
#       baseline VT/PFVC tercile, plus the drift in baseline covariates between
#       the day-2 cohort and the day-7 survivors
#   sevenday_shape_{marker}_{site}.csv     the contrast per day (free), the
#       linear and spline fits, and the AIC comparison, in both samples
#   sevenday_feasibility_{marker}_{site}.pdf
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(here)
                                 library(splines); library(nlme); library(patchwork) })
rm(list = ls())
source("utils/config.R")
site_name  <- config$site_name
output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)

# the 7-day daily grid is the point of the script; a finer grid costs time and
# buys nothing for a question about days 3 to 7
if (!nzchar(Sys.getenv("PBWPFVC_JM_GRID")))    Sys.setenv(PBWPFVC_JM_GRID = "daily")
if (!nzchar(Sys.getenv("PBWPFVC_JM_HORIZON"))) Sys.setenv(PBWPFVC_JM_HORIZON = "7")
source(here("code", "13_biotrauma_grid.R"))

MARKER   <- Sys.getenv("PBWPFVC_INJ_MARKER", "sf")
MIN_CELL <- 10L                       # no reported group smaller than this
BAND     <- c(5, 8)                   # the VT/PBW inclusion band
okabe    <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#000000")
theme_set(theme_minimal(base_size = 10))

y_col  <- c(creatinine = "creatinine", ne_equiv = "ne_equiv_peak", platelets = "platelets",
            bilirubin = "bilirubin", sf = "sf", dp = "dp")[[MARKER]]
y0_col <- c(creatinine = "creatinine_0", ne_equiv = "ne_equiv_0", platelets = "platelet_0",
            bilirubin = "bilirubin_0", sf = "sf_0", dp = "dp_0")[[MARKER]]
offset <- if (MARKER == "ne_equiv") 0.01 else 0
marker_lab <- c(creatinine = "creatinine", platelets = "platelets", sf = "SF ratio",
                ne_equiv = "norepinephrine equivalent", bilirubin = "bilirubin", dp = "driving pressure")[[MARKER]]

panel_path <- file.path(output_dir, paste0("jm_long_", h_suffix, ".parquet"))
if (!file.exists(panel_path))
  stop("no ", h_suffix, " panel: run  PBWPFVC_JM_GRID=daily PBWPFVC_JM_HORIZON=", JM_HORIZON,
       " Rscript code/13_biotrauma_panel.R")
long <- read_parquet(panel_path)
surv <- read_parquet(file.path(output_dir, paste0("jm_surv_", h_suffix, ".parquet")))
message("=== 13_seven_day_feasibility (", MARKER, ", ", h_suffix, " panel): ",
        nrow(surv), " patients, ", nrow(long), " patient-periods ===")

# baseline VT/PFVC tercile: the INDEX day's value, fixed for the whole window, so
# nothing in the grouping is a function of what happened later
surv <- surv %>%
  mutate(vtpfvc_base = vtpfvc_0,
         tercile = if (sum(is.finite(vtpfvc_base)) >= 3 * MIN_CELL)
           cut(vtpfvc_base, breaks = quantile(vtpfvc_base, c(0, 1/3, 2/3, 1), na.rm = TRUE),
               include.lowest = TRUE, labels = c("low VT/PFVC", "middle", "high VT/PFVC")) else NA)

# =============================================================================
# 1. the risk set, the marker rows and the exposure, day by day
# =============================================================================
days <- seq(0, JM_HORIZON, by = STEP)
riskset <- map_dfr(days, function(dd) {
  at_risk <- surv %>% filter(event_day >= dd)   # not yet dead, not yet extubated
  rows    <- long %>% filter(abs(vent_day - dd) < STEP / 2)
  obs     <- rows %>% filter(!is.na(.data[[y_col]]))
  vtpbw   <- rows$vt_ml / surv$pbw[match(rows$hospitalization_id, surv$hospitalization_id)]
  tibble(day = dd,
         n_at_risk        = nrow(at_risk),
         n_deaths_cum     = sum(surv$event == 1L & surv$event_day <= dd),
         n_extub_cum      = sum(surv$event == 2L & surv$event_day <= dd),
         n_marker_rows    = nrow(obs),
         n_marker_patients = n_distinct(obs$hospitalization_id),
         marker_coverage  = if (nrow(at_risk)) n_distinct(obs$hospitalization_id) / nrow(at_risk) else NA_real_,
         vtpbw_q25 = quantile(vtpbw, 0.25, na.rm = TRUE), vtpbw_median = median(vtpbw, na.rm = TRUE),
         vtpbw_q75 = quantile(vtpbw, 0.75, na.rm = TRUE),
         frac_outside_band = mean(vtpbw < BAND[1] | vtpbw > BAND[2], na.rm = TRUE),
         vtpfvc_q25 = quantile(rows$vtpfvc, 0.25, na.rm = TRUE),
         vtpfvc_median = median(rows$vtpfvc, na.rm = TRUE),
         vtpfvc_q75 = quantile(rows$vtpfvc, 0.75, na.rm = TRUE))
}) %>%
  mutate(across(c(vtpbw_q25:vtpfvc_q75), ~ if_else(n_marker_patients >= MIN_CELL, ., NA_real_)),
         marker = MARKER, site = site_name)
write_csv(riskset, file.path(final_dir, paste0("sevenday_riskset_", MARKER, "_", site_name, ".csv")))

# =============================================================================
# 2. differential dropout: cumulative extubation and death by baseline tercile
# =============================================================================
dropout <- if (all(is.na(surv$tercile))) tibble() else
  crossing(day = days, tercile = levels(surv$tercile)) %>%
  rowwise() %>%
  mutate(n_tercile   = sum(surv$tercile == tercile, na.rm = TRUE),
         cum_extub   = sum(surv$tercile == tercile & surv$event == 2L & surv$event_day <= day, na.rm = TRUE),
         cum_death   = sum(surv$tercile == tercile & surv$event == 1L & surv$event_day <= day, na.rm = TRUE),
         n_still_vent = sum(surv$tercile == tercile & surv$event_day >= day, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(pct_extub = 100 * cum_extub / n_tercile, pct_death = 100 * cum_death / n_tercile,
         across(c(cum_extub, cum_death, pct_extub, pct_death),
                ~ if_else(n_still_vent >= MIN_CELL | day == 0, ., NA_real_)),
         marker = MARKER, site = site_name)

# selection drift: the cohort at day 2 against the patients still ventilated at day 7
drift <- {
  vars <- c("age10", "np_sofa", "bmi", "height_cm", "vtpfvc_base", "log_pfvc_sd", "ldisc_sd", "sf_0")
  vars <- intersect(vars, names(surv))
  a <- surv %>% filter(event_day >= 2)
  b <- surv %>% filter(event_day >= JM_HORIZON)
  if (nrow(b) < MIN_CELL) tibble() else
    map_dfr(vars, function(v) {
      x <- a[[v]]; y <- b[[v]]
      s <- sqrt((var(x, na.rm = TRUE) + var(y, na.rm = TRUE)) / 2)
      tibble(variable = v, mean_day2 = mean(x, na.rm = TRUE), mean_day7_survivors = mean(y, na.rm = TRUE),
             smd = if (is.finite(s) && s > 0) (mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE)) / s else NA_real_)
    }) %>% mutate(n_day2 = nrow(a), n_day7 = nrow(b), marker = MARKER, site = site_name)
}
write_csv(bind_rows(dropout %>% mutate(table = "cumulative_dropout"),
                    drift   %>% mutate(table = "selection_drift")),
          file.path(final_dir, paste0("sevenday_dropout_", MARKER, "_", site_name, ".csv")))

# =============================================================================
# 3. the shape of the contrast: free per day, spline, linear; two samples
# =============================================================================
d_all <- long %>%
  filter(period >= 1L, !is.na(.data[[y_col]]), !is.na(l_sf), !is.na(l_pressor)) %>%
  inner_join(surv %>% select(hospitalization_id, np_sofa, bmi, age10, sex_category, race_category,
                             vtpfvc_c, vtpbw_pt_mean, event_day, tercile, all_of(y0_col)),
             by = "hospitalization_id") %>%
  filter(!is.na(.data[[y0_col]]), !is.na(np_sofa), !is.na(vtpfvc_c)) %>%
  mutate(log_y = log(.data[[y_col]] + offset), log_y0 = log(.data[[y0_col]] + offset),
         l_log_sf = log(l_sf), day_f = factor(round(vent_day / STEP) * STEP))

# the covariates the quick LME uses, minus the marker's own lag when it is the marker
lag_terms <- setdiff(c("l_log_sf", "l_pressor"), c(sf = "l_log_sf", ne_equiv = "l_pressor")[MARKER])
base_rhs  <- paste(c("log_y0", "np_sofa", lag_terms, if (MARKER == "dp") "bmi",
                     "vtpbw_pt_mean", "l_vtpbw_within"), collapse = " + ")
ctrl <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200, returnObject = TRUE)

fit_shape <- function(dd, sample_lab) {
  if (n_distinct(dd$hospitalization_id) < 3 * MIN_CELL) return(tibble())
  forms <- list(
    linear = paste("log_y ~ vent_day * vtpfvc_c +", base_rhs),
    spline = paste("log_y ~ ns(vent_day, 3) * vtpfvc_c +", base_rhs),
    free   = paste("log_y ~ day_f * vtpfvc_c +", base_rhs))
  fits <- imap(forms, function(f, nm) {
    try(lme(as.formula(f), random = list(hospitalization_id = pdDiag(~ vent_day)),
            data = dd, control = ctrl, method = "ML"), silent = TRUE)
  })
  ok <- !map_lgl(fits, ~ inherits(.x, "try-error"))
  aic <- tibble(sample = sample_lab, form = names(forms),
                aic = map_dbl(fits, ~ if (inherits(.x, "try-error")) NA_real_ else AIC(.x)),
                converged = ok, kind = "aic")
  # the free fit's contrast per day: the exposure main effect plus that day's interaction
  contrasts <- if (!ok[["free"]]) tibble() else {
    f <- fits[["free"]]; b <- fixef(f); V <- vcov(f)
    lv <- levels(dd$day_f)
    map_dfr(lv, function(l) {
      tm <- c("vtpfvc_c", paste0("day_f", l, ":vtpfvc_c"), paste0("vtpfvc_c:day_f", l))
      tm <- intersect(tm, names(b))
      k <- as.numeric(names(b) %in% tm)
      est <- sum(b[tm]); se <- sqrt(as.numeric(t(k) %*% V %*% k))
      n_d <- n_distinct(dd$hospitalization_id[dd$day_f == l])
      tibble(sample = sample_lab, form = "free", kind = "contrast", day = as.numeric(l),
             estimate = est, se = se, lo = est - 1.96 * se, hi = est + 1.96 * se,
             n_patients = n_d)
    }) %>% mutate(across(c(estimate, se, lo, hi), ~ if_else(n_patients >= MIN_CELL, ., NA_real_)))
  }
  # the linear fit's level and divergence, for the record
  lin <- if (!ok[["linear"]]) tibble() else {
    b <- fixef(fits[["linear"]]); V <- vcov(fits[["linear"]])
    tn <- intersect(c("vent_day:vtpfvc_c", "vtpfvc_c:vent_day"), names(b))
    tibble(sample = sample_lab, form = "linear", kind = "term",
           term = c("level", "divergence_per_day"),
           estimate = c(b[["vtpfvc_c"]], b[[tn[1]]]),
           se = c(sqrt(V["vtpfvc_c", "vtpfvc_c"]), sqrt(V[tn[1], tn[1]]))) %>%
      mutate(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se,
             n_patients = n_distinct(dd$hospitalization_id))
  }
  bind_rows(aic, contrasts, lin)
}

shape <- bind_rows(
  fit_shape(d_all, "all at risk"),
  fit_shape(d_all %>% filter(event_day >= JM_HORIZON), paste0("ventilated through day ", JM_HORIZON))) %>%
  mutate(marker = MARKER, horizon_days = JM_HORIZON, unit = "per point of VT/PFVC, log marker",
         site = site_name)
write_csv(shape, file.path(final_dir, paste0("sevenday_shape_", MARKER, "_", site_name, ".csv")))

# =============================================================================
# 4. the figure
# =============================================================================
p1 <- riskset %>%
  select(day, `still ventilated` = n_at_risk, extubated = n_extub_cum, died = n_deaths_cum) %>%
  pivot_longer(-day) %>%
  ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1) + geom_point() +
  scale_colour_manual(values = okabe[c(5, 3, 1)], name = NULL) +
  labs(title = "Who is left to contribute a trajectory", x = "ventilator day", y = "patients")

p2 <- riskset %>%
  ggplot(aes(day, vtpbw_median)) +
  geom_ribbon(aes(ymin = vtpbw_q25, ymax = vtpbw_q75), alpha = 0.18, fill = okabe[4]) +
  geom_line(colour = okabe[4], linewidth = 1) +
  geom_hline(yintercept = BAND, linetype = 2, colour = "grey55") +
  labs(title = "Is the dose still in the band?", x = "ventilator day",
       y = "VT/PBW (mL/kg), median and IQR")

p3 <- if (!nrow(dropout)) plot_spacer() else
  dropout %>% filter(!is.na(pct_extub)) %>%
  ggplot(aes(day, pct_extub, colour = tercile)) + geom_line(linewidth = 1) + geom_point() +
  scale_colour_manual(values = okabe[c(4, 6, 5)], name = NULL) +
  labs(title = "Differential dropout: cumulative extubation by baseline VT/PFVC",
       subtitle = "terciles that separate here manufacture a crossover in any survivor-only estimate",
       x = "ventilator day", y = "% extubated")

p4 <- {
  cc <- shape %>% filter(kind == "contrast", !is.na(estimate))
  if (!nrow(cc)) plot_spacer() else
    ggplot(cc, aes(day, estimate, colour = sample, fill = sample)) +
    geom_hline(yintercept = 0, colour = "grey55") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) + geom_point() +
    scale_colour_manual(values = okabe[c(1, 2)], name = NULL) +
    scale_fill_manual(values = okabe[c(1, 2)], name = NULL) +
    labs(title = paste0("Marker difference per point of VT/PFVC, free per day (", marker_lab, ")"),
         subtitle = "a gap between the two samples is dropout, not physiology",
         x = "ventilator day", y = "log marker per point")
}

pdf_path <- file.path(final_dir, paste0("sevenday_feasibility_", MARKER, "_", site_name, ".pdf"))
ggsave(pdf_path, (p1 + p2) / (p3 + p4) +
         plot_annotation(title = paste0(site_name, ": can the window go to ", JM_HORIZON, " days, and with what shape?"),
                         subtitle = paste0(marker_lab, "; groups under ", MIN_CELL, " patients suppressed")) &
         theme(legend.position = "top"),
       width = 12, height = 8)

# =============================================================================
message("\nRisk set by day:")
print(as.data.frame(riskset %>% select(day, n_at_risk, n_deaths_cum, n_extub_cum, n_marker_patients,
                                       marker_coverage, vtpbw_median, frac_outside_band, vtpfvc_median) %>%
                      mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
if (nrow(drift)) {
  message("\nSelection drift, day-2 cohort vs patients still ventilated at day ", JM_HORIZON,
          " (|SMD| > 0.1 = the 7-day cohort is a different cohort):")
  print(as.data.frame(drift %>% select(variable, mean_day2, mean_day7_survivors, smd) %>%
                        mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
}
message("\nShape: AIC by functional form (lower is better; 'free' beating 'linear' by more than ~10 means the line misses the shape):")
print(as.data.frame(shape %>% filter(kind == "aic") %>% select(sample, form, aic, converged) %>%
                      mutate(aic = signif(aic, 6))), row.names = FALSE)
message("\nContrast per day (per point of VT/PFVC; the crossover, if it exists, is where this changes sign):")
print(as.data.frame(shape %>% filter(kind == "contrast") %>%
                      select(sample, day, estimate, lo, hi, n_patients) %>%
                      mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
message("\n13_seven_day_feasibility complete -> ", final_dir)
