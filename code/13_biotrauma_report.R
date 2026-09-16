# =============================================================================
# Script 13 (report): Biotrauma joint models -- trajectories, associations, figures
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# Reads the fit bundles 13_biotrauma_fit.R saved and writes the aggregate,
# poolable reads of the three questions plus the figures:
#
#   final/jm_trajectory_grid_{H}d_{site}.csv   Q1 as a picture: the posterior
#       population trajectory of each log marker over days 1..H under a constant
#       previous-day VT/PFVC of 9, 11, 14 and 17% (the ARMA-anchored strain
#       levels), relative to day 1 at 11%, with 95% credible bands, on a fixed grid
#   final/jm_strain_effects_{H}d_{site}.csv    Q1 as numbers: the previous-day
#       strain and cumulative-days coefficients per marker, adjusted and
#       unadjusted, per unit and per SD of the log marker
#   final/jm_association_hr_{H}d_{site}.csv    Q2: hazard ratio for death and for
#       extubation per SD of the current log marker (value) and per unit slope
#   final/jm_heterogeneity_{H}d_{site}.csv     the strain slope at the 10th, 50th
#       and 90th percentile of baseline specific elastance (Ers x PFVC)
#   final/jm_absorption_{H}d_{site}.csv        Q3 (written by the fit script;
#       re-read here for the figure)
#   final/jm_trajectories_{H}d_{site}.pdf, jm_forest_{H}d_{site}.pdf
#
# Population trajectories come from the fixed-effects design matrix times the
# joint posterior draws of the longitudinal betas (the 05e route: predict.jm
# fails on competing-risk fits), so the bands carry the JM's joint uncertainty.
#
# Usage: Rscript code/13_biotrauma_report.R     (PBWPFVC_JM_HORIZON, PBWPFVC_JM_BASELINE as in the fit)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(splines)
  library(nlme)
  library(survival)
  library(JMbayes2)
  library(GLMMadaptive)   # fixef() for the logistic mixed model of the any-pressor fit
  library(patchwork)
})
rm(list = ls())
source("utils/config.R")

site_name  <- config$site_name
output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
source(here("code", "13_biotrauma_grid.R"))   # JM_GRID, STEP, JM_HORIZON, N_PERIODS, h_suffix
BASELINE_FORM <- Sys.getenv("PBWPFVC_JM_BASELINE", "free")
MOD_FORM      <- Sys.getenv("PBWPFVC_JM_MODIFIER", "disc")
stopifnot(MOD_FORM %in% c("disc", "saturated", "none", "pfvc", "disc_level"))
out_tag  <- paste0(if (BASELINE_FORM == "offset") "offset_" else "",
                   if (MOD_FORM != "disc") paste0(MOD_FORM, "_") else "",
                   h_suffix, "_", site_name)
okabe <- c("#009E73", "#56B4E9", "#E69F00", "#D55E00", "#0072B2", "#CC79A7")
DOSE_LEVELS    <- c(6, 8, 10)     # VT/PBW, mL/kg: the LTVV target, its upper bound, conventional
DOSE_REF       <- 6
DISC_PCT       <- c(10, 50, 90)   # PBW/PFVC discordance percentiles for the facets
N_DRAWS <- 1000L
set.seed(20260913)

manifest <- read_csv(file.path(final_dir, paste0("jm_manifest_", out_tag, ".csv")), show_col_types = FALSE)
usable <- manifest %>% filter(status %in% c("converged", "rhat_fail"))
if (nrow(usable) == 0L) stop("No fitted joint models in the manifest for ", out_tag)
message("=== 13_biotrauma_report (", out_tag, "): ", nrow(usable), " fits, of which ",
        sum(usable$status == "converged"), " pass the R-hat gate ===")

# --- posterior draws of the longitudinal fixed effects, stacked across chains
beta_draws <- function(jm, lme_fit) {
  b <- do.call(rbind, jm$mcmc$betas1)
  ref <- names(fixef(lme_fit))
  if (is.null(colnames(b))) colnames(b) <- ref
  stopifnot(all(ref %in% colnames(b)))
  b[, ref, drop = FALSE]
}
# population design row: continuous covariates at their patient-level median,
# factors at the reference level; the strain terms are set by the caller
population_row <- function(ld) {
  pt <- ld %>% distinct(hospitalization_id, .keep_all = TRUE)
  row <- tibble(
    np_sofa = median(pt$np_sofa), bmi = median(pt$bmi), age10 = median(pt$age10),
    age10_c = 0,   # the dose x age interaction is centred at the median age
    sex_category  = factor(levels(factor(ld$sex_category))[1],  levels = levels(factor(ld$sex_category))),
    race_category = factor(levels(factor(ld$race_category))[1], levels = levels(factor(ld$race_category))),
    l_log_sf = median(ld$l_log_sf), l_pressor = 0)
  if ("log_y0" %in% names(ld))     row$log_y0     <- median(pt$log_y0)
  if ("ers_pfvc_0" %in% names(ld)) row$ers_pfvc_0 <- median(pt$ers_pfvc_0, na.rm = TRUE)
  row
}

trajectory_rows <- list(); strain_rows <- list(); assoc_rows <- list(); hetero_rows <- list()
level_rows <- list()
LEVEL_HOURS <- c(24, 48, 72)
traj_plots <- list()

for (i in seq_len(nrow(usable))) {
  u <- usable[i, ]
  tag <- paste(u$marker, u$model, u$adjustment, sep = "_")
  f <- file.path(output_dir, paste0("jm_fit_", tag, "_", BASELINE_FORM,
                                    if (MOD_FORM != "disc") paste0("_", MOD_FORM) else "", "_", h_suffix, ".rds"))
  if (!file.exists(f)) stop("fit bundle missing: ", f)
  b <- readRDS(f); jm <- b$jm; ld <- b$long_data
  draws <- beta_draws(jm, b$lme)
  keep <- sample.int(nrow(draws), min(N_DRAWS, nrow(draws)))
  draws <- draws[keep, , drop = FALSE]
  binary <- isTRUE(b$binary)
  sd_log_y <- if (binary) 1 else sd(ld$log_y)   # binary outcome: report on the log-odds scale
  gate <- u$status == "converged"
  message(sprintf("  %-40s %s", tag, if (gate) "" else "(R-hat gate failed; reported for plumbing only)"))

  # ---- Q1 coefficients, per unit and per SD of the log marker
  for (term in c("l_vtpbw_within", "l_vtpbw_within:ldisc_c", "l_vtpbw_within:age10_c", "vtpbw_pt_mean",
                 "l_vtpbw_within:log_pbw", "l_vtpbw_within:log_pfvc",
                 "mean_prior_vtpfvc", "cum_days_above", "ers_pfvc_0:l_vtpbw_within")) {
    if (!term %in% colnames(draws)) next
    v <- draws[, term]
    strain_rows[[length(strain_rows) + 1L]] <- tibble(
      marker = u$marker, model = u$model, adjustment = u$adjustment, term = term,
      estimate = mean(v), lo = quantile(v, 0.025), hi = quantile(v, 0.975),
      per_sd_estimate = mean(v) / sd_log_y, per_sd_lo = quantile(v, 0.025) / sd_log_y,
      per_sd_hi = quantile(v, 0.975) / sd_log_y, sd_log_marker = sd_log_y,
      n_patients = u$n_patients, rhat_gate = gate)
  }

  # ---- PFVC-level question: marker difference per SD of log PFVC (or of log
  #      PBW/PFVC) at each horizon hour within the grid, level + divergence x time,
  #      from the joint posterior (death before H handled by the shared random effects)
  for (ex in c("log_pfvc_sd", "ldisc_sd")) {
    if (!ex %in% colnames(draws)) next
    b_lev <- draws[, ex]
    tcol <- intersect(c(paste0(ex, ":vent_day"), paste0("vent_day:", ex)), colnames(draws))   # R orders the pair by appearance
    b_tim <- if (length(tcol)) draws[, tcol[1]] else 0
    for (hh in LEVEL_HOURS[LEVEL_HOURS <= JM_HORIZON * 24]) {
      v <- b_lev + b_tim * hh / 24
      level_rows[[length(level_rows) + 1L]] <- tibble(
        marker = u$marker, model = u$model, adjustment = u$adjustment, exposure = ex, horizon_h = hh,
        estimate = mean(v), lo = quantile(v, 0.025), hi = quantile(v, 0.975),
        p_gt0 = mean(v > 0), per_sd_marker = if (binary) NA_real_ else mean(v) / sd_log_y,
        scale = if (binary) "log-odds of any pressor" else "log marker",
        n_patients = u$n_patients, n_deaths = u$n_deaths, rhat_gate = gate)
    }
  }

  # ---- Q2 association: HR per SD of the current log marker (value) and per unit slope
  al <- do.call(rbind, jm$mcmc$alphas)[keep, , drop = FALSE]
  for (cn in colnames(al)) {
    kind  <- if (grepl("value", cn)) "value" else "slope"
    cause <- if (grepl("death", cn)) "death" else "extubation"
    scale <- if (kind == "value") sd_log_y else 1
    v <- al[, cn] * scale
    assoc_rows[[length(assoc_rows) + 1L]] <- tibble(
      marker = u$marker, model = u$model, adjustment = u$adjustment, kind = kind, cause = cause,
      log_hr = mean(v), log_hr_lo = quantile(v, 0.025), log_hr_hi = quantile(v, 0.975),
      hr = exp(mean(v)), hr_lo = exp(quantile(v, 0.025)), hr_hi = exp(quantile(v, 0.975)),
      per = if (binary) "1 logit unit of P(any pressor)" else if (kind == "value") "1 SD of log marker" else "1 log-unit per day",
      n_patients = u$n_patients, n_deaths = u$n_deaths, rhat_gate = gate)
  }

  # ---- heterogeneity: strain slope at Ers x PFVC percentiles
  if (u$model == "hetero" && "ers_pfvc_0:l_vtpbw_within" %in% colnames(draws)) {
    pt <- ld %>% distinct(hospitalization_id, ers_pfvc_0)
    q <- quantile(pt$ers_pfvc_0, c(0.1, 0.5, 0.9))
    for (k in seq_along(q)) {
      v <- draws[, "l_vtpbw_within"] + draws[, "ers_pfvc_0:l_vtpbw_within"] * q[[k]]
      hetero_rows[[length(hetero_rows) + 1L]] <- tibble(
        marker = u$marker, ers_pfvc_pct = c(10, 50, 90)[k], ers_pfvc_value = q[[k]],
        strain_slope = mean(v), lo = quantile(v, 0.025), hi = quantile(v, 0.975),
        n_patients = u$n_patients, rhat_gate = gate)
    }
  }

  # ---- Q1 as a picture: dose-response trajectories by discordance. A patient at
  #      the cohort-median dose level and median age whose previous-day VT/PBW is
  #      held at 6, 8 or 10 mL/kg (the within term), at the 10th, 50th and 90th
  #      percentile of PBW/PFVC discordance (the modifier); every draw is
  #      re-centred at day 1 of the 6 mL/kg curve within its discordance facet, so
  #      the spread between curves in a facet is the dose slope there and the
  #      change in that spread across facets is the interaction.
  if (u$model == "main" && all(c("l_vtpbw_within", "ldisc_c") %in% names(ld))) {
    pt <- ld %>% distinct(hospitalization_id, vtpbw_pt_mean, ldisc_c)
    dose_med <- median(pt$vtpbw_pt_mean)
    disc_q   <- quantile(pt$ldisc_c, DISC_PCT / 100)
    t_grid <- seq(STEP, JM_HORIZON, by = min(STEP, 0.25))
    grid <- expand_grid(vent_day = t_grid, dose = DOSE_LEVELS, disc_pct = DISC_PCT) %>%
      mutate(vtpbw_pt_mean = dose_med, l_vtpbw_within = dose - dose_med,
             ldisc_c = disc_q[match(disc_pct, DISC_PCT)],
             log_pbw = median(ld$log_pbw), log_pfvc = median(ld$log_pfvc),
             log_pfvc_sd = 0, ldisc_sd = 0,
             mean_prior_vtpfvc = median(ld$l_vtpfvc, na.rm = TRUE), cum_days_above = 0)
    grid <- bind_cols(grid, population_row(ld)[rep(1L, nrow(grid)), ])
    tt <- delete.response(b$mf_terms)     # predvars carry the fitted ns() knots
    X <- model.matrix(tt, model.frame(tt, grid))[, colnames(draws), drop = FALSE]
    Y <- X %*% t(draws)                                    # rows = grid, cols = draws
    Yc <- Y
    for (dp in DISC_PCT) {
      ref1 <- which(grid$vent_day == t_grid[1] & grid$dose == DOSE_REF & grid$disc_pct == dp)
      rows <- which(grid$disc_pct == dp)
      Yc[rows, ] <- sweep(Y[rows, , drop = FALSE], 2, Y[ref1, ])
    }
    tg <- grid %>% select(vent_day, dose, disc_pct) %>%
      mutate(ldisc_c = disc_q[match(disc_pct, DISC_PCT)],
             mean = rowMeans(Yc), lo = apply(Yc, 1, quantile, 0.025), hi = apply(Yc, 1, quantile, 0.975),
             marker = u$marker, adjustment = u$adjustment, baseline_form = BASELINE_FORM,
             n_patients = u$n_patients, rhat_gate = gate, site = site_name)
    trajectory_rows[[length(trajectory_rows) + 1L]] <- tg
    if (u$adjustment == "adjusted") {
      traj_plots[[u$marker]] <- ggplot(tg %>% mutate(facet = factor(paste0("PBW/PFVC p", disc_pct), paste0("PBW/PFVC p", DISC_PCT))),
                                       aes(vent_day, mean, colour = factor(dose), fill = factor(dose))) +
        geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
        geom_line(linewidth = 1) +
        facet_wrap(~ facet, nrow = 1) +
        scale_colour_manual(values = okabe[seq_along(DOSE_LEVELS)], name = "VT/PBW (mL/kg)") +
        scale_fill_manual(values = okabe[seq_along(DOSE_LEVELS)], name = "VT/PBW (mL/kg)") +
        labs(title = b$marker$label, x = "Ventilator day",
             y = sprintf("Log marker, relative to the first period at %g mL/kg", DOSE_REF),
             subtitle = sprintf("n = %d patients%s", u$n_patients,
                                if (gate) "" else " (R-hat gate failed)")) +
        theme_minimal(base_size = 11)
    }
  }
}

trajectory_grid <- bind_rows(trajectory_rows)
strain_effects  <- bind_rows(strain_rows)  %>% mutate(grid = JM_GRID, horizon_days = JM_HORIZON, baseline_form = BASELINE_FORM, site = site_name)
association_hr  <- bind_rows(assoc_rows)   %>% mutate(grid = JM_GRID, horizon_days = JM_HORIZON, baseline_form = BASELINE_FORM, site = site_name)
heterogeneity   <- bind_rows(hetero_rows)  %>% mutate(grid = JM_GRID, horizon_days = JM_HORIZON, baseline_form = BASELINE_FORM, site = site_name)
level_contrast  <- bind_rows(level_rows)   %>% mutate(grid = JM_GRID, horizon_days = JM_HORIZON, baseline_form = BASELINE_FORM, site = site_name)
if (nrow(level_contrast)) {
  write_csv(level_contrast, file.path(final_dir, paste0("jm_level_contrast_", out_tag, ".csv")))
  message("--- marker difference per SD of the size exposure at each horizon (log units)")
  print(as.data.frame(level_contrast %>% select(marker, adjustment, exposure, horizon_h, estimate, lo, hi, p_gt0, n_patients, n_deaths) %>%
                        mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
}
write_csv(trajectory_grid, file.path(final_dir, paste0("jm_trajectory_grid_", out_tag, ".csv")))
write_csv(strain_effects,  file.path(final_dir, paste0("jm_strain_effects_",  out_tag, ".csv")))
write_csv(association_hr,  file.path(final_dir, paste0("jm_association_hr_",  out_tag, ".csv")))
if (nrow(heterogeneity)) write_csv(heterogeneity, file.path(final_dir, paste0("jm_heterogeneity_", out_tag, ".csv")))

# ---- figures
if (length(traj_plots)) {
  p <- wrap_plots(traj_plots, ncol = 1, guides = "collect") +
    plot_annotation(title = sprintf("Dose-response trajectories by PBW/PFVC discordance (first %s, %s)",
                                    if (JM_GRID == "6h") paste(as.integer(JM_HORIZON * 24), "hours") else paste(JM_HORIZON, "days"), site_name))
  ggsave(file.path(final_dir, paste0("jm_trajectories_", out_tag, ".pdf")), p,
         width = 11, height = 3.5 * length(traj_plots))
}
if (nrow(association_hr)) {
  fa <- association_hr %>% filter(model == "main") %>%
    mutate(lab = paste(marker, kind), adjustment = factor(adjustment, c("adjusted", "unadjusted")))
  p1 <- ggplot(fa, aes(hr, lab, colour = adjustment)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = hr_lo, xmax = hr_hi), position = position_dodge(width = 0.5)) +
    facet_wrap(~ cause, scales = "free_x") + scale_x_log10() +
    scale_colour_manual(values = okabe[c(5, 3)]) +
    labs(title = "Q2: cause-specific hazard per SD of the current log marker (value) or per unit slope",
         x = "Hazard ratio", y = NULL) + theme_minimal(base_size = 11)
  fs <- strain_effects %>% filter(model == "main", term %in% c("l_vtpbw_within", "l_vtpbw_within:ldisc_c")) %>%
    mutate(marker = paste(marker, if_else(term == "l_vtpbw_within", "dose slope", "x log discordance"))) %>%
    mutate(adjustment = factor(adjustment, c("adjusted", "unadjusted")))
  p2 <- ggplot(fs, aes(per_sd_estimate, marker, colour = adjustment)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = per_sd_lo, xmax = per_sd_hi), position = position_dodge(width = 0.5)) +
    scale_colour_manual(values = okabe[c(5, 3)]) +
    labs(title = "Q1: log marker (SD units) per 1 mL/kg PBW above the patient's own mean, and its modification by log PBW/PFVC",
         x = "SD of log marker per mL/kg PBW (slope) or per mL/kg per log-unit discordance (interaction)", y = NULL) +
    theme_minimal(base_size = 11)
  ggsave(file.path(final_dir, paste0("jm_forest_", out_tag, ".pdf")), p2 / p1, width = 10, height = 9)
}
message("13_biotrauma_report complete: ", nrow(trajectory_grid), " grid rows, ",
        nrow(strain_effects), " strain terms, ", nrow(association_hr), " association terms -> ", final_dir)
