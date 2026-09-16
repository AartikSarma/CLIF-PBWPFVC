# =============================================================================
# Script 13 (figures): the PFVC-level question, drawn from the aggregate tables
# =============================================================================
# Reads only the site's final/ CSVs (no patient rows), so it runs on any site's
# outputs and, with PBWPFVC_FIG_DIR, on a pooled folder. Four figures:
#
#   biotrauma_fig_estimators_{tag}.pdf   the marker difference per unit PFVC (100 mL by default; PBWPFVC_PFVC_EXPO)
#       at each horizon from the three estimators side by side: the joint model
#       (death before H modelled), the quick LME (longitudinal submodel alone),
#       and the fixed-horizon comparator (survivors only), adjusted and unadjusted
#   biotrauma_fig_divergence_{tag}.pdf   the rate term, log marker per day per
#       100 mL PFVC, adjusted beside unadjusted: the part of the effect that
#       the identification problem cannot reach
#   biotrauma_fig_trajectory_{tag}.pdf   the joint model's predicted difference
#       from the median-PFVC patient over the window at PFVC -1 L, median and +1 L
#       (level + divergence x time), with the 95% band
#   biotrauma_fig_pressor_{tag}.pdf      the any-pressor part as odds ratios per
#       100 mL LOWER PFVC at each horizon, joint model beside the comparator
#
# A lower PFVC is the negative of every log-marker estimate; the figures label
# the injury direction per marker so the eye does not have to flip signs.
#
# Usage: PBWPFVC_JM_MODIFIER=pfvc Rscript code/13_biotrauma_figures.R
#        (PBWPFVC_JM_GRID / _HORIZON_H select the tag as for the fit; PBWPFVC_FIG_DIR
#         points at another folder of the same CSVs, e.g. a site's or the pooled one)
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(here); library(patchwork) })
rm(list = ls())
source("utils/config.R")
site_name <- config$site_name
source(here("code", "13_biotrauma_grid.R"))
MOD_FORM  <- Sys.getenv("PBWPFVC_JM_MODIFIER", "pfvc")
fig_dir   <- Sys.getenv("PBWPFVC_FIG_DIR", here("output", paste0(site_name, "_output"), "final"))
tag       <- paste0(if (MOD_FORM != "disc") paste0(MOD_FORM, "_") else "", h_suffix, "_", site_name)
okabe <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9")
theme_set(theme_minimal(base_size = 11))
worse <- c(creatinine = "higher", platelets = "lower", bilirubin = "higher", sf = "lower", dp = "higher",
           ne_equiv_peak = "higher", ne_equiv = "higher", any_pressor = "higher")
lab <- c(creatinine = "Creatinine", platelets = "Platelets", bilirubin = "Bilirubin", sf = "SF ratio",
         dp = "Driving pressure", ne_equiv_peak = "NE-equivalent dose (per kg, flagged)", ne_equiv = "NE-equivalents",
         any_pressor = "Any vasopressor (log-odds)")
read_if <- function(f) if (file.exists(f)) read_csv(f, show_col_types = FALSE) else NULL
marker_label <- function(m) paste0(lab[m], "\n(worse = ", worse[m], ")")

# ---- inputs
lc <- read_if(file.path(fig_dir, paste0("jm_level_contrast_", tag, ".csv")))
es <- read_if(file.path(fig_dir, paste0("jm_estimates_", tag, ".csv")))
if (is.null(lc) || is.null(es)) stop("no joint-model tables for tag ", tag, " in ", fig_dir)
ih <- list.files(fig_dir, "^injury_at_horizon_[a-z_]+_" , full.names = TRUE) %>%
  discard(~ grepl("counts_", .x)) %>% map_dfr(read_if)
ql <- list.files(fig_dir, "^quick_lme_[a-z_]+_", full.names = TRUE) %>% map_dfr(read_if) %>%
  { if (nrow(.) && !"model_horizon_h" %in% names(.)) mutate(., model_horizon_h = horizon_h) else . }

# ---- 1. estimator comparison
est <- bind_rows(
  lc %>% filter(exposure == PFVC_EXPO, model == "main") %>%
    transmute(marker, adjustment, horizon_h, estimate, lo, hi, estimator = "Joint model (death modelled)"),
  if (nrow(ql)) ql %>% filter(exposure == PFVC_EXPO) %>%
    transmute(marker, adjustment, horizon_h = model_horizon_h, estimate, lo, hi, estimator = "Longitudinal model only"),
  if (nrow(ih)) ih %>% filter(exposure == PFVC_EXPO, outcome_type %in% c("log marker at H", "any pressor at H")) %>%
    mutate(marker = if_else(outcome_type == "any pressor at H", "any_pressor", marker)) %>%
    transmute(marker, adjustment, horizon_h, estimate, lo, hi, estimator = "Fixed horizon, survivors only")
) %>%
  filter(marker %in% names(lab)) %>%
  mutate(estimator = factor(estimator, c("Joint model (death modelled)", "Longitudinal model only", "Fixed horizon, survivors only")),
         adjustment = factor(adjustment, c("adjusted", "unadjusted")),
         horizon = factor(paste(horizon_h, "h"), paste(sort(unique(horizon_h)), "h")),
         marker_lab = marker_label(marker))
p1 <- ggplot(est, aes(estimate, estimator, colour = adjustment)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6)) +
  facet_grid(marker_lab ~ horizon, scales = "free_x") +
  scale_colour_manual(values = okabe[c(1, 2)], name = NULL) +
  labs(title = paste("Marker difference", PFVC_UNIT, "at a given VT/PBW"),
       subtitle = paste0(site_name, ": a LOWER PFVC is the negative of each estimate; log-odds for any vasopressor"),
       x = paste("log units", PFVC_UNIT, "(95% interval)"), y = NULL) +
  theme(legend.position = "top", strip.text.y = element_text(angle = 0))
ggsave(file.path(fig_dir, paste0("biotrauma_fig_estimators_", tag, ".pdf")), p1,
       width = 10, height = 2 + 1.6 * n_distinct(est$marker))

# ---- 2. divergence per day (the rate term), adjusted beside unadjusted
div <- es %>% filter(block == "longitudinal", term %in% c(paste0("vent_day:", PFVC_EXPO), paste0(PFVC_EXPO, ":vent_day")),
                     model == "main", marker %in% names(lab)) %>%
  transmute(marker, adjustment = factor(adjustment, c("adjusted", "unadjusted")), estimate, lo, hi, rhat,
            marker_lab = marker_label(marker))
lev <- es %>% filter(block == "longitudinal", term == PFVC_EXPO, model == "main", marker %in% names(lab)) %>%
  transmute(marker, adjustment = factor(adjustment, c("adjusted", "unadjusted")), estimate, lo, hi, rhat,
            marker_lab = marker_label(marker))
p2 <- (ggplot(lev, aes(estimate, marker_lab, colour = adjustment)) +
         geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
         geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5)) +
         scale_colour_manual(values = okabe[c(1, 2)], name = NULL) +
         labs(title = "Level: difference at the start", x = paste("log units", PFVC_UNIT), y = NULL)) /
      (ggplot(div, aes(estimate, marker_lab, colour = adjustment)) +
         geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
         geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5)) +
         scale_colour_manual(values = okabe[c(1, 2)], name = NULL) +
         labs(title = "Divergence: change per day", x = paste("log units per day", PFVC_UNIT), y = NULL)) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "The PFVC effect split into a level and a rate (joint model)",
                  subtitle = paste0(site_name, ": a rate that is the same adjusted and unadjusted is not the age channel")) &
  theme(legend.position = "top")
ggsave(file.path(fig_dir, paste0("biotrauma_fig_divergence_", tag, ".pdf")), p2, width = 8, height = 2.5 + 1.1 * n_distinct(div$marker))

# ---- 3. predicted trajectory difference at PFVC -1 / 0 / +1 step over the window (1 L, or 1 SD of log PFVC)
TRAJ_LABELS <- if (PFVC_EXPO == "pfvc_100") c("-1 L (smaller lung)", "median", "+1 L (larger lung)") else
                                            c("-1 SD (smaller lung)", "median", "+1 SD (larger lung)")
#         (level + divergence x t; the band is the interval of the sum at each t,
#          taken from the posterior intervals of the two terms assuming independence,
#          so it is approximate; the contrast table carries the exact intervals at 24/48/72 h)
traj <- lev %>% transmute(marker, adjustment, level = estimate, level_se = (hi - lo) / 3.92) %>%
  inner_join(div %>% transmute(marker, adjustment, slope = estimate, slope_se = (hi - lo) / 3.92),
             by = c("marker", "adjustment")) %>%
  crossing(t_day = seq(0, JM_HORIZON, by = STEP), pfvc_sd = c(-1, 0, 1)) %>%
  mutate(step = if (PFVC_EXPO == "pfvc_100") 10 else 1,   # +/- 1 L when the unit is 100 mL, else +/- 1 SD
         diff = pfvc_sd * step * (level + slope * t_day),
         se = abs(pfvc_sd) * step * sqrt(level_se^2 + (slope_se * t_day)^2),
         lo = diff - 1.96 * se, hi = diff + 1.96 * se,
         pfvc_lab = factor(TRAJ_LABELS[pfvc_sd + 2], TRAJ_LABELS),
         marker_lab = marker_label(marker))
p3 <- ggplot(traj %>% filter(pfvc_sd != 0), aes(t_day * 24, diff, colour = pfvc_lab, fill = pfvc_lab)) +
  geom_hline(yintercept = 0, colour = "grey55") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  facet_grid(marker_lab ~ adjustment, scales = "free_y") +
  scale_colour_manual(values = okabe[c(4, 3)], name = "PFVC") + scale_fill_manual(values = okabe[c(4, 3)], name = "PFVC") +
  labs(title = "Predicted marker difference from the median-PFVC patient over the window",
       subtitle = paste0(site_name, ": joint model, level + divergence x time; approximate band"),
       x = "Hours from the index", y = "difference in log marker (log-odds for any vasopressor)") +
  theme(legend.position = "top", strip.text.y = element_text(angle = 0))
ggsave(file.path(fig_dir, paste0("biotrauma_fig_trajectory_", tag, ".pdf")), p3, width = 9, height = 2 + 1.5 * n_distinct(traj$marker))

# ---- 4. any vasopressor as odds ratios per 100 mL LOWER PFVC, joint model beside comparator
ap <- est %>% filter(marker == "any_pressor") %>%
  mutate(or = exp(-estimate), or_lo = exp(-hi), or_hi = exp(-lo))
if (nrow(ap)) {
  p4 <- ggplot(ap, aes(or, estimator, colour = adjustment)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
    geom_pointrange(aes(xmin = or_lo, xmax = or_hi), position = position_dodge(width = 0.6)) +
    facet_wrap(~ horizon, nrow = 1) + scale_x_log10() +
    scale_colour_manual(values = okabe[c(1, 2)], name = NULL) +
    labs(title = paste0("Odds of a vasopressor running ", PFVC_UNIT, " LOWER, at a given VT/PBW"),
         subtitle = paste0(site_name, ": the joint model's odds ratio is subject-specific, the comparator's marginal"),
         x = paste0("odds ratio ", PFVC_UNIT, " lower (log scale)"), y = NULL) +
    theme(legend.position = "top")
  ggsave(file.path(fig_dir, paste0("biotrauma_fig_pressor_", tag, ".pdf")), p4, width = 9, height = 3.5)
}
message("13_biotrauma_figures: ", n_distinct(est$marker), " markers, ", n_distinct(est$estimator), " estimators -> ", fig_dir)
