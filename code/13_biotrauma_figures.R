# =============================================================================
# Script 13 (figures): the PFVC-level question, drawn from the aggregate tables
# =============================================================================
# Reads only the site's final/ CSVs (no patient rows), so it runs on any site's
# outputs and, with PBWPFVC_FIG_DIR, on a pooled folder. Five figures:
#
#   biotrauma_fig_level_contrast_{tag}.pdf   THE question: the joint model's marker
#       difference at 24/48/72 h per SD LOWER log PFVC at a given VT/PBW, oriented
#       so that right of zero is more injury for every marker (log units;
#       odds ratio for any vasopressor), adjusted beside unadjusted, with P(harm)
#   biotrauma_fig_estimators_{tag}.pdf   the marker difference per SD of log PFVC
#       at each horizon from the three estimators side by side: the joint model
#       (death before H modelled), the quick LME (longitudinal submodel alone),
#       and the fixed-horizon comparator (survivors only), adjusted and unadjusted
#   biotrauma_fig_divergence_{tag}.pdf   the rate term, log marker per day per
#       SD of log PFVC, adjusted beside unadjusted: the part of the effect that
#       the identification problem cannot reach
#   biotrauma_fig_trajectory_{tag}.pdf   the joint model's predicted difference
#       from the median-PFVC patient over the window at PFVC -1 L, median and +1 L
#       (level + divergence x time), with the 95% band
#   biotrauma_fig_pressor_{tag}.pdf      the any-pressor part as odds ratios per
#       SD LOWER PFVC at each horizon, joint model beside the comparator
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
SIZE_EX   <- switch(MOD_FORM, disc_level = "ldisc_sd", vtpfvc = "vtpfvc_c", "log_pfvc_sd")   # the form's size exposure column
SIZE_LAB  <- switch(MOD_FORM, disc_level = "per SD of log PBW/PFVC (VT/PFVC at a given VT/PBW)",
                    vtpfvc = "per point of VT/PFVC (% of predicted FVC) at a given VT/PBW", "per SD of log PFVC")
TRAJ_STEP <- if (MOD_FORM == "vtpfvc") 2 else 1   # the trajectory figure's contrast: +/- 2 points of VT/PFVC, else +/- 1 SD
FLIP_INJ  <- MOD_FORM %in% c("disc_level", "vtpfvc")   # a HIGHER value of these is the smaller lung
fig_dir   <- Sys.getenv("PBWPFVC_FIG_DIR", here("output", paste0(site_name, "_output"), "final"))
tag       <- paste0(if (MOD_FORM != "disc") paste0(MOD_FORM, "_") else "", h_suffix, "_", site_name)
okabe <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9")
theme_set(theme_minimal(base_size = 11))
worse <- c(creatinine = "higher", platelets = "lower", bilirubin = "higher", sf = "lower", dp = "higher",
           ne_equiv_peak = "higher", ne_equiv = "higher", any_pressor = "higher",
           osi = "higher", oi = "higher")
lab <- c(creatinine = "Creatinine", platelets = "Platelets", bilirubin = "Bilirubin", sf = "SF ratio",
         dp = "Driving pressure", ne_equiv_peak = "NE-equivalent dose (per kg, flagged)", ne_equiv = "NE-equivalents",
         any_pressor = "Any vasopressor (log-odds)",
         osi = "Oxygen saturation index\n(numerator-driven, flagged)", oi = "Oxygenation index\n(numerator-driven, flagged)")
read_if <- function(f) if (file.exists(f)) read_csv(f, show_col_types = FALSE) else NULL
marker_label <- function(m) paste0(lab[m], "\n(worse = ", worse[m], ")")

# ---- inputs
lc <- read_if(file.path(fig_dir, paste0("jm_level_contrast_", tag, ".csv")))
es <- read_if(file.path(fig_dir, paste0("jm_estimates_", tag, ".csv")))
if (is.null(lc) || is.null(es)) stop("no joint-model tables for tag ", tag, " in ", fig_dir)
ih <- list.files(fig_dir, "^injury_at_horizon_[a-z_]+_.*\\.csv$", full.names = TRUE) %>%
  discard(~ grepl("counts_", .x)) %>% map_dfr(read_if)
ql <- list.files(fig_dir, "^quick_lme_[a-z_]+_.*\\.csv$", full.names = TRUE) %>% map_dfr(read_if) %>%
  { if (nrow(.) && !"model_horizon_h" %in% names(.)) mutate(., model_horizon_h = horizon_h) else . }

# ---- channels form: the size effect identified through each GLI input, from the joint model
if (MOD_FORM == "channels") {
  chd <- lc %>% filter(exposure %in% CHANNELS, model == "main", marker %in% names(lab)) %>%
    mutate(channel = factor(sub("^ch_", "", exposure), c("height", "age", "sex", "race")),
           inj = if_else(worse[marker] == "higher", -estimate, estimate),
           inj_lo = if_else(worse[marker] == "higher", -hi, lo), inj_hi = if_else(worse[marker] == "higher", -lo, hi),
           horizon = factor(paste(horizon_h, "h"), paste(sort(unique(horizon_h)), "h")),
           marker_lab = sprintf("%s\n(worse = %s; n = %d, deaths = %d)", lab[marker], worse[marker], n_patients, n_deaths),
           panel = sprintf("%s, p(equal) = %.2g", horizon, p_equal))
  p_ch <- ggplot(chd, aes(inj, channel)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    geom_pointrange(aes(xmin = inj_lo, xmax = inj_hi), colour = okabe[1]) +
    facet_grid(marker_lab ~ panel, scales = "free_x") +
    labs(title = "The size effect identified through each input to PFVC, joint model (death before H modelled)",
         subtitle = paste0(site_name, ": difference toward injury per log unit LOWER of each GLI piece; ",
                           "equal pieces = lung size is the operative quantity; log-odds for any vasopressor"),
         x = "difference in the log marker toward injury per log unit lower (95% interval)", y = NULL) +
    theme(strip.text.y = element_text(angle = 0))
  ggsave(file.path(fig_dir, paste0("biotrauma_fig_channels_", tag, ".pdf")), p_ch,
         width = 10, height = 2 + 1.6 * n_distinct(chd$marker))
  message("13_biotrauma_figures (channels form): ", n_distinct(chd$marker), " markers -> ", fig_dir)
  quit(save = "no")
}

# ---- 0. the level contrasts: the scientific question, one figure
#         Joint-model marker difference at 24/48/72 h per SD LOWER log PFVC at a
#         given VT/PBW, oriented so that right of zero is MORE injury for every
#         marker (sign flipped for markers whose worse direction is lower), in log
#         units; any vasopressor as an odds ratio. Each point is
#         labelled with the posterior probability of harm.
lc0 <- lc %>% filter(exposure == SIZE_EX, model == "main", marker %in% names(lab)) %>%
  mutate(binary = marker == "any_pressor",
         # log-scale contrast per unit LOWER PFVC, in the injury direction
         inj = if_else(worse[marker] == "higher", -estimate, estimate),
         inj_lo = if_else(worse[marker] == "higher", -hi, lo),
         inj_hi = if_else(worse[marker] == "higher", -lo, hi),
         p_harm = if_else(worse[marker] == "higher", 1 - p_gt0, p_gt0),
         adjustment = factor(adjustment, c("adjusted", "unadjusted")),
         horizon = factor(paste(horizon_h, "h"), paste(sort(unique(horizon_h)), "h")),
         marker_lab = paste0(lab[marker], "\n(worse = ", worse[marker], "; n = ", n_patients, ", deaths = ", n_deaths, ")"),
         p_lab = sprintf("P(harm) %.2f", p_harm))
unit_lower <- switch(MOD_FORM, disc_level = "per SD HIGHER log PBW/PFVC (a lung smaller than PBW predicts)",
                     vtpfvc = "per point HIGHER VT/PFVC (% of predicted FVC) at the same VT/PBW", "per SD lower log PFVC")
# for the discordance form a HIGHER discordance is the smaller lung: the injury direction flips
if (FLIP_INJ) lc0 <- lc0 %>% mutate(inj = -inj, i_lo = -inj_hi, inj_hi = -inj_lo, inj_lo = i_lo, p_harm = 1 - p_harm,
                                                    p_lab = sprintf("P(harm) %.2f", p_harm)) %>% select(-i_lo)
p0_cont <- lc0 %>% filter(!binary)
p0_bin  <- lc0 %>% filter(binary)
plot_level <- function(d, xvar, lovar, hivar, ref, xlab, title, log_x = FALSE) {
  p <- ggplot(d, aes(.data[[xvar]], horizon, colour = adjustment)) +
    geom_vline(xintercept = ref, linetype = 2, colour = "grey55") +
    geom_pointrange(aes(xmin = .data[[lovar]], xmax = .data[[hivar]]), position = position_dodge(width = 0.6)) +
    geom_text(aes(label = p_lab, x = .data[[hivar]]), position = position_dodge(width = 0.6), hjust = -0.15, size = 2.8,
              show.legend = FALSE) +
    facet_grid(marker_lab ~ ., scales = "free_x") +
    scale_colour_manual(values = okabe[c(1, 2)], name = NULL) +
    (if (log_x) scale_x_log10(expand = expansion(mult = c(0.05, 0.35))) else
                scale_x_continuous(expand = expansion(mult = c(0.05, 0.35)))) +
    labs(title = title, x = xlab, y = "Horizon") +
    theme(legend.position = "top", strip.text.y = element_text(angle = 0))
  p
}
p0_list <- list()
if (nrow(p0_cont)) p0_list$cont <- plot_level(
  p0_cont, "inj", "inj_lo", "inj_hi", 0,
  xlab = paste0("difference in the log marker toward injury, ", unit_lower, " (95% interval)"),
  title = "Injury markers")
if (nrow(p0_bin)) p0_list$bin <- plot_level(
  p0_bin %>% mutate(or = exp(inj), or_lo = exp(inj_lo), or_hi = exp(inj_hi)), "or", "or_lo", "or_hi", 1,
  xlab = paste0("odds ratio of a vasopressor running, ", unit_lower, " (log scale)"),
  title = "Any vasopressor", log_x = TRUE)
if (!length(p0_list)) stop("no level contrasts for exposure ", SIZE_EX, " in tag ", tag)
p0 <- wrap_plots(p0_list, ncol = 1, heights = c(if (nrow(p0_cont)) n_distinct(p0_cont$marker), if (nrow(p0_bin)) 1)) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("Is a lower PFVC, at a given VT/PBW, associated with more injury? Joint model, ", unit_lower),
    subtitle = paste0(site_name, ": right of the line = more injury; death before the horizon is modelled\n",
                      "P(harm) = posterior probability that the contrast lies in the injury direction")) &
  theme(legend.position = "top")
ggsave(file.path(fig_dir, paste0("biotrauma_fig_level_contrast_", tag, ".pdf")), p0,
       # one row per marker x horizon, so a 7-day run needs the height a 48-hour one did not
       width = 10, height = 2.5 + 0.42 * nrow(distinct(lc0, marker, horizon_h)))

# ---- 0b. the trend across the window (three or more horizons only)
#      The level contrast is a level plus a rate times time, so on a long window
#      it is a curve rather than a set of points, and the curve is the result:
#      across three markers at MIMIC the contrast starts on the wrong side of
#      zero, crosses, and accumulates. The two panels separate the two halves.
#      Left, the contrast at each horizon, with the exact intervals from the
#      posterior (not the approximate band the trajectory figure uses). Right,
#      the rate term alone, adjusted beside unadjusted: a rate that does not move
#      when age, sex and race enter the model is not the age channel, and that
#      invariance is the argument, so it is drawn rather than described.
if (n_distinct(lc0$horizon_h) >= 3) {
  RHAT_GATE <- 1.1
  inj_sign <- function(m) if_else(worse[m] == "higher", -1, 1) * if_else(FLIP_INJ, -1, 1)
  # A marker whose exposure terms did not converge must not be drawn like one
  # that did: the any-vasopressor rate at MIMIC came back with an R-hat of 2.5
  # beside two labs at 1.01, and on the page they looked the same. Non-converged
  # estimates are drawn hollow with a dashed line and their R-hat in the strip.
  rate <- es %>%
    filter(block == "longitudinal", model == "main", marker %in% names(lab),
           term %in% c(paste0(SIZE_EX, ":vent_day"), paste0("vent_day:", SIZE_EX))) %>%
    transmute(marker, adjustment = factor(adjustment, c("adjusted", "unadjusted")),
              s = inj_sign(marker), e = s * estimate, l = pmin(s * lo, s * hi), h = pmax(s * lo, s * hi),
              rhat, ok = is.finite(rhat) & rhat <= RHAT_GATE)
  worst <- rate %>% group_by(marker) %>% summarise(mx = max(rhat, na.rm = TRUE), .groups = "drop")
  strip <- function(m) {
    w <- worst$mx[match(m, worst$marker)]
    paste0(marker_label(m), if_else(is.finite(w) & w > RHAT_GATE, sprintf("\nR-hat %.2f: DID NOT CONVERGE", w), ""))
  }
  lvl <- sort(unique(strip(unique(rate$marker))))   # one order for both panels
  rate  <- rate  %>% mutate(marker_lab = factor(strip(marker), lvl))
  trend <- lc0 %>% mutate(day = horizon_h / 24,
                          ok = marker %in% rate$marker[rate$ok],
                          marker_lab = factor(strip(marker), lvl))
  pt_a <- ggplot(trend, aes(day, inj, colour = adjustment, fill = adjustment)) +
    geom_hline(yintercept = 0, colour = "grey55") +
    geom_ribbon(aes(ymin = inj_lo, ymax = inj_hi), alpha = 0.12, colour = NA) +
    geom_line(aes(linetype = ok), linewidth = 1) +
    geom_point(aes(shape = ok), size = 1.8) +
    facet_wrap(~ marker_lab, scales = "free_y", ncol = 1) +
    scale_colour_manual(values = okabe[c(1, 2)], name = NULL) +
    scale_fill_manual(values = okabe[c(1, 2)], guide = "none") +
    scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"), guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
    labs(title = "The contrast as it accumulates", subtitle = "above zero is more injury",
         x = "days from the index", y = "log units (log-odds for any vasopressor)")
  pt_b <- ggplot(rate, aes(e, adjustment, colour = adjustment)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    geom_pointrange(aes(xmin = l, xmax = h, shape = ok)) +
    facet_wrap(~ marker_lab, scales = "free_x", ncol = 1) +
    scale_colour_manual(values = okabe[c(1, 2)], guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
    labs(title = "The rate alone, per day",
         subtitle = "a rate that adjustment does not move is not the age channel",
         x = paste0("change per day toward injury, ", unit_lower), y = NULL)
  pt <- pt_a + pt_b + plot_layout(widths = c(1.15, 1), guides = "collect") +
    plot_annotation(title = paste0("Marker trends over ", JM_HORIZON, " days (joint model, ", unit_lower, ")"),
                    subtitle = paste0(site_name, ": hollow points and dashed lines did not converge")) &
    theme(legend.position = "top")
  ggsave(file.path(fig_dir, paste0("biotrauma_fig_trend_", tag, ".pdf")), pt,
         width = 13, height = 2.5 + 2.1 * n_distinct(trend$marker))
}

# ---- 1. estimator comparison
est <- bind_rows(
  lc %>% filter(exposure == SIZE_EX, model == "main") %>%
    transmute(marker, adjustment, horizon_h, estimate, lo, hi, estimator = "Joint model (death modelled)"),
  if (nrow(ql)) ql %>% filter(exposure == SIZE_EX) %>%
    transmute(marker, adjustment, horizon_h = model_horizon_h, estimate, lo, hi, estimator = "Longitudinal model only"),
  if (nrow(ih)) ih %>% filter(exposure == SIZE_EX, outcome_type %in% c("log marker at H", "any pressor at H")) %>%
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
  labs(title = paste0("Marker difference ", SIZE_LAB),
       subtitle = paste0(site_name, if (FLIP_INJ) ": a larger lung for the dose is the negative of each estimate" else ": a LOWER PFVC is the negative of each estimate",
                         "; log-odds for any vasopressor"),
       x = paste0("log units ", SIZE_LAB, " (95% interval)"), y = NULL) +
  theme(legend.position = "top", strip.text.y = element_text(angle = 0))
ggsave(file.path(fig_dir, paste0("biotrauma_fig_estimators_", tag, ".pdf")), p1,
       width = 10, height = 2 + 1.6 * n_distinct(est$marker))

# ---- 2. divergence per day (the rate term), adjusted beside unadjusted
div <- es %>% filter(block == "longitudinal", term %in% c(paste0("vent_day:", SIZE_EX), paste0(SIZE_EX, ":vent_day")),
                     model == "main", marker %in% names(lab)) %>%
  transmute(marker, adjustment = factor(adjustment, c("adjusted", "unadjusted")), estimate, lo, hi, rhat,
            marker_lab = marker_label(marker))
lev <- es %>% filter(block == "longitudinal", term == SIZE_EX, model == "main", marker %in% names(lab)) %>%
  transmute(marker, adjustment = factor(adjustment, c("adjusted", "unadjusted")), estimate, lo, hi, rhat,
            marker_lab = marker_label(marker))
p2 <- (ggplot(lev, aes(estimate, marker_lab, colour = adjustment)) +
         geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
         geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5)) +
         scale_colour_manual(values = okabe[c(1, 2)], name = NULL) +
         labs(title = "Level: difference at the start", x = paste0("log units ", SIZE_LAB), y = NULL)) /
      (ggplot(div, aes(estimate, marker_lab, colour = adjustment)) +
         geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
         geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5)) +
         scale_colour_manual(values = okabe[c(1, 2)], name = NULL) +
         labs(title = "Divergence: change per day", x = paste0("log units per day ", SIZE_LAB), y = NULL)) +
  plot_layout(guides = "collect") +
  plot_annotation(title = paste0("The ", switch(MOD_FORM, vtpfvc = "VT/PFVC", disc_level = "PBW/PFVC", "PFVC"), " effect split into a level and a rate (joint model)"),
                  subtitle = paste0(site_name, ": a rate that is the same adjusted and unadjusted is not the age channel")) &
  theme(legend.position = "top")
ggsave(file.path(fig_dir, paste0("biotrauma_fig_divergence_", tag, ".pdf")), p2, width = 8, height = 2.5 + 1.1 * n_distinct(div$marker))

# ---- 3. predicted trajectory difference at PFVC -1 / 0 / +1 SD over the window
TRAJ_LABELS <- if (MOD_FORM == "vtpfvc") c("-2 points (larger lung for the dose)", "median", "+2 points (smaller lung for the dose)") else
               if (FLIP_INJ) c("-1 SD (larger lung for the dose)", "median", "+1 SD (smaller lung for the dose)") else
                              c("-1 SD (smaller lung)", "median", "+1 SD (larger lung)")
#         (level + divergence x t; the band is the interval of the sum at each t,
#          taken from the posterior intervals of the two terms assuming independence,
#          so it is approximate; the contrast table carries the exact intervals at 24/48/72 h)
traj <- lev %>% transmute(marker, adjustment, level = estimate, level_se = (hi - lo) / 3.92) %>%
  inner_join(div %>% transmute(marker, adjustment, slope = estimate, slope_se = (hi - lo) / 3.92),
             by = c("marker", "adjustment")) %>%
  crossing(t_day = seq(0, JM_HORIZON, by = STEP), pfvc_sd = c(-1, 0, 1)) %>%
  mutate(diff = pfvc_sd * TRAJ_STEP * (level + slope * t_day),
         se = abs(pfvc_sd) * TRAJ_STEP * sqrt(level_se^2 + (slope_se * t_day)^2),
         lo = diff - 1.96 * se, hi = diff + 1.96 * se,
         pfvc_lab = factor(TRAJ_LABELS[pfvc_sd + 2], TRAJ_LABELS),
         marker_lab = marker_label(marker))
p3 <- ggplot(traj %>% filter(pfvc_sd != 0), aes(t_day * 24, diff, colour = pfvc_lab, fill = pfvc_lab)) +
  geom_hline(yintercept = 0, colour = "grey55") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  facet_grid(marker_lab ~ adjustment, scales = "free_y") +
  scale_colour_manual(values = if (FLIP_INJ) okabe[c(3, 4)] else okabe[c(4, 3)], name = if (MOD_FORM == "vtpfvc") "VT/PFVC" else if (MOD_FORM == "disc_level") "PBW/PFVC" else "PFVC") +
  scale_fill_manual(values = if (FLIP_INJ) okabe[c(3, 4)] else okabe[c(4, 3)], name = if (MOD_FORM == "vtpfvc") "VT/PFVC" else if (MOD_FORM == "disc_level") "PBW/PFVC" else "PFVC") +
  labs(title = paste0("Predicted marker difference from the median patient over the window (",
                      switch(MOD_FORM, vtpfvc = "VT/PFVC at the same VT/PBW", disc_level = "PBW/PFVC", "PFVC"), ")"),
       subtitle = paste0(site_name, ": joint model, level + divergence x time; approximate band"),
       x = "Hours from the index", y = "difference in log marker (log-odds for any vasopressor)") +
  theme(legend.position = "top", strip.text.y = element_text(angle = 0))
ggsave(file.path(fig_dir, paste0("biotrauma_fig_trajectory_", tag, ".pdf")), p3, width = 9, height = 2 + 1.5 * n_distinct(traj$marker))

# ---- 4. any vasopressor as odds ratios per SD LOWER PFVC, joint model beside comparator
ap <- est %>% filter(marker == "any_pressor") %>%
  mutate(or = exp(-estimate), or_lo = exp(-hi), or_hi = exp(-lo))
if (nrow(ap)) {
  p4 <- ggplot(ap, aes(or, estimator, colour = adjustment)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
    geom_pointrange(aes(xmin = or_lo, xmax = or_hi), position = position_dodge(width = 0.6)) +
    facet_wrap(~ horizon, nrow = 1) + scale_x_log10() +
    scale_colour_manual(values = okabe[c(1, 2)], name = NULL) +
    labs(title = paste0("Odds of a vasopressor running ", unit_lower),
         subtitle = paste0(site_name, ": the joint model's odds ratio is subject-specific, the comparator's marginal"),
         x = paste0("odds ratio ", unit_lower, " (log scale)"), y = NULL) +
    theme(legend.position = "top")
  ggsave(file.path(fig_dir, paste0("biotrauma_fig_pressor_", tag, ".pdf")), p4, width = 9, height = 3.5)
}
message("13_biotrauma_figures: ", n_distinct(est$marker), " markers, ", n_distinct(est$estimator), " estimators -> ", fig_dir)
