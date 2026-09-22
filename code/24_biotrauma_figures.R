# =============================================================================
# Script 24 (figures): the PFVC-level question, drawn from the aggregate tables
# =============================================================================
# Reads only the site's final/ CSVs (no patient rows), so it runs on any site's
# outputs and, with PBWPFVC_FIG_DIR, on a pooled folder.
#
#   biotrauma_fig_main_{tag}.pdf   THE figure, for a window with three or more
#       horizons (the 7-day run): one row per marker, the contrast toward injury
#       from hour 0 to the end of the window, the rate per day adjusted beside
#       unadjusted, and the posterior probability of harm by day. The others are
#       its supporting detail:
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
# Usage: PBWPFVC_JM_MODIFIER=pfvc Rscript code/24_biotrauma_figures.R
#        (PBWPFVC_JM_GRID / _HORIZON_H select the tag as for the fit; PBWPFVC_FIG_DIR
#         points at another folder of the same CSVs, e.g. a site's or the pooled one)
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(here); library(patchwork) })
rm(list = ls())
source("utils/config.R")
site_name <- config$site_name
source(here("code", "20_biotrauma_grid.R"))
MOD_FORM  <- Sys.getenv("PBWPFVC_JM_MODIFIER", "pfvc")
SIZE_EX   <- switch(MOD_FORM, disc_level = "ldisc_sd", vtpfvc = "vtpfvc_c", "log_pfvc_sd")   # the form's size exposure column (pfvc_dose: at the median dose)
SIZE_LAB  <- switch(MOD_FORM, disc_level = "per SD of log PBW/PFVC (VT/PFVC at a given VT/PBW)",
                    vtpfvc = "per point of VT/PFVC (% of predicted FVC) at a given VT/PBW", "per SD of log PFVC")
TRAJ_STEP <- if (MOD_FORM == "vtpfvc") 2 else 1   # the trajectory figure's contrast: +/- 2 points of VT/PFVC, else +/- 1 SD
FLIP_INJ  <- MOD_FORM %in% c("disc_level", "vtpfvc")   # a HIGHER value of these is the smaller lung
fig_dir   <- Sys.getenv("PBWPFVC_FIG_DIR", final_dir_for("injury"))
tag       <- paste0(restrict_tag, if (MOD_FORM != "disc") paste0(MOD_FORM, "_") else "", h_suffix, "_", site_name)
okabe <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9")
theme_set(theme_minimal(base_size = 11))
worse <- c(creatinine = "higher", platelets = "lower", bilirubin = "higher", sf = "lower", dp = "higher",
           ne_equiv_peak = "higher", pressor_dose = "higher", ne_equiv = "higher", any_pressor = "higher",
           osi = "higher", oi = "higher")
lab <- c(creatinine = "Creatinine", platelets = "Platelets", bilirubin = "Bilirubin", sf = "SF ratio",
         dp = "Driving pressure", ne_equiv_peak = "Vasopressor dose\n(NE-equivalents per kg,\nzero days included)",
         pressor_dose = "Vasopressor dose\n(NE-equivalents per kg,\ndays on a pressor)", ne_equiv = "NE-equivalents",
         any_pressor = "Any vasopressor (log-odds)",
         osi = "Oxygen saturation index\n(numerator-driven, flagged)", oi = "Oxygenation index\n(numerator-driven, flagged)")
read_if <- function(f) if (file.exists(f)) read_csv(f, show_col_types = FALSE) else NULL
# Vasopressor dose is per kg, the clinician's dosing scale (heavier patients need
# more drug in absolute terms, which is why it is dosed per kg), so, like VT/PBW,
# it is read on that scale without a size caveat (user, 2026-09-21).
RRT_MARKERS <- character(0)   # markers taken from the dialysis-as-third-cause run
marker_label <- function(m) paste0(lab[m], "\n(worse = ", worse[m], ")",
                                   if_else(m %in% RRT_MARKERS, "\ndialysis modelled as a third cause", ""))

# ---- inputs
lc <- read_if(file.path(fig_dir, paste0("jm_level_contrast_", tag, ".csv")))
es <- read_if(file.path(fig_dir, paste0("jm_estimates_", tag, ".csv")))
# PBWPFVC_JM_WITH_RRT=1 folds the dialysis-as-third-cause run (creatinine) into
# these figures. That fit lives in its own `rrtcause_` tables so a sensitivity
# never overwrites the primary, which also means it is invisible on the shared
# figure unless asked for. It is the better creatinine estimate, because
# dialysis is started BECAUSE creatinine is rising and the two-cause fit drops
# those days without modelling why, so where both exist this one wins and the
# marker is labelled to say so.
if (identical(Sys.getenv("PBWPFVC_JM_WITH_RRT", "0"), "1")) {
  rtag <- paste0("rrtcause_", tag)
  rlc <- read_if(file.path(fig_dir, paste0("jm_level_contrast_", rtag, ".csv")))
  res <- read_if(file.path(fig_dir, paste0("jm_estimates_", rtag, ".csv")))
  if (is.null(rlc) || is.null(res)) {
    message("PBWPFVC_JM_WITH_RRT=1 but no ", rtag, " tables here; ignoring")
  } else {
    RRT_MARKERS <- sort(unique(rlc$marker))
    message("folding in the dialysis-as-third-cause fit for: ", paste(RRT_MARKERS, collapse = ", "))
    as_text <- function(d) d %>% mutate(across(everything(), as.character))
    retype  <- function(d) d %>% type_convert(guess_integer = TRUE, na = c("", "NA"))
    lc <- if (is.null(lc)) rlc else retype(bind_rows(as_text(lc %>% filter(!marker %in% RRT_MARKERS)), as_text(rlc)))
    es <- if (is.null(es)) res else retype(bind_rows(as_text(es %>% filter(!marker %in% RRT_MARKERS)), as_text(res)))
  }
}
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
  message("24_biotrauma_figures (channels form): ", n_distinct(chd$marker), " markers -> ", fig_dir)
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

# ---- 0b. THE figure: one row per marker, three panels (three or more horizons only)
#      A  the contrast as it accumulates: the joint model's marker difference toward
#         injury per unit lower PFVC at each horizon, hour 0 through the end of the
#         window, with its exact posterior interval. On a long window this is a line
#         (a level plus a rate times time), and the line is the result: the level at
#         hour 0 moves with adjustment, the slope does not.
#      B  the rate alone, per day, adjusted beside unadjusted: a rate that does not
#         move when age, sex and race enter the model is not the age channel.
#      C  the posterior probability that the contrast lies in the injury direction,
#         by horizon: where the evidence starts and where it ends up.
#      Rows carry the patient and death counts. A marker whose exposure terms did not
#      converge is drawn hollow and dashed, placed last, and says so in its label:
#      the any-vasopressor rate at MIMIC came back with an R-hat of 2.5 beside labs at
#      1.01, and without the flag they looked the same.
if (n_distinct(lc0$horizon_h) >= 3) {
  RHAT_GATE <- 1.1
  inj_sign <- function(m) if_else(worse[m] == "higher", -1, 1) * if_else(FLIP_INJ, -1, 1)
  rate <- es %>%
    filter(block == "longitudinal", model == "main", marker %in% names(lab),
           term %in% c(paste0(SIZE_EX, ":vent_day"), paste0("vent_day:", SIZE_EX))) %>%
    transmute(marker, adjustment = factor(adjustment, c("adjusted", "unadjusted")),
              s = inj_sign(marker), e = s * estimate, l = pmin(s * lo, s * hi), h = pmax(s * lo, s * hi),
              rhat, ok = is.finite(rhat) & rhat <= RHAT_GATE)
  failed <- rate %>% filter(!ok) %>% group_by(marker) %>%
    summarise(note = paste0("\n", paste(sprintf("%s R-hat %.2f", adjustment, rhat), collapse = "; "), ": not converged"),
              .groups = "drop")
  counts <- es %>% filter(model == "main") %>% distinct(marker, n_patients, n_deaths) %>%
    group_by(marker) %>% slice(1) %>% ungroup()
  # a count suppressed as under 10 (mask_small_counts) reads "<10", not "NA"
  show_count <- function(x) ifelse(is.na(x), "<10", formatC(x, big.mark = ",", format = "d"))
  row_label <- function(m) {
    paste0(marker_label(m),
           sprintf("\nn = %s, deaths = %s", show_count(counts$n_patients[match(m, counts$marker)]),
                   show_count(counts$n_deaths[match(m, counts$marker)])),
           coalesce(failed$note[match(m, failed$marker)], ""))
  }
  # Rows: PBWPFVC_FIG_MARKERS (comma list) picks the markers and their order, e.g.
  # "platelets,bilirubin,creatinine,ne_equiv_peak" to show the vasopressor as dose in
  # place of the unconverged yes/no model. Without it, every marker in the tables:
  # PBWPFVC_FIG_ORDER first (default: the markers where the rate carries the finding),
  # the rest in the order of `lab`. Unconverged markers go last either way.
  fig_markers <- trimws(strsplit(Sys.getenv("PBWPFVC_FIG_MARKERS", ""), ",")[[1]])
  fig_order   <- trimws(strsplit(Sys.getenv("PBWPFVC_FIG_ORDER", "platelets,bilirubin,creatinine"), ",")[[1]])
  present   <- if (length(fig_markers)) intersect(fig_markers, unique(rate$marker)) else
    intersect(unique(c(fig_order, names(lab))), unique(rate$marker))
  if (length(fig_markers) && length(setdiff(fig_markers, present)))
    message("24_biotrauma_figures: PBWPFVC_FIG_MARKERS names markers with no rate in ", tag, ": ",
            paste(setdiff(fig_markers, present), collapse = ", "))
  rate <- rate %>% filter(marker %in% present)
  converged <- setdiff(present, failed$marker)
  row_order <- row_label(c(converged, setdiff(present, converged)))
  rate  <- rate %>% mutate(marker_lab = factor(row_label(marker), row_order))
  trend <- lc0 %>% filter(marker %in% present) %>%
    left_join(rate %>% select(marker, adjustment, ok), by = c("marker", "adjustment")) %>%
    mutate(day = horizon_h / 24, ok = coalesce(ok, FALSE),
           marker_lab = factor(row_label(marker), row_order))
  adj_colours <- okabe[c(1, 2)]
  shared <- list(scale_colour_manual(values = adj_colours, name = NULL),
                 scale_fill_manual(values = adj_colours, guide = "none"),
                 scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"), guide = "none"),
                 scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none"),
                 theme(strip.text.y = element_blank(), plot.title = element_text(face = "bold")))
  day_breaks <- sort(unique(trend$day))
  pm_a <- ggplot(trend, aes(day, inj, colour = adjustment, fill = adjustment)) +
    geom_hline(yintercept = 0, colour = "grey55") +
    geom_ribbon(aes(ymin = inj_lo, ymax = inj_hi), alpha = 0.12, colour = NA) +
    geom_line(aes(linetype = ok), linewidth = 0.9) +
    geom_point(aes(shape = ok), size = 1.6) +
    facet_grid(marker_lab ~ ., scales = "free_y", drop = FALSE) +   # a row for every marker in panel B
    scale_x_continuous(breaks = day_breaks) + shared +
    labs(title = "Ventilated: difference toward injury", subtitle = "above zero = more injury with a smaller lung",
         x = "days from the index", y = "log marker (log-odds for any vasopressor)")
  # ---- panel B: the rate in every arm. The rate is where the argument lives: if small
  #      predicted lungs diverge without a ventilator, the ventilated divergence is not
  #      about the breath; if it grows with baseline hypoxaemia, it tracks the lung.
  #      Arms, left to right:
  #        No support                       the negative control (final/controls/): every
  #                                         patient, the rate read at the ventilated
  #                                         cohort's mean severity anchor (the sevstd_
  #                                         tables; PBWPFVC_JM_SEV_CENTER, 29_run_figure4.sh)
  #        Ventilated, SF <class>           the ventilated cohort by baseline SF
  #        Ventilated, all                  the ventilated cohort (panels A and C)
  #      The noninvasive cohort is deliberately NOT drawn: NIPPV delivers large, unlimited
  #      positive-pressure volumes, so it is a strained group and cannot be a control.
  #      Creatinine takes each arm's dialysis-as-third-cause fit when PBWPFVC_JM_WITH_RRT=1
  #      and that fit exists. Controls are read from final/controls/ of this site, or from
  #      PBWPFVC_FIG_CONTROLS_DIR; with another PBWPFVC_FIG_DIR and none named, no controls.
  ctrl_dir <- Sys.getenv("PBWPFVC_FIG_CONTROLS_DIR",
                         if (nzchar(Sys.getenv("PBWPFVC_FIG_DIR", ""))) "" else file.path(config$final_root, "controls"))
  rate_rows <- function(est) est %>%
    filter(block == "longitudinal", model == "main", marker %in% present,
           term %in% c(paste0(SIZE_EX, ":vent_day"), paste0("vent_day:", SIZE_EX))) %>%
    transmute(marker, adjustment = factor(adjustment, c("adjusted", "unadjusted")),
              s = inj_sign(marker), e = s * estimate, l = pmin(s * lo, s * hi), h = pmax(s * lo, s * hi),
              rhat, ok = is.finite(rhat) & rhat <= RHAT_GATE)
  # one arm's estimates; its dialysis-as-third-cause twin replaces the RRT markers' rows
  read_arm <- function(folder, restriction, site_tag) {
    stub <- paste0(MOD_FORM, "_", h_suffix, "_", site_tag, ".csv")
    est <- read_if(file.path(folder, paste0("jm_estimates_", restriction, stub)))
    if (is.null(est)) return(NULL)
    if (length(RRT_MARKERS)) {
      twin <- read_if(file.path(folder, paste0("jm_estimates_rrtcause_", restriction, stub)))
      if (!is.null(twin)) est <- bind_rows(est %>% filter(!marker %in% RRT_MARKERS) %>% mutate(across(everything(), as.character)),
                                           twin %>% filter(marker %in% RRT_MARKERS) %>% mutate(across(everything(), as.character))) %>%
        type_convert(guess_integer = TRUE, na = c("", "NA"))
    }
    est
  }
  restrictions_in <- function(folder, site_tag) {
    stub <- paste0(MOD_FORM, "_", h_suffix, "_", site_tag, ".csv")
    found <- list.files(folder, pattern = paste0("^jm_estimates_.*", stub, "$"))
    setdiff(sub(paste0(stub, "$"), "", sub("^jm_estimates_", "", found)), NA)
  }
  arms <- list()
  if (nzchar(ctrl_dir) && dir.exists(ctrl_dir)) {
    ctrl_site <- paste0(site_name, "_nosupport")
    # only the severity-standardised control: floors and the unstandardised control are
    # retired designs, and their tables on disk must not reach the figure
    # its rate is per SD of log PFVC in the CONTROL's panel; on the figure it is put on
    # the ventilated cohort's unit (the SDs from 22_biotrauma_fit.R's jm_scale_* tables),
    # and without both SDs the control is not drawn
    scale_v <- read_if(file.path(fig_dir, paste0("jm_scale_", h_suffix, "_", site_name, ".csv")))
    scale_c <- read_if(file.path(ctrl_dir, paste0("jm_scale_", h_suffix, "_", ctrl_site, ".csv")))
    if ("sevstd_" %in% restrictions_in(ctrl_dir, ctrl_site)) {
      if (is.null(scale_v) || is.null(scale_c) || SIZE_EX != "log_pfvc_sd")
        message("24_biotrauma_figures: control arm not drawn: its unit cannot be put on the ventilated cohort's ",
                "(jm_scale_* missing, or the form is not pfvc)")
      else arms[["control"]] <- list(folder = ctrl_dir, restriction = "sevstd_", site = ctrl_site,
                                     label = "No support,\nat ventilated\nseverity", rank = 1,
                                     unit_factor = scale_v$sd_log_pfvc / scale_c$sd_log_pfvc)
    }
  }
  sf_found <- restrictions_in(fig_dir, site_name)
  sf_found <- sf_found[grepl("^sf[0-9.]+to[0-9.]+_$", sf_found)]
  sf_lo <- as.numeric(str_match(sf_found, "^sf([0-9.]+)to")[, 2])
  sf_hi <- as.numeric(str_match(sf_found, "to([0-9.]+)_$")[, 2])
  for (k in order(-sf_lo)) {                                          # mildest hypoxaemia first
    arms[[sf_found[k]]] <- list(folder = fig_dir, restriction = sf_found[k], site = site_name, rank = 3 + match(k, order(-sf_lo)) / 10,
                                label = paste0("Ventilated,\nSF ", if (sf_lo[k] == 0) paste0("<= ", sf_hi[k]) else paste0(sf_lo[k], "-", sf_hi[k])))
  }
  arm_rate <- map_dfr(arms, function(a) {
    est <- read_arm(a$folder, a$restriction, a$site)
    if (is.null(est)) return(NULL)
    f <- if (is.null(a$unit_factor)) 1 else a$unit_factor   # the ventilated cohort's unit
    rate_rows(est %>% mutate(estimate = estimate * f, lo = lo * f, hi = hi * f)) %>%
      mutate(arm = a$label, strain_rank = a$rank)
  })
  all_label <- if (nrow(arm_rate)) "Ventilated,\nall" else "Ventilated"
  rate_arms <- bind_rows(arm_rate, rate %>% select(-marker_lab, -s) %>% mutate(arm = all_label, strain_rank = 9)) %>%
    mutate(arm = factor(arm, unique(arm[order(strain_rank)])),
           marker_lab = factor(row_label(marker), row_order))
  if (nrow(arm_rate)) message("24_biotrauma_figures: arms in the rate panel: ",
                              paste(gsub("\n", " ", levels(rate_arms$arm)), collapse = "; "))
  has_ctrl <- any(grepl("^No support", levels(rate_arms$arm)))
  pm_b <- ggplot(rate_arms, aes(arm, e, colour = adjustment)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
    geom_linerange(aes(ymin = l, ymax = h), linewidth = 0.8, position = position_dodge(width = 0.55)) +
    geom_point(aes(shape = ok), size = 2.2, fill = "white", position = position_dodge(width = 0.55)) +
    facet_grid(marker_lab ~ ., scales = "free_y") + shared +
    guides(colour = "none") +
    labs(title = if (nrow(arm_rate)) "Rate per day, by cohort" else "Rate per day",
         subtitle = paste(c(if (has_ctrl) "controls: no respiratory support",
                            if (length(sf_found)) "ventilated by baseline SF",
                            if (!nrow(arm_rate)) "adjusted vs unadjusted"), collapse = "; "),
         x = NULL, y = "change per day toward injury")
  pm_c <- ggplot(trend, aes(day, p_harm, colour = adjustment)) +
    geom_hline(yintercept = 0.5, colour = "grey55") +
    geom_hline(yintercept = c(0.025, 0.975), linetype = 3, colour = "grey70") +
    geom_line(aes(linetype = ok), linewidth = 0.9) +
    geom_point(aes(shape = ok), size = 1.6) +
    facet_grid(marker_lab ~ ., drop = FALSE) +
    scale_x_continuous(breaks = day_breaks) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) + shared +
    theme(strip.text.y = element_text(angle = 0, hjust = 0)) +
    labs(title = "Ventilated: probability of more injury", subtitle = "posterior; dotted lines at 0.025 and 0.975",
         x = "days from the index", y = "P(contrast in the injury direction)")
  pm <- pm_a + pm_b + pm_c + plot_layout(widths = c(1.35, 0.35 + 0.25 * n_distinct(rate_arms$arm), 1), guides = "collect") +
    plot_annotation(
      tag_levels = "A",
      title = paste0("A smaller predicted lung, at the same VT/PBW, and organ-injury markers over ", JM_HORIZON,
                     " days of ventilation"),
      subtitle = paste0(site_name, ": joint model, death and extubation (in the controls, escalation of support) modelled; ", unit_lower,
                        ".\nA rate unmoved by adjustment for age, sex and race is not the age channel. ",
                        "Row counts are the ventilated cohort's. Hollow points and dashed lines did not converge.",
                        if (nzchar(sev_tag)) "\nSeverity-matched: each marker's floor is on its own anchor, so the rows are different patient subsets." else "")) &
    theme(legend.position = "top")
  ggsave(file.path(fig_dir, paste0("biotrauma_fig_main_", tag, ".pdf")), pm,
         width = 12 + 0.8 * n_distinct(rate_arms$arm), height = 2 + 2.2 * length(present))
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
  # A comparison needs two estimators. Keep the markers this run's joint model covers,
  # at the horizons where a second estimator exists for that marker. That drops hour 0
  # (no comparator can exist at the baseline) and, in a 7-day run, the days beyond the
  # comparators' 24/48/72 h, which the trend figure already shows for the joint model.
  group_by(marker) %>% filter(any(estimator == "Joint model (death modelled)")) %>%
  group_by(marker, horizon_h) %>% filter(n_distinct(estimator) >= 2) %>% ungroup() %>%
  mutate(estimator = factor(estimator, c("Joint model (death modelled)", "Longitudinal model only", "Fixed horizon, survivors only")),
         adjustment = factor(adjustment, c("adjusted", "unadjusted")),
         horizon = factor(paste(horizon_h, "h"), paste(sort(unique(horizon_h)), "h")),
         marker_lab = marker_label(marker))
if (nrow(est)) {
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
         width = max(9, 3 + 2.4 * n_distinct(est$horizon_h)), height = 2 + 1.6 * n_distinct(est$marker))
} else message("24_biotrauma_figures: no horizon has two estimators for a modelled marker; estimator figure skipped")

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

# ---- 5. the check behind figure 4's causal reading (pfvc form, unrestricted run):
#      ventilated vs no support: the ventilated rate, the no-support control read at the
#      ventilated severity, and their difference (27_control_comparison.R), all on the
#      ventilated cohort's unit (per SD of its log PFVC).
# Drawn toward injury: above zero = a smaller predicted lung does worse.
# A marker whose estimates are missing is left out of the panel, not drawn as zero.
if (MOD_FORM == "pfvc" && !nzchar(restrict_tag)) {
  did_tbl <- read_if(file.path(fig_dir, paste0("jm_control_did_pfvc_", h_suffix, "_", site_name, ".csv")))
  toward_injury <- function(m) if_else(worse[m] == "higher", -1, 1)   # a smaller lung is the negative of the per-SD rate
  check_order <- intersect(c("platelets", "bilirubin", "creatinine", "pressor_dose"),
                           unique(did_tbl$marker))
  check_label <- function(m) factor(lab[m], lab[check_order])
  checks <- list()
  if (!is.null(did_tbl) && nrow(did_tbl)) {
    did_rows <- did_tbl %>%
      transmute(marker, adjustment, s = toward_injury(marker),
                v_e = divergence_estimate_ventilated, v_sd = divergence_sd_ventilated, v_ok = divergence_rhat_ventilated <= 1.1,
                c_e = divergence_estimate_control, c_sd = divergence_sd_control, c_ok = divergence_rhat_control <= 1.1,
                d_e = did_estimate, d_lo = did_lo, d_hi = did_hi, d_ok = both_converged) %>%
      { bind_rows(
          transmute(., marker, adjustment, s, arm = "ventilated", e = v_e, l = v_e - 1.96 * v_sd, h = v_e + 1.96 * v_sd, ok = v_ok),
          transmute(., marker, adjustment, s, arm = "no\nsupport", e = c_e, l = c_e - 1.96 * c_sd, h = c_e + 1.96 * c_sd, ok = c_ok),
          transmute(., marker, adjustment, s, arm = "difference", e = d_e, l = d_lo, h = d_hi, ok = d_ok)) } %>%
      mutate(e = s * e, lo = pmin(s * l, s * h), hi = pmax(s * l, s * h), ok = coalesce(ok, FALSE),
             arm = factor(arm, c("ventilated", "no\nsupport", "difference")),
             marker_lab = check_label(marker))
    checks$did <- ggplot(did_rows, aes(arm, e, colour = adjustment)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
      geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.8, position = position_dodge(width = 0.5)) +
      geom_point(aes(shape = ok), size = 2.2, fill = "white", position = position_dodge(width = 0.5)) +
      facet_wrap(~ marker_lab, nrow = 1, scales = "free_y") +
      scale_colour_manual(values = okabe[c(1, 2)], name = NULL) +
      scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
      labs(title = "With and without a ventilator (difference-in-differences)",
           subtitle = "the control keeps every patient and is read at the ventilated cohort's severity and unit",
           x = NULL, y = "change per day toward injury\nper SD of log PFVC")
  }
  if (length(checks)) {
    ggsave(file.path(fig_dir, paste0("biotrauma_fig_checks_", tag, ".pdf")),
           wrap_plots(checks, ncol = 1) +
             plot_annotation(tag_levels = "A",
                             title = paste0(site_name, ": is the lung-size divergence the ventilator's?"),
                             subtitle = "95% intervals; hollow points did not converge (R-hat > 1.1)") &
             theme(legend.position = "top", axis.text.x = element_text(size = 8)),
           width = 3 + 2.6 * max(1, length(check_order)), height = 1.5 + 3.6 * length(checks))
    message("24_biotrauma_figures: checks figure (", paste(names(checks), collapse = ", "), ") -> biotrauma_fig_checks_", tag, ".pdf")
  } else message("24_biotrauma_figures: no DiD table for ", site_name, "; checks figure skipped")
}
message("24_biotrauma_figures: ", n_distinct(est$marker), " markers, ", n_distinct(est$estimator), " estimators -> ", fig_dir)
