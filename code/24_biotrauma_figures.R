# =============================================================================
# Script 24 (figures): the PFVC-level question, drawn from the aggregate tables
# =============================================================================
# Reads only the site's final/ CSVs (no patient rows), so it runs on any site's
# outputs and, with PBWPFVC_FIG_DIR, on a pooled folder.
#
#   biotrauma_fig_main_{tag}.pdf   THE figure, for a window with three or more
#       horizons (the 7-day run): one row per marker, the contrast toward injury
#       from hour 0 to the end of the window, the rate per day adjusted beside
#       unadjusted, and the posterior probability of harm by day
#   biotrauma_fig_checks_{tag}.pdf   the pfvc form's unrestricted run: every
#       outcome against each control (every no-support patient, and the no-support
#       patients hypoxemic on the index day): the markers' difference-in-differences,
#       and 60-day death (supplement/xsec_pfvc_age_control.R)
#   biotrauma_fig_channels_{tag}.pdf   the channels form only: the contrast per
#       GLI piece
#
# The level-contrast, divergence, trajectory, estimator and pressor figures were
# cut on 2026-09-24 (docs/output_manifest.md): the main figure carries the same
# estimates, and the estimator and pressor figures needed the 48-hour runner's
# tables.
#
# A lower PFVC is the negative of every log-marker estimate; the figures label
# the injury direction per marker so the eye does not have to flip signs.
#
# Usage: PBWPFVC_JM_MODIFIER=pfvc uvr run code/24_biotrauma_figures.R
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
FLIP_INJ  <- MOD_FORM %in% c("disc_level", "vtpfvc")   # a HIGHER value of these is the smaller lung
fig_dir   <- Sys.getenv("PBWPFVC_FIG_DIR", final_dir_for("injury"))
tag       <- paste0(restrict_tag, if (MOD_FORM != "disc") paste0(MOD_FORM, "_") else "", h_suffix, "_", site_name)
okabe <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9")
theme_set(theme_minimal(base_size = 11))
worse <- c(creatinine = "higher", platelets = "lower", bilirubin = "higher", sf = "lower", dp = "higher",
           ne_equiv_peak = "higher", pressor_dose = "higher", ne_equiv = "higher", any_pressor = "higher",
           osi = "higher", oi = "higher")
lab <- c(creatinine = "Creatinine", platelets = "Platelets", bilirubin = "Bilirubin", sf = "SF ratio\n(positive control:\nrecruitment, not injury)",
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
if (!nrow(lc0)) stop("no level contrasts for exposure ", SIZE_EX, " in tag ", tag)

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
  show_count <- function(x) formatC(x, big.mark = ",", format = "d")
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
  #                                         tables; PBWPFVC_JM_SEV_CENTER, 29_run_figure4.R)
  #        No support, hypoxemic            the same control restricted to index-day
  #                                         SF < 315, the ventilated cohort's own gate
  #                                         (the sevstd_sf0to315_ tables), so the two
  #                                         arms differ in ventilation, not hypoxaemia
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
    # only the severity-standardised controls: other control tables in the folder are ignored
    # its rate is per SD of log PFVC in the CONTROL's panel; on the figure it is put on
    # the ventilated cohort's unit (the SDs from 22_biotrauma_fit.R's jm_scale_* tables),
    # and without both SDs the control is not drawn
    scale_v <- read_if(file.path(fig_dir, paste0("jm_scale_", h_suffix, "_", site_name, ".csv")))
    scale_c <- read_if(file.path(ctrl_dir, paste0("jm_scale_", h_suffix, "_", ctrl_site, ".csv")))
    # each control arm: its restriction tag, its label and its place on the axis. Both
    # are read at the ventilated severity and put on the ventilated cohort's unit with
    # the whole control's SD of log PFVC (the hypoxemic fit keeps it, 22_biotrauma_fit.R)
    control_arms <- list(
      control           = list(restriction = "sevstd_", rank = 1,
                               label = "No support,\nat ventilated\nseverity"),
      control_hypoxemic = list(restriction = "sevstd_sf0to315_", rank = 2,
                               label = "No support,\nhypoxemic,\nat ventilated\nseverity"))
    for (arm_name in names(control_arms)) {
      control_arm <- control_arms[[arm_name]]
      if (!control_arm$restriction %in% restrictions_in(ctrl_dir, ctrl_site)) next
      if (is.null(scale_v) || is.null(scale_c) || SIZE_EX != "log_pfvc_sd") {
        message("24_biotrauma_figures: ", arm_name, " arm not drawn: its unit cannot be put on the ventilated cohort's ",
                "(jm_scale_* missing, or the form is not pfvc)")
        next
      }
      arms[[arm_name]] <- list(folder = ctrl_dir, restriction = control_arm$restriction, site = ctrl_site,
                               label = control_arm$label, rank = control_arm$rank,
                               unit_factor = scale_v$sd_log_pfvc / scale_c$sd_log_pfvc)
    }
  }
  sf_found <- restrictions_in(fig_dir, site_name)
  sf_found <- sf_found[grepl("^sf[0-9.]+to[0-9.]+_$", sf_found)]
  sf_lo <- as.numeric(str_match(sf_found, "^sf([0-9.]+)to")[, 2])
  sf_hi <- as.numeric(str_match(sf_found, "to([0-9.]+)_$")[, 2])
  for (k in order(-sf_lo)) {                                          # mildest hypoxaemia first
    arms[[sf_found[k]]] <- list(folder = fig_dir, restriction = sf_found[k], site = site_name, rank = 4 + match(k, order(-sf_lo)) / 10,
                                label = paste0("Ventilated,\nSF ", if (sf_lo[k] == 0) paste0("< ", sf_hi[k]) else paste0(sf_lo[k], "-", sf_hi[k])))
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

# ---- 5. the check behind figure 4's causal reading (pfvc form, unrestricted run):
#      the ventilated rate, a control's rate read at the ventilated severity, and their
#      difference (27_control_comparison.R), all on the ventilated cohort's unit (per SD
#      of its log PFVC). One panel per control:
#        A  every no-support patient (jm_control_did_*): the larger sample, but most of
#           the control is not hypoxemic, so the difference contrasts ventilation and
#           hypoxaemia together
#        B  the no-support patients hypoxemic on the index day, SF < 315
#           (jm_hypoxemic_control_did_*): the arms differ in ventilation alone
#      and beneath each, 60-day all-cause death against the same control: the hazard
#      ratio per SD lower log PFVC in each cohort and their difference, read at the
#      ventilated severity (pfvc_age_control_contrast_*, final/supplement/). The
#      oxygen saturation index has no control:
#      it needs a mean airway pressure.
# Drawn toward injury: above zero = a smaller predicted lung does worse.
# A marker whose estimates are missing is left out of the panel, not drawn as zero.
if (MOD_FORM == "pfvc" && !nzchar(restrict_tag)) {
  did_stub <- paste0("pfvc_", h_suffix, "_", site_name, ".csv")
  did_tables <- list(
    all = list(table = read_if(file.path(fig_dir, paste0("jm_control_did_", did_stub))),
               control = "no support,\nall",
               title = "Against every no-support patient",
               subtitle = "read at the ventilated cohort's severity and unit"),
    hypoxemic = list(table = read_if(file.path(fig_dir, paste0("jm_hypoxemic_control_did_", did_stub))),
                     control = "no support,\nhypoxemic",
                     title = "Against no-support patients hypoxemic on the index day",
                     subtitle = "SF < 315, so the arms differ in ventilation, not hypoxaemia"))
  toward_injury <- function(m) if_else(worse[m] == "higher", -1, 1)   # a smaller lung is the negative of the per-SD rate
  check_order <- intersect(c("platelets", "bilirubin", "creatinine", "pressor_dose"),
                           unique(unlist(map(did_tables, ~ .x$table$marker))))
  check_label <- function(m) factor(lab[m], lab[check_order])
  did_panel <- function(did_tbl, control_label, title, subtitle) {
    did_rows <- did_tbl %>%
      transmute(marker, adjustment, s = toward_injury(marker),
                v_e = divergence_estimate_ventilated, v_sd = divergence_sd_ventilated, v_ok = divergence_rhat_ventilated <= 1.1,
                c_e = divergence_estimate_control, c_sd = divergence_sd_control, c_ok = divergence_rhat_control <= 1.1,
                d_e = did_estimate, d_lo = did_lo, d_hi = did_hi, d_ok = both_converged) %>%
      { bind_rows(
          transmute(., marker, adjustment, s, arm = "ventilated", e = v_e, l = v_e - 1.96 * v_sd, h = v_e + 1.96 * v_sd, ok = v_ok),
          transmute(., marker, adjustment, s, arm = control_label, e = c_e, l = c_e - 1.96 * c_sd, h = c_e + 1.96 * c_sd, ok = c_ok),
          transmute(., marker, adjustment, s, arm = "difference", e = d_e, l = d_lo, h = d_hi, ok = d_ok)) } %>%
      mutate(e = s * e, lo = pmin(s * l, s * h), hi = pmax(s * l, s * h), ok = coalesce(ok, FALSE),
             arm = factor(arm, c("ventilated", control_label, "difference")),
             marker_lab = check_label(marker))
    ggplot(did_rows, aes(arm, e, colour = adjustment)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
      geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.8, position = position_dodge(width = 0.5)) +
      geom_point(aes(shape = ok), size = 2.2, fill = "white", position = position_dodge(width = 0.5)) +
      facet_wrap(~ marker_lab, nrow = 1, scales = "free_y", drop = FALSE) +
      scale_colour_manual(values = okabe[c(1, 2)], name = NULL) +
      scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
      labs(title = title, subtitle = subtitle, x = NULL, y = "change per day toward injury\nper SD of log PFVC")
  }
  # 60-day death against each control, from the supplement's contrast table: log hazard
  # ratio per SD of log PFVC (the ventilated cohort's SD) in each cohort, and the
  # difference, turned toward harm (a smaller lung) and to ventilated minus control
  death_contrast <- if (nzchar(Sys.getenv("PBWPFVC_FIG_DIR", ""))) NULL else
    read_if(file.path(config$final_root, "supplement", paste0("pfvc_age_control_contrast_", site_name, ".csv")))
  DEATH_POPULATIONS <- c(all = "everyone", hypoxemic = "hypoxemic at the index (SF < 315)")
  death_panel <- function(population, control_label) {
    if (is.null(death_contrast)) return(NULL)
    rows <- death_contrast %>%
      filter(population == !!population, outcome == "60-day death, all",
             severity == "standardised to ventilated severity",
             quantity %in% c("Ventilated", "No support", "no support minus ventilated"))
    if (!nrow(rows)) return(NULL)
    rows <- rows %>%
      transmute(adjustment,
                arm = factor(recode(quantity, Ventilated = "ventilated", `No support` = control_label,
                                    `no support minus ventilated` = "difference"),
                             c("ventilated", control_label, "difference")),
                # toward harm: per SD LOWER log PFVC; the difference row is already
                # no support minus ventilated, which is ventilated minus control toward harm
                e = if_else(quantity == "no support minus ventilated", log_ratio, -log_ratio),
                lo = e - 1.96 * se, hi = e + 1.96 * se, marker_lab = "60-day death")
    ggplot(rows, aes(arm, e, colour = adjustment)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
      geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.8, position = position_dodge(width = 0.5)) +
      geom_point(size = 2.2, position = position_dodge(width = 0.5)) +
      facet_wrap(~ marker_lab) +
      scale_colour_manual(values = okabe[c(1, 2)], name = NULL) +
      labs(title = NULL, subtitle = "60-day all-cause death, at the ventilated severity (Cox)",
           x = NULL, y = "log hazard ratio toward harm\nper SD lower log PFVC")
  }
  checks <- list()
  for (control_name in names(did_tables)) {
    did_entry <- did_tables[[control_name]]
    if (!is.null(did_entry$table) && nrow(did_entry$table))
      checks[[control_name]] <- did_panel(did_entry$table, did_entry$control, did_entry$title, did_entry$subtitle)
    death <- death_panel(DEATH_POPULATIONS[[control_name]], did_entry$control)
    if (!is.null(death)) checks[[paste0(control_name, "_death")]] <- death
  }
  if (length(checks)) {
    ggsave(file.path(fig_dir, paste0("biotrauma_fig_checks_", tag, ".pdf")),
           wrap_plots(checks, ncol = 1) +
             plot_annotation(tag_levels = "A",
                             title = paste0(site_name, ": is the lung-size divergence the ventilator's?"),
                             subtitle = "95% intervals; hollow points did not converge (R-hat > 1.1)") &
             theme(legend.position = "top", axis.text.x = element_text(size = 8)),
           width = max(8, 3 + 2.6 * length(check_order)), height = 1.5 + 3.6 * length(checks))
    message("24_biotrauma_figures: checks figure (", paste(names(checks), collapse = ", "), ") -> biotrauma_fig_checks_", tag, ".pdf")
  } else message("24_biotrauma_figures: no DiD table for ", site_name, "; checks figure skipped")
}
message("24_biotrauma_figures: ", n_distinct(lc0$marker), " markers -> ", fig_dir)
