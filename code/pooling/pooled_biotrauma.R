# =============================================================================
# Pooled biotrauma results across sites (run centrally, like pooled_estimates.R)
# =============================================================================
# Discovers each site's aggregate biotrauma outputs under a results root that
# holds one subfolder per site (each site's final/ renamed to the site name,
# PBWPFVC_RESULTS_ROOT, default results/). The pooled tables and forests go to
# the All sites/ subfolder, which is excluded from discovery.
#
# How the pooling is done (2026-09-22):
#   units   log_pfvc_sd is standardised inside each site's own panel, so a per-SD
#           estimate means a different lung-size difference at every site. Every
#           per-SD estimate is converted to per 0.1 log units of PFVC using that
#           site's jm_scale_{h}_{site}.csv, and the unit is a grouping key, so a
#           site without a scale table is never pooled with a converted one.
#           VT/PFVC (vtpfvc_c) is per percentage point everywhere and is pooled
#           as it stands.
#   method  common-effect inverse variance, because this project has two or three
#           sites: a random-effects tau2 from k = 2 is not estimable in any
#           useful sense. Heterogeneity (Q p, I2, tau2) is reported beside every
#           pooled row, and a REML + Knapp-Hartung estimate is added as a
#           sensitivity from k = 3. A failed random-effects fit is reported in
#           re_status, never silently replaced.
#   gate    a site enters a pool only if its own chain converged for that term
#           (rhat <= 1.1); the difference-in-differences uses its both_converged
#           flag. The hazard blocks of the 7-day fits do not converge, so the
#           association hazard ratios are pooled without a gate and read as
#           descriptive.
#
# What is pooled (each row is estimate + standard error per site):
#   jm_level_contrast_*    the joint model's marker difference per SD of log PFVC
#                          (or log PBW/PFVC) at 24/48/72 h; SE = (hi - lo) / 3.92
#                          from the posterior interval; the log-odds scale for the
#                          any-pressor part
#   jm_association_hr_*    log hazard ratio of death / extubation per SD of the
#                          current marker (value) or per unit slope
#   jm_estimates_*         the key longitudinal terms (dose slope, dose x
#                          discordance, dose x age) with posterior SD as SE
#   injury_at_horizon_*    the fixed-horizon comparator (estimate, se)
#   quick_lme_*            the longitudinal-submodel-only contrast, SE from the CI
#   injury_channels_*, quick_channels_*            the channel decomposition
#   injury_dose_channels_*, injury_sf_channels_*,  the channel supports
#   injury_negctrl_*, quick_dose_channels_*,       (coefficients only; the
#   quick_sf_channels_*                             nested LR tables are not pooled)
#   jm_control_did_*       figure 4's difference-in-differences: the ventilated
#                          divergence minus the no-support control's, per day
# Joint-model tables are keyed by the panel horizon in the file tag (24h/48h/72h)
# as well as the contrast horizon, so one site's three panels are never pooled as
# three sites. Every pooled row carries k, I2, tau2 and the per-site estimates it
# was built from. Sites are anonymised with utils/site_anonymization.R when present.
#
# Usage: PBWPFVC_RESULTS_ROOT=/path/to/results Rscript code/pooling/pooled_biotrauma.R
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(here); library(metafor) })
rm(list = ls())

root <- Sys.getenv("PBWPFVC_RESULTS_ROOT", here("results"))
if (!dir.exists(root)) stop("results root not found: ", root)
out_dir <- file.path(root, "All sites"); dir.create(out_dir, showWarnings = FALSE)
sites <- setdiff(list.dirs(root, recursive = FALSE, full.names = FALSE), "All sites")
sites <- sites[!startsWith(sites, ".")]
message("Sites: ", paste(sites, collapse = ", "))
# site labels: the project's anonymization (built from each site's cohort size)
# when it can be built, otherwise the folder names
anon <- identity
site_levels <- sort(sites)          # largest cohort first once the alias table is built
if (file.exists(here("utils", "site_anonymization.R"))) {
  source(here("utils", "site_anonymization.R"))
  alias_tbl <- tryCatch(build_site_aliases(file.path(root, sites)), error = function(e) NULL)
  if (!is.null(alias_tbl)) {
    aliases <- alias_tbl$aliases
    anon <- function(x) anonymize_site(x, aliases)
    site_levels <- alias_tbl$table$site_label
  } else message("site anonymization unavailable for this root (", "no cohort sizes); using folder names")
}
# the cross-sectional pooling's palette and ordering (pooled_estimates.R): Okabe-Ito
# by cohort, largest first, so a cohort keeps its colour across every pooled figure;
# the pooled estimate is a black diamond on a row at the foot of each panel
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
POOLED_LABEL <- "Pooled (CE)"        # common effect (see pool_one)
site_colors <- if (length(site_levels) > length(okabe_ito))
  setNames(grDevices::colorRampPalette(okabe_ito)(length(site_levels)), site_levels) else
  setNames(okabe_ito[seq_along(site_levels)], site_levels)
forest_palette <- c(site_colors, setNames("#000000", POOLED_LABEL))
okabe <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#000000")

# Injury is drawn upward in every pooled figure: a marker's estimate is per 0.1 log
# units MORE predicted lung, so it is negated for a marker that is worse when higher,
# and kept for one that is worse when lower. The result reads "change toward injury
# per 10% smaller predicted lung" for every marker, and the strip names the marker
# in words rather than by its column name.
MARKER_LABELS <- c(platelets = "Platelets", creatinine = "Creatinine", bilirubin = "Bilirubin",
                   pressor_dose = "Vasopressor dose\n(NE-equivalents per kg)", any_pressor = "Any vasopressor\n(log-odds)",
                   ne_equiv_peak = "Peak vasopressor\n(NE-equivalents per kg)", osi = "Oxygen saturation index",
                   sf = "SpO2 / FiO2", dp = "Driving pressure")
MARKER_WORSE  <- c(platelets = "lower", creatinine = "higher", bilirubin = "higher", pressor_dose = "higher",
                   any_pressor = "higher", ne_equiv_peak = "higher", osi = "higher", sf = "lower", dp = "higher")
toward_injury <- function(d) {
  unknown <- setdiff(unique(d$marker), names(MARKER_WORSE))
  if (length(unknown)) stop("no injury direction for marker(s): ", paste(unknown, collapse = ", "))
  d %>% mutate(sign = if_else(MARKER_WORSE[marker] == "higher", -1, 1),
               estimate = sign * estimate, lo_new = pmin(sign * lo, sign * hi), hi = pmax(sign * lo, sign * hi), lo = lo_new,
               marker_label = factor(MARKER_LABELS[marker], levels = MARKER_LABELS[names(MARKER_LABELS) %in% marker])) %>%
    select(-sign, -lo_new)
}
# One forest in the cross-sectional style: sites as coloured points with capped
# intervals, the pooled estimate as a black diamond at the foot, one row per marker
# and one column per `column`, each column on its own x scale.
draw_forest <- function(d, title, subtitle, x_label) {
  d <- d %>% mutate(site = factor(site, levels = c(POOLED_LABEL, rev(site_levels))),
                    kind = if_else(site == POOLED_LABEL, "Pooled", "Site"))
  ggplot(d, aes(x = estimate, y = site, colour = site)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, orientation = "y") +
    geom_point(aes(size = kind, shape = kind)) +
    facet_grid(marker_label ~ column, scales = "free_x", drop = FALSE) +
    scale_colour_manual(values = forest_palette, guide = "none") +
    scale_shape_manual(values = c(Site = 16, Pooled = 18), guide = "none") +
    scale_size_manual(values = c(Site = 2.3, Pooled = 3.4), guide = "none") +
    labs(title = title, subtitle = subtitle, x = x_label, y = "Cohort") +
    theme_minimal(base_size = 11) +
    theme(strip.text.x = element_text(face = "bold"),
          strip.text.y = element_text(face = "bold", angle = 0),
          panel.spacing = unit(0.6, "lines"),
          panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.5))
}
forest_height <- function(d) 2 + 0.28 * n_distinct(d$marker_label) * (n_distinct(d$site) + 1.5)

# read one file family from every site, tagging the site; tolerant of absent files
read_family <- function(pattern) {
  map_dfr(sites, function(s) {
    # a returned final/ is sorted by block; the tables pooled here are in injury/. controls/
    # is never listed, so a control cohort cannot enter a pool of the ventilated one. The
    # site folder itself is listed too, for a flat (older) return.
    fs <- list.files(file.path(root, s, c("", "injury")), pattern = pattern, full.names = TRUE)
    if (!length(fs)) return(NULL)
    map_dfr(fs, function(f) read_csv(f, show_col_types = FALSE, guess_max = 1e5) %>%
              mutate(site = s, file = basename(f), .before = 1))
  })
}

# Pool one estimate/SE set into one row. The primary pool is common-effect
# (inverse-variance): with the two or three sites this project has, a
# random-effects tau2 is estimated from too few studies to mean anything, and
# Knapp-Hartung on k = 2 gives a t(1) interval so wide it says nothing. The
# random-effects fit is reported beside it as a sensitivity from k = 3, and its
# heterogeneity (Q, I2, tau2) is reported at every k so a reader can see when the
# sites disagree. A failed random-effects fit is recorded in re_status, never
# swapped silently for another estimator.
pool_one <- function(d) {
  d <- d %>% filter(is.finite(estimate), is.finite(se), se > 0)
  k <- nrow(d)
  if (k == 0) return(NULL)
  one_site <- tibble(k = 1L, method = "single site", pooled = d$estimate[1], se = d$se[1],
                     lo = d$estimate[1] - 1.96 * d$se[1], hi = d$estimate[1] + 1.96 * d$se[1],
                     q_p = NA_real_, i2 = NA_real_, tau2 = NA_real_,
                     re_pooled = NA_real_, re_lo = NA_real_, re_hi = NA_real_, re_status = "k < 3")
  if (k == 1) return(one_site %>% mutate(sites = d$site, site_estimates = as.character(signif(d$estimate, 4))))
  fe <- rma(yi = d$estimate, sei = d$se, method = "FE")
  re <- tryCatch(rma(yi = d$estimate, sei = d$se, method = "REML", test = "knha"),
                 error = function(e) conditionMessage(e))
  re_ok <- inherits(re, "rma")
  tibble(k = k, method = "common effect", pooled = as.numeric(fe$b), se = fe$se,
         lo = fe$ci.lb, hi = fe$ci.ub,
         q_p = fe$QEp, i2 = if (re_ok) re$I2 else NA_real_, tau2 = if (re_ok) re$tau2 else NA_real_,
         re_pooled = if (re_ok && k >= 3) as.numeric(re$b) else NA_real_,
         re_lo = if (re_ok && k >= 3) re$ci.lb else NA_real_,
         re_hi = if (re_ok && k >= 3) re$ci.ub else NA_real_,
         re_status = if (!re_ok) paste("random-effects fit failed:", re) else if (k < 3) "k < 3" else "REML, Knapp-Hartung",
         sites = paste(d$site, collapse = ";"), site_estimates = paste(signif(d$estimate, 4), collapse = ";"))
}
pool_by <- function(d, ...) d %>% group_by(...) %>% group_modify(~ pool_one(.x)) %>% ungroup()
# the modifier form and the panel horizon from a joint-model file name
jm_form  <- function(file, family, default) str_match(file, paste0("^", family, "_(?:(\\w+?)_)?(\\d+[hd])_"))[, 2] %>% replace_na(default)
jm_panel <- function(file, family) str_match(file, paste0("^", family, "_(?:(\\w+?)_)?(\\d+[hd])_"))[, 3]
CHANNELS <- c("ch_height", "ch_age", "ch_sex", "ch_race")
# an interaction is written in whichever order the model formula produced it
# (log_pfvc_sd:vent_day at one site, vent_day:log_pfvc_sd at another): sort the
# components so the same term from two sites lands in one pool
canonical_term <- function(term) map_chr(str_split(term, ":"), ~ paste(sort(.x), collapse = ":"))
RHAT_MAX <- 1.1   # a site's estimate enters a pool only if its own chain converged

# --- units. log_pfvc_sd is standardised inside each site's own panel, so a
# per-SD estimate means a different lung-size difference at every site and the
# raw numbers must not be pooled. jm_scale_{h}_{site}.csv carries each site's
# SD of log PFVC; every per-SD estimate is converted to PER_LOG_PFVC log units
# (0.1 log units is about a 10% smaller predicted lung). VT/PFVC needs no
# conversion: vtpfvc_c is per percentage point of predicted FVC everywhere.
PER_LOG_PFVC <- 0.1
sc <- read_family("^jm_scale_.*\\.csv$")
pfvc_sd <- if (nrow(sc)) sc %>% filter(cohort == "imv") %>%
  transmute(site, panel_h = paste0(horizon_days, "d"), sd_log_pfvc) %>% distinct() else
  tibble(site = character(), panel_h = character(), sd_log_pfvc = double())
# Scale the estimates that are per SD of log PFVC (`per_sd`, one value per row)
# to per PER_LOG_PFVC log units, leaving every other row as it is with the unit
# the caller names. A site with no scale table keeps its per-SD estimate and is
# labelled, so it never joins a pool of converted ones: `unit` is a grouping key
# in every pool below.
to_log_units <- function(d, per_sd, other_unit) {
  d %>% mutate(.per_sd = per_sd, .other_unit = other_unit) %>%
    left_join(pfvc_sd, by = c("site", "panel_h")) %>%
    mutate(scale_factor = if_else(.per_sd & !is.na(sd_log_pfvc), PER_LOG_PFVC / sd_log_pfvc, 1),
           unit = case_when(!.per_sd ~ .other_unit,
                            !is.na(sd_log_pfvc) ~ paste0("per ", PER_LOG_PFVC, " log PFVC"),
                            TRUE ~ "per site SD of log PFVC (not harmonised)"),
           across(any_of(c("estimate", "se", "lo", "hi")), ~ .x * scale_factor)) %>%
    select(-.per_sd, -.other_unit)
}

pooled <- list()

# --- 1. joint-model level contrasts (the PFVC-level question)
lc <- read_family("^jm_level_contrast_.*\\.csv$")
if (nrow(lc)) {
  lc <- lc %>% mutate(se = (hi - lo) / 3.92,
                      grid = if ("grid" %in% names(lc)) grid else NA_character_,
                      form = jm_form(file, "jm_level_contrast", "pfvc"), panel_h = jm_panel(file, "jm_level_contrast")) %>%
    to_log_units(per_sd = .$exposure == "log_pfvc_sd",
                 other_unit = if_else(.$exposure == "vtpfvc_c", "per point of VT/PFVC", "per site unit of the exposure"))
  pooled$level_contrast <- pool_by(lc, marker, model, adjustment, exposure, unit, panel_h, horizon_h, grid, form) %>%
    mutate(scale = "log marker; log-odds for any_pressor")
}

# --- 2. association hazard ratios (Q2)
ah <- read_family("^jm_association_hr_.*\\.csv$")
if (nrow(ah)) {
  ah <- ah %>% mutate(estimate = log_hr, se = (log_hr_hi - log_hr_lo) / 3.92,
                      grid = if ("grid" %in% names(ah)) grid else NA_character_,
                      form = jm_form(file, "jm_association_hr", "disc"), panel_h = jm_panel(file, "jm_association_hr"))
  pooled$association <- pool_by(ah, marker, model, adjustment, kind, cause, panel_h, grid, form) %>%
    mutate(hr = exp(pooled), hr_lo = exp(lo), hr_hi = exp(hi))
}

# --- 3. key longitudinal terms from the estimates tables
es <- read_family("^jm_estimates_.*\\.csv$")
if (nrow(es)) {
  # written with the interaction components in either order; canonical_term sorts them
  key <- canonical_term(c("l_vtpbw_within", "l_vtpbw_within:ldisc_c", "l_vtpbw_within:age10_c",
           "log_pfvc_sd", "log_pfvc_sd:vent_day", "ldisc_sd", "ldisc_sd:vent_day",
           "vtpfvc_c", "vtpfvc_c:vent_day",   # the VT/PFVC companion (figure 4, step 7)
           "ers_pfvc_0:l_vtpbw_within", "vtpbw_pt_mean",
           CHANNELS, paste0(CHANNELS, ":vent_day")))
  es <- es %>% mutate(term = canonical_term(term)) %>% filter(block == "longitudinal", term %in% key) %>%
    mutate(se = sd, grid = if ("grid" %in% names(es)) grid else NA_character_,
           form = jm_form(file, "jm_estimates", "disc"), panel_h = jm_panel(file, "jm_estimates"))
  # a site enters only if its own chain converged for that term; the dropped rows
  # are named, not counted, so a reader can see which marker lost which site
  dropped <- es %>% filter(!is.na(rhat), rhat > RHAT_MAX) %>%
    transmute(what = paste0(site, " ", marker, " ", adjustment, " ", term, " (", form, " ", panel_h,
                            ", rhat ", round(rhat, 2), ")"))
  if (nrow(dropped)) message("longitudinal terms dropped for rhat > ", RHAT_MAX, ":\n  ",
                             paste(dropped$what, collapse = "\n  "))
  es <- es %>% filter(is.na(rhat) | rhat <= RHAT_MAX) %>%
    to_log_units(per_sd = str_detect(.$term, "log_pfvc_sd"),
                 other_unit = if_else(str_detect(.$term, "vtpfvc_c"), "per point of VT/PFVC", "per site unit of the term"))
  pooled$longitudinal_terms <- pool_by(es, marker, model, adjustment, term, unit, panel_h, grid, form)
}

# --- 4. fixed-horizon comparator
ih <- read_family("^injury_at_horizon_.*\\.csv$")   # list.files takes POSIX regex: no lookahead
if (nrow(ih)) ih <- ih %>% filter(!grepl("^injury_at_horizon_counts_", file))
if (nrow(ih)) pooled$injury_at_horizon <- pool_by(ih, marker, horizon_h, outcome_type, exposure, adjustment)

# --- 5. quick LME (no death correction); older horizon-tagged files are stale copies
ql <- read_family("^quick_lme_.*\\.csv$")
if (nrow(ql)) ql <- ql %>% filter(!grepl("^quick_lme_[a-z_]+_\\d+h_", file))
if (nrow(ql)) {
  ql <- ql %>% mutate(se = (hi - lo) / 3.92,
                      model_horizon_h = if ("model_horizon_h" %in% names(ql)) model_horizon_h else horizon_h)
  pooled$quick_lme <- pool_by(ql, marker, exposure, adjustment, model_horizon_h)
}

# --- 6. channel decomposition and its supports (coefficients with a standard error)
ic <- read_family("^injury_channels_.*\\.csv$")
if (nrow(ic)) pooled$injury_channels <- pool_by(ic, marker, horizon_h, outcome_type, exposure, channel)
qc <- read_family("^quick_channels_.*\\.csv$")
if (nrow(qc)) pooled$quick_channels <- pool_by(qc, marker, exposure, channel, model_horizon_h)
for (fam in c("injury_dose_channels", "injury_sf_channels", "injury_negctrl")) {
  x <- read_family(paste0("^", fam, "_.*\\.csv$"))
  if (nrow(x) && "estimate" %in% names(x)) {
    x <- x %>% filter(!is.na(estimate))
    pooled[[fam]] <- if (fam == "injury_negctrl") pool_by(x, marker, outcome, model, term) else
      pool_by(x, marker, horizon_h, outcome_type, model, term)
  }
}
for (fam in c("quick_dose_channels", "quick_sf_channels")) {
  x <- read_family(paste0("^", fam, "_.*\\.csv$"))
  if (nrow(x) && "estimate" %in% names(x))
    pooled[[fam]] <- pool_by(x %>% filter(!is.na(estimate)), marker, exposure, model_horizon_h, model, term)
}

# --- 7. the figure-4 causal support: the difference-in-differences against the
#        no-support control, per SD of log PFVC per day in the ventilated cohort's
#        units, converted by to_log_units.
did <- read_family("^jm_control_did_.*\\.csv$")
if (nrow(did)) {
  did <- did %>% transmute(site, marker, adjustment, form = jm_form(file, "jm_control_did", "pfvc"),
                           panel_h = jm_panel(file, "jm_control_did"),
                           estimate = did_estimate, se = did_sd, both_converged)
  if (any(!did$both_converged))
    message("difference in differences dropped, an arm did not converge:\n  ",
            did %>% filter(!both_converged) %>%
              transmute(what = paste(site, marker, adjustment)) %>% pull(what) %>% paste(collapse = "\n  "))
  did <- did %>% filter(both_converged) %>%
    to_log_units(per_sd = TRUE, other_unit = NA_character_)
  if (nrow(did)) pooled$control_did <- pool_by(did, marker, adjustment, unit, panel_h, form) %>%
    mutate(scale = "ventilated minus no-support divergence, log marker per day")
}

# --- write
for (nm in names(pooled)) {
  write_csv(pooled[[nm]], file.path(out_dir, paste0("pooled_biotrauma_", nm, ".csv")))
  message(nm, ": ", nrow(pooled[[nm]]), " pooled rows")
}

# --- figure 4 pooled: the divergence and the control contrast that tests it, site by
#     site and pooled, one file per adjustment. Main model, the daily panel, PFVC
#     units only: the VT/PFVC divergence is on its own scale and is not drawn here.
DIVERGENCE <- c(pfvc = "log_pfvc_sd:vent_day", vtpfvc = canonical_term("vtpfvc_c:vent_day"))
FIG4_COLUMNS <- c(divergence = "Ventilated cohort:\ndivergence by predicted lung size",
                  did = "Ventilated minus no-support control\n(difference-in-differences)")
# the column label is read from the calling environment (.env) so that no table column
# of the same name can shadow it
fig4_rows <- function(d, column_key, keep = TRUE) {
  if (is.null(d) || !nrow(d)) return(NULL)
  d %>% filter(keep, grepl("d$", panel_h)) %>%
    transmute(marker, adjustment, unit, column = FIG4_COLUMNS[[.env$column_key]], site = anon(site),
              estimate, lo = estimate - 1.96 * se, hi = estimate + 1.96 * se)
}
fig4_pooled <- function(d, column_key, keep = TRUE) {
  if (is.null(d) || !nrow(d)) return(NULL)
  d %>% filter(keep, grepl("d$", panel_h)) %>%
    transmute(marker, adjustment, unit, column = FIG4_COLUMNS[[.env$column_key]], site = POOLED_LABEL,
              estimate = pooled, lo, hi)
}
fd4 <- bind_rows(
  fig4_rows(es, "divergence", es$term == DIVERGENCE[["pfvc"]] & es$model == "main"),
  fig4_pooled(pooled$longitudinal_terms, "divergence",
              pooled$longitudinal_terms$term == DIVERGENCE[["pfvc"]] & pooled$longitudinal_terms$model == "main"),
  fig4_rows(did, "did"), fig4_pooled(pooled$control_did, "did"))
if (!is.null(fd4) && nrow(fd4)) {
  fd4 <- fd4 %>% filter(startsWith(unit, "per ")) %>%
    mutate(column = factor(column, levels = unname(FIG4_COLUMNS))) %>%
    toward_injury()
  for (adj in c("adjusted", "unadjusted")) {
    d <- fd4 %>% filter(adjustment == adj)
    if (!nrow(d)) next
    p4 <- draw_forest(d,
      title = paste0("Figure 4 pooled: divergence by predicted lung size over 7 days of ventilation (", adj, ")"),
      subtitle = paste0("log marker per day per 10% smaller predicted lung (", PER_LOG_PFVC, " log units of PFVC), 95% CI; ",
                        "injury upward for every marker;\nblack diamond = common-effect pooled estimate; dashed line = null (0)"),
      x_label = "change per day toward injury per 10% smaller predicted lung (95% CI)")
    ggsave(file.path(out_dir, paste0("pooled_biotrauma_figure4", if (adj == "unadjusted") "_unadjusted" else "", ".pdf")),
           p4, width = 11, height = forest_height(d), limitsize = FALSE)
  }
}

# --- forests: the PFVC-level contrast per marker and horizon, per site and pooled,
#     one file per adjustment
if (nrow(lc) && any(lc$form == "pfvc" & lc$panel_h == paste0(lc$horizon_h, "h"))) {
  # the primary read: the pfvc form's contrast at each contrast horizon, from the panel of the same length
  fd <- lc %>% filter(exposure == "log_pfvc_sd", model == "main", form == "pfvc", panel_h == paste0(horizon_h, "h")) %>%
    transmute(marker, adjustment, horizon_h, unit, site = anon(site), estimate, lo, hi) %>%
    bind_rows(pooled$level_contrast %>% filter(exposure == "log_pfvc_sd", model == "main", form == "pfvc", panel_h == paste0(horizon_h, "h")) %>%
                transmute(marker, adjustment, horizon_h, unit, site = POOLED_LABEL, estimate = pooled, lo, hi)) %>%
    filter(startsWith(unit, "per ")) %>%
    mutate(column = factor(paste0("Marker difference at ", horizon_h, " h"),
                           levels = paste0("Marker difference at ", sort(unique(horizon_h)), " h"))) %>%
    toward_injury()
  for (adj in c("adjusted", "unadjusted")) {
    d <- fd %>% filter(adjustment == adj)
    if (!nrow(d)) next
    p <- draw_forest(d,
      title = paste0("Marker difference by predicted lung size at the horizon, joint model (", adj, ")"),
      subtitle = paste0("log marker per 10% smaller predicted lung (", PER_LOG_PFVC, " log units of PFVC), 95% CI, death before the horizon modelled; ",
                        "injury upward for every marker;\nblack diamond = common-effect pooled estimate; dashed line = null (0); ",
                        "any-vasopressor rows are log-odds"),
      x_label = "difference toward injury per 10% smaller predicted lung (95% CI)")
    ggsave(file.path(out_dir, paste0("pooled_biotrauma_level_contrast", if (adj == "unadjusted") "_unadjusted" else "", ".pdf")),
           p, width = 4 + 3.5 * n_distinct(d$horizon_h), height = forest_height(d), limitsize = FALSE)
  }
}
message("pooled_biotrauma complete -> ", out_dir)
