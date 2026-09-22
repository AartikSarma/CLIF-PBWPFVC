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
#   jm_pre_placebo_*       the pre-intubation placebo: the divergence in the days
#                          before the first IMV record, beside the joint model's
#                          post-intubation rate, pooled as separate periods
#   jm_pre_stacked_*       the stacked within-patient placebo: before, after and
#                          the change at intubation in the same patients
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
if (file.exists(here("utils", "site_anonymization.R"))) {
  source(here("utils", "site_anonymization.R"))
  aliases <- tryCatch(build_site_aliases(file.path(root, sites))$aliases, error = function(e) NULL)
  if (!is.null(aliases)) anon <- function(x) anonymize_site(x, aliases) else
    message("site anonymization unavailable for this root (", "no cohort sizes); using folder names")
}
okabe <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#000000")

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
           across(any_of(c("estimate", "se")), ~ .x * scale_factor)) %>%
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
  # a site enters only if its own chain converged for that term
  dropped <- sum(!is.na(es$rhat) & es$rhat > RHAT_MAX)
  if (dropped) message("longitudinal terms: ", dropped, " site-rows dropped for rhat > ", RHAT_MAX)
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

# --- 7. the figure-4 causal supports: the difference-in-differences against the
#        no-support control, and the pre-intubation placebo (two-sample and
#        stacked within-patient). All three are per SD of log PFVC per day in the
#        ventilated cohort's units, so all three are converted by to_log_units.
did <- read_family("^jm_control_did_.*\\.csv$")
if (nrow(did)) {
  did <- did %>% transmute(site, marker, adjustment, form = jm_form(file, "jm_control_did", "pfvc"),
                           panel_h = jm_panel(file, "jm_control_did"),
                           estimate = did_estimate, se = did_sd, both_converged)
  if (any(!did$both_converged))
    message("difference in differences: ", sum(!did$both_converged), " site-rows dropped, an arm did not converge")
  did <- did %>% filter(both_converged) %>%
    to_log_units(per_sd = TRUE, other_unit = NA_character_)
  if (nrow(did)) pooled$control_did <- pool_by(did, marker, adjustment, unit, panel_h, form) %>%
    mutate(scale = "ventilated minus no-support divergence, log marker per day")
}

pl <- read_family("^jm_pre_placebo_.*\\.csv$")
if (nrow(pl)) {
  pl <- pl %>% filter(status == "fitted") %>%
    mutate(panel_h = jm_panel(file, "jm_pre_placebo")) %>%
    # the pre-intubation rate and the joint model's post-intubation rate, pooled apart
    { bind_rows(transmute(., site, marker, adjustment, panel_h, period = "pre-intubation",
                          estimate = pre_estimate, se = pre_se, rhat = NA_real_),
                transmute(., site, marker, adjustment, panel_h, period = "post-intubation",
                          estimate = post_estimate, se = (post_hi - post_lo) / 3.92, rhat = post_rhat)) }
  if (any(!is.na(pl$rhat) & pl$rhat > RHAT_MAX))
    message("placebo: ", sum(!is.na(pl$rhat) & pl$rhat > RHAT_MAX),
            " post-intubation site-rows dropped for rhat > ", RHAT_MAX)
  pl <- pl %>% filter(is.na(rhat) | rhat <= RHAT_MAX) %>%
    to_log_units(per_sd = TRUE, other_unit = NA_character_)
  if (nrow(pl)) pooled$pre_placebo <- pool_by(pl, marker, adjustment, period, unit, panel_h) %>%
    mutate(scale = "divergence by lung size, log marker per day")
}

st <- read_family("^jm_pre_stacked_.*\\.csv$")
if (nrow(st)) {
  st <- st %>% filter(status == "fitted") %>%
    transmute(site, marker, adjustment, quantity, estimate, se,
              panel_h = jm_panel(file, "jm_pre_stacked")) %>%
    to_log_units(per_sd = TRUE, other_unit = NA_character_)
  if (nrow(st)) pooled$pre_stacked <- pool_by(st, marker, adjustment, quantity, unit, panel_h) %>%
    mutate(scale = "within-patient divergence, log marker per day; change = after - before")
}

# --- write
for (nm in names(pooled)) {
  write_csv(pooled[[nm]], file.path(out_dir, paste0("pooled_biotrauma_", nm, ".csv")))
  message(nm, ": ", nrow(pooled[[nm]]), " pooled rows")
}

# --- figure 4 pooled: the divergence and the three supports that test it, site by
#     site and pooled. Adjusted, main model, the daily panel, PFVC units only: the
#     VT/PFVC divergence is on its own scale and gets its own row of panels.
DIVERGENCE <- c(pfvc = "log_pfvc_sd:vent_day", vtpfvc = canonical_term("vtpfvc_c:vent_day"))
fig4_rows <- function(d, quantity, keep = TRUE) {
  if (is.null(d) || !nrow(d)) return(NULL)
  d %>% filter(keep, adjustment == "adjusted", grepl("d$", panel_h)) %>%
    transmute(marker, unit, quantity, site = anon(site), estimate, lo = estimate - 1.96 * se,
              hi = estimate + 1.96 * se, is_pooled = FALSE)
}
fig4_pooled <- function(d, quantity, keep = TRUE) {
  if (is.null(d) || !nrow(d)) return(NULL)
  d %>% filter(keep, adjustment == "adjusted", grepl("d$", panel_h)) %>%
    transmute(marker, unit, quantity, site = "Pooled", estimate = pooled, lo, hi, is_pooled = TRUE)
}
fd4 <- bind_rows(
  fig4_rows(es, "divergence, ventilated", es$term == DIVERGENCE[["pfvc"]] & es$model == "main"),
  fig4_pooled(pooled$longitudinal_terms, "divergence, ventilated",
              pooled$longitudinal_terms$term == DIVERGENCE[["pfvc"]] & pooled$longitudinal_terms$model == "main"),
  fig4_rows(did, "difference in differences"), fig4_pooled(pooled$control_did, "difference in differences"),
  fig4_rows(pl, "placebo, before intubation", pl$period == "pre-intubation"),
  fig4_pooled(pooled$pre_placebo, "placebo, before intubation", pooled$pre_placebo$period == "pre-intubation"),
  fig4_rows(st, "stacked, change at intubation", st$quantity == "change"),
  fig4_pooled(pooled$pre_stacked, "stacked, change at intubation", pooled$pre_stacked$quantity == "change"))
if (!is.null(fd4) && nrow(fd4)) {
  fd4 <- fd4 %>% filter(startsWith(unit, "per ")) %>%
    mutate(site = factor(site, levels = c(sort(unique(setdiff(site, "Pooled"))), "Pooled")),
           quantity = factor(quantity, c("divergence, ventilated", "difference in differences",
                                         "placebo, before intubation", "stacked, change at intubation")))
  worse <- c(creatinine = "higher", platelets = "lower", bilirubin = "higher", pressor_dose = "higher",
             osi = "higher", any_pressor = "higher", sf = "lower", dp = "higher")
  p4 <- ggplot(fd4, aes(estimate, site, shape = is_pooled)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = lo, xmax = hi), colour = okabe[1]) +
    facet_grid(paste0(marker, "\n(worse = ", worse[marker], ")") ~ quantity, scales = "free_x") +
    scale_shape_manual(values = c(16, 18), guide = "none") +
    labs(title = "Figure 4 pooled: divergence by predicted lung size over 7 days of ventilation",
         subtitle = paste0("log marker per day per ", PER_LOG_PFVC,
                           " log units of PFVC (about a 10% smaller predicted lung); adjusted, common-effect pool"),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10)
  ggsave(file.path(out_dir, "pooled_biotrauma_figure4.pdf"), p4,
         width = 13, height = 2 + 1.1 * n_distinct(fd4$marker) * (1 + n_distinct(fd4$site) / 6), limitsize = FALSE)
}

# --- forests: the PFVC-level contrast per marker and horizon, per site and pooled
if (nrow(lc) && any(lc$form == "pfvc" & lc$panel_h == paste0(lc$horizon_h, "h"))) {
  # the primary read: the pfvc form's contrast at each contrast horizon, from the panel of the same length
  fd <- lc %>% filter(exposure == "log_pfvc_sd", model == "main", form == "pfvc", panel_h == paste0(horizon_h, "h")) %>%
    transmute(marker, adjustment, horizon_h, site = anon(site), estimate, lo, hi, pooled = FALSE) %>%
    bind_rows(pooled$level_contrast %>% filter(exposure == "log_pfvc_sd", model == "main", form == "pfvc", panel_h == paste0(horizon_h, "h")) %>%
                transmute(marker, adjustment, horizon_h, site = "Pooled", estimate = pooled, lo, hi, pooled = TRUE)) %>%
    mutate(site = factor(site, levels = c(sort(unique(setdiff(site, "Pooled"))), "Pooled")),
           adjustment = factor(adjustment, c("adjusted", "unadjusted")))
  # one panel per marker x horizon with its own x scale (facet_grid shares x down a
  # column, so the any-pressor log-odds would set the scale for the labs); the strip
  # names the injury direction so the sign reads without flipping
  worse <- c(creatinine = "higher", platelets = "lower", bilirubin = "higher", sf = "lower", dp = "higher",
             ne_equiv_peak = "higher", any_pressor = "higher")
  fd <- fd %>% mutate(panel = factor(paste0(marker, " (worse = ", worse[marker], ")\n", horizon_h, " h"),
                                     levels = unique(paste0(marker, " (worse = ", worse[marker], ")\n", horizon_h, " h")[order(marker, horizon_h)])))
  p <- ggplot(fd, aes(estimate, site, colour = adjustment, shape = pooled)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5)) +
    facet_wrap(~ panel, scales = "free_x", ncol = n_distinct(fd$horizon_h), dir = "h") +
    scale_colour_manual(values = okabe[1:2]) + scale_shape_manual(values = c(16, 18), guide = "none") +
    labs(title = "Marker difference per SD of log PFVC at the horizon, joint model (death before H modelled)",
         subtitle = "a lower PFVC is the negative of the estimate; log-odds scale for any_pressor",
         x = "per SD of log PFVC", y = NULL) +
    theme_minimal(base_size = 10)
  ggsave(file.path(out_dir, "pooled_biotrauma_level_contrast.pdf"), p, width = 4 + 3.5 * n_distinct(fd$horizon_h),
         height = 2 + 1.2 * n_distinct(fd$marker) * (1 + n_distinct(fd$site) / 6))
}
message("pooled_biotrauma complete -> ", out_dir)
