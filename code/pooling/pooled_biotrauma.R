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
#   jm_control_did_*       figure 4's difference-in-differences: the ventilated
#                          divergence minus the no-support control's, per day
#   jm_hypoxemic_control_did_*  the same against the hypoxemic control (SF <= 315)
#   pfvc_age_control_channels_*, _channel_vcov_*, _contrast_*   the mortality control
#                          contrast and its GLI channel breakdown (supplement/)
#   crs_channels_estimates_*, crs_channels_tests_*   the compliance channels and the
#                          compliance head-to-head (supplement/); the channel results
#                          are drawn in All sites/pooled_biotrauma_channels.pdf
#   fingerprint_*, fingerprint_did_*   the height fingerprint (28): the rate per log
#                          unit of PBW/PFVC moved by height within sex, the whole
#                          ratio's rate it is read against, and the DiD
# Joint-model tables are keyed by the panel horizon in the file tag (24h/48h/72h)
# as well as the contrast horizon, so one site's three panels are never pooled as
# three sites. Every pooled row carries k, I2, tau2 and the per-site estimates it
# was built from. Sites are anonymised with utils/site_anonymization.R when present.
#
# Usage: PBWPFVC_RESULTS_ROOT=/path/to/results uvr run code/pooling/pooled_biotrauma.R
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
                   sf = "SpO2 / FiO2\n(positive control)", dp = "Driving pressure")
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
read_family <- function(pattern, folders = c("", "injury")) {
  map_dfr(sites, function(s) {
    # a returned final/ is sorted by block; most tables pooled here are in injury/, the
    # supplement's in supplement/ (folders = c("", "supplement")). controls/ is never
    # listed, so a control cohort cannot enter a pool of the ventilated one. The site
    # folder itself is listed too, for a flat (older) return.
    fs <- list.files(file.path(root, s, folders), pattern = pattern, full.names = TRUE)
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
  # the panel tag comes from the file name (jm_scale_48h_*, jm_scale_7d_*), as for every
  # other family: horizon_days would turn a 48-hour panel into "2d" and miss its estimates
  transmute(site, panel_h = str_match(file, "^jm_scale_(\\d+[hd])_")[, 2], sd_log_pfvc) %>% distinct() else
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

# --- 7b. the same difference against the hypoxemic control (index-day SF <= 315,
#         27_control_comparison.R, jm_hypoxemic_control_did_*): the arms then differ in
#         ventilation and not in hypoxemia. Its own file family, so it never mixes with
#         figure 4's DiD; not drawn in the pooled figure 4.
hdid <- read_family("^jm_hypoxemic_control_did_.*\\.csv$")
if (nrow(hdid)) {
  hdid <- hdid %>% transmute(site, marker, adjustment, form = jm_form(file, "jm_hypoxemic_control_did", "pfvc"),
                             panel_h = jm_panel(file, "jm_hypoxemic_control_did"),
                             estimate = did_estimate, se = did_sd, both_converged)
  if (any(!hdid$both_converged))
    message("hypoxemic-control difference in differences dropped, an arm did not converge:\n  ",
            hdid %>% filter(!both_converged) %>% transmute(what = paste(site, marker, adjustment)) %>%
              pull(what) %>% paste(collapse = "\n  "))
  hdid <- hdid %>% filter(both_converged) %>% to_log_units(per_sd = TRUE, other_unit = NA_character_)
  if (nrow(hdid)) pooled$hypoxemic_control_did <- pool_by(hdid, marker, adjustment, unit, panel_h, form) %>%
    mutate(scale = "ventilated minus hypoxemic no-support divergence, log marker per day")
}

# --- 8. the height fingerprint (28_height_fingerprint.R): the rate per log unit of
#        PBW/PFVC moved by height within sex, beside the whole ratio's rate (the
#        predicted value), and the difference in differences against the no-support
#        control. Log units of the ratio are the same at every site, so nothing is
#        converted. The pooled minimum detectable effect (80% power) says whether
#        the pool can see the predicted value.
fp <- read_family("^fingerprint_.*\\.csv$")
if (nrow(fp)) {
  fp_est <- fp %>% filter(grepl("^fingerprint_[a-z]+_", file), !grepl("^fingerprint_(ladder|curves|did)_", file))
  if (nrow(fp_est)) {
    predicted <- fp_est %>% filter(model == "whole ratio", quantity == "rate") %>%
      pool_by(marker, adjustment) %>% select(marker, adjustment, predicted_rate = pooled)
    pooled$fingerprint <- pool_by(fp_est, marker, model, quantity, shared_df, adjustment) %>%
      left_join(predicted, by = c("marker", "adjustment")) %>%
      mutate(predicted_rate = if_else(model == "fingerprint" & quantity == "rate", predicted_rate, NA_real_),
             mde_80 = 2.80 * se, detectable = abs(predicted_rate) >= mde_80,
             scale = "log marker (per day for the rate) per log unit of PBW/PFVC")
  }
  fp_did <- fp %>% filter(grepl("^fingerprint_did_", file))
  if (nrow(fp_did))
    pooled$fingerprint_did <- fp_did %>% transmute(site, marker, shared_df, adjustment, estimate = did_estimate, se = did_se) %>%
      pool_by(marker, shared_df, adjustment) %>%
      mutate(scale = "ventilated minus no-support fingerprint rate, log marker per day per log unit of PBW/PFVC")
}

# --- 9. the supplement's contrasts (2026-09-24): the mortality control contrast and its
#        channel breakdown (supplement/xsec_pfvc_age_control.R), and the compliance
#        channels (supplement/xsec_crs_channels.R). Read from each site's supplement/.
#   channels     each GLI piece per cohort and ventilated minus no support, per log
#                unit of the piece: the same unit at every site, pooled as it stands.
#                Each site reports the pieces on two exposure scales (log PFVC, and
#                log PBW/PFVC, the strain error), so every pool and test keys on exposure.
#                The pooled test that the pieces' differences agree needs each site's
#                covariance between pieces (pfvc_age_control_channel_vcov_*), pooled by
#                multivariate common-effect inverse variance; a site without that table
#                enters the per-piece pools but not the test, and the test row names the
#                sites it used.
#   contrast     per SD of log PFVC in each site's own ventilated cohort; converted to
#                per 0.1 log units with the site's exported SD, and kept per site SD,
#                labelled and never pooled with converted rows, where the SD is missing
#   compliance   Crs exponents (unitless) pooled as they stand; the head-to-head AIC
#                differences summed across sites (AIC is additive over independent
#                samples), with the count of sites and how many favoured each side
SUPPLEMENT <- c("", "supplement")
equal_test_p <- function(b, V, contrast_matrix) {
  d <- contrast_matrix %*% b
  as.numeric(pchisq(t(d) %*% solve(contrast_matrix %*% V %*% t(contrast_matrix)) %*% d,
                    df = nrow(contrast_matrix), lower.tail = FALSE))
}
PIECES <- c("height", "age", "sex", "race")

channels_tbl <- read_family("^pfvc_age_control_channels_.*\\.csv$", SUPPLEMENT)
if (nrow(channels_tbl)) {
  pooled$age_control_channels <- channels_tbl %>% filter(!is.na(log_ratio), !is.na(se)) %>%
    transmute(site, exposure, population, outcome, ratio_type, quantity, piece, estimate = log_ratio, se) %>%
    pool_by(exposure, population, outcome, ratio_type, quantity, piece) %>%
    mutate(ratio_per_0.1 = exp(0.1 * pooled), lo_per_0.1 = exp(0.1 * lo), hi_per_0.1 = exp(0.1 * hi),
           scale = paste0("log ratio per log unit of the GLI piece, ", exposure, " scale; ratio_per_0.1 per 0.1 log units"))
  vcov_tbl <- read_family("^pfvc_age_control_channel_vcov_.*\\.csv$", SUPPLEMENT)
  if (nrow(vcov_tbl)) {
    differences <- channels_tbl %>% filter(quantity == "ventilated minus no support", piece %in% PIECES)
    pooled$age_control_channel_tests <- vcov_tbl %>% distinct(exposure, population, outcome) %>%
      pmap_dfr(function(exposure, population, outcome) {
        per_site <- map(unique(vcov_tbl$site), function(s) {
          b <- differences %>% filter(site == s, .data$exposure == .env$exposure, .data$population == .env$population,
                                      .data$outcome == .env$outcome)
          V <- vcov_tbl %>% filter(site == s, .data$exposure == .env$exposure, .data$population == .env$population,
                                   .data$outcome == .env$outcome)
          if (nrow(b) != 4 || nrow(V) != 16 || anyNA(b$log_ratio)) return(NULL)
          b_vec <- setNames(b$log_ratio, b$piece)[PIECES]
          V_mat <- matrix(NA_real_, 4, 4, dimnames = list(PIECES, PIECES))
          V_mat[cbind(V$piece_row, V$piece_col)] <- V$covariance
          list(site = s, b = b_vec, W = solve(V_mat))
        }) %>% compact()
        if (!length(per_site)) return(NULL)
        W_sum <- Reduce(`+`, map(per_site, "W"))
        V_pooled <- solve(W_sum)
        b_pooled <- as.numeric(V_pooled %*% Reduce(`+`, map(per_site, ~ .x$W %*% .x$b)))
        names(b_pooled) <- PIECES
        four_equal <- rbind(c(1, -1, 0, 0), c(1, 0, -1, 0), c(1, 0, 0, -1))
        size_equal <- rbind(c(1, 0, -1, 0), c(1, 0, 0, -1))    # height = sex = race, age left out
        tibble(exposure = exposure, population = population, outcome = outcome, k = length(per_site),
               sites = paste(map_chr(per_site, "site"), collapse = ";"),
               test = c("the four differences are equal (3 df)", "height = sex = race differences (2 df)"),
               p = c(equal_test_p(b_pooled, V_pooled, four_equal), equal_test_p(b_pooled, V_pooled, size_equal)))
      })
  }
}

contrast_tbl <- read_family("^pfvc_age_control_contrast_.*\\.csv$", SUPPLEMENT)
if (nrow(contrast_tbl)) {
  if (!"ventilated_log_pfvc_sd" %in% names(contrast_tbl)) contrast_tbl$ventilated_log_pfvc_sd <- NA_real_
  pooled$age_control_contrast <- contrast_tbl %>% filter(!is.na(log_ratio), !is.na(se)) %>%
    mutate(converted = !is.na(ventilated_log_pfvc_sd) & !grepl("x anchor", quantity),
           scale_factor = if_else(converted, PER_LOG_PFVC / ventilated_log_pfvc_sd, 1),
           unit = case_when(grepl("x anchor", quantity) ~ "PFVC x anchor term, per site SD (not harmonised)",
                            converted ~ paste0("per ", PER_LOG_PFVC, " log PFVC"),
                            TRUE ~ "per site SD of log PFVC (not harmonised)"),
           estimate = log_ratio * scale_factor, se = se * scale_factor) %>%
    select(site, population, outcome, ratio_type, adjustment, severity, quantity, unit, estimate, se) %>%
    pool_by(population, outcome, ratio_type, adjustment, severity, quantity, unit) %>%
    mutate(ratio = exp(pooled), ratio_lo = exp(lo), ratio_hi = exp(hi))
}

# each site labels its short-women subgroup with its own median female height ("women
# shorter than 163 cm"); the subgroup is the same definition everywhere, so the label is
# made common before pooling, or each site would pool alone
common_short_women <- function(label) sub("women shorter than [0-9.]+ cm", "women below the median female height", label)
crs_tbl <- read_family("^crs_channels_estimates_.*\\.csv$", SUPPLEMENT)
if (nrow(crs_tbl)) {
  pooled$crs_channels <- crs_tbl %>% filter(!is.na(estimate), !is.na(se)) %>%
    mutate(model = common_short_women(model)) %>%
    select(site, sample, model, term, estimate, se) %>%
    pool_by(sample, model, term) %>%
    mutate(scale = "exponent of log Crs (1 = proportional scaling)")
  crs_tests <- read_family("^crs_channels_tests_.*\\.csv$", SUPPLEMENT)
  if (nrow(crs_tests)) pooled$crs_channel_aic <- crs_tests %>% filter(grepl("AIC", test)) %>%
    mutate(test = common_short_women(test)) %>%
    group_by(sample, test) %>%
    summarise(k = n(), summed_delta_aic = sum(statistic), sites_below_zero = sum(statistic < 0),
              sites = paste(site, collapse = ";"), site_delta_aic = paste(round(statistic, 1), collapse = ";"),
              .groups = "drop") %>%
    mutate(note = "AIC differences summed across sites (additive over independent samples); below 0 favours the non-PBW exposure")
}

# --- 10. figures for the channel results (2026-09-24): one PDF, three pages.
#   page 1  the candidate figure 5: for each GLI piece, the ventilated-minus-no-support
#           mortality difference by site and pooled, for the three outcomes, beside the
#           Crs exponent through the same piece (the inputs that move measured compliance
#           are the ones strain predicts will carry a ventilator-specific association)
#   page 2  the same differences across populations (everyone, full code, hypoxemic),
#           in-hospital death
#   page 3  each cohort's own piece coefficients, pooled: where each difference comes from
# Ratios per 0.1 log units of the piece (about 10% of PFVC); below 1 = a larger predicted
# lung through that input goes with less death. Pooled rows are diamonds. The pages are
# drawn on the log PFVC scale; the strain-error scale is in the tables only.
FIGURE_EXPOSURE <- "log PFVC"
PIECE_ORDER <- c("height", "sex", "race", "age", "all four (one beta)")
PIECE_COLOURS <- setNames(okabe[c(1, 3, 5, 2, 8)], PIECE_ORDER)
DIFFERENCE <- "ventilated minus no support"
if (exists("channels_tbl") && nrow(channels_tbl) && !is.null(pooled$age_control_channels)) {
  per_0.1 <- function(d) d %>% mutate(ratio = exp(0.1 * log_ratio), lo = exp(0.1 * (log_ratio - 1.96 * se)),
                                      hi = exp(0.1 * (log_ratio + 1.96 * se)))
  channel_rows <- function(quantities, populations, outcomes) {
    sites_part <- channels_tbl %>%
      filter(exposure == FIGURE_EXPOSURE, quantity %in% quantities, population %in% populations, outcome %in% outcomes, !is.na(log_ratio)) %>%
      per_0.1() %>% transmute(population, outcome, quantity, piece, site = anon(site), ratio, lo, hi, is_pooled = FALSE)
    pooled_part <- pooled$age_control_channels %>%
      filter(exposure == FIGURE_EXPOSURE, quantity %in% quantities, population %in% populations, outcome %in% outcomes) %>%
      transmute(population, outcome, quantity, piece, site = "Pooled", ratio = ratio_per_0.1, lo = lo_per_0.1,
                hi = hi_per_0.1, is_pooled = TRUE)
    bind_rows(sites_part, pooled_part) %>%
      mutate(piece = factor(piece, levels = PIECE_ORDER),
             site = factor(site, levels = c(sort(unique(setdiff(site, "Pooled"))), "Pooled")))
  }
  channel_theme <- theme_minimal(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
  forest <- function(d, facet_formula, title, subtitle) {
    ggplot(d, aes(ratio, fct_rev(site), colour = piece, shape = is_pooled)) +
      geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
      geom_pointrange(aes(xmin = lo, xmax = hi, size = is_pooled)) +
      scale_size_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.6), guide = "none") +
      facet_grid(facet_formula, switch = "y") + scale_x_log10() +
      scale_colour_manual(values = PIECE_COLOURS, guide = "none") +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
      labs(title = title, subtitle = subtitle, x = "ratio per 0.1 log units of the piece (log scale)", y = NULL) +
      channel_theme + theme(strip.text.y.left = element_text(angle = 0), strip.placement = "outside")
  }
  outcomes_main <- c("in-hospital death (logistic)", "60-day death, all", "60-day death, before invasive ventilation")
  page_1a <- forest(channel_rows(DIFFERENCE, "everyone", outcomes_main) %>% mutate(outcome = factor(outcome, outcomes_main)),
                    piece ~ outcome, "A. Ventilated minus no support, by GLI input",
                    "mortality OR (in-hospital) or cause-specific HR; below 1 = more protective under ventilation than without it")
  # the Crs exponent through each piece, by site and pooled (supplement/xsec_crs_channels.R)
  page_1b <- NULL
  if (exists("crs_tbl") && nrow(crs_tbl) && !is.null(pooled$crs_channels)) {
    crs_rows <- bind_rows(
      crs_tbl %>% filter(sample == "all plateau-measured", model == "channels") %>%
        transmute(piece = term, site = anon(site), estimate, lo = estimate - 1.96 * se, hi = estimate + 1.96 * se, is_pooled = FALSE),
      pooled$crs_channels %>% filter(sample == "all plateau-measured", model == "channels") %>%
        transmute(piece = term, site = "Pooled", estimate = pooled, lo, hi, is_pooled = TRUE)) %>%
      mutate(piece = factor(piece, levels = PIECE_ORDER),
             site = factor(site, levels = c(sort(unique(setdiff(site, "Pooled"))), "Pooled")))
    page_1b <- ggplot(crs_rows, aes(estimate, fct_rev(site), colour = piece, shape = is_pooled)) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
      geom_vline(xintercept = 1, colour = "grey30") +
      geom_pointrange(aes(xmin = lo, xmax = hi, size = is_pooled)) +
      scale_size_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.6), guide = "none") +
      facet_grid(piece ~ ., switch = "y") +
      scale_colour_manual(values = PIECE_COLOURS, guide = "none") +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
      labs(title = "B. Crs exponent, same input", subtitle = "1 = proportional to PFVC",
           x = "d log Crs / d log PFVC (piece)", y = NULL) +
      channel_theme + theme(strip.text.y.left = element_text(angle = 0), strip.placement = "outside")
  }
  page_1 <- if (is.null(page_1b)) page_1a else patchwork::wrap_plots(page_1a, page_1b, widths = c(3, 1))
  populations_shown <- intersect(c("everyone", "full code at the index", "hypoxemic at the index (SF < 315)"),
                                 unique(channels_tbl$population))
  page_2 <- forest(channel_rows(DIFFERENCE, populations_shown, "in-hospital death (logistic)") %>%
                     mutate(population = factor(population, populations_shown)),
                   piece ~ population, "Ventilated minus no support, by GLI input and population",
                   "in-hospital death, OR per 0.1 log units of the piece; hypoxemic = index SF < 315 in both cohorts")
  cohort_rows <- pooled$age_control_channels %>%
    filter(exposure == FIGURE_EXPOSURE, population == "everyone", outcome %in% outcomes_main, quantity %in% c("Ventilated", "No support")) %>%
    transmute(outcome = factor(outcome, outcomes_main), quantity, piece = factor(piece, levels = PIECE_ORDER),
              ratio = ratio_per_0.1, lo = lo_per_0.1, hi = hi_per_0.1)
  page_3 <- ggplot(cohort_rows, aes(ratio, fct_rev(piece), colour = quantity)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5), size = 0.35) +
    facet_wrap(~ outcome) + scale_x_log10() +
    scale_colour_manual(values = c(Ventilated = okabe[[1]], `No support` = okabe[[2]]), breaks = c("Ventilated", "No support"), name = NULL) +
    labs(title = "Each cohort's own association through each GLI input (pooled)",
         subtitle = "everyone; the difference in figure A is the gap between the two points of each input",
         x = "ratio per 0.1 log units of the piece (log scale)", y = NULL) +
    channel_theme
  pdf(file.path(out_dir, "pooled_biotrauma_channels.pdf"), width = 13, height = 8.5)
  print(page_1); print(page_2); print(page_3)
  invisible(dev.off())
  message("channel figures -> ", file.path(out_dir, "pooled_biotrauma_channels.pdf"))
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
