# =============================================================================
# Pooled biotrauma results across sites (run centrally, like pooled_estimates.R)
# =============================================================================
# Discovers each site's aggregate biotrauma outputs under a results root that
# holds one subfolder per site (each site's final/ renamed to the site name,
# PBWPFVC_RESULTS_ROOT, default results/) and pools them by random-effects
# meta-analysis (metafor::rma, REML, as pooled_estimates.R does). The pooled
# tables and forests go to the All sites/ subfolder, which is excluded from
# discovery.
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

# random-effects pool of one estimate/SE set; returns one row
pool_one <- function(d) {
  d <- d %>% filter(is.finite(estimate), is.finite(se), se > 0)
  k <- nrow(d)
  if (k == 0) return(NULL)
  if (k == 1) return(tibble(k = 1L, pooled = d$estimate, se = d$se, lo = d$estimate - 1.96 * d$se,
                            hi = d$estimate + 1.96 * d$se, i2 = NA_real_, tau2 = NA_real_,
                            sites = d$site, site_estimates = as.character(signif(d$estimate, 4))))
  fit <- tryCatch(rma(yi = d$estimate, sei = d$se, method = "REML"),
                  error = function(e) rma(yi = d$estimate, sei = d$se, method = "DL"))
  tibble(k = k, pooled = as.numeric(fit$b), se = fit$se, lo = fit$ci.lb, hi = fit$ci.ub,
         i2 = fit$I2, tau2 = fit$tau2,
         sites = paste(d$site, collapse = ";"), site_estimates = paste(signif(d$estimate, 4), collapse = ";"))
}
pool_by <- function(d, ...) d %>% group_by(...) %>% group_modify(~ pool_one(.x)) %>% ungroup()
# the modifier form and the panel horizon from a joint-model file name
jm_form  <- function(file, family, default) str_match(file, paste0("^", family, "_(?:(\\w+?)_)?(\\d+[hd])_"))[, 2] %>% replace_na(default)
jm_panel <- function(file, family) str_match(file, paste0("^", family, "_(?:(\\w+?)_)?(\\d+[hd])_"))[, 3]
CHANNELS <- c("ch_height", "ch_age", "ch_sex", "ch_race")

pooled <- list()

# --- 1. joint-model level contrasts (the PFVC-level question)
lc <- read_family("^jm_level_contrast_.*\\.csv$")
if (nrow(lc)) {
  lc <- lc %>% mutate(se = (hi - lo) / 3.92,
                      grid = if ("grid" %in% names(lc)) grid else NA_character_,
                      form = jm_form(file, "jm_level_contrast", "pfvc"), panel_h = jm_panel(file, "jm_level_contrast"))
  pooled$level_contrast <- pool_by(lc, marker, model, adjustment, exposure, panel_h, horizon_h, grid, form) %>%
    mutate(scale = "log marker (log-odds for any_pressor) per SD of the size exposure (per log unit for ch_* pieces)")
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
  key <- c("l_vtpbw_within", "l_vtpbw_within:ldisc_c", "l_vtpbw_within:age10_c", "log_pfvc_sd",
           "vent_day:log_pfvc_sd", "log_pfvc_sd:vent_day", "ldisc_sd", "vent_day:ldisc_sd", "ldisc_sd:vent_day",
           "ers_pfvc_0:l_vtpbw_within", "vtpbw_pt_mean",
           CHANNELS, paste0("vent_day:", CHANNELS), paste0(CHANNELS, ":vent_day"))
  es <- es %>% filter(block == "longitudinal", term %in% key) %>%
    mutate(se = sd, grid = if ("grid" %in% names(es)) grid else NA_character_,
           form = jm_form(file, "jm_estimates", "disc"), panel_h = jm_panel(file, "jm_estimates"))
  pooled$longitudinal_terms <- pool_by(es, marker, model, adjustment, term, panel_h, grid, form)
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

# --- write
for (nm in names(pooled)) {
  write_csv(pooled[[nm]], file.path(out_dir, paste0("pooled_biotrauma_", nm, ".csv")))
  message(nm, ": ", nrow(pooled[[nm]]), " pooled rows")
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
