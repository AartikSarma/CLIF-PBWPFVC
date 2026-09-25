# =============================================================================
# Script 23 (report): Biotrauma joint models -- contrasts, movement, associations
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# Reads the fit bundles 22_biotrauma_fit.R saved and writes the aggregate,
# poolable tables, all tagged {tag} as the fit's own files:
#
#   final/jm_level_contrast_{tag}.csv   figure 4A: the marker difference per SD of
#       log PFVC (or per log unit of each GLI piece) at each horizon, level plus
#       divergence x time, from the joint posterior draws
#   final/jm_movement_{tag}.csv         the observed mean change from baseline by
#       day (descriptive), read by 27 before a control's divergence
#   final/jm_lme_check_{tag}.csv        the level and divergence terms from the joint
#       model beside the same terms from the longitudinal submodel fitted alone (no
#       death or extubation correction): whether the unconverged hazard and
#       association blocks move the answer
#   final/jm_association_hr_{tag}.csv   hazard ratio for death, extubation (and RRT)
#       per SD of the current log marker (value) and per unit slope
#
# The within-patient dose terms are not reported: the manuscript does not read the
# dose inside the band, which is confounded by indication.
#
# Switches (environment), each as set for the fit being reported:
#   PBWPFVC_JM_MODIFIER    size form (default pfvc)
#   PBWPFVC_JM_RRT_EVENT   1 = the creatinine fit with dialysis as a third competing
#                          cause, written as rrtcause_ tables (default 0)
#   PBWPFVC_JM_MARKERS     comma list; report only these markers and merge their rows
#                          into the existing tables (default: every marker fitted)
#   PBWPFVC_JM_SF_BAND     "lo,hi": the baseline SF class (default none)
#   PBWPFVC_JM_SEV_CENTER  a control read at the ventilated severity (default none)
#   PBWPFVC_JM_GRID, PBWPFVC_JM_HORIZON   grid and window (default daily, 7 days)
#   PBWPFVC_JM_BASELINE    baseline form (default free)
# The last five are read through 20_biotrauma_grid.R, which builds the tag from them.
# Figure 4 = PBWPFVC_JM_MODIFIER=pfvc, daily grid, 7 days, as set by 29_run_figure4.R.
#
# Usage: uvr run code/23_biotrauma_report.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(splines)
  library(nlme)
  library(survival)
  library(JMbayes2)
  library(GLMMadaptive)   # fixef() for the logistic mixed model of the any-pressor fit
})
rm(list = ls())
source("utils/config.R")

site_name  <- config$site_name
output_dir <- config$output_dir
final_dir  <- final_dir_for("injury")
source(here("code", "20_biotrauma_grid.R"))   # JM_GRID, STEP, JM_HORIZON, N_PERIODS, h_suffix
BASELINE_FORM <- Sys.getenv("PBWPFVC_JM_BASELINE", "free")
MOD_FORM      <- Sys.getenv("PBWPFVC_JM_MODIFIER", "pfvc")   # figure 4
# RRT as a third competing cause (creatinine only): its own tables and its own
# bundles, so the two-cause primary is never overwritten by the sensitivity
RRT_EVENT     <- identical(Sys.getenv("PBWPFVC_JM_RRT_EVENT", "0"), "1")
stopifnot(MOD_FORM %in% c("disc", "saturated", "none", "pfvc", "disc_level", "channels", "vtpfvc", "pfvc_dose"))
out_tag  <- paste0(if (RRT_EVENT) "rrtcause_" else "", restrict_tag,
                   if (BASELINE_FORM == "offset") "offset_" else "",
                   if (MOD_FORM != "disc") paste0(MOD_FORM, "_") else "",
                   h_suffix, "_", site_name)
N_DRAWS <- 1000L
set.seed(20260913)

manifest <- read_csv(file.path(final_dir, paste0("jm_manifest_", out_tag, ".csv")), show_col_types = FALSE)
# the gate by block, from the estimates table (so every fit is gated the same
# way, whatever its manifest carries): longitudinal for the trajectory
# contrasts and Q1, association for the death correction and Q2, hazard for Q3
RHAT_GATE <- 1.1   # the standard convergence threshold
# why each marker has no fit, from the manifest (skipped for too few patients or
# deaths, or failed), for the messages below
why_not <- function(m) {
  r <- manifest %>% filter(marker %in% m) %>% distinct(marker, status, reason)
  paste(paste0(r$marker, " ", r$status, ifelse(is.na(r$reason), "", paste0(" (", r$reason, ")"))), collapse = "; ")
}
# An arm where no fit ran (a small SF class, a thin control) is a skip, not an error:
# say so and stop cleanly, so a runner over many arms carries on.
if (!any(manifest$status %in% c("converged", "rhat_fail"))) {
  message("nothing to report for ", out_tag, ": no fit ran: ", why_not(unique(manifest$marker)))
  quit(save = "no", status = 0)
}
est_tbl <- read_csv(file.path(final_dir, paste0("jm_estimates_", out_tag, ".csv")), show_col_types = FALSE)
# The exposure gate: the longitudinal block also carries nuisance terms (the
# baseline marker, the age spline) whose chains mix worse than the exposure
# terms; the trajectory contrasts need only the exposure terms and their
# time interactions, so those are gated separately (`exposure_rhat`).
SIZE_EXPOS <- c("log_pfvc_sd", "ldisc_sd", "vtpfvc_c", CHANNELS)
DOSE_MOD   <- c("log_pfvc_sd:vtpbw_c", "vtpbw_c:log_pfvc_sd",
                "log_pfvc_sd:vtpbw_c:vent_day", "vent_day:log_pfvc_sd:vtpbw_c",
                "log_pfvc_sd:vent_day:vtpbw_c", "vtpbw_c:vent_day")
EXPO_TERMS <- c("l_vtpbw_within", "vtpbw_pt_mean", SIZE_EXPOS, paste0(SIZE_EXPOS, ":vent_day"),
                paste0("vent_day:", SIZE_EXPOS), DOSE_MOD)
block_gates <- est_tbl %>% group_by(marker, model, adjustment) %>%
  summarise(longitudinal_rhat = suppressWarnings(max(rhat[block == "longitudinal"], na.rm = TRUE)),
            exposure_rhat     = suppressWarnings(max(rhat[block == "longitudinal" & term %in% EXPO_TERMS], na.rm = TRUE)),
            association_rhat  = suppressWarnings(max(rhat[block == "association"],  na.rm = TRUE)),
            hazard_rhat       = suppressWarnings(max(rhat[block == "survival"],     na.rm = TRUE)), .groups = "drop") %>%
  mutate(across(ends_with("_rhat"), ~ if_else(is.finite(.), ., NA_real_)))
manifest <- manifest %>% select(-any_of(c("longitudinal_rhat", "exposure_rhat", "association_rhat", "hazard_rhat"))) %>%
  left_join(block_gates, by = c("marker", "model", "adjustment"))
usable <- manifest %>% filter(status %in% c("converged", "rhat_fail"))
# PBWPFVC_JM_MARKERS restricts the report to those markers, as it does the fit:
# re-reading every bundle costs minutes per marker, and a run that added one
# marker should not have to redo the rest. Each table is then merged on write,
# replacing only its own marker rows (`report_write` below). The figures are
# drawn from the tables afterwards, so they still show every marker.
want_markers <- trimws(strsplit(Sys.getenv("PBWPFVC_JM_MARKERS", ""), ",")[[1]])
if (length(want_markers) && nzchar(want_markers[1])) {
  # a requested marker without a usable fit (skipped or failed in this arm) is named
  # and left out; the others are reported
  missing <- setdiff(want_markers, unique(usable$marker))
  if (length(missing)) message("no usable fit in ", out_tag, " for: ", why_not(missing),
                               if (!length(intersect(missing, manifest$marker))) " (not in the manifest)" else "")
  usable <- usable %>% filter(marker %in% want_markers)
  if (!nrow(usable)) {
    message("nothing to report for ", out_tag, " among ", paste(want_markers, collapse = ", "))
    quit(save = "no", status = 0)
  }
  message("restricted to markers: ", paste(want_markers, collapse = ", "))
}
RESTRICTED <- length(want_markers) && nzchar(want_markers[1])
# merge on write: keep the rows of markers this run did not refit
report_write <- function(new, name) {
  path <- file.path(final_dir, paste0("jm_", name, "_", out_tag, ".csv"))
  if (RESTRICTED && file.exists(path) && nrow(new)) {
    old <- read_csv(path, show_col_types = FALSE) %>% filter(!marker %in% want_markers)
    as_text <- function(d) d %>% mutate(across(everything(), as.character))
    if (nrow(old)) {
      message("  ", name, ": kept ", nrow(old), " rows from other markers")
      new <- bind_rows(as_text(old), as_text(new))
    }
  }
  if (nrow(new)) write_csv(new, path)
}
if (nrow(usable) == 0L) stop("No fitted joint models in the manifest for ", out_tag)
message("=== 23_biotrauma_report (", out_tag, "): ", nrow(usable), " fits, of which ",
        sum(usable$status == "converged"), " pass the R-hat gate ===")

# --- posterior draws of the longitudinal fixed effects, stacked across chains
beta_draws <- function(jm, lme_fit) {
  b <- do.call(rbind, jm$mcmc$betas1)
  ref <- names(fixef(lme_fit))
  if (is.null(colnames(b))) colnames(b) <- ref
  stopifnot(all(ref %in% colnames(b)))
  b[, ref, drop = FALSE]
}
assoc_rows <- list()
level_rows <- list(); movement_rows <- list(); lme_check_rows <- list()
# Horizons for the level contrast: one per day out to the run's own endpoint,
# plus the endpoint itself when the horizon is not a whole number of days. Even
# spacing, because the contrast is a level plus a rate times time and a reader
# comparing rows is reading a slope off the page.
# t = 0 is included so the trend panel starts at the index rather than at day 1.
# The contrast there is the level term alone, which is an extrapolation: no
# marker row exists at day 0 (each trajectory starts the period after that
# patient's baseline draw), so read it as the model's anchor, not as data.
LEVEL_HOURS <- sort(unique(c(0, seq(24, floor(JM_HORIZON) * 24, by = 24), JM_HORIZON * 24)))

for (i in seq_len(nrow(usable))) {
  u <- usable[i, ]
  tag <- paste(u$marker, u$model, u$adjustment, sep = "_")
  f <- file.path(output_dir, paste0("jm_fit_", tag, "_", BASELINE_FORM,
                                    if (MOD_FORM != "disc") paste0("_", MOD_FORM) else "",
                                    if (RRT_EVENT && u$marker == "creatinine") "_rrtcause" else "",
                                    restrict_sfx_for(u$marker), "_", h_suffix, ".rds"))
  if (!file.exists(f)) stop("fit bundle missing: ", f)
  b <- readRDS(f); jm <- b$jm; ld <- b$long_data
  draws <- beta_draws(jm, b$lme)
  keep <- sample.int(nrow(draws), min(N_DRAWS, nrow(draws)))
  draws <- draws[keep, , drop = FALSE]
  binary <- isTRUE(b$binary)
  # How much the marker moves at all, by day: the observed change from baseline among
  # patients still observed (descriptive, survivor-selected, not a model quantity). A
  # control cohort whose marker does not move cannot show a divergence by lung size,
  # so this is read BEFORE its divergence term.
  # (the dose part has a change from baseline only for patients on a pressor at day 0)
  movement_rows[[length(movement_rows) + 1L]] <- ld %>%
    filter(if ("on_y0" %in% names(ld)) on_y0 == 1 else TRUE) %>%
    mutate(change = if (binary) NA_real_ else if (BASELINE_FORM == "offset") log_y else log_y - log_y0,
           day = floor(vent_day + 1e-9)) %>%
    filter(day >= 1) %>%
    group_by(hospitalization_id, day) %>%
    summarise(change = mean(change), level = mean(log_y), .groups = "drop") %>%
    group_by(day) %>%
    summarise(n_patients = n(), mean_change = mean(change), sd_change = sd(change),
              mean_abs_change = mean(abs(change)), mean_level = mean(level), .groups = "drop") %>%
    mutate(marker = u$marker, model = u$model, adjustment = u$adjustment, binary = binary, .before = 1)
  sd_log_y <- if (binary) 1 else sd(ld$log_y)   # binary outcome: report on the log-odds scale
  gate <- u$status == "converged"
  gate_long  <- isTRUE(is.finite(u$longitudinal_rhat) && u$longitudinal_rhat <= RHAT_GATE)
  gate_assoc <- isTRUE(is.finite(u$association_rhat)  && u$association_rhat  <= RHAT_GATE)
  gate_expo  <- isTRUE(is.finite(u$exposure_rhat)     && u$exposure_rhat     <= RHAT_GATE)
  message(sprintf("  %-40s longitudinal %s (exposure terms %s), association %s, hazard %s", tag,
                  if (gate_long) "pass" else sprintf("FAIL (%.2f)", u$longitudinal_rhat),
                  if (gate_expo) "pass" else sprintf("FAIL (%.2f)", u$exposure_rhat),
                  if (gate_assoc) "pass" else sprintf("FAIL (%.2f)", u$association_rhat),
                  if (isTRUE(is.finite(u$hazard_rhat) && u$hazard_rhat <= RHAT_GATE)) "pass" else sprintf("FAIL (%.2f)", u$hazard_rhat)))

  # ---- PFVC-level question: marker difference per SD of log PFVC (or of log
  #      PBW/PFVC) at each horizon hour within the grid, level + divergence x time,
  #      from the joint posterior (death before H handled by the shared random effects)
  #      Channels form: one contrast per GLI piece (per log unit of the piece), and
  #      a Wald test on the posterior mean and covariance of the four contrasts
  #      that they are equal (equal = lung size is the operative quantity).
  contrast_draws <- function(ex, hh) {
    tcol <- intersect(c(paste0(ex, ":vent_day"), paste0("vent_day:", ex)), colnames(draws))   # R orders the pair by appearance
    draws[, ex] + (if (length(tcol)) draws[, tcol[1]] else 0) * hh / 24
  }
  for (hh in LEVEL_HOURS[LEVEL_HOURS <= JM_HORIZON * 24]) {
    exs <- intersect(c("log_pfvc_sd", "ldisc_sd", "vtpfvc_c", CHANNELS), colnames(draws))
    if (!length(exs)) next
    V <- sapply(exs, contrast_draws, hh = hh)                       # draws x exposures
    p_equal <- if (all(CHANNELS %in% exs)) channels_equal_p(colMeans(V[, CHANNELS]), cov(V[, CHANNELS])) else NA_real_
    for (ex in exs) {
      v <- V[, ex]
      level_rows[[length(level_rows) + 1L]] <- tibble(
        marker = u$marker, model = u$model, adjustment = u$adjustment, exposure = ex, horizon_h = hh,
        estimate = mean(v), lo = quantile(v, 0.025), hi = quantile(v, 0.975),
        p_gt0 = mean(v > 0), per_sd_marker = if (binary) NA_real_ else mean(v) / sd_log_y,
        unit = if (ex %in% CHANNELS) "per log unit of the piece" else if (ex == "vtpfvc_c") "per point of VT/PFVC (% of predicted FVC)" else "per SD of the exposure",
        p_equal = if (ex %in% CHANNELS) p_equal else NA_real_,
        scale = if (binary) "log-odds of any pressor" else "log marker",
        n_patients = u$n_patients, n_deaths = u$n_deaths,
        rhat_gate = gate_long, rhat_gate_exposure = gate_expo, longitudinal_rhat = u$longitudinal_rhat,
        exposure_rhat = u$exposure_rhat, association_rhat = u$association_rhat)
    }
  }

  # ---- does the death correction move the answer? The fit bundle carries the
  #      longitudinal submodel as fitted alone (nlme, before the joint model), so the
  #      size terms can be read with and without the linkage to death and extubation.
  #      If the hazard or association block fails R-hat but the joint model's level
  #      and divergence agree with the LME's, the unconverged blocks are not what the
  #      estimate rests on. Continuous markers only
  #      (the any-pressor part is a GLMMadaptive fit on another scale).
  if (!binary && inherits(b$lme, "lme")) {
    fe <- nlme::fixef(b$lme); fe_v <- vcov(b$lme)
    for (ex in intersect(c("log_pfvc_sd", "ldisc_sd", "vtpfvc_c", CHANNELS), colnames(draws))) {
      rate <- intersect(c(paste0(ex, ":vent_day"), paste0("vent_day:", ex)), colnames(draws))
      for (tm in c(ex, rate)) {
        if (!tm %in% names(fe)) next
        v <- draws[, tm]
        lme_check_rows[[length(lme_check_rows) + 1L]] <- tibble(
          marker = u$marker, model = u$model, adjustment = u$adjustment, exposure = ex,
          term = if (tm == ex) "level" else "divergence per day",
          jm_estimate = mean(v), jm_lo = quantile(v, 0.025), jm_hi = quantile(v, 0.975),
          lme_estimate = unname(fe[[tm]]), lme_se = sqrt(fe_v[tm, tm]),
          lme_lo = lme_estimate - 1.96 * lme_se, lme_hi = lme_estimate + 1.96 * lme_se,
          jm_minus_lme_in_lme_se = (jm_estimate - lme_estimate) / lme_se,
          exposure_rhat = u$exposure_rhat, association_rhat = u$association_rhat, hazard_rhat = u$hazard_rhat,
          n_patients = u$n_patients, n_deaths = u$n_deaths)
      }
    }
  }

  # ---- Q2 association: HR per SD of the current log marker (value) and per unit slope
  al <- do.call(rbind, jm$mcmc$alphas)[keep, , drop = FALSE]
  for (cn in colnames(al)) {
    kind  <- if (grepl("value", cn)) "value" else "slope"
    # a third cause (RRT, creatinine only) must not be silently labelled extubation
    cause <- if (grepl("death", cn)) "death" else if (grepl("rrt", cn)) "rrt" else "extubation"
    scale <- if (kind == "value") sd_log_y else 1
    v <- al[, cn] * scale
    assoc_rows[[length(assoc_rows) + 1L]] <- tibble(
      marker = u$marker, model = u$model, adjustment = u$adjustment, kind = kind, cause = cause,
      log_hr = mean(v), log_hr_lo = quantile(v, 0.025), log_hr_hi = quantile(v, 0.975),
      hr = exp(mean(v)), hr_lo = exp(quantile(v, 0.025)), hr_hi = exp(quantile(v, 0.975)),
      per = if (binary) "1 logit unit of P(any pressor)" else if (kind == "value") "1 SD of log marker" else "1 log-unit per day",
      n_patients = u$n_patients, n_deaths = u$n_deaths, rhat_gate = gate_assoc, association_rhat = u$association_rhat)
  }
}

association_hr  <- bind_rows(assoc_rows)   %>% mutate(grid = JM_GRID, horizon_days = JM_HORIZON, baseline_form = BASELINE_FORM, site = site_name)
movement        <- bind_rows(movement_rows) %>% mutate(grid = JM_GRID, horizon_days = JM_HORIZON, baseline_form = BASELINE_FORM, site = site_name)
level_contrast  <- bind_rows(level_rows)   %>% mutate(grid = JM_GRID, horizon_days = JM_HORIZON, baseline_form = BASELINE_FORM, site = site_name)
if (nrow(level_contrast)) {
  report_write(level_contrast, "level_contrast")
  message("--- marker difference per unit of the size exposure at each horizon (log units; ",
          "per SD for log_pfvc_sd / ldisc_sd, per log unit for the ch_* pieces; p_equal tests the four pieces equal)")
  print(as.data.frame(level_contrast %>% select(marker, adjustment, exposure, horizon_h, estimate, lo, hi, p_gt0, p_equal, n_patients, n_deaths) %>%
                        mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
}
if (nrow(movement)) report_write(movement, "movement")
lme_check <- bind_rows(lme_check_rows) %>% mutate(grid = JM_GRID, horizon_days = JM_HORIZON, site = site_name)
if (nrow(lme_check)) {
  report_write(lme_check, "lme_check")
  message("--- the size terms with and without the death correction (joint model against the LME alone)")
  print(as.data.frame(lme_check %>% select(marker, adjustment, exposure, term, jm_estimate, lme_estimate, jm_minus_lme_in_lme_se,
                                           association_rhat, hazard_rhat) %>% mutate(across(where(is.numeric), ~ signif(., 3)))),
        row.names = FALSE)
}
report_write(association_hr,  "association_hr")

message("23_biotrauma_report complete: ", nrow(level_contrast), " contrast rows, ",
        nrow(association_hr), " association terms -> ", final_dir)
