# =============================================================================
# Script 20 (grid): the time grid shared by the biotrauma panel, fit and report
# =============================================================================
# Sourced by 21_biotrauma_panel.R, 22_biotrauma_fit.R, 23_biotrauma_report.R,
# 24_biotrauma_figures.R, 27_control_comparison.R, 28_height_fingerprint.R and
# pooling/pooled_biotrauma.R (and
# by supplement/xsec_crs_channels.R, supplement/xsec_pfvc_age_control.R and
# tools/creatinine_positive_control.R for pfvc_channels() or h_suffix), so they
# agree on the grid, the horizon, the output suffix and the cohort restrictions.
#
#   PBWPFVC_JM_GRID       "daily" (figure 4, the default) or "6h" (six-hour periods;
#                         not in the paper)
#   PBWPFVC_JM_HORIZON_H  horizon in hours for the 6h grid (48)
#   PBWPFVC_JM_HORIZON    horizon in days for the daily grid (7)
#   PBWPFVC_JM_ICU_DAY0   0 | 1   keep patients on the cohort's support at ICU admission
#                         (icu_day0, script 03): the ventilated arm of every
#                         ventilated-vs-control comparison; file tag "day0_"
#   PBWPFVC_JM_NO_LAGS    0 | 1   drop the previous-day SF and pressor terms from the
#                         longitudinal submodel (sensitivity); file tag "nolag_"
#
# Defines: JM_GRID, STEP_H (hours per period), STEP (days per period),
# JM_HORIZON (days), N_PERIODS (last period index), h_suffix ("48h" / "7d"), and the
# convergence gate every script applies (RHAT_GATE, size_terms_for(), fit_convergence()).
# Time in every model is `vent_day` in days (period x STEP), so coefficients on
# time and the random slope have the same units on both grids.
# =============================================================================
JM_GRID <- Sys.getenv("PBWPFVC_JM_GRID", "daily")
stopifnot(JM_GRID %in% c("6h", "daily"))
if (JM_GRID == "6h") {
  STEP_H     <- 6
  horizon_h  <- as.numeric(Sys.getenv("PBWPFVC_JM_HORIZON_H", "48"))
  stopifnot(is.finite(horizon_h), horizon_h >= 12, horizon_h %% STEP_H == 0)
  JM_HORIZON <- horizon_h / 24
  h_suffix   <- paste0(as.integer(horizon_h), "h")
} else {
  STEP_H     <- 24
  JM_HORIZON <- as.integer(Sys.getenv("PBWPFVC_JM_HORIZON", "7"))
  stopifnot(is.finite(JM_HORIZON), JM_HORIZON >= 2L)
  h_suffix   <- paste0(JM_HORIZON, "d")
}
STEP      <- STEP_H / 24
N_PERIODS <- as.integer(round(JM_HORIZON / STEP))

# =============================================================================
# Channel decomposition of the size exposure (the channels form of 22 and 23)
# =============================================================================
# GLI-2012 log PFVC is, to a small remainder, a sum of four pieces: a height
# term, an age curve, a sex shift and a race shift. Each piece is computed by
# moving one input while the other three sit at reference values (cohort median
# height and age, male, white), so every piece is in log-PFVC units. A model
# that gives each piece its own coefficient is the saturated demographic model
# with GLI-shaped functional forms; if lung size is the operative quantity the
# four coefficients are equal, and the model collapses to one beta on log PFVC.
# The same pieces are built for log PBW/PFVC (VT/PFVC at a given VT/PBW): the
# PBW pieces (Devine: height and sex) minus the PFVC pieces.
#
# pfvc_channels(d, expo) takes height_cm, age10, sex_category, race_category and
# returns ch_height, ch_age, ch_sex, ch_race (log units), ch_sum, and
# ch_remainder = exposure (centred at the reference patient) - ch_sum, which is
# non-zero only because GLI's height and age coefficients differ slightly by sex.
CHANNELS <- c("ch_height", "ch_age", "ch_sex", "ch_race")
pfvc_channels <- function(d, expo = c("log_pfvc", "ldisc")) {
  expo <- match.arg(expo)
  stopifnot(all(c("height_cm", "age10", "sex_category", "race_category") %in% names(d)))
  sex  <- if_else(d$sex_category == "Female", 2L, 1L)
  race <- case_when(d$race_category == "WHITE" ~ 1L, d$race_category == "BLACK" ~ 2L, TRUE ~ 5L)   # GLI codes, as in 03
  age  <- d$age10 * 10
  h    <- d$height_cm
  n    <- nrow(d)
  ref  <- list(h = median(h, na.rm = TRUE), age = median(age, na.rm = TRUE), sex = 1L, race = 1L)
  gli  <- function(age, h, sex, race) log(rspiro::pred_GLI(age = age, height = h / 100, gender = sex, ethnicity = race, param = "FVC"))
  dev  <- function(h, sex) log(if_else(sex == 1L, 50, 45.5) + 2.3 * (h / 2.54 - 60))       # Devine PBW, kg
  l0   <- gli(ref$age, ref$h, ref$sex, ref$race)
  out  <- tibble(
    ch_height = gli(rep(ref$age, n), h, rep(ref$sex, n), rep(ref$race, n)) - l0,
    ch_age    = gli(age, rep(ref$h, n), rep(ref$sex, n), rep(ref$race, n)) - l0,
    ch_sex    = gli(rep(ref$age, n), rep(ref$h, n), sex, rep(ref$race, n)) - l0,
    ch_race   = gli(rep(ref$age, n), rep(ref$h, n), rep(ref$sex, n), race) - l0)
  exact <- gli(age, h, sex, race) - l0
  if (expo == "ldisc") {
    p0 <- dev(ref$h, ref$sex)
    out <- out %>% mutate(ch_height = (dev(h, ref$sex) - p0) - ch_height,
                          ch_sex    = (dev(ref$h, sex) - p0) - ch_sex,
                          ch_age    = -ch_age, ch_race = -ch_race)
    exact <- (dev(h, sex) - p0) - exact
  }
  out %>% mutate(ch_sum = ch_height + ch_age + ch_sex + ch_race, ch_remainder = exact - ch_sum)
}
# The height fingerprint of log PBW/PFVC (28_height_fingerprint.R): the ratio's
# height piece at each patient's OWN sex, with age and race at fixed reference
# values. ch_height above holds sex at male, which erases what this term is for:
# Devine PBW is a line with an intercept and GLI FVC a power law in height, so the
# ratio is an inverted U in height that peaks near 166 cm in men and 183 cm in
# women. GLI's log-height coefficient does not depend on age or race, so the
# reference age (60 years) and race (white) move each sex's curve by a constant and
# nothing else, so any fixed choice gives the same fingerprint up to that constant,
# which the sex term of any model that uses the fingerprint absorbs.
# Returned in log units, relative to a 170 cm man.
height_fingerprint <- function(height_cm, sex_category) {
  bad_sex <- setdiff(unique(sex_category), c("Male", "Female"))
  if (length(bad_sex)) stop("height_fingerprint(): GLI has two sex equations; got sex_category ",
                            paste(bad_sex, collapse = ", "))
  sex <- if_else(sex_category == "Female", 2L, 1L)
  log_ratio <- function(h, sex) {
    n <- length(h)
    log(if_else(sex == 1L, 50, 45.5) + 2.3 * (h / 2.54 - 60)) -
      log(rspiro::pred_GLI(age = rep(60, n), height = h / 100, gender = sex, ethnicity = rep(1L, n), param = "FVC"))
  }
  log_ratio(height_cm, sex) - log_ratio(170, 1L)
}

# ---- The convergence gate ------------------------------------------------------
# One gate for every script that reads a joint model (22, 23, 24, 27 and the
# pooling): a fit counts as converged when the R-hat of its LUNG-SIZE terms, the level
# and the divergence the figure reads, is at or below RHAT_GATE. The hazard links (the
# survival submodel and the association of the marker with each cause) are reported
# beside it as hazard_rhat and read with the longitudinal-only comparison
# (jm_lme_check_*, 23_biotrauma_report.R), which shows whether they move the estimate.
RHAT_GATE <- 1.1   # the standard convergence threshold
# the size terms per modifier form (22_biotrauma_fit.R, mod_terms): each form's size
# level and its divergence (the level x vent_day interaction); an interaction is
# matched in either order, as R writes the pair by appearance.
#   pfvc, pfvc_dose  log_pfvc_sd, log_pfvc_sd:vent_day (pfvc_dose's dose-modified
#                    terms are not size terms)
#   disc_level       ldisc_sd, ldisc_sd:vent_day
#   vtpfvc           vtpfvc_c, vtpfvc_c:vent_day
#   channels         ch_height, ch_age, ch_sex, ch_race and each x vent_day
# The dose-modification forms (not in the paper) have no divergence; their size terms are the
# ones they read:
#   disc             ldisc_c, l_vtpbw_within:ldisc_c
#   saturated        log_pbw, log_pfvc, l_vtpbw_within:log_pbw, l_vtpbw_within:log_pfvc
#   none             no size term; the dose slope l_vtpbw_within is what it reads
SIZE_TERMS <- list(
  pfvc       = c("log_pfvc_sd", "log_pfvc_sd:vent_day"),
  pfvc_dose  = c("log_pfvc_sd", "log_pfvc_sd:vent_day"),
  disc_level = c("ldisc_sd", "ldisc_sd:vent_day"),
  vtpfvc     = c("vtpfvc_c", "vtpfvc_c:vent_day"),
  channels   = c(CHANNELS, paste0(CHANNELS, ":vent_day")),
  disc       = c("ldisc_c", "l_vtpbw_within:ldisc_c"),
  saturated  = c("log_pbw", "log_pfvc", "l_vtpbw_within:log_pbw", "l_vtpbw_within:log_pfvc"),
  none       = "l_vtpbw_within")
# a term with its components sorted, so vent_day:log_pfvc_sd matches log_pfvc_sd:vent_day
sorted_term <- function(term) vapply(strsplit(term, ":"), function(p) paste(sort(p), collapse = ":"), character(1))
size_terms_for <- function(form) {
  if (!form %in% names(SIZE_TERMS)) stop("unknown modifier form '", form, "': no size terms")
  SIZE_TERMS[[form]]
}
# One fit's size-term R-hat: the maximum over the size terms of the longitudinal block.
# `term`, `rhat` and `block` are the columns of one fit's rows of jm_estimates_*
# (use inside summarise() grouped by fit). A fit whose table has none of the form's
# size terms is an error, not a pass.
size_rhat_of <- function(term, rhat, block, form) {
  keep <- block == "longitudinal" & sorted_term(term) %in% sorted_term(size_terms_for(form))
  if (!any(keep)) stop("no ", form, "-form size term among the fit's estimates (",
                       paste(head(unique(term), 8), collapse = ", "), ", ...)")
  max(rhat[keep])
}
# One fit's hazard R-hat: the maximum over the survival submodel and the association
# parameters (value(log_y):strata...)
hazard_rhat_of <- function(rhat, block) {
  keep <- block %in% c("survival", "association")
  if (any(keep)) max(rhat[keep]) else NA_real_
}
passes_rhat_gate <- function(rhat) !is.na(rhat) & is.finite(rhat) & rhat <= RHAT_GATE
# The gate for every fit in an estimates table: one row per fit (the `by` columns) with
# size_terms_rhat, size_gate, hazard_rhat and hazard_gate
fit_convergence <- function(est, form, by = c("marker", "model", "adjustment")) {
  est %>% group_by(across(all_of(by))) %>%
    summarise(size_terms_rhat = size_rhat_of(term, rhat, block, form),
              hazard_rhat     = hazard_rhat_of(rhat, block), .groups = "drop") %>%
    mutate(size_gate = passes_rhat_gate(size_terms_rhat), hazard_gate = passes_rhat_gate(hazard_rhat))
}

# Wald test that the four channel coefficients (or contrasts) are equal:
# est is a length-4 vector, V its covariance; returns the chi-square p-value on 3 df
channels_equal_p <- function(est, V) {
  L <- rbind(c(1, -1, 0, 0), c(1, 0, -1, 0), c(1, 0, 0, -1))
  d <- L %*% est
  as.numeric(pchisq(t(d) %*% solve(L %*% V %*% t(L)) %*% d, df = 3, lower.tail = FALSE))
}

# ---- Cohort restrictions shared by the fit, the report and the figures ----------
# The severity anchor, PER MARKER: the sum of the index-day cardiovascular,
# coagulation, liver and renal SOFA components, leaving out the marker's OWN
# component, so the anchor never contains the outcome. Two components are never in an
# anchor: respiratory (collinear with SF) and neurological (on the day of intubation
# the worst GCS is sedation). The oxygenation and mechanics markers have no component
# to drop.
ANCHOR_POOL <- c("sofa_cv_97", "sofa_coag", "sofa_liver", "sofa_renal")
ANCHOR_DROP <- c(creatinine = "sofa_renal", platelets = "sofa_coag", bilirubin = "sofa_liver",
                 ne_equiv_peak = "sofa_cv_97", any_pressor = "sofa_cv_97", pressor_dose = "sofa_cv_97")
anchor_components <- function(marker) setdiff(ANCHOR_POOL, unname(ANCHOR_DROP[marker]))
anchor_label <- function(marker)
  paste(sub("_97", "", sub("sofa_", "", anchor_components(marker))), collapse = " + ")
# Severity standardisation of the control (PBWPFVC_JM_SEV_CENTER = "creatinine=2.61,
# platelets=3.05,..."; figure 4 uses it). The
# objection to an unmatched control is effect modification, not confounding: PFVC is
# fixed by height, age, sex and race, so illness cannot move it, but a smaller lung
# might show in the trajectory only under physiological stress, and a healthier
# control could be null for that reason alone. So the control keeps EVERY patient and
# its divergence is allowed to vary with the marker's own anchor, which is centred at
# the VENTILATED cohort's mean anchor (22_biotrauma_fit.R, section 22f):
#   log_pfvc_sd:vent_day             the control's rate at the ventilated severity;
#                                    linear in the anchor, so this is also the rate
#                                    averaged over the ventilated anchor distribution
#   log_pfvc_sd:vent_day:sev_anchor_c  does the rate grow with severity? (the test of
#                                    the objection itself)
# The centre enters the cache names, so a new ventilated cohort refits the control.
SEV_CENTER_SPEC <- trimws(Sys.getenv("PBWPFVC_JM_SEV_CENTER", ""))
sev_centers <- if (nzchar(SEV_CENTER_SPEC)) {
  pairs <- strsplit(trimws(strsplit(SEV_CENTER_SPEC, ",")[[1]]), "=")
  centers <- setNames(suppressWarnings(as.numeric(vapply(pairs, `[`, "", 2))), trimws(vapply(pairs, `[`, "", 1)))
  if (anyNA(centers) || any(!nzchar(names(centers)))) stop("PBWPFVC_JM_SEV_CENTER must be 'marker=number,...'; got '", SEV_CENTER_SPEC, "'")
  centers
} else numeric()
sev_center_for <- function(marker) {
  if (!nzchar(SEV_CENTER_SPEC)) return(NA_real_)
  if (!marker %in% names(sev_centers)) stop("PBWPFVC_JM_SEV_CENTER has no centre for ", marker)
  sev_centers[[marker]]
}
sev_center_sfx_for <- function(marker) { v <- sev_center_for(marker); if (is.na(v)) "" else sprintf("_sevstd%.3f", v) }
sev_center_tag <- if (nzchar(SEV_CENTER_SPEC)) "sevstd_" else ""

# Status at ICU admission (PBWPFVC_JM_ICU_DAY0 = 1): keep icu_day0 patients. In the
# ventilated cohort these are the patients on invasive ventilation at ICU admission,
# the ventilated arm of every comparison with the no-support control (whose patients
# are all indexed at ICU admission, so the filter keeps all of them).
ICU_DAY0 <- identical(Sys.getenv("PBWPFVC_JM_ICU_DAY0", "0"), "1")
icu_day0_tag <- if (ICU_DAY0) "day0_" else ""
icu_day0_sfx <- if (ICU_DAY0) "_day0" else ""

# Baseline SF band (PBWPFVC_JM_SF_BAND = "lo,hi"): keep patients with lo <= SF < hi at
# the index timepoint (sf_index, the SF ratio of script 03's index row), the same gate
# in every cohort. The strata in use are "235,315", "115,235" and "0,115" (315 and 235
# are the Rice 2007 SF equivalents of P/F 300 and 200).
# The upper bound is strict so that "0,315" is the ventilated cohort's own gate, SF < 315.
SF_BAND_RULE <- "lo <= SF < hi"   # stored with each fit: a fit made under another rule is refitted
SF_BAND <- trimws(Sys.getenv("PBWPFVC_JM_SF_BAND", ""))
sf_band_limits <- if (nzchar(SF_BAND)) suppressWarnings(as.numeric(strsplit(SF_BAND, ",")[[1]])) else c(NA_real_, NA_real_)
if (nzchar(SF_BAND) && (length(sf_band_limits) != 2L || anyNA(sf_band_limits) || sf_band_limits[1] >= sf_band_limits[2]))
  stop("PBWPFVC_JM_SF_BAND must be 'lo,hi' with lo < hi; got '", SF_BAND, "'")
sf_sfx <- if (nzchar(SF_BAND)) paste0("_sf", sf_band_limits[1], "to", sf_band_limits[2]) else ""
sf_tag <- if (nzchar(SF_BAND)) paste0("sf", sf_band_limits[1], "to", sf_band_limits[2], "_") else ""

# Sensitivity without the previous-day SF and pressor terms (PBWPFVC_JM_NO_LAGS = 1):
# both are measured after the previous day's dose, so they may carry part of its
# effect. The rows are those of the main fit; only the model changes.
NO_LAGS <- identical(Sys.getenv("PBWPFVC_JM_NO_LAGS", "0"), "1")
no_lags_tag <- if (NO_LAGS) "nolag_" else ""
no_lags_sfx <- if (NO_LAGS) "_nolag" else ""

# every restriction and variant, in the order every script uses: ICU day 0, severity
# standardisation, SF band, no lags
restrict_tag <- paste0(icu_day0_tag, sev_center_tag, sf_tag, no_lags_tag)
restrict_sfx_for <- function(marker) paste0(icu_day0_sfx, sev_center_sfx_for(marker), sf_sfx, no_lags_sfx)
