# =============================================================================
# Script 20 (grid): the time grid shared by the biotrauma panel, fit and report
# =============================================================================
# Sourced by 21_biotrauma_panel.R, 22_biotrauma_fit.R and 23_biotrauma_report.R
# so the three agree on the grid, the horizon and the output suffix.
#
#   PBWPFVC_JM_GRID       "daily" (figure 4, the default) or "6h" (six-hour periods;
#                         not in the paper)
#   PBWPFVC_JM_HORIZON_H  horizon in hours for the 6h grid (48)
#   PBWPFVC_JM_HORIZON    horizon in days for the daily grid (7)
#
# Defines: JM_GRID, STEP_H (hours per period), STEP (days per period),
# JM_HORIZON (days), N_PERIODS (last period index), h_suffix ("48h" / "7d").
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
# reference age and race move each sex's curve by a constant and nothing else;
# the sex term of any model that uses the fingerprint absorbs that constant.
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

# Wald test that the four channel coefficients (or contrasts) are equal:
# est is a length-4 vector, V its covariance; returns the chi-square p-value on 3 df
channels_equal_p <- function(est, V) {
  L <- rbind(c(1, -1, 0, 0), c(1, 0, -1, 0), c(1, 0, 0, -1))
  d <- L %*% est
  as.numeric(pchisq(t(d) %*% solve(L %*% V %*% t(L)) %*% d, df = 3, lower.tail = FALSE))
}

# ---- Cohort restrictions shared by the fit, the report and the figures ----------
# Severity floor (PBWPFVC_JM_SEV_MIN): keep patients whose baseline severity anchor
# is at or above a value, to severity-match the no-support control to the ventilated
# cohort. The anchor is PER MARKER: the sum of the index-day cardiovascular,
# coagulation, liver and renal SOFA components, leaving out the marker's OWN
# component. A floor on a score containing the outcome selects extreme baselines,
# which bends the average time trend through regression to the mean. Two components
# are never in an anchor: respiratory (collinear with SF) and neurological (on the
# day of intubation the worst GCS is sedation, so matching on it matches on the
# treatment). The oxygenation and mechanics markers have no component to drop.
ANCHOR_POOL <- c("sofa_cv_97", "sofa_coag", "sofa_liver", "sofa_renal")
ANCHOR_DROP <- c(creatinine = "sofa_renal", platelets = "sofa_coag", bilirubin = "sofa_liver",
                 ne_equiv_peak = "sofa_cv_97", any_pressor = "sofa_cv_97", pressor_dose = "sofa_cv_97")
anchor_components <- function(marker) setdiff(ANCHOR_POOL, unname(ANCHOR_DROP[marker]))
anchor_label <- function(marker)
  paste(sub("_97", "", sub("sofa_", "", anchor_components(marker))), collapse = " + ")
# The knob is one number (each marker's own anchor gets the same floor) or a list,
# "platelets=2,bilirubin=1". The floor enters every cache name and the output tag, so
# a restricted run never reuses an unrestricted fit.
SEV_SPEC <- trimws(Sys.getenv("PBWPFVC_JM_SEV_MIN", ""))
SEV_BY_MARKER <- grepl("=", SEV_SPEC)
sev_floors <- if (SEV_BY_MARKER) {
  pairs <- strsplit(trimws(strsplit(SEV_SPEC, ",")[[1]]), "=")
  floors <- setNames(suppressWarnings(as.numeric(vapply(pairs, `[`, "", 2))), trimws(vapply(pairs, `[`, "", 1)))
  if (anyNA(floors) || any(!nzchar(names(floors)))) stop("PBWPFVC_JM_SEV_MIN must be a number or 'marker=number,...'; got '", SEV_SPEC, "'")
  floors[order(names(floors))]          # sorted, so the same floors always give the same tag
} else numeric()
sev_floor_for <- function(marker) {
  if (!nzchar(SEV_SPEC)) return(NA_real_)
  if (!SEV_BY_MARKER) return(as.numeric(SEV_SPEC))
  if (!marker %in% names(sev_floors)) stop("PBWPFVC_JM_SEV_MIN lists floors by marker but has none for ", marker)
  sev_floors[[marker]]
}
sev_sfx_for <- function(marker) { v <- sev_floor_for(marker); if (is.na(v)) "" else paste0("_sev", v) }
# The output tag carries the floors themselves ("sev2_", or "sev_bilirubin1_platelets2_"):
# a second choice of floor writes its own tables instead of replacing the first's.
sev_tag <- if (!nzchar(SEV_SPEC)) "" else if (SEV_BY_MARKER)
  paste0("sev_", paste0(names(sev_floors), sev_floors, collapse = "_"), "_") else paste0("sev", SEV_SPEC, "_")

# Severity standardisation of the control (PBWPFVC_JM_SEV_CENTER = "creatinine=2.61,
# platelets=3.05,..."; 2026-09-21, replacing the floor above in figure 4). The
# objection to an unmatched control is effect modification, not confounding: PFVC is
# fixed by height, age, sex and race, so illness cannot move it, but a smaller lung
# might show in the trajectory only under physiological stress, and a healthier
# control could be null for that reason alone. So the control keeps EVERY patient and
# its divergence is allowed to vary with the marker's own anchor, which is centred at
# the VENTILATED cohort's mean anchor (22_biotrauma_fit.R, section 13f):
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
if (nzchar(SEV_CENTER_SPEC) && nzchar(SEV_SPEC)) stop("PBWPFVC_JM_SEV_CENTER and PBWPFVC_JM_SEV_MIN are alternatives: set one")
sev_center_for <- function(marker) {
  if (!nzchar(SEV_CENTER_SPEC)) return(NA_real_)
  if (!marker %in% names(sev_centers)) stop("PBWPFVC_JM_SEV_CENTER has no centre for ", marker)
  sev_centers[[marker]]
}
sev_center_sfx_for <- function(marker) { v <- sev_center_for(marker); if (is.na(v)) "" else sprintf("_sevstd%.3f", v) }
sev_center_tag <- if (nzchar(SEV_CENTER_SPEC)) "sevstd_" else ""

# Baseline SF band (PBWPFVC_JM_SF_BAND = "lo,hi"): keep patients with lo <= SF < hi on
# the index day. The strata in use are "235,315", "115,235" and "0,115" (user-specified,
# 2026-09-18; 315 and 235 are the Rice 2007 SF equivalents of P/F 300 and 200).
# The upper bound is strict so that "0,315" is the ventilated cohort's own gate, SF < 315.
SF_BAND_RULE <- "lo <= SF < hi"   # stored with each fit: a fit made under another rule is refitted
SF_BAND <- trimws(Sys.getenv("PBWPFVC_JM_SF_BAND", ""))
sf_band_limits <- if (nzchar(SF_BAND)) suppressWarnings(as.numeric(strsplit(SF_BAND, ",")[[1]])) else c(NA_real_, NA_real_)
if (nzchar(SF_BAND) && (length(sf_band_limits) != 2L || anyNA(sf_band_limits) || sf_band_limits[1] >= sf_band_limits[2]))
  stop("PBWPFVC_JM_SF_BAND must be 'lo,hi' with lo < hi; got '", SF_BAND, "'")
sf_sfx <- if (nzchar(SF_BAND)) paste0("_sf", sf_band_limits[1], "to", sf_band_limits[2]) else ""
sf_tag <- if (nzchar(SF_BAND)) paste0("sf", sf_band_limits[1], "to", sf_band_limits[2], "_") else ""
# both restrictions, in the order every script uses
restrict_tag <- paste0(sev_tag, sev_center_tag, sf_tag)
restrict_sfx_for <- function(marker) paste0(sev_sfx_for(marker), sev_center_sfx_for(marker), sf_sfx)
