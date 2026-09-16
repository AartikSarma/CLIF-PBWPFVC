# =============================================================================
# Script 13 (grid): the time grid shared by the biotrauma panel, fit and report
# =============================================================================
# Sourced by 13_biotrauma_panel.R, 13_biotrauma_fit.R and 13_biotrauma_report.R
# so the three agree on the grid, the horizon and the output suffix.
#
#   PBWPFVC_JM_GRID       "6h" (PRIMARY) or "daily" (sensitivity)
#   PBWPFVC_JM_HORIZON_H  horizon in hours for the 6h grid (48)
#   PBWPFVC_JM_HORIZON    horizon in days for the daily grid (7)
#
# Defines: JM_GRID, STEP_H (hours per period), STEP (days per period),
# JM_HORIZON (days), N_PERIODS (last period index), h_suffix ("48h" / "7d").
# Time in every model is `vent_day` in days (period x STEP), so coefficients on
# time and the random slope have the same units on both grids.
# =============================================================================
JM_GRID <- Sys.getenv("PBWPFVC_JM_GRID", "6h")
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
# Channel decomposition of the size exposure (13_injury_at_horizon.R, 13_quick_lme.R)
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
# Wald test that the four channel coefficients (or contrasts) are equal:
# est is a length-4 vector, V its covariance; returns the chi-square p-value on 3 df
channels_equal_p <- function(est, V) {
  L <- rbind(c(1, -1, 0, 0), c(1, 0, -1, 0), c(1, 0, 0, -1))
  d <- L %*% est
  as.numeric(pchisq(t(d) %*% solve(L %*% V %*% t(L)) %*% d, df = 3, lower.tail = FALSE))
}
