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
