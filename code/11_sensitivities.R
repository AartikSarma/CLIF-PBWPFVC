# =============================================================================
# 11_sensitivities: run-all wrapper for the split TTE sensitivity pieces (lead-site / supplement)
# =============================================================================
# NOT part of the per-site primary bundle (NOT in 11_run_all). Runs EVERY sensitivity in one shot.
# Each sweep also runs STANDALONE -- source the individual 11_sens_*.R file to rerun just that test
# without rebuilding the others (the engine's primary design is built once and reused via the
# guard in 11_sens_common.R). All pieces recompute the 11.X discordance-gradient statistics with
# the SAME pooled, SOFA-adjusted, IPC-weighted standardized estimator and write the SAME output
# filenames the monolithic version did (so the report + pooling are unaffected).
#
# Pieces (each independently runnable):
#   11_sens_thresholds.R       1. strain-threshold sweep (C_LOW x C_HIGH)
#   11_sens_censoring.R        2. common-support trim x day-weight cap
#   11_sens_weightmodel.R      3-4. richer S/F deviation model + weight-timing diagnostic
#   11_sens_deadspace.R        5. ventilatory-ratio (dead-space) confounder, OUT vs IN the weights
#   11_sens_severity_ladder.R  6. time-varying severity-confounder ladder (daily SOFA / VR / DP)
# Future: fold in the remaining archived 11.C-W sensitivities (see code/SENSITIVITY_INVENTORY.md).
# =============================================================================
library(here)
source(here::here("code", "11_sens_common.R"))   # build the primary design once; helpers in scope

for (piece in c("11_sens_thresholds", "11_sens_censoring", "11_sens_weightmodel",
                "11_sens_deadspace", "11_sens_severity_ladder")) {
  message("\n================  ", piece, "  ================")
  source(here::here("code", paste0(piece, ".R")))
}
message("\n11_sensitivities: all pieces complete.")
