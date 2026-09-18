# =============================================================================
# tte_sensitivities: run-all wrapper for the split TTE sensitivity pieces (lead-site / supplement)
# =============================================================================
# NOT part of the per-site primary bundle (NOT in 32_tte_run_all). Runs EVERY sensitivity in one shot.
# Each sweep also runs STANDALONE -- source the individual tte_sens_*.R file to rerun just that test
# without rebuilding the others (the engine's primary design is built once and reused via the
# guard in tte_sens_common.R). All pieces recompute the 37_tte_discordance_benefit discordance-gradient statistics with
# the SAME pooled, SOFA-adjusted, IPC-weighted standardized estimator and write the SAME output
# filenames the monolithic version did (so the report + pooling are unaffected).
#
# Pieces (each independently runnable):
#   tte_sens_thresholds.R       1. strain-threshold sweep (C_LOW x C_HIGH)
#   tte_sens_censoring.R        2. common-support trim x day-weight cap
#   tte_sens_weightmodel.R      3-4. richer S/F deviation model + weight-timing diagnostic
#   tte_sens_deadspace.R        5. ventilatory-ratio (dead-space) confounder, OUT vs IN the weights
#   tte_sens_severity_ladder.R  6. time-varying severity-confounder ladder (daily SOFA / VR / DP)
# Future: fold in the remaining archived 11.C-W sensitivities (see code/supplement/README.md).
# =============================================================================
library(here)
source(here::here("code", "supplement", "tte_sens_common.R"))   # build the primary design once; helpers in scope

for (piece in c("tte_sens_thresholds", "tte_sens_censoring", "tte_sens_weightmodel",
                "tte_sens_deadspace", "tte_sens_severity_ladder")) {
  message("\n================  ", piece, "  ================")
  source(here::here("code", "supplement", paste0(piece, ".R")))
}
message("\ntte_sensitivities: all pieces complete.")
