# =============================================================================
# Script 11_run_all: PRIMARY TTE refresh -- source common ONCE, then the primary leaves.
# =============================================================================
# Driver for the federated per-site TTE deliverable. Sources the shared engine once (so
# the expensive baseline/panel/primary-design build runs a single time) then the PRIMARY
# leaves: 11.A (primary TTE), 11.B (diagnostics), 11.X (discordance-HTE primary). Each
# leaf carries an `if (!exists("build_design")) source(common)` guard so it also runs
# standalone; here that guard short-circuits because the engine is already in .GlobalEnv.
#
# The ~20 TTE sensitivities/diagnostics that used to run here were moved to the gitignored
# code/archive/ (see code/SENSITIVITY_INVENTORY.md). They are to be rebuilt as a single
# consolidated sensitivity script once the primary thresholds (ARMA-derived) are finalized.
#
# RNG note: the only RNG consumers are the synthetic survival sim (§10a, in common) and the
# cluster bootstraps (11.A and 11.X), each on a fixed clusterSetRNGStream; leaf order does
# not change any output.
# =============================================================================
library(here)
source(here::here("code", "10_tte_engine.R"))

leaves <- c(
  "11.A_primary.R",
  "11.B_diagnostics.R",
  "11.X_discordance_benefit.R")

for (leaf in leaves) {
  message("\n>>> running ", leaf, " ...")
  source(here::here("code", leaf))
}
message("\n11_run_all: completed all ", length(leaves), " primary TTE leaves.")
