# =============================================================================
# Script 32_tte_run_all: PRIMARY TTE refresh -- source the engine ONCE, then the primary leaves.
# =============================================================================
# Driver for the federated per-site TTE deliverable. Sources the shared engine once (so
# the expensive baseline/panel/primary-design build runs a single time) then the leaves:
#
#   33_tte_ceiling.R          PRIMARY: the normalizer head-to-head on tidal volume --
#                               PFVC-anchored vs PBW-anchored VT ceiling (bite-matched), the
#                               cap-on-top secondary, and the weight-refit bootstrap that is
#                               the primary interval (PBWPFVC_TTE_EXPO_FAMILY=vt)
#   34_tte_titration.R        CO-PRIMARY: one-step titration toward the PFVC-anchored target
#                               (bounded modified treatment policy, LMTP); the plain additive
#                               shift is its sensitivity (PBWPFVC_VT_POLICY=additive)
#   35_tte_primary.R              strain-limiting (VT/PFVC <= 11%) vs permissive (<= 16%) ceiling
#   36_tte_diagnostics.R          its diagnostics
#   37_tte_discordance_benefit.R  its discordance-HTE reads
#
# Every leaf sources 31_tte_engine.R unconditionally; the engine restores the cached build for
# this site, so the repeat sources are cheap. The mechanical-power family of 33_tte_ceiling and the
# sensitivity suite (tte_sensitivities.R) are lead-site supplements, not part of this bundle.
#
# RNG note: the only RNG consumers are the synthetic survival sim (§10a, in the engine) and
# the cluster bootstraps, each on a fixed clusterSetRNGStream; leaf order does not change any
# output. The weight-refit bootstrap rebuilds both arms per replicate and is the slow piece;
# at a large site set PBWPFVC_CEIL_REFIT_BOOT (e.g. 200) rather than skip it.
# =============================================================================
library(here)
source(here::here("code", "31_tte_engine.R"))
Sys.setenv(PBWPFVC_TTE_EXPO_FAMILY = "vt")

leaves <- c(
  "33_tte_ceiling.R",
  "34_tte_titration.R",
  "35_tte_primary.R",
  "36_tte_diagnostics.R",
  "37_tte_discordance_benefit.R")

for (leaf in leaves) {
  message("\n>>> running ", leaf, " ...")
  source(here::here("code", leaf))
}
message("\n32_tte_run_all: completed all ", length(leaves), " primary TTE leaves.")
