# =============================================================================
# Script 11_run_all: full TTE refresh -- source common ONCE, then every 11.* leaf.
# =============================================================================
# Driver for the split TTE pipeline (10_tte_common.R + 11.*_*.R). Sourcing
# common once (so the expensive baseline/panel/primary-design build runs a single time)
# and then each analysis leaf in order produces the full result set. The 11.*
# leaves carry an `if (!exists("build_design")) source(common)` guard so each also runs
# standalone; here that guard short-circuits because common is already in .GlobalEnv.
#
# RNG note: the only RNG consumers are the synthetic survival sim (§10a, in common) and
# the primary cluster bootstrap (§10e, in 11.A, on a fixed clusterSetRNGStream). The
# leaf order does not change any output because the sensitivities/diagnostics draw no
# main-process RNG. 11.N (refit calibration) reuses `overall` populated by 11.A.
# =============================================================================
library(here)
source(here::here("code", "10_tte_engine.R"))

leaves <- c(
  "11.A_primary.R",
  "11.B_diagnostics.R",
  "11.C_sens_weightcap.R",
  "11.D_sens_ceiling_grace.R",
  "11.E_sens_rule.R",
  "11.F_sens_aggregation.R",
  "11.G_sens_cumweight.R",
  "11.H_sens_ph.R",
  "11.I_sens_dp.R",
  "11.J_sens_pf.R",
  "11.K_sens_numerator.R",
  "11.L_sens_trim.R",
  "11.M_sens_mtp.R",
  "11.N_refit_boot.R",
  "11.O_pfvc_overlap.R",
  "11.P_sens_pfvc_overlap.R",
  "11.Q_sens_dailysofa.R",
  "11.R_age_strata_balance.R",
  "11.S_sofa_trajectory.R",
  "11.T_demand_selection.R",
  "11.U_balance_reference.R",
  "11.V_grace_compare.R",
  "11.X_discordance_benefit.R")

for (leaf in leaves) {
  message("\n>>> running ", leaf, " ...")
  source(here::here("code", leaf))
}
message("\n11_run_all: completed all ", length(leaves), " TTE analysis leaves.")
