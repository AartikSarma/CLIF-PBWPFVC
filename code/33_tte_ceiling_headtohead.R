# =============================================================================
# 33_tte_ceiling_headtohead: TTE of a PFVC-anchored vs a PBW-anchored ceiling (bite-matched)
# =============================================================================
# One piece of the split 33_tte_ceiling (see 33_tte_ceiling_common.R for the design + the shared estimator).
# Runs ONLY the "ceiling" design's fixed-weight bootstrap and writes tte_<fam>_ceiling_*. The
# primary interval comes from the weight-refit bootstrap in 33_tte_ceiling_diagnostics.R.
# Family: PBWPFVC_TTE_EXPO_FAMILY (vt = primary, mp = secondary).
# =============================================================================
library(here)
source(here::here("code", "33_tte_ceiling_common.R"))
res_ceiling <- analyse_design("ceiling")
message("Wrote ", PFX, "ceiling_{overall,disc_hte,disc_gradient,disc_cate_curve,disc_cate_slope,overlap_*,disc_dose}_",
        site_name, ".csv + cuminc / disc_benefit .pdf to ", final_dir)
