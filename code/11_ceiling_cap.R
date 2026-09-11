# =============================================================================
# 11_ceiling_cap: TTE of [PBW ceiling + PFVC safety cap] vs [PBW ceiling alone] (cap-on-top)
# =============================================================================
# One piece of the split 11.M (see 11_ceiling_common.R for the design + the shared estimator).
# Runs ONLY the "cap" design's fixed-weight bootstrap and writes tte_<fam>_cap_*. The Concordant
# tertile is the design's negative control; 11_ceiling_diagnostics.R tests whether its RD is a
# weight artifact and supplies the primary (weight-refit) interval.
# Family: PBWPFVC_TTE_EXPO_FAMILY (vt = primary, mp = secondary).
# =============================================================================
library(here)
source(here::here("code", "11_ceiling_common.R"))
res_cap <- analyse_design("cap")
message("Wrote ", PFX, "cap_{overall,disc_hte,disc_gradient,disc_cate_curve,disc_cate_slope,overlap_*,disc_dose}_",
        site_name, ".csv + cuminc / disc_benefit .pdf to ", final_dir)
