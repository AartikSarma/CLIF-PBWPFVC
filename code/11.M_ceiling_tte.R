# =============================================================================
# 11.M: run-all wrapper for the normalizer head-to-head TTE pieces
# =============================================================================
# The analysis is split so each test runs on its own without redoing the others:
#   11_ceiling_common.R        shared build (cheap; guarded) + thresholds
#   11_ceiling_headtohead.R    PFVC-anchored vs PBW-anchored ceiling, bite-matched  (tte_<fam>_ceiling_*)
#   11_ceiling_cap.R           PBW ceiling + PFVC cap vs PBW ceiling                (tte_<fam>_cap_*)
#   11_ceiling_diagnostics.R   weight-REFIT bootstrap (the PRIMARY interval) + placebo contrast
#                              for every design                                     (tte_<fam>_*_overall_refit_*, tte_<fam>_diag_*)
# PBWPFVC_TTE_EXPO_FAMILY selects the exposure family (vt = PRIMARY, mp = secondary);
# PBWPFVC_CEIL_PIECES selects which pieces run (default: all three). Run a single piece
# (Rscript code/11_ceiling_cap.R) to rerun just that test.
# =============================================================================
library(here)
PIECES <- strsplit(Sys.getenv("PBWPFVC_CEIL_PIECES", "headtohead,cap,diagnostics"), ",")[[1]]
stopifnot(all(PIECES %in% c("headtohead", "cap", "diagnostics")))
for (p in PIECES) source(here::here("code", paste0("11_ceiling_", p, ".R")))
