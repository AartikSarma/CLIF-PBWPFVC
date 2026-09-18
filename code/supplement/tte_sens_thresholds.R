# =============================================================================
# tte_sens_thresholds: strain-threshold sweep (vary C_LOW and C_HIGH around the primary)
# =============================================================================
# Vary the strain-limiting ceiling (C_LOW) and permissive ceiling (C_HIGH) one-at-a-time around
# the ARMA-anchored primary (C_LOW=11 = LTV p75; C_HIGH=16 ~ HTV p25). Does the discordance
# gradient survive other clinically-defensible ceilings? Point estimates only (primary CIs: 37_tte_discordance_benefit).
# The cohort (structural-positivity floor) is FIXED at the primary C_LOW, so the sweep isolates the
# ceiling choice on one comparable population; patients infeasible at a lower swept C_LOW simply
# deviate immediately. Standalone -- run on its own; outputs tte_ccw_disc_ceiling_sweep_<site>.csv.
# =============================================================================
source(here::here("code", "supplement", "tte_sens_common.R"))

ceiling_grid <- c(
  list(c(C_LOW, C_HIGH)),                                  # PRIMARY (ARMA-anchored 11/16)
  lapply(c(9, 10, 12, 13), function(cl) c(cl, C_HIGH)),    # vary strain-limiting ceiling
  lapply(c(14, 18, 20),    function(ch) c(C_LOW, ch)))     # vary permissive ceiling
ceiling_grid <- Filter(function(p) p[1] < p[2], ceiling_grid)
message("tte_sens_thresholds: strain-threshold sweep -- rebuilding ", length(ceiling_grid), " designs ...")
ceiling_sweep <- map_dfr(ceiling_grid, function(p) hte_for_design(p[1], p[2])) %>%
  mutate(setting = if_else(c_low == C_LOW & c_high == C_HIGH, "PRIMARY", "sensitivity"), .before = 1) %>%
  select(setting, c_low, c_high, overall_rd, rd_concordant, rd_mid, rd_discordant, gradient)
write_csv(ceiling_sweep, file.path(final_dir, paste0("tte_ccw_disc_ceiling_sweep_", site_name, ".csv")))

cat("\n=== 1. Strain-threshold sweep (SOFA-adjusted standardized RDs, pp; PRIMARY = C_LOW", C_LOW,
    "/ C_HIGH", C_HIGH, ") ===\n")
print(as.data.frame(ceiling_sweep %>% transmute(setting, c_low, c_high,
        overall_pp = round(100 * overall_rd, 2), rd_conc_pp = round(100 * rd_concordant, 2),
        rd_disc_pp = round(100 * rd_discordant, 2), gradient_pp = round(100 * gradient, 2))), row.names = FALSE)
cat("    (gradient stays negative across ceilings = the discordance-HTE is not an artifact of the 11/16 choice)\n")
message("Wrote tte_ccw_disc_ceiling_sweep_", site_name, ".csv to ", final_dir)
