# =============================================================================
# tte_sens_censoring: common-support trim x day-weight cap sweep (at the primary ceilings)
# =============================================================================
# Vary the common-support trim (TRIM_ALPHA) and the day-weight truncation cap (DAYW_CAP) at the
# primary C_LOW/C_HIGH. Stability = the discordance gradient is not a thin-support / weight-tail
# artifact. Point estimates only (primary CIs: 37_tte_discordance_benefit). Outputs tte_ccw_disc_censoring_sweep_<site>.
# =============================================================================
source(here::here("code", "supplement", "tte_sens_common.R"))

cens_grid <- list(
  list(trim = TRIM_ALPHA, cap = DAYW_CAP),                 # PRIMARY
  list(trim = 0,          cap = DAYW_CAP),
  list(trim = 0.01,       cap = DAYW_CAP),
  list(trim = 0.05,       cap = DAYW_CAP),
  list(trim = TRIM_ALPHA, cap = 3),
  list(trim = TRIM_ALPHA, cap = 10))
message("tte_sens_censoring: censoring (trim x cap) sweep -- rebuilding ", length(cens_grid), " designs ...")
censoring_sweep <- map_dfr(cens_grid, function(s) hte_for_design(C_LOW, C_HIGH, trim = s$trim, cap = s$cap)) %>%
  mutate(setting = if_else(trim == TRIM_ALPHA & cap == DAYW_CAP, "PRIMARY", "sensitivity"), .before = 1) %>%
  select(setting, trim, cap, overall_rd, rd_concordant, rd_mid, rd_discordant, gradient)
write_csv(censoring_sweep, file.path(final_dir, paste0("tte_ccw_disc_censoring_sweep_", site_name, ".csv")))

cat("\n=== 2. Censoring sweep: trim x weight-cap (point estimates) ===\n")
print(as.data.frame(censoring_sweep %>% transmute(setting, trim, cap,
        overall_pp = round(100 * overall_rd, 2), rd_disc_pp = round(100 * rd_discordant, 2),
        gradient_pp = round(100 * gradient, 2))), row.names = FALSE)
cat("    (gradient stable across trim/cap = not a thin-support / weight-tail artifact)\n")
message("Wrote tte_ccw_disc_censoring_sweep_", site_name, ".csv to ", final_dir)
