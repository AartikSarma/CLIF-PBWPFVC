# =============================================================================
# 11_sens_deadspace: ventilatory-ratio (dead-space) confounder sensitivity
# =============================================================================
# VR = (VE x PaCO2) / (PBW x 100 x 37.5) indexes dead space. Clinicians often escalate ventilation
# (-> strain -> deviation) off CO2 clearance, so VR is a plausible deviation confounder the primary
# weight model OMITS -- and l_sf may be a noisy proxy for it. Add lagged VR to the deviation/weight
# model (engine `conf` hook) and test, on the gas-documented subset: (a) HTE robustness -- RD/
# gradient with VR IN vs NOT-IN the weights (same day-set); (b) does VR ABSORB the residual l_sf
# imbalance. CAVEAT: PaCO2 is an arterial gas (MNAR/sparse) + conf-eligibility (lagged VR recorded)
# restricts further, so this is a PARTIAL test on a selected subset of a full-cohort residual.
# Outputs tte_ccw_disc_vr_{hte,sfbalance}_<site>.csv.
# =============================================================================
source(here::here("code", "11_sens_common.R"))
library(sandwich); library(lmtest)

vr       <- make_vr_panel()
vr_daily <- vr$vr_daily; panel_vr <- vr$panel_vr
cat(sprintf("\n=== 5. Ventilatory-ratio confounder sensitivity (gas subset: %d patients, %d patient-days with VR) ===\n",
            n_distinct(vr_daily$hospitalization_id), nrow(vr_daily)))

# same day-set (lagged-VR-recorded); VR NOT-in vs IN the weight model (engine conf_in_model toggle)
des_vr_base <- tryCatch(build_design(C_LOW, C_HIGH, pnl = panel_vr, conf = "vr", conf_in_model = FALSE),
                        error = function(e) { message("  VR base design failed: ", conditionMessage(e)); NULL })
des_vr_adj  <- tryCatch(build_design(C_LOW, C_HIGH, pnl = panel_vr, conf = "vr", conf_in_model = TRUE),
                        error = function(e) { message("  VR adj design failed: ", conditionMessage(e)); NULL })
vr_row <- function(dsg, lab) {
  if (is.null(dsg)) return(tibble(model = lab, overall_rd = NA_real_, rd_discordant = NA_real_, gradient = NA_real_))
  v <- hte_from_design(dsg)
  tibble(model = lab, overall_rd = unname(v["All"]), rd_discordant = unname(v["Discordant"]),
         gradient = unname(v["Discordant"] - v["Concordant"]))
}
vr_hte <- bind_rows(vr_row(des, "primary (full cohort)"),
                    vr_row(des_vr_base, "VR subset, VR NOT in weights"),
                    vr_row(des_vr_adj,  "VR subset, VR IN weights"))
write_csv(vr_hte, file.path(final_dir, paste0("tte_ccw_disc_vr_hte_", site_name, ".csv")))
cat("--- 5a. HTE robustness to dead-space (VR) confounding ---\n")
print(as.data.frame(vr_hte %>% transmute(model, overall_pp = round(100 * overall_rd, 2),
        rd_disc_pp = round(100 * rd_discordant, 2), gradient_pp = round(100 * gradient, 2))), row.names = FALSE)
cat("    (Discordant RD/gradient stable VR-not-in vs VR-in => robust to dead-space confounding)\n")

# (b) does VR absorb the residual l_sf imbalance? l_sf weighted log-OR on the VR subset, base vs adj
sf_on_vrsubset <- function(dsg, lab) {
  if (is.null(dsg)) return(NULL)
  de <- panel_vr %>% select(hospitalization_id, vent_day, vtpfvc, sf, map, on_pressor, vr, disc_grp, death_day) %>%
    group_by(hospitalization_id) %>% arrange(vent_day) %>%
    mutate(viol = vent_day > GRACE & vtpfvc > C_LOW, prior_dev = lag(cumsum(viol), default = 0) > 0,
           l_sf = lag(sf), l_map = lag(map), l_pressor = lag(on_pressor), l_vr = lag(vr),
           alive = is.na(death_day) | vent_day <= death_day) %>% ungroup() %>%
    filter(!prior_dev, vent_day > GRACE, vent_day %in% 2:7, alive,
           !is.na(l_sf), !is.na(l_map), !is.na(l_pressor), !is.na(l_vr)) %>%
    left_join(dsg$bl$wday, by = c("hospitalization_id", "vent_day")) %>%
    mutate(cumw = trunc_w(coalesce(cumw, 1)), z_l_sf = as.numeric(scale(l_sf)),
           z_l_map = as.numeric(scale(l_map)), z_l_pressor = as.numeric(scale(l_pressor)))
  map_dfr(DISC_LEVELS, function(g) {
    d <- de %>% filter(as.character(disc_grp) == g); nd <- nrow(d); ndev <- sum(d$viol)
    if (nd < 30 || ndev < 5 || ndev > nd - 5)
      return(tibble(model = lab, disc_grp = g, n = nd, deviation_rate = round(ndev / max(nd, 1), 3),
                    sf_logOR = NA_real_, lo = NA_real_, hi = NA_real_))
    m  <- suppressWarnings(glm(viol ~ z_l_sf + z_l_map + z_l_pressor, data = d, family = binomial, weights = cumw))
    ct <- coeftest(m, vcov = sandwich::vcovCL(m, cluster = d$hospitalization_id))
    tibble(model = lab, disc_grp = g, n = nd, deviation_rate = round(ndev / nd, 3),
           sf_logOR = unname(ct["z_l_sf", "Estimate"]), lo = unname(ct["z_l_sf", "Estimate"] - 1.96 * ct["z_l_sf", "Std. Error"]),
           hi = unname(ct["z_l_sf", "Estimate"] + 1.96 * ct["z_l_sf", "Std. Error"]))
  })
}
vr_sfbal <- bind_rows(sf_on_vrsubset(des_vr_base, "VR NOT in weights"),
                      sf_on_vrsubset(des_vr_adj,  "VR IN weights")) %>%
  mutate(disc_grp = factor(disc_grp, DISC_LEVELS)) %>% arrange(disc_grp, model)
write_csv(vr_sfbal, file.path(final_dir, paste0("tte_ccw_disc_vr_sfbalance_", site_name, ".csv")))
cat("\n--- 5b. Residual l_sf log-OR on the VR subset: does adding VR to the weights shrink it? ---\n")
print(as.data.frame(vr_sfbal %>% transmute(disc_grp, model, n, deviation_rate,
        sf_logOR = round(sf_logOR, 3), ci = sprintf("[%.2f, %.2f]", lo, hi))), row.names = FALSE)
cat("    (smaller l_sf log-OR under 'VR IN weights' => the residual S/F was dead-space confounding that VR absorbs)\n")
message("Wrote tte_ccw_disc_vr_{hte,sfbalance}_", site_name, ".csv to ", final_dir)
