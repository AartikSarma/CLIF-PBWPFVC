# =============================================================================
# 11_sens_weightmodel: S/F weight-model robustness (richer deviation model + weight-timing)
# =============================================================================
# Two reads on the residual lagged-S/F imbalance 11.X's time-varying balance flagged:
#   3. RICHER deviation/weight model -- enrich the S/F term to ns(l_sf,3)+l_sf:disc_grp and check
#      (a) the weighted-deviation S/F balance and (b) the HTE move (stable/larger Discordant RD =
#      the linear-model under-correction was conservative).
#   4. WEIGHT-TIMING -- is the "imbalance" an artifact of weighting the balance check by cumw(t)
#      (which includes today's own deviation factor) vs the principled prior-day cumw(t-1)?
# Point estimates only. Outputs tte_ccw_disc_sf_weightmodel_{balance,hte}_ + sf_weighttiming_<site>.
# =============================================================================
source(here::here("code", "11_sens_common.R"))
library(sandwich); library(lmtest)

# --- 3. Richer weight-model sensitivity --------------------------------------------------------
# 11.X's time-varying balance found lagged S/F still predicts deviation AFTER weighting in the
# Mid/Discordant strata -- the pooled LINEAR l_sf term under-corrects oxygenation where censoring
# is heaviest. Enrich the deviation/weight model's S/F term to ns(l_sf,3) + l_sf:disc_grp and
# (a) re-run the weighted-deviation balance check, (b) recompute the HTE. Expectation: the weighted
# S/F log-OR/SD shrinks toward 0, and the Discordant RD is STABLE or slightly LARGER (linear
# under-correction biased the strain arm sicker => toward null), confirming the headline is
# conservative, not inflated.
RICH_SF <- "ns(l_sf, 3) + l_sf:disc_grp"

# weighted deviation ~ lagged-confounder S/F balance for one design (strain arm, days 2-7) -- 11.X D2
tv_sf_check <- function(dsg, label) {
  de <- panel %>% select(hospitalization_id, vent_day, vtpfvc, sf, map, on_pressor, disc_grp, death_day) %>%
    group_by(hospitalization_id) %>% arrange(vent_day) %>%
    mutate(viol = vent_day > GRACE & vtpfvc > C_LOW, prior_dev = lag(cumsum(viol), default = 0) > 0,
           l_sf = lag(sf), l_map = lag(map), l_pressor = lag(on_pressor),
           alive = is.na(death_day) | vent_day <= death_day) %>% ungroup() %>%
    filter(!prior_dev, vent_day > GRACE, vent_day %in% 2:7, alive,
           !is.na(l_sf), !is.na(l_map), !is.na(l_pressor)) %>%
    left_join(dsg$bl$wday, by = c("hospitalization_id", "vent_day")) %>%
    mutate(cumw = trunc_w(coalesce(cumw, 1)), z_l_sf = as.numeric(scale(l_sf)),
           z_l_map = as.numeric(scale(l_map)), z_l_pressor = as.numeric(scale(l_pressor)))
  map_dfr(DISC_LEVELS, function(g) {
    d <- de %>% filter(as.character(disc_grp) == g); nd <- nrow(d); ndev <- sum(d$viol)
    if (nd < 50 || ndev < 10 || ndev > nd - 10)
      return(tibble(model = label, disc_grp = g, deviation_rate = round(ndev / nd, 3),
                    sf_wt_logOR = NA_real_, lo = NA_real_, hi = NA_real_))
    m  <- suppressWarnings(glm(viol ~ z_l_sf + z_l_map + z_l_pressor, data = d, family = binomial, weights = cumw))
    ct <- coeftest(m, vcov = sandwich::vcovCL(m, cluster = d$hospitalization_id))
    tibble(model = label, disc_grp = g, deviation_rate = round(ndev / nd, 3),
           sf_wt_logOR = unname(ct["z_l_sf", "Estimate"]),
           lo = unname(ct["z_l_sf", "Estimate"] - 1.96 * ct["z_l_sf", "Std. Error"]),
           hi = unname(ct["z_l_sf", "Estimate"] + 1.96 * ct["z_l_sf", "Std. Error"]))
  })
}
message("11_sens_weightmodel: richer weight-model (S/F) sensitivity -- rebuilding rich design ...")
des_rich <- tryCatch(build_design(C_LOW, C_HIGH, GRACE, DAYW_CAP, "simple", sf_term = RICH_SF),
                     error = function(e) { message("  rich design failed: ", conditionMessage(e)); NULL })
sf_balance <- bind_rows(tv_sf_check(des, "primary (linear l_sf)"),
                        if (!is.null(des_rich)) tv_sf_check(des_rich, "rich (ns + l_sf:disc_grp)")) %>%
  mutate(disc_grp = factor(disc_grp, DISC_LEVELS)) %>% arrange(disc_grp, model)
write_csv(sf_balance, file.path(final_dir, paste0("tte_ccw_disc_sf_weightmodel_balance_", site_name, ".csv")))
sf_hte <- bind_rows(hte_for_design(C_LOW, C_HIGH) %>% mutate(model = "primary (linear l_sf)", .before = 1),
                    hte_for_design(C_LOW, C_HIGH, sf_term = RICH_SF) %>% mutate(model = "rich (ns + l_sf:disc_grp)", .before = 1)) %>%
  select(model, overall_rd, rd_concordant, rd_mid, rd_discordant, gradient)
write_csv(sf_hte, file.path(final_dir, paste0("tte_ccw_disc_sf_weightmodel_hte_", site_name, ".csv")))

cat("\n=== 3a. Residual lagged-S/F imbalance: linear vs richer weight model (strain arm, days 2-7) ===\n")
print(as.data.frame(sf_balance %>% transmute(disc_grp, model, deviation_rate,
        sf_wt_logOR = round(sf_wt_logOR, 3), ci = sprintf("[%.2f, %.2f]", lo, hi))), row.names = FALSE)
cat("    (richer model should pull the Discordant/Mid weighted S/F log-OR toward 0 = oxygenation now balanced)\n")
cat("\n=== 3b. HTE under linear vs richer weight model (does the Discordant RD move?) ===\n")
print(as.data.frame(sf_hte %>% transmute(model, overall_pp = round(100 * overall_rd, 2),
        rd_disc_pp = round(100 * rd_discordant, 2), gradient_pp = round(100 * gradient, 2))), row.names = FALSE)
cat("    (Discordant RD stable or LARGER under the richer model = the linear-model imbalance was conservative)\n")

# --- 4. Weight-timing diagnostic ---------------------------------------------------------------
# 11.X D2 weighted each eligible day t by cumw(t) -- which includes day t's OWN deviation factor
# (1-p_den(t)), large precisely for high-l_sf DEVIATORS (high oxygenation -> high modeled deviation
# prob -> small denominator -> big weight). So cumw(t) up-weights high-l_sf deviators and can
# manufacture a positive viol~l_sf association (unweighted->weighted INCREASED at UCSF). The
# principled weight is the PRIOR-day cumulative weight cumw(t-1) (history through t-1, independent
# of today's decision). Recompute the strain-arm S/F balance (primary design) under both. If
# prior-day collapses Mid/Discordant toward 0, the imbalance was a check artifact, not confounding.
wlag <- des$bl$wday %>% group_by(hospitalization_id) %>% arrange(vent_day) %>%
  mutate(cumw_prior = lag(cumw, default = 1)) %>% ungroup()
de_t <- panel %>% select(hospitalization_id, vent_day, vtpfvc, sf, map, on_pressor, disc_grp, death_day) %>%
  group_by(hospitalization_id) %>% arrange(vent_day) %>%
  mutate(viol = vent_day > GRACE & vtpfvc > C_LOW, prior_dev = lag(cumsum(viol), default = 0) > 0,
         l_sf = lag(sf), l_map = lag(map), l_pressor = lag(on_pressor),
         alive = is.na(death_day) | vent_day <= death_day) %>% ungroup() %>%
  filter(!prior_dev, vent_day > GRACE, vent_day %in% 2:7, alive, !is.na(l_sf), !is.na(l_map), !is.na(l_pressor)) %>%
  left_join(wlag %>% select(hospitalization_id, vent_day, cumw, cumw_prior), by = c("hospitalization_id", "vent_day")) %>%
  mutate(w_current = trunc_w(coalesce(cumw, 1)), w_prior = trunc_w(coalesce(cumw_prior, 1)),
         z_l_sf = as.numeric(scale(l_sf)), z_l_map = as.numeric(scale(l_map)), z_l_pressor = as.numeric(scale(l_pressor)))
sf_logor <- function(d, wv) {
  m  <- suppressWarnings(glm(viol ~ z_l_sf + z_l_map + z_l_pressor, data = d, family = binomial, weights = wv))
  ct <- coeftest(m, vcov = sandwich::vcovCL(m, cluster = d$hospitalization_id))
  c(est = unname(ct["z_l_sf", "Estimate"]), se = unname(ct["z_l_sf", "Std. Error"]))
}
sf_timing <- map_dfr(DISC_LEVELS, function(g) {
  d <- de_t %>% filter(as.character(disc_grp) == g); nd <- nrow(d); ndev <- sum(d$viol)
  if (nd < 50 || ndev < 10 || ndev > nd - 10)
    return(tibble(disc_grp = g, deviation_rate = round(ndev / nd, 3), sf_current = NA_real_,
                  lo_cur = NA_real_, hi_cur = NA_real_, sf_prior = NA_real_, lo_pri = NA_real_, hi_pri = NA_real_))
  cw <- sf_logor(d, d$w_current); pw <- sf_logor(d, d$w_prior)
  tibble(disc_grp = g, deviation_rate = round(ndev / nd, 3),
         sf_current = cw["est"], lo_cur = cw["est"] - 1.96 * cw["se"], hi_cur = cw["est"] + 1.96 * cw["se"],
         sf_prior  = pw["est"], lo_pri = pw["est"] - 1.96 * pw["se"], hi_pri = pw["est"] + 1.96 * pw["se"])
}) %>% mutate(disc_grp = factor(disc_grp, DISC_LEVELS)) %>% arrange(disc_grp)
write_csv(sf_timing, file.path(final_dir, paste0("tte_ccw_disc_sf_weighttiming_", site_name, ".csv")))
cat("\n=== 4. S/F balance: current-day vs prior-day IPCW weighting (primary design, strain arm) ===\n")
print(as.data.frame(sf_timing %>% transmute(disc_grp, deviation_rate,
        current = sprintf("%.2f [%.2f, %.2f]", sf_current, lo_cur, hi_cur),
        prior   = sprintf("%.2f [%.2f, %.2f]", sf_prior, lo_pri, hi_pri))), row.names = FALSE)
cat("    (prior-day collapsing Mid/Discordant toward 0 [CI incl 0] => the imbalance was a weight-timing artifact, not residual confounding)\n")
message("Wrote tte_ccw_disc_sf_weightmodel_{balance,hte}_ + sf_weighttiming_", site_name, ".csv to ", final_dir)
