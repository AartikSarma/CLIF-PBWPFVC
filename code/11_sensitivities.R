# =============================================================================
# Script 11_sensitivities: consolidated TTE sensitivity sweeps (lead-site / supplement)
# =============================================================================
# NOT part of the per-site primary bundle (NOT in 11_run_all) -- run once (lead site or on
# pooled data) for the supplement. Rebuilds the design under alternative analytic choices and
# recomputes the 11.X primary HTE statistics (overall RD, per-tertile RDs, Discordant-
# Concordant gradient) via the SAME pooled, SOFA-adjusted, IPC-weighted standardized estimator
# as 11.X. POINT estimates only: each setting rebuilds the design + refits the deviation models
# and the outcome MSM; no per-setting bootstrap (the primary CIs live in 11.X).
#
# Sweeps:
#   1. STRAIN THRESHOLDS -- vary C_LOW (strain-limiting ceiling) and C_HIGH (permissive),
#      one-at-a-time around the ARMA-anchored primary (C_LOW=11 = LTV p75; C_HIGH=16 ~ HTV p25).
#      Does the discordance gradient survive other clinically-defensible ceilings?
#   2. CENSORING (trim x weight-cap) -- moved here from 11.X to keep the every-site leaf lean.
#
# The cohort (structural-positivity feasibility floor) is FIXED at the primary C_LOW, so each
# sweep isolates the ceiling/censoring choice on one comparable population; patients infeasible
# at a lower swept C_LOW simply deviate immediately (handled by the deviation/trim mechanism).
# Future: fold in the archived 11.C-W sensitivities (see code/SENSITIVITY_INVENTORY.md).
# =============================================================================
library(here)
source(here::here("code", "10_tte_engine.R"))

DISC_LEVELS <- c("Concordant", "Mid", "Discordant")
FORM  <- died ~ arm * disc_grp + arm * ns(day, 4) + disc_grp * ns(day, 4) + sofa_total
arm_f <- function() factor(c("permissive", "strain_limiting"), c("permissive", "strain_limiting"))

# standardized RD in tertile g (g="All" => overall), from a fitted pooled MSM -- identical to 11.X
std_rd <- function(fit, prof, g) {
  pp <- if (g == "All") prof else prof %>% filter(as.character(disc_grp) == g)
  if (nrow(pp) < 50) return(NA_real_)
  cells <- pp %>% count(disc_grp, sofa_total, name = "wt")
  grid  <- tidyr::crossing(cells, day = 1:HORIZON, arm = arm_f())
  grid$haz <- predict(fit, grid, type = "response")
  ci <- grid %>% group_by(arm, day) %>% summarise(h = weighted.mean(haz, wt), .groups = "drop") %>%
    group_by(arm) %>% arrange(day) %>% summarise(cif = 1 - prod(1 - h), .groups = "drop")
  ci$cif[ci$arm == "strain_limiting"] - ci$cif[ci$arm == "permissive"]
}

# one design -> overall RD + per-tertile RDs + Discordant-Concordant gradient. Fault-tolerant:
# an extreme ceiling can collapse the eligible deviation-model set (degenerate ns(vent_day)),
# esp. on the small synthetic cohort -> that setting degrades to NA instead of halting the sweep.
hte_for_design <- function(c_low, c_high, trim = TRIM_ALPHA, cap = DAYW_CAP, sf_term = "l_sf") {
  r <- tryCatch({
    dd <- build_design(c_low, c_high, GRACE, cap, "simple", trim = trim, sf_term = sf_term)
    lj <- dd$long %>% left_join(base %>% select(hospitalization_id, sofa_total), by = "hospitalization_id")
    pj <- lj %>% distinct(hospitalization_id, disc_grp, sofa_total)
    f  <- suppressWarnings(glm(FORM, data = lj, family = binomial, weights = ipcw))
    c(setNames(vapply(DISC_LEVELS, function(g) std_rd(f, pj, g), numeric(1)), DISC_LEVELS),
      All = std_rd(f, pj, "All"))
  }, error = function(e) { message("  setting C_LOW=", c_low, " C_HIGH=", c_high,
      " trim=", trim, " cap=", cap, " failed: ", conditionMessage(e));
      setNames(rep(NA_real_, 4), c(DISC_LEVELS, "All")) })
  tibble(c_low = c_low, c_high = c_high, trim = trim, cap = cap,
         overall_rd = unname(r["All"]), rd_concordant = unname(r["Concordant"]),
         rd_mid = unname(r["Mid"]), rd_discordant = unname(r["Discordant"]),
         gradient = unname(r["Discordant"] - r["Concordant"]))
}

# =============================================================================
# 1. Strain-threshold sweep: vary C_LOW and C_HIGH one-at-a-time around the primary
# =============================================================================
ceiling_grid <- c(
  list(c(C_LOW, C_HIGH)),                                  # PRIMARY (ARMA-anchored 11/16)
  lapply(c(9, 10, 12, 13), function(cl) c(cl, C_HIGH)),    # vary strain-limiting ceiling
  lapply(c(14, 18, 20),    function(ch) c(C_LOW, ch)))     # vary permissive ceiling
ceiling_grid <- Filter(function(p) p[1] < p[2], ceiling_grid)
message("11_sensitivities: strain-threshold sweep -- rebuilding ", length(ceiling_grid), " designs ...")
ceiling_sweep <- map_dfr(ceiling_grid, function(p) hte_for_design(p[1], p[2])) %>%
  mutate(setting = if_else(c_low == C_LOW & c_high == C_HIGH, "PRIMARY", "sensitivity"), .before = 1) %>%
  select(setting, c_low, c_high, overall_rd, rd_concordant, rd_mid, rd_discordant, gradient)
write_csv(ceiling_sweep, file.path(final_dir, paste0("tte_ccw_disc_ceiling_sweep_", site_name, ".csv")))

# =============================================================================
# 2. Censoring sweep: common-support trim x day-weight cap (at the primary ceilings)
# =============================================================================
cens_grid <- list(
  list(trim = TRIM_ALPHA, cap = DAYW_CAP),                 # PRIMARY
  list(trim = 0,          cap = DAYW_CAP),
  list(trim = 0.01,       cap = DAYW_CAP),
  list(trim = 0.05,       cap = DAYW_CAP),
  list(trim = TRIM_ALPHA, cap = 3),
  list(trim = TRIM_ALPHA, cap = 10))
message("11_sensitivities: censoring (trim x cap) sweep -- rebuilding ", length(cens_grid), " designs ...")
censoring_sweep <- map_dfr(cens_grid, function(s) hte_for_design(C_LOW, C_HIGH, trim = s$trim, cap = s$cap)) %>%
  mutate(setting = if_else(trim == TRIM_ALPHA & cap == DAYW_CAP, "PRIMARY", "sensitivity"), .before = 1) %>%
  select(setting, trim, cap, overall_rd, rd_concordant, rd_mid, rd_discordant, gradient)
write_csv(censoring_sweep, file.path(final_dir, paste0("tte_ccw_disc_censoring_sweep_", site_name, ".csv")))

# =============================================================================
# Console
# =============================================================================
cat("\n=== 1. Strain-threshold sweep (SOFA-adjusted standardized RDs, pp; PRIMARY = C_LOW", C_LOW,
    "/ C_HIGH", C_HIGH, ") ===\n")
print(as.data.frame(ceiling_sweep %>% transmute(setting, c_low, c_high,
        overall_pp = round(100 * overall_rd, 2), rd_conc_pp = round(100 * rd_concordant, 2),
        rd_disc_pp = round(100 * rd_discordant, 2), gradient_pp = round(100 * gradient, 2))), row.names = FALSE)
cat("    (gradient stays negative across ceilings = the discordance-HTE is not an artifact of the 11/16 choice)\n")

cat("\n=== 2. Censoring sweep: trim x weight-cap (point estimates) ===\n")
print(as.data.frame(censoring_sweep %>% transmute(setting, trim, cap,
        overall_pp = round(100 * overall_rd, 2), rd_disc_pp = round(100 * rd_discordant, 2),
        gradient_pp = round(100 * gradient, 2))), row.names = FALSE)
cat("    (gradient stable across trim/cap = not a thin-support / weight-tail artifact)\n")

# =============================================================================
# 3. Richer weight-model sensitivity: address the residual lagged-S/F imbalance (11.X TV-balance)
# =============================================================================
# 11.X's time-varying balance found lagged S/F still predicts deviation AFTER weighting in the
# Mid/Discordant strata -- the pooled LINEAR l_sf term under-corrects oxygenation where censoring
# is heaviest. Here the deviation/weight model's S/F term is enriched to
# ns(l_sf,3) + l_sf:disc_grp (nonlinear + stratum-specific slope), and we (a) re-run the same
# weighted-deviation-regression balance check on the richer weights and (b) recompute the HTE.
# Expectation: the weighted S/F log-OR/SD shrinks toward 0, and the Discordant RD is STABLE or
# slightly LARGER (the linear-model under-correction biased the strain arm sicker => toward null),
# confirming the headline is conservative, not inflated.
library(sandwich); library(lmtest)
RICH_SF <- "ns(l_sf, 3) + l_sf:disc_grp"

# weighted deviation ~ lagged-confounder S/F balance for one design (strain arm, days 2-7) -- as 11.X D2
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
message("11_sensitivities: richer weight-model (S/F) sensitivity -- rebuilding rich design ...")
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

# =============================================================================
# 4. Weight-timing diagnostic: is the S/F "imbalance" an artifact of weighting the balance
#    check by the cumulative IPCW THROUGH the current day?
# =============================================================================
# 11.X D2 weighted each eligible day t by cumw(t) -- which includes day t's OWN deviation factor
# (1-p_den(t)), large precisely for high-l_sf DEVIATORS (high oxygenation -> high modeled
# deviation prob -> small denominator -> big weight). So cumw(t) up-weights high-l_sf deviators
# and can manufacture a positive viol~l_sf association (note unweighted->weighted INCREASED at
# UCSF). The principled weight is the PRIOR-day cumulative weight cumw(t-1) (history through t-1,
# independent of today's decision). Here we recompute the strain-arm S/F balance (primary design)
# under both. If prior-day collapses Mid/Discordant toward 0, the imbalance was a check artifact
# (not residual confounding) and D2 should switch to the prior-day weight.
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
message("Wrote tte_ccw_disc_{ceiling,censoring}_sweep_ + sf_weightmodel_{balance,hte}_ + sf_weighttiming_", site_name, ".csv to ", final_dir)
