# =============================================================================
# Script 36.B: Diagnostics -- positivity / weights / balance (§10g, 10g2, 10g3, 10g4)
#              + the two write-free objects computed in common (panel drops, struct excl)
# =============================================================================
# Sources the shared common file (30_tte_common.R). Writes the diagnostic CSVs/PDF and
# the two CSVs whose computation lives in common (panel_drop_summary, struct_excl).
# =============================================================================
library(here)
source(here::here("code", "31_tte_engine.R"))

# --- CSV writes pulled out of common (objects computed in 10b / 10a) ----------
# [T11] missing-covariate drop diagnostic (panel_drop_summary built in common §10b)
write_csv(panel_drop_summary, file.path(final_dir, paste0("tte_ccw_panel_drops_", site_name, ".csv")))
# [T9b] structural-positivity exclusion (struct_excl built in common §10a, pre-restriction)
write_csv(struct_excl, file.path(final_dir, paste0("tte_ccw_structural_excluded_", site_name, ".csv")))

# =============================================================================
# 10g. Diagnostics: deviation, weights, per-arm positivity ([T7])
# =============================================================================
diag <- bind_rows(des$bl$idsum %>% mutate(arm = "strain_limiting"),
                  des$bh$idsum %>% mutate(arm = "permissive")) %>%
  group_by(arm) %>%
  summarise(n_patients = n(), frac_deviated = mean(is.finite(dev_day)),
            wt_p99 = quantile(trunc_w(ipcw_term), 0.99), wt_max = max(ipcw_term),
            ess_frac = ess_frac(ipcw_term), .groups = "drop")
write_csv(diag, file.path(final_dir, paste0("tte_ccw_diagnostics_", site_name, ".csv")))
cat("\n=== CCW diagnostics (per-arm positivity / weights) ===\n")
print(as.data.frame(diag %>% mutate(across(where(is.numeric), ~ round(., 3)))), row.names = FALSE)

# =============================================================================
# 10g2. Positivity scan ([T7]): can the strain-limiting target even be reached,
#       and is it reached often enough for IPCW to stand on? Two complementary
#       reads, both split by AGE TERTILE (older -> lower PFVC -> the cell where
#       deviation concentrates and weights blow up):
#   structural -- at the LTVV VT floor (VT_FLOOR_MLKG x PBW), is the implied
#                 VT/PFVC already above the ceiling? If so adherence is PHYSICALLY
#                 impossible: P(adhere)=0 by construction, not a small probability,
#                 and no weight repairs a stratum with zero adherers.
#   empirical  -- distribution of the weight model's predicted P(adhere)=1-p_den on
#                 eligible days. A mass near 0 is a near positivity violation: the
#                 rare adherers get enormous weights and dominate that arm/age cell.
# A non-trivial structural-infeasible fraction in the oldest tertile is the signal
# that the static <=C_LOW estimand is not identified there, motivating the feasible-
# shift (modified treatment policy) reframing rather than leaning harder on weights.
# =============================================================================
ceil_by_arm <- c(strain_limiting = C_LOW, permissive = C_HIGH)
struct <- base %>% mutate(min_vtpfvc = VT_FLOOR_MLKG * pbw / pfvc * 0.1)
pos_structural <- imap_dfr(ceil_by_arm, function(ceil, arm_lab) {
  by_age <- struct %>% group_by(age_grp) %>%
    summarise(n_patients = n(), frac_infeasible = mean(min_vtpfvc > ceil),
              median_min_vtpfvc = median(min_vtpfvc), .groups = "drop") %>%
    mutate(age_grp = as.character(age_grp))
  overall <- struct %>%
    summarise(age_grp = "All", n_patients = n(),
              frac_infeasible = mean(min_vtpfvc > ceil),
              median_min_vtpfvc = median(min_vtpfvc))
  bind_rows(by_age, overall) %>% mutate(arm = arm_lab, ceiling = ceil)
}) %>% select(arm, ceiling, age_grp, n_patients, frac_infeasible, median_min_vtpfvc)
write_csv(pos_structural, file.path(final_dir, paste0("tte_ccw_positivity_structural_", site_name, ".csv")))

pday <- bind_rows(des$bl$pday %>% mutate(arm = "strain_limiting"),
                  des$bh$pday %>% mutate(arm = "permissive"))
emp_one <- function(g) summarise(g,
  n_eligible_days = n(), frac_padhere_lt05 = mean(p_adhere < 0.05),
  frac_padhere_lt02 = mean(p_adhere < 0.02), p01_padhere = quantile(p_adhere, 0.01),
  median_padhere = median(p_adhere), .groups = "drop")
pos_empirical <- bind_rows(
  pday %>% group_by(arm, age_grp) %>% emp_one() %>% mutate(age_grp = as.character(age_grp)),
  pday %>% group_by(arm) %>% emp_one() %>% mutate(age_grp = "All")) %>%
  arrange(arm, age_grp)
write_csv(pos_empirical, file.path(final_dir, paste0("tte_ccw_positivity_empirical_", site_name, ".csv")))
# Same empirical positivity scan, but by PBW/PFVC DISCORDANCE tertile -- the strain-limiting
# contrast lives in the Discordant band, so this shows whether positivity is thinnest exactly
# where the effect concentrates (the TTE analog of the MP discordance work).
pos_empirical_disc <- bind_rows(
  pday %>% group_by(arm, disc_grp) %>% emp_one() %>% mutate(disc_grp = as.character(disc_grp)),
  pday %>% group_by(arm) %>% emp_one() %>% mutate(disc_grp = "All")) %>% arrange(arm, disc_grp)
write_csv(pos_empirical_disc, file.path(final_dir, paste0("tte_ccw_positivity_empirical_disc_", site_name, ".csv")))
cat("\n=== positivity: structural (infeasible at ", VT_FLOOR_MLKG, " mL/kg PBW floor) ===\n", sep = "")
print(as.data.frame(pos_structural %>% mutate(across(where(is.numeric), ~ round(., 3)))), row.names = FALSE)
cat("\n=== positivity: empirical (predicted P(adhere) on eligible days) ===\n")
print(as.data.frame(pos_empirical %>% mutate(across(where(is.numeric), ~ round(., 3)))), row.names = FALSE)

# =============================================================================
# 10g3. Weighted covariate balance over follow-up ([T7]): the figure that shows
#       the weights doing their job. At each horizon day, take the clones still
#       uncensored-by-deviation and still in follow-up, and compare their BASELINE
#       covariates to the full arm cohort as a standardized mean difference (SMD),
#       both unweighted and IPCW-weighted. Informative censoring makes the
#       unweighted at-risk set drift (the strain arm loses its older/lower-PFVC
#       clones, so its surviving set gets younger); if IPCW is doing its job the
#       weighted SMD stays flat near 0 while the unweighted drifts past +/-0.1.
#       NOTE: this uses the PRIMARY (time-only-numerator) weights, which target V
#       balance -- so weighted SMD -> 0 is the right expectation here. (The old
#       baseline-numerator stabilized weights left V in by design, so weighted ~=
#       unweighted under them; that is why V had to move out of the numerator. If
#       even these weights cannot pull V's SMD to 0, the cap is binding and residual
#       confounding remains -- a real finding, not a stabilization artifact.)
# =============================================================================
# Early days included deliberately: most of the LTVV effect is generated in the first week
# (the SOFA trajectory separates at day 2), and early is also where survivorship contamination
# is smallest -- so the early balance is both the more relevant and the cleaner read. The late
# days are kept to show how survivorship and size-selection grow over follow-up.
bal_days <- c(1L, 2L, 3L, 5L, 7L, 14L, 21L, 28L)
bal_covs <- c(age10 = "Age (per 10y)", pfvc = "PFVC (L)", sofa_total = "SOFA")
balance_for_arm <- function(b, arm_lab) {
  bb  <- b$idsum %>% left_join(base %>% select(hospitalization_id, pfvc, sofa_total),
                               by = "hospitalization_id")
  ref <- bb %>% summarise(across(all_of(names(bal_covs)),
                                 list(m = ~mean(.), s = ~sd(.)), .names = "{.col}__{.fn}"))
  map_dfr(bal_days, function(d) {
    atrisk <- bb %>% filter(dev_day > d, trim_day > d, is.na(death_day) | death_day >= d)
    if (nrow(atrisk) < 20) return(NULL)
    w_at_d <- b$wday %>% filter(vent_day <= d) %>%
      group_by(hospitalization_id) %>% arrange(vent_day) %>%
      summarise(w = last(cumw), .groups = "drop")
    atrisk <- atrisk %>% left_join(w_at_d, by = "hospitalization_id") %>%
      mutate(w = trunc_w(coalesce(w, 1)))
    map_dfr(names(bal_covs), function(cv) {
      x <- atrisk[[cv]]
      tibble(arm = arm_lab, day = d, covariate = unname(bal_covs[[cv]]),
             smd_unweighted = (mean(x) - ref[[paste0(cv, "__m")]]) / ref[[paste0(cv, "__s")]],
             smd_weighted   = (weighted.mean(x, atrisk$w) - ref[[paste0(cv, "__m")]]) /
                              ref[[paste0(cv, "__s")]],
             n_atrisk   = nrow(atrisk),
             ess_atrisk = (sum(atrisk$w)^2 / sum(atrisk$w^2)))
    })
  })
}
balance <- bind_rows(balance_for_arm(des$bl, "strain_limiting"),
                     balance_for_arm(des$bh, "permissive"))
write_csv(balance, file.path(final_dir, paste0("tte_ccw_balance_", site_name, ".csv")))
bal_long <- balance %>%
  pivot_longer(c(smd_unweighted, smd_weighted), names_to = "weighting",
               values_to = "smd", names_prefix = "smd_")
pb <- ggplot(bal_long, aes(day, smd, colour = weighting, linetype = arm, shape = arm)) +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dotted", colour = "grey60") +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_line() + geom_point(size = 2) +
  facet_wrap(~ covariate, scales = "free_y") +
  scale_colour_manual(values = c(unweighted = okabe[5], weighted = okabe[4]), name = NULL) +
  scale_linetype_discrete(name = NULL) + scale_shape_discrete(name = NULL) +
  labs(x = "Days from index ventilation", y = "SMD vs full-cohort baseline",
       title = "IPCW covariate balance over follow-up (at-risk clones vs baseline)",
       subtitle = paste0(site_name, if (is_synthetic) " (SYNTHETIC - plumbing only)" else "",
         " - weighted SMD near 0 = informative censoring corrected; dotted = +/-0.1")) +
  theme_minimal(base_size = 10) + theme(legend.position = "top")
ggsave(file.path(final_dir, paste0("tte_ccw_balance_", site_name, ".pdf")), pb, width = 9, height = 4)
cat("\n=== IPCW balance over follow-up (SMD vs baseline; |SMD|<0.1 = balanced) ===\n")
print(as.data.frame(balance %>% mutate(across(where(is.numeric), ~ round(., 3)))), row.names = FALSE)

# =============================================================================
# 10g4. Weights by AGE TERTILE x arm ([T7]): the per-arm diagnostic above hides
#       that the damage is local. ESS collapse and extreme weights concentrate in
#       the oldest tertile of the strain-limiting arm; mean stabilized weight should
#       sit near 1 (a systematic departure flags model misspecification or a
#       positivity problem). A near-zero ESS here means the old-subgroup RD rests on
#       a handful of effective patients -- say so before a reviewer does.
# =============================================================================
diag_age <- bind_rows(des$bl$idsum %>% mutate(arm = "strain_limiting"),
                      des$bh$idsum %>% mutate(arm = "permissive")) %>%
  group_by(arm, age_grp) %>%
  summarise(n_patients = n(), frac_deviated = mean(is.finite(dev_day)),
            mean_stab_w = mean(ipcw_term),
            wt_p99 = quantile(trunc_w(ipcw_term), 0.99), wt_max = max(ipcw_term),
            ess_frac = ess_frac(ipcw_term), .groups = "drop")
write_csv(diag_age, file.path(final_dir, paste0("tte_ccw_weights_by_age_", site_name, ".csv")))
cat("\n=== weights by age tertile x arm (mean stab. wt ~1; watch ESS in Old x strain) ===\n")
print(as.data.frame(diag_age %>% mutate(across(where(is.numeric), ~ round(., 3)))), row.names = FALSE)

# Same, by PBW/PFVC discordance tertile: the strain arm's effective sample should collapse in
# the Discordant band (small lungs structurally cannot reach the ceiling) -- the positivity
# caveat made visible on the thesis's own axis, alongside the discordance subgroup RD (35_tte_primary).
diag_disc <- bind_rows(des$bl$idsum %>% mutate(arm = "strain_limiting"),
                       des$bh$idsum %>% mutate(arm = "permissive")) %>%
  group_by(arm, disc_grp) %>%
  summarise(n_patients = n(), frac_deviated = mean(is.finite(dev_day)),
            mean_stab_w = mean(ipcw_term),
            wt_p99 = quantile(trunc_w(ipcw_term), 0.99), wt_max = max(ipcw_term),
            ess_frac = ess_frac(ipcw_term), .groups = "drop")
write_csv(diag_disc, file.path(final_dir, paste0("tte_ccw_weights_by_disc_", site_name, ".csv")))
cat("\n=== weights by PBW/PFVC discordance tertile x arm (watch ESS in Discordant x strain) ===\n")
print(as.data.frame(diag_disc %>% mutate(across(where(is.numeric), ~ round(., 3)))), row.names = FALSE)
