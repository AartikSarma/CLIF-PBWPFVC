# =============================================================================
# 11_sens_severity_ladder: time-varying severity-confounder ladder (gradient robustness)
# =============================================================================
# The sharpest gradient worry: time-varying severity drives BOTH adherence (deviation/censoring)
# AND death, DIFFERENTIALLY across discordance, so an IPCW fit on a respiratory/cardiovascular core
# (l_sf,l_fio2,l_peep,l_rr,l_map,l_pressor) + BASELINE sofa could leave a gradient that is really
# severity. Each rung adds ONE lagged time-varying severity confounder to the WEIGHT model (engine
# conf hook) on ITS recorded subset -- the g-methods-correct place for a tv confounder -- and reports
# the discordance gradient with that confounder OUT vs IN the weights on the SAME day-set. We do NOT
# stack them (each recorded subset differs; intersecting collapses the sample); the claim is that the
# gradient is stable to adding ANY of them. Outputs tte_ccw_disc_severity_ladder_<site>.csv.
# =============================================================================
source(here::here("code", "11_sens_common.R"))

cat("\n=== 6. Time-varying severity-confounder ladder (discordance gradient OUT vs IN the IPCW) ===\n")
ladder_row <- function(d_out, d_in, label, n) {
  one <- function(d) if (is.null(d)) c(NA_real_, NA_real_) else {
    v <- hte_from_design(d); c(unname(v["Discordant"] - v["Concordant"]), unname(v["All"])) }
  o <- one(d_out); i <- one(d_in)
  tibble(rung = label, subset_n = n, gradient_conf_out = o[1], gradient_conf_in = i[1],
         overall_rd_in = i[2], change = i[1] - o[1])
}
# rung: daily SOFA (CNS/renal/hepatic/coag organs the core panel misses) -- realign 03's
# admission-indexed sofa_daily to vent_day (the 11.Q construction), LOCF, lag into the weights.
sofa_ladder <- tryCatch({
  adm    <- read_parquet(file.path(output_dir, "cohort_demographics.parquet")) %>%
    select(hospitalization_id, admission_dttm)
  offset <- base %>% left_join(adm, by = "hospitalization_id") %>%
    transmute(hospitalization_id, index_day = floor(as.numeric(difftime(t0, admission_dttm, units = "days"))))
  sofa_d <- read_parquet(file.path(output_dir, "sofa_daily.parquet")) %>%
    select(hospitalization_id, sofa_day, sofa_total) %>%
    inner_join(offset, by = "hospitalization_id") %>%
    transmute(hospitalization_id, vent_day = sofa_day - index_day, sofa = sofa_total) %>%
    filter(vent_day >= 0, vent_day <= MAX_VENT_DAY)
  panel_sofa <- panel %>% filter(hospitalization_id %in% unique(sofa_d$hospitalization_id)) %>%
    left_join(sofa_d, by = c("hospitalization_id", "vent_day")) %>%
    group_by(hospitalization_id) %>% arrange(vent_day) %>% fill(sofa, .direction = "down") %>% ungroup()
  n_s <- n_distinct(panel_sofa %>% filter(!is.na(sofa)) %>% pull(hospitalization_id))
  if (n_s < 100) NULL else list(n = n_s,
    base = build_design(C_LOW, C_HIGH, pnl = panel_sofa, conf = "sofa", conf_in_model = FALSE),
    adj  = build_design(C_LOW, C_HIGH, pnl = panel_sofa, conf = "sofa", conf_in_model = TRUE))
}, error = function(e) { message("  daily-SOFA rung failed: ", conditionMessage(e)); NULL })
# rung: ventilatory ratio (dead space) -- shared VR panel from common; build its own designs.
vr_ladder <- tryCatch({
  vrp <- make_vr_panel()
  list(n = n_distinct(vrp$vr_daily$hospitalization_id),
       base = build_design(C_LOW, C_HIGH, pnl = vrp$panel_vr, conf = "vr", conf_in_model = FALSE),
       adj  = build_design(C_LOW, C_HIGH, pnl = vrp$panel_vr, conf = "vr", conf_in_model = TRUE))
}, error = function(e) { message("  VR rung failed: ", conditionMessage(e)); NULL })
# rung: driving pressure (plateau-recorded only, NOT forward-filled) -- engine dp_daily.
dp_ladder <- tryCatch({
  panel_dp <- panel %>% filter(hospitalization_id %in% unique(dp_daily$hospitalization_id)) %>%
    left_join(dp_daily, by = c("hospitalization_id", "vent_day"))
  n_d <- n_distinct(panel_dp %>% filter(!is.na(dp)) %>% pull(hospitalization_id))
  if (n_d < 100) NULL else list(n = n_d,
    base = build_design(C_LOW, C_HIGH, pnl = panel_dp, conf = "dp", conf_in_model = FALSE),
    adj  = build_design(C_LOW, C_HIGH, pnl = panel_dp, conf = "dp", conf_in_model = TRUE))
}, error = function(e) { message("  driving-pressure rung failed: ", conditionMessage(e)); NULL })

sev_ladder <- bind_rows(
  ladder_row(des, des, "baseline SOFA only (primary, full cohort)", length(ids)),
  if (!is.null(sofa_ladder)) ladder_row(sofa_ladder$base, sofa_ladder$adj, "+ daily SOFA (CNS/renal/hepatic/coag)", sofa_ladder$n),
  if (!is.null(vr_ladder))   ladder_row(vr_ladder$base,   vr_ladder$adj,   "+ ventilatory ratio (dead space)", vr_ladder$n),
  if (!is.null(dp_ladder))   ladder_row(dp_ladder$base,   dp_ladder$adj,   "+ driving pressure (plateau-recorded)", dp_ladder$n))
write_csv(sev_ladder, file.path(final_dir, paste0("tte_ccw_disc_severity_ladder_", site_name, ".csv")))
print(as.data.frame(sev_ladder %>% mutate(across(c(gradient_conf_out, gradient_conf_in, overall_rd_in, change),
                                                 ~ round(100 * ., 2)))), row.names = FALSE)
cat("    (gradient_conf_in ~ gradient_conf_out at every rung => the discordance gradient is NOT generated by\n")
cat("     time-varying severity; small |change|, same sign = stable to g-methods-correct tv-severity control)\n")
message("Wrote tte_ccw_disc_severity_ladder_", site_name, ".csv to ", final_dir)
