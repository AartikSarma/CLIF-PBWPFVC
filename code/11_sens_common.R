# =============================================================================
# 11_sens_common: shared substrate for the split TTE sensitivity pieces
# =============================================================================
# Sourced by each 11_sens_*.R topic file (and by the 11_sensitivities.R run-all wrapper) so a new
# test can be run on its own WITHOUT rebuilding every other sweep. Defines the shared discordance-
# gradient estimator (identical to 11.X) + a lazy VR-panel builder. Sourcing the engine (the
# expensive primary-design build) is GUARDED: it runs once per session, so re-sourcing common from
# several pieces (e.g. via the run-all wrapper) is cheap.
#
# Each piece writes the SAME output filenames as the old monolithic 11_sensitivities.R, so the
# report (figures/make_tte_report.py) and pooling are unaffected by the split.
# =============================================================================
library(here)
if (!exists("des") || !exists("long_all") || !exists("base"))
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

# overall + per-tertile standardized RDs from a PRE-BUILT design (the shared estimator).
hte_from_design <- function(dsg) {
  lj <- dsg$long %>% left_join(base %>% select(hospitalization_id, sofa_total), by = "hospitalization_id")
  pj <- lj %>% distinct(hospitalization_id, disc_grp, sofa_total)
  f  <- suppressWarnings(glm(FORM, data = lj, family = binomial, weights = ipcw))
  c(All = std_rd(f, pj, "All"),
    setNames(vapply(DISC_LEVELS, function(g) std_rd(f, pj, g), numeric(1)), DISC_LEVELS))
}
# one design (built from ceiling/trim/cap/sf_term) -> RD row. Fault-tolerant: an extreme ceiling
# can collapse the eligible deviation-model set (degenerate ns(vent_day)) on the small synthetic
# cohort -> that setting degrades to NA instead of halting the sweep.
hte_for_design <- function(c_low, c_high, trim = TRIM_ALPHA, cap = DAYW_CAP, sf_term = "l_sf") {
  r <- tryCatch(
    hte_from_design(build_design(c_low, c_high, GRACE, cap, "simple", trim = trim, sf_term = sf_term)),
    error = function(e) { message("  setting C_LOW=", c_low, " C_HIGH=", c_high, " trim=", trim,
      " cap=", cap, " failed: ", conditionMessage(e)); setNames(rep(NA_real_, 4), c("All", DISC_LEVELS)) })
  tibble(c_low = c_low, c_high = c_high, trim = trim, cap = cap,
         overall_rd = unname(r["All"]), rd_concordant = unname(r["Concordant"]),
         rd_mid = unname(r["Mid"]), rd_discordant = unname(r["Discordant"]),
         gradient = unname(r["Discordant"] - r["Concordant"]))
}

# Lazy VR-panel builder (the gas-subset dead-space confounder), used by 11_sens_deadspace and the
# severity ladder. VR = (VE x PaCO2) / (PBW x 100 x 37.5). Returns daily VR + the VR-joined panel;
# callers build their own designs (conf = "vr") from panel_vr so each piece stays self-contained.
make_vr_panel <- function() {
  ve_daily <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet")) %>%
    select(hospitalization_id, recorded_dttm, tidal_volume_set, resp_rate_set) %>%
    filter(!is.na(tidal_volume_set), tidal_volume_set > 0, !is.na(resp_rate_set), resp_rate_set > 0) %>%
    inner_join(base %>% select(hospitalization_id, t0, pbw), by = "hospitalization_id") %>%
    mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
    filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
    group_by(hospitalization_id, vent_day, pbw) %>%
    summarise(ve = median(tidal_volume_set * resp_rate_set), .groups = "drop")   # mL/min
  paco2_daily <- read_parquet(file.path(output_dir, "cohort_labs_clean.parquet")) %>%
    filter(lab_category == "pco2_arterial", !is.na(lab_value_numeric)) %>%
    inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
    mutate(vent_day = floor(as.numeric(difftime(lab_result_dttm, t0, units = "days")))) %>%
    filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
    group_by(hospitalization_id, vent_day) %>% summarise(paco2 = max(lab_value_numeric), .groups = "drop")
  vr_daily <- ve_daily %>% inner_join(paco2_daily, by = c("hospitalization_id", "vent_day")) %>%
    mutate(vr = (ve * paco2) / (pbw * 100 * 37.5)) %>%
    filter(is.finite(vr), between(vr, 0.2, 10)) %>% select(hospitalization_id, vent_day, vr)
  list(vr_daily = vr_daily,
       panel_vr = panel %>% left_join(vr_daily, by = c("hospitalization_id", "vent_day")))
}
