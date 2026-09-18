# =============================================================================
# 33_tte_ceiling_common: shared substrate for the normalizer head-to-head TTE pieces (33_tte_ceiling_*.R)
# =============================================================================
# The target trial that asks the paper's question directly: does dosing to PREDICTED LUNG SIZE
# (PFVC) rather than PREDICTED BODY WEIGHT (PBW) change outcomes? Two prescribable protocols are
# emulated with the clone-censor-weight engine of 35_tte_primary / 37_tte_discordance_benefit, on one of two EXPOSURE FAMILIES:
#
#   PBWPFVC_TTE_EXPO_FAMILY = "vt"   (PRIMARY)  tidal volume: VT/PFVC (% pred FVC) vs VT/PBW (mL/kg)
#                           = "mp"   (secondary) mechanical power: MP/PFVC (J/min/L) vs MP/PBW (J/min/kg)
#                           = "work" (secondary) work per breath: W/PFVC (J/breath/L) vs W/PBW -- power
#                             with the respiratory rate removed, so the shared numerator varies less
#                             and the normalizer contrast is correspondingly larger (see section 1)
#
# Two designs per family (section 3):
#   HEAD-TO-HEAD ("ceiling"):  X/PFVC <= TAU_PFVC   vs   X/PBW <= TAU_PBW, BITE-MATCHED -- each ceiling
#     is the (1 - BITE) quantile of its own exposure over post-grace patient-days, so both arms
#     restrict the same share of ventilator-days and differ only in WHICH days. In per-lung terms
#     the PBW ceiling allows TAU_PBW x (PBW/PFVC) per litre -- more the more PBW oversizes the lung
#     -- while the PFVC ceiling allows the same TAU_PFVC to every lung. The policies coincide at the
#     crossover discordance; above it the PFVC ceiling is the stricter one. The BITE is symmetric by
#     construction (each ceiling is a quantile of its own distribution); OVERLAP of the day-level
#     adherence decision is a separate, empirical property -- clinicians set VT early and rarely
#     change it, so P(adhere) can still collapse in the tails -- and is CHECKED by read B, not assumed.
#     Primary design because both arms are prescribable and the contrast is the normalizer question
#     itself. RD = PFVC arm - PBW arm; negative favours PFVC-anchored dosing.
#   CAP-ON-TOP ("cap"):  [X/PBW <= TAU_PBW AND X/PFVC <= TAU_CAP]  vs  X/PBW <= TAU_PBW. Usual
#     PBW-guided care plus a PFVC safety cap, vs usual care alone -- the trial a clinician would run.
#     TAU_CAP is anchored so the cap COINCIDES with the PBW ceiling at the anchor discordance
#     (default the Concordant tertile's median PBW/PFVC), so below the anchor the cap never binds
#     beyond the PBW ceiling: the Concordant tertile is a NEGATIVE CONTROL and all contrast
#     accumulates in the discordant tail. The bite is NOT symmetric here (the cap arm is a strict
#     superset of restrictions, binding in the discordant, older tail), so it is the secondary,
#     clinical-translation design and is read against its per-tertile overlap tables.
#
# Both arms' deviation models condition on the SAME lagged history (S/F, and every normalization of
# yesterday's exposure; plus lagged peak pressure for the MP family), so the arms differ only in the
# deviation indicator, not in weight-model specification. Each design reports the 37_tte_discordance_benefit reads:
# (A) tertile RDs + Discordant-Concordant gradient, (A2) continuous CATE by discordance,
# (B) positivity by tertile x arm, (C) the per-lung correction each arm delivers.
#
# UNCERTAINTY. The design pieces (33_tte_ceiling_headtohead.R / 33_tte_ceiling_cap.R) carry a fixed-weight
# cluster bootstrap, which for two arms that share nearly all their clones conditions on one weight
# fit and is too narrow. 33_tte_ceiling_diagnostics.R therefore runs a weight-REFIT bootstrap (both arms
# rebuilt inside every replicate) and writes tte_<fam>_<design>_overall_refit_*: THAT is the primary
# interval for the manuscript; the fixed-weight CI is kept beside it.
#
# The MP family uses the tte_mppbw_additive_hte formula (0.098 x RR x VT[L] x PIP, daily median) so tte_mppbw_additive_hte and the MP
# pieces share one exposure; the VT family uses the engine's daily-median VT/PFVC directly. The
# cohort is the engine's `base`, which carries the VT/PFVC structural-positivity exclusion
# (min feasible strain > C_LOW). Synthetic survival is simulated => plumbing only.
# Env: PBWPFVC_TTE_EXPO_FAMILY (vt|mp), PBWPFVC_CEIL_BITE (default 0.25), PBWPFVC_CEIL_TAU_PBW /
# _TAU_PFVC (explicit ceilings), PBWPFVC_CEIL_CAP_ANCHOR ("concordant_median" or a number),
# PBWPFVC_DISC_CATE_DF, PBWPFVC_NBOOT, PBWPFVC_CORES.
# =============================================================================
library(here); library(patchwork)
source(here::here("code", "31_tte_engine.R"))
FAM <- Sys.getenv("PBWPFVC_TTE_EXPO_FAMILY", "vt")
stopifnot(FAM %in% c("vt", "mp", "work"))
.ceil_key <- paste(site_name, FAM)
# guard: rebuild only if never built this session, or the engine/family changed
if (!exists(".ceil_loaded_key") || !identical(.ceil_loaded_key, .ceil_key)) {

  BITE        <- as.numeric(Sys.getenv("PBWPFVC_CEIL_BITE", "0.25"))      # share of post-grace days each ceiling binds
  TAU_PBW_ENV <- suppressWarnings(as.numeric(Sys.getenv("PBWPFVC_CEIL_TAU_PBW",  unset = NA)))
  TAU_PFVC_ENV<- suppressWarnings(as.numeric(Sys.getenv("PBWPFVC_CEIL_TAU_PFVC", unset = NA)))
  stopifnot(is.finite(BITE), BITE > 0, BITE < 1)
  DISC_LEVELS <- c("Concordant", "Mid", "Discordant")
  # x_pfvc = x_pbw x (PBW/PFVC) x K_CONV: VT/PFVC is in % of predicted FVC (x 0.1), MP/PFVC in J/min/L
  K_CONV <- if (FAM == "vt") 0.1 else 1
  UNITS  <- switch(FAM,
    vt   = c(raw = "VT (mL)", pbw = "VT/PBW (mL/kg)", pfvc = "VT/PFVC (% pred FVC)"),
    mp   = c(raw = "MP (J/min)", pbw = "MP/PBW (J/min/kg)", pfvc = "MP/PFVC (J/min/L)"),
    work = c(raw = "W (J/breath)", pbw = "W/PBW (J/breath/kg)", pfvc = "W/PFVC (J/breath/L)"))
  FAM_LAB <- switch(FAM, vt = "tidal-volume", mp = "mechanical-power", work = "work-per-breath")
  PFX     <- paste0("tte_", FAM, "_")

  # =============================================================================
  # 1. Daily exposure panel: x_raw / x_pbw / x_pfvc on the engine panel
  # =============================================================================
  if (FAM == "vt") {
    # The engine panel already carries the daily-median VT/PFVC; VT/PBW follows from the identity.
    panel_x <- panel %>%
      mutate(x_pfvc = vtpfvc, x_raw = vtpfvc * 10 * pfvc, x_pbw = x_raw / pbw) %>%
      filter(is.finite(x_pbw), is.finite(x_pfvc))
  } else {
    # Raw `wf` is cache-blocklisted, so the waterfall is re-read. Filters are tte_mppbw_additive_hte's. Daily MEDIAN
    # of the per-record exposure + daily median PIP (a lagged confounder for these families).
    #   mp   = 0.098 x RR x VT[L] x PIP   (J/min)     -- power, carries the respiratory rate
    #   work = 0.098 x      VT[L] x PIP   (J/breath)  -- the same thing per CYCLE, rate removed
    # Work per breath is the energy of one stress-strain event, and stripping the rate removes
    # the component with the widest spread. Because the two ceilings of a head-to-head differ
    # ONLY by the (demographic) normalizer ratio, their disagreement is largest when the shared
    # numerator varies least: VT (protocolized 6-8 mL/kg) gave ~20% of days, MP ~7%, and work per
    # breath should sit between them. That is the point of running this family.
    KC <- 0.098
    x_daily <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet")) %>%
      select(hospitalization_id, recorded_dttm, tidal_volume_set, resp_rate_set, peep_set,
             peak_inspiratory_pressure_obs) %>%
      inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
      mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
      filter(vent_day >= 0, vent_day <= MAX_VENT_DAY,
             !is.na(tidal_volume_set), between(tidal_volume_set, 50, 2000),
             !is.na(resp_rate_set), between(resp_rate_set, 4, 60),
             !is.na(peep_set), !is.na(peak_inspiratory_pressure_obs),
             between(peak_inspiratory_pressure_obs, 5, 80)) %>%
      mutate(rate_term = if (FAM == "mp") resp_rate_set else 1,
             x_rec = KC * rate_term * (tidal_volume_set / 1000) * peak_inspiratory_pressure_obs) %>%
      group_by(hospitalization_id, vent_day) %>%
      summarise(x_raw = median(x_rec), pip = median(peak_inspiratory_pressure_obs), .groups = "drop")
    # INNER join: an NA exposure would make `above` NA and corrupt the cumulative deviation flag
    # inside arm_build, so panel-days without a recorded peak pressure are dropped and REPORTED.
    panel_x <- panel %>%
      inner_join(x_daily, by = c("hospitalization_id", "vent_day")) %>%
      mutate(x_pbw = x_raw / pbw, x_pfvc = x_raw / pfvc) %>%
      filter(is.finite(x_pbw), is.finite(x_pfvc))
  }
  x_drop <- tibble(
    family = FAM, panel_days = nrow(panel), x_days = nrow(panel_x),
    frac_days_dropped = 1 - nrow(panel_x) / nrow(panel),
    panel_patients = n_distinct(panel$hospitalization_id),
    x_patients = n_distinct(panel_x$hospitalization_id),
    patients_lost = n_distinct(panel$hospitalization_id) - n_distinct(panel_x$hospitalization_id))
  message(sprintf("33_tte_ceiling [%s]: %d of %d patient-days carry the exposure (%.1f%% dropped); %d of %d patients",
                  FAM, x_drop$x_days, x_drop$panel_days, 100 * x_drop$frac_days_dropped,
                  x_drop$x_patients, x_drop$panel_patients))
  if (x_drop$x_patients < 100) stop("33_tte_ceiling: fewer than 100 patients carry the exposure panel -- not estimable here.")

  # =============================================================================
  # 2. Bite-matched ceilings + where the two policies disagree
  # =============================================================================
  post <- panel_x %>% filter(vent_day > GRACE)
  tau_pbw  <- if (is.finite(TAU_PBW_ENV))  TAU_PBW_ENV  else unname(quantile(post$x_pbw,  1 - BITE))
  tau_pfvc <- if (is.finite(TAU_PFVC_ENV)) TAU_PFVC_ENV else unname(quantile(post$x_pfvc, 1 - BITE))
  post <- post %>% mutate(above_pbw = x_pbw > tau_pbw, above_pfvc = x_pfvc > tau_pfvc)
  thresholds <- tibble(
    family = FAM, bite_target = BITE,
    anchoring = if (is.finite(TAU_PBW_ENV) || is.finite(TAU_PFVC_ENV)) "explicit_env" else "bite_matched_quantile",
    unit_pbw = UNITS["pbw"], unit_pfvc = UNITS["pfvc"],
    tau_pbw = tau_pbw, tau_pfvc = tau_pfvc,
    implied_abs_pbw_arm  = tau_pbw  * median(base$pbw),            # raw units at the median PBW
    implied_abs_pfvc_arm = tau_pfvc * median(base$pfvc) / K_CONV,  # raw units at the median PFVC
    crossover_discordance = tau_pfvc / (tau_pbw * K_CONV),         # PBW/PFVC above which the PFVC ceiling is stricter
    frac_days_above_pbw  = mean(post$above_pbw), frac_days_above_pfvc = mean(post$above_pfvc),
    frac_days_disagree   = mean(post$above_pbw != post$above_pfvc),
    frac_days_pfvc_only  = mean(post$above_pfvc & !post$above_pbw),
    frac_days_pbw_only   = mean(post$above_pbw & !post$above_pfvc),
    n_post_grace_days = nrow(post))
  if (thresholds$frac_days_disagree == 0)
    stop("33_tte_ceiling: the two ceilings bind on exactly the same days -- no head-to-head contrast exists at this bite.")
  write_csv(bind_cols(thresholds, x_drop %>% select(-family)),
            file.path(final_dir, paste0(PFX, "ceiling_thresholds_", site_name, ".csv")))
  message(sprintf("33_tte_ceiling [%s]: %s <= %.3g (~%.0f raw), %s <= %.3g (~%.0f raw); crossover discordance %.1f; %.1f%% of post-grace days disagree",
                  FAM, UNITS["pbw"], tau_pbw, thresholds$implied_abs_pbw_arm, UNITS["pfvc"], tau_pfvc,
                  thresholds$implied_abs_pfvc_arm, thresholds$crossover_discordance, 100 * thresholds$frac_days_disagree))

  # =============================================================================
  # 3. The two contrasts (see header) -- cap anchor, cap driver column, shared lagged history
  # =============================================================================
  # The cap arm's exposure is the day's value as a RATIO to its tighter ceiling,
  # max(x_pbw / TAU_PBW, x_pfvc / TAU_CAP), with ceiling 1 -- one column drives the deviation rule.
  anchor_env <- Sys.getenv("PBWPFVC_CEIL_CAP_ANCHOR", "concordant_median")
  disc_anchor <- if (identical(anchor_env, "concordant_median")) {
    base %>% filter(disc_grp == "Concordant") %>% summarise(m = median(pbw / pfvc)) %>% pull(m)
  } else as.numeric(anchor_env)
  stopifnot(is.finite(disc_anchor), disc_anchor > 0)
  tau_cap <- tau_pbw * disc_anchor * K_CONV
  panel_x <- panel_x %>% mutate(x_cap = pmax(x_pbw / tau_pbw, x_pfvc / tau_cap))
  # SHARED lagged history for the deviation (IPCW denominator) models. arm_build adds the arm's own
  # lagged exposure as l_expo; every OTHER lagged exposure (+ lagged peak pressure, MP family) is
  # passed through sf_term so both arms of a design condition on the same column space. Without
  # this each arm's weights came from a differently-specified model and the Concordant negative
  # control of the cap design sat at -0.6 pp for a 1% difference in binding days. Lags are row-based
  # within patient, exactly as arm_build computes l_expo, so they are NA on the same rows.
  panel_x <- panel_x %>% group_by(hospitalization_id) %>% arrange(vent_day) %>%
    mutate(l_x_pbw = lag(x_pbw), l_x_pfvc = lag(x_pfvc), l_x_cap = lag(x_cap),
           l_pip = if (FAM == "mp") lag(pip) else NA_real_) %>%
    ungroup()
  # The cap ratio's lag enters only the cap design's models: the head-to-head is the primary
  # contrast and must not depend on the cap anchor, which plays no role in it.
  SHARED_LAGS <- c("l_sf", if (FAM == "mp") "l_pip", "l_x_pbw", "l_x_pfvc")
  sf_for <- function(expo, design = "ceiling") {
    lags <- if (identical(design, "cap")) c(SHARED_LAGS, "l_x_cap") else SHARED_LAGS
    paste(setdiff(lags, paste0("l_", expo)), collapse = " + ")
  }
  post <- post %>% mutate(above_cap = x_pbw > tau_pbw | x_pfvc > tau_cap,
                          cap_only  = above_cap & !above_pbw)
  cap_thresholds <- tibble(
    family = FAM, anchoring = anchor_env, anchor_discordance = disc_anchor,
    unit_pbw = UNITS["pbw"], unit_pfvc = UNITS["pfvc"],
    tau_pbw = tau_pbw, tau_cap = tau_cap,
    implied_abs_cap_at_median_pfvc = tau_cap * median(base$pfvc) / K_CONV,
    frac_days_above_pbw = mean(post$above_pbw), frac_days_above_cap_arm = mean(post$above_cap),
    frac_days_cap_only  = mean(post$cap_only),                 # = frac_days_disagree (cap arm is a superset)
    n_post_grace_days = nrow(post))
  cap_by_tertile <- post %>% group_by(disc_grp) %>%
    summarise(frac_days_above_pbw = mean(above_pbw), frac_days_above_cap_arm = mean(above_cap),
              frac_days_cap_only = mean(cap_only),
              frac_patients_any_cap_only = n_distinct(hospitalization_id[cap_only]) / n_distinct(hospitalization_id),
              .groups = "drop") %>% mutate(family = FAM, disc_grp = factor(disc_grp, DISC_LEVELS)) %>% arrange(disc_grp)
  if (cap_thresholds$frac_days_cap_only == 0)
    stop("33_tte_ceiling: the PFVC cap never binds beyond the PBW ceiling -- no cap-on-top contrast at this anchor.")
  write_csv(cap_thresholds, file.path(final_dir, paste0(PFX, "cap_thresholds_", site_name, ".csv")))
  write_csv(cap_by_tertile, file.path(final_dir, paste0(PFX, "cap_bite_by_tertile_", site_name, ".csv")))
  message(sprintf("33_tte_ceiling [%s] cap: anchor discordance %.2f -> %s cap %.3g; cap binds beyond the PBW ceiling on %.1f%% of post-grace days (Concordant %.1f%% / Mid %.1f%% / Discordant %.1f%%)",
                  FAM, disc_anchor, UNITS["pfvc"], tau_cap, 100 * cap_thresholds$frac_days_cap_only,
                  100 * cap_by_tertile$frac_days_cap_only[1], 100 * cap_by_tertile$frac_days_cap_only[2],
                  100 * cap_by_tertile$frac_days_cap_only[3]))

# =============================================================================
# 4. One full analysis per design: design build -> HTE (37_tte_discordance_benefit read A) -> bootstrap -> tables -> figures
# =============================================================================
# --- shared estimator pieces (top level so the diagnostics piece can reuse them) ----------------
# Engine slot "strain_limiting" = the PFVC-informed arm, "permissive" = the PBW-ceiling arm in
# both designs; RD = PFVC-informed arm - PBW arm, NEGATIVE favours the PFVC-informed policy.
build_ceiling_design <- function(design, pnl = panel_x, keep_pday = FALSE) {
  is_cap  <- identical(design, "cap")
  expo_sl <- if (is_cap) "x_cap" else "x_pfvc"
  build_design(if (is_cap) 1 else tau_pfvc, tau_pbw, GRACE, DAYW_CAP, "simple",
               pnl = pnl, conf = NULL, keep_pday = keep_pday,
               sf_term = c(sf_for(expo_sl, design), sf_for("x_pbw", design)), expo = c(expo_sl, "x_pbw"))
}
FORM   <- died ~ arm * disc_grp + arm * ns(day, 4) + disc_grp * ns(day, 4) + sofa_total
arm_f  <- function() factor(c("permissive", "strain_limiting"), c("permissive", "strain_limiting"))
std_rd <- function(fit, prof, g) {
  pp <- if (g == "All") prof else prof %>% filter(as.character(disc_grp) == g)
  if (nrow(pp) < 50) return(c(rd = NA_real_, risk_sl = NA_real_, risk_pm = NA_real_))
  cells <- pp %>% count(disc_grp, sofa_total, name = "wt")
  grid  <- tidyr::crossing(cells, day = 1:HORIZON, arm = arm_f())
  grid$haz <- predict(fit, grid, type = "response")
  ci <- grid %>% group_by(arm, day) %>% summarise(h = weighted.mean(haz, wt), .groups = "drop") %>%
    group_by(arm) %>% arrange(day) %>% summarise(ci = 1 - prod(1 - h), .groups = "drop")
  rs <- ci$ci[ci$arm == "strain_limiting"]; rp <- ci$ci[ci$arm == "permissive"]
  c(rd = rs - rp, risk_sl = rs, risk_pm = rp)
}
design_prefix <- function(design) paste0(PFX, if (identical(design, "cap")) "cap_" else "ceiling_")

analyse_design <- function(design) {
  is_cap <- identical(design, "cap")
  prefix <- design_prefix(design)
  tag    <- function(x) paste0(prefix, x, "_", site_name)
  ARM    <- if (is_cap) c(strain_limiting = "pbw_plus_pfvc_cap", permissive = "pbw_ceiling")
            else        c(strain_limiting = "pfvc_ceiling",       permissive = "pbw_ceiling")
  lab_sl <- if (is_cap) sprintf("%s <= %.3g AND %s <= %.3g", UNITS["pbw"], tau_pbw, UNITS["pfvc"], tau_cap)
            else        sprintf("%s <= %.3g", UNITS["pfvc"], tau_pfvc)
  lab_pm <- sprintf("%s <= %.3g", UNITS["pbw"], tau_pbw)
  message("\n=== 33_tte_ceiling [", FAM, "] design '", design, "': ", ARM["strain_limiting"], " vs ", ARM["permissive"], " ===")

  des <- build_ceiling_design(design, keep_pday = TRUE)
  long <- des$long; lib <- des$lib; ids <- unique(long$hospitalization_id)
  point <- rd_from(long); lib_pt <- lib_diff(lib)

  # --- HTE: pooled SOFA-adjusted IPC-weighted MSM, tertiles + continuous (37_tte_discordance_benefit read A) ---------
  long_s   <- long %>% left_join(base %>% select(hospitalization_id, sofa_total), by = "hospitalization_id")
  prof_all <- long_s %>% distinct(hospitalization_id, disc_grp, sofa_total)
  n_disc   <- prof_all %>% count(disc_grp) %>% mutate(disc_grp = as.character(disc_grp))
  fit0  <- suppressWarnings(glm(FORM, data = long_s, family = binomial, weights = ipcw))
  pt_l  <- setNames(lapply(DISC_LEVELS, function(g) std_rd(fit0, prof_all, g)), DISC_LEVELS)
  pt    <- vapply(pt_l, `[`, numeric(1), "rd")
  overall_rd_adj <- unname(std_rd(fit0, prof_all, "All")["rd"])
  grad_pt <- unname(pt["Discordant"] - pt["Concordant"])

  KD      <- as.integer(Sys.getenv("PBWPFVC_DISC_CATE_DF", "3"))
  disc_c  <- base %>% transmute(hospitalization_id, ldisc = log(pbw / pfvc))
  long_c  <- long_s %>% left_join(disc_c, by = "hospitalization_id")
  prof_c  <- long_c %>% distinct(hospitalization_id, ldisc, sofa_total)
  Bspl    <- ns(prof_c$ldisc, df = KD)
  zc      <- paste0("z", seq_len(KD))
  add_z   <- function(d) { m <- predict(Bspl, d$ldisc); for (j in seq_len(KD)) d[[zc[j]]] <- m[, j]; d }
  long_cz <- add_z(long_c)
  FORM_C  <- as.formula(paste0("died ~ arm*(", paste(zc, collapse = "+"), ") + arm*ns(day,4) + (",
                               paste(zc, collapse = "+"), ")*ns(day,4) + sofa_total"))
  DGRID   <- seq(quantile(prof_c$ldisc, .025), quantile(prof_c$ldisc, .975), length.out = 40)
  qd      <- quantile(prof_c$ldisc, c(.10, .90))
  EVAL    <- sort(unique(c(DGRID, qd)))
  sofa_cells <- prof_c %>% count(sofa_total, name = "wt")
  grid_c  <- tidyr::crossing(ldisc = EVAL, sofa_cells, day = 1:HORIZON, arm = arm_f()) %>% add_z()
  cate_from_fit <- function(fit) {
    g <- grid_c; g$haz <- predict(fit, g, type = "response")
    g %>% group_by(ldisc, arm, day) %>% summarise(h = weighted.mean(haz, wt), .groups = "drop") %>%
      group_by(ldisc, arm) %>% arrange(day) %>% summarise(cif = 1 - prod(1 - h), .groups = "drop") %>%
      group_by(ldisc) %>%
      summarise(rd = cif[arm == "strain_limiting"] - cif[arm == "permissive"], .groups = "drop") %>%
      arrange(ldisc)
  }
  fitc0 <- suppressWarnings(glm(FORM_C, data = long_cz, family = binomial, weights = ipcw))
  cate0 <- cate_from_fit(fitc0)
  i10 <- which.min(abs(cate0$ldisc - qd[1])); i90 <- which.min(abs(cate0$ldisc - qd[2]))
  slope_pt <- cate0$rd[i90] - cate0$rd[i10]

  # --- ONE fixed-weight cluster bootstrap: overall + liberation + tertiles + gradient + CATE curve --
  # (the weight-REFIT bootstrap in 33_tte_ceiling_diagnostics.R supplies the primary interval)
  NCURVE <- length(EVAL)
  boot_template <- c(overall = NA_real_, lib = NA_real_, setNames(rep(NA_real_, 3), DISC_LEVELS),
                     gradient = NA_real_, setNames(rep(NA_real_, NCURVE), paste0("c", seq_len(NCURVE))))
  boot_one <- function() {
    samp <- tibble(hospitalization_id = sample(ids, replace = TRUE))
    out  <- boot_template
    bl   <- long %>% inner_join(samp, by = "hospitalization_id", relationship = "many-to-many")
    blib <- lib  %>% inner_join(samp, by = "hospitalization_id", relationship = "many-to-many")
    out["overall"] <- tryCatch(unname(rd_from(bl)["rd"]), error = function(e) NA_real_)
    out["lib"]     <- tryCatch(lib_diff(blib),            error = function(e) NA_real_)
    lb  <- long_cz %>% inner_join(samp, by = "hospitalization_id", relationship = "many-to-many")
    pb  <- samp %>% left_join(prof_all, by = "hospitalization_id")
    fb  <- suppressWarnings(tryCatch(glm(FORM, data = lb, family = binomial, weights = ipcw), error = function(e) NULL))
    if (!is.null(fb)) {
      r <- vapply(DISC_LEVELS, function(g) unname(std_rd(fb, pb, g)["rd"]), numeric(1))
      out[DISC_LEVELS] <- r; out["gradient"] <- unname(r["Discordant"] - r["Concordant"])
    }
    fc  <- suppressWarnings(tryCatch(glm(FORM_C, data = lb, family = binomial, weights = ipcw), error = function(e) NULL))
    if (!is.null(fc)) {
      cv <- tryCatch(cate_from_fit(fc)$rd, error = function(e) rep(NA_real_, NCURVE))
      if (length(cv) == NCURVE) out[paste0("c", seq_len(NCURVE))] <- cv
    }
    out
  }
  n_cores_used <- min(N_CORES, N_BOOT)
  message("33_tte_ceiling '", design, "' fixed-weight cluster bootstrap (", N_BOOT, " reps across ", n_cores_used, " core(s)) ...")
  boot_t0 <- Sys.time()
  chunks   <- split(seq_len(N_BOOT), cut(seq_len(N_BOOT), min(20L, N_BOOT), labels = FALSE))
  bts_list <- vector("list", N_BOOT); done <- 0L
  if (n_cores_used > 1) {
    cl <- makeCluster(n_cores_used, type = "PSOCK")
    clusterEvalQ(cl, { library(tidyverse); library(splines); library(survival) })
    clusterExport(cl, envir = .GlobalEnv, varlist = c(
      "HORIZON", "WT_TRUNC", "DISC_LEVELS", "rd_from", "ci_curve", "cif_lib", "lib_diff", "trunc_w",
      "FORM", "std_rd", "arm_f"))
    clusterExport(cl, envir = environment(), varlist = c(
      "ids", "long", "lib", "long_cz", "prof_all", "boot_template", "boot_one", "FORM_C",
      "cate_from_fit", "grid_c", "NCURVE"))
    clusterSetRNGStream(cl, if (is_cap) 20260905L else 20260904L)
    tryCatch(
      for (ch in chunks) {
        bts_list[ch] <- parLapply(cl, ch, function(bb) boot_one())
        done <- done + length(ch)
        message(sprintf("  bootstrap %d/%d (%2d%%) | elapsed %4.0fs", done, N_BOOT,
                        as.integer(round(100 * done / N_BOOT)),
                        as.numeric(difftime(Sys.time(), boot_t0, units = "secs"))))
      }, finally = stopCluster(cl))
  } else {
    for (ch in chunks) {
      for (b in ch) bts_list[[b]] <- boot_one()
      done <- done + length(ch); message(sprintf("  bootstrap %d/%d", done, N_BOOT))
    }
  }
  bts <- do.call(rbind, bts_list)
  ci  <- function(col) quantile(bts[, col], c(.025, .975), na.rm = TRUE)
  n_bad <- sum(apply(bts[, paste0("c", seq_len(NCURVE)), drop = FALSE], 1, function(r) all(is.na(r))))
  if (n_bad > 0) message(sprintf("  CATE curve: %d/%d bootstrap reps failed -> dropped via na.rm%s",
          n_bad, N_BOOT, if (n_bad > N_BOOT / 2) " (WARNING: >half failed -- band unreliable)" else ""))

  # --- tables -----------------------------------------------------------------------------------
  evalue <- function(rr) { rr <- if (rr >= 1) rr else 1 / rr; rr + sqrt(rr * (rr - 1)) }
  rr_point <- unname(point["risk_sl"] / point["risk_pm"])
  rr_lo <- unname((point["risk_pm"] + ci("overall")[1]) / point["risk_pm"])
  rr_hi <- unname((point["risk_pm"] + ci("overall")[2]) / point["risk_pm"])
  rr_near <- if (rr_point >= 1) min(rr_lo, rr_hi) else max(rr_lo, rr_hi)
  rr_near <- if ((rr_point >= 1) != (rr_near >= 1)) 1 else rr_near
  overall <- tibble(
    family = FAM, design = design, contrast = paste0(ARM["strain_limiting"], "_minus_", ARM["permissive"]),
    unit_pbw = UNITS["pbw"], unit_pfvc = UNITS["pfvc"],
    tau_pbw = tau_pbw, tau_pfvc = if (is_cap) NA_real_ else tau_pfvc, tau_cap = if (is_cap) tau_cap else NA_real_,
    risk_pfvc_arm = unname(point["risk_sl"]), risk_pbw_arm = unname(point["risk_pm"]),
    rd = unname(point["rd"]), rd_lo = ci("overall")[1], rd_hi = ci("overall")[2],
    ci_type = "fixed_weight_bootstrap",
    rd_sofa_adjusted = overall_rd_adj,
    lib_diff = lib_pt, lib_lo = ci("lib")[1], lib_hi = ci("lib")[2],
    evalue_point = evalue(rr_point), evalue_ci = evalue(rr_near),
    n_patients = length(ids),
    frac_trimmed_pfvc_arm = mean(is.finite(des$bl$idsum$trim_day)),
    frac_trimmed_pbw_arm  = mean(is.finite(des$bh$idsum$trim_day)),
    frac_deviated_pfvc_arm = mean(is.finite(des$bl$idsum$dev_day)),
    frac_deviated_pbw_arm  = mean(is.finite(des$bh$idsum$dev_day)))
  write_csv(overall, file.path(final_dir, paste0(tag("overall"), ".csv")))

  hte <- tibble(disc_grp = DISC_LEVELS, rd = pt[DISC_LEVELS],
                rd_lo = vapply(DISC_LEVELS, function(k) ci(k)[1], numeric(1)),
                rd_hi = vapply(DISC_LEVELS, function(k) ci(k)[2], numeric(1)),
                risk_pfvc_arm = vapply(pt_l, `[`, numeric(1), "risk_sl"),
                risk_pbw_arm  = vapply(pt_l, `[`, numeric(1), "risk_pm"),
                adjustment = "sofa_baseline", ci_type = "fixed_weight_bootstrap",
                family = FAM, design = design) %>%
    left_join(n_disc, by = "disc_grp") %>% mutate(disc_grp = factor(disc_grp, DISC_LEVELS))
  gradient <- tibble(family = FAM, design = design, statistic = "RD(Discordant) - RD(Concordant)",
                     adjustment = "sofa_baseline", ci_type = "fixed_weight_bootstrap",
                     estimate = grad_pt, lo = ci("gradient")[1], hi = ci("gradient")[2])
  write_csv(hte,      file.path(final_dir, paste0(tag("disc_hte"), ".csv")))
  write_csv(gradient, file.path(final_dir, paste0(tag("disc_gradient"), ".csv")))

  cmat <- bts[, paste0("c", seq_len(NCURVE)), drop = FALSE]
  qlo  <- apply(cmat, 2, quantile, .025, na.rm = TRUE); qhi <- apply(cmat, 2, quantile, .975, na.rm = TRUE)
  slope_ci <- quantile(cmat[, i90] - cmat[, i10], c(.025, .975), na.rm = TRUE)
  curve_tbl <- tibble(family = FAM, design = design, ldisc = cate0$ldisc, discordance = exp(cate0$ldisc),
                      rd = cate0$rd, rd_lo = qlo, rd_hi = qhi) %>% filter(ldisc %in% DGRID)
  slope_tbl <- tibble(family = FAM, design = design, statistic = "RD(p90 discordance) - RD(p10 discordance)",
                      adjustment = "sofa_baseline", ci_type = "fixed_weight_bootstrap",
                      disc_p10 = exp(qd[1]), disc_p90 = exp(qd[2]),
                      reference_discordance = if (is_cap) disc_anchor else thresholds$crossover_discordance,
                      rd_p10 = cate0$rd[i10], rd_p90 = cate0$rd[i90],
                      estimate = slope_pt, lo = slope_ci[1], hi = slope_ci[2])
  write_csv(curve_tbl, file.path(final_dir, paste0(tag("disc_cate_curve"), ".csv")))
  write_csv(slope_tbl, file.path(final_dir, paste0(tag("disc_cate_slope"), ".csv")))

  # --- read B: positivity by tertile x arm ---------------------------------------------------------
  pday <- bind_rows(des$bl$pday %>% mutate(arm = unname(ARM["strain_limiting"])),
                    des$bh$pday %>% mutate(arm = unname(ARM["permissive"])))
  emp1 <- function(g) summarise(g, n_eligible_days = n(),
    frac_padhere_lt05 = mean(p_adhere < 0.05), frac_padhere_lt02 = mean(p_adhere < 0.02),
    p01_padhere = quantile(p_adhere, 0.01), median_padhere = median(p_adhere), .groups = "drop")
  overlap_emp <- bind_rows(
    pday %>% group_by(arm, disc_grp) %>% emp1() %>% mutate(disc_grp = as.character(disc_grp)),
    pday %>% group_by(arm) %>% emp1() %>% mutate(disc_grp = "All")) %>%
    mutate(family = FAM, design = design) %>% arrange(arm, disc_grp)
  support_for_arm <- function(b, arm_lab) {
    map_dfr(c(7L, 14L, 28L), function(d) {
      w_at_d <- b$wday %>% filter(vent_day <= d) %>% group_by(hospitalization_id) %>%
        arrange(vent_day) %>% summarise(w = last(cumw), .groups = "drop")
      b$idsum %>% filter(dev_day > d, trim_day > d, is.na(death_day) | death_day >= d) %>%
        left_join(w_at_d, by = "hospitalization_id") %>% mutate(w = trunc_w(coalesce(w, 1))) %>%
        group_by(disc_grp) %>%
        summarise(arm = arm_lab, day = d, n_atrisk = n(),
                  ess_atrisk = if (n() > 0) (sum(w)^2 / sum(w^2)) else 0, .groups = "drop")
    })
  }
  overlap_support <- bind_rows(support_for_arm(des$bl, unname(ARM["strain_limiting"])),
                               support_for_arm(des$bh, unname(ARM["permissive"]))) %>%
    mutate(family = FAM, design = design) %>% arrange(arm, day, disc_grp)
  overlap_loss <- bind_rows(des$bl$idsum %>% mutate(arm = unname(ARM["strain_limiting"])),
                            des$bh$idsum %>% mutate(arm = unname(ARM["permissive"]))) %>%
    group_by(arm, disc_grp) %>%
    summarise(n = n(), frac_deviated = mean(is.finite(dev_day)), frac_trimmed = mean(is.finite(trim_day)),
              mean_stab_w = mean(ipcw_term), .groups = "drop") %>%
    mutate(family = FAM, design = design) %>% arrange(arm, disc_grp)
  write_csv(overlap_emp,     file.path(final_dir, paste0(tag("overlap_empirical"), ".csv")))
  write_csv(overlap_support, file.path(final_dir, paste0(tag("overlap_support"), ".csv")))
  write_csv(overlap_loss,    file.path(final_dir, paste0(tag("overlap_loss"), ".csv")))

  # --- read C: what each arm delivers per tertile (per-lung cut on binding days, PFVC units) -------
  # PFVC-informed arm ceiling in per-lung terms: head-to-head = tau_pfvc; cap = min(PBW ceiling in
  # per-lung terms, tau_cap). The PBW ceiling in per-lung terms is tau_pbw x discordance x K_CONV.
  dose_tbl <- post %>%
    mutate(disc = pbw / pfvc,
           ceil_pbw_arm_per_lung = tau_pbw * disc * K_CONV,
           ceil_pfvc_arm_per_lung = if (is_cap) pmin(ceil_pbw_arm_per_lung, tau_cap) else tau_pfvc,
           above_pfvc_arm  = x_pfvc > ceil_pfvc_arm_per_lung,
           cut_pfvc_arm = pmax(x_pfvc - ceil_pfvc_arm_per_lung, 0),
           cut_pbw_arm  = pmax(x_pfvc - ceil_pbw_arm_per_lung, 0)) %>%
    group_by(disc_grp) %>%
    summarise(family = FAM, design = design, n_patients = n_distinct(hospitalization_id), n_pdays = n(),
              median_x_raw = median(x_raw), median_x_pbw = median(x_pbw), median_x_pfvc = median(x_pfvc),
              median_allowed_per_lung_pbw_arm  = median(ceil_pbw_arm_per_lung),
              median_allowed_per_lung_pfvc_arm = median(ceil_pfvc_arm_per_lung),
              frac_days_above_pbw_arm = mean(above_pbw), frac_days_above_pfvc_arm = mean(above_pfvc_arm),
              frac_days_pfvc_arm_only = mean(above_pfvc_arm & !above_pbw),
              cut_per_lung_pbw_arm  = median(cut_pbw_arm[above_pbw]),
              cut_per_lung_pfvc_arm = median(cut_pfvc_arm[above_pfvc_arm]),
              .groups = "drop") %>%
    mutate(disc_grp = factor(disc_grp, DISC_LEVELS)) %>% arrange(disc_grp)
  write_csv(dose_tbl, file.path(final_dir, paste0(tag("disc_dose"), ".csv")))

  # --- figures ------------------------------------------------------------------------------------
  cc <- ci_curve(long) %>%
    pivot_longer(c(strain_limiting, permissive), names_to = "arm", values_to = "cuminc") %>%
    mutate(arm = unname(ARM[arm]))
  cols <- setNames(c(okabe[3], okabe[1]), unname(ARM))
  labs_arm <- setNames(c(lab_sl, lab_pm), unname(ARM))
  p1 <- ggplot(cc, aes(day, 100 * cuminc, colour = arm)) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = cols, labels = labs_arm, name = NULL) +
    labs(x = "Days from index ventilation", y = "Per-protocol cumulative mortality (%)",
         title = if (is_cap) paste0("TTE: PBW ceiling + PFVC safety cap vs PBW ceiling alone (", FAM_LAB, ")")
                 else        paste0("Head-to-head TTE: PFVC-anchored vs PBW-anchored ", FAM_LAB, " ceiling"),
         subtitle = sprintf("%s%s - RD %+.2f pp [%.2f, %.2f] (fixed-weight CI; refit CI in %s); trim P(adhere) >= %g",
                            site_name, if (is_synthetic) " (SYNTHETIC - plumbing only)" else "",
                            100 * overall$rd, 100 * overall$rd_lo, 100 * overall$rd_hi,
                            paste0(prefix, "overall_refit"), TRIM_ALPHA)) +
    theme_minimal(base_size = 10) + theme(legend.position = "top")
  ggsave(file.path(final_dir, paste0(tag("cuminc"), ".pdf")), p1, width = 8, height = 5)

  disc_med <- base %>% mutate(dc = pbw / pfvc) %>% group_by(disc_grp) %>%
    summarise(disc_med = median(dc), .groups = "drop") %>% mutate(disc_grp = factor(disc_grp, DISC_LEVELS))
  hte_ov <- hte %>% left_join(disc_med, by = "disc_grp")
  ref_x  <- if (is_cap) disc_anchor else thresholds$crossover_discordance
  pA <- ggplot(curve_tbl, aes(discordance, 100 * rd)) +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_hline(yintercept = 100 * overall_rd_adj, linetype = "dashed", colour = "grey55") +
    geom_vline(xintercept = ref_x, linetype = "dotted", colour = "grey40") +
    geom_ribbon(aes(ymin = 100 * rd_lo, ymax = 100 * rd_hi), alpha = 0.15, fill = "#0072B2") +
    geom_line(colour = "#0072B2", linewidth = 1) +
    geom_pointrange(data = hte_ov, aes(x = disc_med, y = 100 * rd, ymin = 100 * rd_lo, ymax = 100 * rd_hi,
                                       colour = disc_grp), inherit.aes = FALSE) +
    scale_colour_manual(values = setNames(okabe[c(2, 1, 6)], DISC_LEVELS), name = "tertile") +
    labs(x = "PBW/PFVC discordance (higher = PBW oversizes -> misdosed)",
         y = paste0("CATE: 28-d mortality RD (pp)\n", ARM["strain_limiting"], " - ", ARM["permissive"]),
         title = "A. Continuous CATE by discordance (spline; SOFA-adjusted)",
         subtitle = sprintf("RD(p90)-RD(p10): %+.2f pp [%.2f, %.2f]\ndotted = %s %.1f (%s)",
                            100 * slope_pt, 100 * slope_ci[1], 100 * slope_ci[2],
                            if (is_cap) "cap anchor" else "crossover discordance", ref_x,
                            if (is_cap) "cap binds only to the right; Concordant = negative control"
                            else "PFVC ceiling stricter to the right")) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  dose_long <- dose_tbl %>% transmute(disc_grp, pbw = cut_per_lung_pbw_arm, pfvc = cut_per_lung_pfvc_arm) %>%
    pivot_longer(-disc_grp, names_to = "arm", values_to = "cut") %>%
    mutate(arm = unname(ARM[ifelse(arm == "pfvc", "strain_limiting", "permissive")]))
  pB <- ggplot(dose_long, aes(disc_grp, cut, fill = arm)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    scale_fill_manual(values = cols, name = NULL) +
    labs(x = NULL, y = paste0("Median cut on binding days\n(", UNITS["pfvc"], ")"),
         title = "C. Per-lung correction each arm delivers",
         subtitle = "Rising with discordance in the PFVC-informed arm only\n= the mechanical engine of A") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  fig <- pA + pB + patchwork::plot_annotation(
    title = paste0(if (is_cap) "PBW ceiling + PFVC cap vs PBW ceiling" else "PFVC-anchored vs PBW-anchored ceiling",
                   " (", FAM_LAB, ") by PBW/PFVC discordance - ", site_name,
                   if (is_synthetic) " (SYNTHETIC - plumbing only)" else ""),
    subtitle = paste0("Negative = the PFVC-informed arm lowers mortality more. Read with ", prefix,
                      "overlap_*: a collapsed Discordant ESS = extrapolation, not a targeting win."))
  ggsave(file.path(final_dir, paste0(tag("disc_benefit"), ".pdf")), fig, width = 12, height = 5)

  # --- console ------------------------------------------------------------------------------------
  cat(sprintf("\n=== 33_tte_ceiling [%s] '%s': %s vs %s ===\n", FAM, design, ARM["strain_limiting"], ARM["permissive"]))
  cat(sprintf("    28-d mortality RD: %+.2f pp [%.2f, %.2f] fixed-weight  (risks %.1f%% vs %.1f%%; SOFA-adjusted %+.2f pp)\n",
              100 * overall$rd, 100 * overall$rd_lo, 100 * overall$rd_hi,
              100 * overall$risk_pfvc_arm, 100 * overall$risk_pbw_arm, 100 * overall$rd_sofa_adjusted))
  cat(sprintf("    liberation CIF diff: %+.3f [%.3f, %.3f]; E-value point %.2f / CI %.2f; trimmed %.1f%% / %.1f%%; deviated %.1f%% / %.1f%%\n",
              overall$lib_diff, overall$lib_lo, overall$lib_hi, overall$evalue_point, overall$evalue_ci,
              100 * overall$frac_trimmed_pfvc_arm, 100 * overall$frac_trimmed_pbw_arm,
              100 * overall$frac_deviated_pfvc_arm, 100 * overall$frac_deviated_pbw_arm))
  cat("--- RD by discordance tertile (SOFA-adjusted, standardized) ---\n")
  print(as.data.frame(hte %>% transmute(disc_grp, n, rd_pp = round(100 * rd, 2),
          ci = sprintf("[%.2f, %.2f]", 100 * rd_lo, 100 * rd_hi))), row.names = FALSE)
  cat(sprintf("    GRADIENT RD(Discordant)-RD(Concordant): %+.2f pp [%.2f, %.2f]  (negative => PFVC-informed arm helps the misdosed more)\n",
              100 * grad_pt, 100 * gradient$lo, 100 * gradient$hi))
  cat(sprintf("    continuous RD(p90)-RD(p10): %+.2f pp [%.2f, %.2f]\n", 100 * slope_pt, 100 * slope_ci[1], 100 * slope_ci[2]))
  cat("--- Positivity by tertile (at-risk ESS, day 28) ---\n")
  print(as.data.frame(overlap_support %>% filter(day == 28L) %>%
          transmute(arm, disc_grp, n_atrisk, ess_atrisk = round(ess_atrisk, 1))), row.names = FALSE)
  cat("--- What each arm delivers, by tertile ---\n")
  print(as.data.frame(dose_tbl %>% transmute(disc_grp, n_patients,
          x_pbw = signif(median_x_pbw, 3), x_pfvc = signif(median_x_pfvc, 3),
          pct_above_pbw = round(100 * frac_days_above_pbw_arm), pct_above_pfvc_arm = round(100 * frac_days_above_pfvc_arm),
          pct_pfvc_arm_only = round(100 * frac_days_pfvc_arm_only),
          cut_pbw_arm = signif(cut_per_lung_pbw_arm, 3), cut_pfvc_arm = signif(cut_per_lung_pfvc_arm, 3))), row.names = FALSE)
  invisible(list(overall = overall, hte = hte, gradient = gradient, slope = slope_tbl, dose = dose_tbl))
}

.ceil_loaded_key <- .ceil_key
cat(sprintf("\n=== 33_tte_ceiling [%s] thresholds ===\n", FAM))
print(as.data.frame(thresholds %>% mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
print(as.data.frame(cap_by_tertile %>% mutate(across(where(is.numeric), ~ round(., 3)))), row.names = FALSE)
} else message("33_tte_ceiling_common: already built this session for ", .ceil_key, " -- reusing.")
