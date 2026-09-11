# =============================================================================
# Script 10 (common): Longitudinal target trial emulation -- SHARED setup + compute
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# Shared engine for the longitudinal TTE; sourced by every code/11.*_*.R leaf.
# This file holds everything the downstream 11.* analysis scripts share: setup,
# libraries, design knobs, the seeded RNG state, the baseline + panel + engine
# functions, and the PRIMARY design build (up to but NOT including the bootstrap).
# It performs NO file writes (pure compute + message()); each 11.* leaf writes its
# own CSV/PDF outputs. See the monolith for the full design rationale ([T1]-[T12]).
#
# After sourcing this file a caller has, in .GlobalEnv:
#   data:      base, panel, panel_full, daily, dp_daily, ph_daily, ph_daily_art,
#              pao2_daily, panel_drop_summary (object), struct_excl (object)
#   engine:    trunc_w, ess_frac, arm_build, make_long, build_design, ci_curve,
#              rd_from, cif_lib, lib_diff, sub_vars, sg_rd
#   primary:   des, long_all, lib_all, point, lib_pt, sg_point, ids
#   knobs:     C_LOW, C_HIGH, GRACE, VT_FLOOR_MLKG, TRIM_ALPHA, DEESC_FRAC,
#              DAYW_CAP, HORIZON, MAX_VENT_DAY, WT_TRUNC, is_synthetic, N_BOOT, N_CORES
#   misc:      okabe, final_dir, output_dir, site_name
# =============================================================================

# Pin BLAS to one thread per process so the PSOCK bootstrap workers don't each
# run all-core multithreaded BLAS (N workers x cores threads -> oversubscription).
# Must precede any BLAS use; inherited by the workers spawned later. (cf. script 06.)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")

library(data.table)   # rolling SpO2->FiO2 join for the daily worst SF ratio
library(tidyverse)    # loaded AFTER data.table so dplyr first()/last()/between() win
library(arrow)
library(here)
library(splines)
library(survival)
library(parallel)

source("utils/config.R")
site_name  <- config$site_name
output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
RNGkind("L'Ecuyer-CMRG"); set.seed(20260617)

# --- design knobs (provisional; see [T3]) ------------------------------------
C_LOW        <- as.numeric(Sys.getenv("PBWPFVC_TTE_CLOW", "11"))   # strain-limiting ceiling, VT/PFVC %
C_HIGH       <- as.numeric(Sys.getenv("PBWPFVC_TTE_CHIGH", "16"))  # permissive ceiling (~ usual care)
                        # ARMA-ANCHORED: on the PFVC strain scale the proven-protective LTV arm
                        # (6 mL/kg PBW) delivered VT/PFVC p75 = 11.2 -> C_LOW = 11 is "the strain
                        # ~75% of the protective arm met" (an achievable ceiling, not the over-
                        # stringent LTV median ~10 which also tightens positivity). The harmful HTV
                        # arm (12 mL/kg) sat at VT/PFVC p25 = 17.4 -> C_HIGH = 16 is just under HTV
                        # territory. Alternatives are swept in 11_sensitivities.R (strain-threshold
                        # sweep). env-overridable so the FVC_age25 normalizer arm can re-calibrate
                        # (FVC_age25 > PFVC => same numeric ceiling is more lenient there; ARMA
                        # FVC_age25 anchors = LTV p75 9.4 / HTV p25 16.9). C_LOW < C_HIGH enforced below.
stopifnot(is.finite(C_LOW), is.finite(C_HIGH), C_LOW < C_HIGH)
GRACE        <- 1L      # days allowed above ceiling before deviation. 1d (one calendar
                        # day to titrate to protective settings, then enforce) mirrors
                        # ARMA, which reduced VT to 6 mL/kg over hours and checked plateau
                        # >=4h -- a 2d amnesty was more permissive than the trial standard.
                        # 11.D brackets it with grace 0 (strict) / 2 / 3 (lenient); RD is
                        # monotone in grace (stricter = larger) but stable in direction.
VT_FLOOR_MLKG <- 4      # lowest VT any LTVV protocol delivers (mL/kg PBW); the
                        # structural-positivity floor: if even VT_FLOOR_MLKG x PBW
                        # gives VT/PFVC > ceiling, adherence is physically impossible
TRIM_ALPHA   <- 0.02    # common-support trim ([T9], PRIMARY is trimmed): censor a
                        # clone at the first eligible day where modeled P(adhere) =
                        # 1 - p_den < TRIM_ALPHA (deviation near-certain -> the
                        # positivity-violation tail where weights explode). Lower-tail
                        # only; restricts the estimand to the overlap region. trim=0
                        # disables. Decisive at UCSF (MARGINAL overlap on this tail), ~free at MIMIC.
DEESC_FRAC   <- 0.05    # de-escalation MTP ([T10] sensitivity): a clone above the
                        # ceiling must cut next-day VT by >= DEESC_FRAC (relative) to
                        # stay adherent -- the dynamic de-escalation decision (the one with day-level overlap).
DAYW_CAP     <- 5       # day-weight truncation (IPCW clamped to [1/cap, cap]); see weight-cap sensitivity
HORIZON      <- 28L     # days, primary outcome (28-d: bulk of ICU mortality, ~= adherence window)
MAX_VENT_DAY <- 27L     # ventilation/adherence window
WT_TRUNC     <- c(0.01, 0.99)
is_synthetic <- identical(site_name, "synthetic_clif")
# cluster bootstrap reps; override on ANY site via PBWPFVC_NBOOT (e.g. =25 to prototype fast --
# each rep refits the MSM on a resampled long panel, so this is the dominant cost at real scale).
# Unset -> 100 (synthetic) / 500 (real).
N_BOOT       <- suppressWarnings(as.integer(Sys.getenv("PBWPFVC_NBOOT", unset = NA)))
if (is.na(N_BOOT)) N_BOOT <- if (is_synthetic) 100L else 500L
message("TTE: N_BOOT = ", N_BOOT, "  (PBWPFVC_NBOOT env = '", Sys.getenv("PBWPFVC_NBOOT"), "')")
# bootstrap worker count; override with PBWPFVC_CORES (default leaves 1 core free)
N_CORES      <- suppressWarnings(as.integer(Sys.getenv("PBWPFVC_CORES", unset = NA)))
if (is.na(N_CORES)) N_CORES <- max(1L, detectCores() - 1L)
message("Site: ", site_name, " | ceilings ", C_LOW, "/", C_HIGH, " | grace ", GRACE, "d")

# =============================================================================
# 10a. Baseline (PFVC, demographics, t0, death day) -- mortality handling
#      identical to scripts 06-08 incl. synthetic-only workaround.
# =============================================================================
# Cohort source. DEFAULT = the primary cross-sectional cohort (index selected at VT/PBW
# 6-8 mL/kg). Set PBWPFVC_TTE_VTPBW="lo,hi" to RE-SELECT the index over a wider VT/PBW band
# from the pre-gate per-timepoint dataset (analysis_all_eligible_timepoints) -- a positivity
# sensitivity that lets smaller lungs be dosed below 6 mL/kg and actually reach the strain
# ceiling. Same two-tier index rule and SF<315 hypoxemia gate as script 03; only the VT/PBW
# band changes. Outputs go to a final/vtpbw_<lo>_<hi>/ subfolder so they never collide with
# the primary cohort's outputs. Every 11.* leaf inherits this final_dir unchanged.
vtpbw_band <- Sys.getenv("PBWPFVC_TTE_VTPBW", "")
if (nzchar(vtpbw_band)) {
  bw <- suppressWarnings(as.numeric(strsplit(vtpbw_band, ",")[[1]]))
  stopifnot(length(bw) == 2, !any(is.na(bw)), bw[1] < bw[2])
  message("*** BROADER COHORT: re-selecting the index over VT/PBW [", bw[1], ", ", bw[2],
          "] from analysis_all_eligible_timepoints (positivity sensitivity) ***")
  ae <- read_parquet(file.path(output_dir, "analysis_all_eligible_timepoints.parquet"))
  qual <- ae %>% filter(has_all_data, vtpbw >= bw[1], vtpbw <= bw[2], sf_ratio < 315)
  imv_start <- ae %>% group_by(hospitalization_id) %>%
    summarise(imv_start_dttm = min(recorded_dttm), .groups = "drop")
  t1 <- qual %>% left_join(imv_start, by = "hospitalization_id") %>%
    filter(!is.na(dp), recorded_dttm <= imv_start_dttm + lubridate::hours(6)) %>%
    group_by(hospitalization_id) %>% slice_min(recorded_dttm, n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(-imv_start_dttm)
  t2 <- qual %>% filter(!hospitalization_id %in% t1$hospitalization_id) %>%
    group_by(hospitalization_id) %>% slice_min(recorded_dttm, n = 1, with_ties = FALSE) %>% ungroup()
  cs <- bind_rows(t1, t2)
  final_dir <- file.path(final_dir, paste0("vtpbw_", bw[1], "_", bw[2]))
  dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
  message("  broader cohort: ", nrow(cs), " patients; outputs -> ", final_dir)
} else {
  cs <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"))
}

# --- normalizer switch (PBWPFVC_TTE_NORM = pfvc | pfvc_age25) -----------------
# The CEILING estimand is normalizer-DEPENDENT (unlike the 12 multiplicative shift, where the
# normalizer cancels). Default = pfvc (the published age+sex+race+height correction). pfvc_age25
# = the STRUCTURAL-only normalizer (GLI sex/race/height with age pinned to 25, derived in 03):
# it KEEPS the size correction validated by the race-invariance of specific elastance (05c E)
# and DROPS the age channel, which is both prognostically empty (05c D) AND the positivity-
# breaker here -- min feasible VT/normalizer = VT_FLOOR x (pbw/normalizer) x 0.1 scales with the
# discordance, ~77% of which is age. So pfvc_age25 is the "validated-correction, fixed-positivity"
# robustness arm. Outputs go to final/norm_pfvc_age25/ so they never collide with the primary;
# every 11.* leaf inherits final_dir unchanged. base ALWAYS carries pfvc_age25 (the 12 secondary
# CATE + the 11.W positivity diagnostic need it regardless of which normalizer drives the ceiling).
tte_norm <- Sys.getenv("PBWPFVC_TTE_NORM", "pfvc")
stopifnot(tte_norm %in% c("pfvc", "pfvc_age25"))
if (!"pfvc_age25" %in% names(cs))
  stop("cross_sectional lacks pfvc_age25 -- re-run script 03 (it now derives the structural normalizer).")
if (!identical(tte_norm, "pfvc"))
  message("*** TTE NORMALIZER = ", tte_norm, " (structural / age-stripped strain ceiling; positivity sensitivity) ***")
# Isolate outputs for ANY non-default normalizer or ceiling so sensitivity runs never collide
# with the primary (or each other). Mirrors the cache-file suffix in 10_tte_engine.R.
.dir_tag <- paste(c(
  if (!identical(tte_norm, "pfvc")) paste0("norm_", tte_norm),
  if (C_LOW != 11 || C_HIGH != 16) sprintf("c%g_%g", C_LOW, C_HIGH)), collapse = "_")
if (nzchar(.dir_tag)) {
  final_dir <- file.path(final_dir, .dir_tag)
  dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
  message("  sensitivity outputs -> ", final_dir)
}
rtrunc_lnorm <- function(n_needed, meanlog, sdlog, lo, hi) {
  acc <- numeric(0)
  while (length(acc) < n_needed) {
    cand <- rlnorm(max(n_needed * 2L, 1000L), meanlog, sdlog)
    cand <- cand[cand > lo & cand <= hi]; acc <- c(acc, cand) }
  acc[seq_len(n_needed)]
}
if (is_synthetic) {
  message("*** SYNTHETIC SITE: simulated survival (plumbing only; synthetic CLIF mortality is unreliable). ***")
  set.seed(20260615); n <- nrow(cs); died_h <- rbinom(n, 1L, 0.35)
  tte <- rep(NA_real_, n); tte[died_h == 1L] <- rtrunc_lnorm(sum(died_h), log(9), 0.95, 0.04, HORIZON)
  cs <- cs %>% mutate(death_day = if_else(died_h == 1L, floor(tte), NA_real_))
} else {
  cs <- cs %>% mutate(idx = as.numeric(difftime(death_dttm, recorded_dttm, units = "days")),
                      death_day = if_else(!is.na(idx) & idx >= 0 & idx <= HORIZON, floor(idx), NA_real_))
}
age_breaks <- quantile(cs$age_at_admission, c(1/3, 2/3), na.rm = TRUE)
base <- cs %>%
  filter(!is.na(pfvc), pfvc > 0, !is.na(pfvc_age25), pfvc_age25 > 0,
         !is.na(age_at_admission), !is.na(sex_category),
         !is.na(race_category), !is.na(sofa_total), !is.na(height_cm), !is.na(pbw), pbw > 0) %>%
  group_by(sex_category) %>% mutate(height_z = as.numeric(scale(height_cm))) %>% ungroup() %>%
  transmute(hospitalization_id, t0 = recorded_dttm,
            pfvc = .data[[tte_norm]],   # the CEILING normalizer (switch); downstream stays normalizer-agnostic
            pfvc_age25,                 # ALWAYS carry the structural normalizer (12 secondary CATE + positivity)
            pbw, death_day,
            age10 = age_at_admission / 10, sex_category, race_category, sofa_total,
            age_grp = cut(age_at_admission, c(-Inf, age_breaks, Inf),
                          labels = c("Young", "Middle", "Old")),
            height_grp = cut(height_z, c(-Inf, quantile(height_z, c(1/3, 2/3), na.rm = TRUE), Inf),
                             labels = c("Short", "Middle", "Tall")),
            # PBW/PFVC discordance: at fixed VT/PBW the strain-limiting contrast lives ENTIRELY
            # in the high-discordance (small-lung) patients, so this tertile is the TTE's own
            # version of the MP discordance HTE (12.E) -- the effect should concentrate in the
            # Discordant band, where positivity is also thinnest.
            disc_grp = cut(pbw / pfvc, c(-Inf, quantile(pbw / pfvc, c(1/3, 2/3), na.rm = TRUE), Inf),
                           labels = c("Concordant", "Mid", "Discordant")))

# [T9b] Structural-positivity restriction (enforces the diagnostic at 10g2 in the
# estimand). A patient whose lung is so small that even the LTVV floor
# (VT_FLOOR_MLKG mL/kg PBW) implies VT/PFVC > C_LOW can NEVER adhere to the strain-
# limiting ceiling: P(adhere)=0 by construction, a structural positivity violation no
# weight repairs. The empirical TRIM_ALPHA trim is model-based and may miss these if
# the weight model is misspecified, so exclude them directly and REPORT the share
# removed (overall + by age tertile -- the violation concentrates in older/lower-PFVC
# patients). pos_structural (10g2) then re-confirms ~0 infeasible on the analytic set.
# NOTE: the struct_excl CSV write is performed by 11.B_diagnostics.R (this file is
# write-free); the object is computed here so it reflects the pre-restriction `base`.
base <- base %>% mutate(min_vtpfvc = VT_FLOOR_MLKG * pbw / pfvc * 0.1)
struct_excl <- bind_rows(
  base %>% group_by(age_grp) %>%
    summarise(n = n(), n_infeasible = sum(min_vtpfvc > C_LOW),
              frac_infeasible = mean(min_vtpfvc > C_LOW), .groups = "drop") %>%
    mutate(age_grp = as.character(age_grp)),
  base %>% summarise(age_grp = "All", n = n(), n_infeasible = sum(min_vtpfvc > C_LOW),
                     frac_infeasible = mean(min_vtpfvc > C_LOW)))
message(sprintf("Structural positivity: excluding %d of %d patients (%.1f%%) infeasible for the strain ceiling at the %g mL/kg PBW floor.",
        sum(base$min_vtpfvc > C_LOW), nrow(base), 100 * mean(base$min_vtpfvc > C_LOW), VT_FLOOR_MLKG))
base <- base %>% filter(min_vtpfvc <= C_LOW) %>% select(-min_vtpfvc)

# =============================================================================
# 10b. Daily exposure + time-varying confounder panel
# =============================================================================
wf <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet")) %>%
  select(hospitalization_id, recorded_dttm, tidal_volume_set, fio2_set, peep_set,
         resp_rate_set, plateau_pressure_obs) %>%
  filter(!is.na(tidal_volume_set), tidal_volume_set > 0) %>%
  inner_join(base %>% select(hospitalization_id, t0, pfvc), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  mutate(vtpfvc = tidal_volume_set / pfvc * 0.1)
daily <- wf %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(vtpfvc = median(vtpfvc, na.rm = TRUE), vtpfvc_max = max(vtpfvc, na.rm = TRUE),
            fio2 = median(fio2_set, na.rm = TRUE),
            peep = median(peep_set, na.rm = TRUE), rr = median(resp_rate_set, na.rm = TRUE),
            .groups = "drop")   # vtpfvc_max = within-day PEAK strain, for the Q9 aggregation sensitivity

# Daily WORST (max) driving pressure for the [T5c] sensitivity: DP = plateau - PEEP
# from RECORDED plateaus only (plateau_pressure_obs is never forward-filled), taking
# the worst value in each 24h vent-day. Clinicians titrate VT to plateau (ARMA) and
# driving pressure (Amato NEJM 2015), so lagged worst-DP is a behaviorally-real
# driver of the deviation decision -- a time-varying confounder, run as a sensitivity.
dp_daily <- wf %>%
  filter(!is.na(plateau_pressure_obs), !is.na(peep_set),
         plateau_pressure_obs - peep_set > 0) %>%
  mutate(dp = plateau_pressure_obs - peep_set) %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(dp = max(dp, na.rm = TRUE), .groups = "drop")
message("DP panel: ", nrow(dp_daily), " patient-days with a recorded plateau, ",
        n_distinct(dp_daily$hospitalization_id), " patients")
# MAP: daily median (typical) from vitals.
vit <- read_parquet(file.path(output_dir, "cohort_vitals_clean.parquet")) %>%
  filter(vital_category == "map") %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(map = median(vital_value, na.rm = TRUE), .groups = "drop")

# SF ratio computed PER SpO2 measurement -- each SpO2 matched to the most recent FiO2
# within 4h (the canonical rolling join from script 03, fio2 across ALL modes) -- then
# reduced to the daily WORST (lowest SF = worst oxygenation, the value most likely to
# drive a tidal-volume decision). SpO2 is CLAMPED to [80,97] (the linear part of the
# oxyhemoglobin dissociation curve) rather than filtered, so a fully-oxygenated day
# keeps a high (good) SF instead of being dropped, and off-curve readings are bounded.
fio2_dt <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet")) %>%
  filter(!is.na(fio2_set)) %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  transmute(hospitalization_id, fio2_set, t = as.numeric(recorded_dttm)) %>%
  as.data.table()
setkey(fio2_dt, hospitalization_id, t)
spo2_dt <- read_parquet(file.path(output_dir, "cohort_vitals_clean.parquet")) %>%
  filter(vital_category == "spo2", !is.na(vital_value)) %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  transmute(hospitalization_id, vent_day,
            spo2_clamped = pmin(pmax(vital_value, 80), 97), t = as.numeric(recorded_dttm)) %>%
  as.data.table()
setkey(spo2_dt, hospitalization_id, t)
sf_daily <- fio2_dt[spo2_dt, roll = 4 * 3600, on = .(hospitalization_id, t)] %>%
  as_tibble() %>%
  filter(!is.na(fio2_set)) %>%
  mutate(fio2_frac = if_else(fio2_set > 1.5, fio2_set / 100, fio2_set),
         sf_pt = spo2_clamped / fio2_frac) %>%
  filter(is.finite(sf_pt)) %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(sf = min(sf_pt), .groups = "drop")   # daily WORST (lowest) SF
message("Worst-SF panel: ", nrow(sf_daily), " patient-days, ",
        n_distinct(sf_daily$hospitalization_id), " patients")
med <- read_parquet(file.path(output_dir, "cohort_meds.parquet")) %>%
  filter(med_group == "vasoactives") %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(admin_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  distinct(hospitalization_id, vent_day) %>% mutate(on_pressor = 1L)

# True extubation from the FULL IMV course (any ventilator mode), independent of
# tidal_volume_set. The old proxy (last volume-targeted day + 1) fires at the END of
# volume-control ventilation, but patients are routinely weaned onto pressure support
# before extubation -- at our sites ~28-38% of IMV patients have their last IMV day on a
# non-volume mode, so the old proxy truncated ventilation/liberation early. The IPCW
# weight still freezes at the last VOLUME-TARGETED day (no exposure info accrues on non-
# volume days); only the competing-risk LIBERATION time uses this true extubation.
imv_extub <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet")) %>%
  select(hospitalization_id, recorded_dttm, device_category) %>%
  filter(tolower(device_category) == "imv") %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  group_by(hospitalization_id) %>%
  summarise(imv_extub_day = max(vent_day) + 1L, .groups = "drop")
base <- base %>% left_join(imv_extub, by = "hospitalization_id")
message("IMV-course extubation derived for ", sum(!is.na(base$imv_extub_day)), " of ",
        nrow(base), " patients (every index-IMV patient should resolve).")

panel_full <- daily %>%
  left_join(vit, by = c("hospitalization_id", "vent_day")) %>%          # map (daily median)
  left_join(sf_daily, by = c("hospitalization_id", "vent_day")) %>%      # sf (daily worst)
  left_join(med, by = c("hospitalization_id", "vent_day")) %>%
  left_join(base, by = "hospitalization_id") %>%
  mutate(on_pressor = coalesce(on_pressor, 0L),
         keep = is.finite(vtpfvc) & is.finite(fio2) & is.finite(peep) & is.finite(rr) &
                is.finite(sf) & is.finite(map))
panel <- panel_full %>% filter(keep) %>% select(-keep)
# Extubation for the liberation endpoint = the IMV-course day (imv_extub_day, derived
# above, carried in via base); the IPCW weight-freeze point stays at the last volume-
# targeted day (arm idsum last_vent). The old panel-based extub proxy is retired ([T12]).
message("Panel: ", nrow(panel), " patient-days, ", n_distinct(panel$hospitalization_id), " patients")

# [T11] Missing-covariate diagnostic (the drop above must not be silent). ICU charting
# should make vtpfvc/fio2/peep/rr/sf/map dense, but a dropped vent-day (a) breaks the
# lag-1 spacing of the time-varying confounders (lag() is row-based, so a gap makes
# "yesterday" actually 2-3 days back) and (b) if it falls at the END of a stay, pulls
# the extubation proxy (= last observed vent-day + 1) earlier, biasing ventilation
# duration and the liberation competing-risk endpoint. Report the magnitude so it is
# auditable; a non-trivial share here is a signal to LOCF-fill rather than list-delete.
# NOTE: the panel_drop_summary CSV write is performed by 11.B_diagnostics.R (this file
# is write-free); the object is computed here.
drop_diag <- panel_full %>% group_by(hospitalization_id) %>%
  summarise(days_total = n(), days_kept = sum(keep),
            last_full = max(vent_day),
            last_kept = { k <- vent_day[keep]; if (length(k)) max(k) else NA_real_ },
            .groups = "drop")
panel_drop_summary <- tibble(
  patient_days_total     = sum(drop_diag$days_total),
  patient_days_dropped   = sum(drop_diag$days_total - drop_diag$days_kept),
  frac_days_dropped      = sum(drop_diag$days_total - drop_diag$days_kept) / sum(drop_diag$days_total),
  patients_any_drop      = sum(drop_diag$days_kept < drop_diag$days_total),
  patients_extub_shifted = sum(!is.na(drop_diag$last_kept) & drop_diag$last_kept < drop_diag$last_full),
  patients_lost_entirely = sum(drop_diag$days_kept == 0))
message(sprintf(paste0("Missing-covariate drops: %.1f%% of patient-days (%d of %d); ",
        "%d patients lose >=1 day, %d have extubation pulled earlier, %d lost entirely."),
        100 * panel_drop_summary$frac_days_dropped, panel_drop_summary$patient_days_dropped,
        panel_drop_summary$patient_days_total, panel_drop_summary$patients_any_drop,
        panel_drop_summary$patients_extub_shifted, panel_drop_summary$patients_lost_entirely))

# Daily WORST (lowest) pH for the [T5b] sensitivity -- the acidosis nadir is the gas
# most likely to drive a tidal-volume increase (permissive hypercapnia -> deviation
# from a low-VT arm). Read both gases once, then derive two daily worst-pH series:
# (1) POOLED arterial-equivalent = arterial + venous imputed as venous +
# 0.05 (venous pH runs ~0.03-0.05 below arterial); (2) ARTERIAL-ONLY (drops the
# imputed venous values, to confirm the +0.05 imputation isn't doing the work).
# Kept SEPARATE from the core panel: gas sampling is indication-driven (MNAR).
gas <- read_parquet(file.path(output_dir, "cohort_labs_clean.parquet")) %>%
  filter(lab_category %in% c("ph_arterial", "ph_venous"), !is.na(lab_value_numeric)) %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(lab_result_dttm, t0, units = "days"))),
         ph_art_eq = if_else(lab_category == "ph_venous",
                             lab_value_numeric + 0.05, lab_value_numeric)) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY)
ph_daily <- gas %>% group_by(hospitalization_id, vent_day) %>%        # pooled (art + venous+0.05)
  summarise(ph = min(ph_art_eq, na.rm = TRUE), .groups = "drop")      # worst (lowest) pH
ph_daily_art <- gas %>% filter(lab_category == "ph_arterial") %>%     # arterial only
  group_by(hospitalization_id, vent_day) %>%
  summarise(ph = min(lab_value_numeric, na.rm = TRUE), .groups = "drop")  # worst (lowest) pH
message("pH panel: ", nrow(ph_daily), " patient-days (pooled), ",
        n_distinct(ph_daily$hospitalization_id), " patients; ",
        n_distinct(ph_daily_art$hospitalization_id), " with an arterial gas")

# Daily worst (lowest) PaO2 for the PF-ratio sensitivity (10i2). Oxygenation drives the
# FiO2/PEEP titration that co-determines whether a low-VT arm can be held, so lagged
# P/F is a behaviorally-real confounder of the deviation decision -- run as a SENSITIVITY
# (arterial gas is indication-driven / MNAR, like pH). PF is formed in the sensitivity
# block by joining this PaO2 to the panel's daily FiO2 (same fraction/percent handling).
pao2_daily <- read_parquet(file.path(output_dir, "cohort_labs_clean.parquet")) %>%
  filter(lab_category == "po2_arterial", !is.na(lab_value_numeric)) %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(lab_result_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(pao2 = min(lab_value_numeric, na.rm = TRUE), .groups = "drop")   # worst-of-day
message("PaO2 panel: ", nrow(pao2_daily), " patient-days, ",
        n_distinct(pao2_daily$hospitalization_id), " patients with an arterial PaO2")

# =============================================================================
# 10c. Parameterized arm builder: IPCW day-weights (lagged-confounder model)
# =============================================================================
# Deviation hazard from LAGGED confounders only (never the concurrent vtpfvc that
# DEFINES deviation). Parameterized by ceiling, grace, day-weight cap, and rule so
# the sensitivities below are one-liners.
#   rule = "simple":    any exceedance after grace is a deviation.
#   rule = "corrected": an exceedance is a deviation ONLY if NOT brought back under
#                       the ceiling by the next assessment ([T4] -- transient
#                       excursions that the clinician corrects are not deviations).
trunc_w <- function(w) { q <- quantile(w, WT_TRUNC, na.rm = TRUE); pmin(pmax(w, q[1]), q[2]) }
ess_frac <- function(w) { w <- trunc_w(w); (sum(w)^2 / sum(w^2)) / length(w) }

# pnl: the panel to build on (default = global `panel`). conf: name of an extra daily
#      column to add as a LAGGED confounder; eligibility is restricted to days whose
#      lagged value is recorded ([T5b] pH, [T5c] driving pressure). conf_in_model
#      toggles whether it enters the IPCW denominator -- FALSE gives the same-day-set
#      base (so base vs adjusted differ ONLY by the confounder term, not the day-set).
arm_build <- function(ceiling, grace = GRACE, cap = DAYW_CAP, rule = "simple",
                      pnl = panel, conf = NULL, conf_in_model = TRUE,
                      keep_pday = FALSE, num_spec = "time_only",
                      trim = TRIM_ALPHA, deesc_frac = DEESC_FRAC, sf_term = "l_sf",
                      expo = "vtpfvc") {
  # sf_term: the lagged-oxygenation term in the deviation DENOMINATOR model. Default "l_sf"
  # (linear) = the primary. 11_sensitivities.R passes a richer spec (ns(l_sf,3)+l_sf:disc_grp)
  # to address the residual stratum-varying l_sf imbalance the 11.X TV-balance found. Additive:
  # default reproduces the primary den_rhs exactly.
  # expo: the daily exposure column the ceiling is defined on. Default "vtpfvc" (the primary
  # VT/PFVC strain ceiling). 11.M passes a mechanical-power column ("mp_pfvc" / "mp_pbw") so
  # the same clone-censor-weight machinery emulates a power-ceiling trial; its lag enters the
  # deviation model as l_expo in place of the lagged strain.
  if (!expo %in% names(pnl)) stop("arm_build: exposure column '", expo, "' not in the panel")
  p <- pnl %>% group_by(hospitalization_id) %>% arrange(vent_day) %>%
    mutate(above = vent_day > grace & .data[[expo]] > ceiling, lead_vt = lead(.data[[expo]]),
           # deviation rule. "simple": any post-grace exceedance. "corrected":
           # exceedance not brought back under by next assessment ([T4]).
           # "deescalate": the de-escalation MTP ([T10]) -- above-ceiling day where
           # the clinician did NOT cut next-day VT by >= deesc_frac (no next day =
           # no decision = not a deviation). This is the dynamic de-escalation decision (day-level overlap).
           viol = if (rule == "corrected") {
                    above & (is.na(lead_vt) | lead_vt > ceiling)
                  } else if (rule == "deescalate") {
                    above & !is.na(lead_vt) & lead_vt > .data[[expo]] * (1 - deesc_frac)
                  } else above,
           prior_dev = lag(cumsum(viol), default = 0) > 0,
           l_expo = lag(.data[[expo]]), l_fio2 = lag(fio2), l_peep = lag(peep),
           l_rr = lag(rr), l_sf = lag(sf), l_map = lag(map), l_pressor = lag(on_pressor)) %>%
    ungroup()
  if (!is.null(conf)) p <- p %>% group_by(hospitalization_id) %>% arrange(vent_day) %>%
    mutate(l_conf = lag(.data[[conf]])) %>% ungroup()
  fr <- p %>% filter(!prior_dev, vent_day > grace, !is.na(l_expo), !is.na(l_sf), !is.na(l_map))
  if (!is.null(conf)) fr <- fr %>% filter(!is.na(l_conf))
  # Stabilizing numerator. PRIMARY is "time_only": V (age/sex/race/SOFA) is in the
  # DENOMINATOR only, so the weights balance V and the marginal outcome MSM (no V) is
  # valid by construction -- this is the fix for the V-in-numerator/not-in-MSM mismatch
  # the balance diagnostic exposed (stabilized "baseline" weights leave V imbalanced,
  # which then requires V in the MSM). "baseline" (V in the numerator) is retained for
  # the 10k consistency check. ESS is near-identical between the two (10k), so time_only
  # costs ~nothing in efficiency while keeping subgroups/bootstrap on the cheap marginal MSM.
  num_rhs <- if (identical(num_spec, "time_only")) "ns(vent_day, 3)" else
    "ns(vent_day, 3) + age10 + sex_category + race_category + sofa_total"
  num <- glm(as.formula(paste("viol ~", num_rhs)), data = fr, family = binomial)
  den_rhs <- paste("ns(vent_day, 3) + l_expo + l_fio2 + l_peep + l_rr +", sf_term, "+ l_map +",
                   "l_pressor + age10 + sex_category + race_category + sofa_total",
                   if (!is.null(conf) && conf_in_model) "+ l_conf" else "")
  den <- glm(as.formula(paste("viol ~", den_rhs)), data = fr, family = binomial)
  p <- p %>% mutate(
    elig  = !prior_dev & vent_day > grace & !is.na(l_expo) & !is.na(l_sf) & !is.na(l_map) &
            (if (!is.null(conf)) !is.na(l_conf) else TRUE),
    p_num = predict(num, newdata = ., type = "response"),
    p_den = predict(den, newdata = ., type = "response"),
    day_w = if_else(elig, pmin(pmax((1 - p_num) / (1 - p_den), 1 / cap), cap), 1),
    # common-support trim ([T9]): eligible day in the near-certain-deviation tail
    # (modeled P(adhere) = 1 - p_den below `trim`). The clone is censored at the
    # first such day, restricting the primary estimand to the overlap region.
    out_support = elig & (1 - p_den) < trim) %>%
    group_by(hospitalization_id) %>% arrange(vent_day) %>%
    mutate(cumw = cumprod(day_w)) %>% ungroup()
  idsum <- p %>% group_by(hospitalization_id) %>% arrange(vent_day) %>%
    summarise(dev_day = { w <- which(viol); if (length(w)) min(vent_day[w]) else Inf },
              trim_day = { w <- which(out_support); if (length(w)) min(vent_day[w]) else Inf },
              last_vent = max(vent_day), death_day = first(death_day), ipcw_term = last(cumw),
              age10 = first(age10), sex_category = first(sex_category),
              race_category = first(race_category), age_grp = first(age_grp),
              height_grp = first(height_grp), disc_grp = first(disc_grp),
              imv_extub_day = first(imv_extub_day),
              .groups = "drop")
  wday <- p %>% select(hospitalization_id, vent_day, cumw) %>%
    group_by(hospitalization_id) %>%
    complete(vent_day = full_seq(c(0, max(vent_day)), 1)) %>%
    fill(cumw, .direction = "down") %>% mutate(cumw = coalesce(cumw, 1)) %>% ungroup()
  out <- list(idsum = idsum, wday = wday)
  # Per-day predicted adherence for the empirical positivity scan ([T7]): on
  # eligible days, p_adhere = 1 - p_den is the model's probability the clone stays
  # on protocol; a mass of near-zero values = a near positivity violation. Kept
  # only for the primary design (avoids bloating the bootstrap/sensitivity calls).
  if (keep_pday) out$pday <- p %>%
    filter(elig) %>%
    transmute(hospitalization_id, vent_day, age_grp, height_grp, disc_grp,
              sex_category, race_category, p_adhere = 1 - p_den)
  out
}

# =============================================================================
# 10d. Long panel + IPC-weighted MSM + competing-risk liberation, by design
# =============================================================================
make_long <- function(b, arm_lab) {
  b$idsum %>%
    mutate(mort_cens = pmin(dev_day, trim_day, HORIZON),
           event_id  = is.finite(death_day) & death_day <= mort_cens,
           fu_end    = pmax(pmin(ifelse(event_id, death_day, mort_cens), HORIZON), 1)) %>%
    transmute(hospitalization_id, fu_end, event_id, last_vent, age_grp, height_grp, disc_grp,
              sex_category, race_category, arm = arm_lab) %>%
    cross_join(tibble(day = 1:HORIZON)) %>% filter(day <= fu_end) %>%
    mutate(died = as.integer(event_id & day == fu_end), vent_day = pmin(day, last_vent)) %>%
    left_join(b$wday, by = c("hospitalization_id", "vent_day")) %>%
    mutate(ipcw = coalesce(cumw, 1)) %>% select(-cumw)
}
build_design <- function(c_low, c_high, grace = GRACE, cap = DAYW_CAP, rule = "simple",
                         pnl = panel, conf = NULL, conf_in_model = TRUE,
                         keep_pday = FALSE, num_spec = "time_only",
                         trim = TRIM_ALPHA, deesc_frac = DEESC_FRAC, sf_term = "l_sf",
                         expo = "vtpfvc") {
  # expo: length 1 (both arms on one exposure column, different ceilings -- the primary) or
  # length 2 (c(low-arm column, high-arm column) -- a head-to-head of ceilings defined on two
  # different normalizations of the same quantity, as in 11.M's MP/PFVC vs MP/PBW).
  # sf_term likewise length 1 (shared) or 2 (per arm): 11.M gives each arm the SAME lagged-history
  # covariate set minus the arm's own lagged exposure (which arm_build adds as l_expo), so the two
  # arms' deviation models span one column space and differ only through the deviation indicator.
  expo <- rep_len(expo, 2); sf_term <- rep_len(sf_term, 2)
  bl <- arm_build(c_low, grace, cap, rule, pnl, conf, conf_in_model, keep_pday, num_spec, trim, deesc_frac, sf_term[1], expo = expo[1])
  bh <- arm_build(c_high, grace, cap, rule, pnl, conf, conf_in_model, keep_pday, num_spec, trim, deesc_frac, sf_term[2], expo = expo[2])
  long <- bind_rows(make_long(bl, "strain_limiting"), make_long(bh, "permissive")) %>%
    mutate(arm = factor(arm, levels = c("permissive", "strain_limiting")))
  lib <- bind_rows(bl$idsum %>% mutate(arm = "strain_limiting"),
                   bh$idsum %>% mutate(arm = "permissive")) %>%
    mutate(dd = ifelse(is.finite(death_day), death_day, Inf),
           # liberation = true IMV-course extubation; pmax guards the rare clone whose
           # last volume-targeted day extends past its last imv-labelled day (set VT under
           # a non-imv device label) so extubation is never before the weighted course.
           extub_day = pmax(imv_extub_day, last_vent + 1L),
           cr_end = pmin(dd, extub_day, dev_day, trim_day, HORIZON),
           cr_status = case_when(dd == cr_end ~ 1L, extub_day == cr_end & extub_day < dd ~ 2L, TRUE ~ 0L),
           cr_time = pmax(cr_end, 1))
  list(long = long, lib = lib, bl = bl, bh = bh)
}
ci_curve <- function(dat) {
  fit <- suppressWarnings(glm(died ~ arm * ns(day, 4), data = dat, family = binomial, weights = ipcw))
  haz <- function(a) predict(fit, tibble(arm = factor(a, c("permissive", "strain_limiting")),
                                         day = 1:HORIZON), type = "response")
  tibble(day = 1:HORIZON, strain_limiting = 1 - cumprod(1 - haz("strain_limiting")),
         permissive = 1 - cumprod(1 - haz("permissive")))
}
rd_from <- function(dat) {
  cc <- ci_curve(dat); rs <- cc$strain_limiting[HORIZON]; rp <- cc$permissive[HORIZON]
  c(risk_sl = unname(rs), risk_pm = unname(rp), rd = unname(rs - rp))
}
cif_lib <- function(df, day = 28L) {
  if (sum(df$cr_status == 2L) < 5) return(NA_real_)
  df <- df %>% mutate(w = trunc_w(ipcw_term),
                      st = factor(cr_status, c(0L, 1L, 2L), c("censored", "death", "liberation")))
  fit <- survfit(Surv(cr_time, st) ~ 1, data = df, weights = w)
  if (!"liberation" %in% colnames(fit$pstate)) return(NA_real_)
  fit$pstate[max(which(fit$time <= day)), "liberation"]
}
lib_diff <- function(df) cif_lib(df %>% filter(arm == "strain_limiting")) -
                         cif_lib(df %>% filter(arm == "permissive"))

# per-(var, level) subgroup RD, used for point estimates AND inside the bootstrap
sub_vars <- c(age_grp = "Age tertile", height_grp = "Height tertile (w/in sex)",
              disc_grp = "PBW/PFVC discordance tertile",
              sex_category = "Sex", race_category = "Race")
sg_rd <- function(long) {
  imap_dfr(sub_vars, function(lab, v) {
    map_dfr(as.character(na.omit(unique(long[[v]]))), function(g) {
      d <- long %>% filter(as.character(.data[[v]]) == g)
      n <- n_distinct(d$hospitalization_id); if (n < 100) return(NULL)
      tibble(subgroup = lab, level = g, key = paste(lab, g, sep = "||"),
             rd = unname(rd_from(d)["rd"]), n = n)
    })
  })
}

# =============================================================================
# 10e (design-build portion). PRIMARY design (the bootstrap itself lives in 11.A).
# =============================================================================
des <- build_design(C_LOW, C_HIGH, GRACE, DAYW_CAP, "simple", keep_pday = TRUE)
long_all <- des$long; lib_all <- des$lib
point   <- rd_from(long_all); lib_pt <- lib_diff(lib_all)
sg_point <- sg_rd(long_all)
ids <- unique(long_all$hospitalization_id)
