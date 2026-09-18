# =============================================================================
# Script 15: first-stage check for the tidal-volume rounding ("sawtooth") instrument
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Standalone and exploratory; NOT in 00_run_pipeline.R. Reads the pre-gate per-timepoint
# table script 03 saves (analysis_all_eligible_timepoints). Writes aggregates only.
#
# THE IDEA. Clinicians set tidal volume in round numbers (400, 450, 500 mL), while PBW and
# PFVC are smooth in height. If a clinician aims at tau mL/kg PBW and rounds the volume to a
# grid of m mL, the delivered dose is
#     Z(PBW) = round(tau x PBW / m) x m / PBW,
# a sawtooth in PBW, and so in height within sex. Two patients a few centimetres apart, on
# either side of a rounding boundary, receive doses that differ by up to m / (tau x PBW)
# (about 10% for m = 50 mL at 450 mL) with the same age, sex, and race. After a smooth
# function of height, what remains of Z is variation in strain that owes nothing to
# demographics or to the overlap of any treatment decision. Any confounder that is smooth in
# height is orthogonal to it.
#
# WHAT THIS SCRIPT DECIDES: whether that variation exists and is strong enough to carry an
# instrumental-variable analysis. It does NOT fit the second stage, the reduced form, or any
# mortality model. Reads, in order:
#   A. Heaping of set VT: the share of volumes on multiples of 50, 25, 10, and 5 mL against
#      the share expected with no rounding, and the most common volumes.
#   B. Heaping of HEIGHT, the main threat. Z is a function of PBW, so if height is recorded
#      coarsely the sawtooth has few distinct phases and a spline in height can absorb it.
#      The read is the number of distinct observed heights inside one sawtooth period
#      (period = m / tau kg PBW = m / (tau x 0.9055) cm of height, Devine slope). Rule of
#      thumb: >= 3 distinct heights per period survives (whole inches at m = 50); ~2 is
#      marginal (5-cm recording); <= 1 is dead (10-cm recording).
#   C. Concentration of the mL/kg target. The boundaries are sharp only if one target
#      dominates. If 6, 7, and 8 mL/kg are chosen about equally often, the sawtooths
#      superimpose out of phase and the instrument blurs.
#   D. First stage, over a grid of tau in {6, 6.5, 7, 7.5, 8} and m in {50, 25, 10}:
#        VT/PBW  ~ Z + smooth(height) + sex + ns(age, 4) + race        (HC1; one row/patient)
#        VT/PFVC ~ the same                                            (the eventual exposure)
#      Smooth arms: ns(height, df) for df 3..6, and ns(height, 4) x sex. If the partial F
#      dies as the smooth gets richer, the sawtooth is not separable from height.
#      Strength is reported in strain units, not F alone (with n ~ 10^4, F clears 10 on a
#      partial R^2 of 0.001): the SD of the Z-driven component of VT/PFVC, in % of predicted
#      FVC, beside the SD of VT/PFVC explained by demographics. That ratio decides power.
#   E. Placebo grids m in {37, 43} mL, volumes no clinician rounds to. A real first stage is
#      much weaker there. A large F on a placebo grid means the height smooth is underfitting,
#      not that the instrument works.
#   F. Balance: each baseline covariate (age, SOFA, SF ratio, sex, race) regressed on Z with
#      the same smooth, in covariate-SD units per SD of residual Z. Z must not predict who
#      arrives sicker. The sex row is uninformative by construction: Devine shifts PBW by
#      4.5 kg between sexes, about half a 50-mL period, so with sex dropped from the controls
#      residual Z carries sex mechanically. Read age, SOFA, SF ratio, and race.
#
# COHORT. The pre-gate table, NOT the analytic cohort: the 6-8 mL/kg gate selects on the
# endogenous variable and truncates the sawtooth. One row per patient, the first complete
# hypoxemic timepoint (SF < 315; the band scan's index rule) with no VT/PBW gate. PRIMARY:
# volume-control timepoints only (mode_category = "assist control-volume control", script 03's
# VCV class), where set VT is the clinician's rounded choice. SENSITIVITY: all modes.
#
# WHAT A STRONG FIRST STAGE WOULD AND WOULD NOT BUY.
#   * It identifies the effect of strain WITHIN demographic strata. It does not by itself
#     identify the normalizer contrast (PBW vs PFVC): that needs the per-mL effect compared
#     across strata against 1/PFVC versus 1/PBW, which will need pooling.
#   * The exclusion threat specific to this design: units may differ in which grid they round
#     to, so a grid-specific Z can pick up unit-level mortality. The eventual 2SLS needs unit
#     fixed effects; this check does not. The same concern motivates a second instrument form,
#     snapping tau x PBW to the empirically most common volumes rather than to one grid; it is
#     left for the second stage, because choosing the heaps from the VT distribution uses the
#     exposure itself.
#   * Synthetic CLIF does not heap VT (checked 2026-09-18: share on multiples of 50 mL is 0.84x
#     the no-rounding expectation, heights continuous, no dominant target), so on synthetic data
#     this script is a NEGATIVE CONTROL. It behaved as one: across all 500 fitted cells the
#     maximum partial F was 7.2, with 22 cells above 3.84 and 5 above 6.63, against 25 and 5
#     expected under the null. Its binned-dose panel falls steeply with height because synthetic
#     VT is roughly constant in mL; real data dosed at 6-8 mL/kg PBW should be much flatter.
#
# READING RULE. Under no instrument the partial F is chi-square(1): about 1 in 20 cells exceeds
# 3.84 by chance, so across the 25-cell grid one or two "significant" cells are expected (the
# synthetic 37 mL placebo reached 4.2). Read the instrument as present only if (i) the 50 or 25 mL
# rounding grids at the modal target sit far above every placebo cell (F well above 10, not above
# 3.84), (ii) F survives the ns(height, 6) and height x sex smooths, and (iii) the implied SE
# of the eventual 2SLS coefficient is small enough to matter. The script prints it:
#     SE per point of VT/PFVC ~ sd(death) / (SD of the Z-driven component x sqrt(n)).
# At MIMIC scale (n ~ 12,000, mortality ~ 0.27) a 0.5-point SE needs a Z-driven SD near
# 0.8 % of predicted FVC; synthetic gave 0.03 to 0.06. A single strong
# cell off the modal target, or a placebo cell comparable to the rounding grids, means no instrument.
#
# Outputs (aggregates only; any row describing fewer than MIN_CELL patients is dropped):
#   final/sawtooth_heaping_{site}.csv       VT and height heaping, target concentration
#   final/sawtooth_first_stage_{site}.csv   one row per cohort x tau x m x smooth x outcome
#   final/sawtooth_balance_{site}.csv       covariate balance on Z
#   final/sawtooth_{site}.pdf               VT heaping, VT/PBW target, F heatmap, binned sawtooth
# Usage: Rscript code/15_sawtooth_first_stage.R
#        (PBWPFVC_SAW_TAUS="6,6.5,7,7.5,8"; PBWPFVC_SAW_GRIDS="50,25,10";
#         PBWPFVC_SAW_PLACEBO="37,43")
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(arrow); library(here); library(splines)
  library(sandwich); library(patchwork)
})
rm(list = ls())
source("utils/config.R")
site_name  <- config$site_name
output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

parse_num_list <- function(env_name, default) as.numeric(strsplit(Sys.getenv(env_name, default), ",")[[1]])
TARGETS_ML_PER_KG <- parse_num_list("PBWPFVC_SAW_TAUS", "6,6.5,7,7.5,8")
ROUNDING_GRIDS_ML <- parse_num_list("PBWPFVC_SAW_GRIDS", "50,25,10")
PLACEBO_GRIDS_ML  <- parse_num_list("PBWPFVC_SAW_PLACEBO", "37,43")
SF_HYPOXEMIA_THRESHOLD <- 315
MIN_CELL <- 10
DEVINE_CM_PER_KG <- 2.54 / 2.3     # Devine: 2.3 kg PBW per inch of height -> 1.104 cm per kg
VCV_MODES <- c("assist control-volume control")   # script 03's mp_mode_class == "vcv"
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

# =============================================================================
# 1. Cohort: one row per patient, first complete hypoxemic timepoint, no VT/PBW gate
# =============================================================================
timepoints <- read_parquet(file.path(output_dir, "analysis_all_eligible_timepoints.parquet")) %>%
  filter(has_all_data, sf_ratio < SF_HYPOXEMIA_THRESHOLD, is.finite(tidal_volume_set),
         tidal_volume_set > 0, is.finite(height_cm), is.finite(pbw), is.finite(pfvc)) %>%
  mutate(mode_lower = tolower(mode_category),
         is_volume_control = mode_lower %in% VCV_MODES)

first_timepoint_per_patient <- function(tp) {
  tp %>% group_by(hospitalization_id) %>%
    slice_min(recorded_dttm, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(hospitalization_id, tidal_volume_set, vtpbw, vtpfvc, height_cm, pbw, pfvc,
              pbwpfvc = pbw / pfvc, age = age_at_admission,
              sex_category = factor(sex_category), race_category = factor(race_category),
              sofa_total, sf_ratio, deceased)
}
cohorts <- list(
  volume_control = first_timepoint_per_patient(timepoints %>% filter(is_volume_control)),
  all_modes      = first_timepoint_per_patient(timepoints))
cohort_sizes <- map_int(cohorts, nrow)
message(sprintf("Index patients: volume control %d, all modes %d (%.0f%% of all-mode index rows are VCV)",
                cohort_sizes[["volume_control"]], cohort_sizes[["all_modes"]],
                100 * mean(timepoints %>% group_by(hospitalization_id) %>%
                             slice_min(recorded_dttm, n = 1, with_ties = FALSE) %>% pull(is_volume_control))))
if (cohort_sizes[["volume_control"]] < 300)
  stop("Fewer than 300 volume-control index patients: the first stage is not estimable here. ",
       "Check the mode_category strings against CLIF mCIDE before concluding anything.")

# =============================================================================
# 2. Heaping reads (A, B, C)
# =============================================================================
on_multiple <- function(x, step, tol = 0.5) abs(x - round(x / step) * step) < tol

heaping_reads <- function(cohort, cohort_name) {
  vt <- cohort$tidal_volume_set
  vt_heaping <- tibble(grid_ml = c(50, 25, 10, 5)) %>%
    mutate(share_on_grid = map_dbl(grid_ml, ~ mean(on_multiple(vt, .x))),
           share_expected_without_rounding = 1 / grid_ml,
           excess_ratio = share_on_grid / share_expected_without_rounding,
           read = "A_vt_on_grid", value_label = paste0("multiple of ", grid_ml, " mL"))

  top_volumes <- cohort %>% count(tidal_volume_set, name = "n_patients") %>%
    filter(n_patients >= MIN_CELL) %>% arrange(desc(n_patients)) %>% slice_head(n = 15) %>%
    mutate(share = n_patients / nrow(cohort), cumulative_share = cumsum(share),
           read = "A_top_volumes", value_label = paste0(signif(tidal_volume_set, 4), " mL"))

  height_heaping <- tibble(
    read = "B_height_recording",
    value_label = c("distinct heights (0.1 cm)", "share on whole cm", "share on whole inch",
                    "share on multiple of 5 cm"),
    share_on_grid = c(n_distinct(round(cohort$height_cm, 1)),
                      mean(on_multiple(cohort$height_cm, 1, 0.05)),
                      mean(on_multiple(cohort$height_cm / 2.54, 1, 0.02)),
                      mean(on_multiple(cohort$height_cm, 5, 0.05))))

  # distinct observed heights inside one sawtooth period, centred on each sex's median height
  heights_per_period <- expand_grid(target_ml_per_kg = TARGETS_ML_PER_KG,
                                    grid_ml = c(ROUNDING_GRIDS_ML, PLACEBO_GRIDS_ML),
                                    sex_category = levels(cohort$sex_category)) %>%
    mutate(period_cm = grid_ml / target_ml_per_kg * DEVINE_CM_PER_KG,
           distinct_heights_per_period = pmap_dbl(
             list(period_cm, sex_category), function(period_cm, sex_level) {
               sex_heights <- round(cohort$height_cm[cohort$sex_category == sex_level], 1)
               centre <- median(sex_heights)
               n_distinct(sex_heights[abs(sex_heights - centre) <= period_cm / 2])
             }),
           n_patients_sex = map_int(sex_category, ~ sum(cohort$sex_category == .x)),
           read = "B_heights_per_period",
           value_label = sprintf("tau %g, m %g, %s", target_ml_per_kg, grid_ml, sex_category)) %>%
    filter(n_patients_sex >= MIN_CELL)

  target_concentration <- tibble(target_ml_per_kg = TARGETS_ML_PER_KG) %>%
    mutate(share_on_grid = map_dbl(target_ml_per_kg, ~ mean(abs(cohort$vtpbw - .x) <= 0.25)),
           n_patients = map_int(target_ml_per_kg, ~ sum(abs(cohort$vtpbw - .x) <= 0.25)),
           read = "C_share_within_0.25_of_target",
           value_label = paste0(target_ml_per_kg, " mL/kg")) %>%
    filter(n_patients >= MIN_CELL)

  bind_rows(vt_heaping, top_volumes, height_heaping, heights_per_period, target_concentration) %>%
    mutate(cohort = cohort_name, n_cohort = nrow(cohort), site = site_name, .before = 1)
}
heaping_tbl <- imap_dfr(cohorts, heaping_reads)
write_csv(heaping_tbl, file.path(final_dir, paste0("sawtooth_heaping_", site_name, ".csv")))

modal_target <- cohorts$volume_control %>%
  mutate(nearest_half = round(vtpbw * 2) / 2) %>% count(nearest_half) %>%
  filter(nearest_half %in% TARGETS_ML_PER_KG) %>% slice_max(n, n = 1, with_ties = FALSE) %>%
  pull(nearest_half)
message("Modal VT/PBW target among volume-control index patients (0.5 mL/kg grid): ", modal_target)

# =============================================================================
# 3. First stage (D) and placebo grids (E)
# =============================================================================
sawtooth_instrument <- function(pbw, target_ml_per_kg, grid_ml)
  round(target_ml_per_kg * pbw / grid_ml) * grid_ml / pbw

SMOOTH_SPECS <- c(
  height_ns3        = "ns(height_cm, df = 3) + sex_category",
  height_ns4        = "ns(height_cm, df = 4) + sex_category",
  height_ns5        = "ns(height_cm, df = 5) + sex_category",
  height_ns6        = "ns(height_cm, df = 6) + sex_category",
  height_ns4_by_sex = "ns(height_cm, df = 4) * sex_category")
DEMOGRAPHIC_CONTROLS <- "ns(age, df = 4) + race_category"
OUTCOMES <- c(vtpbw = "VT/PBW (mL/kg)", vtpfvc = "VT/PFVC (% pred FVC)")

# One regressor of interest -> partial F = (coef / HC1 SE)^2. Partial R^2 by Frisch-Waugh:
# the squared correlation of the outcome and the instrument, each residualized on the controls.
# The strain-unit read is the SD of the instrument-driven component, coef x residual Z.
first_stage_fit <- function(cohort, instrument, outcome, controls_formula) {
  model_data <- cohort %>% mutate(instrument_z = instrument)
  full_fit   <- lm(as.formula(paste(outcome, "~ instrument_z +", controls_formula)), data = model_data)
  robust_se  <- sqrt(diag(vcovHC(full_fit, type = "HC1")))[["instrument_z"]]
  coef_z     <- coef(full_fit)[["instrument_z"]]
  z_residual <- resid(lm(as.formula(paste("instrument_z ~", controls_formula)), data = model_data))
  y_residual <- resid(lm(as.formula(paste(outcome, "~", controls_formula)), data = model_data))
  tibble(coef = coef_z, se_hc1 = robust_se, partial_f = (coef_z / robust_se)^2,
         partial_r2 = cor(z_residual, y_residual)^2,
         sd_residual_instrument = sd(z_residual),
         sd_instrument_driven_component = abs(coef_z) * sd(z_residual))
}

# the benchmark the instrument-driven SD is read against: how much VT/PFVC the demographics
# (and height) explain, on the same patients
demographic_gradient_sd <- function(cohort, outcome) {
  sd(fitted(lm(as.formula(paste(outcome, "~ ns(height_cm, df = 4) * sex_category +",
                                DEMOGRAPHIC_CONTROLS)), data = cohort)))
}

first_stage_grid <- expand_grid(cohort = names(cohorts), target_ml_per_kg = TARGETS_ML_PER_KG,
                                grid_ml = c(ROUNDING_GRIDS_ML, PLACEBO_GRIDS_ML),
                                smooth = names(SMOOTH_SPECS), outcome = names(OUTCOMES))
message("Fitting ", nrow(first_stage_grid), " first-stage regressions")
first_stage_tbl <- first_stage_grid %>%
  mutate(fit = pmap(list(cohort, target_ml_per_kg, grid_ml, smooth, outcome),
                    function(cohort, target_ml_per_kg, grid_ml, smooth, outcome) {
                      d <- cohorts[[cohort]]
                      first_stage_fit(d, sawtooth_instrument(d$pbw, target_ml_per_kg, grid_ml), outcome,
                                      paste(SMOOTH_SPECS[[smooth]], "+", DEMOGRAPHIC_CONTROLS))
                    })) %>%
  unnest(fit) %>%
  group_by(cohort, outcome) %>%
  mutate(sd_outcome = sd(cohorts[[first(cohort)]][[first(outcome)]]),
         sd_demographic_gradient = demographic_gradient_sd(cohorts[[first(cohort)]], first(outcome)),
         sd_death = sd(cohorts[[first(cohort)]]$deceased)) %>%
  ungroup() %>%
  mutate(instrument_to_demographic_sd_ratio = sd_instrument_driven_component / sd_demographic_gradient,
         # approximate SE of the eventual 2SLS mortality coefficient, per unit of the outcome
         # column (per point of VT/PFVC for vtpfvc): sd(death) / (SD of Z-driven part x sqrt(n))
         implied_2sls_se = sd_death / (sd_instrument_driven_component * sqrt(cohort_sizes[cohort])),
         grid_type = if_else(grid_ml %in% PLACEBO_GRIDS_ML, "placebo", "rounding"),
         is_modal_target = target_ml_per_kg == modal_target,
         is_primary_cell = cohort == "volume_control" & is_modal_target & grid_ml == 50 &
           smooth == "height_ns4_by_sex",
         outcome_label = unname(OUTCOMES[outcome]),
         period_cm = grid_ml / target_ml_per_kg * DEVINE_CM_PER_KG,
         n_patients = cohort_sizes[cohort], site = site_name)
write_csv(first_stage_tbl, file.path(final_dir, paste0("sawtooth_first_stage_", site_name, ".csv")))

# the question in one line per grid: does the richest smooth still leave a first stage, and does
# the placebo grid have none?
headline <- first_stage_tbl %>%
  filter(cohort == "volume_control", is_modal_target, smooth %in% c("height_ns4_by_sex", "height_ns6")) %>%
  select(outcome, grid_ml, grid_type, smooth, partial_f, partial_r2, sd_instrument_driven_component,
         instrument_to_demographic_sd_ratio, implied_2sls_se) %>%
  arrange(outcome, grid_type, grid_ml, smooth)
message(sprintf("\nFirst stage at the modal target (%g mL/kg), volume-control patients:", modal_target))
print(as.data.frame(headline %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)

# =============================================================================
# 4. Balance (F): Z must not predict baseline severity or demographics
# =============================================================================
BALANCE_COVARIATES <- c(age = "age", sofa_total = "SOFA", sf_ratio = "SF ratio",
                        male = "male sex", race_black = "Black race")
balance_fit <- function(cohort, target_ml_per_kg, grid_ml) {
  d <- cohort %>%
    mutate(instrument_z = sawtooth_instrument(pbw, target_ml_per_kg, grid_ml),
           male = as.numeric(tolower(sex_category) == "male"),
           race_black = as.numeric(toupper(race_category) %in% c("BLACK", "BLACK OR AFRICAN AMERICAN")))
  map_dfr(names(BALANCE_COVARIATES), function(covariate) {
    # a covariate cannot sit on both sides: drop age's spline, sex, or race from the controls
    controls <- c(height = "ns(height_cm, df = 4)",
                  sex    = if (covariate != "male") "sex_category",
                  age    = if (covariate != "age") "ns(age, df = 4)",
                  race   = if (covariate != "race_black") "race_category")
    controls <- paste(controls[lengths(controls) > 0], collapse = " + ")
    if (sd(d[[covariate]], na.rm = TRUE) == 0)
      return(tibble(covariate = covariate, std_coef = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
    fit <- lm(as.formula(paste(covariate, "~ instrument_z +", controls)), data = d)
    z_sd <- sd(resid(lm(as.formula(paste("instrument_z ~", controls)), data = d)))
    coef_z <- coef(fit)[["instrument_z"]]
    se_z <- sqrt(diag(vcovHC(fit, type = "HC1")))[["instrument_z"]]
    scale <- z_sd / sd(d[[covariate]], na.rm = TRUE)   # covariate SDs per SD of residual Z
    tibble(covariate = covariate, std_coef = coef_z * scale,
           lo = (coef_z - 1.96 * se_z) * scale, hi = (coef_z + 1.96 * se_z) * scale,
           p = 2 * pnorm(-abs(coef_z / se_z)))
  })
}
balance_tbl <- expand_grid(cohort = names(cohorts), target_ml_per_kg = TARGETS_ML_PER_KG,
                           grid_ml = ROUNDING_GRIDS_ML) %>%
  mutate(fit = pmap(list(cohort, target_ml_per_kg, grid_ml),
                    function(cohort, target_ml_per_kg, grid_ml)
                      balance_fit(cohorts[[cohort]], target_ml_per_kg, grid_ml))) %>%
  unnest(fit) %>%
  mutate(covariate_label = unname(BALANCE_COVARIATES[covariate]),
         is_modal_target = target_ml_per_kg == modal_target,
         n_patients = cohort_sizes[cohort], site = site_name)
write_csv(balance_tbl, file.path(final_dir, paste0("sawtooth_balance_", site_name, ".csv")))

# =============================================================================
# 5. Figure
# =============================================================================
vc <- cohorts$volume_control
theme_set(theme_minimal(base_size = 10))

vt_residue <- vc %>% mutate(residue_mod_50 = round(tidal_volume_set) %% 50) %>%
  count(residue_mod_50) %>% filter(n >= MIN_CELL)
panel_vt <- ggplot(vt_residue, aes(residue_mod_50, n)) +
  geom_col(fill = okabe[5], width = 0.9) +
  labs(title = "A. Set VT modulo 50 mL", subtitle = "a spike at 0 is rounding to 50 mL",
       x = "VT mod 50 (mL)", y = "patients")

target_hist <- vc %>% mutate(bin = floor(vtpbw * 4) / 4) %>% count(bin) %>%
  filter(n >= MIN_CELL, bin >= 3, bin <= 12)
panel_target <- ggplot(target_hist, aes(bin + 0.125, n)) +
  geom_col(fill = okabe[3], width = 0.24) +
  geom_vline(xintercept = modal_target, linetype = 2) +
  labs(title = "C. Delivered VT/PBW", subtitle = sprintf("modal target %g mL/kg (dashed)", modal_target),
       x = "VT/PBW (mL/kg)", y = "patients")

heat <- first_stage_tbl %>%
  filter(cohort == "volume_control", smooth == "height_ns4_by_sex", outcome == "vtpfvc") %>%
  mutate(grid_label = factor(paste0(grid_ml, " mL", if_else(grid_type == "placebo", " (placebo)", "")),
                             levels = paste0(c(ROUNDING_GRIDS_ML, PLACEBO_GRIDS_ML), " mL",
                                             c(rep("", length(ROUNDING_GRIDS_ML)),
                                               rep(" (placebo)", length(PLACEBO_GRIDS_ML))))))
panel_heat <- ggplot(heat, aes(factor(target_ml_per_kg), grid_label, fill = log10(pmax(partial_f, 0.01)))) +
  geom_tile(colour = "white") +
  geom_text(aes(label = signif(partial_f, 2)), size = 3, colour = "white") +
  scale_fill_viridis_c(name = "log10 partial F") +
  labs(title = "D. First-stage partial F on VT/PFVC",
       subtitle = "ns(height, 4) x sex + ns(age, 4) + race; placebo rows should be dark",
       x = "target tau (mL/kg PBW)", y = "rounding grid m")

# binned sawtooth: mean delivered VT/PBW by 1-cm height bin within sex, against the modal-target
# Z curve. If rounding drives the dose, the bin means trace the teeth.
binned <- vc %>% mutate(height_bin = floor(height_cm) + 0.5) %>%
  group_by(sex_category, height_bin) %>%
  summarise(n = n(), mean_vtpbw = mean(vtpbw), mean_pbw = mean(pbw), .groups = "drop") %>%
  filter(n >= MIN_CELL) %>%
  mutate(z_curve = sawtooth_instrument(mean_pbw, modal_target, 50))
panel_saw <- ggplot(binned, aes(height_bin)) +
  geom_line(aes(y = z_curve, colour = "Z: rounded to 50 mL"), linewidth = 0.4) +
  geom_point(aes(y = mean_vtpbw, colour = "observed mean", size = n), alpha = 0.8) +
  facet_wrap(~ sex_category, scales = "free_x") +
  scale_colour_manual(values = c("Z: rounded to 50 mL" = okabe[6], "observed mean" = okabe[5]), name = NULL) +
  scale_size_area(max_size = 3, guide = "none") +
  labs(title = sprintf("E. Delivered dose by height against the %g mL/kg sawtooth", modal_target),
       subtitle = "1-cm bins with >= 10 patients", x = "height (cm)", y = "VT/PBW (mL/kg)")

fig <- (panel_vt + panel_target) / panel_heat / panel_saw +
  plot_layout(heights = c(1, 1.1, 1.1)) +
  plot_annotation(title = sprintf("Rounding-sawtooth instrument: first-stage check (%s)", site_name),
                  subtitle = sprintf("volume-control index patients, n = %d; pre-gate cohort (no VT/PBW band)",
                                     cohort_sizes[["volume_control"]]))
ggsave(file.path(final_dir, paste0("sawtooth_", site_name, ".pdf")), fig, width = 11, height = 13)
message("15_sawtooth_first_stage complete -> ", final_dir)
