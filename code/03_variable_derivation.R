# =============================================================================
# Script 03: Variable Derivation
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Derives every analysis variable, applies the analytic inclusion gates, and picks
# one index timepoint per patient: the cross-sectional cohort of figures 1-3 and the
# baseline of figure 4.
#
# Population (PBWPFVC_COHORT, environment variable only; the cohort script 01 built):
#   imv (default; the paper's cohort)  patients with height 150-210 cm, known sex, and
#       at least one ventilator timepoint that has VT/PBW, VT/PFVC and an SF ratio,
#       is lung-protective (VT/PBW 6-8 mL/kg) AND hypoxemic (SF < 315, the Rice 2007
#       equivalent of PF 300). A patient qualifies if ANY such timepoint exists. The
#       index is the first qualifying timepoint with a driving pressure within
#       INDEX_WINDOW_HOURS of the first IMV row, else the first qualifying timepoint.
#   nosupport (the negative control)   the first room-air or nasal-cannula row with an
#       SF ratio within CONTROL_ICU_WINDOW_H of an ICU admission, with no advanced
#       support before it; patients escalated to any support within
#       ESCALATION_LANDMARK_H of that index are removed.
#   niv (built on request only)        the first HFNC / NIPPV / CPAP row with SF < 315.
# SOFA is scored once per patient by clifR's compute_sofa(), from the worst value of
# each input over the SOFA_WINDOW_H hours from the index (section 3e2).
#
# Sections: 3a PBW and PFVC (Devine; race-specific GLI-2012), 3b SF and PF ratios,
# 3d per-timepoint dose, mechanics, vasopressor and lab variables, 3e the index,
# 3e2 SOFA, 3f ventilator-free days, 3h-3i saved tables and the attrition log,
# 3k the non-hypoxemic negative-control cohorts, 3l federated PBW/PFVC distributions.
#
# Inputs (config$output_dir): the _clean tables of script 02, and cohort_demographics,
#   cohort_weights, cohort_meds, cohort_assessments, nc_cohort, cohort_icu_stays and
#   attrition_log_partial.csv from script 01.
# Outputs, patient-level (config$output_dir, never shared):
#   analysis_cross_sectional    one row per patient at the index -> 04, 05, 10, 29, supplement/
#   analysis_all_timepoints     every ventilator timepoint of those patients -> 04
#   analysis_all_eligible_timepoints  every timepoint before the VT/PBW gate -> 10
#   analysis_broad_pfvc         everyone with PBW and PFVC, before ventilation gates -> 04
#   analysis_negative_control   the non-hypoxemic cohorts -> 04 (4j)
#   ne_equiv_admin              norepinephrine-equivalent dose per administration -> 10, 21
# Outputs, aggregate (final/cross_sectional/, returned by the site):
#   attrition_log_{site}.csv    the 7-step cohort funnel -> 04 (CONSORT), pooling
#   dist_histograms_{site}.csv, dist_quantiles_{site}.csv  PBW/PFVC by demographic group
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(data.table)
library(lubridate)
library(rspiro)

source("utils/config.R")
source("utils/standardize_pressor_dose.R")
source("utils/attrition_log.R")

site_name <- config$site_name
output_dir <- config$output_dir
final_dir <- final_dir_for("cross_sectional")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

# Fixed bin edges for the federated PBW:PFVC distribution exports (must match
# across sites so per-site histograms can be summed into a pooled distribution).
PBWPFVC_BIN_EDGES <- seq(0, 50, by = 1)
AGE_BIN_EDGES     <- c(18, 30, 40, 50, 60, 70, 80, Inf)
HEIGHT_BIN_EDGES  <- c(150, 155, 160, 165, 170, 175, 180, 185, 190, 210)

# External (ARMA) constants for the DIFFERENCE parameterization of PBW vs PFVC (3d).
# Fixed, not cohort-derived, so vt_excess_ml means the same thing at every site.
VT_PER_KG_ARMA   <- 6    # mL/kg PBW, the low-VT arm's prescribed dose
VT_PCT_PFVC_ARMA <- 11   # % of predicted FVC: about the 75th percentile of VT/PFVC in the ARMA low-VT arm

# Time windows, in hours, chosen to capture early ventilation.
# FIO2_LOOKBACK_H (4 h; utils/config.R): an SpO2 or PaO2 takes the most recent FiO2
#   recorded up to this long before it.
MEASUREMENT_CARRY_H   <- 4   # an SF/PF ratio, lab value or vasopressor rate is carried
                             # forward to a ventilator timepoint for up to this long
INDEX_WINDOW_HOURS    <- 6   # ventilated index: prefer a timepoint with driving pressure
                             # within this long of the first IMV row
CONTROL_ICU_WINDOW_H  <- 6   # no-support index: first qualifying row within this long
                             # of an ICU admission
ESCALATION_LANDMARK_H <- 24  # no-support control: escalation to any support within this
                             # long of the index removes the patient
SOFA_WINDOW_H         <- 24  # SOFA from the worst values over this long from the index

# =============================================================================
# Load cleaned intermediate data
# =============================================================================

cohort_ids <- readRDS(file.path(output_dir, "cohort_hospitalization_ids.rds"))
resp_waterfall <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet")) %>%
  estimate_fio2_nosupport()   # no-support control only: FiO2 on room air / cannula
cohort_demographics <- read_parquet(file.path(output_dir, "cohort_demographics.parquet"))
cohort_vitals <- read_parquet(file.path(output_dir, "cohort_vitals_clean.parquet"))
cohort_labs <- read_parquet(file.path(output_dir, "cohort_labs_clean.parquet"))
cohort_heights <- read_parquet(file.path(output_dir, "cohort_heights_clean.parquet"))
cohort_weights <- read_parquet(file.path(output_dir, "cohort_weights.parquet"))
cohort_meds <- read_parquet(file.path(output_dir, "cohort_meds.parquet"))
cohort_assessments <- read_parquet(file.path(output_dir, "cohort_assessments.parquet"))

message("Loaded all cleaned intermediate data")

# =============================================================================
# 3a. PBW and PFVC calculation
# =============================================================================

# Join demographics with heights
pbw_pfvc_data <- cohort_demographics %>%
  select(hospitalization_id, patient_id, age_at_admission,
         sex_category, race_category, deceased) %>%
  inner_join(cohort_heights, by = "hospitalization_id") %>%
  # Height inclusion range from the original analysis (PBWvsFVC): 150-210 cm.
  # The GLI-2012 reference equations are only valid within this range, and the
  # published cohort excluded heights outside it. Restricting here cascades to all
  # downstream analyses, since they inner-join pbw_pfvc_data.
  filter(!is.na(height_cm), height_cm >= 150, height_cm <= 210) %>%
  mutate(
    # Map sex to numeric for pred_GLI: 1=male, 2=female
    sex_numeric = case_when(
      sex_category == "Male"   ~ 1L,
      sex_category == "Female" ~ 2L,
      TRUE                     ~ NA_integer_
    ),
    # GLI ethnicity: 1=Caucasian, 2=African-American, 5=Other/Mixed
    race_numeric = case_when(
      race_category == "WHITE" ~ 1L,
      race_category == "BLACK" ~ 2L,
      TRUE                     ~ 5L
    ),
    # PBW via Devine formula (kg)
    pbw = case_when(
      sex_numeric == 1 ~ 50.0 + 2.3 * (height_cm / 2.54 - 60),
      sex_numeric == 2 ~ 45.5 + 2.3 * (height_cm / 2.54 - 60)
    ),
    # PFVC via GLI-2012 (litres)
    pfvc = pred_GLI(
      age       = age_at_admission,
      height    = height_cm / 100,
      gender    = sex_numeric,
      ethnicity = race_numeric,
      param     = "FVC"
    ),
    # FVC_age25: GLI prediction with age pinned to 25 (peak adult), keeping the
    # patient's actual height/sex/race -- an age-STRIPPED, structural-size surrogate.
    # It is the "structural-only" point on the age-correction spectrum:
    #   PBW = no age correction, FVC_age25 = GLI structure but age-flat, PFVC = full.
    # VT/FVC_age25 vs VT/PFVC isolates whether PFVC's age slope earns its keep; the
    # triple is read as a BRACKET, not a truth -- the strain denominator's age behaviour
    # in the old/critically ill is unidentified (statistically, physiologically,
    # mechanically), so the analyses sweep the correction rather than pick one.
    pfvc_age25 = pred_GLI(
      age       = rep(25, n()),
      height    = height_cm / 100,
      gender    = sex_numeric,
      ethnicity = race_numeric,
      param     = "FVC"
    )
  ) %>%
  filter(!is.na(pbw), !is.na(pfvc), pfvc > 0) %>%
  # Reference categories: Male and White (first factor level). Releveling here
  # cascades to every downstream model (analysis_data inner-joins these columns),
  # so all regressions report effects relative to male / white patients.
  mutate(
    sex_category  = factor(sex_category,  levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER"))
  )

message("PBW/PFVC calculated: ", nrow(pbw_pfvc_data), " hospitalizations")

# Broad PFVC cohort for the PFVC-vs-PBW table: all eligible patients with
# height/age/sex/race/PFVC, before any ventilation-based inclusion. Saved for
# script 04.
write_parquet(pbw_pfvc_data, file.path(output_dir, "analysis_broad_pfvc.parquet"))

# =============================================================================
# 3b. Compute SF and PF ratios at measurement time using concurrent FiO2
# =============================================================================
# Ratios are computed at the time SpO2/PaO2 were measured, using the FiO2
# that was active at that moment (most recent prior FiO2 within FIO2_LOOKBACK_H).
# These pre-computed ratios are then carried forward to IMV timepoints.

# --- FiO2 from waterfall (for matching to SpO2/PaO2 measurement times) ---
fio2_data <- resp_waterfall %>%
  filter(!is.na(fio2_set)) %>%
  select(hospitalization_id, recorded_dttm, fio2_set) %>%
  mutate(fio2_dttm = as.numeric(recorded_dttm))

fio2_dt <- as.data.table(fio2_data)
setkey(fio2_dt, hospitalization_id, fio2_dttm)

# --- SF ratio: SpO2 / FiO2 at time of SpO2 measurement ---
spo2_data <- cohort_vitals %>%
  filter(vital_category == "spo2", !is.na(vital_value)) %>%
  # The SF ratio uses SpO2 80-97 only, the original analysis's (PBWvsFVC) range.
  # Above 97 the SpO2-PaO2 relation is too flat (the curve saturates) and below 80
  # too steep for SF to track PF. Both bounds are applied here: this is the only
  # place the analytic cohort's SF is bounded, and script 01 keeps every SpO2.
  filter(vital_value >= 80, vital_value <= 97) %>%
  select(hospitalization_id, recorded_dttm, vital_value) %>%
  rename(spo2_value = vital_value, spo2_dttm = recorded_dttm) %>%
  mutate(spo2_dttm_num = as.numeric(spo2_dttm))

spo2_dt <- as.data.table(spo2_data)
setkey(spo2_dt, hospitalization_id, spo2_dttm_num)

# Rolling join: for each SpO2, find most recent prior FiO2 within FIO2_LOOKBACK_H.
# The joined frame's fio2_dttm column holds the SpO2 time (a rolling join keeps the
# i-side key), so sf_dttm is when the SpO2 was measured, not when FiO2 was set.
sf_joined <- fio2_dt[spo2_dt,
                      roll = FIO2_LOOKBACK_H * 3600,
                      on = .(hospitalization_id, fio2_dttm = spo2_dttm_num)]

sf_joined <- as_tibble(sf_joined) %>%
  filter(!is.na(fio2_set), !is.na(spo2_value)) %>%
  mutate(sf_ratio = spo2_value / fio2_set,
         sf_dttm = fio2_dttm)  # timestamp = when SpO2 was measured

message("SF ratio computed at ", nrow(sf_joined), " SpO2 measurement times")

# --- PF ratio: PaO2 / FiO2 at time of PaO2 measurement ---
pao2_data <- cohort_labs %>%
  filter(lab_category == "po2_arterial", !is.na(lab_value_numeric)) %>%
  select(hospitalization_id, lab_result_dttm, lab_value_numeric) %>%
  rename(pao2_value = lab_value_numeric, pao2_dttm = lab_result_dttm) %>%
  mutate(pao2_dttm_num = as.numeric(pao2_dttm))

pao2_dt <- as.data.table(pao2_data)
setkey(pao2_dt, hospitalization_id, pao2_dttm_num)

pf_joined <- fio2_dt[pao2_dt,
                      roll = FIO2_LOOKBACK_H * 3600,
                      on = .(hospitalization_id, fio2_dttm = pao2_dttm_num)]

pf_joined <- as_tibble(pf_joined) %>%
  filter(!is.na(fio2_set), !is.na(pao2_value)) %>%
  mutate(pf_ratio = pao2_value / fio2_set,
         pf_dttm = fio2_dttm)  # timestamp = when PaO2 was measured

message("PF ratio computed at ", nrow(pf_joined), " PaO2 measurement times")

# =============================================================================
# 3c. SOFA is scored after the index is chosen (section 3e2): the worst values over
# the 24 hours from each patient's index, by clifR's definition.
# =============================================================================

# =============================================================================
# 3d. Derive analysis variables
# =============================================================================



# Build per-timepoint dataset from the waterfall output (IMV rows only; for the
# nosupport and niv cohorts, their own device rows, with every ventilator setting
# blanked so no dose or pressure variable is derived)
spine_devices <- switch(config$cohort, imv = "imv", niv = NIV_DEVICES, nosupport = NOSUPPORT_DEVICES)
imv_timepoints <- resp_waterfall %>%
  filter(tolower(device_category) %in% spine_devices) %>%
  select(hospitalization_id, recorded_dttm, device_category, mode_category,
         fio2_set, tidal_volume_set, peep_set, plateau_pressure_obs,
         resp_rate_set, peak_inspiratory_pressure_obs,
         any_of("minute_vent_obs"))
if (config$cohort != "imv")
  imv_timepoints <- imv_timepoints %>%
    mutate(tidal_volume_set = NA_real_, plateau_pressure_obs = NA_real_, peak_inspiratory_pressure_obs = NA_real_)

# Join with PBW/PFVC data
analysis_data <- imv_timepoints %>%
  inner_join(
    pbw_pfvc_data %>% select(hospitalization_id, age_at_admission,
                              sex_category, race_category, sex_numeric,
                              race_numeric, height_cm, pbw, pfvc, pfvc_age25, deceased),
    by = "hospitalization_id"
  ) %>%
  mutate(
    pbwpfvc = pbw / pfvc,
    pbwpfvc_age25 = pbw / pfvc_age25,          # discordance vs the structural (age-flat) size
    # --- the same disagreement as a DIFFERENCE, not a ratio ---------------------
    # PBW is kg and PFVC is litres, so they cannot be subtracted directly, and scaling
    # one to the other WITHIN a cohort would make the exposure depend on the cohort.
    # Both constants below are EXTERNAL (ARMA), so the difference is comparable across
    # sites and interpretable at the bedside:
    #   VT_PER_KG_ARMA  6 mL/kg PBW  -- the protective dose the trial prescribed
    #   VT_PCT_PFVC_ARMA 11% of PFVC -- the VT/PFVC that same arm actually delivered (p75)
    # vt_excess_ml = the millilitres of tidal volume the PBW rule prescribes OVER the
    # PFVC rule at the protective dose. Positive = PBW over-doses this patient. It is the
    # ratio's difference-scale twin, but weighted differently: the ratio treats a 10%
    # mis-sizing the same in a 40 kg and a 90 kg patient, the difference does not, and
    # the difference is the quantity a clinician can act on (mL on the ventilator).
    vt_excess_ml   = VT_PER_KG_ARMA * pbw - (VT_PCT_PFVC_ARMA / 100) * pfvc * 1000,
    vtpbw = tidal_volume_set / pbw,                # mL/kg PBW
    vtpfvc = tidal_volume_set / pfvc * 0.1,        # % of predicted FVC: VT (mL) / PFVC (L) x 0.1
    vtpfvc_age25 = tidal_volume_set / pfvc_age25 * 0.1,   # structural-only normalizer (age stripped)
    # QC AT SOURCE: plateau <= PEEP is non-physiologic (a ventilated patient
    # receiving a tidal volume cannot have zero/negative driving pressure) and
    # reflects a charting/measurement error -> driving pressure is undefined (NA).
    # Doing this here means every downstream script gets clean dp (and everything
    # derived from it: crs, ers, mechanical power), instead of each having to remember
    # to filter dp > 0.
    dp = if_else(
      !is.na(plateau_pressure_obs) & !is.na(peep_set) & plateau_pressure_obs > peep_set,
      plateau_pressure_obs - peep_set,
      NA_real_
    ),
    crs = if_else(!is.na(dp) & dp > 0, tidal_volume_set / dp, NA_real_),
    # Elastance in cmH2O/L. crs is mL/cmH2O, so 1/crs is cmH2O/mL; the x1000
    # converts to cmH2O/L so downstream tables don't round elastance to 0.0.
    # ers_pbw / ers_pfvc (elastance normalized to PBW / PFVC) inherit this scale.
    ers = if_else(!is.na(crs) & crs > 0, 1000 / crs, NA_real_),
    ers_pbw = ers * pbw,
    ers_pfvc = ers * pfvc,
    ers_pfvc_age25 = ers * pfvc_age25,         # specific-elastance bracket: structural-size scaling
    # --- Mechanical power (J/min), simplified equation -------------------------
    # The simplified power equation differs by inspiratory flow pattern, which is
    # set by the ventilator mode:
    #   * Volume control (constant / square flow): a resistive correction is
    #     applied via the driving pressure --
    #       MP = 0.098 * RR * VT(L) * (Ppeak - 1/2 * (Pplat - PEEP))   [Gattinoni 2016]
    #   * Pressure-targeted modes (pressure control, PRVC; decelerating flow): the
    #     pressure-volume loop is rectangular, so there is no -1/2*dP correction --
    #       MP = 0.098 * RR * VT(L) * Ppeak                            [Becher 2019]
    # Spontaneous / unclassifiable modes (pressure support/CPAP, SIMV, etc.) are
    # left undefined (MP = NA). RR uses the set rate, consistent with the set
    # tidal volume and PEEP used elsewhere; Ppeak and Pplat are observed values
    # (never forward-filled), so MP is computable only where they were recorded.
    # NOTE: the VCV/PCV mode_category sets below are the CLIF mCIDE strings as used
    # elsewhere in the pipeline; extend them if a site reports other controlled modes.
    mp_mode_class = case_when(
      mode_category %in% c("assist control-volume control")                        ~ "vcv",
      mode_category %in% c("pressure control", "pressure-regulated volume control") ~ "pcv",
      TRUE ~ NA_character_
    ),
    mechanical_power = case_when(
      mp_mode_class == "vcv" & !is.na(resp_rate_set) & !is.na(tidal_volume_set) &
        !is.na(peak_inspiratory_pressure_obs) & !is.na(dp) & dp > 0 ~
        0.098 * resp_rate_set * (tidal_volume_set / 1000) *
          (peak_inspiratory_pressure_obs - 0.5 * dp),
      mp_mode_class == "pcv" & !is.na(resp_rate_set) & !is.na(tidal_volume_set) &
        !is.na(peak_inspiratory_pressure_obs) ~
        0.098 * resp_rate_set * (tidal_volume_set / 1000) * peak_inspiratory_pressure_obs,
      TRUE ~ NA_real_
    ),
    # --- Normalizations of mechanical power ---------------------------------------
    # Unnormalized power is meaningless across lung sizes: the same J/min would shred
    # a small lung and under-ventilate a large one. Elastic power per breath is
    # 1/2 * VT^2 * Ers; writing VT = strain * V0 and Ers = E_spec / V0 (V0 = resting
    # lung volume, E_spec = specific elastance) gives
    #     power / V0 = RR * 1/2 * strain^2 * E_spec,
    # which depends only on strain, tissue property and rate. Dividing power by a
    # PREDICTED lung volume is therefore size-invariant by construction, and PFVC is
    # the available proxy for V0. MP/PBW = MP/PFVC x (PFVC/PBW): it carries the PBW/PFVC
    # discordance and under-reads power in exactly the patients PBW over-sizes.
    #   mp_pfvc  primary size normalization (J/min per L predicted FVC)
    #   mp_pbw   the conventional (Gattinoni) normalization, the biased comparator
    #   mp_crs   power per unit MEASURED compliance = MP x Ers (cmH2O/mL): energy per
    #            aerated lung, the functional-lung comparator; ~ DP^2, the most
    #            recoil-laden of the three (04 and supplement/xsec_mortality_prediction.R
    #            use it as such)
    mp_pbw  = mechanical_power / pbw,
    mp_pfvc = mechanical_power / pfvc,
    mp_pfvc_age25 = mechanical_power / pfvc_age25,
    mp_crs  = mechanical_power * ers / 1000,
    # Elastic TIDAL component of power (J/min): the triangle 1/2 * VT * dP per breath,
    # the energy stored in the lung by the tidal breath, independent of flow pattern
    # (the mode split above concerns only the resistive/PEEP components of the total).
    # Divided by PFVC it is the specific elastic power, a sensitivity analysis for the
    # total peak-pressure power. Defined wherever dP is.
    mp_elastic      = if_else(!is.na(dp) & dp > 0 & !is.na(resp_rate_set) & !is.na(tidal_volume_set),
                              0.098 * resp_rate_set * (tidal_volume_set / 1000) * 0.5 * dp,
                              NA_real_),
    mp_elastic_pbw  = mp_elastic / pbw,
    mp_elastic_pfvc = mp_elastic / pfvc
  )

message("Analysis timepoints: ", nrow(analysis_data), " rows, ",
        n_distinct(analysis_data$hospitalization_id), " hospitalizations")

# QC diagnostic: how many timepoints had plateau <= PEEP (driving pressure NA'd).
n_dp_bad <- sum(!is.na(analysis_data$plateau_pressure_obs) &
                !is.na(analysis_data$peep_set) &
                analysis_data$plateau_pressure_obs <= analysis_data$peep_set)
if (n_dp_bad > 0) {
  n_dp_obs <- sum(!is.na(analysis_data$plateau_pressure_obs) & !is.na(analysis_data$peep_set))
  message("QC: ", n_dp_bad, " of ", n_dp_obs, " plateau-measured timepoint(s) had ",
          "plateau <= PEEP (", round(100 * n_dp_bad / n_dp_obs, 1),
          "%); driving pressure set to NA. If this fraction is large or systematic, ",
          "inspect the site's CLIF respiratory_support mapping (plateau/PEEP swap, units) ",
          "before assuming sporadic charting error.")
}

# --- NE equivalents: compute from cohort_meds and rolling-join to IMV timepoints ---
# All catecholamine doses are first standardized to mcg/kg/min; vasopressin is
# handled separately because it is dosed in units/min (units/hr), not mcg-based.
# The conversion factors below (and vasopressin's 2.5 per unit/min) are published
# norepinephrine-equivalence factors (citation in the Methods).
catecholamines <- c("norepinephrine", "epinephrine", "dopamine",
                    "phenylephrine", "dobutamine")

ne_equiv_cat <- cohort_meds %>%
  filter(med_category %in% catecholamines, !is.na(med_dose)) %>%
  standardize_pressor_dose(
    weights = cohort_weights,
    out_col = "dose_mcg_kg_min",
    label = "NE-equiv catecholamines"
  ) %>%
  filter(!is.na(dose_mcg_kg_min)) %>%
  mutate(
    ne_equiv = case_when(
      med_category == "norepinephrine" ~ dose_mcg_kg_min,
      med_category == "epinephrine"    ~ dose_mcg_kg_min,
      med_category == "dopamine"       ~ dose_mcg_kg_min / 100,
      med_category == "phenylephrine"  ~ dose_mcg_kg_min / 10,
      med_category == "dobutamine"     ~ 0
    )
  ) %>%
  select(hospitalization_id, admin_dttm, ne_equiv)

# Vasopressin: units/min -> NE-equiv uses the raw rate (standard practice),
# so we pass through without weight-based conversion. Units/hr entries are
# normalized to units/min.
ne_equiv_vaso <- cohort_meds %>%
  filter(med_category == "vasopressin", !is.na(med_dose)) %>%
  mutate(
    med_dose = as.numeric(med_dose),
    unit_clean = str_to_lower(str_replace_all(med_dose_unit, "\\s+", "")),
    units_per_min = case_when(
      unit_clean %in% c("units/min", "u/min")   ~ med_dose,
      unit_clean %in% c("units/hr", "u/hr",
                        "units/h",  "u/h")      ~ med_dose / 60,
      TRUE                                       ~ NA_real_
    ),
    ne_equiv = units_per_min * 2.5
  ) %>%
  filter(!is.na(ne_equiv)) %>%
  select(hospitalization_id, admin_dttm, ne_equiv)

ne_equiv <- bind_rows(ne_equiv_cat, ne_equiv_vaso) %>%
  filter(ne_equiv > 0) %>%
  group_by(hospitalization_id, admin_dttm) %>%
  summarise(ne_equiv_total = sum(ne_equiv, na.rm = TRUE), .groups = "drop")

# Persist the per-administration NE-equivalent table: the daily panel
# (10_panel_common.R) reduces it to a daily peak dose for the biotrauma joint
# models, and the rolling join below keeps only the value in force at each
# IMV timepoint.
write_parquet(ne_equiv, file.path(output_dir, "ne_equiv_admin.parquet"))

ne_dt <- as.data.table(ne_equiv)
ne_dt[, join_dttm := as.numeric(admin_dttm)]
setkey(ne_dt, hospitalization_id, join_dttm)

analysis_ne_dt <- as.data.table(analysis_data)
analysis_ne_dt[, join_dttm := as.numeric(recorded_dttm)]
setkey(analysis_ne_dt, hospitalization_id, join_dttm)

analysis_ne_joined <- ne_dt[, .(hospitalization_id, join_dttm, ne_equiv_total)][
  analysis_ne_dt, roll = MEASUREMENT_CARRY_H * 3600, on = .(hospitalization_id, join_dttm)
]
analysis_ne_joined[, join_dttm := NULL]

analysis_data <- as_tibble(analysis_ne_joined) %>%
  mutate(ne_equiv_total = coalesce(ne_equiv_total, 0))

message("NE equivalents joined: ", sum(analysis_data$ne_equiv_total > 0),
        " timepoints with vasopressors")

# --- Join pre-computed SF ratio to IMV timepoints (most recent prior within MEASUREMENT_CARRY_H) ---
sf_for_join <- as.data.table(sf_joined %>% select(hospitalization_id, sf_dttm, sf_ratio))
setkey(sf_for_join, hospitalization_id, sf_dttm)

analysis_dt <- as.data.table(analysis_data)
analysis_dt[, recorded_dttm_num := as.numeric(recorded_dttm)]
setkey(analysis_dt, hospitalization_id, recorded_dttm_num)

analysis_with_sf <- sf_for_join[analysis_dt,
                                 roll = MEASUREMENT_CARRY_H * 3600,
                                 on = .(hospitalization_id, sf_dttm = recorded_dttm_num)]
analysis_with_sf[, sf_dttm := NULL]

# --- Join pre-computed PF ratio to IMV timepoints (most recent prior within MEASUREMENT_CARRY_H) ---
pf_for_join <- as.data.table(pf_joined %>% select(hospitalization_id, pf_dttm, pf_ratio))
setkey(pf_for_join, hospitalization_id, pf_dttm)

# Need numeric key for the join
analysis_with_sf[, recorded_dttm_num := as.numeric(recorded_dttm)]
setkey(analysis_with_sf, hospitalization_id, recorded_dttm_num)

analysis_with_sf_pf <- pf_for_join[analysis_with_sf,
                                     roll = MEASUREMENT_CARRY_H * 3600,
                                     on = .(hospitalization_id, pf_dttm = recorded_dttm_num)]
analysis_with_sf_pf[, pf_dttm := NULL]

analysis_with_sf <- as_tibble(analysis_with_sf_pf)

message("SF ratio available at ", sum(!is.na(analysis_with_sf$sf_ratio)),
        " / ", nrow(analysis_with_sf), " timepoints")
message("PF ratio available at ", sum(!is.na(analysis_with_sf$pf_ratio)),
        " / ", nrow(analysis_with_sf), " timepoints")

# --- Join creatinine to IMV timepoints (most recent within MEASUREMENT_CARRY_H) ---
creat_labs <- cohort_labs %>%
  filter(lab_category == "creatinine", !is.na(lab_value_numeric)) %>%
  select(hospitalization_id, creat_dttm = lab_result_dttm,
         creatinine = lab_value_numeric)

creat_for_join <- as.data.table(creat_labs)
setkey(creat_for_join, hospitalization_id, creat_dttm)

analysis_creat_dt <- as.data.table(analysis_with_sf)
analysis_creat_dt[, recorded_dttm_num := as.numeric(recorded_dttm)]
setkey(analysis_creat_dt, hospitalization_id, recorded_dttm_num)

analysis_with_creat <- creat_for_join[analysis_creat_dt,
                                       roll = MEASUREMENT_CARRY_H * 3600,
                                       on = .(hospitalization_id, creat_dttm = recorded_dttm_num)]
analysis_with_creat[, creat_dttm := NULL]

# Compute baseline creatinine (first value per patient) and delta
analysis_with_creat[, creatinine_baseline := creatinine[which(!is.na(creatinine))[1]],
                     by = hospitalization_id]
analysis_with_creat[, delta_creatinine := creatinine - creatinine_baseline]

analysis_with_sf <- as_tibble(analysis_with_creat)

message("Creatinine available at ", sum(!is.na(analysis_with_sf$creatinine)),
        " / ", nrow(analysis_with_sf), " timepoints")
message("Delta creatinine available at ", sum(!is.na(analysis_with_sf$delta_creatinine)),
        " / ", nrow(analysis_with_sf), " timepoints")

# --- Join platelet count to IMV timepoints (most recent within MEASUREMENT_CARRY_H) ---
platelet_labs <- cohort_labs %>%
  filter(lab_category == "platelet_count", !is.na(lab_value_numeric)) %>%
  select(hospitalization_id, platelet_dttm = lab_result_dttm,
         platelet_count = lab_value_numeric)

platelet_for_join <- as.data.table(platelet_labs)
setkey(platelet_for_join, hospitalization_id, platelet_dttm)

analysis_platelet_dt <- as.data.table(analysis_with_sf)
analysis_platelet_dt[, recorded_dttm_num := as.numeric(recorded_dttm)]
setkey(analysis_platelet_dt, hospitalization_id, recorded_dttm_num)

analysis_with_platelet <- platelet_for_join[analysis_platelet_dt,
                                             roll = MEASUREMENT_CARRY_H * 3600,
                                             on = .(hospitalization_id,
                                                    platelet_dttm = recorded_dttm_num)]
analysis_with_platelet[, platelet_dttm := NULL]

analysis_with_platelet[, platelet_baseline := platelet_count[which(!is.na(platelet_count))[1]],
                        by = hospitalization_id]
analysis_with_platelet[, delta_platelet := platelet_count - platelet_baseline]

analysis_with_sf <- as_tibble(analysis_with_platelet)

message("Platelet count available at ", sum(!is.na(analysis_with_sf$platelet_count)),
        " / ", nrow(analysis_with_sf), " timepoints")
message("Delta platelet available at ", sum(!is.na(analysis_with_sf$delta_platelet)),
        " / ", nrow(analysis_with_sf), " timepoints")

analysis_with_sf <- analysis_with_sf %>%
  left_join(
    cohort_demographics %>% select(hospitalization_id, admission_dttm),
    by = "hospitalization_id"
  )

# Join demographics for survival time
analysis_with_sf <- analysis_with_sf %>%
  left_join(
    cohort_demographics %>%
      select(hospitalization_id, discharge_dttm, death_dttm),
    by = "hospitalization_id"
  ) %>%
  mutate(
    # Time from admission to death (any cause), using the patient-level death_dttm,
    # which captures out-of-hospital (post-discharge) deaths as well as in-hospital
    # deaths. discharge_dttm is deliberately NOT used as the survival endpoint:
    # censoring survivors at hospital discharge would discard known post-discharge
    # vital status.
    death_day = as.numeric(difftime(death_dttm, admission_dttm, units = "days")),
    # All-cause mortality within the 60-day horizon. Patients with no recorded death
    # (or a death after day 60) are alive at the horizon and censored at day 60.
    # This assumes complete vital-status ascertainment to 60 days from the death
    # registry linkage (out-of-hospital deaths captured => no competing risk).
    mortality_event_60 = if_else(
      !is.na(death_day) & death_day >= 0 & death_day <= 60, 1L, 0L
    ),
    surv_time = if_else(mortality_event_60 == 1L, death_day, 60)
  )

# =============================================================================
# 3e. Cross-sectional cohort selection
# =============================================================================
# Cross-sectional design: one index timepoint per patient. A patient is included
# if ANY IMV timepoint with all core variables observed is lung-protective
# (VT/PBW 6-8 mL/kg) AND hypoxemic (SF ratio < 315); the index is then chosen
# among those qualifying timepoints (two-tier rule below).

# SF 315 is the Rice 2007 equivalent of PF 300 (mild ARDS / hypoxemia threshold).
SF_HYPOXEMIA_THRESHOLD <- 315

# A timepoint has "all available data" when the core analysis variables used in
# the cross-sectional models are all observed: VT/PBW, VT/PFVC and the SF ratio.
# SOFA is scored after the index is chosen (3e2), over the 24 hours from it, so it
# is not part of this gate. Driving-pressure-derived measures (dp, crs) are intentionally
# NOT required here, since plateau pressure is frequently unrecorded and would
# otherwise shrink the cohort dramatically.
analysis_with_completeness <- analysis_with_sf %>%
  mutate(
    # the control has no tidal volume: complete = SF observed
    has_all_data = (config$cohort != "imv" | (!is.na(vtpbw) & !is.na(vtpfvc))) &
      !is.na(sf_ratio)
  )

# Save the pre-gate per-timepoint dataset (all eligible IMV timepoints with derived
# variables, BEFORE the lung-protective VT/PBW gate). This lets exploratory scripts
# re-derive the cohort under a different (e.g. liberalized) VT/PBW band without
# re-running the pipeline. Not a consortium deliverable.
write_parquet(analysis_with_completeness,
              file.path(output_dir, "analysis_all_eligible_timepoints.parquet"))

# Apply the inclusion predicates BEFORE reducing to one row per patient, matching
# the original analysis (filter qualifying timepoints, then distinct(subject_id)).
# A patient is included if ANY complete-data IMV timepoint is simultaneously
# lung-protective (VT/PBW 6-8) and hypoxemic (SF < threshold); the index timepoint
# for that patient is then chosen by the two-tier rule below (preferring a
# pressure-complete timepoint early in ventilation).
# The niv cohort keeps the hypoxemia gate and has no lung-protective band to apply.
# The no-support control has no hypoxemia gate either (its patients are, by and
# large, not hypoxemic): a qualifying row is a room-air / nasal-cannula row with SF
# observed, and the index is set at ICU admission (below).
qualifying_timepoints <- analysis_with_completeness %>%
  filter(has_all_data, config$cohort != "imv" | (vtpbw >= 6 & vtpbw <= 8),
         config$cohort == "nosupport" | sf_ratio < SF_HYPOXEMIA_THRESHOLD)

# ---- the no-support control is indexed at ICU ADMISSION
# The ventilated cohort's clock starts at intubation. The control's starts at ICU
# entry: its index is the first qualifying row within CONTROL_ICU_WINDOW_H hours
# after an ICU admission, with no advanced support (HFNC, NIV, IMV) before it. A
# qualifying row earlier in the stay (the ED or the ward on arrival) is not used, so
# a week of follow-up means a week from the same kind of event in both arms. The
# escalation landmark below then runs from this index.
if (config$cohort == "nosupport") {
  icu_stays <- read_parquet(file.path(output_dir, "cohort_icu_stays.parquet"))
  first_advanced <- resp_waterfall %>%
    filter(tolower(device_category) %in% SUPPORT_DEVICES | (!is.na(tidal_volume_set) & tidal_volume_set > 0)) %>%
    group_by(hospitalization_id) %>% summarise(first_adv = min(recorded_dttm), .groups = "drop")
  icu_index <- qualifying_timepoints %>%
    inner_join(icu_stays %>% select(hospitalization_id, icu_in = in_dttm), by = "hospitalization_id",
               relationship = "many-to-many") %>%
    filter(recorded_dttm >= icu_in, recorded_dttm <= icu_in + lubridate::hours(CONTROL_ICU_WINDOW_H)) %>%
    left_join(first_advanced, by = "hospitalization_id") %>%
    filter(is.na(first_adv) | first_adv > recorded_dttm) %>%      # no advanced support before the index
    group_by(hospitalization_id) %>% slice_min(recorded_dttm, n = 1, with_ties = FALSE) %>% ungroup() %>%
    select(-icu_in, -first_adv)
  message("No-support index: ", nrow(icu_index), " patients have a qualifying row within ",
          CONTROL_ICU_WINDOW_H, " h of an ICU admission with no advanced support before it")
  qualifying_timepoints <- icu_index
}

message("Patients with >=1 complete-data timepoint: ",
        n_distinct(analysis_with_completeness$hospitalization_id[analysis_with_completeness$has_all_data]))

# Index timepoint selection (two-tier, to reduce missing pressures in the
# driving-pressure / elastance / mechanical-power analyses):
#   1. PREFERRED: the patient's first qualifying timepoint that ALSO has recorded
#      pressures (driving pressure computable, dp non-missing) within the first
#      INDEX_WINDOW_HOURS of ventilation -- so dp/crs/ers (and MP, when peak
#      pressure is also present) are observed at the index.
#   2. FALLBACK: if no qualifying timepoint within that window has pressures, use
#      the patient's first qualifying timepoint.
# The cohort (set of included patients) is unchanged -- every patient with >= 1
# qualifying timepoint is still included; only which timepoint represents them
# changes, and the chosen index always satisfies the lung-protective + hypoxemic
# eligibility criteria.

# imv_start_dttm = each patient's first IMV timepoint; the index window runs from it.
imv_start <- analysis_with_completeness %>%
  group_by(hospitalization_id) %>%
  summarise(imv_start_dttm = min(recorded_dttm), .groups = "drop")

# Tier 1: first qualifying + pressure-complete timepoint within the window.
index_tier1 <- qualifying_timepoints %>%
  left_join(imv_start, by = "hospitalization_id") %>%
  filter(!is.na(dp),
         recorded_dttm <= imv_start_dttm + lubridate::hours(INDEX_WINDOW_HOURS)) %>%
  group_by(hospitalization_id) %>%
  slice_min(recorded_dttm, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-imv_start_dttm)

# Tier 2 (fallback): first qualifying timepoint for patients not covered by tier 1.
index_tier2 <- qualifying_timepoints %>%
  filter(!hospitalization_id %in% index_tier1$hospitalization_id) %>%
  group_by(hospitalization_id) %>%
  slice_min(recorded_dttm, n = 1, with_ties = FALSE) %>%
  ungroup()

cross_sectional <- bind_rows(index_tier1, index_tier2) %>%
  # Attach weight and BMI (needed as a covariate in the driving-pressure /
  # elastance demographic-bias models in script 04). Weight is extracted and
  # imputed in script 01 but is not otherwise carried downstream.
  left_join(cohort_weights, by = "hospitalization_id") %>%
  mutate(bmi = if_else(!is.na(weight_kg) & height_cm > 0,
                       weight_kg / (height_cm / 100)^2, NA_real_))

eligible_patients <- cross_sectional$hospitalization_id
if (config$cohort != "imv") {
  # escalation: the first row of the next level of support at or after the index
  # (invasive ventilation for the niv cohort; any advanced support for the no-support
  # control), the biotrauma suite's competing event. Here t0 is the index time (not
  # the first IMV row).
  esc_devices <- if (config$cohort == "niv") "imv" else SUPPORT_DEVICES
  esc <- resp_waterfall %>%
    filter(tolower(device_category) %in% esc_devices | (!is.na(tidal_volume_set) & tidal_volume_set > 0)) %>%
    inner_join(cross_sectional %>% select(hospitalization_id, t0 = recorded_dttm), by = "hospitalization_id") %>%
    filter(recorded_dttm >= t0) %>%
    group_by(hospitalization_id) %>% summarise(escalation_dttm = min(recorded_dttm), .groups = "drop")
  cross_sectional <- cross_sectional %>% left_join(esc, by = "hospitalization_id")
  if (config$cohort == "nosupport") {
    # escalation landmark: a patient escalated within ESCALATION_LANDMARK_H of the
    # index is the pre-support stub of a supported course, not an unsupported patient
    early <- !is.na(cross_sectional$escalation_dttm) &
      cross_sectional$escalation_dttm < cross_sectional$recorded_dttm + lubridate::hours(ESCALATION_LANDMARK_H)
    message("No-support control: ", sum(early), " patients escalated within ", ESCALATION_LANDMARK_H,
            " h of the index removed")
    cross_sectional <- cross_sectional[!early, ]
    eligible_patients <- cross_sectional$hospitalization_id
  }
  message("Control cohort (", config$cohort, "): ", nrow(cross_sectional), " patients at the index; ",
          sum(!is.na(cross_sectional$escalation_dttm)), " later escalated (recorded as a competing event)")
}

if (config$cohort == "imv") {
  message("Index timepoint: ", nrow(index_tier1), " patients used a pressure-complete ",
          "qualifying timepoint within ", INDEX_WINDOW_HOURS, " h of ventilation; ",
          nrow(index_tier2), " fell back to the first qualifying timepoint.")
  message("Index timepoints with driving pressure observed: ",
          sum(!is.na(cross_sectional$dp)), " of ", nrow(cross_sectional))
  message("Included patients (VT/PBW 6-8 AND SF<", SF_HYPOXEMIA_THRESHOLD,
          " at any complete-data timepoint): ", length(eligible_patients))
}

# =============================================================================
# 3e2. SOFA: the worst values over the 24 hours from the index
# =============================================================================
# One SOFA score per patient, from the worst value of each input between the index
# timepoint and SOFA_WINDOW_H hours after it: the ventilator row chosen in 3e for
# the ventilated cohort, the first qualifying row after ICU admission for the
# no-support control. The scoring
# is clifR's compute_sofa() (github.com/AartikSarma/clifR, the commit pinned in
# uvr.toml), a port of clifpy's, so the definition is the consortium's:
#   cardiovascular  dopamine > 15, or epinephrine or norepinephrine > 0.1 mcg/kg/min: 4;
#                   dopamine > 5, or epinephrine or norepinephrine at or below 0.1: 3;
#                   dopamine at or below 5, or any dobutamine: 2; MAP < 70: 1
#   respiratory     measured PaO2 / FiO2 (the worst of each over the window); below
#                   200 scores 3-4 only on IMV, NIPPV or CPAP
#   coagulation, liver, renal, CNS   platelets, bilirubin, creatinine, GCS on the
#                   standard SOFA cut points
#   A component with no data in the window scores 0.
# Vasoactive doses are converted to mcg/kg/min first (utils/standardize_pressor_dose.R).
SOFA_COMPONENTS <- c("sofa_cv_97", "sofa_coag", "sofa_liver", "sofa_resp", "sofa_cns", "sofa_renal")
sofa_window <- cross_sectional %>%
  transmute(hospitalization_id, start_time = recorded_dttm,
            end_time = recorded_dttm + lubridate::hours(SOFA_WINDOW_H))
sofa_ids <- sofa_window$hospitalization_id
# clifR's device names (IMV, High Flow NC, ...); the waterfall holds them in lower case
clifr_device_names <- setNames(names(clifR::DEVICE_RANK_DICT), tolower(names(clifR::DEVICE_RANK_DICT)))
sofa_pressor_events <- function(drug) {
  cohort_meds %>%
    filter(med_category == drug, hospitalization_id %in% sofa_ids) %>%
    standardize_pressor_dose(weights = cohort_weights, out_col = "dose_mcg_kg_min",
                             label = paste0(drug, " (SOFA)")) %>%
    filter(!is.na(dose_mcg_kg_min)) %>%
    transmute(hospitalization_id, event_time = admin_dttm, variable = paste0(drug, "_mcg_kg_min"),
              value = dose_mcg_kg_min)
}
# every numeric input as one row per measurement: patient, time, variable, value
sofa_numeric_events <- bind_rows(
  cohort_labs %>%
    filter(hospitalization_id %in% sofa_ids, !is.na(lab_value_numeric),
           lab_category %in% c("creatinine", "bilirubin_total", "platelet_count", "po2_arterial")) %>%
    transmute(hospitalization_id, event_time = lab_result_dttm, variable = lab_category, value = lab_value_numeric),
  cohort_vitals %>%
    filter(hospitalization_id %in% sofa_ids, !is.na(vital_value), vital_category %in% c("map", "spo2")) %>%
    transmute(hospitalization_id, event_time = recorded_dttm, variable = vital_category, value = vital_value),
  cohort_assessments %>%
    filter(hospitalization_id %in% sofa_ids, assessment_category == "gcs_total") %>%
    transmute(hospitalization_id, event_time = recorded_dttm, variable = "gcs_total",
              value = as.numeric(numerical_value)) %>%
    filter(!is.na(value)),
  resp_waterfall %>%
    filter(hospitalization_id %in% sofa_ids, !is.na(fio2_set)) %>%
    transmute(hospitalization_id, event_time = recorded_dttm, variable = "fio2_set", value = fio2_set),
  map_dfr(c("norepinephrine", "epinephrine", "dopamine", "dobutamine"), sofa_pressor_events)
)
# compute_sofa() reads a wide table. One row per measurement keeps every value; its
# worst-value aggregation per patient ignores the empty cells.
sofa_wide <- bind_rows(
  sofa_numeric_events %>%
    mutate(measurement = row_number()) %>%
    pivot_wider(id_cols = c(hospitalization_id, event_time, measurement),
                names_from = variable, values_from = value) %>%
    select(-measurement),
  resp_waterfall %>%
    filter(hospitalization_id %in% sofa_ids, !is.na(device_category)) %>%
    transmute(hospitalization_id, event_time = recorded_dttm,
              device_category = unname(clifr_device_names[tolower(device_category)]))
)
for (input_column in c(clifR::MAX_ITEMS, setdiff(clifR::MIN_ITEMS, "pao2_imputed")))
  if (!input_column %in% names(sofa_wide)) sofa_wide[[input_column]] <- NA_real_
sofa_index <- clifR::compute_sofa(sofa_wide, cohort_df = sofa_window, id_name = "hospitalization_id") %>%
  select(hospitalization_id, sofa_total, all_of(SOFA_COMPONENTS))
cross_sectional <- cross_sectional %>% left_join(sofa_index, by = "hospitalization_id")
# every index row is itself a measurement inside its window, so a missing score is a bug
if (any(is.na(cross_sectional$sofa_total)))
  stop(sum(is.na(cross_sectional$sofa_total)), " patients have no SOFA over the 24 h from the index")
message("SOFA over the 24 h from the index (clifR): median ", median(cross_sectional$sofa_total),
        " (IQR ", paste(quantile(cross_sectional$sofa_total, c(0.25, 0.75)), collapse = "-"), "), ",
        nrow(cross_sectional), " patients; mean component scores: ",
        paste(sprintf("%s %.2f", sub("sofa_", "", SOFA_COMPONENTS), colMeans(cross_sectional[SOFA_COMPONENTS])),
              collapse = ", "))

# All IMV timepoints for the included patients, retained for descriptive
# summaries; no per-timepoint eligibility flags are computed here.
analysis_all <- analysis_with_sf %>%
  filter(hospitalization_id %in% eligible_patients)

# =============================================================================
# 3f. 28-day ventilator-free days (VFDs)
# =============================================================================

# Count total hours on IMV per patient within first 28 days
imv_hours_28d <- resp_waterfall %>%
  filter(hospitalization_id %in% eligible_patients,
         tolower(device_category) == "imv") %>%
  inner_join(
    cohort_demographics %>% select(hospitalization_id, admission_dttm),
    by = "hospitalization_id"
  ) %>%
  mutate(t_days = as.numeric(difftime(recorded_dttm, admission_dttm, units = "days"))) %>%
  filter(t_days >= 0, t_days <= 28) %>%
  group_by(hospitalization_id) %>%
  summarise(imv_hours = n(), .groups = "drop") %>%
  mutate(imv_days = imv_hours / 24)

# VFD-28: 28 minus ventilator days, deaths within 28 days get 0.
# Death is all-cause (in- or out-of-hospital) within 28 days, from the patient-level
# death_dttm — consistent with the 60-day survival endpoint. An in-hospital-only
# condition (deceased == 1) would credit ventilator-free days to patients who died out
# of hospital within 28 days.
vfd_data <- cohort_demographics %>%
  filter(hospitalization_id %in% eligible_patients) %>%
  select(hospitalization_id, deceased, death_dttm, admission_dttm) %>%
  left_join(imv_hours_28d, by = "hospitalization_id") %>%
  mutate(
    imv_days = replace_na(imv_days, 0),
    death_day = as.numeric(difftime(death_dttm, admission_dttm, units = "days")),
    died_within_28 = !is.na(death_day) & death_day >= 0 & death_day <= 28,
    vfd_28 = if_else(died_within_28, 0, pmax(28 - imv_days, 0)),
    # --- Competing-risks VFD outcome (Yehya & Harhay, AJRCCM 2019) -------------
    # VFDs are best analyzed as a competing-risks outcome rather than a continuous
    # value: event of interest = extubation (liberation from ventilation),
    # competing risk = death within 28 days, censored at day 28 if still
    # ventilated. Time to liberation is the cumulative ventilator days within 28 d
    # (the duration component of the VFD construct); any death within 28 d is the
    # competing event (consistent with the VFD convention that death -> 0 VFDs).
    # vfd_status: 0 = censored (still ventilated at 28 d), 1 = extubation,
    # 2 = death. vfd_time is floored at ~1 h so the survival models have time > 0.
    vent_days_28 = pmin(imv_days, 28),
    vfd_status = case_when(
      died_within_28     ~ 2L,
      vent_days_28 >= 28 ~ 0L,
      TRUE               ~ 1L
    ),
    vfd_time = pmax(
      case_when(
        vfd_status == 2L ~ death_day,
        vfd_status == 0L ~ 28,
        TRUE             ~ vent_days_28
      ),
      1 / 24
    )
  ) %>%
  select(hospitalization_id, vfd_28, vfd_time, vfd_status)

# Attach VFDs to the cross-sectional cohort (one row per included patient)
cross_sectional <- cross_sectional %>%
  left_join(vfd_data, by = "hospitalization_id")

message("28-day VFDs computed. Median VFD-28: ",
        round(median(vfd_data$vfd_28, na.rm = TRUE), 1))

# =============================================================================
# 3h. Save outputs
# =============================================================================

write_parquet(analysis_all, file.path(output_dir, "analysis_all_timepoints.parquet"))
write_parquet(cross_sectional, file.path(output_dir, "analysis_cross_sectional.parquet"))

# =============================================================================
# 3i. Complete the attrition log (steps 4-7) and write the full CONSORT table
# =============================================================================
# Patient-level, monotonic counts continuing the funnel from script 01.
n_step4 <- n_distinct(pbw_pfvc_data$hospitalization_id)
n_step5 <- n_distinct(
  analysis_with_completeness$hospitalization_id[analysis_with_completeness$has_all_data]
)
n_step6 <- analysis_with_completeness %>%
  filter(has_all_data, config$cohort != "imv" | (vtpbw >= 6 & vtpbw <= 8)) %>%
  summarise(n = n_distinct(hospitalization_id)) %>% pull(n)
n_step7 <- length(eligible_patients)

partial_path <- file.path(output_dir, "attrition_log_partial.csv")
if (!file.exists(partial_path)) {
  stop("Missing attrition_log_partial.csv from script 01: ", partial_path)
}
# Step 4 also loses patients with missing or unknown sex (PBW is undefined) and
# those whose PFVC cannot be computed, so its reason names them.
attrition_steps <- attrition_steps_for(config$cohort)   # step labels (utils/attrition_log.R)
attrition <- read_csv(partial_path, show_col_types = FALSE) %>%
  attrition_add(attrition_steps[4], n_step4,
                exclusion_reason = "Height missing or outside 150-210 cm, or sex, PBW or PFVC missing") %>%
  attrition_add(attrition_steps[5], n_step5,
                exclusion_reason = if (config$cohort != "imv") "Incomplete index data (SF)" else
                  "Incomplete index data (VT/PBW, VT/PFVC, SF)") %>%
  attrition_add(attrition_steps[6], n_step6,
                exclusion_reason = if (config$cohort != "imv") "(no tidal-volume band applies)" else
                  "Not lung-protective (VT/PBW outside 6-8)") %>%
  attrition_add(attrition_steps[7], n_step7,
                exclusion_reason = if (config$cohort == "nosupport")
                  paste0("No qualifying row within ", CONTROL_ICU_WINDOW_H, " h of ICU admission, or escalated within ",
                         ESCALATION_LANDMARK_H, " h of the index") else "Not hypoxemic (SF ratio >= 315)") %>%
  mutate(site = site_name, .before = 1)

write_csv(attrition, file.path(final_dir, paste0("attrition_log_", site_name, ".csv")))
message("Attrition log written (7 steps): ",
        paste(attrition$n_remaining, collapse = " -> "))

# =============================================================================
# 3k. Negative-control cohorts: PBW/PFVC + survival fields (from script 01)
# =============================================================================
# Non-hypoxemic patients, ventilated and not, in which no lung-protective dosing
# decision was made (see script 01). The same Devine / GLI-2012 derivation and the
# same height window as the analytic cohort; the same 60-day all-cause survival
# fields. Script 04 (4j) fits PBW/PFVC, PFVC and height against mortality in each.
# Ventilated control: non-hypoxemic over the whole stay AND classifiably non-hypoxemic
# over the entire ventilated period (script 01: SF >= 315 at every SpO2, PaO2/FiO2 >= 300,
# FiO2 never > 0.40, FiO2 documented). Non-ventilated control: non-hypoxemic over the stay.
nc_cohort <- read_parquet(file.path(output_dir, "nc_cohort.parquet")) %>%
  filter(!hypoxemic, !is.na(height_cm), height_cm >= 150, height_cm <= 210,
         !is.na(age_at_admission), sex_category %in% c("Male", "Female"),
         !imv_set_vt | (!is.na(hypoxemic_during_vent) & !hypoxemic_during_vent)) %>%
  mutate(
    # Labels say what each cohort controls for: the ventilated non-hypoxemic group
    # still receives PBW-dosed tidal volumes (a positive control for the dosing
    # pathway, a negative control for hypoxemia); the non-ventilated group is the
    # true no-dose control.
    nc_cohort = if_else(imv_set_vt, "Ventilated, non-hypoxemic (dosed, uninjured lung)",
                        "Not ventilated (no tidal volume)"),
    sex_numeric  = if_else(sex_category == "Male", 1L, 2L),
    race_numeric = case_when(race_category == "WHITE" ~ 1L, race_category == "BLACK" ~ 2L, TRUE ~ 5L),
    pbw = if_else(sex_numeric == 1L, 50.0 + 2.3 * (height_cm / 2.54 - 60), 45.5 + 2.3 * (height_cm / 2.54 - 60)),
    pfvc = pred_GLI(age = age_at_admission, height = height_cm / 100, gender = sex_numeric,
                    ethnicity = race_numeric, param = "FVC"),
    pfvc_age25 = pred_GLI(age = rep(25, n()), height = height_cm / 100, gender = sex_numeric,
                          ethnicity = race_numeric, param = "FVC"),
    pbwpfvc = pbw / pfvc,
    vt_excess_ml = VT_PER_KG_ARMA * pbw - (VT_PCT_PFVC_ARMA / 100) * pfvc * 1000,
    sex_category  = factor(sex_category,  levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
    death_day = as.numeric(difftime(death_dttm, admission_dttm, units = "days")),
    mortality_event_60 = if_else(!is.na(death_day) & death_day >= 0 & death_day <= 60, 1L, 0L),
    surv_time = if_else(mortality_event_60 == 1L, death_day, 60)
  ) %>%
  filter(!is.na(pfvc), pfvc > 0)
# Delivered dose for the ventilated control: these patients are in the script-01 IMV
# cohort (the hypoxemia gate is applied here, in 3e), so their volume-targeted
# timepoints are in the pre-gate frame. Per-hospitalization median VT/PBW and VT/PFVC
# over all IMV timepoints with a set tidal volume.
nc_dose <- analysis_with_completeness %>%
  filter(hospitalization_id %in% nc_cohort$hospitalization_id, !is.na(vtpbw), !is.na(vtpfvc)) %>%
  group_by(hospitalization_id) %>%
  summarise(vtpbw = median(vtpbw), vtpfvc = median(vtpfvc), n_vt_timepoints = n(), .groups = "drop")
nc_cohort <- nc_cohort %>% left_join(nc_dose, by = "hospitalization_id")
write_parquet(nc_cohort, file.path(output_dir, "analysis_negative_control.parquet"))
message("Negative-control cohorts: ", paste(capture.output(print(table(nc_cohort$nc_cohort))), collapse = " | "))

# =============================================================================
# 3l. Federated PBW:PFVC distribution exports (site-specific; poolable)
# =============================================================================
# Aggregated summaries of the PBW:PFVC ratio by demographic group only -- no
# row-level data leaves the site. (a) histograms on fixed bins (summable across
# sites into a pooled distribution); (b) quantile summaries for exact per-site
# boxplots.
dist_groups <- bind_rows(
  cross_sectional %>% transmute(group_type = "sex",
    group_value = as.character(sex_category), value = pbwpfvc),
  cross_sectional %>% transmute(group_type = "race",
    group_value = as.character(race_category), value = pbwpfvc),
  cross_sectional %>%
    mutate(gv = cut(age_at_admission, AGE_BIN_EDGES, right = FALSE)) %>%
    transmute(group_type = "age_bin", group_value = as.character(gv), value = pbwpfvc),
  cross_sectional %>%
    mutate(gv = cut(height_cm, HEIGHT_BIN_EDGES, right = FALSE)) %>%
    transmute(group_type = "height_bin", group_value = as.character(gv), value = pbwpfvc)
) %>%
  filter(!is.na(value), !is.na(group_value), group_value != "NA")


# (a) histograms with tails clamped into the edge bins so per-group totals == N
.bin_lo <- min(PBWPFVC_BIN_EDGES); .bin_hi <- max(PBWPFVC_BIN_EDGES)
dist_histograms <- dist_groups %>%
  mutate(value = pmin(pmax(value, .bin_lo), .bin_hi - 1e-9),
         bin_i = findInterval(value, PBWPFVC_BIN_EDGES, rightmost.closed = TRUE)) %>%
  count(group_type, group_value, bin_i, name = "count") %>%
  mutate(bin_left  = PBWPFVC_BIN_EDGES[bin_i],
         bin_right = PBWPFVC_BIN_EDGES[bin_i + 1]) %>%
  transmute(site = site_name, group_type, group_value, bin_left, bin_right, count) %>%
  arrange(group_type, group_value, bin_left)

# (b) per-group quantile summaries
dist_quantiles <- dist_groups %>%
  group_by(group_type, group_value) %>%
  summarise(
    n = n(), median = median(value),
    q1 = quantile(value, 0.25), q3 = quantile(value, 0.75),
    p10 = quantile(value, 0.10), p20 = quantile(value, 0.20),
    p30 = quantile(value, 0.30), p40 = quantile(value, 0.40),
    p50 = quantile(value, 0.50), p60 = quantile(value, 0.60),
    p70 = quantile(value, 0.70), p80 = quantile(value, 0.80),
    p90 = quantile(value, 0.90),
    .groups = "drop"
  ) %>%
  mutate(site = site_name, .before = 1)

write_csv(dist_histograms, file.path(final_dir, paste0("dist_histograms_", site_name, ".csv")))
write_csv(dist_quantiles, file.path(final_dir, paste0("dist_quantiles_", site_name, ".csv")))
message("Federated distribution exports written for ",
        n_distinct(dist_groups$group_value), " demographic groups")

message("Script 03 complete.")
message("  Cross-sectional (1 row/included patient): ", nrow(cross_sectional))
message("  All timepoints (included patients): ", nrow(analysis_all))
