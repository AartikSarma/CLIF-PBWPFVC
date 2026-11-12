# =============================================================================
# causal_common: the data and the within-person and ladder designs. Sourced by
# causal_triangulation.R, causal_testbed.R and causal_tests.R; never run; writes
# nothing. Also here: the outcomes, the time course and the decision rule. The planted-truth generator that tests the designs is causal_generator.R.
# =============================================================================
#
# The hypothesis: at a fixed VT/PBW, a high VT/PFVC causes biotrauma. At a fixed
# VT/PBW, log VT/PFVC moves only with lambda = log(PBW/PFVC), and lambda is a fixed
# function of age, sex, race and height. So a model that adjusts for age, sex and race
# leaves lambda only its height variation within sex (a few percent of its variance):
# the adjusted contrast is close to a positivity violation, and the unadjusted one
# cannot tell strain from anything else those demographics do to a sick patient.
#
# No single design answers that. Each design here gets the lambda gradient of a
# marker's slope (change in log marker per day, per log unit of PBW/PFVC) from a
# different source of variation, and so fails under a different assumption:
#
#   unadjusted        ventilated, on the ventilator; severity only.
#                     Fails if demographics move the marker for any other reason.
#   adjusted          the same with ns(age, 4), sex and race: height within sex is all
#                     that identifies it. Valid under any demographic confounding that
#                     those terms hold, and nearly uninformative (see partial_r2).
#   fingerprint       28_height_fingerprint.R's contrast: adjusted, with a height smooth
#                     the sexes share; only the sex-reversed height curve is left.
#   control_did       ventilated minus no support (figure 4's design): the same
#                     demographics, no tidal volume. Fails if demographics move the
#                     marker differently in the two cohorts for a reason other than
#                     strain (sicker patients, or ventilation itself).
#   within_onoff      NEW. Each ventilated patient before and after intubation:
#                     the lambda gradient on the ventilator minus the same patients'
#                     gradient in the days before it. A patient's demographics do not
#                     change at intubation, so nothing is adjusted and positivity is
#                     not at stake. Fails if the demographic gradient of an illness
#                     steepens when the illness does (intubation marks deterioration).
#   within_triple     NEW. within_onoff minus the same before/after change in the
#                     no-support controls around their ICU admission, which removes a
#                     demographic gradient that steepens with any ICU-level illness.
#                     Fails only if demographics act through ventilation itself by a
#                     path other than strain.
#   off_reversal      NEW. The gradient on the ventilator minus the gradient in the
#                     days after extubation, in the patients liberated alive. Strain
#                     stops at extubation; age does not. Selected on survival to
#                     extubation, so it is read as support, not as the estimate.
#   free_<input>      NEW. The ladder: the unadjusted model with ONE demographic input
#                     set free (age as ns(age, 4), sex, race, or a height smooth), so
#                     lambda is identified by the other three. If one input's own
#                     effect is the confounder, the rung that frees it drops to zero
#                     and the others do not; under strain all four agree. Each rung
#                     keeps most of lambda's variance, unlike the adjusted model.
#   free_<input>_did  the same rung as a ventilated-minus-control contrast.
#
# The estimate is always the same quantity, so the designs can be laid side by side
# and each can be scored against a planted truth (causal_generator.R, causal_testbed.R).
#
# The model. Every design is one regression of log marker on time with a separate
# intercept for every patient-phase (so levels, and everything fixed about a patient,
# drop out), and slopes that depend on patient covariates:
#     log y_ij = a_(i,phase) + t_ij * (x_(i,phase)' gamma) + e_ij
# fitted by least squares after centring within patient-phase, with standard errors
# clustered on the patient (which allows each patient a slope of their own). It is the
# longitudinal model alone: like 28_height_fingerprint.R it is conditional on staying
# under observation, with no joint model for death or liberation.
#
# Phases, on a clock in days from t0:
#   pre    [-PRE_DAYS, 0)                    before intubation (control: before the index)
#   on     [0, min(t_off, ON_DAYS))          on the ventilator (control: after the index,
#                                            to escalation or ON_DAYS)
#   post   [t_off, t_off + POST_DAYS)        after liberation (ventilated only)
# t0 is the first invasive-ventilation record of the stay (the on switch; script 03's
# index is the first qualifying timepoint, which can be hours later), and the control's
# index. t_off is the last invasive-ventilation record.
#
# Needs: utils/config.R and code/20_biotrauma_grid.R sourced (pfvc_channels(),
# height_fingerprint()), and scripts 01-03 run for the ventilated cohort and for
# PBWPFVC_COHORT=nosupport.
# =============================================================================

PRE_DAYS  <- as.numeric(Sys.getenv("PBWPFVC_CAUSAL_PRE_DAYS", "3"))
ON_DAYS   <- as.numeric(Sys.getenv("PBWPFVC_CAUSAL_ON_DAYS", "7"))
POST_DAYS <- as.numeric(Sys.getenv("PBWPFVC_CAUSAL_POST_DAYS", "3"))
MIN_SPAN_D <- 0.25          # a patient-phase needs two measurements at least 6 h apart
MIN_PATIENTS <- 50L         # a design with fewer patients is reported as not estimable
# the ventilated arm's VT/PBW slope term (PBWPFVC_DID_VTPBW, as in the other
# ventilated-versus-control supplements): "1" holds the index VT/PBW fixed, "0" drops it
CAUSAL_VTPBW <- identical(Sys.getenv("PBWPFVC_DID_VTPBW", "1"), "1")

# The outcomes. source says where the measurements come from, family how the slope
# model reads them:
#   gaussian   least squares on y (the log of a lab or of OSI; a 0/1 indicator, which
#              makes the model a linear probability model)
#   poisson    a multiplicative model of the MEAN of a non-negative outcome with its
#              zeros kept: log E[y] = a_(i,phase) + t * (x' gamma), fitted by Poisson
#              quasi-likelihood with the patient-phase intercepts conditioned out. It
#              needs only that mean to be right (not Poisson variation), has no added
#              constant, and its gamma is in the same units as a log marker's.
# lab markers come from cohort_labs_clean (script 02). The vasopressor dose is the
# norepinephrine-equivalent rate of script 02's ne_equiv_admin, as the highest rate
# charted in each PRESSOR_BIN_H-hour bin a patient was under observation, zero when
# none was charted. OSI is 21_biotrauma_panel.R's definition (FiO2 x mean airway
# pressure / SpO2, SpO2 clamped to 80-97), as the median of each bin; it exists only
# on the ventilator, so only the ventilated on-phase designs can be fitted to it.
CAUSAL_MARKERS <- tibble::tribble(
  ~marker,        ~source,   ~lab_category,    ~injury_sign, ~family,
  "platelets",    "lab",     "platelet_count", -1,           "gaussian",
  "creatinine",   "lab",     "creatinine",     +1,           "gaussian",
  "pressor_dose", "pressor", NA,               +1,           "poisson",
  "pressor_any",  "pressor", NA,               +1,           "gaussian",
  "osi",          "osi",     NA,               +1,           "gaussian")
PRESSOR_BIN_H <- as.numeric(Sys.getenv("PBWPFVC_CAUSAL_BIN_H", "6"))
# the on-ventilator phase cut into segments, for the time course (days from t0)
CAUSAL_KNOTS <- as.numeric(strsplit(Sys.getenv("PBWPFVC_CAUSAL_KNOTS", "2,4"), ",")[[1]])

# The strain denominator. lambda = log(PBW / predicted volume); PFVC is the primary one.
# The others are script 03's reference volumes, missing above their age limit (80 for
# GLI 2021), and have no race term.
CAUSAL_DENOMINATORS <- c(pfvc = "ldisc", pfvc_age25 = "ldisc25", tlc = "ldisc_tlc",
                         alveolar_volume = "ldisc_va", frc = "ldisc_frc")

# ---- Data -----------------------------------------------------------------------
# One row per patient, both arms. Columns: id, arm, vent (1/0), t0, t_off_d (days from
# t0; NA when never off), the demographics, lambda (ldisc) and its pieces, severity.
load_causal_patients <- function(config) {
  site_root <- dirname(config$output_dir)
  if (config$cohort != "imv") stop("the causal designs read both cohorts themselves; unset PBWPFVC_COHORT")
  vent_file <- file.path(config$output_dir, paste0("analysis_cross_sectional", config$cs_suffix, ".parquet"))
  ctrl_dir  <- file.path(config$output_dir, "controls", "nosupport")
  ctrl_file <- file.path(ctrl_dir, "analysis_cross_sectional.parquet")
  if (!file.exists(vent_file)) stop("missing ", vent_file, ": run scripts 01-03 first")
  if (!file.exists(ctrl_file)) stop("no no-support cohort: run scripts 01-03 with PBWPFVC_COHORT=nosupport first")
  keep <- c("hospitalization_id", "patient_id", "index_dttm", "age_at_admission", "sex_category",
            "race_category", "height_cm", "pbw", "pfvc", "pfvc_age25", "pva", "ptlc", "pfrc", "sf_ratio", "sofa_total",
            "admission_dttm", "discharge_dttm", "death_dttm")
  vent <- arrow::read_parquet(vent_file) %>%
    select(all_of(keep), vtpbw, last_imv_dttm, any_of(c("crs", "peep_set", "bmi"))) %>% mutate(arm = "ventilated")
  for (col in c("crs", "peep_set", "bmi")) if (!col %in% names(vent)) vent[[col]] <- NA_real_
  ctrl <- arrow::read_parquet(ctrl_file) %>%
    select(all_of(keep), escalation_dttm, any_of("bmi")) %>% mutate(arm = "no support", vtpbw = NA_real_, crs = NA_real_, peep_set = NA_real_)
  if (!"bmi" %in% names(ctrl)) ctrl$bmi <- NA_real_
  # the on switch: the stay's first invasive-ventilation record
  wf_file <- file.path(config$output_dir, "resp_support_waterfall_clean.parquet")
  if (!file.exists(wf_file)) stop("missing ", wf_file, ": run script 02 first")
  first_imv <- arrow::read_parquet(wf_file, col_select = c("hospitalization_id", "recorded_dttm", "device_category")) %>%
    filter(tolower(device_category) == "imv") %>%
    group_by(hospitalization_id) %>% summarise(first_imv_dttm = min(recorded_dttm), .groups = "drop")
  vent <- vent %>% left_join(first_imv, by = "hospitalization_id") %>%
    mutate(t0 = pmin(coalesce(first_imv_dttm, index_dttm), index_dttm),
           t_off = last_imv_dttm) %>% select(-first_imv_dttm, -last_imv_dttm)
  ctrl <- ctrl %>% filter(!patient_id %in% vent$patient_id) %>%
    mutate(t0 = index_dttm, t_off = escalation_dttm) %>% select(-escalation_dttm)
  days <- function(a, b) as.numeric(difftime(a, b, units = "days"))
  bind_rows(vent, ctrl) %>%
    filter(sex_category %in% c("Male", "Female"), !is.na(height_cm), !is.na(pfvc), !is.na(pbw),
           !is.na(sf_ratio), !is.na(sofa_total), is.na(vtpbw) | vtpbw > 0) %>%
    mutate(id = hospitalization_id, vent = as.numeric(arm == "ventilated"),
           t_off_d = days(t_off, t0), stay_end_d = days(discharge_dttm, t0),
           # hospital time before t0 (what a pre phase can be drawn from) and the day of death
           pre_avail_d = pmin(pmax(days(t0, admission_dttm), 0), PRE_DAYS),
           death_d = days(death_dttm, t0),
           index_lag_h = 24 * days(index_dttm, t0),     # hours from the start of ventilation to the index
           age10 = age_at_admission / 10,
           sex_category = factor(as.character(sex_category), levels = c("Male", "Female")),
           race_category = factor(as.character(race_category), levels = c("WHITE", "BLACK", "OTHER")),
           ldisc = log(pbw / pfvc),
           lpbw = log(pbw),             # log PBW alone, for the designs with PBW and PFVC as separate terms
           # aeration: index compliance over predicted lung size (the fraction still
           # aerated, up to a constant); 0 in the control, which has no compliance.
           # NA where a ventilated patient has no compliance at the index.
           laer0 = if_else(vent == 1, if_else(!is.na(crs) & crs > 0, log(crs / pfvc), NA_real_), 0),
           laer0 = laer0 - mean(laer0[vent == 1], na.rm = TRUE) * vent,
           lbmi = if_else(!is.na(bmi) & bmi > 0, log(bmi), NA_real_), lbmi = lbmi - mean(lbmi[vent == 1], na.rm = TRUE),
           fp = height_fingerprint(height_cm, as.character(sex_category)),
           # log VT/PBW about the ventilated mean; 0 in the control, which has no tidal volume
           lv0 = if_else(vent == 1, log(vtpbw) - mean(log(vtpbw[vent == 1]), na.rm = TRUE), 0),
           sofa_z = as.numeric(scale(sofa_total)), sf_z = as.numeric(scale(sf_ratio)),
           age_z = as.numeric(scale(age10))) %>%
    mutate(# the same ratio on two other denominators, for the testbed's wrong-denominator
           # scenarios: GLI at age 25 (no age term) and GLI without its race term
           ldisc25 = log(pbw / pfvc_age25),
           ldisc_norace = ldisc + pfvc_channels(pick(height_cm, age10, sex_category, race_category), "log_pfvc")$ch_race,
           # lambda's age part (minus GLI's age term in log PFVC, other inputs at reference)
           # and the rest of lambda, for the age-mediation check
           lage = -pfvc_channels(pick(height_cm, age10, sex_category, race_category), "log_pfvc")$ch_age,
           lrest = ldisc - lage,
           # the other reference volumes (NA above their age limit)
           ldisc_tlc = log(pbw / ptlc), ldisc_va = log(pbw / pva), ldisc_frc = log(pbw / pfrc)) %>%
    mutate(across(c(ldisc, ldisc25, ldisc_norace, ldisc_tlc, ldisc_va, ldisc_frc, lpbw, lage, lrest),
                  ~ .x - mean(.x[vent == 1], na.rm = TRUE))) %>%
    select(id, patient_id, arm, vent, t0, t_off_d, stay_end_d, pre_avail_d, death_d, age10, age_z,
           sex_category, race_category, height_cm, pbw, pfvc, ptlc, pva, pfrc, pfvc_age25, crs, peep_set, index_lag_h, sf_ratio_raw = sf_ratio, laer0, lbmi, lpbw, lage, lrest, ldisc, ldisc25, ldisc_norace, ldisc_tlc, ldisc_va, ldisc_frc,
           fp, lv0, sofa_z, sf_z)
}

# The patients with lambda on another denominator. `common` keeps only the patients in
# whom every denominator is defined, so that denominators are compared on one sample.
# ldiff is PFVC's lambda minus this one: what PFVC adds to the denominator (mostly
# GLI's age decline in FVC, and its race term).
use_denominator <- function(patients, denominator, common = TRUE) {
  col <- CAUSAL_DENOMINATORS[[denominator]]
  ok <- if (common) stats::complete.cases(patients[unname(CAUSAL_DENOMINATORS)]) else !is.na(patients[[col]])
  out <- patients[ok, , drop = FALSE]
  pfvc <- out$ldisc
  out$ldisc <- out[[col]] - mean(out[[col]][out$vent == 1])
  out$ldiff <- (pfvc - mean(pfvc[out$vent == 1])) - out$ldisc
  out
}

# the phase of a time t (days from t0)
phase_of <- function(t, vent, t_off_d) {
  on_end <- pmin(dplyr::coalesce(t_off_d, Inf), ON_DAYS)
  dplyr::case_when(t >= -PRE_DAYS & t < 0 ~ "pre",
                   t >= 0 & t < on_end ~ "on",
                   vent == 1 & !is.na(t_off_d) & t_off_d < ON_DAYS + POST_DAYS &
                     t >= t_off_d & t < t_off_d + POST_DAYS ~ "post",
                   TRUE ~ NA_character_)
}
both_cohorts <- function(config, file) {
  dirs <- c(config$output_dir, file.path(config$output_dir, "controls", "nosupport"))
  lapply(dirs, function(dir) {
    f <- file.path(dir, file)
    if (!file.exists(f)) stop("missing ", f, ": run script 02 first")
    f
  })
}

# Every measurement of one lab in the window, with its phase. t is days from t0.
# complete = FALSE keeps the patient-phases with a single measurement too (the testbed
# calibrates its measurement rate on them).
load_causal_panel <- function(config, patients, lab_category, complete = TRUE) {
  read_labs <- function(f) {
    labs <- arrow::read_parquet(f)
    value <- intersect(c("lab_value_numeric", "lab_value"), names(labs))[1]
    labs %>% filter(lab_category == !!lab_category) %>%
      transmute(id = hospitalization_id, lab_dttm = lab_result_dttm, value = suppressWarnings(as.numeric(.data[[value]])))
  }
  purrr::map_dfr(both_cohorts(config, "cohort_labs_clean.parquet"), read_labs) %>%
    inner_join(patients %>% select(id, vent, t0, t_off_d), by = "id") %>%
    mutate(t = as.numeric(difftime(lab_dttm, t0, units = "days")), phase = phase_of(t, vent, t_off_d)) %>%
    filter(!is.na(phase), is.finite(value), value > 0) %>%
    transmute(id, phase, t, y = log(value)) %>%
    distinct() %>%
    complete_phases(complete) %>%
    arrange(id, phase, t)
}

# The vasopressor dose on a regular grid: one row per patient and bin the patient was
# under observation for (in hospital, alive, inside one phase), y the highest
# norepinephrine-equivalent rate charted in the bin and 0 when none was. A continuous
# infusion is charted when it changes, so a bin with no record is read as no infusion,
# as scripts 10 and 21 read a day with none.
load_pressor_panel <- function(config, patients, complete = TRUE) {
  bin <- PRESSOR_BIN_H / 24
  dose <- purrr::map_dfr(both_cohorts(config, "ne_equiv_admin.parquet"), function(f)
    arrow::read_parquet(f) %>% transmute(id = hospitalization_id, dttm = admin_dttm, dose = ne_equiv_total)) %>%
    inner_join(patients %>% select(id, t0), by = "id") %>%
    mutate(k = floor(as.numeric(difftime(dttm, t0, units = "days")) / bin)) %>%
    filter(is.finite(dose), dose >= 0) %>%
    group_by(id, k) %>% summarise(dose = max(dose), .groups = "drop")
  tidyr::expand_grid(id = patients$id, k = seq(-PRE_DAYS / bin, (ON_DAYS + POST_DAYS) / bin - 1)) %>%
    inner_join(patients %>% select(id, vent, t_off_d, stay_end_d, death_d, pre_avail_d), by = "id") %>%
    mutate(start = k * bin, t = start + bin / 2,
           observed_to = pmin(dplyr::coalesce(stay_end_d, Inf), dplyr::coalesce(death_d, Inf)),
           phase = phase_of(start, vent, t_off_d)) %>%
    # the whole bin under observation and inside one phase
    filter(!is.na(phase), start >= -pre_avail_d, start + bin <= observed_to,
           phase == phase_of(start + bin - 1e-9, vent, t_off_d)) %>%
    left_join(dose, by = c("id", "k")) %>%
    transmute(id, phase, t, y = dplyr::coalesce(dose, 0)) %>%
    complete_phases(complete) %>%
    arrange(id, phase, t)
}

# The oxygen saturation index, log scale, as the median of each bin on the ventilator.
# 21_biotrauma_panel.R's definition: at each SpO2 (clamped to 80-97), the last FiO2
# within 4 hours and the nearest recorded mean airway pressure within 1 hour.
load_osi_panel <- function(config, patients, complete = TRUE) {
  bin <- PRESSOR_BIN_H / 24
  vent <- patients %>% filter(vent == 1) %>% select(id, t0, t_off_d)
  wf <- arrow::read_parquet(file.path(config$output_dir, "resp_support_waterfall_clean.parquet"),
                            col_select = c("hospitalization_id", "recorded_dttm", "fio2_set", "mean_airway_pressure_obs")) %>%
    rename(id = hospitalization_id) %>% filter(id %in% vent$id) %>% mutate(s = as.numeric(recorded_dttm))
  fio2 <- data.table::as.data.table(wf %>% filter(!is.na(fio2_set), fio2_set > 0) %>% select(id, s, fio2_set))
  maw  <- data.table::as.data.table(wf %>% filter(!is.na(mean_airway_pressure_obs), mean_airway_pressure_obs > 0) %>%
                                      transmute(id, s, maw = mean_airway_pressure_obs, s_maw = s))
  spo2 <- arrow::read_parquet(file.path(config$output_dir, "cohort_vitals_clean.parquet")) %>%
    filter(vital_category == "spo2", hospitalization_id %in% vent$id, !is.na(vital_value), vital_value >= 50, vital_value <= 100) %>%
    transmute(id = hospitalization_id, s = as.numeric(recorded_dttm), spo2 = pmin(pmax(vital_value, 80), 97))
  spo2 <- data.table::as.data.table(spo2)
  data.table::setkeyv(fio2, c("id", "s")); data.table::setkeyv(maw, c("id", "s")); data.table::setkeyv(spo2, c("id", "s"))
  pts <- fio2[spo2, roll = 4 * 3600, on = c("id", "s")]
  pts <- maw[pts, roll = "nearest", on = c("id", "s")]
  tibble::as_tibble(pts) %>%
    filter(!is.na(fio2_set), !is.na(maw), abs(s - s_maw) <= 3600) %>%
    mutate(fio2 = if_else(fio2_set > 1.5, fio2_set / 100, fio2_set), osi = 100 * maw * fio2 / spo2) %>%
    inner_join(vent, by = "id") %>%
    mutate(t = (s - as.numeric(t0)) / 86400, phase = phase_of(t, 1, t_off_d)) %>%
    filter(phase == "on", is.finite(osi), osi > 0) %>%
    group_by(id, phase, k = floor(t / bin)) %>%
    summarise(t = (k[1] + 0.5) * bin, y = log(median(osi)), .groups = "drop") %>%
    select(id, phase, t, y) %>%
    complete_phases(complete) %>%
    arrange(id, phase, t)
}

# one outcome's panel, by its name in CAUSAL_MARKERS
load_causal_outcome <- function(config, patients, marker, complete = TRUE) {
  m <- CAUSAL_MARKERS[CAUSAL_MARKERS$marker == marker, ]
  if (nrow(m) != 1) stop("unknown marker ", marker)
  panel <- switch(m$source,
                  lab = load_causal_panel(config, patients, m$lab_category, complete),
                  pressor = load_pressor_panel(config, patients, complete),
                  osi = load_osi_panel(config, patients, complete))
  if (marker == "pressor_any") panel$y <- as.numeric(panel$y > 0)
  panel
}
# a patient-phase enters a design with two measurements at least MIN_SPAN_D apart
complete_phases <- function(panel, complete = TRUE) {
  if (!complete || !nrow(panel)) return(panel)
  key <- paste(panel$id, panel$phase)
  ok <- tapply(panel$t, key, function(t) length(t) >= 2 && max(t) - min(t) >= MIN_SPAN_D)
  panel[ok[key], , drop = FALSE]
}

# ---- The designs ----------------------------------------------------------------
# A design is: the arms and phases it reads, whether a patient must have every phase,
# the slope formula over patient-phase covariates, and the coefficient that is the
# lambda gradient. `on` is 1 on the ventilator (after the index, in the control) and 0
# in the comparison phase.
causal_designs <- function(vtpbw = CAUSAL_VTPBW) {
  dose <- if (vtpbw) " + lv_on" else ""          # index VT/PBW acts only while ventilated
  sev  <- "sofa_z + sf_z"
  demo <- "splines::ns(age10, 4) + sex_category + race_category"
  free <- c(age = "splines::ns(age10, 4)", sex = "sex_category", race = "race_category",
            height = "splines::ns(height_cm, 3)")
  d <- function(design, arms, phases, formula, target, sign = 1, fails_if)
    tibble(design = design, arms = list(arms), phases = list(phases), formula = formula,
           target = target, sign = sign, fails_if = fails_if)
  V <- "ventilated"; C <- "no support"
  bind_rows(
    d("unadjusted", V, "on", paste0("~ ldisc + ", sev, dose), "ldisc",
      fails_if = "demographics move the marker by any path other than strain"),
    d("adjusted", V, "on", paste0("~ ldisc + ", sev, " + ", demo, dose), "ldisc",
      fails_if = "height acts on the marker directly; identified by height within sex only (positivity)"),
    d("fingerprint", V, "on", paste0("~ fp + splines::ns(height_cm, 3) + ", sev, " + ", demo, dose), "fp",
      fails_if = "a height effect that differs by sex; identified by the sex-reversed height curve only"),
    d("control_did", c(V, C), "on", paste0("~ vent * (", sev, ") + ldisc + vent:ldisc", dose), "vent:ldisc",
      fails_if = "demographics move the marker differently in ventilated and control patients, apart from strain"),
    d("within_onoff", V, c("pre", "on"), paste0("~ on * (", sev, ") + ldisc + on:ldisc", dose), "on:ldisc",
      fails_if = "the demographic gradient of the illness steepens at intubation for a reason other than strain"),
    d("within_triple", c(V, C), c("pre", "on"),
      paste0("~ on * vent * (", sev, ") + on * vent * ldisc", dose), "on:vent:ldisc",
      fails_if = "demographics act through ventilation itself by a path other than strain"),
    d("off_reversal", V, c("post", "on"), paste0("~ on * (", sev, ") + ldisc + on:ldisc", dose), "on:ldisc",
      fails_if = "selection on liberation alive; a demographic gradient that ends at extubation for another reason"),
    purrr::imap_dfr(free, function(term, input)
      bind_rows(
        d(paste0("free_", input), V, "on", paste0("~ ldisc + ", term, " + ", sev, dose), "ldisc",
          fails_if = paste0("an input other than ", input, " moves the marker apart from strain")),
        d(paste0("free_", input, "_did"), c(V, C), "on",
          paste0("~ vent * (", sev, " + ", term, ") + ldisc + vent:ldisc", dose), "vent:ldisc",
          fails_if = paste0("an input other than ", input, " moves the marker differently in the two cohorts")))))
}

# Prepare one design on one panel: everything that does not depend on the outcome.
# Not estimable (ok = FALSE) when an arm has too few patients with the phases. `rows`
# are the panel rows the design reads, in the order of its outcome vector.
# raw = TRUE keeps the uncentred slope columns, which the Poisson model needs.
prepare_design <- function(design, patients, panel, raw = FALSE) {
  arms <- design$arms[[1]]; phases <- design$phases[[1]]
  pid <- match(panel$id, patients$id)             # the patient row of each measurement
  phase_code <- match(panel$phase, c("pre", "on", "post"))
  rows <- which(panel$phase %in% phases & (patients$arm %in% arms)[pid])
  if (length(phases) > 1) {                      # a within-person design needs every phase
    n_phases <- tabulate(unique(pid[rows] * 4L + phase_code[rows]) %/% 4L, nbins = nrow(patients))
    rows <- rows[n_phases[pid[rows]] == length(phases)]
  }
  i <- pid[rows]
  patient_rows <- unique(i)
  n_by_arm <- table(factor(patients$arm[patient_rows], levels = arms))
  empty <- list(ok = FALSE, design = design$design, n_patients = length(patient_rows), n_obs = length(rows))
  if (any(n_by_arm < MIN_PATIENTS)) return(empty)
  # the design's rows of the patient table, built column by column (row-subsetting a
  # data frame is the slow step when this runs thousands of times)
  vars <- all.vars(as.formula(design$formula))
  extra_vars <- c("lv0", if ("laer_on" %in% vars) "laer0", if ("lbmi_on" %in% vars) c("lbmi", "vent"))
  pan <- lapply(patients[vars %>% intersect(names(patients)) %>% union(extra_vars)], `[`, i)
  pan$t <- panel$t[rows]
  pan$on <- as.numeric(panel$phase[rows] == "on")
  pan$lv_on <- pan$lv0 * pan$on
  if ("laer_on" %in% vars) pan$laer_on <- pan$laer0 * pan$on
  if ("lbmi_on" %in% vars) pan$lbmi_on <- pan$lbmi * pan$on * pan$vent
  pan <- structure(pan, class = "data.frame", row.names = .set_row_names(length(rows)))
  X <- model.matrix(as.formula(design$formula), data = pan)
  if (!design$target %in% colnames(X)) {         # R may write an interaction in another order
    alt <- colnames(X)[sorted_term(colnames(X)) == sorted_term(design$target)]
    if (length(alt) != 1) stop("design ", design$design, ": no coefficient ", design$target)
    colnames(X)[colnames(X) == alt] <- design$target
  }
  key <- i * 4L + phase_code[rows]
  g <- match(key, unique(key)); n_g <- tabulate(g)   # the patient-phase of each row
  centre <- function(m) m - (rowsum(m, g, reorder = TRUE) / n_g)[g, , drop = FALSE]   # within patient-phase
  Xt <- X * pan$t                                 # slopes: every column times time
  Z <- centre(Xt)
  keep <- which(colSums(Z^2) > 1e-10)
  qrz <- qr(Z[, keep, drop = FALSE])
  keep <- keep[qrz$pivot[seq_len(qrz$rank)]]
  if (!design$target %in% colnames(Z)[keep]) return(empty)
  Z <- Z[, keep, drop = FALSE]
  ZtZ_inv <- solve(crossprod(Z))
  k <- which(colnames(Z) == design$target)
  # how much of the target's own variation is left beside the other terms
  partial_r2 <- if (ncol(Z) > 1) 1 - 1 / (ZtZ_inv[k, k] * sum(Z[, k]^2)) else 0
  list(ok = TRUE, design = design$design, sign = design$sign, Z = Z, ZtZ_inv = ZtZ_inv, k = k,
       centre = centre, cluster = match(i, patient_rows), rows = rows, group = g,
       patient_ids = patients$id[patient_rows], Zraw = if (raw) Xt[, keep, drop = FALSE] else NULL,
       n_patients = length(patient_rows), n_obs = length(rows), n_by_arm = n_by_arm,
       partial_r2 = partial_r2, df_resid = nrow(Z) - ncol(Z) - max(g))
}

# The coefficients of a prepared design for one outcome vector: beta, each patient's
# influence on it (one row per patient who contributes, `patient` their index in
# prep$patient_ids) and the cluster-robust variance those influences give. NULL when
# the model cannot be fitted.
fit_coefs <- function(prep, y, family = "gaussian") {
  if (family == "gaussian") {
    yc <- prep$centre(matrix(y, ncol = 1))[, 1]
    beta <- prep$ZtZ_inv %*% crossprod(prep$Z, yc)
    resid <- yc - prep$Z %*% beta
    scores <- rowsum(prep$Z * as.numeric(resid), prep$cluster, reorder = TRUE)
    influence <- scores %*% prep$ZtZ_inv
    patient <- as.integer(rownames(scores))
    resid_sd <- sqrt(sum(resid^2) / max(prep$df_resid, 1)); n_obs <- length(y)
  } else if (family == "poisson") {
    # log E[y] = a_group + Z beta. Conditioning on a group's total removes its
    # intercept: within a group the outcomes are shares p = exp(Z beta) / sum of it.
    # Groups with a total of zero say nothing about beta and drop out.
    if (is.null(prep$Zraw)) stop("prepare_design(raw = TRUE) is needed for the Poisson model")
    if (any(y < 0)) stop("the Poisson model needs a non-negative outcome")
    use <- (rowsum(y, prep$group, reorder = TRUE)[, 1] > 0)[prep$group]
    Z <- prep$Zraw[use, , drop = FALSE]; y <- y[use]
    g <- match(prep$group[use], unique(prep$group[use])); cl <- prep$cluster[use]
    total <- rowsum(y, g, reorder = TRUE)[g, 1]
    shares <- function(beta) {
      eta <- as.numeric(Z %*% beta); eta <- eta - (rowsum(eta, g, reorder = TRUE) / tabulate(g))[g, 1]
      m <- exp(eta); m / rowsum(m, g, reorder = TRUE)[g, 1]
    }
    loglik <- function(pr) sum(y[y > 0] * log(pr[y > 0]))
    beta <- rep(0, ncol(Z)); pr <- shares(beta); ll <- loglik(pr); converged <- FALSE
    for (iter in 1:100) {
      mu <- total * pr
      Zc <- Z - rowsum(Z * pr, g, reorder = TRUE)[g, , drop = FALSE]      # centred at the group's fitted mean
      H <- crossprod(Zc * mu, Zc)
      step <- tryCatch(solve(H, crossprod(Zc, y - mu)), error = function(e) NULL)
      if (is.null(step)) return(NULL)
      for (half in 0:20) {                         # the likelihood is concave; halve a step that overshoots
        cand <- beta + step[, 1] / 2^half; pr_c <- shares(cand); ll_c <- loglik(pr_c)
        if (is.finite(ll_c) && ll_c >= ll - 1e-10) break
      }
      done <- max(abs(cand - beta)) < 1e-8
      beta <- cand; pr <- pr_c; ll <- ll_c
      if (done) { converged <- TRUE; break }
    }
    if (!converged) return(NULL)
    mu <- total * pr
    Zc <- Z - rowsum(Z * pr, g, reorder = TRUE)[g, , drop = FALSE]
    H_inv <- tryCatch(solve(crossprod(Zc * mu, Zc)), error = function(e) NULL)
    if (is.null(H_inv)) return(NULL)
    scores <- rowsum(Zc * (y - mu), cl, reorder = TRUE)
    influence <- scores %*% H_inv
    patient <- as.integer(rownames(scores))
    beta <- matrix(beta, ncol = 1)
    resid_sd <- sqrt(sum((y - mu)^2 / pmax(mu, 1e-12)) / max(length(y) - ncol(Z) - max(g), 1))   # Pearson dispersion
    n_obs <- length(y)
  } else stop("unknown family ", family)
  G <- nrow(influence)
  rownames(beta) <- colnames(prep$Z); colnames(influence) <- colnames(prep$Z)
  list(beta = beta[, 1], V = crossprod(influence) * G / (G - 1), influence = influence, patient = patient,
       resid_sd = resid_sd, n_obs = n_obs, n_patients = G)
}

# Fit a prepared design to an outcome vector aligned with the prepared rows: the
# lambda gradient, and beside it the gradient on the index log VT/PBW where the
# design has that term. Under strain, with dose = log VT - log size, the two are the
# same number (vt_gradient = estimate); p_vt_equal tests that. influence = TRUE
# attaches each patient's influence on the estimate (named by patient id).
NO_FIT <- tibble::tibble(estimate = NA_real_, se = NA_real_, resid_sd = NA_real_, vt_gradient = NA_real_,
                         vt_gradient_se = NA_real_, p_vt_equal = NA_real_)
fit_design <- function(prep, y, family = "gaussian", influence = FALSE) {
  fit <- if (prep$ok) fit_coefs(prep, y, family) else NULL
  if (is.null(fit)) return(dplyr::bind_cols(tibble::tibble(design = prep$design), NO_FIT))
  k <- prep$k; v <- which(names(fit$beta) == "lv_on")
  out <- tibble::tibble(design = prep$design, estimate = prep$sign * fit$beta[[k]], se = sqrt(fit$V[k, k]),
                        resid_sd = fit$resid_sd, vt_gradient = NA_real_, vt_gradient_se = NA_real_, p_vt_equal = NA_real_)
  if (length(v) == 1) {
    difference <- fit$beta[[k]] - fit$beta[[v]]
    out$vt_gradient <- fit$beta[[v]]; out$vt_gradient_se <- sqrt(fit$V[v, v])
    out$p_vt_equal <- 2 * pnorm(-abs(difference) / sqrt(fit$V[k, k] + fit$V[v, v] - 2 * fit$V[k, v]))
  }
  if (influence) attr(out, "influence") <- setNames(prep$sign * fit$influence[, k], prep$patient_ids[fit$patient])
  out
}
# the outcome of a panel in a prepared design's row order
align_outcome <- function(prep, panel, column = "y") panel[[column]][prep$rows]

# Fit every design to one panel: one row per design (NA where not estimable)
fit_all_designs <- function(designs, patients, panel, column = "y", family = "gaussian") {
  out <- vector("list", nrow(designs))
  for (i in seq_len(nrow(designs))) {
    prep <- prepare_design(designs[i, ], patients, panel, raw = family == "poisson")
    fit <- fit_design(prep, if (prep$ok) align_outcome(prep, panel, column) else NULL, family)
    fit$n_patients <- prep$n_patients; fit$n_obs <- prep$n_obs
    fit$partial_r2 <- if (prep$ok) prep$partial_r2 else NA_real_
    out[[i]] <- fit
  }
  dplyr::bind_rows(out)
}

# ---- PBW and PFVC as separate terms ----------------------------------------------
# With VT/PBW dropped, lambda alone forces log PBW and log PFVC to carry equal and
# opposite coefficients. Under strain, slope = b (log VT - log PFVC); if tidal volume
# was set as log VT = rho log PBW + other variation, the slope is
# b rho log PBW - b log PFVC: the two are equal and opposite only when rho = 1, VT at
# the same multiple of PBW in everyone. Otherwise the lambda coefficient mixes the
# PFVC effect with a diluted tidal-volume effect. Adding log PBW beside lambda frees
# them: the lambda coefficient is then minus the log PFVC coefficient, the gradient
# along PFVC at a fixed PBW (and so roughly a fixed VT) whatever rho is, and the
# coefficient on log PBW is b (rho - 1), zero when the ratio form is adequate. It also
# takes any direct effect of height or sex on the marker, so a coefficient of zero
# is not proof of rho = 1, nor a non-zero one proof of rho != 1: under strain it
# should be close to estimate x (rho - 1), with rho from vt_follows_pbw(). PBW is
# almost entirely height within sex, so at a fixed PBW lambda moves mainly with the age
# and race terms of PFVC: these designs are less precise and lean on age.

# A design with log PBW freed beside lambda wherever lambda enters (main effect and
# interactions alike). NULL for a design without lambda (the fingerprint) or with the
# index VT/PBW in it (there tidal volume is in the model already).
with_pbw_separate <- function(design) {
  if (!grepl("\\bldisc\\b", design$formula) || grepl("lv_on", design$formula)) return(NULL)
  design$formula <- gsub("\\bldisc\\b", "(ldisc + lpbw)", design$formula)
  design$extra <- gsub("\\bldisc\\b", "lpbw", design$target)
  design
}

# A design's target and one further coefficient of it (design$extra) from the same fit:
# both estimates, their standard errors and covariance, and the test that the extra is 0.
fit_with_extra <- function(design, patients, panel, family = "gaussian") {
  prep <- prepare_design(design, patients, panel, raw = family == "poisson")
  fit <- if (prep$ok) fit_coefs(prep, align_outcome(prep, panel), family) else NULL
  out <- tibble::tibble(design = design$design, estimate = NA_real_, se = NA_real_, extra = NA_real_, extra_se = NA_real_,
                        cov_extra = NA_real_, p_extra_zero = NA_real_, n_patients = prep$n_patients, n_obs = prep$n_obs)
  if (is.null(fit)) return(out)
  e <- which(sorted_term(names(fit$beta)) == sorted_term(design$extra))
  k <- prep$k; s <- prep$sign
  out$estimate <- s * fit$beta[[k]]; out$se <- sqrt(fit$V[k, k])
  if (length(e) != 1) return(out)                 # the extra term dropped as collinear
  out$extra <- s * fit$beta[[e]]; out$extra_se <- sqrt(fit$V[e, e]); out$cov_extra <- fit$V[k, e]
  out$p_extra_zero <- 2 * pnorm(-abs(fit$beta[[e]]) / sqrt(fit$V[e, e]))
  out
}

# Every lambda design of a set refitted with log PBW freed beside lambda. With rho (from
# vt_follows_pbw), strain predicts pbw_extra = (rho - 1) x estimate: pbw_extra_strain is
# that prediction and p_strain_consistent tests the difference (rho taken as known; its
# standard error is small beside the designs'). Where rho is not 1, a non-zero
# pbw_extra is what strain predicts, so p_ratio_adequate alone misleads there.
fit_pbw_separate <- function(designs, patients, panel, family = "gaussian", rho = NA_real_) {
  rows <- lapply(seq_len(nrow(designs)), function(i) with_pbw_separate(designs[i, ]))
  rows <- dplyr::bind_rows(rows)
  if (!nrow(rows)) return(tibble::tibble())
  purrr::map_dfr(seq_len(nrow(rows)), function(i) fit_with_extra(rows[i, ], patients, panel, family)) %>%
    dplyr::mutate(rho = rho, pbw_extra_strain = (rho - 1) * estimate,
                  p_strain_consistent = 2 * stats::pnorm(-abs(extra - pbw_extra_strain) /
                    sqrt(extra_se^2 + (rho - 1)^2 * se^2 - 2 * (rho - 1) * cov_extra))) %>%
    dplyr::rename(pbw_extra = extra, pbw_extra_se = extra_se, p_ratio_adequate = p_extra_zero) %>%
    dplyr::select(-cov_extra)
}

# rho: how closely the index tidal volume followed PBW. The slope of log VT on log PBW
# across ventilated patients, one row per patient (log VT = log VT/PBW + log PBW, so
# rho is 1 plus the slope of log VT/PBW on log PBW). 1 is VT at the same multiple of
# PBW in everyone; below 1, smaller patients got relatively more. Overall and by sex
# (PBW, and so the protocol's tidal volume, differs by sex at the same height).
# In the default cohort, gated on VT/PBW 6-8 mL/kg at the index, the gate truncates
# log VT/PBW and pulls rho toward 1: there rho describes the analysed sample, not how
# tidal volume was set. Run with PBWPFVC_VTPBW_GATE=0 for that.
vt_follows_pbw <- function(patients) {
  v <- patients[patients$vent == 1 & is.finite(patients$lv0) & is.finite(patients$lpbw), , drop = FALSE]
  one <- function(d, group) {
    if (nrow(d) < 10 || stats::var(d$lpbw) == 0) return(tibble::tibble(group = group, n_patients = nrow(d)))
    f <- stats::lm(lv0 ~ lpbw, data = d); b <- stats::coef(f)[["lpbw"]]; se <- sqrt(stats::vcov(f)["lpbw", "lpbw"])
    tibble::tibble(group = group, n_patients = nrow(d), rho = 1 + b, se = se,
                   conf_low = 1 + b - stats::qnorm(0.975) * se, conf_high = 1 + b + stats::qnorm(0.975) * se,
                   p_rho_one = 2 * stats::pnorm(-abs(b) / se), sd_log_pbw = stats::sd(d$lpbw), sd_log_vtpbw = stats::sd(d$lv0))
  }
  dplyr::bind_rows(one(v, "all"),
                   purrr::map_dfr(levels(droplevels(factor(v$sex_category))), function(x) one(v[v$sex_category == x, ], as.character(x))))
}

# ---- Does PFVC mediate the age gradient? ------------------------------------------
# PFVC is a fixed function of age, height, sex and race, so given those it does not
# vary: the share of age's association with the marker that runs through PFVC cannot
# be read off the data alone. It can under one assumption, the one strain itself makes:
# a log unit of PFVC acts the same whatever made it small. lambda is split into its age
# part (lage: GLI's age term) and the rest (lrest: height, sex, race and Devine PBW),
# each with its own coefficient. b_rest is the PFVC gradient learnt from non-age
# variation; b_age is the gradient along lambda's age part, which also carries any
# effect of age that does not run through PFVC. Under the assumption, the share of the
# age gradient that PFVC accounts for is b_rest / b_age (1: all of it, as strain alone
# would have; 0: none), and p_full_mediation tests b_rest = b_age. lrest is open to
# direct effects of sex, race and height as the ladder rungs are; read the control
# contrast, where effects common to both cohorts cancel.
age_mediation_design <- function(design) {
  if (!grepl("\\bldisc\\b", design$formula) || grepl("ns\\(age10", design$formula)) return(NULL)
  design$formula <- gsub("\\bldisc\\b", "(lrest + lage)", design$formula)
  design$extra <- gsub("\\bldisc\\b", "lage", design$target)
  design$target <- gsub("\\bldisc\\b", "lrest", design$target)
  design
}
fit_age_mediation <- function(designs, patients, panel, family = "gaussian") {
  rows <- dplyr::bind_rows(lapply(seq_len(nrow(designs)), function(i) age_mediation_design(designs[i, ])))
  if (!nrow(rows)) return(tibble::tibble())
  v <- patients$vent == 1
  per_decade <- unname(stats::coef(stats::lm(patients$lage[v] ~ patients$age10[v]))[2])   # lage per decade of age
  purrr::map_dfr(seq_len(nrow(rows)), function(i) fit_with_extra(rows[i, ], patients, panel, family)) %>%
    dplyr::transmute(design, b_rest = estimate, b_rest_se = se, b_age = extra, b_age_se = extra_se,
                     p_full_mediation = 2 * stats::pnorm(-abs(b_rest - b_age) / sqrt(se^2 + extra_se^2 - 2 * cov_extra)),
                     share_mediated = b_rest / b_age,
                     share_mediated_se = sqrt(se^2 / b_age^2 + b_rest^2 * extra_se^2 / b_age^4 - 2 * b_rest * cov_extra / b_age^3),
                     lage_per_decade = per_decade,
                     age_gradient_per_decade = b_age * per_decade, through_pfvc_per_decade = b_rest * per_decade,
                     n_patients, n_obs)
}

# ---- Aeration: compliance as a measure of the lung that takes the breath -----------
# With specific lung elastance (the pressure that stretches the lung by a given
# fraction of its resting volume) taken as one value, compliance is proportional to the
# aerated end-expiratory volume, and the driving pressure splits exactly into
#     log dP = log VT/PBW + lambda - log f + constant,   log f = log(Crs / PFVC) + c,
# f the fraction of the predicted lung still aerated. Strain gives the three equal
# weight: the gradients on log VT/PBW and lambda equal, and the one on log(Crs/PFVC)
# their negative (p_strain_restriction, a 2-df Wald test). Because Crs and VT/PBW come
# from the same set tidal volume at the index, log(Crs/PFVC) = lambda + log VT/PBW -
# log dP exactly, so that test asks: given the driving pressure, do PBW/PFVC and VT/PBW
# add anything? It can also reject under strain: log(Crs/PFVC) measures aeration with
# error (below), which attenuates its coefficient. The quantity to read is the shift of
# the VT/PBW gradient from base_vt_gradient to vt_gradient: if tidal volume was lowered
# for stiff lungs, conditioning on aeration moves it toward the lambda gradient.
# Limits, read with the numbers:
#   - specific elastance is similar on AVERAGE across groups (Chiumello et al.), not
#     constant: in that cohort it spans about two-fold from low to high in every
#     group. log(Crs / PFVC) estimates log f with an error of that size (a two-fold
#     range is log 2 = 0.69; if it covers about four SDs, an SD near 0.17 log units:
#     a derived figure, not one reported), so the adjustment is incomplete and the aeration
#     gradient is attenuated; and if specific elastance moves with age or sex, the
#     term carries a demographic part too.
#   - Crs is the respiratory system's, not the lung's: the chest wall (obesity,
#     abdominal pressure) lowers it without loss of aeration. log BMI enters as a
#     covariate and, while ventilated, beside the aeration term.
#   - plateau pressure needs an inspiratory hold in a volume-controlled mode, so Crs
#     exists for a subset; every design is also fitted without the aeration terms on
#     that same subset (base_estimate, base_vt_gradient) so the two are comparable.
#   - the index is not before ventilation: it falls at or after its start (t0), within
#     hours for most patients (median_index_lag_h, and 90th centile, in the table), so
#     Crs is measured early in the outcome window, not before it. A larger VT for the
#     size can overdistend the lung and lower Crs, and the first hours of ventilation
#     can already have changed it: a covariate measured early, not a baseline.
aeration_design <- function(design) {
  if (!grepl("lv_on", design$formula) || !grepl("\\bldisc\\b", design$formula)) return(NULL)
  design$formula <- paste0(design$formula, " + laer_on + lbmi + lbmi_on")
  design
}
# The aeration designs on the patients with what they need (compliance in the
# ventilated, BMI in both), each beside the same design without the aeration terms.
fit_aeration <- function(designs, patients, panel, family = "gaussian") {
  pats <- patients[!is.na(patients$lbmi) & !is.na(patients$laer0), , drop = FALSE]
  pan <- panel[panel$id %in% pats$id, , drop = FALSE]
  coef_of <- function(fit, term) { j <- which(sorted_term(names(fit$beta)) == sorted_term(term)); if (length(j) == 1) j else NA_integer_ }
  purrr::map_dfr(seq_len(nrow(designs)), function(r) {
    d <- aeration_design(designs[r, ])
    if (is.null(d)) return(NULL)
    out <- tibble::tibble(design = d$design, estimate = NA_real_, se = NA_real_, vt_gradient = NA_real_, vt_gradient_se = NA_real_,
                          aeration_gradient = NA_real_, aeration_gradient_se = NA_real_, p_strain_restriction = NA_real_,
                          base_estimate = NA_real_, base_se = NA_real_, base_vt_gradient = NA_real_,
                          n_patients = NA_integer_, n_with_compliance = sum(pats$vent == 1),
                          median_index_lag_h = stats::median(pats$index_lag_h[pats$vent == 1], na.rm = TRUE),
                          p90_index_lag_h = stats::quantile(pats$index_lag_h[pats$vent == 1], 0.9, na.rm = TRUE, names = FALSE))
    base <- fit_all_designs(designs[r, ], pats, pan, family = family)
    out$base_estimate <- base$estimate; out$base_se <- base$se; out$base_vt_gradient <- base$vt_gradient
    prep <- prepare_design(d, pats, pan, raw = family == "poisson")
    out$n_patients <- prep$n_patients
    fit <- if (prep$ok) fit_coefs(prep, align_outcome(prep, pan), family) else NULL
    if (is.null(fit)) return(out)
    k <- prep$k; v <- coef_of(fit, "lv_on"); a <- coef_of(fit, "laer_on"); s <- prep$sign
    out$estimate <- s * fit$beta[[k]]; out$se <- sqrt(fit$V[k, k])
    if (!is.na(v)) { out$vt_gradient <- s * fit$beta[[v]]; out$vt_gradient_se <- sqrt(fit$V[v, v]) }
    if (!is.na(a)) { out$aeration_gradient <- s * fit$beta[[a]]; out$aeration_gradient_se <- sqrt(fit$V[a, a]) }
    if (!is.na(v) && !is.na(a)) {
      # strain: beta_lambda = beta_VT and beta_aeration = -beta_lambda
      L <- matrix(0, 2, length(fit$beta)); L[1, k] <- 1; L[1, v] <- -1; L[2, k] <- 1; L[2, a] <- 1
      dlt <- L %*% fit$beta; W <- tryCatch(as.numeric(t(dlt) %*% solve(L %*% fit$V %*% t(L)) %*% dlt), error = function(e) NA_real_)
      out$p_strain_restriction <- stats::pchisq(W, df = 2, lower.tail = FALSE)
    }
    out
  })
}

# Which predicted volume does compliance scale with? log Crs on log PBW and the log of
# each predicted volume, with log BMI for the chest wall, in all ventilated patients
# with compliance and in the least injured (S/F at or above the median and PEEP at or
# below 8), aged 80 or under so that every volume is defined on one sample. A coefficient near 1 on the volume and near 0 on PBW says compliance scales
# with that volume and not with PBW. The individual spread of specific elastance is in
# the residual SD; what is estimated here is how the AVERAGE scales.
compliance_scaling <- function(patients) {
  vols <- c(pfvc = "pfvc", pfvc_age25 = "pfvc_age25", tlc = "ptlc", alveolar_volume = "pva", frc = "pfrc")
  # one sample for every volume (TLC and FRC are defined to age 80), so that the AICs compare
  ok <- patients$vent == 1 & !is.na(patients$crs) & patients$crs > 0 & !is.na(patients$lbmi) &
    stats::complete.cases(patients[unname(vols)])
  v <- patients[ok, , drop = FALSE]
  if (nrow(v) < 30) return(tibble::tibble())
  v$least_injured <- v$sf_ratio_raw >= stats::median(v$sf_ratio_raw, na.rm = TRUE) & !is.na(v$peep_set) & v$peep_set <= 8
  purrr::map_dfr(c("all", "least injured"), function(subset) {
    d <- if (subset == "all") v else v[v$least_injured %in% TRUE, , drop = FALSE]
    purrr::imap_dfr(vols, function(col, name) {
      dd <- d[is.finite(log(d[[col]])), , drop = FALSE]
      if (nrow(dd) < 30) return(NULL)
      f <- stats::lm(log(crs) ~ log(pbw) + log(dd[[col]]) + lbmi, data = dd)
      cf <- stats::coef(summary(f))
      tibble::tibble(subset = subset, volume = name, n_patients = nrow(dd),
                     volume_coef = cf[3, 1], volume_se = cf[3, 2], pbw_coef = cf[2, 1], pbw_se = cf[2, 2],
                     bmi_coef = cf[4, 1], bmi_se = cf[4, 2], residual_sd = summary(f)$sigma, r2 = summary(f)$r.squared,
                     aic = stats::AIC(f))
    })
  })
}

# ---- The time course ------------------------------------------------------------
# The designs above summarise the ventilator days with one straight line. A marker
# that responds within a day or two and then plateaus, or one that starts late, is
# poorly described by it, and markers differ. fit_time_course() makes no assumption
# about the shape: it cuts the on-ventilator phase at CAUSAL_KNOTS, fits the design
# with the on phase restricted to each segment in turn (a segment has its own
# intercept; the comparison phases and arms are unchanged), and adds the segments up:
#   segment      the lambda gradient of the slope within the segment
#   cumulative   sum of (segment gradient x segment length) up to each cut: the lambda
#                gradient of the marker's total change by that day. Its standard error
#                comes from the patients' influences added across segments, so it
#                allows for one patient contributing to several.
# A fast marker shows its effect in the first segment and nothing after; a slow one
# the reverse; the cumulative effect at the end of the window is comparable between
# them. Later segments hold only the patients still ventilated and measured then.
fit_time_course <- function(design, patients, panel, family = "gaussian", knots = CAUSAL_KNOTS) {
  cuts <- c(0, knots[knots > 0 & knots < ON_DAYS], ON_DAYS)
  on <- panel$phase == "on"
  running <- NULL; estimable <- TRUE; rows <- list()
  for (s in seq_len(length(cuts) - 1)) {
    a <- cuts[s]; b <- cuts[s + 1]
    seg <- complete_phases(panel[!on | (panel$t >= a & panel$t < b), , drop = FALSE])
    prep <- prepare_design(design, patients, seg, raw = family == "poisson")
    fit <- fit_design(prep, if (prep$ok) align_outcome(prep, seg) else NULL, family, influence = TRUE)
    window <- sprintf("days %g to %g", a, b)
    rows[[length(rows) + 1]] <- tibble::tibble(design = design$design, quantity = "segment", window = window, from = a, to = b,
                                               estimate = fit$estimate, se = fit$se, n_patients = prep$n_patients, n_obs = prep$n_obs)
    estimable <- estimable && !is.na(fit$estimate)
    if (estimable) {
      infl <- attr(fit, "influence") * (b - a)
      ids <- union(names(running$influence), names(infl))
      total <- setNames(numeric(length(ids)), ids)
      total[names(running$influence)] <- running$influence
      total[names(infl)] <- total[names(infl)] + infl
      running <- list(estimate = sum(running$estimate, fit$estimate * (b - a)), influence = total)
    }
    G <- length(running$influence)
    rows[[length(rows) + 1]] <- tibble::tibble(design = design$design, quantity = "cumulative", window = sprintf("by day %g", b), from = 0, to = b,
                                               estimate = if (estimable) running$estimate else NA_real_,
                                               se = if (estimable) sqrt(sum(running$influence^2) * G / (G - 1)) else NA_real_,
                                               n_patients = G, n_obs = NA_integer_)
  }
  dplyr::bind_rows(rows)
}

# ---- The decision rule --------------------------------------------------------------
# One verdict from one site's fifteen estimates, fixed before any data are read
# (CAUSAL_TESTBED_PROTOCOL.md, section 6; scored by causal_testbed.R). The data support the hypothesis only if
# ALL of these hold, in the direction of injury at two-sided 5%:
#   1. between cohorts    control_did is significant
#   2. within patients    off_reversal or within_onoff is significant, and neither is
#                         significant the other way
#   3. not age, not sex   free_age_did and free_sex_did are each significant
#   4. the rungs agree    all four free_*_did rungs have the sign of injury and lie
#                         within a factor of two of control_did
# A design that is not estimable counts as not significant. The rule is scored as a
# whole: how often it says "supports" under each truth.
RULE_Z <- qnorm(0.975)
RULE_BAND <- 2
triangulation_supports <- function(fits, sign) {
  z <- setNames(sign * fits$estimate / fits$se, fits$design)
  est <- setNames(sign * fits$estimate, fits$design)
  up   <- function(d) isTRUE(z[[d]] > RULE_Z)
  down <- function(d) isTRUE(z[[d]] < -RULE_Z)
  rungs <- paste0("free_", c("age", "sex", "race", "height"), "_did")
  ratio <- est[rungs] / est[["control_did"]]
  c(rule = up("control_did") &&
      (up("off_reversal") || up("within_onoff")) && !down("off_reversal") && !down("within_onoff") &&
      up("free_age_did") && up("free_sex_did") &&
      all(!is.na(ratio)) && all(est[rungs] > 0) && all(ratio >= 1 / RULE_BAND & ratio <= RULE_BAND),
    # the comparators: what a reader concludes from one familiar analysis alone
    unadjusted_alone = up("unadjusted"), control_did_alone = up("control_did"))
}
# The rule with the index VT/PBW held and with it dropped. Holding VT/PBW is what "at a
# fixed VT/PBW" means, but where tidal volume is chosen by how sick the lung is, VT/PBW
# is a consequence of both lambda and severity, and holding it lets severity into the
# lambda gradient. Dropping it avoids that and changes the estimand to the gradient as
# practice delivers it. Neither is right under every way of setting tidal volume, so
# both are reported and the hypothesis is supported only if the rule holds both ways.
triangulation_supports_both <- function(fits_held, fits_dropped, sign) {
  held <- triangulation_supports(fits_held, sign); dropped <- triangulation_supports(fits_dropped, sign)
  c(rule_both = held[["rule"]] && dropped[["rule"]], rule_held = held[["rule"]], rule_dropped = dropped[["rule"]],
    unadjusted_alone = held[["unadjusted_alone"]], control_did_alone = held[["control_did_alone"]])
}

# ---- Injury and death as states (causal_illness_death.R) ---------------------------
# A threshold turns a marker into a state a patient enters once; death is then the
# state after it, or the state reached without it. See causal_illness_death.R for the
# reading. Settings:
HORIZON_DAYS   <- as.numeric(Sys.getenv("PBWPFVC_CAUSAL_HORIZON", "28"))    # death is followed this long
INJURY_GRACE_H <- as.numeric(Sys.getenv("PBWPFVC_CAUSAL_GRACE_H", "6"))     # labs up to here are baseline
AKI_STAGE      <- as.integer(Sys.getenv("PBWPFVC_CAUSAL_AKI_STAGE", "2"))   # KDIGO stage counted as injury: 2 or 1
PLT_BELOW      <- as.numeric(Sys.getenv("PBWPFVC_CAUSAL_PLT_BELOW", "100"))
PLT_FALL       <- as.numeric(Sys.getenv("PBWPFVC_CAUSAL_PLT_FALL", "0.5"))
PRESSOR_MIN    <- as.numeric(Sys.getenv("PBWPFVC_CAUSAL_PRESSOR_MIN", "0"))
MIN_EVENTS     <- 20L           # a transition with fewer events in an arm is not estimated
if (!AKI_STAGE %in% 1:2) stop("PBWPFVC_CAUSAL_AKI_STAGE must be 1 or 2")
INJURY_MARKERS <- c(
  creatinine = if (AKI_STAGE == 2) "KDIGO stage 2 or worse: creatinine at least 2x baseline, or renal replacement started" else
    "KDIGO stage 1 or worse: creatinine at least 1.5x baseline or up 0.3 mg/dL within 48 h, or renal replacement started",
  platelets = paste0("platelet count below ", PLT_BELOW, ", or a fall of at least ", 100 * PLT_FALL, "% from baseline"),
  pressor = if (PRESSOR_MIN > 0) paste0("norepinephrine-equivalent rate of at least ", PRESSOR_MIN, " mcg/kg/min after t0") else
    "any vasopressor after t0")
ILLNESS_DEATH_TRANSITIONS <- c("injury", "death without injury", "death after injury", "death, any")
CONTRAST <- "contrast: injury minus death without injury"

# When each patient is first injured: one row per patient with
#   eligible   has a baseline and is not already past the threshold at the end of it
#   reason     why not, otherwise "eligible"
#   t_injury   days from t0 of the first crossing in the on phase; NA if none
# `values` (id, t, value) may be given in place of the site's files (the testbed does).
load_injury_times <- function(config, patients, marker, values = NULL, rrt = NULL) {
  grace <- INJURY_GRACE_H / 24
  who <- patients %>% transmute(id, t0, on_end = pmin(dplyr::coalesce(t_off_d, Inf), ON_DAYS))
  days_from_t0 <- function(d, time) d %>% inner_join(who %>% select(id, t0), by = "id") %>%
    mutate(t = as.numeric(difftime(.data[[time]], t0, units = "days"))) %>% select(-t0)
  if (marker == "pressor") {
    if (is.null(values)) values <- purrr::map_dfr(both_cohorts(config, "ne_equiv_admin.parquet"), function(f)
      arrow::read_parquet(f) %>% transmute(id = hospitalization_id, dttm = admin_dttm, value = ne_equiv_total)) %>% days_from_t0("dttm")
    running <- values %>% filter(is.finite(value), value > 0, value >= PRESSOR_MIN)
    out <- who %>%
      left_join(running %>% filter(t >= -1, t < 0) %>% distinct(id) %>% mutate(before = TRUE), by = "id") %>%
      left_join(running %>% filter(t >= 0) %>% group_by(id) %>% summarise(first = min(t), .groups = "drop"), by = "id") %>%
      mutate(eligible = is.na(before), reason = if_else(eligible, "eligible", "already injured at baseline"))
  } else {
    lab <- c(creatinine = "creatinine", platelets = "platelet_count")[[marker]]
    if (is.null(values)) values <- purrr::map_dfr(both_cohorts(config, "cohort_labs_clean.parquet"), function(f) {
      labs <- arrow::read_parquet(f)
      col <- intersect(c("lab_value_numeric", "lab_value"), names(labs))[1]
      labs %>% filter(lab_category == lab) %>%
        transmute(id = hospitalization_id, dttm = lab_result_dttm, value = suppressWarnings(as.numeric(.data[[col]])))
    }) %>% days_from_t0("dttm")
    values <- values %>% filter(is.finite(value), value > 0) %>% arrange(id, t)
    best <- if (marker == "creatinine") min else max
    base <- values %>% filter(t <= grace) %>% group_by(id) %>%
      summarise(baseline = best(value), last = dplyr::last(value), .groups = "drop")
    later <- values %>% filter(t > grace) %>% inner_join(base, by = "id")
    if (marker == "creatinine") {
      ratio <- if (AKI_STAGE == 2) 2 else 1.5
      base <- base %>% mutate(already = last >= ratio * baseline)
      later <- later %>% mutate(crossed = value >= ratio * baseline)
      if (AKI_STAGE == 1) {                       # or up 0.3 mg/dL on the lowest value of the 48 hours before
        low48 <- function(t, v) vapply(seq_along(t), function(k) { w <- t < t[k] & t >= t[k] - 2; if (any(w)) min(v[w]) else NA_real_ }, numeric(1))
        rise <- values %>% group_by(id) %>% mutate(low = low48(t, value)) %>% ungroup() %>%
          filter(!is.na(low), value - low >= 0.3)
        base <- base %>% mutate(already = already | id %in% rise$id[rise$t <= grace])
        later <- later %>% left_join(rise %>% transmute(id, t, rose = TRUE), by = c("id", "t")) %>%
          mutate(crossed = crossed | dplyr::coalesce(rose, FALSE))
      }
      # renal replacement: script 01's CRRT and dialysis records, where the site has them
      if (is.null(rrt) && !is.null(config)) rrt <- purrr::map_dfr(c("cohort_crrt.parquet", "cohort_dialysis.parquet"), function(file)
        purrr::map_dfr(c(config$output_dir, file.path(config$output_dir, "controls", "nosupport")), function(dir) {
          f <- file.path(dir, file)
          if (!file.exists(f)) return(NULL)
          arrow::read_parquet(f) %>% transmute(id = hospitalization_id, dttm = recorded_dttm)
        })) %>% { if (nrow(.)) days_from_t0(., "dttm") else tibble(id = character(), t = numeric()) }
      if (!is.null(rrt) && nrow(rrt)) {
        base <- base %>% mutate(already = already | id %in% rrt$id[rrt$t <= grace])
        later <- bind_rows(later, rrt %>% filter(t > grace, id %in% base$id) %>% mutate(crossed = TRUE))
      }
    } else {
      base <- base %>% mutate(already = last < PLT_BELOW | last <= (1 - PLT_FALL) * baseline)
      later <- later %>% mutate(crossed = value < PLT_BELOW | value <= (1 - PLT_FALL) * baseline)
    }
    out <- who %>% left_join(base %>% select(id, baseline, already), by = "id") %>%
      left_join(later %>% filter(crossed) %>% group_by(id) %>% summarise(first = min(t), .groups = "drop"), by = "id") %>%
      mutate(eligible = !is.na(baseline) & !already,
             reason = case_when(is.na(baseline) ~ "no baseline value", already ~ "already injured at baseline", TRUE ~ "eligible"))
  }
  out %>% transmute(id, eligible, reason, t_injury = if_else(eligible & !is.na(first) & first < on_end, first, NA_real_))
}

# One row per patient with the three states' times, days from t0:
#   t_injury   first injury, if before death and the end of follow-up
#   t_death    death within the horizon and the follow-up; NA otherwise
#   t_end      end of follow-up: death, the horizon, or (control) escalation of support
illness_death_states <- function(patients, injury, horizon = HORIZON_DAYS) {
  patients %>% inner_join(injury, by = "id") %>%
    mutate(follow = if_else(vent == 1, horizon, pmin(horizon, dplyr::coalesce(t_off_d, Inf))),
           t_death = if_else(!is.na(death_d) & death_d <= follow, pmax(death_d, 1e-3), NA_real_),
           t_end = pmin(follow, dplyr::coalesce(t_death, Inf)),
           t_injury = if_else(!is.na(t_injury) & t_injury < t_end, pmax(t_injury, 1e-4), NA_real_)) %>%
    filter(t_end > 0)
}
illness_death_counts <- function(states) {
  states %>% group_by(arm) %>%
    summarise(n_patients = n(), no_baseline = sum(reason == "no baseline value"),
              already_injured = sum(reason == "already injured at baseline"), eligible = sum(eligible),
              injured = sum(!is.na(t_injury[reason == "eligible"])),
              died_without_injury = sum(reason == "eligible" & is.na(t_injury) & !is.na(t_death)),
              died_after_injury = sum(reason == "eligible" & !is.na(t_injury) & !is.na(t_death)), .groups = "drop")
}
# the counting-process rows of one transition, eligible patients only
transition_rows <- function(states, transition) {
  s <- states[states$eligible, , drop = FALSE]
  on_end <- pmin(dplyr::coalesce(s$t_off_d, Inf), ON_DAYS)
  injured <- !is.na(s$t_injury); died <- !is.na(s$t_death)
  s$start <- 0
  if (transition == "injury") {
    s$stop <- pmin(dplyr::coalesce(s$t_injury, Inf), on_end, s$t_end); s$event <- as.numeric(injured)
  } else if (transition == "death without injury") {
    s$stop <- pmin(dplyr::coalesce(s$t_injury, Inf), s$t_end); s$event <- as.numeric(died & !injured)
  } else if (transition == "death after injury") {
    s <- s[injured, , drop = FALSE]; s$start <- s$t_injury; s$stop <- s$t_end; s$event <- as.numeric(!is.na(s$t_death))
  } else if (transition == "death, any") {
    s$stop <- s$t_end; s$event <- as.numeric(died)
  } else stop("unknown transition ", transition)
  s[s$stop > s$start, , drop = FALSE]
}
# The lambda gradient of every transition's hazard under every single-phase design: a
# Cox model stratified by cohort, robust standard errors from each patient's influence.
fit_illness_death <- function(designs, states) {
  designs <- designs[vapply(designs$phases, function(p) identical(p, "on"), logical(1)), ]
  states$lv_on <- states$lv0; states$on <- 1
  one <- function(design, transition) {
    arms <- design$arms[[1]]
    d <- transition_rows(states[states$arm %in% arms, , drop = FALSE], transition)
    out <- tibble::tibble(design = design$design, transition = transition, estimate = NA_real_, se = NA_real_,
                          n_patients = nrow(d), n_events = sum(d$event))
    events <- table(factor(d$arm[d$event == 1], levels = arms))
    if (any(events < MIN_EVENTS)) return(list(row = out, influence = NULL))
    # the cohort is a stratum, so its own term goes; everything else is the design's formula
    rhs <- update(as.formula(design$formula), ~ . - vent)
    form <- as.formula(paste("survival::Surv(start, stop, event) ~", paste(deparse(rhs[[2]], width.cutoff = 500), collapse = " "),
                             "+ survival::strata(vent)"))
    fit <- tryCatch(suppressWarnings(survival::coxph(form, data = d, ties = "efron")), error = function(e) NULL)
    k <- if (is.null(fit)) integer() else which(sorted_term(names(coef(fit))) == sorted_term(design$target))
    if (length(k) != 1 || is.na(coef(fit)[[k]])) return(list(row = out, influence = NULL))
    infl <- residuals(fit, type = "dfbeta"); if (is.null(dim(infl))) infl <- matrix(infl, ncol = 1)
    out$estimate <- coef(fit)[[k]]; out$se <- sqrt(sum(infl[, k]^2))
    list(row = out, influence = setNames(infl[, k], d$id))
  }
  dplyr::bind_rows(lapply(seq_len(nrow(designs)), function(i) {
    fits <- lapply(setNames(ILLNESS_DEATH_TRANSITIONS, ILLNESS_DEATH_TRANSITIONS), function(tr) one(designs[i, ], tr))
    rows <- dplyr::bind_rows(lapply(fits, `[[`, "row"))
    a <- fits[["injury"]]; b <- fits[["death without injury"]]
    contrast <- tibble::tibble(design = designs$design[i], transition = CONTRAST, estimate = NA_real_, se = NA_real_,
                               n_patients = a$row$n_patients, n_events = NA_integer_)
    if (!is.null(a$influence) && !is.null(b$influence)) {
      ids <- union(names(a$influence), names(b$influence)); diff <- setNames(numeric(length(ids)), ids)
      diff[names(a$influence)] <- a$influence; diff[names(b$influence)] <- diff[names(b$influence)] - b$influence
      contrast$estimate <- a$row$estimate - b$row$estimate; contrast$se <- sqrt(sum(diff^2))
    }
    dplyr::bind_rows(rows, contrast)
  }))
}
