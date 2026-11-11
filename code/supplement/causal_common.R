# =============================================================================
# causal_common: the within-person and ladder designs, and the planted-truth
# generator that tests them. Sourced by causal_triangulation.R and
# causal_testbed.R; never run; writes nothing.
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
# and each can be scored against a planted truth (causal_testbed.R).
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

# marker -> lab_category in cohort_labs_clean (script 02), and the direction of injury
CAUSAL_MARKERS <- tibble::tribble(
  ~marker,      ~lab_category,    ~injury_sign,
  "platelets",  "platelet_count", -1,
  "creatinine", "creatinine",     +1)

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
            "race_category", "height_cm", "pbw", "pfvc", "sf_ratio", "sofa_total", "discharge_dttm")
  vent <- arrow::read_parquet(vent_file) %>%
    select(all_of(keep), vtpbw, last_imv_dttm) %>% mutate(arm = "ventilated")
  ctrl <- arrow::read_parquet(ctrl_file) %>%
    select(all_of(keep), escalation_dttm) %>% mutate(arm = "no support", vtpbw = NA_real_)
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
           age10 = age_at_admission / 10,
           sex_category = factor(as.character(sex_category), levels = c("Male", "Female")),
           race_category = factor(as.character(race_category), levels = c("WHITE", "BLACK", "OTHER")),
           ldisc = log(pbw / pfvc),
           fp = height_fingerprint(height_cm, as.character(sex_category)),
           # log VT/PBW about the ventilated mean; 0 in the control, which has no tidal volume
           lv0 = if_else(vent == 1, log(vtpbw) - mean(log(vtpbw[vent == 1]), na.rm = TRUE), 0),
           sofa_z = as.numeric(scale(sofa_total)), sf_z = as.numeric(scale(sf_ratio)),
           age_z = as.numeric(scale(age10))) %>%
    mutate(ldisc = ldisc - mean(ldisc[vent == 1])) %>%
    select(id, patient_id, arm, vent, t0, t_off_d, stay_end_d, age10, age_z, sex_category, race_category,
           height_cm, ldisc, fp, lv0, sofa_z, sf_z)
}

# Every measurement of one marker in the window, with its phase. t is days from t0.
load_causal_panel <- function(config, patients, lab_category) {
  read_labs <- function(dir) {
    f <- file.path(dir, "cohort_labs_clean.parquet")
    if (!file.exists(f)) stop("missing ", f, ": run script 02 first")
    labs <- arrow::read_parquet(f)
    value <- intersect(c("lab_value_numeric", "lab_value"), names(labs))[1]
    labs %>% filter(lab_category == !!lab_category) %>%
      transmute(id = hospitalization_id, lab_dttm = lab_result_dttm, value = suppressWarnings(as.numeric(.data[[value]])))
  }
  labs <- bind_rows(read_labs(config$output_dir), read_labs(file.path(config$output_dir, "controls", "nosupport")))
  labs %>% inner_join(patients %>% select(id, vent, t0, t_off_d), by = "id") %>%
    mutate(t = as.numeric(difftime(lab_dttm, t0, units = "days")),
           on_end = pmin(coalesce(t_off_d, Inf), ON_DAYS),
           phase = case_when(t >= -PRE_DAYS & t < 0 ~ "pre",
                             t >= 0 & t < on_end ~ "on",
                             vent == 1 & !is.na(t_off_d) & t_off_d < ON_DAYS + POST_DAYS &
                               t >= t_off_d & t < t_off_d + POST_DAYS ~ "post",
                             TRUE ~ NA_character_)) %>%
    filter(!is.na(phase), is.finite(value), value > 0) %>%
    transmute(id, phase, t, y = log(value)) %>%
    distinct() %>%
    group_by(id, phase) %>% filter(n() >= 2, max(t) - min(t) >= MIN_SPAN_D) %>% ungroup() %>%
    arrange(id, phase, t)
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
# Returns NULL-estimable (ok = FALSE) when too few patients have the phases.
prepare_design <- function(design, patients, panel) {
  arms <- design$arms[[1]]; phases <- design$phases[[1]]
  pat <- patients %>% filter(arm %in% arms)
  pan <- panel %>% filter(phase %in% phases, id %in% pat$id)
  if (length(phases) > 1) {                      # a within-person design needs both phases
    both <- pan %>% distinct(id, phase) %>% count(id) %>% filter(n == length(phases))
    pan <- pan %>% filter(id %in% both$id)
  }
  n_by_arm <- pan %>% distinct(id) %>% inner_join(pat %>% select(id, arm), by = "id") %>% count(arm)
  empty <- list(ok = FALSE, design = design$design, n_patients = n_distinct(pan$id), n_obs = nrow(pan))
  if (nrow(n_by_arm) < length(arms) || any(n_by_arm$n < MIN_PATIENTS)) return(empty)
  pan <- pan %>% inner_join(pat, by = "id") %>%
    mutate(on = as.numeric(phase == "on"), lv_on = lv0 * on,
           group = paste(id, phase))
  X <- model.matrix(as.formula(design$formula), data = pan)
  if (!design$target %in% colnames(X)) {         # R may write an interaction in another order
    alt <- colnames(X)[sorted_term(colnames(X)) == sorted_term(design$target)]
    if (length(alt) != 1) stop("design ", design$design, ": no coefficient ", design$target)
    colnames(X)[colnames(X) == alt] <- design$target
  }
  g <- as.integer(factor(pan$group)); n_g <- tabulate(g)
  centre <- function(m) m - (rowsum(m, g, reorder = TRUE) / n_g)[g, , drop = FALSE]   # within patient-phase
  Z <- centre(X * pan$t)                          # slopes: every column times time
  keep <- which(apply(Z, 2, function(z) sum(z^2)) > 1e-10)
  qrz <- qr(Z[, keep, drop = FALSE])
  keep <- keep[qrz$pivot[seq_len(qrz$rank)]]
  if (!design$target %in% colnames(Z)[keep]) return(empty)
  Z <- Z[, keep, drop = FALSE]
  ZtZ_inv <- solve(crossprod(Z))
  k <- which(colnames(Z) == design$target)
  # how much of the target's own variation is left beside the other terms
  partial_r2 <- if (ncol(Z) > 1) 1 - 1 / (ZtZ_inv[k, k] * sum(Z[, k]^2)) else 0
  list(ok = TRUE, design = design$design, sign = design$sign, Z = Z, ZtZ_inv = ZtZ_inv, k = k,
       centre = centre, cluster = as.integer(factor(pan$id)), row_key = paste(pan$id, pan$phase, pan$t),
       n_patients = n_distinct(pan$id), n_obs = nrow(pan), n_by_arm = n_by_arm,
       partial_r2 = partial_r2, df_resid = nrow(Z) - ncol(Z) - max(g))
}

# Fit a prepared design to an outcome vector aligned with the prepared rows.
fit_design <- function(prep, y) {
  if (!prep$ok) return(tibble(design = prep$design, estimate = NA_real_, se = NA_real_))
  yc <- prep$centre(matrix(y, ncol = 1))[, 1]
  beta <- prep$ZtZ_inv %*% crossprod(prep$Z, yc)
  resid <- yc - prep$Z %*% beta
  scores <- rowsum(prep$Z * as.numeric(resid), prep$cluster, reorder = FALSE)
  G <- nrow(scores)
  V <- prep$ZtZ_inv %*% crossprod(scores) %*% prep$ZtZ_inv * G / (G - 1)
  tibble(design = prep$design, estimate = prep$sign * beta[prep$k], se = sqrt(V[prep$k, prep$k]),
         resid_sd = sqrt(sum(resid^2) / max(prep$df_resid, 1)))
}
# the outcome of a panel in a prepared design's row order
align_outcome <- function(prep, panel, column = "y")
  panel[[column]][match(prep$row_key, paste(panel$id, panel$phase, panel$t))]

# ---- The planted truth ----------------------------------------------------------
# A semisynthetic outcome on the site's own patients, phases and measurement times:
# only the marker values are replaced. Slopes are in log units per day; A is the
# confounding demographic index (standardised age by default). The trajectory is
# continuous and piecewise linear, with these slopes:
#   every phase   illness + sev * severity + d_all * A * (arm_scale in the control)
#   after t0      + acute + d_acute * A                        (both arms)
#   ventilated    + d_vent * A + beta * (lambda + log VT/PBW)  (on the ventilator only)
#   after t_off   the after-t0 terms stay; the ventilated terms stop
# beta is the truth every design is scored against: the lambda gradient of the slope
# while ventilated, at a fixed VT/PBW.
plant_outcome <- function(panel, patients, truth, noise, rng_seed) {
  set.seed(rng_seed)
  A <- switch(truth$index, age = patients$age_z,
              sex = as.numeric(patients$sex_category == "Female") - mean(patients$sex_category == "Female"),
              stop("unknown confounding index ", truth$index))
  n <- nrow(patients); vent <- patients$vent
  u_level <- rnorm(n, 0, noise$sd_level)
  s_all  <- noise$illness + noise$sev * patients$sofa_z + rnorm(n, 0, noise$sd_slope) +
    truth$d_all * A * ifelse(vent == 1, 1, truth$arm_scale)
  s_t0   <- noise$acute + truth$d_acute * A
  s_vent <- vent * (truth$d_vent * A + truth$beta * (patients$ldisc + patients$lv0))
  i <- match(panel$id, patients$id)
  t_on  <- pmax(panel$t, 0)                                           # time since t0
  t_vnt <- pmin(t_on, dplyr::coalesce(patients$t_off_d[i], Inf))      # time on the ventilator
  panel$y_sim <- u_level[i] + s_all[i] * panel$t + s_t0[i] * t_on + s_vent[i] * t_vnt +
    rnorm(nrow(panel), 0, noise$sd_resid)
  panel
}
# The slope per unit of the confounding index that gives the unadjusted model the same
# lambda gradient as a strain effect of size beta: beta * var(lambda) / cov(lambda, A),
# in the ventilated patients.
confounding_equal_to <- function(beta, patients, index) {
  v <- patients[patients$vent == 1, ]
  A <- switch(index, age = v$age_z, sex = as.numeric(v$sex_category == "Female") - mean(v$sex_category == "Female"))
  beta * var(v$ldisc) / cov(v$ldisc, A)
}

# Noise for the planted outcome, read off the site's observed marker: the residual SD
# around each patient-phase line and the between-patient SDs of level and slope.
observed_noise <- function(panel) {
  per <- panel %>% group_by(id, phase) %>%
    summarise(n = n(), sxx = sum((t - mean(t))^2), level = mean(y),
              slope = sum((t - mean(t)) * (y - mean(y))) / sxx,
              rss = sum((y - mean(y) - slope * (t - mean(t)))^2), .groups = "drop")
  sd_resid <- sqrt(sum(per$rss) / max(sum(per$n - 2), 1))
  on <- per %>% filter(phase == "on")
  # slope variance beyond what measurement noise alone would give a fitted line
  sd_slope <- sqrt(max(weighted.mean((on$slope - weighted.mean(on$slope, on$sxx))^2 - sd_resid^2 / on$sxx, on$sxx), 0))
  list(sd_resid = sd_resid, sd_slope = sd_slope, sd_level = sd(on$level),
       illness = weighted.mean(on$slope, on$sxx), sev = 0, acute = 0)
}
