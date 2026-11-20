# =============================================================================
# causal_generator: the planted-truth data generator, its calibration to a site,
# and the realism summaries (the decision rule is in causal_common.R). Sourced by
# causal_testbed.R and causal_tests.R after causal_common.R; never run; writes nothing.
# =============================================================================
#
# simulate_site() keeps a site's patients (their demographics, severity, cohort and
# hospital time before t0) and generates everything else: the marker's trajectory,
# when it is measured, and when a patient is liberated, escalated or dies. The
# protocol (CAUSAL_TESTBED_PROTOCOL.md) says why each piece is there; in short, a
# generator that is the analysis model run backwards can only show that the code
# works. The pieces, each switched by one scenario parameter:
#
#   Truth              beta, the strain effect; the confounding index and where it acts
#                      (d_all every phase, d_acute after t0 in both cohorts, d_vent on
#                      the ventilator only; arm_scale the control's share of d_all)
#   Trajectory shape   curvature: an acute fall after t0 that recovers, with a
#                      patient-specific bend. time_shape: how a ventilation effect
#                      accrues with ventilator time (linear, saturating, delayed).
#                      dose_shape: linear in log strain, or a threshold
#   Errors             ar_share, ar_rho: a share of the residual variance is a
#                      stationary autocorrelated process (a patient's own drifting
#                      severity), the rest is measurement error with error_df degrees of
#                      freedom (Inf is Gaussian). u_slope: an unmeasured severity moves
#                      the slope
#   Observation        info_measure: the marker is drawn more often when it is worse.
#                      info_off, info_death: liberation is slower, and escalation and
#                      death faster, when it is worse
#   Selection          het_pre: the strain effect differs in the patients who have
#                      hospital time before intubation (the only ones a within-person
#                      design can use). pre_exposure: a share of patients are under
#                      positive pressure before the recorded intubation
#   VT/PBW             vt_model "indication": the tidal volume follows predicted size
#                      and the unmeasured severity; gate: only VT/PBW 6-8 enters
#   Exposure error     sd_lung: true lung size varies about GLI's prediction.
#                      height_error_cm: charted height is the true height with error.
#                      denominator: strain follows GLI (gli), GLI at age 25 (age25) or
#                      GLI without the race term (norace), while the analysis uses GLI
#
# The truth a design is scored against is the same in every scenario: the
# least-squares gradient, on the analysis lambda at a fixed VT/PBW, of each
# ventilated patient's planted strain effect per day over the ON_DAYS window. With a
# linear effect on GLI's lambda it is beta; with a threshold, a saturating time course
# or another denominator it is the linear summary the designs can at best return.
# =============================================================================

GRID_STEP <- 0.25               # days: the resolution of measurement and event times

# ---- Scenario parameters ----------------------------------------------------------
# every parameter with its matched-model value: the generator reduced to the analysis
# model (straight lines, independent Gaussian errors, nothing informative)
SCENARIO_DEFAULTS <- list(
  beta = 0, index = "age", d_all = 0, d_acute = 0, d_vent = 0, arm_scale = 1,
  curvature = 0, time_shape = "linear", dose_shape = "linear",
  ar_share = 0, ar_rho = 0.7, error_df = Inf, u_slope = 0,
  info_measure = 0, info_off = 0, info_death = 0,
  # death by routes other than the marker: log hazard per log unit of strain dose in the
  # ventilated (a direct path from strain to death), and per SD of age in both cohorts
  death_direct = 0, death_age = 0,
  het_pre = 1, pre_exposure = 0,
  vt_model = "observed", vt_lambda = 0.5, vt_severity = 0.08, vt_sd = 0.12, gate = FALSE,
  sd_lung = 0, height_error_cm = 0, denominator = "gli")
# the realistic base: every piece of realism at a moderate value. The one-at-a-time
# stress scenarios move a single parameter from here.
REALISTIC_BASE <- list(curvature = 1, ar_share = 0.5, error_df = 5, u_slope = 0.5,
                       info_measure = 0.3, info_off = 0.3, info_death = 0.5, sd_lung = 0.13)
scenario <- function(..., base = c("realistic", "matched")) {
  base <- match.arg(base)
  sc <- SCENARIO_DEFAULTS
  if (base == "realistic") sc[names(REALISTIC_BASE)] <- REALISTIC_BASE
  over <- list(...)
  unknown <- setdiff(names(over), names(sc))
  if (length(unknown)) stop("unknown scenario parameter: ", paste(unknown, collapse = ", "))
  sc[names(over)] <- over
  sc$base <- base
  sc
}

# the demographic index a confounder acts through, standardised in the ventilated
confounding_index <- function(patients, index) {
  v <- patients$vent == 1
  z <- function(x) (x - mean(x[v])) / sd(x[v])
  female <- as.numeric(patients$sex_category == "Female"); black <- as.numeric(patients$race_category == "BLACK")
  switch(index,
         age       = z(patients$age10),
         age_hinge = z(pmax(patients$age10 - 6.5, 0)),                 # nothing before 65
         sex       = female - mean(female[v]),
         # several inputs at once, in proportions that are not GLI's
         composite = z(0.6 * z(patients$age10) + 0.8 * female - 0.6 * black + 0.5 * z(patients$height_cm)),
         # the same demographics in exactly GLI's proportions: lambda itself
         gli       = z(patients$ldisc),
         stop("unknown confounding index ", index))
}
# The slope per unit of the index that gives the unadjusted model the same lambda
# gradient as a strain effect of size beta: beta * var(lambda) / cov(lambda, index)
confounding_equal_to <- function(beta, patients, index) {
  v <- patients$vent == 1
  beta * var(patients$ldisc[v]) / cov(patients$ldisc[v], confounding_index(patients, index)[v])
}

# ---- Calibration to the site ------------------------------------------------------
# Noise read off the site's observed marker, by cohort and phase: the residual SD about
# each patient-phase line and the mean slope; and, from the on phase, the
# between-patient SDs of level and slope.
PHASES <- c("pre", "on", "post"); ARMS <- c("ventilated", "no support")
# One least-squares line per patient-phase, computed for all of them at once: the
# count, span, level, slope, residual sum of squares, and the lag-1 products of the
# residuals (for phases with four or more measurements).
phase_lines <- function(panel) {
  o <- order(panel$id, panel$phase, panel$t); panel <- panel[o, , drop = FALSE]
  key <- paste(panel$id, panel$phase); first <- !duplicated(key)
  g <- cumsum(first); n <- tabulate(g)
  mean_g <- function(x) (rowsum(x, g, reorder = FALSE) / n)[g]
  tc <- panel$t - mean_g(panel$t); yc <- panel$y - mean_g(panel$y)
  sxx <- as.numeric(rowsum(tc^2, g, reorder = FALSE))
  slope <- as.numeric(rowsum(tc * yc, g, reorder = FALSE)) / sxx
  e <- yc - slope[g] * tc
  pair <- c(FALSE, g[-1] == g[-length(g)]) & (n >= 4)[g]            # a residual and the one before it
  lag <- function(x) { out <- numeric(length(n)); s <- rowsum(x[pair], g[pair], reorder = FALSE); out[as.integer(rownames(s))] <- s; out }
  data.frame(id = panel$id[first], phase = panel$phase[first], n = n,
             span = panel$t[cumsum(n)] - panel$t[cumsum(n) - n + 1],
             sxx = sxx, level = as.numeric(rowsum(panel$y, g, reorder = FALSE)) / n, slope = slope,
             rss = as.numeric(rowsum(e^2, g, reorder = FALSE)),
             lag_xy = lag(c(0, e[-1] * e[-length(e)])), lag_xx = ifelse(n >= 4, as.numeric(rowsum(e^2, g, reorder = FALSE)), 0),
             stringsAsFactors = FALSE)
}
observed_noise <- function(panel, patients) {
  per <- phase_lines(panel)
  per$arm <- patients$arm[match(per$id, patients$id)]
  cell <- function(f) {
    m <- matrix(NA_real_, 2, 3, dimnames = list(ARMS, PHASES))
    for (a in ARMS) for (p in PHASES) { d <- per[per$arm == a & per$phase == p, ]; if (nrow(d)) m[a, p] <- f(d) }
    m
  }
  sd_resid <- cell(function(d) sqrt(sum(d$rss) / max(sum(d$n - 2), 1)))
  drift    <- cell(function(d) weighted.mean(d$slope, d$sxx))
  pooled   <- sqrt(sum(per$rss) / max(sum(per$n - 2), 1))
  sd_resid[is.na(sd_resid) | sd_resid == 0] <- pooled; drift[is.na(drift)] <- 0
  on <- per[per$phase == "on", ]
  resid_on <- sd_resid[cbind(on$arm, "on")]; drift_on <- drift[cbind(on$arm, "on")]
  # slope variance beyond what measurement noise alone gives a fitted line
  sd_slope <- sqrt(max(weighted.mean((on$slope - drift_on)^2 - resid_on^2 / on$sxx, on$sxx), 0))
  list(sd_resid = sd_resid, drift = drift, sd_resid_pooled = pooled, sd_slope = sd_slope,
       sd_level = sd(on$level), mean_level = mean(on$level))
}

# Person-time at risk in each phase, and how each on phase ended, from the site's
# own patients. Used for the measurement rate and the event hazards.
phase_exposure <- function(patients) {
  death <- ifelse(is.na(patients$death_d), Inf, pmax(patients$death_d, 0))
  off   <- ifelse(is.na(patients$t_off_d), Inf, pmax(patients$t_off_d, 0))
  stay  <- ifelse(is.na(patients$stay_end_d), Inf, pmax(patients$stay_end_d, 0))
  # a death within an hour of the last ventilator record is a death on the ventilator
  died_on <- death <= off + 1 / 24
  on_end  <- pmin(off, death, stay, ON_DAYS)
  ended   <- ifelse(on_end >= ON_DAYS, "window", ifelse(died_on & death <= on_end, "death",
                                                       ifelse(off <= on_end, "off", "discharge")))
  post <- ifelse(patients$vent == 1 & ended == "off", pmax(pmin(POST_DAYS, stay - off, death - off), 0), 0)
  data.frame(pre = patients$pre_avail_d, on = on_end, post = post, ended = ended)
}

# calibrate_site(): everything the generator takes from the site's observed marker.
#   panel_raw   every measurement in the window (load_causal_panel(complete = FALSE))
#   panel       the analysis panel (complete patient-phases)
calibrate_site <- function(patients, panel_raw, panel) {
  expo <- phase_exposure(patients)
  arm_f <- factor(patients$arm, ARMS)
  hazard <- function(event) tapply(expo$ended == event, arm_f, sum) / tapply(expo$on, arm_f, sum)
  shape <- phase_shape(panel, patients)
  # starting values for the two-part measurement model, tuned by tune_generator():
  # the share of patients followed in a phase, and their draws per day
  monitor <- pmin(pmax(shape$share * 1.2, 0.01), 1)
  time <- sapply(PHASES, function(p) tapply(expo[[p]], arm_f, mean))
  rate <- shape$per_phase / pmax(time, 0.25)
  rate[!is.finite(rate)] <- 0
  list(noise = observed_noise(panel, patients), shape = shape, monitor = monitor, rate = rate,
       off_hazard = hazard("off"), death_hazard = hazard("death"),
       # corrections tune_generator() sets: multipliers for the residual SD and the slope
       # SD, and an added drift
       scale = list(resid = rate * 0 + 1, slope = 1, drift = rate * 0))
}
# what the analysis sees of a cohort's sampling: the share of its patients with a
# complete phase, and the measurements in one
phase_shape <- function(panel, patients) {
  arm <- patients$arm[match(panel$id, patients$id)]
  per <- unique(data.frame(id = panel$id, arm = arm, phase = panel$phase))
  n_phase <- table(factor(per$arm, ARMS), factor(per$phase, PHASES))
  n_meas  <- table(factor(arm, ARMS), factor(panel$phase, PHASES))
  as_matrix <- function(m) matrix(as.numeric(m), 2, 3, dimnames = list(ARMS, PHASES))
  list(share = as_matrix(n_phase) / as.numeric(table(factor(patients$arm, ARMS))),
       per_phase = as_matrix(ifelse(n_phase > 0, n_meas / n_phase, 0)))
}

# ---- The generator ----------------------------------------------------------------
# simulate_site(patients, calib, sc, sign, seed) -> list(patients, panel, truth)
#   patients  the analysis patients: the site's, with the simulated t_off_d and (when
#             the scenario simulates the tidal volume) lv0 and the gate applied
#   panel     id, phase, t, y for complete patient-phases
#   truth     the lambda gradient the designs are scored against
# sign is the marker's direction of injury (-1 platelets, +1 creatinine).
simulate_site <- function(patients, calib, sc, sign, seed) {
  set.seed(seed)
  n <- nrow(patients); vent <- patients$vent == 1
  noise <- calib$noise; sd_slope <- noise$sd_slope
  sd_resid <- noise$sd_resid * calib$scale$resid    # cohort x phase
  drift <- noise$drift + calib$scale$drift          # mean slope, cohort x phase
  female <- patients$sex_category == "Female"

  # --- exposure: the lambda strain follows, and the lambda the analysis sees ---
  ld_strain <- switch(sc$denominator, gli = patients$ldisc, age25 = patients$ldisc25, norace = patients$ldisc_norace,
                      stop("unknown denominator ", sc$denominator))
  if (sc$height_error_cm > 0) {
    if (sc$denominator != "gli") stop("height error is implemented for the GLI denominator only")
    # classical error: charted = true + error. Given the charted height, the true one
    # is shrunk towards its sex's mean by the reliability, plus what the error left.
    true_height <- patients$height_cm
    for (s in unique(female)) {
      r <- female == s; m <- mean(patients$height_cm[r])
      reliability <- max(1 - sc$height_error_cm^2 / var(patients$height_cm[r]), 0.05)
      true_height[r] <- m + reliability * (patients$height_cm[r] - m) + rnorm(sum(r), 0, sqrt(reliability) * sc$height_error_cm)
    }
    log_ratio <- function(h) {
      pbw <- ifelse(female, 45.5, 50) + 2.3 * (h / 2.54 - 60)
      race <- ifelse(patients$race_category == "WHITE", 1L, ifelse(patients$race_category == "BLACK", 2L, 5L))
      log(pbw / rspiro::pred_GLI(age = patients$age10 * 10, height = h / 100, gender = ifelse(female, 2L, 1L),
                                 ethnicity = race, param = "FVC"))
    }
    ld_strain <- patients$ldisc + log_ratio(true_height) - log_ratio(patients$height_cm)
  }
  lung_dev <- rnorm(n, 0, sc$sd_lung)             # true lung size about its prediction
  severity <- rnorm(n)                            # unmeasured, independent of the measured scores

  # --- tidal volume: the site's own, or set by predicted size and severity ---
  keep <- rep(TRUE, n); lv <- patients$lv0
  if (sc$vt_model == "indication") {
    log_vtpbw <- log(7) + sc$vt_lambda * patients$ldisc - sc$vt_severity * severity + rnorm(n, 0, sc$vt_sd)
    if (sc$gate) keep <- !vent | (log_vtpbw >= log(6) & log_vtpbw <= log(8))
    lv <- ifelse(vent, log_vtpbw - mean(log_vtpbw[vent & keep]), 0)
  } else if (sc$vt_model != "observed") stop("unknown vt_model ", sc$vt_model)
  strain <- ld_strain + lv - lung_dev             # log VT over true lung size, centred
  dose <- switch(sc$dose_shape, linear = strain,
                 threshold = 2 * pmax(strain - median(strain[vent & keep]), 0),   # nothing below the median
                 stop("unknown dose_shape ", sc$dose_shape))
  # the strain effect in the patients with hospital time before intubation, relative
  # to the rest, the cohort's mean effect held at beta
  has_pre <- patients$pre_avail_d >= MIN_SPAN_D
  share <- mean(has_pre[vent & keep])
  het <- ifelse(has_pre, sc$het_pre, (1 - share * sc$het_pre) / (1 - share))
  beta_i <- sc$beta * het

  # --- slopes, log units per day ---
  A <- confounding_index(patients, sc$index)
  level   <- noise$mean_level + rnorm(n, 0, noise$sd_level)
  s_all   <- rnorm(n, 0, sd_slope * calib$scale$slope) + sign * sc$u_slope * sd_slope * severity +
    sc$d_all * A * ifelse(vent, 1, sc$arm_scale)
  s_t0    <- sc$curvature * sign * sd_slope + sc$d_acute * A
  bend    <- sc$curvature * (-sign * sd_slope / ON_DAYS + rnorm(n, 0, 0.5 * sd_slope / ON_DAYS))
  s_vent  <- ifelse(vent, sc$d_vent * A + beta_i * dose, 0)
  accrue  <- switch(sc$time_shape, linear = function(v) v,
                    saturating = function(v) 3 * (1 - exp(-v / 3)),
                    delayed = function(v) pmax(v - 1, 0),
                    stop("unknown time_shape ", sc$time_shape))
  # positive pressure before the recorded intubation, in a share of the ventilated
  lead <- ifelse(vent & runif(n) < sc$pre_exposure, runif(n, 0.25, 1.5), 0)
  t_off <- rep(Inf, n)                            # ventilator time stops here
  vent_time <- function(t) pmax(pmin(t, t_off) + lead, 0) * (t > -lead)
  arm_i <- ifelse(vent, 1L, 2L)                   # the cohort's row of the calibration tables
  # the cohort's mean path: its observed mean slope in each phase, joined at t0 and t_off
  mean_path <- function(t) drift[arm_i, "pre"] * pmin(t, 0) + drift[arm_i, "on"] * pmin(pmax(t, 0), t_off) +
    drift[arm_i, "post"] * pmax(t - t_off, 0)
  # a patient's departure from it: their own slope, the acute bend, the confounder, strain
  departure <- function(t) {
    t_on <- pmax(t, 0)
    s_all * t + s_t0 * t_on + bend * t_on^2 + s_vent * accrue(vent_time(t))
  }

  # --- errors: an autocorrelated process and measurement error ---
  phi <- sc$ar_rho^GRID_STEP
  x <- rnorm(n)                                   # unit variance; scaled by the phase's SD
  meas_error <- function(k) if (is.finite(sc$error_df))
    rt(k, sc$error_df) / sqrt(sc$error_df / (sc$error_df - 2)) else rnorm(k)
  sd_dev <- sqrt(noise$sd_resid_pooled^2 + (3 * sd_slope)^2)   # the scale of "worse than expected"
  # measurement is two-part: whether the marker is being followed in a phase at all
  # (many patients have none before intubation), and then a schedule: labs come at
  # roughly regular intervals, not at random, and sooner when the marker is worse
  monitored <- matrix(runif(n * 3) < calib$monitor[cbind(rep(arm_i, 3), rep(1:3, each = n))], n, 3)
  next_draw <- rep(NA_real_, n); last_phase <- rep(0L, n)

  # --- the clock: measurements and events, one grid step at a time ---
  rate <- calib$rate
  ended <- rep(FALSE, n)                          # dead, or a control escalated
  t_death <- rep(NA_real_, n)
  extubated <- rep(FALSE, n)
  grid <- seq(-PRE_DAYS, ON_DAYS + POST_DAYS - GRID_STEP, by = GRID_STEP)
  out_id <- out_phase <- out_t <- out_y <- vector("list", length(grid))
  for (j in seq_along(grid)) {
    t <- grid[j]
    phase <- rep(NA_integer_, n)
    if (t < 0) phase[t >= -patients$pre_avail_d] <- 1L else {
      phase[!ended & !extubated & t < ON_DAYS] <- 2L
      phase[vent & extubated & t < t_off + POST_DAYS] <- 3L
    }
    at_risk <- which(!is.na(phase))
    cell <- cbind(arm_i[at_risk], phase[at_risk])
    sd_now <- rep(NA_real_, n); sd_now[at_risk] <- sd_resid[cell]
    ar <- sqrt(sc$ar_share) * sd_now * x
    # how much worse than expected the marker is now, in the direction of injury
    worse <- sign * (departure(t) + ifelse(is.na(ar), 0, ar)) / sd_dev
    gap <- rep(NA_real_, n)
    gap[at_risk] <- exp(-sc$info_measure * worse[at_risk]) / pmax(rate[cell], 1e-6)
    entered <- at_risk[phase[at_risk] != last_phase[at_risk]]         # a new phase: its first draw comes early
    next_draw[entered] <- t + runif(length(entered)) * 0.5 * gap[entered]
    last_phase[at_risk] <- phase[at_risk]
    drawn <- at_risk[monitored[cbind(at_risk, phase[at_risk])] & next_draw[at_risk] < t + GRID_STEP]
    if (length(drawn)) {
      tm <- rep(t, n); tm[drawn] <- pmax(next_draw[drawn], t)
      next_draw[drawn] <- tm[drawn] + gap[drawn] * runif(length(drawn), 0.6, 1.4)
      out_id[[j]] <- drawn; out_phase[[j]] <- phase[drawn]; out_t[[j]] <- tm[drawn]
      out_y[[j]] <- (level + mean_path(tm) + departure(tm) + ar)[drawn] +
        sqrt(1 - sc$ar_share) * sd_now[drawn] * meas_error(length(drawn))
    }
    if (t >= 0) {                                 # events take effect at the end of the step
      on <- which(phase == 2L)
      p <- function(h, a) 1 - exp(-h * GRID_STEP * exp(a * worse[on]))
      other_routes <- exp(sc$death_direct * strain[on] * vent[on] + sc$death_age * patients$age_z[on])
      dies <- runif(length(on)) < p(calib$death_hazard[arm_i[on]] * other_routes, sc$info_death)
      # liberation is slower when the marker is worse; a control's escalation is faster
      off  <- runif(length(on)) < p(calib$off_hazard[arm_i[on]], ifelse(vent[on], -sc$info_off, sc$info_off))
      ended[on[dies]] <- TRUE; t_death[on[dies]] <- t + GRID_STEP
      left <- on[off & !dies]
      t_off[left] <- t + GRID_STEP
      extubated[left[vent[left]]] <- TRUE
      ended[left[!vent[left]]] <- TRUE
    }
    x <- phi * x + sqrt(1 - phi^2) * rnorm(n)
  }
  panel <- data.frame(id = patients$id[unlist(out_id)], phase = c("pre", "on", "post")[unlist(out_phase)],
                      t = unlist(out_t), y = unlist(out_y), stringsAsFactors = FALSE)
  panel <- panel[patients$id[unlist(out_id)] %in% patients$id[keep], , drop = FALSE]
  panel <- complete_phases(panel)

  sim <- patients[keep, , drop = FALSE]
  sim$lv0 <- lv[keep]; sim$t_off_d <- ifelse(is.finite(t_off[keep]), t_off[keep], NA_real_)
  sim$death_d <- t_death[keep]; sim$level0 <- level[keep]      # the simulated death, and the marker at t0
  # the truth: the planted strain effect per day over the window, as a gradient on
  # the lambda the analysis uses, at a fixed VT/PBW
  v <- vent & keep
  per_day <- beta_i[v] * dose[v] * accrue(ON_DAYS) / ON_DAYS
  truth <- if (sc$beta == 0) 0 else unname(coef(lm(per_day ~ patients$ldisc[v] + lv[v]))[2])
  # the same gradient with VT/PBW not held (what the designs without that term estimate),
  # and of the total effect accrued by each cut of the time course
  truth_dropped <- if (sc$beta == 0) 0 else unname(coef(lm(per_day ~ patients$ldisc[v]))[2])
  cuts <- c(CAUSAL_KNOTS[CAUSAL_KNOTS > 0 & CAUSAL_KNOTS < ON_DAYS], ON_DAYS)
  truth_by_day <- setNames(vapply(cuts, function(b) if (sc$beta == 0) 0 else
    unname(coef(lm(I(beta_i[v] * dose[v] * accrue(b)) ~ patients$ldisc[v] + lv[v]))[2]), numeric(1)), cuts)
  list(patients = sim, panel = panel, truth = truth, truth_dropped = truth_dropped, truth_by_day = truth_by_day,
       ended = c(off = mean(is.finite(t_off[v])), dead = mean(ended[v])))
}

# tune_generator(): the generator's own noise and measurement rate, set so that what
# the ANALYSIS would read off a simulated site matches what it reads off the real one.
# The observed residual SD and slope SD are estimates made from short, sparse lines; a
# generator that plugs them in as if they were the truth produces data that are too
# noisy or too smooth once curvature and autocorrelation are added. A few passes of a
# correction under the scenario's own base, at no effect.
tune_generator <- function(patients, calib, base, sign, seed, passes = 8) {
  target <- calib$noise; shape <- calib$shape
  n_arm <- as.numeric(table(factor(patients$arm, ARMS)))
  bounded <- function(x, lo, hi) pmin(pmax(x, lo), hi)
  for (pass in seq_len(passes)) {
    sim <- simulate_site(patients, calib, scenario(base = base), sign, seed + pass)
    got <- observed_noise(sim$panel, sim$patients); got_shape <- phase_shape(sim$panel, sim$patients)
    wanted <- shape$share > 0
    # a cell is tuned on its noise only when the simulated site has enough of it to read
    enough <- wanted & got_shape$share * n_arm >= 30
    calib$scale$resid[enough] <- bounded(calib$scale$resid[enough] * target$sd_resid[enough] / got$sd_resid[enough], 0.2, 5)
    calib$scale$drift[enough] <- calib$scale$drift[enough] + bounded(target$drift[enough] - got$drift[enough], -0.05, 0.05)
    if (got$sd_slope > 0 && target$sd_slope > 0)
      calib$scale$slope <- bounded(calib$scale$slope * target$sd_slope / got$sd_slope, 0.05, 5)
    # sampling: half-steps, so that the share followed and the rate settle together
    share_ratio <- ifelse(got_shape$share > 0, shape$share / got_shape$share, 2)
    count_ratio <- ifelse(got_shape$share > 0, shape$per_phase / got_shape$per_phase, 1.5)
    calib$monitor[wanted] <- bounded(calib$monitor[wanted] * bounded(share_ratio[wanted], 0.5, 2)^0.5, 0.005, 1)
    calib$rate[wanted]    <- calib$rate[wanted] * bounded(count_ratio[wanted], 0.5, 2)^0.5
    calib$monitor[!wanted] <- 0
  }
  calib
}

# ---- Realism: the summaries a simulated site is compared with the real one on -----
panel_summaries <- function(panel, patients) {
  per <- phase_lines(panel)
  per$arm <- patients$arm[match(per$id, patients$id)]
  n_arm <- table(patients$arm)
  do.call(rbind, lapply(split(per, list(per$arm, per$phase), drop = TRUE), function(d) {
    data.frame(arm = d$arm[1], phase = d$phase[1],
               statistic = c("share of patients with the phase", "measurements per patient-phase", "span, days (median)",
                             "slope per day (10th centile)", "slope per day (median)", "slope per day (90th centile)",
                             "residual SD", "lag-1 residual autocorrelation"),
               value = c(nrow(d) / n_arm[[d$arm[1]]], mean(d$n), median(d$span),
                         quantile(d$slope, c(0.1, 0.5, 0.9), names = FALSE),
                         sqrt(sum(d$rss) / max(sum(d$n - 2), 1)),
                         if (sum(d$lag_xx) > 0) sum(d$lag_xy) / sum(d$lag_xx) else NA_real_))
  }))
}
