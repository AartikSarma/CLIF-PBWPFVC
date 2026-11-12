# =============================================================================
# Unit tests for the causal supplements (causal_common.R, causal_generator.R)
# =============================================================================
#
# Checks, on a small fabricated cohort (no CLIF data, no site configuration):
#   1. the within-patient-phase estimator equals least squares with a dummy for
#      every patient-phase, coefficient and cluster-robust standard error
#   2. every design returns the planted effect exactly when there is no noise
#   3. a planted age confounder is removed exactly by the designs that free age and,
#      when it acts in every phase, nearly by the control and within-person designs;
#      not by the unadjusted model
#   4. the generator keeps its own rules (phase windows, complete phases, the gate,
#      reproducibility by seed) and reports the planted truth
#   5. the decision rule returns what its definition says, one way and both ways
#   6. the Poisson model for an outcome with zeros equals glm(poisson) with a dummy
#      for every patient-phase, coefficient and cluster-robust standard error
#   7. the time course: every segment returns the planted effect without noise, the
#      cumulative effect is the planted effect times the days, and its standard error
#      from added influences equals the direct one where the two can be compared
#   8. the VT/PBW gradient equals the lambda gradient when strain is planted
#   9. injury and death as states: the thresholds, who is left out, which transition
#      each patient's events are counted in, and the gradients of a simulated site
# Stops at the first failure; prints "all tests passed" otherwise.
#
# Usage: uvr run code/supplement/causal_tests.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
rm(list = ls())
# 20_biotrauma_grid.R needs a site configuration; the tests need only this helper of it
sorted_term <- function(term) vapply(strsplit(term, ":"), function(p) paste(sort(p), collapse = ":"), character(1))
source(file.path("code", "supplement", "causal_common.R"))
source(file.path("code", "supplement", "causal_generator.R"))

check <- function(ok, what) { if (!isTRUE(ok)) stop("FAILED: ", what, call. = FALSE); message("  ok  ", what) }

# ---- a fabricated cohort -----------------------------------------------------------
set.seed(1)
n <- 800
fab <- data.frame(id = sprintf("p%04d", seq_len(n)), stringsAsFactors = FALSE)
fab$arm <- rep(c("ventilated", "no support"), each = n / 2); fab$vent <- as.numeric(fab$arm == "ventilated")
fab$age10 <- runif(n, 2, 9); fab$height_cm <- rnorm(n, 170, 9)
fab$sex_category <- factor(sample(c("Male", "Female"), n, TRUE), levels = c("Male", "Female"))
fab$race_category <- factor(sample(c("WHITE", "BLACK", "OTHER"), n, TRUE, c(0.6, 0.25, 0.15)), levels = c("WHITE", "BLACK", "OTHER"))
# a height piece that differs by sex (the fingerprint), and a lambda that is that piece
# plus terms linear in age and shifts by sex and race: the structure of GLI's ratio
fab$fp <- ifelse(fab$sex_category == "Female", -1, 1) * 0.004 * (fab$height_cm - 170) - 0.00015 * (fab$height_cm - 170)^2
fab$ldisc <- fab$fp + 0.045 * fab$age10 + 0.05 * (fab$sex_category == "Female") + 0.12 * (fab$race_category == "BLACK")
fab$ldisc <- fab$ldisc - mean(fab$ldisc[fab$vent == 1])
fab$ldisc25 <- fab$ldisc - 0.045 * (fab$age10 - mean(fab$age10)); fab$ldisc_norace <- fab$ldisc - 0.12 * (fab$race_category == "BLACK")
fab$pbw <- 60; fab$age_z <- as.numeric(scale(fab$age10))
fab$lv0 <- ifelse(fab$vent == 1, rnorm(n, 0, 0.08), 0); fab$lv0[fab$vent == 1] <- fab$lv0[fab$vent == 1] - mean(fab$lv0[fab$vent == 1])
fab$sofa_z <- rnorm(n); fab$sf_z <- rnorm(n)
fab$pre_avail_d <- ifelse(runif(n) < 0.7, runif(n, 1, PRE_DAYS), 0)
fab$t_off_d <- NA_real_; fab$stay_end_d <- 30; fab$death_d <- NA_real_
ones <- matrix(1, 2, 3, dimnames = list(ARMS, PHASES))
fab_calib <- function(noise_sd) list(
  noise = list(sd_resid = ones * noise_sd, drift = ones * 0, sd_resid_pooled = max(noise_sd, 1e-6),
               sd_slope = noise_sd / 4, sd_level = 0.5, mean_level = 5),
  monitor = ones, rate = ones * 1.5, off_hazard = c(0.15, 0.05), death_hazard = c(0.01, 0.005),
  scale = list(resid = ones, slope = 1, drift = ones * 0))
designs <- causal_designs(vtpbw = TRUE)

# ---- 1. the estimator against dummy-variable least squares ----------------------------
message("1. The within estimator")
sim <- simulate_site(fab, fab_calib(0.1), scenario(beta = -0.15, base = "matched"), sign = -1, seed = 11)
small_ids <- c(head(fab$id[fab$vent == 1], 150))
small <- sim$panel[sim$panel$id %in% small_ids & sim$panel$phase == "on", ]
prep <- prepare_design(designs[designs$design == "unadjusted", ], sim$patients, small)
fit  <- fit_design(prep, align_outcome(prep, small))
d <- cbind(small, sim$patients[match(small$id, sim$patients$id), c("ldisc", "sofa_z", "sf_z", "lv0")])
dummy <- lm(y ~ 0 + factor(id) + t + t:ldisc + t:sofa_z + t:sf_z + t:lv0, data = d)
check(abs(fit$estimate - coef(dummy)[["t:ldisc"]]) < 1e-9, "coefficient equals the dummy-variable regression")
X <- model.matrix(dummy); e <- resid(dummy); bread <- solve(crossprod(X))
scores <- rowsum(X * e, d$id); G <- nrow(scores)
se_dummy <- sqrt((bread %*% crossprod(scores) %*% bread * G / (G - 1))["t:ldisc", "t:ldisc"])
check(abs(fit$se - se_dummy) < 1e-9, "cluster-robust standard error equals the dummy-variable sandwich")
check(max(abs(rowsum(prep$Z, paste(small$id, small$phase)[prep$rows]))) < 1e-9, "centred columns sum to zero within a patient-phase")

per <- phase_lines(sim$panel)
agree <- vapply(sample(nrow(per), 25), function(k) {
  d <- sim$panel[sim$panel$id == per$id[k] & sim$panel$phase == per$phase[k], ]; f <- lm(y ~ t, d); e <- resid(f)[order(d$t)]
  abs(coef(f)[[2]] - per$slope[k]) < 1e-9 && abs(sum(e^2) - per$rss[k]) < 1e-9 && nrow(d) == per$n[k] &&
    abs(diff(range(d$t)) - per$span[k]) < 1e-12 && (nrow(d) < 4 || abs(sum(e[-1] * e[-length(e)]) - per$lag_xy[k]) < 1e-9)
}, logical(1))
check(all(agree), "the per-phase lines used for calibration equal lm, phase by phase")

# ---- 2. no noise: every design returns the planted effect ------------------------------
message("2. Exact recovery without noise")
exact <- simulate_site(fab, fab_calib(0), scenario(beta = -0.15, base = "matched"), sign = -1, seed = 12)
fits <- fit_all_designs(designs, exact$patients, exact$panel)
check(all(!is.na(fits$estimate)), "all fifteen designs are estimable")
check(max(abs(fits$estimate - (-0.15))) < 1e-6, paste0("every design returns -0.15 (largest error ", signif(max(abs(fits$estimate + 0.15)), 2), ")"))
check(abs(exact$truth - (-0.15)) < 1e-10, "the reported truth is the planted beta")

# ---- 3. a planted confounder ----------------------------------------------------------
message("3. A planted age confounder, no strain effect")
d_age <- confounding_equal_to(-0.15, fab, "age")
conf <- simulate_site(fab, fab_calib(0), scenario(d_vent = d_age, base = "matched"), sign = -1, seed = 13)
fits <- fit_all_designs(designs, conf$patients, conf$panel)
est <- setNames(fits$estimate, fits$design)
check(conf$truth == 0, "the truth is zero")
check(all(abs(est[c("free_age", "free_age_did", "adjusted")]) < 1e-6), "the designs that free age return zero")
check(abs(est[["unadjusted"]]) > 0.1 && abs(est[["control_did"]]) > 0.1, "the unadjusted model and the control contrast do not")
everywhere <- simulate_site(fab, fab_calib(0), scenario(d_all = d_age, base = "matched"), sign = -1, seed = 14)
est <- with(fit_all_designs(designs, everywhere$patients, everywhere$panel), setNames(estimate, design))
# It cancels in expectation, not identically: the two phases weight a patient by the
# spread of their measurement times, and the two cohorts are different patients.
check(all(abs(est[c("control_did", "within_onoff", "within_triple", "off_reversal")]) < 0.2 * abs(est[["unadjusted"]])),
      "a confounder in every phase and both cohorts nearly cancels in the control and within-person designs")

# ---- 4. the generator's own rules -----------------------------------------------------
message("4. Generator invariants")
g <- simulate_site(fab, fab_calib(0.1), scenario(beta = -0.15), sign = -1, seed = 21)
p <- g$panel; pat <- g$patients[match(p$id, g$patients$id), ]
check(all(p$t[p$phase == "pre"] < 0 & p$t[p$phase == "pre"] >= -pat$pre_avail_d[p$phase == "pre"]), "pre measurements lie in the patient's time before t0")
check(all(p$t[p$phase == "on"] >= 0 & p$t[p$phase == "on"] < pmin(ifelse(is.na(pat$t_off_d), Inf, pat$t_off_d), ON_DAYS)[p$phase == "on"]),
      "on measurements lie between t0 and liberation or the end of the window")
post <- p$phase == "post"
check(all(pat$vent[post] == 1) && all(p$t[post] >= pat$t_off_d[post] & p$t[post] < pat$t_off_d[post] + POST_DAYS), "post measurements follow a ventilated patient's liberation")
spans <- tapply(p$t, paste(p$id, p$phase), function(t) c(length(t), max(t) - min(t)))
check(all(vapply(spans, function(s) s[1] >= 2 && s[2] >= MIN_SPAN_D, logical(1))), "every patient-phase is complete")
again <- simulate_site(fab, fab_calib(0.1), scenario(beta = -0.15), sign = -1, seed = 21)
other <- simulate_site(fab, fab_calib(0.1), scenario(beta = -0.15), sign = -1, seed = 22)
check(identical(g$panel, again$panel), "the same seed gives the same site")
check(!identical(g$panel$y, other$panel$y), "another seed gives another")
gated <- simulate_site(fab, fab_calib(0.1), scenario(vt_model = "indication", gate = TRUE), sign = -1, seed = 23)
vt <- exp(gated$patients$lv0[gated$patients$vent == 1])
check(nrow(gated$patients) < n && max(vt) / min(vt) <= 8 / 6 + 1e-9 && sum(gated$patients$vent == 0) == n / 2,
      "the gate keeps ventilated patients within VT/PBW 6-8 and every control")
check(simulate_site(fab, fab_calib(0.1), scenario(), sign = -1, seed = 24)$truth == 0, "no planted effect, zero truth")
thr <- simulate_site(fab, fab_calib(0.1), scenario(beta = -0.15, dose_shape = "threshold"), sign = -1, seed = 25)$truth
check(thr < -0.075 && thr > -0.3, "a threshold effect has a linear summary of beta's sign and order")
check(inherits(try(scenario(not_a_parameter = 1), silent = TRUE), "try-error"), "an unknown scenario parameter is refused")

# ---- 5. the decision rule ---------------------------------------------------------------
message("5. The decision rule")
verdict <- function(z, est = NULL) {
  f <- tibble(design = designs$design, estimate = -0.15, se = 0.15 / 4)       # all significant, injury sign -1
  for (dn in names(z)) f$se[f$design == dn] <- abs(f$estimate[f$design == dn] / z[[dn]])
  for (dn in names(z)) if (z[[dn]] < 0) f$estimate[f$design == dn] <- 0.15
  for (dn in names(est)) f$estimate[f$design == dn] <- est[[dn]]
  triangulation_supports(f, sign = -1)[["rule"]]
}
check(verdict(list()), "every design significant in the direction of injury: supports")
check(!verdict(list(free_age_did = 1)), "the age rung not significant: does not support")
check(!verdict(list(control_did = 1)), "the control contrast not significant: does not support")
check(verdict(list(within_onoff = 1)), "one within-person design significant is enough")
check(!verdict(list(within_onoff = 1, off_reversal = 1)), "neither within-person design significant: does not support")
check(!verdict(list(within_onoff = -3)), "a within-person design significant the other way: does not support")
check(!verdict(list(), est = c(free_sex_did = -0.45)), "a rung three times the control contrast: does not support")
na <- tibble(design = designs$design, estimate = -0.15, se = 0.15 / 4); na$estimate[na$design == "within_onoff"] <- NA
check(triangulation_supports(na, -1)[["rule"]], "a design that is not estimable counts as not significant, no more")

good <- tibble(design = designs$design, estimate = -0.15, se = 0.15 / 4)
bad <- good; bad$se[bad$design == "control_did"] <- 0.15
check(triangulation_supports_both(good, good, -1)[["rule_both"]], "supported with VT/PBW held and dropped: supports")
check(!triangulation_supports_both(good, bad, -1)[["rule_both"]] && !triangulation_supports_both(bad, good, -1)[["rule_both"]],
      "supported one way only: does not support")

# ---- 6. the Poisson model -----------------------------------------------------------------
message("6. The Poisson model for an outcome with zeros")
set.seed(31)
pp <- sim$panel[sim$panel$id %in% small_ids & sim$panel$phase == "on", ]
pat <- sim$patients[match(pp$id, sim$patients$id), ]
pp$y <- rpois(nrow(pp), exp(-0.5 + 0.2 * pp$t * pat$ldisc * 5 + 0.1 * pp$t * pat$sofa_z)) * rgamma(nrow(pp), 2, 2)   # zeros, not Poisson
prep <- prepare_design(designs[designs$design == "unadjusted", ], sim$patients, pp, raw = TRUE)
fit  <- fit_design(prep, align_outcome(prep, pp), family = "poisson")
d <- cbind(pp, pat[c("ldisc", "sofa_z", "sf_z", "lv0")])
d <- d[d$id %in% names(which(tapply(d$y, d$id, sum) > 0)), ]
dummy <- suppressWarnings(glm(y ~ 0 + factor(id) + t + t:ldisc + t:sofa_z + t:sf_z + t:lv0, family = quasipoisson, data = d))
check(mean(pp$y == 0) > 0.2, "the fabricated outcome has zeros")
check(abs(fit$estimate - coef(dummy)[["t:ldisc"]]) < 1e-6, "coefficient equals the dummy-variable Poisson regression")
X <- model.matrix(dummy); mu <- fitted(dummy); bread <- solve(crossprod(X * sqrt(mu)))
scores <- rowsum(X * (d$y - mu), d$id); G <- nrow(scores)
se_dummy <- sqrt((bread %*% crossprod(scores) %*% bread * G / (G - 1))["t:ldisc", "t:ldisc"])
check(abs(fit$se / se_dummy - 1) < 1e-4, "cluster-robust standard error equals the dummy-variable sandwich")

# ---- 7. the time course ---------------------------------------------------------------------
message("7. The time course")
tc <- fit_time_course(designs[designs$design == "control_did", ], exact$patients, exact$panel)
seg <- tc[tc$quantity == "segment", ]; cum <- tc[tc$quantity == "cumulative", ]
check(max(abs(seg$estimate + 0.15)) < 1e-6, "every segment returns the planted effect without noise")
check(max(abs(cum$estimate + 0.15 * cum$to)) < 1e-5, "the cumulative effect is the planted effect times the days")
check(all(abs(exact$truth_by_day + 0.15 * as.numeric(names(exact$truth_by_day))) < 1e-9), "the reported cumulative truth agrees")
noisy <- fit_time_course(designs[designs$design == "unadjusted", ], sim$patients, sim$panel)
first <- noisy[noisy$quantity == "segment", ][1, ]; first_cum <- noisy[noisy$quantity == "cumulative", ][1, ]
check(abs(first_cum$se - first$se * first$to) < 1e-12, "over one segment the added influences give the segment's own standard error")
sat <- simulate_site(fab, fab_calib(0), scenario(beta = -0.15, time_shape = "saturating", base = "matched"), sign = -1, seed = 15)
tc <- fit_time_course(designs[designs$design == "unadjusted", ], sat$patients, sat$panel)
cum <- tc[tc$quantity == "cumulative", ]
check(all(abs(cum$estimate / sat$truth_by_day - 1) < 0.1), "a saturating effect: the cumulative estimate is within 10% of the accrued truth at every cut")
check(abs(tc$estimate[tc$quantity == "segment"][1]) > 1.5 * abs(tc$estimate[tc$quantity == "segment"][3]), "and the first segment is steeper than the last")

# ---- 8. the VT/PBW gradient -------------------------------------------------------------------
message("8. The gradient on VT/PBW")
f <- fit_design(prepare_design(designs[designs$design == "unadjusted", ], exact$patients, exact$panel),
                exact$panel$y[exact$panel$phase == "on" & exact$panel$id %in% exact$patients$id[exact$patients$vent == 1]])
check(abs(f$vt_gradient + 0.15) < 1e-6, "with strain planted on log VT - log size, the VT/PBW gradient equals the lambda gradient")

# ---- 9. injury and death as states ---------------------------------------------------------------
message("9. Injury and death as states")
six <- data.frame(id = paste0("q", 1:6), t0 = 0, t_off_d = c(NA, NA, NA, NA, 3, NA), stringsAsFactors = FALSE)
creat <- data.frame(id = c("q1", "q1", "q1", "q2", "q2", "q3", "q4", "q4", "q5", "q5", "q6", "q6"),
                    t  = c(-1, 0.1, 2,    -1, 0.2,   1,    0, 5,      0, 4,       0, 8),
                    value = c(1, 1.1, 2.3, 1, 2.5,   3,    1, 1.9,    1, 2.5,     1, 3))
inj <- load_injury_times(NULL, six, "creatinine", values = creat)
check(identical(inj$reason, c("eligible", "already injured at baseline", "no baseline value", "eligible", "eligible", "eligible")),
      "creatinine: eligible, already injured at baseline, and no baseline are told apart")
check(isTRUE(all.equal(inj$t_injury, c(2, NA, NA, NA, NA, NA))),
      "creatinine doubles from the lowest baseline value: injured then; 1.9x, after liberation, or after the window: not")
plt <- data.frame(id = c("q1", "q1", "q2", "q2", "q3", "q3", "q4", "q4"), t = c(0, 1, 0, 2, 0, 1, 0, 1),
                  value = c(200, 90, 400, 190, 80, 70, 200, 150))
inj_p <- load_injury_times(NULL, six[1:4, ], "platelets", values = plt)
check(isTRUE(all.equal(inj_p$t_injury, c(1, 2, NA, NA))) && inj_p$reason[3] == "already injured at baseline",
      "platelets: below 100, or half of baseline; a count already below 100 is left out")
ne <- data.frame(id = c("q1", "q2", "q2", "q3"), t = c(1.5, -0.5, 1, 0.2), value = c(0.1, 0.05, 0.1, 0))
inj_v <- load_injury_times(NULL, six[1:4, ], "pressor", values = ne)
check(isTRUE(all.equal(inj_v$t_injury, c(1.5, NA, NA, NA))) && inj_v$reason[2] == "already injured at baseline",
      "pressor: the first rate above zero after t0; one running the day before is left out")
st <- data.frame(id = paste0("r", 1:5), arm = "ventilated", vent = 1, t_off_d = NA_real_, death_d = c(NA, 5, 5, 40, 2),
                 eligible = TRUE, reason = "eligible", t_injury = c(NA, NA, 2, 3, 4), stringsAsFactors = FALSE)
st <- illness_death_states(st[setdiff(names(st), c("eligible", "reason", "t_injury"))], st[c("id", "eligible", "reason", "t_injury")], horizon = 28)
ev <- function(tr) { r <- transition_rows(st, tr); setNames(r$event, r$id) }
check(identical(unname(ev("injury")), c(0, 0, 1, 1, 0)), "an injury after death is not an injury")
check(identical(unname(ev("death without injury")), c(0, 1, 0, 0, 1)), "death without injury: the patients who died uninjured")
check(identical(ev("death after injury"), c(r3 = 1, r4 = 0)), "death after injury: the injured only, a death past the horizon not counted")
check(sum(ev("death, any")) == sum(ev("death without injury")) + sum(ev("death after injury")), "every death is counted in exactly one route")
big <- simulate_site(fab, fab_calib(0.1), scenario(beta = 1.5, base = "matched", info_death = 0.3), sign = 1, seed = 41)
vals <- rbind(data.frame(id = big$patients$id, t = 0, value = exp(big$patients$level0)),
              data.frame(id = big$panel$id, t = big$panel$t, value = exp(big$panel$y)))
big$patients$t0 <- 0
sts <- illness_death_states(big$patients, load_injury_times(NULL, big$patients, "creatinine", values = vals), horizon = ON_DAYS)
idm <- fit_illness_death(designs[designs$design == "unadjusted", ], sts)
check(idm$estimate[idm$transition == "injury"] / idm$se[idm$transition == "injury"] > 3,
      "a simulated site with a large strain effect: the gradient on injury is clearly positive")
one_fit <- survival::coxph(survival::Surv(start, stop, event) ~ ldisc + sofa_z + sf_z + lv0, data = transition_rows(sts[sts$vent == 1, ], "injury"), cluster = id)
check(abs(idm$estimate[idm$transition == "injury"] - coef(one_fit)[["ldisc"]]) < 1e-8 &&
        abs(idm$se[idm$transition == "injury"] - sqrt(one_fit$var[1, 1])) < 1e-8,
      "the injury gradient and its robust standard error equal coxph's")

# ---- 10. PBW and PFVC as separate terms ---------------------------------------------------------
message("10. PBW freed beside lambda, VT/PBW dropped")
# strain on log VT - log PFVC, with tidal volume set as log VT = rho log PBW: the slope is
# b (lambda + log VT/PBW) = b lambda + b (rho - 1) log PBW
sep_pat <- exact$patients; set.seed(31)
sep_pat$lpbw <- 0.5 * sep_pat$ldisc + rnorm(nrow(sep_pat), 0, 0.1); rho <- 0.6; b <- -0.15
sep_pat$lv0 <- ifelse(sep_pat$vent == 1, (rho - 1) * sep_pat$lpbw, 0)
sep_pan <- exact$panel; i <- match(sep_pan$id, sep_pat$id)
on_vent <- sep_pan$phase == "on" & sep_pat$vent[i] == 1
sep_pan$y <- 3 + 0.1 * (i %% 7) + ifelse(on_vent, b * (sep_pat$ldisc[i] + sep_pat$lv0[i]), 0) * sep_pan$t
dropped <- causal_designs(vtpbw = FALSE)
sep <- fit_pbw_separate(dropped[dropped$design %in% c("unadjusted", "control_did"), ], sep_pat, sep_pan)
check(max(abs(sep$estimate - b)) < 1e-6, "with log PBW freed, the lambda coefficient is the PFVC gradient whatever rho is")
check(max(abs(sep$pbw_extra - b * (rho - 1))) < 1e-6, "and the log PBW coefficient is b (rho - 1)")
sep_rho <- fit_pbw_separate(dropped[dropped$design == "unadjusted", ], sep_pat, sep_pan, rho = rho)
check(abs(sep_rho$pbw_extra_strain - sep_rho$pbw_extra) < 1e-6, "given rho, the strain prediction for the log PBW coefficient is met")
sep_pan2 <- sep_pan; sep_pan2$y <- sep_pan2$y + ifelse(on_vent, 0.05 * sep_pat$lpbw[i], 0) * sep_pan2$t + rnorm(nrow(sep_pan2), 0, 0.01)
off <- fit_pbw_separate(dropped[dropped$design == "unadjusted", ], sep_pat, sep_pan2, rho = rho)
check(off$p_strain_consistent < 0.001, "a direct PBW effect on top of strain is flagged as inconsistent with it")
ratio_only <- fit_all_designs(dropped[dropped$design == "unadjusted", ], sep_pat, sep_pan)
check(abs(ratio_only$estimate - b) > 0.01, "lambda alone misses it when rho is not 1 and PBW moves with lambda")
check(is.null(with_pbw_separate(dropped[dropped$design == "fingerprint", ])) &&
        is.null(with_pbw_separate(causal_designs(TRUE)[1, ])), "no separate-terms version of the fingerprint or of a design with VT/PBW")
r <- suppressWarnings(vt_follows_pbw(sep_pat))   # an exact fit, by construction
check(abs(r$rho[r$group == "all"] - rho) < 1e-8, "rho is recovered from log VT/PBW and log PBW")

# ---- 11. the age-mediation split -----------------------------------------------------------------
message("11. lambda split into its age part and the rest")
med_pat <- exact$patients; set.seed(41)
med_pat$lage <- 0.6 * (med_pat$age10 - mean(med_pat$age10)) / 10; med_pat$lrest <- med_pat$ldisc - med_pat$lage
med_pan <- exact$panel; i <- match(med_pan$id, med_pat$id); on_v <- med_pan$phase == "on" & med_pat$vent[i] == 1
make_y <- function(b_rest, b_age) 3 + 0.1 * (i %% 7) + ifelse(on_v, b_rest * med_pat$lrest[i] + b_age * med_pat$lage[i], 0) * med_pan$t
med_pan$y <- make_y(-0.15, -0.15)
full <- fit_age_mediation(designs[designs$design %in% c("unadjusted", "control_did"), ], med_pat, med_pan)
check(max(abs(full$share_mediated - 1)) < 1e-6, "age acting only through lambda: share mediated 1")
med_pan$y <- make_y(-0.15, -0.45)
part <- fit_age_mediation(designs[designs$design %in% c("unadjusted", "control_did"), ], med_pat, med_pan)
check(max(abs(part$share_mediated - 1 / 3)) < 1e-6 && max(abs(part$b_age + 0.45)) < 1e-6, "an extra age effect: share mediated b_rest / b_age")
check(is.null(age_mediation_design(designs[designs$design == "free_age", ])), "no split for a design that frees age")

# ---- 12. aeration ------------------------------------------------------------------------------
message("12. Compliance as aerated volume")
# strain on log VT - log(aerated volume), aerated = f x PFVC: slope = b (lambda + log VT/PBW - log f)
aer_pat <- exact$patients; set.seed(51)
aer_pat$lbmi <- rnorm(nrow(aer_pat), 0, 0.15)
aer_pat$laer0 <- ifelse(aer_pat$vent == 1, rnorm(nrow(aer_pat), 0, 0.3), 0)
aer_pat$lv0 <- ifelse(aer_pat$vent == 1, -0.5 * aer_pat$laer0 + rnorm(nrow(aer_pat), 0, 0.05), 0)   # VT lowered for stiff lungs
aer_pan <- exact$panel; i <- match(aer_pan$id, aer_pat$id); on_v <- aer_pan$phase == "on" & aer_pat$vent[i] == 1
aer_pan$y <- 3 + 0.1 * (i %% 7) + ifelse(on_v, b * (aer_pat$ldisc[i] + aer_pat$lv0[i] - aer_pat$laer0[i]), 0) * aer_pan$t
aer <- fit_aeration(designs[designs$design %in% c("unadjusted", "control_did"), ], aer_pat, aer_pan)
check(max(abs(c(aer$estimate, aer$vt_gradient) - b)) < 1e-6 && max(abs(aer$aeration_gradient + b)) < 1e-6,
      "with aeration in the model, the lambda, VT/PBW and aeration gradients are b, b and -b")
aer_pan$y <- aer_pan$y + rnorm(nrow(aer_pan), 0, 0.05)    # the Wald test needs noise
aer <- fit_aeration(designs[designs$design %in% c("unadjusted", "control_did"), ], aer_pat, aer_pan)
check(min(aer$p_strain_restriction) > 0.01, "and, with noise, the strain restriction is not rejected")
check(max(abs(aer$base_vt_gradient - b)) > 0.05, "without it, tidal volume lowered for stiff lungs biases the VT/PBW gradient")

# ---- 13. longitudinal: tidal volume as it changes within a patient ---------------------------------
message("13. Within-patient changes in tidal volume")
# step integrals against a direct sum
seg <- step_segments(pid = c(1, 1, 1, 2), s = c(0, 1, 2.5, 0.5), v = c(2, NA, 3, 5), hold = 1, upto = c(3, 10))
check(identical(round(seg$end, 10), c(1, 3, 1.5)) , "a value holds to the next record, at most `hold`, not past the end")
check(abs(step_integral(seg, c(1, 1, 2), c(0.5, 3, 1)) - c(1, 2 + 0.5 * 3, 2.5)) |> max() < 1e-12, "the step integral matches a hand sum")
# a fabricated course: patients with their own level, slope and lambda, tidal volume changed
# by the clinician and in response to a severity that also lowers the marker
set.seed(61)
np <- 300; times <- seq(0.25, 6, by = 0.5)
lf <- expand.grid(t = times, pid = seq_len(np)); lf$id <- paste0("q", lf$pid)
lam <- rnorm(np, 0, 0.15); lf$ldisc <- lam[lf$pid]
# the regressors as integrals of piecewise-constant series sampled on a fine grid
fine <- expand.grid(s = seq(0, 6, by = 1 / 24), pid = seq_len(np))
sev <- ave(rnorm(nrow(fine), 0, 0.3), fine$pid, FUN = function(x) stats::filter(x, 0.97, method = "recursive"))
fine$z <- sev; fine$x <- -0.5 * c(0, sev[-length(sev)]) + ave(rnorm(nrow(fine), 0, 0.02), fine$pid, FUN = cumsum)
fine$x <- fine$x - ave(fine$x, fine$pid)
cumint <- function(v) ave(v / 24, fine$pid, FUN = cumsum)
fine$Cx <- cumint(fine$x); fine$Cz <- cumint(fine$z)
at <- match(paste(lf$pid, lf$t), paste(fine$pid, round(fine$s, 6)))
lf$C <- fine$Cx[at]; lf$Z_sev <- fine$Cz[at]; lf$K <- lf$t
b_true <- -0.06
lf$y <- 5 + rnorm(np)[lf$pid] + (rnorm(np, 0, 0.05) + 0.3 * lam)[lf$pid] * lf$t +
  as.numeric(splines::ns(lf$t, df = 3) %*% c(-0.1, -0.2, 0.05)) + b_true * lf$C - 0.1 * lf$Z_sev   # a common curve in the model's span
fit_a <- fit_longitudinal(lf, c("C", "K", "Z_sev"))
check(abs(fit_a$beta[["C"]] - b_true) < 1e-8, "a planted effect of tidal volume is recovered")
small <- lf[lf$pid <= 40, ]
fit_s <- fit_longitudinal(small, c("C", "Z_sev"))
tb <- splines::ns(small$t, df = 3)
ref <- lm(small$y ~ small$C + small$Z_sev + tb + factor(small$pid) + factor(small$pid):small$t)
check(abs(fit_s$beta[["C"]] - coef(ref)[["small$C"]]) < 1e-8, "the estimate equals least squares with a patient intercept and slope")
lf2 <- lf; lf2$y <- lf2$y + 0.5 * lf2$ldisc * lf2$t
check(abs(fit_longitudinal(lf2, c("C", "K", "Z_sev"))$beta[["C"]] - b_true) < 1e-8, "a lambda x time effect leaves it unchanged")
check(abs(fit_longitudinal(lf, c("C", "K"))$beta[["C"]] - b_true) > 0.02, "tidal volume lowered as severity rises biases it without the severity term")

# end to end: records in, the frame's integrals and the estimate out
set.seed(71)
ne <- 120; hours <- seq(0, 7 * 24 - 1)
recs <- expand.grid(h = hours, k = seq_len(ne)); recs$id <- paste0("e", recs$k); recs$s <- recs$h / 24
step_vt <- function(k) { n_ch <- sample(0:3, 1); at <- sort(sample(hours, n_ch)); lv <- log(450) + cumsum(c(0, rnorm(n_ch, 0, 0.15)))
                         lv[findInterval(hours, c(-1, at))] }
recs$lvt <- unlist(lapply(seq_len(ne), step_vt))
recs$mode <- "assist control-volume control"; recs$vt <- exp(recs$lvt); recs$fio2 <- 0.5; recs$peep <- 8
e_pat <- data.frame(id = paste0("e", seq_len(ne)), vent = 1, t_off_d = NA_real_, ldisc = rnorm(ne, 0, 0.1))
xbar <- tapply(recs$lvt, recs$k, mean)
truth_C <- function(k, t) { h <- hours[hours / 24 < t]; sum((recs$lvt[recs$k == k][h + 1] - xbar[k]) * pmin(1 / 24, t - h / 24)) }
e_pan <- expand.grid(t = seq(0.3, 6.8, by = 0.5), k = seq_len(ne)); e_pan$id <- paste0("e", e_pan$k); e_pan$phase <- "on"
e_pan$C_true <- mapply(truth_C, e_pan$k, e_pan$t)
e_pan$y <- 5 + rnorm(ne)[e_pan$k] + rnorm(ne, 0, 0.05)[e_pan$k] * e_pan$t - 0.07 * e_pan$C_true
e_frame <- longitudinal_frame(e_pat, e_pan, list(vent = recs[, c("id", "s", "mode", "vt", "fio2", "peep")],
                                                pressor = data.frame(id = character(), s = numeric(), dose = numeric())))
check(max(abs(e_frame$C - e_pan$C_true[match(paste(e_frame$id, e_frame$t), paste(e_pan$id, e_pan$t))])) < 1e-9,
      "the frame's C is the time integral of log VT about the patient's own mean")
e_fit <- fit_longitudinal(e_frame, c("C", "C_lead", "K"))
check(abs(e_fit$beta[["C"]] + 0.07) < 1e-6 && abs(e_fit$beta[["C_lead"]]) < 1e-6, "records to estimate: the planted effect, and none on future tidal volume")

message("\nall tests passed")
