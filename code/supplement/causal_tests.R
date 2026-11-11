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

message("\nall tests passed")
