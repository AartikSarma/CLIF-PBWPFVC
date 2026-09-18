# =============================================================================
# Script 34_tte_titration: CO-PRIMARY TTE -- a one-step titration TOWARD PFVC dosing (LMTP)
# =============================================================================
# The head-to-head ceiling TTE (33_tte_ceiling_*) is the primary: two prescribable protocols. This
# leaf is its co-primary and the bedside-intuitive version of the same question, estimated by
# doubly-robust LMTP as a MODIFIED TREATMENT POLICY on the yardstick clinicians use (VT/PBW):
#
#   step_to_pfvc (PRIMARY policy): on each ventilated day, if the delivered VT/PBW is above the
#       patient's PFVC-anchored target, turn it down by at most DELTA mL/kg PBW toward that
#       target; patients already at or below their target are untouched. The target is the
#       bite-matched PFVC ceiling of the head-to-head (TAU_PFVC % of predicted FVC, the same
#       quantile rule as 33_tte_ceiling_common; PBWPFVC_VT_TARGET_PFVC overrides, e.g. 11 = ARMA),
#       expressed per patient on the VT/PBW scale: target_i = TAU_PFVC x 10 x PFVC_i / PBW_i.
#       This IS "dose toward PFVC" as a feasible shift: it acts only where PBW over-doses, by an
#       amount bounded by DELTA so the shifted value stays inside observed support (the MTP
#       positivity condition), and it leaves concordant patients alone.
#   additive (sensitivity): VT/PBW -> max(VT/PBW - DELTA, FLOOR) for everyone -- the plain
#       "half a click down" shift. It answers "does a small VT reduction help?", with the
#       normalizer entering only through the discordance gradient; kept as the sensitivity
#       that separates "toward PFVC" from "less VT for all".
#
# WHY A SHIFT AND NOT A STATIC TARGET: VT/PFVC is near-deterministic in demographics, so a static
# "everyone <= 11%" contrast has no overlap (the cross-sectional positivity wall). A bounded step
# stays in support by construction; the density-ratio trim and the ATE-width warning below are
# the positivity diagnostics. Per tte_within_demographic_hte the discordance gradient is demographically patterned
# (who is re-dosed), not a mechanism claim.
#
# Estimator: lmtp_tmle (binomial, mtp = TRUE; shift bites only on on-vent-and-alive days; 28-day
# all-cause mortality); CATE by discordance and PFVC via the DR-learner (Kennedy): per-patient DR
# pseudo-outcome, ITE = pseudo(shift) - pseudo(natural), spline of the log-modifier -> continuous
# CATE, per-log slope and p90-p10 gradient (bootstrap CI). Age enters every nuisance model as
# ns(age10, 4). In 32_tte_run_all as the co-primary. Synthetic mortality is simulated => plumbing only.
# Env: PBWPFVC_VT_POLICY (step_to_pfvc|additive), PBWPFVC_VT_TARGET_PFVC, PBWPFVC_CEIL_BITE,
# PBWPFVC_VT_{DELTA,FLOOR,HORIZON,FOLDS,LEARNERS,BOOT}.
# =============================================================================
library(here); library(lmtp); library(splines); library(patchwork)
# lmtp reports progress via progressr, OFF by default -> enable a handler so the per-fold/timepoint
# bar actually renders (PBWPFVC_PROGRESS=0 to silence). cli if available, else base txtprogressbar.
if (!identical(Sys.getenv("PBWPFVC_PROGRESS", "1"), "0") && requireNamespace("progressr", quietly = TRUE)) {
  progressr::handlers(global = TRUE)
  progressr::handlers(if (requireNamespace("cli", quietly = TRUE)) "cli" else "txtprogressbar")
}
source(here::here("code", "31_tte_engine.R"))

POLICY <- Sys.getenv("PBWPFVC_VT_POLICY", "step_to_pfvc"); stopifnot(POLICY %in% c("step_to_pfvc", "additive"))
DELTA <- as.numeric(Sys.getenv("PBWPFVC_VT_DELTA", "0.5"))    # maximum step, mL/kg PBW (half a click)
FLOOR <- as.numeric(Sys.getenv("PBWPFVC_VT_FLOOR", "4.0"))    # never reduce VT/PBW below this (feasible MTP)
K     <- as.integer(Sys.getenv("PBWPFVC_VT_HORIZON", "14"))
FOLDS <- as.integer(Sys.getenv("PBWPFVC_VT_FOLDS", "5"))
BGRAD <- as.integer(Sys.getenv("PBWPFVC_VT_BOOT", "1000"))    # bootstrap reps for the CATE-gradient CI
# Richer nuisance library: flexible additive (gam) + forest (ranger) so a uniform confounding bias
# moves the LEVEL but the gradient is read off a well-specified fit. Any learner whose backend pkg
# is absent is dropped (glmnet often is); SL.glm + SL.mean are always retained.
# SuperLearner's SL.ranger defaults to num.threads=1 (single-core!) -- the reason the forest-heavy
# LMTP crawls. This wrapper lets each ranger fit use all cores (override via PBWPFVC_CORES). Use
# THIS *or* a future multisession plan, not both -- combining oversubscribes (folds x all-core ranger).
NTHREAD <- suppressWarnings(as.integer(Sys.getenv("PBWPFVC_CORES", unset = NA)))
if (is.na(NTHREAD)) NTHREAD <- max(1L, parallel::detectCores())
SL.ranger.mc <- function(...) SuperLearner::SL.ranger(..., num.threads = NTHREAD)
LRN   <- strsplit(Sys.getenv("PBWPFVC_VT_LEARNERS", "SL.glm,SL.gam,SL.ranger.mc,SL.mean"), ",")[[1]]
.lrn_pkg <- c(SL.gam = "gam", SL.ranger = "ranger", SL.ranger.mc = "ranger", SL.glmnet = "glmnet",
              SL.earth = "earth", SL.xgboost = "xgboost", SL.nnet = "nnet")
.avail <- vapply(LRN, function(l) { p <- unname(.lrn_pkg[l]); is.na(p) || requireNamespace(p, quietly = TRUE) }, logical(1))
if (any(!.avail)) message("34_tte_titration: dropping unavailable learners: ", paste(LRN[!.avail], collapse = ", "))
LRN <- unique(c(LRN[.avail], "SL.glm", "SL.mean"))
# STABLE treatment-density learner: log(VT/PBW) is near-constant (protocolized ~6-8), so asking a
# forest/gam to estimate its conditional DENSITY (lmtp's classification trick) explodes the density
# ratio -> the EIF variance blows up (ATE CI spanned +-100pp). A simple parametric trt model is far
# more stable here; the rich library stays on the OUTCOME side. Plus density-ratio trim + ITE winsor.
LRN_TRT <- strsplit(Sys.getenv("PBWPFVC_VT_TRT_LEARNERS", "SL.glm,SL.mean"), ",")[[1]]
TRIM    <- as.numeric(Sys.getenv("PBWPFVC_VT_TRIM", "0.99"))   # lmtp_control(.trim): bound the density ratio
WINSOR  <- as.numeric(Sys.getenv("PBWPFVC_VT_WINSOR", "1"))    # clamp ITEs to +-WINSOR before the CATE spline
cat(sprintf("    LMTP learners: outcome={%s} trt={%s}; ranger num.threads=%d; trim=%.3f, winsor=+-%.1f\n",
            paste(LRN, collapse = ","), paste(LRN_TRT, collapse = ","), NTHREAD, TRIM, WINSOR))

# --- daily VT/PBW (mL/kg) + survival nodes (binomial, on-vent-only) ---------------------------
vtd <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet")) %>%
  inner_join(base %>% select(hospitalization_id, t0, pbw), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 1, vent_day <= K, !is.na(tidal_volume_set), tidal_volume_set >= 50,
         tidal_volume_set <= 2000) %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(vt = median(tidal_volume_set), pbw = first(pbw), .groups = "drop") %>%
  transmute(hospitalization_id, vent_day, A = log(vt / pbw)) %>% filter(is.finite(A))   # log(mL/kg)

# PFVC-anchored target: the head-to-head's bite-matched ceiling (same rule as 33_tte_ceiling_common:
# the (1 - BITE) quantile of daily VT/PFVC over post-grace patient-days), unless overridden.
BITE <- as.numeric(Sys.getenv("PBWPFVC_CEIL_BITE", "0.25"))
TAU_ENV <- suppressWarnings(as.numeric(Sys.getenv("PBWPFVC_VT_TARGET_PFVC", unset = NA)))
TAU_PFVC <- if (is.finite(TAU_ENV)) TAU_ENV else unname(quantile(panel$vtpfvc[panel$vent_day > GRACE], 1 - BITE))
state <- panel %>% select(hospitalization_id, vent_day, sf, map, on_pressor)
demo  <- base %>% select(hospitalization_id, age10, sex_category, race_category,
                         height_grp, sofa_total, death_day, imv_extub_day)
ids <- unique(vtd$hospitalization_id)
long <- tidyr::expand_grid(hospitalization_id = ids, vent_day = 1:K) %>%
  left_join(vtd, by = c("hospitalization_id", "vent_day")) %>%
  left_join(state, by = c("hospitalization_id", "vent_day")) %>%
  left_join(demo, by = "hospitalization_id") %>%
  group_by(hospitalization_id) %>% arrange(vent_day) %>%
  mutate(alive = is.na(death_day) | vent_day <= death_day, exposed = as.integer(!is.na(A) & alive)) %>%
  fill(A, sf, map, .direction = "downup") %>% mutate(on_pressor = coalesce(on_pressor, 0L)) %>% ungroup()
bad <- long %>% filter(is.na(A) | is.na(sf) | is.na(map)) %>% pull(hospitalization_id) %>% unique()
long <- long %>% filter(!hospitalization_id %in% bad)

mods <- base %>% transmute(hospitalization_id, discord = pbw / pfvc, pfvc = pfvc,
                           # the patient's PFVC-anchored target on the VT/PBW scale (mL/kg):
                           # VT/PFVC% = VT/PFVC x 0.1  =>  VT/PBW = VT/PFVC% x 10 x PFVC/PBW
                           vt_target = TAU_PFVC * 10 * pfvc / pbw,
                           Y = as.integer(!is.na(death_day) & death_day <= 28)) %>%
  filter(!is.na(discord), !is.na(pfvc))
wide <- long %>%
  select(hospitalization_id, vent_day, A, sf, map, on_pressor, exposed,
         age10, sex_category, race_category, height_grp, sofa_total) %>%
  pivot_wider(id_cols = c(hospitalization_id, age10, sex_category, race_category, height_grp, sofa_total),
              names_from = vent_day, values_from = c(A, sf, map, on_pressor, exposed), names_sep = "_") %>%
  inner_join(mods, by = "hospitalization_id") %>% as.data.frame()
cat(sprintf("34_tte_titration [%s]: %d patients; PBW/PFVC discordance median %.2f, PFVC median %.2f L; PFVC target %.1f%% (VT/PBW-scale target median %.2f mL/kg); step %.2f mL/kg, floor %.1f, K %d, folds %d\n",
            POLICY, nrow(wide), median(wide$discord), median(wide$pfvc), TAU_PFVC, median(wide$vt_target), DELTA, FLOOR, K, FOLDS))

# age enters the nuisance models as a natural cubic spline: lmtp takes column names, so the
# basis is materialized as columns (fixed knots from the analytic cohort)
age_basis <- splines::ns(wide$age10, 4); for (j in 1:4) wide[[paste0("age_ns", j)]] <- age_basis[, j]
# vt_target is a baseline covariate: the shift depends on it, and lmtp hands the shift
# function only the model columns (a non-model column comes back NULL and empties the shift)
base_cov <- c(paste0("age_ns", 1:4), "sex_category", "race_category", "height_grp", "sofa_total", "vt_target")
tv <- lapply(1:K, function(t) paste0(c("sf", "map", "on_pressor", "exposed"), "_", t))
trt <- paste0("A_", 1:K)
# The policy on the VT/PBW scale (mL/kg), applied on the log scale the density learners see.
#   step_to_pfvc: A' = max(FLOOR, max(A - DELTA, min(A, target)))  -- down by at most DELTA, never
#                 below the patient's PFVC target, untouched if already at or below it
#   additive:     A' = max(A - DELTA, FLOOR) for everyone
apply_policy <- function(a, target) {
  if (POLICY == "step_to_pfvc") pmax(FLOOR, pmax(a - DELTA, pmin(a, target))) else pmax(a - DELTA, FLOOR)
}
shift_fun <- function(data, t) { i <- sub("^A_", "", t); ex <- data[[paste0("exposed_", i)]]
  out <- data[[t]]; reduced <- log(apply_policy(exp(out), data[["vt_target"]]))
  out[ex == 1L] <- reduced[ex == 1L]; out }
# realized intervention intensity (audit the feasible-MTP bite)
.exp_days <- long %>% filter(exposed == 1L) %>%
  inner_join(mods %>% select(hospitalization_id, vt_target, discord), by = "hospitalization_id") %>%
  mutate(raw = exp(A), shifted = apply_policy(raw, vt_target), binds = shifted < raw - 1e-9,
         reaches_target = POLICY == "step_to_pfvc" & binds & shifted <= vt_target + 1e-9)
bite_tbl <- .exp_days %>%
  mutate(disc_grp = cut(discord, quantile(discord, c(0, 1/3, 2/3, 1)), include.lowest = TRUE,
                        labels = c("Concordant", "Mid", "Discordant"))) %>%
  group_by(disc_grp) %>%
  summarise(n_days = n(), frac_days_shifted = mean(binds), mean_cut_mlkg = mean(raw - shifted),
            frac_shifted_days_reaching_target = if (any(binds)) mean(reaches_target[binds]) else NA_real_,
            median_vtpbw = median(raw), median_target = median(vt_target), .groups = "drop") %>%
  mutate(policy = POLICY, tau_pfvc = TAU_PFVC, delta = DELTA, site = site_name)
write_csv(bite_tbl, file.path(final_dir, paste0("vtpbw_titration_bite_", site_name, ".csv")))
cat(sprintf("    realized cut: %.1f%% of exposed days shifted; median VT/PBW %.2f -> %.2f mL/kg; mean cut %.2f mL/kg; %.1f%% of shifted days reach the PFVC target in one step\n",
            100 * mean(.exp_days$binds), median(.exp_days$raw), median(.exp_days$shifted),
            mean(.exp_days$raw - .exp_days$shifted),
            100 * mean(.exp_days$reaches_target[.exp_days$binds])))
print(as.data.frame(bite_tbl %>% transmute(disc_grp, n_days, pct_shifted = round(100 * frac_days_shifted), mean_cut = round(mean_cut_mlkg, 2),
        pct_reach_target = round(100 * frac_shifted_days_reaching_target))), row.names = FALSE)

args <- list(data = wide, trt = trt, outcome = "Y", baseline = base_cov, time_vary = tv,
             mtp = TRUE, outcome_type = "binomial", learners_outcome = LRN, learners_trt = LRN_TRT,
             folds = FOLDS, control = lmtp_control(.trim = TRIM))
.t <- Sys.time()
message("[", format(Sys.time(), "%H:%M:%S"), "] fitting SHIFT-policy LMTP (slow: ", K, " timepoints x ", FOLDS, " folds x 2 nuisance models) ...")
fs <- do.call(lmtp_tmle, c(args, list(shift = shift_fun)))
message("[", format(Sys.time(), "%H:%M:%S"), "] shift done (", round(difftime(Sys.time(), .t, units = "mins"), 1), " min); fitting NATURAL policy ...")
.t <- Sys.time()
fo <- do.call(lmtp_tmle, c(args, list(shift = NULL)))
message("[", format(Sys.time(), "%H:%M:%S"), "] natural done (", round(difftime(Sys.time(), .t, units = "mins"), 1), " min)")
ate <- lmtp_contrast(fs, ref = fo, type = "additive")$estimates

# DR pseudo-outcomes for EACH policy (uncentered = estimate + EIF). The difference is the
# per-patient ITE (absolute scale); keeping them separate also allows the RELATIVE scale below.
ps_dr <- fs$estimate@x + fs$estimate@eif      # pseudo-outcome under the shifted policy
pn_dr <- fo$estimate@x + fo$estimate@eif      # pseudo-outcome under the natural course
ite <- ps_dr - pn_dr
stopifnot(length(ite) == nrow(wide))
cat(sprintf("    ATE [%s, step <= %.2f mL/kg]: risk_natural %.1f%% (CANARY ~27-32 real), RD %+.2f pp [%.2f, %.2f]\n",
            POLICY, DELTA, 100 * ate$ref, 100 * ate$estimate, 100 * ate$conf.low, 100 * ate$conf.high))
write_csv(tibble(policy = POLICY, tau_pfvc = TAU_PFVC, delta = DELTA, floor = FLOOR, horizon_days = K,
                 n_patients = nrow(wide), risk_natural = ate$ref, risk_policy = ate$ref + ate$estimate,
                 rd = ate$estimate, rd_lo = ate$conf.low, rd_hi = ate$conf.high, rd_se = ate$std.error,
                 frac_exposed_days_shifted = mean(.exp_days$binds), site = site_name),
          file.path(final_dir, paste0("vtpbw_titration_ate_", site_name, ".csv")))
if (abs(ate$conf.high - ate$conf.low) > 0.40)
  message("  WARNING: ATE CI width ", round(100 * (ate$conf.high - ate$conf.low)), "pp -- VT/PBW shift still poorly identified ",
          "(positivity: protocolized treatment). Read as 'no estimable PBW-scale contrast', not a result.")
# winsorize EIF tails to tame the CATE SHAPE variance, THEN recenter to the unbiased TMLE ATE so the
# curve's level is honest (winsorizing the heavy negative tail otherwise inflates the level; the
# gradient/shape is a contrast and is unchanged by the recentering).
n_wins <- sum(abs(ite) > WINSOR)
itew <- pmin(pmax(ite, -WINSOR), WINSOR)
itew <- itew - mean(itew) + ate$estimate
cat(sprintf("    winsorized %d/%d ITEs to +-%.1f, recentered CATE to TMLE ATE %.2f pp (gradient/shape unaffected)\n",
            n_wins, length(ite), WINSOR, 100 * ate$estimate))

# --- CATE by each modifier: spline of ITE on log(modifier); slope + p90-p10 gradient -----------
# The CATE is a risk difference (shift - natural), so MORE NEGATIVE = MORE BENEFIT:
# discordance: signal expected as a DOWN-slope (more benefit where PBW oversizes the lung).
# pfvc:        signal expected as an UP-slope (more benefit at smaller absolute lung).
# Both scales. The absolute CATE (risk difference) is modified by baseline risk: a constant
# RELATIVE effect in a higher-risk group mechanically produces a larger RD, and the modifiers
# here (high discordance, small PFVC) mark the highest-risk patients. So the relative curve is
# the check that decides whether a sloped RD curve is targeting or baseline risk: RR(modifier) =
# E[pseudo_shift | modifier] / E[pseudo_natural | modifier], each fitted on the SAME spline basis.
# A flat RR with a sloped RD => baseline risk. A sloped RR => genuine effect modification.
cate_one <- function(mod, lab) {
  dd  <- tibble(lx = log(wide[[mod]]), ite = itew, ps = ps_dr, pn = pn_dr) %>% filter(is.finite(lx))
  kn  <- attr(ns(dd$lx, 3), "knots"); bd <- attr(ns(dd$lx, 3), "Boundary.knots")  # FIX basis
  fit_spl <- function(d) lm(ite ~ ns(lx, knots = kn, Boundary.knots = bd), data = d)
  spl <- fit_spl(dd)
  # relative scale: fit each policy's pseudo-outcome mean, take the ratio of the fitted curves.
  # Pseudo-outcomes are risks smeared by the influence function, so individual values can fall
  # outside [0,1]; the fitted MEANS are the estimable objects, and the ratio is reported only
  # where the fitted natural-course risk exceeds RR_FLOOR (dividing by a near-zero risk is noise).
  RR_FLOOR <- 0.05
  fit_rr <- function(d) {
    a <- lm(ps ~ ns(lx, knots = kn, Boundary.knots = bd), data = d)
    b <- lm(pn ~ ns(lx, knots = kn, Boundary.knots = bd), data = d)
    function(newd) { num <- predict(a, newd); den <- predict(b, newd)
                     ifelse(den > RR_FLOOR, num / den, NA_real_) }
  }
  rr_f <- fit_rr(dd)
  tr  <- summary(lm(ite ~ lx, data = dd))$coefficients["lx", ]   # per-log-unit linear slope
  g   <- tibble(lx = quantile(dd$lx, seq(0.1, 0.9, 0.1)))
  pr  <- predict(spl, g, se.fit = TRUE)
  tbl <- g %>% transmute(modifier = lab, pctile = seq(10, 90, 10), value = round(exp(lx), 2),
                         cate_pp = round(100 * pr$fit, 2),
                         lo = round(100 * (pr$fit - 1.96 * pr$se.fit), 2),
                         hi = round(100 * (pr$fit + 1.96 * pr$se.fit), 2))
  qg  <- quantile(dd$lx, c(0.1, 0.9))
  grad_pt <- diff(as.numeric(predict(spl, tibble(lx = qg))))
  rr_pt <- as.numeric(rr_f(tibble(lx = qg)))          # RR at p10 and p90
  rr_ratio_pt <- rr_pt[2] / rr_pt[1]                  # >1 or <1 = the relative effect differs
  bsm <- vapply(seq_len(BGRAD), function(b) {
    d2 <- dd[sample.int(nrow(dd), replace = TRUE), ]
    fb <- tryCatch(fit_spl(d2), error = function(e) NULL)
    g  <- if (is.null(fb)) NA_real_ else diff(as.numeric(predict(fb, tibble(lx = qg))))
    rb <- tryCatch({ r <- as.numeric(fit_rr(d2)(tibble(lx = qg))); r[2] / r[1] }, error = function(e) NA_real_)
    c(g, rb)
  }, numeric(2))
  bs <- bsm[1, ][is.finite(bsm[1, ])]
  bs_rr <- bsm[2, ][is.finite(bsm[2, ])]
  slope <- tibble(modifier = lab,
    statistic   = c("per_log_unit_slope", "grad_p90_minus_p10"),
    estimate_pp = round(100 * c(tr["Estimate"], grad_pt), 2),
    lo = round(100 * c(tr["Estimate"] - 1.96 * tr["Std. Error"], quantile(bs, 0.025)), 2),
    hi = round(100 * c(tr["Estimate"] + 1.96 * tr["Std. Error"], quantile(bs, 0.975)), 2),
    p  = c(signif(tr["Pr(>|t|)"], 2), NA_real_), n_boot = c(NA_integer_, length(bs)))
  # the relative-scale companion to grad_p90_minus_p10: the RATIO of risk ratios at p90 vs p10
  slope <- bind_rows(slope, tibble(modifier = lab, statistic = "rr_ratio_p90_over_p10",
    estimate_pp = round(rr_ratio_pt, 3),
    lo = round(unname(quantile(bs_rr, 0.025)), 3), hi = round(unname(quantile(bs_rr, 0.975)), 3),
    p = NA_real_, n_boot = length(bs_rr)))
  fg  <- tibble(lx = seq(quantile(dd$lx, 0.02), quantile(dd$lx, 0.98), length.out = 60))
  fp  <- predict(spl, fg, se.fit = TRUE)
  fig <- tibble(modifier = lab, x = exp(fg$lx), cate = 100 * fp$fit,
                lo = 100 * (fp$fit - 1.96 * fp$se.fit), hi = 100 * (fp$fit + 1.96 * fp$se.fit))
  # relative curve + its bootstrap band, on the same grid
  rr_grid <- as.numeric(rr_f(fg))
  rr_bs <- vapply(seq_len(min(BGRAD, 400)), function(b) {
    d2 <- dd[sample.int(nrow(dd), replace = TRUE), ]
    tryCatch(as.numeric(fit_rr(d2)(fg)), error = function(e) rep(NA_real_, nrow(fg)))
  }, numeric(nrow(fg)))
  fig_rr <- tibble(modifier = lab, x = exp(fg$lx), rr = rr_grid,
                   lo = apply(rr_bs, 1, quantile, 0.025, na.rm = TRUE),
                   hi = apply(rr_bs, 1, quantile, 0.975, na.rm = TRUE))
  # rug: the modifier's observed distribution, so the reader sees where the spline is constrained
  rug <- tibble(modifier = lab, x = exp(dd$lx))
  list(tbl = tbl, trend = tr, fig = fig, fig_rr = fig_rr, rug = rug, slope = slope)
}
ct_disc <- cate_one("discord", "PBW/PFVC discordance")
ct_pfvc <- cate_one("pfvc",    "PFVC (predicted size, L)")
cate_tbl <- bind_rows(ct_disc$tbl, ct_pfvc$tbl) %>% mutate(policy = POLICY, site = site_name)
write_csv(cate_tbl, file.path(final_dir, paste0("vtpbw_titration_cate_", site_name, ".csv")))
slope_tbl <- bind_rows(ct_disc$slope, ct_pfvc$slope) %>% mutate(policy = POLICY, site = site_name)
write_csv(slope_tbl, file.path(final_dir, paste0("vtpbw_titration_slope_", site_name, ".csv")))

# --- figure: one panel per modifier, ATE reference line ---------------------------------------
ate_pp <- 100 * ate$estimate   # official TMLE ATE as the figure reference (not the winsorized mean)
POLICY_LAB <- sprintf("%s (step <= %.2f mL/kg PBW)", POLICY, DELTA)
mk_panel <- function(ct, xlab) {
  ggplot(ct$fig, aes(x, cate)) +
    geom_hline(yintercept = ate_pp, linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, fill = "#0072B2") +
    geom_line(colour = "#0072B2", linewidth = 1) +
    geom_rug(data = ct$rug %>% slice_sample(n = min(nrow(ct$rug), 2000)), aes(x = x),
             inherit.aes = FALSE, alpha = 0.06, length = unit(0.03, "npc")) +
    labs(x = xlab, y = sprintf("CATE: 28-d mortality RD of the %s policy (pp)", POLICY), title = ct$fig$modifier[1]) +
    theme_minimal(base_size = 11)
}
mk_panel_rr <- function(ct, xlab) {
  ggplot(ct$fig_rr, aes(x, rr)) +
    geom_hline(yintercept = 1, colour = "grey80") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, fill = "#D55E00") +
    geom_line(colour = "#D55E00", linewidth = 1) +
    geom_rug(data = ct$rug %>% slice_sample(n = min(nrow(ct$rug), 2000)), aes(x = x),
             inherit.aes = FALSE, alpha = 0.06, length = unit(0.03, "npc")) +
    scale_y_log10() +
    labs(x = xlab, y = "Relative CATE: risk ratio (policy / natural)", title = NULL) +
    theme_minimal(base_size = 11)
}
rr_tbl <- bind_rows(ct_disc$fig_rr, ct_pfvc$fig_rr) %>% mutate(site = site_name)
write_csv(rr_tbl, file.path(final_dir, paste0("vtpbw_titration_cate_rr_", site_name, ".csv")))
# 2 x 2: absolute risk difference (top) and relative risk ratio (bottom) for each modifier,
# with a rug of the modifier's distribution under every panel. Read the two rows together: an
# absolute curve that bends where the relative curve is flat is baseline risk, not targeting.
fig <- (mk_panel(ct_disc, "PBW/PFVC discordance - higher = PBW oversizes (PFVC says smaller)") | mk_panel(ct_pfvc, "PFVC (L) - lower = smaller predicted lung")) /
       (mk_panel_rr(ct_disc, "PBW/PFVC discordance - higher = PBW oversizes (PFVC says smaller)") | mk_panel_rr(ct_pfvc, "PFVC (L) - lower = smaller predicted lung")) +
  patchwork::plot_annotation(
    title = paste0("CATE of the ", POLICY, " titration (step <= ", DELTA, " mL/kg PBW), by PFVC-derived modifiers - ", site_name, if (is_synthetic) " (SYNTHETIC)" else ""),
    subtitle = paste0("Top: absolute CATE (risk difference); dashed = ATE; negative = benefit. ",
                      "Bottom: relative CATE (risk ratio, policy / natural course); 1 = no effect. ",
                      "Rug = the modifier's distribution. A sloped absolute curve with a flat relative curve is ",
                      "baseline risk (the modifiers mark the highest-risk patients), not effect modification."))
ggsave(file.path(final_dir, paste0("vtpbw_titration_cate_", site_name, ".pdf")), fig, width = 12, height = 9)

cat(sprintf("\n=== 34_tte_titration [%s] continuous CATE by PFVC-derived modifiers ===\n", POLICY))
cat("--- effect-modification trends (read the SHAPE, not the level: confounding-by-severity inflates the ATE) ---\n")
print(as.data.frame(slope_tbl), row.names = FALSE)
cat("    (discordance grad_p90_minus_p10 MORE NEGATIVE => the misdosed gain more from a fixed bedside VT/PBW cut;\n")
cat("     per_log slope CI = analytic regression; gradient CI = patient bootstrap -- both conditional on the DR pseudo-outcomes)\n")
print(as.data.frame(cate_tbl), row.names = FALSE)
message("Wrote vtpbw_titration_{cate,slope}_", site_name, ".csv + cate .pdf to ", final_dir)
