# =============================================================================
# Script 11.Y: CATE of an ADDITIVE MP/PBW power-reduction MTP, by PBW/PFVC discordance + PFVC
# =============================================================================
# The mechanical-power, LMTP-estimated SIBLING to 11.X. 11.X is a clone-censor-weight TTE of a
# VT/PFVC strain-limiting CEILING and shows the benefit concentrates in the misdosed (high
# PBW/PFVC). This leaf asks the SAME targeting question on a different exposure (total
# mechanical power, not just VT) with a different estimator (doubly-robust LMTP, not CCW), and
# under the PBW-anchored yardstick clinicians actually dose on. Two PFVC-derived modifiers:
#   1. PBW/PFVC discordance -- higher = PBW oversizes the lung (PFVC says it is smaller);
#   2. PFVC (absolute predicted size, L) -- lower = smaller lung regardless of the PBW gap.
# CATE is a mortality risk difference (shift - natural), so MORE NEGATIVE = more benefit. A
# DOWN-sloping CATE in discordance and/or an UP-sloping CATE in PFVC => power-reduction helps
# most in exactly the patients PBW over-sizes -- the MP analogue of the 11.X discordance gradient.
#
# WHY ADDITIVE (the normalizer-DEPENDENT choice, deliberately): a *multiplicative* "MP/PBW x
# (1-delta)" is normalizer-INVARIANT -- identical to "MP x (1-delta)" and to "MP/PFVC x
# (1-delta)", so the policy itself carries no PBW-vs-PFVC content and only the slicing of
# heterogeneity differs (that was the old 12.H framing; rejected). An *additive* policy
# "MP/PBW -> MP/PBW - DELTA" removes an absolute power of DELTA x PBW, i.e. in true-power-per-
# lung terms DELTA x (PBW/PFVC) = DELTA x discordance. So a fixed PBW-anchored decrement cuts
# the MISDOSED hardest per unit of actual lung -- the policy is normalizer-DEPENDENT on the
# discordance axis exactly as the 11.X VT/PFVC ceiling is, and the resulting discordance
# down-slope (more benefit at higher discordance) is the genuine MP parallel (NOT a relabelling artifact). DELTA is in (J/min)/kg PBW;
# default 0.05 ~ a 3 J/min absolute cut at PBW 60 kg. Feasible MTP: never reduced below FLOOR.
#
# Design = the validated 12.D/F estimator: doubly-robust binomial LMTP (lmtp, mtp=TRUE), shift
# bites only on on-vent-and-alive days, 28-day all-cause mortality. CATE via the DR-learner
# (Kennedy): fit shift + natural on the FULL cohort, form each patient's DR pseudo-outcome
# (estimate + influence value), ITE_i = pseudo(shift)_i - pseudo(natural)_i (mean = ATE/RD),
# then regress ITE on a spline of the log-modifier -> CATE(modifier); slope = the trend. CIs
# are the regression's (approximate), as in 12.F.
#
# HONEST CAVEAT (inherited from 12.E/F): both modifiers are deterministic in age/sex/race/
# height, so HTE-by-them is entangled with HTE-by-demographics -- this informs TARGETING, it
# cannot prove "strain beyond age". By the bedside logic it need not: if power-reduction helps
# the PFVC-flagged group more, reading PFVC is the better tool whatever the mechanism.
# PEAK-MP/PBW exposure. Age enters every nuisance model as ns(age10, 4), matching the engine and
# the VT titration. This is the POWER-axis analogue of the VT co-primary (11_vtpbw_titration): an
# additive PBW-anchored shift whose per-lung bite scales with discordance. Read it beside the MP
# head-to-head (tte_mp_ceiling_*), where the two normalizers' ceilings disagreed on only ~7% of
# days versus ~20% on the VT axis, because rate and pressure are shared by both normalizers and
# dilute the denominator difference -- so a flat CATE here is a small-contrast null, not evidence
# against the dosing hypothesis. Standalone supplement -- NOT in 11_run_all. Synthetic mortality
# is simulated => plumbing only. Env: PBWPFVC_MP_{DELTA,FLOOR,HORIZON,FOLDS,LEARNERS}.
# =============================================================================
library(here); library(lmtp); library(splines); library(patchwork)
# lmtp reports progress via progressr, OFF by default -> enable a handler so the per-fold/timepoint
# bar actually renders (PBWPFVC_PROGRESS=0 to silence). cli if available, else base txtprogressbar.
if (!identical(Sys.getenv("PBWPFVC_PROGRESS", "1"), "0") && requireNamespace("progressr", quietly = TRUE)) {
  progressr::handlers(global = TRUE)
  progressr::handlers(if (requireNamespace("cli", quietly = TRUE)) "cli" else "txtprogressbar")
}
source(here::here("code", "10_tte_engine.R"))

DELTA <- as.numeric(Sys.getenv("PBWPFVC_MP_DELTA", "0.05"))   # ADDITIVE cut, (J/min)/kg PBW (~3 J/min at 60 kg)
FLOOR <- as.numeric(Sys.getenv("PBWPFVC_MP_FLOOR", "0.10"))   # never reduce MP/PBW below this (feasible MTP)
K     <- as.integer(Sys.getenv("PBWPFVC_MP_HORIZON", "14"))
FOLDS <- as.integer(Sys.getenv("PBWPFVC_MP_FOLDS", "5"))      # richer cross-fitting than the 12.x prototypes
BGRAD <- as.integer(Sys.getenv("PBWPFVC_MP_BOOT", "1000"))    # bootstrap reps for the CATE-gradient CI
# Richer default nuisance library than the 12.x prototypes (SL.glm,SL.mean): a flexible additive
# (gam) + a forest (ranger) so a uniform confounding bias moves the LEVEL but the gradient is read
# off a well-specified fit. Any learner whose backend pkg is absent is dropped (glmnet often is);
# SL.glm + SL.mean are always retained so the library is never empty.
# SuperLearner's SL.ranger defaults to num.threads=1 (single-core!) -- the reason the forest-heavy
# LMTP crawls. This wrapper lets each ranger fit use all cores (override via PBWPFVC_CORES). Use
# THIS *or* a future multisession plan, not both -- combining oversubscribes (folds x all-core ranger).
NTHREAD <- suppressWarnings(as.integer(Sys.getenv("PBWPFVC_CORES", unset = NA)))
if (is.na(NTHREAD)) NTHREAD <- max(1L, parallel::detectCores())
SL.ranger.mc <- function(...) SuperLearner::SL.ranger(..., num.threads = NTHREAD)
LRN   <- strsplit(Sys.getenv("PBWPFVC_MP_LEARNERS", "SL.glm,SL.gam,SL.ranger.mc,SL.mean"), ",")[[1]]
.lrn_pkg <- c(SL.gam = "gam", SL.ranger = "ranger", SL.ranger.mc = "ranger", SL.glmnet = "glmnet",
              SL.earth = "earth", SL.xgboost = "xgboost", SL.nnet = "nnet")
.avail <- vapply(LRN, function(l) { p <- unname(.lrn_pkg[l]); is.na(p) || requireNamespace(p, quietly = TRUE) }, logical(1))
if (any(!.avail)) message("11.Y: dropping unavailable learners: ", paste(LRN[!.avail], collapse = ", "))
LRN <- unique(c(LRN[.avail], "SL.glm", "SL.mean"))
# STABLE treatment-density learner (matches 11_vtpbw): a forest/gam estimating the conditional
# DENSITY of the (peaked) treatment via lmtp's classification trick is unstable -> the density ratio
# and EIF variance blow up. Use a simple parametric trt model; keep the rich library on the OUTCOME
# side. Plus density-ratio trim + ITE winsor/recenter for a defensible, stable CATE.
LRN_TRT <- strsplit(Sys.getenv("PBWPFVC_MP_TRT_LEARNERS", "SL.glm,SL.mean"), ",")[[1]]
TRIM    <- as.numeric(Sys.getenv("PBWPFVC_MP_TRIM", "0.99"))   # lmtp_control(.trim): bound the density ratio
WINSOR  <- as.numeric(Sys.getenv("PBWPFVC_MP_WINSOR", "1"))    # clamp ITEs to +-WINSOR before the CATE spline
cat(sprintf("    LMTP learners: outcome={%s} trt={%s}; ranger num.threads=%d; trim=%.3f, winsor=+-%.1f\n",
            paste(LRN, collapse = ","), paste(LRN_TRT, collapse = ","), NTHREAD, TRIM, WINSOR))
KC <- 0.098

# --- daily peak-MP/PBW + survival nodes (binomial, on-vent-only) -- as in 12.D/F --------
mpd <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet")) %>%
  inner_join(base %>% select(hospitalization_id, t0, pbw), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 1, vent_day <= K, !is.na(tidal_volume_set), tidal_volume_set >= 50,
         tidal_volume_set <= 2000, !is.na(resp_rate_set), resp_rate_set >= 4, resp_rate_set <= 60,
         !is.na(peep_set), !is.na(peak_inspiratory_pressure_obs),
         between(peak_inspiratory_pressure_obs, 5, 80)) %>%
  mutate(mp_raw = KC * resp_rate_set * (tidal_volume_set / 1000) * peak_inspiratory_pressure_obs) %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(mp_raw = median(mp_raw), pbw = first(pbw), .groups = "drop") %>%
  transmute(hospitalization_id, vent_day, A = log(mp_raw / pbw)) %>% filter(is.finite(A))

state <- panel %>% select(hospitalization_id, vent_day, sf, map, on_pressor)
demo  <- base %>% select(hospitalization_id, age10, sex_category, race_category,
                         height_grp, sofa_total, death_day, imv_extub_day)
ids <- unique(mpd$hospitalization_id)
long <- tidyr::expand_grid(hospitalization_id = ids, vent_day = 1:K) %>%
  left_join(mpd, by = c("hospitalization_id", "vent_day")) %>%
  left_join(state, by = c("hospitalization_id", "vent_day")) %>%
  left_join(demo, by = "hospitalization_id") %>%
  group_by(hospitalization_id) %>% arrange(vent_day) %>%
  mutate(alive = is.na(death_day) | vent_day <= death_day, exposed = as.integer(!is.na(A) & alive)) %>%
  fill(A, sf, map, .direction = "downup") %>% mutate(on_pressor = coalesce(on_pressor, 0L)) %>% ungroup()
bad <- long %>% filter(is.na(A) | is.na(sf) | is.na(map)) %>% pull(hospitalization_id) %>% unique()
long <- long %>% filter(!hospitalization_id %in% bad)

# per-patient effect modifiers (the two the user asked for) + outcome
mods <- base %>% transmute(hospitalization_id, discord = pbw / pfvc, pfvc = pfvc,
                           Y = as.integer(!is.na(death_day) & death_day <= 28)) %>%
  filter(!is.na(discord), !is.na(pfvc))
wide <- long %>%
  select(hospitalization_id, vent_day, A, sf, map, on_pressor, exposed,
         age10, sex_category, race_category, height_grp, sofa_total) %>%
  pivot_wider(id_cols = c(hospitalization_id, age10, sex_category, race_category, height_grp, sofa_total),
              names_from = vent_day, values_from = c(A, sf, map, on_pressor, exposed), names_sep = "_") %>%
  inner_join(mods, by = "hospitalization_id") %>% as.data.frame()
cat(sprintf("11.Y: %d patients; PBW/PFVC discordance median %.2f, PFVC median %.2f L; ADDITIVE delta %.3f (J/min)/kg, floor %.2f, K %d, folds %d\n",
            nrow(wide), median(wide$discord), median(wide$pfvc), DELTA, FLOOR, K, FOLDS))

# age enters the nuisance models as a natural cubic spline, matching the engine's weight
# models and the titration: the exposures are deterministic in age through the GLI equations,
# so a linear age term leaves a curved residual correlated with any PFVC-derived modifier.
# lmtp takes column names, so the basis is materialized as columns.
age_basis <- splines::ns(wide$age10, 4); for (j in 1:4) wide[[paste0("age_ns", j)]] <- age_basis[, j]
base_cov <- c(paste0("age_ns", 1:4), "sex_category", "race_category", "height_grp", "sofa_total")
tv <- lapply(1:K, function(t) paste0(c("sf", "map", "on_pressor", "exposed"), "_", t))
trt <- paste0("A_", 1:K)
# ADDITIVE (normalizer-DEPENDENT) MTP: MP/PBW -> max(MP/PBW - DELTA, FLOOR), implemented on the
# log scale the treatment density learners see. A fixed PBW-normalized decrement = a larger true
# (MP/PFVC) cut in the misdosed, by DELTA x (PBW/PFVC). Floored so already-low days are untouched.
shift_fun <- function(data, t) { i <- sub("^A_", "", t); ex <- data[[paste0("exposed_", i)]]
  out <- data[[t]]; reduced <- log(pmax(exp(out) - DELTA, FLOOR))
  out[ex == 1L] <- reduced[ex == 1L]; out }
# realized intervention intensity (audit the feasible-MTP bite): mean cut + fraction floor-bound
.exp_days <- long %>% filter(exposed == 1L) %>% mutate(raw = exp(A))
cat(sprintf("    realized cut: median MP/PBW %.3f -> %.3f (J/min)/kg; mean reduction %.1f%%; %.1f%% of exposed days hit the floor\n",
            median(.exp_days$raw), median(pmax(.exp_days$raw - DELTA, FLOOR)),
            100 * mean(1 - pmax(.exp_days$raw - DELTA, FLOOR) / .exp_days$raw),
            100 * mean(.exp_days$raw - DELTA < FLOOR)))

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

# DR pseudo-outcomes (uncentered = estimate + EIF); ITE = shift - natural (mean = ATE/RD)
# DR pseudo-outcomes for EACH policy (uncentered = estimate + EIF). The difference is the
# per-patient ITE (absolute scale); keeping them separate also allows the RELATIVE scale below.
ps_dr <- fs$estimate@x + fs$estimate@eif      # pseudo-outcome under the shifted policy
pn_dr <- fo$estimate@x + fo$estimate@eif      # pseudo-outcome under the natural course
ite <- ps_dr - pn_dr
stopifnot(length(ite) == nrow(wide))
cat(sprintf("    ATE (MP/PBW -%.3f (J/min)/kg additive, PBW-anchored): risk_natural %.1f%% (CANARY ~27-32 real), RD %+.2f pp [%.2f, %.2f]\n",
            DELTA, 100 * ate$ref, 100 * ate$estimate, 100 * ate$conf.low, 100 * ate$conf.high))
if (abs(ate$conf.high - ate$conf.low) > 0.40)
  message("  WARNING: ATE CI width ", round(100 * (ate$conf.high - ate$conf.low)), "pp -- MP/PBW shift poorly identified; read the gradient, not the level.")
# ATE table (same columns as the titration's, so pooled_tte.R pools the two policies alike)
write_csv(tibble(policy = "mppbw_additive", delta = DELTA, floor = FLOOR, horizon_days = K,
                 n_patients = nrow(wide), risk_natural = ate$ref, risk_policy = ate$ref + ate$estimate,
                 rd = ate$estimate, rd_lo = ate$conf.low, rd_hi = ate$conf.high, rd_se = ate$std.error,
                 frac_exposed_days_shifted = mean(.exp_days$raw - DELTA > FLOOR), site = site_name),
          file.path(final_dir, paste0("mppbw_additive_ate_", site_name, ".csv")))
# winsorize EIF tails to tame the CATE SHAPE variance, THEN recenter to the unbiased TMLE ATE (the
# gradient/shape is a contrast and is unchanged by recentering; the level is the honest ATE).
n_wins <- sum(abs(ite) > WINSOR)
itew <- pmin(pmax(ite, -WINSOR), WINSOR)
itew <- itew - mean(itew) + ate$estimate
cat(sprintf("    winsorized %d/%d ITEs to +-%.1f, recentered CATE to TMLE ATE %.2f pp (gradient/shape unaffected)\n",
            n_wins, length(ite), WINSOR, 100 * ate$estimate))

# --- CATE by each modifier: spline of ITE on log(modifier); slope on a linear fit -------
# discordance: signal expected as a DOWN-slope (RD more negative where PBW oversizes).
# pfvc:        signal expected as an UP-slope (RD more negative at smaller absolute lung).
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
  # CATE gradient = CATE(p90) - CATE(p10); CI by patient bootstrap on (modifier, pseudo-outcome)
  # with FIXED knots (the regression's sampling variability, conditional on the DR pseudo-outcomes
  # -- same approximation as the per-log slope CI; does NOT propagate first-stage LMTP uncertainty).
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
cate_tbl <- bind_rows(ct_disc$tbl, ct_pfvc$tbl) %>% mutate(policy = "mppbw_additive", site = site_name)
write_csv(cate_tbl, file.path(final_dir, paste0("mppbw_additive_cate_", site_name, ".csv")))
slope_tbl <- bind_rows(ct_disc$slope, ct_pfvc$slope) %>% mutate(policy = "mppbw_additive", site = site_name)
write_csv(slope_tbl, file.path(final_dir, paste0("mppbw_additive_slope_", site_name, ".csv")))

# --- figure: one panel per modifier, ATE reference line --------------------------------
ate_pp <- 100 * ate$estimate   # official TMLE ATE as the figure reference (not the winsorized mean)
POLICY_LAB <- sprintf("MP/PBW -%.3f (J/min)/kg", DELTA)
mk_panel <- function(ct, xlab) {
  ggplot(ct$fig, aes(x, cate)) +
    geom_hline(yintercept = ate_pp, linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, fill = "#0072B2") +
    geom_line(colour = "#0072B2", linewidth = 1) +
    geom_rug(data = ct$rug %>% slice_sample(n = min(nrow(ct$rug), 2000)), aes(x = x),
             inherit.aes = FALSE, alpha = 0.06, length = unit(0.03, "npc")) +
    labs(x = xlab, y = sprintf("CATE: 28-d mortality RD of MP/PBW -%.3f (pp)", DELTA), title = ct$fig$modifier[1]) +
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
write_csv(rr_tbl, file.path(final_dir, paste0("mppbw_additive_cate_rr_", site_name, ".csv")))
# 2 x 2: absolute risk difference (top) and relative risk ratio (bottom) for each modifier,
# with a rug of the modifier's distribution under every panel. Read the two rows together: an
# absolute curve that bends where the relative curve is flat is baseline risk, not targeting.
fig <- (mk_panel(ct_disc, "PBW/PFVC discordance - higher = PBW oversizes (PFVC says smaller)") | mk_panel(ct_pfvc, "PFVC (L) - lower = smaller predicted lung")) /
       (mk_panel_rr(ct_disc, "PBW/PFVC discordance - higher = PBW oversizes (PFVC says smaller)") | mk_panel_rr(ct_pfvc, "PFVC (L) - lower = smaller predicted lung")) +
  patchwork::plot_annotation(
    title = paste0("CATE of an ADDITIVE MP/PBW power-reduction MTP, by PFVC-derived modifiers - ", site_name, if (is_synthetic) " (SYNTHETIC)" else ""),
    subtitle = paste0("Top: absolute CATE (risk difference); dashed = ATE; negative = benefit. ",
                      "Bottom: relative CATE (risk ratio, policy / natural course); 1 = no effect. ",
                      "Rug = the modifier's distribution. A sloped absolute curve with a flat relative curve is ",
                      "baseline risk (the modifiers mark the highest-risk patients), not effect modification."))
ggsave(file.path(final_dir, paste0("mppbw_additive_cate_", site_name, ".pdf")), fig, width = 12, height = 9)

cat("\n=== 11.Y continuous CATE of an ADDITIVE MP/PBW reduction, by PFVC-derived modifiers ===\n")
cat("--- effect-modification trends (read the SHAPE, not the level: confounding-by-severity inflates the ATE) ---\n")
print(as.data.frame(slope_tbl), row.names = FALSE)
cat("    (discordance grad_p90_minus_p10 MORE NEGATIVE => the misdosed gain more from a fixed PBW-anchored power cut;\n")
cat("     per_log slope CI = analytic regression; gradient CI = patient bootstrap -- both conditional on the DR pseudo-outcomes)\n")
print(as.data.frame(cate_tbl), row.names = FALSE)
message("Wrote mppbw_additive_{cate,slope}_", site_name, ".csv + cate .pdf to ", final_dir)
