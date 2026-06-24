# =============================================================================
# Script 11.X: Does the strain-limiting (PFVC-guided) policy benefit the MISDOSED most?
#              Discordance-HTE + positivity-by-discordance + dose-correction decomposition
# =============================================================================
# The causal payoff of the whole TTE for the PBW-vs-PFVC thesis. The primary contrast
# (strain_limiting VT/PFVC <= C_LOW vs permissive <= C_HIGH ~ usual PBW-care) is a
# PFVC-GUIDED dosing policy vs standard care, and -- unlike the 12.x multiplicative MP
# shift, where the normalizer cancels -- the CEILING estimand is normalizer-DEPENDENT.
# So the targeting question 12.E-H could not answer is answerable HERE: does holding a
# PFVC strain ceiling help most in the patients PBW oversizes (high PBW/PFVC discordance =
# the "misdosed")? This leaf makes that case in four reads (A continuous + A-D by tertile):
#
#   A. HETEROGENEITY (the claim). 28-d mortality RD within Concordant/Mid/Discordant
#      tertiles, PLUS the formal Discordant - Concordant GRADIENT with a bootstrap CI --
#      the actual "benefits the misdosed MORE" statistic. 11.A bootstraps each subgroup RD
#      but never their contrast; that contrast is the interaction test, computed here. From
#      ONE pooled, SOFA-adjusted, IPC-weighted MSM (arm x disc_grp interaction on a shared
#      day-spline), standardized per tertile -- stable CIs + corrects the residual severity
#      imbalance read D exposed. See the §A block header for the full rationale.
#
#   B. POSITIVITY (the honest caveat). The misdosed tertile is exactly where the policy
#      diverges most from observed care, so it is also where overlap is THINNEST. 11.O
#      runs this by PFVC and age; here it is indexed by the discordance tertile that
#      indexes read A, so the HTE and its support are read on the same axis. Empirical
#      P(adhere), weighted at-risk ESS, and deviation/trim per tertile x arm: if the
#      Discordant strain arm's ESS collapses, the Discordant RD is extrapolation and is
#      reported as such, NOT as a clean targeting win.
#
#   C. MECHANISM (why it works -- mechanical, not biological). The benefit need not come
#      from the discordant lung RESPONDING differently; it comes from the algorithm
#      DELIVERING A LARGER CORRECTION there. Per tertile: observed VT/PBW (mL/kg -- these
#      patients look LTVV-compliant on the PBW yardstick), observed VT/PFVC (% -- high
#      strain the PBW yardstick hides), the share of vent-days above the C_LOW ceiling,
#      and the VT cut the ceiling implies (dVT = observed VT - ceiling VT, in mL and
#      mL/kg PBW). A larger dVT in the Discordant tertile IS the mechanical engine of read A.
#
# Reuses the PRIMARY design (des/long_all from 10_tte_common via the engine loader); the
# only new fit is a focused cluster bootstrap for the per-tertile RDs + the gradient.
# Synthetic survival is simulated => plumbing only. Registered at the end of 11_run_all.
# =============================================================================
library(here); library(patchwork)
source(here::here("code", "10_tte_engine.R"))

DISC_LEVELS <- c("Concordant", "Mid", "Discordant")   # ascending PBW/PFVC; "Discordant" = misdosed

# =============================================================================
# A. Discordance HTE: ONE pooled, SOFA-adjusted, IPC-weighted MSM -> standardized per-tertile
#    28-d mortality RDs + the Discordant-Concordant gradient (stable bootstrap)
# =============================================================================
# REBUILD of the per-subset estimator. The old version refit `arm * ns(day,4)` separately
# WITHIN each discordance tertile; on resampled single tertiles that glm hit separation /
# extreme fitted risk, blowing up the MIMIC bootstrap CIs (Mid upper bound +24 pp). Two fixes
# in one model:
#   (1) STABILITY -- fit ONE pooled weighted MSM on the full clone panel with an arm x disc_grp
#       interaction (the HTE) on a SHARED day-spline. The day-dynamics are estimated from the
#       whole sample, so no single-tertile cell can separate; per-tertile RDs come out by
#       g-standardization, not subset refits.
#   (2) SOFA CORRECTION -- the 11.X balance read showed the IPCW leaves a residual baseline-
#       SOFA gap (strain arm ~0.10-0.15 SD LESS sick, NOT removed by the weights; the strongest
#       in the Discordant stratum). Adding baseline sofa_total to the OUTCOME model and
#       standardizing over each tertile's SOFA distribution mops up that severity tilt -- a
#       doubly-robust touch: weights handle the (already-clean) size axis, the outcome model
#       handles the residual severity. The Discordant RD is now SOFA-adjusted; compare it to the
#       unadjusted -8.7/-12.6 to read how much of the benefit was severity selection.
# RD per tertile = standardized E[Y^strain] - E[Y^permissive] over that tertile's empirical
# (disc_grp, sofa_total) profiles (SOFA is integer -> a handful of cells -> exact + fast).
long_s   <- long_all %>% left_join(base %>% select(hospitalization_id, sofa_total),
                                   by = "hospitalization_id")
prof_all <- long_s %>% distinct(hospitalization_id, disc_grp, sofa_total)
FORM     <- died ~ arm * disc_grp + arm * ns(day, 4) + disc_grp * ns(day, 4) + sofa_total
n_disc   <- prof_all %>% count(disc_grp) %>% mutate(disc_grp = as.character(disc_grp))

# Standardized RD in tertile g (g = "All" => marginal over the whole cohort), from a fitted
# pooled model: average the daily hazard over the target population's profiles per arm,
# cumulate to HORIZON, contrast arms. prof = patient-level (disc_grp, sofa_total) rows (with
# bootstrap multiplicity); count() collapses to unique cells weighted by frequency.
arm_f  <- function() factor(c("permissive", "strain_limiting"), c("permissive", "strain_limiting"))
std_rd <- function(fit, prof, g) {
  pp <- if (g == "All") prof else prof %>% filter(as.character(disc_grp) == g)
  if (nrow(pp) < 50) return(c(rd = NA_real_, risk_sl = NA_real_, risk_pm = NA_real_))
  cells <- pp %>% count(disc_grp, sofa_total, name = "wt")
  grid  <- tidyr::crossing(cells, day = 1:HORIZON, arm = arm_f())
  grid$haz <- predict(fit, grid, type = "response")
  ci <- grid %>% group_by(arm, day) %>% summarise(h = weighted.mean(haz, wt), .groups = "drop") %>%
    group_by(arm) %>% arrange(day) %>% summarise(ci = 1 - prod(1 - h), .groups = "drop")
  rs <- ci$ci[ci$arm == "strain_limiting"]; rp <- ci$ci[ci$arm == "permissive"]
  c(rd = rs - rp, risk_sl = rs, risk_pm = rp)
}

fit0 <- suppressWarnings(glm(FORM, data = long_s, family = binomial, weights = ipcw))
pt_l <- setNames(lapply(DISC_LEVELS, function(g) std_rd(fit0, prof_all, g)), DISC_LEVELS)
pt   <- vapply(pt_l, `[`, numeric(1), "rd")
overall_rd <- unname(std_rd(fit0, prof_all, "All")["rd"])
grad_pt    <- unname(pt["Discordant"] - pt["Concordant"])   # negative = helps the misdosed MORE

# stable bootstrap: resample patients, refit the POOLED model, re-standardize (no subset refit)
boot_one <- function() {
  samp   <- tibble(hospitalization_id = sample(ids, replace = TRUE))
  long_b <- long_s  %>% inner_join(samp, by = "hospitalization_id", relationship = "many-to-many")
  prof_b <- samp %>% left_join(prof_all, by = "hospitalization_id")   # one row per draw (keeps multiplicity)
  fit_b  <- suppressWarnings(tryCatch(glm(FORM, data = long_b, family = binomial, weights = ipcw),
                                      error = function(e) NULL))
  if (is.null(fit_b)) return(setNames(rep(NA_real_, 4), c(DISC_LEVELS, "gradient")))
  r <- vapply(DISC_LEVELS, function(g) unname(std_rd(fit_b, prof_b, g)["rd"]), numeric(1))
  c(r, gradient = unname(r["Discordant"] - r["Concordant"]))
}
n_cores_used <- min(N_CORES, N_BOOT)
message("11.X discordance-HTE bootstrap (", N_BOOT, " reps across ", n_cores_used,
        " core(s); pooled SOFA-adjusted MSM, standardized) ...")
boot_t0 <- Sys.time()
chunks   <- split(seq_len(N_BOOT), cut(seq_len(N_BOOT), min(20L, N_BOOT), labels = FALSE))
bts_list <- vector("list", N_BOOT); done <- 0L
if (n_cores_used > 1) {
  cl <- makeCluster(n_cores_used, type = "PSOCK")
  clusterEvalQ(cl, { library(tidyverse); library(splines) })
  clusterExport(cl, envir = .GlobalEnv, varlist = c(
    "ids", "long_s", "prof_all", "FORM", "std_rd", "arm_f", "DISC_LEVELS", "HORIZON", "boot_one"))
  clusterSetRNGStream(cl, 20260624)
  tryCatch(
    for (ch in chunks) {
      bts_list[ch] <- parLapply(cl, ch, function(bb) boot_one())
      done <- done + length(ch)
      message(sprintf("  bootstrap %d/%d (%2d%%) | elapsed %4.0fs", done, N_BOOT,
                      as.integer(round(100 * done / N_BOOT)),
                      as.numeric(difftime(Sys.time(), boot_t0, units = "secs"))))
    }, finally = stopCluster(cl))
} else {
  for (ch in chunks) {
    for (b in ch) bts_list[[b]] <- boot_one()
    done <- done + length(ch)
    message(sprintf("  bootstrap %d/%d", done, N_BOOT))
  }
}
bts <- do.call(rbind, bts_list)
ci  <- function(col) quantile(bts[, col], c(.025, .975), na.rm = TRUE)

hte <- tibble(disc_grp = DISC_LEVELS, rd = pt[DISC_LEVELS],
              rd_lo = vapply(DISC_LEVELS, function(k) ci(k)[1], numeric(1)),
              rd_hi = vapply(DISC_LEVELS, function(k) ci(k)[2], numeric(1)),
              risk_sl = vapply(pt_l, `[`, numeric(1), "risk_sl"),
              risk_pm = vapply(pt_l, `[`, numeric(1), "risk_pm"),
              adjustment = "sofa_baseline") %>%
  left_join(n_disc, by = "disc_grp") %>%
  mutate(disc_grp = factor(disc_grp, DISC_LEVELS))
gradient <- tibble(statistic = "RD(Discordant) - RD(Concordant)", adjustment = "sofa_baseline",
                   estimate = grad_pt, lo = ci("gradient")[1], hi = ci("gradient")[2])
write_csv(hte, file.path(final_dir, paste0("tte_ccw_disc_hte_", site_name, ".csv")))
write_csv(gradient, file.path(final_dir, paste0("tte_ccw_disc_gradient_", site_name, ".csv")))

# --- discordance-STRATIFIED E-value (11.A's E-value is overall only) -------------------
# The Discordant RD rests on the heaviest deviation/trim censoring AND the thinnest overlap,
# so it is the MOST vulnerable to unmeasured confounding of the adherence-censoring -- and the
# overall E-value understates that. Per tertile: the minimum risk-ratio-scale association an
# unmeasured time-varying confounder would need with BOTH deviation and mortality to explain
# the stratum RD away (point + the CI bound nearest the null). Same approximate-bound caveat as
# 11.A (textbook E-value is for a point exposure; here it indexes censoring-confounding). If the
# Discordant E-value is NOT larger than the Concordant, the EXTRA benefit in the misdosed would
# need an unmeasured confounder concentrated there to be explained away.
eval_fn <- function(rr) { rr <- if (rr >= 1) rr else 1 / rr; rr + sqrt(rr * (rr - 1)) }
disc_eval <- hte %>% transmute(disc_grp, risk_sl, risk_pm, rd, rd_lo, rd_hi) %>% rowwise() %>%
  mutate(rr_point = risk_sl / risk_pm,
         rr_lo = (risk_pm + rd_lo) / risk_pm, rr_hi = (risk_pm + rd_hi) / risk_pm,
         rr_ci_near = if (rr_point >= 1) min(rr_lo, rr_hi) else max(rr_lo, rr_hi),
         rr_ci_near = if ((rr_point >= 1) != (rr_ci_near >= 1)) 1 else rr_ci_near,  # CI crosses null
         evalue_point = eval_fn(rr_point), evalue_ci = eval_fn(rr_ci_near)) %>%
  ungroup() %>% select(disc_grp, rr_point, evalue_point, rr_ci_bound = rr_ci_near, evalue_ci)
write_csv(disc_eval, file.path(final_dir, paste0("tte_ccw_disc_evalue_", site_name, ".csv")))

# =============================================================================
# A2. CONTINUOUS CATE(discordance): spline effect-modification (threshold-capable) -- PRIMARY
# =============================================================================
# The tertile gradient is a 2-point summary of a continuous relationship, and the MIMIC
# Mid-null hinted the effect is NON-monotone -- plausibly a THRESHOLD (strain-limiting matters
# only once PBW/PFVC discordance is extreme). This is the primary effect-modification read: the
# SAME pooled, SOFA-adjusted, IPC-weighted MSM, but discordance enters CONTINUOUSLY as a natural
# spline interacted with arm (so the CATE can bend), standardized along a discordance grid -> a
# smooth CATE curve with a bootstrap band. Flat-then-falling = threshold; straight tilt =
# linear gradient. Parallels 12.F/12.H's continuous CATE for the MP shift (FLAT -- normalizer
# cancels); a SLOPE here = the ceiling policy's normalizer-dependence on the SAME axis. SOFA is
# standardized over the marginal (full-cohort) distribution at every discordance (a partial-
# effect curve), so the tertile points (population estimands, local SOFA) can sit slightly off.
# Fixed spline basis (knots from full data) so the bootstrap band reflects coefficient
# uncertainty only, not knot wobble. KD = ns df (>=1 interior knot -> threshold-capable).
KD      <- as.integer(Sys.getenv("PBWPFVC_DISC_CATE_DF", "3"))
disc_c  <- base %>% transmute(hospitalization_id, ldisc = log(pbw / pfvc))
long_c  <- long_s %>% left_join(disc_c, by = "hospitalization_id")
prof_c  <- long_c %>% distinct(hospitalization_id, ldisc, sofa_total)
Bspl    <- ns(prof_c$ldisc, df = KD)                     # FIXED basis -> stable band
zc      <- paste0("z", seq_len(KD))
add_z   <- function(d) { m <- predict(Bspl, d$ldisc); for (j in seq_len(KD)) d[[zc[j]]] <- m[, j]; d }
long_cz <- add_z(long_c)
FORM_C  <- as.formula(paste0("died ~ arm*(", paste(zc, collapse = "+"), ") + arm*ns(day,4) + (",
                             paste(zc, collapse = "+"), ")*ns(day,4) + sofa_total"))
# fixed eval grid (2.5-97.5th pct of log-discordance) + p10/p90 threshold-contrast anchors
DGRID   <- seq(quantile(prof_c$ldisc, .025), quantile(prof_c$ldisc, .975), length.out = 40)
qd      <- quantile(prof_c$ldisc, c(.10, .90))
EVAL    <- sort(unique(c(DGRID, qd)))
sofa_cells <- prof_c %>% count(sofa_total, name = "wt")  # marginal SOFA reference (fixed)
grid_c  <- tidyr::crossing(ldisc = EVAL, sofa_cells, day = 1:HORIZON, arm = arm_f()) %>% add_z()
cate_from_fit <- function(fit) {
  g <- grid_c; g$haz <- predict(fit, g, type = "response")
  g %>% group_by(ldisc, arm, day) %>% summarise(h = weighted.mean(haz, wt), .groups = "drop") %>%
    group_by(ldisc, arm) %>% arrange(day) %>% summarise(cif = 1 - prod(1 - h), .groups = "drop") %>%
    group_by(ldisc) %>%
    summarise(rd = cif[arm == "strain_limiting"] - cif[arm == "permissive"], .groups = "drop") %>%
    arrange(ldisc)
}
fitc0 <- suppressWarnings(glm(FORM_C, data = long_cz, family = binomial, weights = ipcw))
cate0 <- cate_from_fit(fitc0)                            # rows ordered by EVAL
i10 <- which.min(abs(cate0$ldisc - qd[1])); i90 <- which.min(abs(cate0$ldisc - qd[2]))
slope_pt <- cate0$rd[i90] - cate0$rd[i10]                # RD(p90 disc) - RD(p10 disc)

boot_c <- function() {
  samp <- tibble(hospitalization_id = sample(ids, replace = TRUE))
  lb   <- long_cz %>% inner_join(samp, by = "hospitalization_id", relationship = "many-to-many")
  fb   <- suppressWarnings(tryCatch(glm(FORM_C, data = lb, family = binomial, weights = ipcw),
                                    error = function(e) NULL))
  if (is.null(fb)) return(rep(NA_real_, length(EVAL)))
  # a rank-deficient resample can make predict/standardize throw or drop rows; never let one
  # replicate halt the parallel run -- degrade it to NA so the pointwise band uses na.rm.
  r <- tryCatch(cate_from_fit(fb)$rd, error = function(e) rep(NA_real_, length(EVAL)))
  if (length(r) != length(EVAL)) rep(NA_real_, length(EVAL)) else r
}
message("11.X continuous CATE bootstrap (", N_BOOT, " reps across ", n_cores_used, " core(s)) ...")
ct_list <- vector("list", N_BOOT); done <- 0L
if (n_cores_used > 1) {
  cl <- makeCluster(n_cores_used, type = "PSOCK")
  clusterEvalQ(cl, { library(tidyverse); library(splines) })
  clusterExport(cl, envir = .GlobalEnv, varlist = c("ids", "long_cz", "grid_c", "FORM_C",
    "cate_from_fit", "Bspl", "add_z", "zc", "boot_c", "HORIZON", "EVAL"))
  clusterSetRNGStream(cl, 20260625)
  tryCatch(for (ch in chunks) {
    ct_list[ch] <- parLapply(cl, ch, function(bb) boot_c())
    done <- done + length(ch); message(sprintf("  cate boot %d/%d", done, N_BOOT))
  }, finally = stopCluster(cl))
} else for (ch in chunks) {
  for (b in ch) ct_list[[b]] <- boot_c()
  done <- done + length(ch); message(sprintf("  cate boot %d/%d", done, N_BOOT))
}
ctm <- do.call(rbind, ct_list)                           # reps x length(EVAL)
n_bad <- sum(apply(ctm, 1, function(r) all(is.na(r))))   # whole-rep failures (rank-deficient resamples)
if (n_bad > 0) message(sprintf("  continuous CATE: %d/%d bootstrap reps failed -> dropped via na.rm%s",
        n_bad, N_BOOT, if (n_bad > N_BOOT / 2) " (WARNING: >half failed -- band unreliable)" else ""))
qlo <- apply(ctm, 2, quantile, .025, na.rm = TRUE); qhi <- apply(ctm, 2, quantile, .975, na.rm = TRUE)
slope_ci <- quantile(ctm[, i90] - ctm[, i10], c(.025, .975), na.rm = TRUE)

curve_tbl <- tibble(ldisc = cate0$ldisc, discordance = exp(cate0$ldisc),
                    rd = cate0$rd, rd_lo = qlo, rd_hi = qhi) %>% filter(ldisc %in% DGRID)
write_csv(curve_tbl, file.path(final_dir, paste0("tte_ccw_disc_cate_curve_", site_name, ".csv")))
slope_tbl <- tibble(statistic = "RD(p90 discordance) - RD(p10 discordance)", adjustment = "sofa_baseline",
                    disc_p10 = exp(qd[1]), disc_p90 = exp(qd[2]),
                    rd_p10 = cate0$rd[i10], rd_p90 = cate0$rd[i90],
                    estimate = slope_pt, lo = slope_ci[1], hi = slope_ci[2])
write_csv(slope_tbl, file.path(final_dir, paste0("tte_ccw_disc_cate_slope_", site_name, ".csv")))

# =============================================================================
# B. Positivity by discordance tertile (the support behind read A) -- mirrors 11.O on disc_grp
# =============================================================================
disc_key <- base %>% transmute(hospitalization_id, disc = pbw / pfvc, disc_grp)

# (1) empirical P(adhere) on eligible days (des carries pday: keep_pday=TRUE for primary)
pday <- bind_rows(des$bl$pday %>% mutate(arm = "strain_limiting"),
                  des$bh$pday %>% mutate(arm = "permissive"))   # pday already carries disc_grp
emp1 <- function(g) summarise(g, n_eligible_days = n(),
  frac_padhere_lt05 = mean(p_adhere < 0.05), frac_padhere_lt02 = mean(p_adhere < 0.02),
  p01_padhere = quantile(p_adhere, 0.01), median_padhere = median(p_adhere), .groups = "drop")
disc_empirical <- bind_rows(
  pday %>% group_by(arm, disc_grp) %>% emp1() %>% mutate(disc_grp = as.character(disc_grp)),
  pday %>% group_by(arm) %>% emp1() %>% mutate(disc_grp = "All")) %>% arrange(arm, disc_grp)
write_csv(disc_empirical, file.path(final_dir, paste0("tte_ccw_disc_overlap_empirical_", site_name, ".csv")))

# (2) weighted at-risk ESS by tertile x arm over follow-up
support_days <- c(7L, 14L, 28L)
support_for_arm <- function(b, arm_lab) {
  bb <- b$idsum   # idsum already carries disc_grp (first(disc_grp) in the engine)
  map_dfr(support_days, function(d) {
    w_at_d <- b$wday %>% filter(vent_day <= d) %>% group_by(hospitalization_id) %>%
      arrange(vent_day) %>% summarise(w = last(cumw), .groups = "drop")
    bb %>% filter(dev_day > d, trim_day > d, is.na(death_day) | death_day >= d) %>%
      left_join(w_at_d, by = "hospitalization_id") %>% mutate(w = trunc_w(coalesce(w, 1))) %>%
      group_by(disc_grp) %>%
      summarise(arm = arm_lab, day = d, n_atrisk = n(),
                ess_atrisk = if (n() > 0) (sum(w)^2 / sum(w^2)) else 0, .groups = "drop")
  })
}
disc_support <- bind_rows(support_for_arm(des$bl, "strain_limiting"),
                          support_for_arm(des$bh, "permissive")) %>% arrange(arm, day, disc_grp)
write_csv(disc_support, file.path(final_dir, paste0("tte_ccw_disc_overlap_support_", site_name, ".csv")))

# (3) deviation / trim by tertile x arm -- the mechanism behind any support gap
disc_loss <- bind_rows(des$bl$idsum %>% mutate(arm = "strain_limiting"),
                       des$bh$idsum %>% mutate(arm = "permissive")) %>%
  left_join(disc_key %>% select(hospitalization_id, disc), by = "hospitalization_id") %>%  # disc_grp already in idsum
  group_by(arm, disc_grp) %>%
  summarise(n = n(), median_disc = median(disc),
            frac_deviated = mean(is.finite(dev_day)), frac_trimmed = mean(is.finite(trim_day)),
            mean_stab_w = mean(ipcw_term), .groups = "drop") %>% arrange(arm, disc_grp)
write_csv(disc_loss, file.path(final_dir, paste0("tte_ccw_disc_overlap_loss_", site_name, ".csv")))

# =============================================================================
# C. Dose-correction decomposition: how big a VT cut the C_LOW ceiling delivers, by tertile
# =============================================================================
# Read the waterfall directly (raw `wf` is cache-blocklisted, so not guaranteed in memory).
# pfvc here is the ceiling normalizer carried in base (pfvc or pfvc_age25 per PBWPFVC_TTE_NORM),
# so the correction is computed against the SAME normalizer the ceiling is defined on.
dose_day <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet")) %>%
  select(hospitalization_id, recorded_dttm, tidal_volume_set) %>%
  filter(!is.na(tidal_volume_set), tidal_volume_set > 0) %>%
  inner_join(base %>% select(hospitalization_id, t0, pfvc, pbw, disc_grp), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  group_by(hospitalization_id, vent_day, pfvc, pbw, disc_grp) %>%
  summarise(vt = median(tidal_volume_set), .groups = "drop") %>%
  mutate(vtpbw  = vt / pbw,                 # mL/kg PBW (LTVV yardstick: 6-8 = "compliant")
         vtpfvc = vt / pfvc * 0.1,          # % predicted FVC (the strain the PBW yardstick hides)
         vt_ceiling = C_LOW * pfvc * 10,    # VT (mL) that exactly meets the C_LOW strain ceiling
         above  = vtpfvc > C_LOW,           # day the PFVC ceiling would require a cut
         dvt_ml   = pmax(vt - vt_ceiling, 0),
         dvt_mlkg = dvt_ml / pbw)           # the VT cut, per kg PBW
dose_tbl <- dose_day %>% group_by(disc_grp) %>%
  summarise(n_patients = n_distinct(hospitalization_id), n_pdays = n(),
            median_vtpbw  = median(vtpbw),                       # ~6-8 across tertiles = looks LTVV-OK
            median_vtpfvc = median(vtpfvc),                      # rises with discordance = hidden strain
            frac_days_above_clow = mean(above),                  # the misdosed exceed the ceiling more
            dvt_ml_above   = median(dvt_ml[above]),              # correction the algorithm delivers...
            dvt_mlkg_above = median(dvt_mlkg[above]),            # ...larger in the Discordant tertile
            .groups = "drop") %>%
  mutate(disc_grp = factor(disc_grp, DISC_LEVELS)) %>% arrange(disc_grp)
write_csv(dose_tbl, file.path(final_dir, paste0("tte_ccw_disc_dose_correction_", site_name, ".csv")))

# =============================================================================
# D. Within-tertile weighting balance: does the IPCW balance the arms INSIDE each
#    discordance stratum (where the per-tertile RDs, esp. the Discordant one, live)?
# =============================================================================
# The primary IPC weights are fit ONCE on the full cohort; 11.X then estimates RDs within
# tertiles using those weights. So the subgroup-relevant balance question is whether those
# full-cohort weights still balance the arms WITHIN each tertile. Metric = the 11.U
# convention: weighted strain - permissive mean difference / full-cohort SD (survivorship
# cancels in the difference -> the pure between-arm confounding gap), reported UNWEIGHTED vs
# WEIGHTED so the weights' work is visible. |SMD| < 0.1 = well balanced. NOTE: standardized
# by the FULL-cohort SD (as in 11.U, for cross-comparability); within the Discordant tertile
# PFVC's range is compressed, so its weighted SMD here is if anything CONSERVATIVE (small).
# "All" reproduces 11.U's strain_minus_permissive row as a cross-check.
BAL_COVS <- c("pfvc", "sofa_total", "age10")   # size axis (the known +SMD concern) + severity + age
BAL_DAYS <- c(2L, 7L, 14L, 28L)                # day 2 = first deviation (cleanest); 7-28 = arms diverge
ref_bal  <- base %>% summarise(across(all_of(BAL_COVS),
              list(m = ~mean(., na.rm = TRUE), s = ~sd(., na.rm = TRUE)), .names = "{.col}__{.fn}"))
# join only covariates NOT already on idsum (age10 is carried by idsum via first(age10);
# joining it again would create age10.x/.y and break the lookup). The rest come from idsum.
covdat   <- base %>% select(hospitalization_id, any_of(setdiff(BAL_COVS, names(des$bl$idsum))))
# weighted + unweighted arm means in tertile g (g = "All" => no tertile filter) at day d,
# on the at-risk (uncensored, alive) set -- exactly the population the day-d MSM hazard uses.
arm_means <- function(b, g, d) {
  wd <- b$wday %>% filter(vent_day <= d) %>% group_by(hospitalization_id) %>%
    arrange(vent_day) %>% summarise(w = last(cumw), .groups = "drop")
  ar <- b$idsum %>%
    filter(g == "All" | as.character(disc_grp) == g,
           dev_day > d, trim_day > d, is.na(death_day) | death_day >= d) %>%
    left_join(covdat, by = "hospitalization_id") %>%
    left_join(wd, by = "hospitalization_id") %>% mutate(w = trunc_w(coalesce(w, 1)))
  list(n = nrow(ar),
       mu = vapply(BAL_COVS, function(cv) mean(ar[[cv]], na.rm = TRUE), numeric(1)),
       mw = vapply(BAL_COVS, function(cv) weighted.mean(ar[[cv]], ar$w, na.rm = TRUE), numeric(1)))
}
disc_balance <- map_dfr(c("All", DISC_LEVELS), function(g) map_dfr(BAL_DAYS, function(d) {
  s <- arm_means(des$bl, g, d); p <- arm_means(des$bh, g, d)
  map_dfr(BAL_COVS, function(cv) tibble(
    disc_grp = g, day = d, covariate = cv, n_strain = s$n, n_perm = p$n,
    smd_unweighted = unname((s$mu[cv] - p$mu[cv]) / ref_bal[[paste0(cv, "__s")]]),
    smd_weighted   = unname((s$mw[cv] - p$mw[cv]) / ref_bal[[paste0(cv, "__s")]])))
})) %>% mutate(disc_grp = factor(disc_grp, c("All", DISC_LEVELS))) %>%
  arrange(disc_grp, day, covariate)
write_csv(disc_balance, file.path(final_dir, paste0("tte_ccw_disc_balance_", site_name, ".csv")))
max_wsmd <- disc_balance %>% group_by(disc_grp) %>%
  summarise(worst_abs_weighted_smd = max(abs(smd_weighted)), .groups = "drop")

# --- D2. TIME-VARYING confounder balance within each discordance stratum --------------
# Read D checked BASELINE covariates; but the IPCW exists to balance the TIME-VARYING drivers
# of the deviation decision (concurrent strain VT/PFVC, oxygenation S/F, MAP, pressor use).
# Whether THOSE balance between arms within the Discordant stratum -- where censoring is
# heaviest -- is the direct test that the informative censoring is actually handled there. For
# each arm we take the at-risk (uncensored, alive) clones at day d, their concurrent panel
# confounders, weighted by the cumulative IPCW; the strain - permissive weighted SMD
# (standardized by the full-panel SD) is the residual imbalance the weights leave. |SMD| < 0.1
# in the Discordant arm-contrast = the deviation-censoring is well-corrected where it bites most.
TV_COVS <- c("vtpfvc", "sf", "map", "on_pressor")
tv_ref  <- panel %>% summarise(across(all_of(TV_COVS), ~ sd(., na.rm = TRUE), .names = "{.col}__s"))
tv_means <- function(b, d) {        # weighted at-risk concurrent-confounder means by stratum, one arm
  w_at_d <- b$wday %>% filter(vent_day <= d) %>% group_by(hospitalization_id) %>%
    arrange(vent_day) %>% summarise(w = last(cumw), .groups = "drop")
  b$idsum %>% filter(dev_day > d, trim_day > d, is.na(death_day) | death_day >= d) %>%
    select(hospitalization_id, disc_grp) %>%
    inner_join(panel %>% filter(vent_day == d) %>% select(hospitalization_id, all_of(TV_COVS)),
               by = "hospitalization_id") %>%
    left_join(w_at_d, by = "hospitalization_id") %>% mutate(w = trunc_w(coalesce(w, 1))) %>%
    group_by(disc_grp) %>%
    summarise(across(all_of(TV_COVS), ~ weighted.mean(., w, na.rm = TRUE)), n = n(), .groups = "drop")
}
tv_balance <- map_dfr(c(2L, 7L, 14L), function(d) {
  j <- inner_join(tv_means(des$bl, d), tv_means(des$bh, d), by = "disc_grp", suffix = c("_s", "_p"))
  out <- tibble(disc_grp = j$disc_grp, day = d, n_strain = j$n_s, n_perm = j$n_p)
  for (cv in TV_COVS)
    out[[paste0("smd_", cv)]] <- (j[[paste0(cv, "_s")]] - j[[paste0(cv, "_p")]]) / tv_ref[[paste0(cv, "__s")]]
  out
}) %>% mutate(disc_grp = factor(disc_grp, DISC_LEVELS)) %>% arrange(disc_grp, day)
write_csv(tv_balance, file.path(final_dir, paste0("tte_ccw_disc_tvbalance_", site_name, ".csv")))
tv_worst <- tv_balance %>% rowwise() %>%
  mutate(rmax = max(abs(c_across(starts_with("smd_"))), na.rm = TRUE)) %>% ungroup() %>%
  group_by(disc_grp) %>% summarise(worst_abs_tv_smd = max(rmax), .groups = "drop")

# =============================================================================
# E. Gradient robustness to common-support trim + day-weight cap (informative-censoring stress)
# =============================================================================
# Does the Discordant-Concordant gradient survive stricter/looser positivity handling? The
# Discordant RD leans on the thin-support, heavily-censored tail; if the gradient flips or
# collapses under a stricter trim, it is a tail artifact. POINT estimates only (rebuilds the
# design per setting; no per-setting bootstrap), via the same pooled SOFA-adjusted standardizer.
grad_for_design <- function(trim_a, cap_a) {
  dd <- build_design(C_LOW, C_HIGH, GRACE, cap_a, "simple", trim = trim_a)
  lj <- dd$long %>% left_join(base %>% select(hospitalization_id, sofa_total), by = "hospitalization_id")
  pj <- lj %>% distinct(hospitalization_id, disc_grp, sofa_total)
  f  <- suppressWarnings(glm(FORM, data = lj, family = binomial, weights = ipcw))
  r  <- vapply(DISC_LEVELS, function(g) unname(std_rd(f, pj, g)["rd"]), numeric(1))
  tibble(trim = trim_a, cap = cap_a, rd_concordant = r["Concordant"], rd_discordant = r["Discordant"],
         gradient = unname(r["Discordant"] - r["Concordant"]))
}
sweep_grid <- list(c(TRIM_ALPHA, DAYW_CAP),                       # primary (reference row)
                   c(0, DAYW_CAP), c(0.01, DAYW_CAP), c(0.05, DAYW_CAP),  # trim sweep
                   c(TRIM_ALPHA, 3), c(TRIM_ALPHA, 10))                   # weight-cap sweep
message("11.X gradient robustness: rebuilding ", length(sweep_grid), " designs (trim x cap) ...")
grad_sweep <- map_dfr(sweep_grid, function(s) grad_for_design(s[1], s[2])) %>%
  mutate(setting = if_else(trim == TRIM_ALPHA & cap == DAYW_CAP, "PRIMARY", "sensitivity"), .before = 1)
write_csv(grad_sweep, file.path(final_dir, paste0("tte_ccw_disc_gradient_sweep_", site_name, ".csv")))

# =============================================================================
# Figure: RD by discordance tertile (with gradient) + dose correction the ceiling delivers
# =============================================================================
# continuous CATE curve (primary) with the tertile RDs overlaid at each tertile's median discordance
disc_med <- base %>% mutate(dc = pbw / pfvc) %>% group_by(disc_grp) %>%
  summarise(disc_med = median(dc), .groups = "drop") %>% mutate(disc_grp = factor(disc_grp, DISC_LEVELS))
hte_ov <- hte %>% left_join(disc_med, by = "disc_grp")
pA <- ggplot(curve_tbl, aes(discordance, 100 * rd)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_hline(yintercept = 100 * overall_rd, linetype = "dashed", colour = "grey55") +
  geom_ribbon(aes(ymin = 100 * rd_lo, ymax = 100 * rd_hi), alpha = 0.15, fill = "#0072B2") +
  geom_line(colour = "#0072B2", linewidth = 1) +
  geom_pointrange(data = hte_ov, aes(x = disc_med, y = 100 * rd, ymin = 100 * rd_lo, ymax = 100 * rd_hi,
                                     colour = disc_grp), inherit.aes = FALSE) +
  scale_colour_manual(values = setNames(okabe[c(2, 1, 6)], DISC_LEVELS), name = "tertile") +
  labs(x = "PBW/PFVC discordance (higher = PBW oversizes -> misdosed)",
       y = "CATE: 28-d mortality RD (pp)\nstrain-limiting - permissive",
       title = "A. Continuous CATE by discordance (spline; SOFA-adjusted)",
       subtitle = sprintf("RD(p90)-RD(p10) discordance: %+.2f pp [%.2f, %.2f]  (curve flat-then-down = threshold; points = tertiles)",
                          100 * slope_pt, 100 * slope_ci[1], 100 * slope_ci[2])) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
pB <- ggplot(dose_tbl, aes(disc_grp, dvt_mlkg_above, fill = disc_grp)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = setNames(okabe[c(2, 1, 6)], DISC_LEVELS), guide = "none") +
  labs(x = NULL, y = "Median VT cut the ceiling delivers\n(mL/kg PBW, above-ceiling days)",
       title = "C. Dose correction delivered",
       subtitle = "Larger in the misdosed = the mechanical engine of A") +
  theme_minimal(base_size = 11)
fig <- pA + pB + patchwork::plot_annotation(
  title = paste0("Strain-limiting (PFVC-guided) dosing benefits the misdosed - ", site_name,
                 if (is_synthetic) " (SYNTHETIC - plumbing only)" else ""),
  subtitle = "Read with positivity-by-discordance (tte_ccw_disc_overlap_*): a collapsed Discordant strain-arm ESS = the RD there is extrapolation.")
ggsave(file.path(final_dir, paste0("tte_ccw_disc_benefit_", site_name, ".pdf")), fig, width = 11, height = 5)

# =============================================================================
# Console summary
# =============================================================================
cat("\n=== 11.X strain-limiting benefit by PBW/PFVC discordance (misdosed = Discordant; SOFA-adjusted, standardized) ===\n")
cat(sprintf("    Overall RD: %+.2f pp\n", 100 * overall_rd))
print(as.data.frame(hte %>% transmute(disc_grp, n,
        rd_pp = round(100 * rd, 2), ci = sprintf("[%.2f, %.2f]", 100 * rd_lo, 100 * rd_hi))),
      row.names = FALSE)
cat(sprintf("    GRADIENT RD(Discordant)-RD(Concordant): %+.2f pp [%.2f, %.2f]  (negative => helps misdosed more)\n",
            100 * grad_pt, 100 * gradient$lo, 100 * gradient$hi))
cat("\n--- Continuous CATE(discordance), SOFA-adjusted spline (PRIMARY effect-modification) ---\n")
cat(sprintf("    RD at low discordance (p10 = %.1f): %+.2f pp;   at high (p90 = %.1f): %+.2f pp\n",
            exp(qd[1]), 100 * cate0$rd[i10], exp(qd[2]), 100 * cate0$rd[i90]))
cat(sprintf("    Threshold contrast RD(p90)-RD(p10): %+.2f pp [%.2f, %.2f]  (flat low / steep high = threshold, not linear)\n",
            100 * slope_pt, 100 * slope_ci[1], 100 * slope_ci[2]))
cat("\n--- Positivity by discordance tertile (strain arm; watch the Discordant ESS) ---\n")
print(as.data.frame(disc_support %>% filter(arm == "strain_limiting", day == 28L) %>%
        transmute(disc_grp, day, n_atrisk, ess_atrisk = round(ess_atrisk, 1))), row.names = FALSE)
cat("\n--- Dose correction the C_LOW =", C_LOW, "% ceiling delivers, by tertile ---\n")
print(as.data.frame(dose_tbl %>% transmute(disc_grp, n_patients,
        median_vtpbw = round(median_vtpbw, 1), median_vtpfvc = round(median_vtpfvc, 1),
        pct_days_above = round(100 * frac_days_above_clow),
        dvt_mlkg_above = round(dvt_mlkg_above, 2))), row.names = FALSE)
cat("    (vtpbw ~6-8 across tertiles = all look LTVV-compliant; vtpfvc + dose-cut rise with discordance = the hidden strain PFVC corrects)\n")

cat("\n--- Within-tertile IPCW balance (weighted strain - permissive SMD; |SMD|<0.1 = balanced) ---\n")
print(as.data.frame(disc_balance %>% filter(day %in% c(2L, 28L)) %>%
        transmute(disc_grp, day, covariate, n_strain, n_perm,
                  smd_unw = round(smd_unweighted, 3), smd_wt = round(smd_weighted, 3))), row.names = FALSE)
cat("    Worst |weighted SMD| (baseline) by stratum:\n")
print(as.data.frame(max_wsmd %>% mutate(worst_abs_weighted_smd = round(worst_abs_weighted_smd, 3))), row.names = FALSE)
cat("    Worst |weighted SMD| (TIME-VARYING confounders the IPCW targets) by stratum:\n")
print(as.data.frame(tv_worst %>% mutate(worst_abs_tv_smd = round(worst_abs_tv_smd, 3))), row.names = FALSE)
cat("    (Discordant-stratum baseline AND time-varying SMD < 0.1 = informative censoring corrected where it bites most)\n")

cat("\n--- Discordance-stratified E-value (robustness to unmeasured censoring-confounding) ---\n")
print(as.data.frame(disc_eval %>% mutate(across(where(is.numeric), ~ round(., 2)))), row.names = FALSE)
cat("    (Discordant E-value >= Concordant => the EXTRA benefit needs an unmeasured confounder concentrated in the misdosed)\n")

cat("\n--- Gradient robustness to trim x weight-cap (point estimates; PRIMARY = trim", TRIM_ALPHA, "cap", DAYW_CAP, ") ---\n")
print(as.data.frame(grad_sweep %>% transmute(setting, trim, cap,
        rd_disc_pp = round(100 * rd_discordant, 2), rd_conc_pp = round(100 * rd_concordant, 2),
        gradient_pp = round(100 * gradient, 2))), row.names = FALSE)
cat("    (gradient stays negative across trim/cap = not a thin-support / weight-tail artifact)\n")
message("Wrote tte_ccw_disc_{hte,gradient,evalue,cate_curve,cate_slope,overlap_*,dose_correction,balance,tvbalance,gradient_sweep}_",
        site_name, ".csv + benefit .pdf to ", final_dir)
