# =============================================================================
# Script 11_vtpbw_titration: CATE of an ADDITIVE VT/PBW titration MTP, by discordance + PFVC
# =============================================================================
# The bedside-DIRECT titration sibling to 11.X. 11.X is a clone-censor-weight TTE of a VT/PFVC
# strain-limiting CEILING. This leaf studies the intervention a clinician actually performs on the
# yardstick they actually use: "turn the tidal volume down half a click" -- an additive
# MODIFIED TREATMENT POLICY, VT/PBW -> VT/PBW - DELTA (default 0.5 mL/kg PBW), estimated by
# doubly-robust LMTP, with CATE by two PFVC-derived modifiers (PBW/PFVC discordance; PFVC).
#
# WHY THIS POLICY, FOR THE TITRATION OBJECTIVE:
#  * IDENTIFIABLE where a VT/PBW CEILING is not. VT/PBW is protocolized to ~6-8 in the data (no
#    static dose contrast), but a -DELTA *shift* stays inside observed support (~5.5-7.5), so the
#    feasible-MTP positivity holds. This is the modified-treatment-policy positivity workaround.
#  * NORMALIZER-DEPENDENT geometry => it targets the misdosed for free. A fixed -DELTA mL/kg PBW
#    removes DELTA x (PBW/PFVC) = DELTA x discordance of true (PFVC-normalized) strain, so the
#    discordant get the bigger strain cut -- the same gradient mechanism as the 11.X ceiling.
#  * LESS severity-confounded than the MP sibling (11.Y): VT/PBW is protocolized, so it carries
#    little severity-driven variation, unlike MP (which carries RR + airway pressures).
#
# READING IT (titration is the objective; the PHYSIOLOGIC case is made separately in 05). Per 11.Z
# discordance is ~deterministic in demographics (R^2~0.99), so "the misdosed" = the demographic
# groups PBW over-doses -- the targeting is demographically PATTERNED (the equity reading), NOT a
# claim of orthogonal physiology. Confounding-by-severity inflates the ATE *level*, so read the
# discordance GRADIENT (which a uniform severity bias does not generate), not the level.
#
# Design = the validated 12.D/F LMTP estimator (binomial, mtp=TRUE; shift bites only on
# on-vent-and-alive days; 28-day all-cause mortality). CATE via the DR-learner (Kennedy): DR
# pseudo-outcome per patient, ITE = pseudo(shift) - pseudo(natural) (mean = ATE/RD), regress on a
# spline of the log-modifier -> CATE(modifier) + per-log slope + p90-p10 gradient (bootstrap CI).
# Standalone supplement -- NOT in 11_run_all. Synthetic mortality is simulated => plumbing only.
# Env: PBWPFVC_VT_{DELTA,FLOOR,HORIZON,FOLDS,LEARNERS,BOOT}.
# =============================================================================
library(here); library(lmtp); library(splines); library(patchwork)
# lmtp reports progress via progressr, OFF by default -> enable a handler so the per-fold/timepoint
# bar actually renders (PBWPFVC_PROGRESS=0 to silence). cli if available, else base txtprogressbar.
if (!identical(Sys.getenv("PBWPFVC_PROGRESS", "1"), "0") && requireNamespace("progressr", quietly = TRUE)) {
  progressr::handlers(global = TRUE)
  progressr::handlers(if (requireNamespace("cli", quietly = TRUE)) "cli" else "txtprogressbar")
}
source(here::here("code", "10_tte_engine.R"))

DELTA <- as.numeric(Sys.getenv("PBWPFVC_VT_DELTA", "0.5"))    # ADDITIVE cut, mL/kg PBW (half a click)
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
if (any(!.avail)) message("11_vtpbw: dropping unavailable learners: ", paste(LRN[!.avail], collapse = ", "))
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
                           Y = as.integer(!is.na(death_day) & death_day <= 28)) %>%
  filter(!is.na(discord), !is.na(pfvc))
wide <- long %>%
  select(hospitalization_id, vent_day, A, sf, map, on_pressor, exposed,
         age10, sex_category, race_category, height_grp, sofa_total) %>%
  pivot_wider(id_cols = c(hospitalization_id, age10, sex_category, race_category, height_grp, sofa_total),
              names_from = vent_day, values_from = c(A, sf, map, on_pressor, exposed), names_sep = "_") %>%
  inner_join(mods, by = "hospitalization_id") %>% as.data.frame()
cat(sprintf("11_vtpbw: %d patients; PBW/PFVC discordance median %.2f, PFVC median %.2f L; ADDITIVE delta %.2f mL/kg, floor %.1f, K %d, folds %d\n",
            nrow(wide), median(wide$discord), median(wide$pfvc), DELTA, FLOOR, K, FOLDS))

base_cov <- c("age10", "sex_category", "race_category", "height_grp", "sofa_total")
tv <- lapply(1:K, function(t) paste0(c("sf", "map", "on_pressor", "exposed"), "_", t))
trt <- paste0("A_", 1:K)
# ADDITIVE (normalizer-DEPENDENT) MTP: VT/PBW -> max(VT/PBW - DELTA, FLOOR), on the log scale the
# treatment density learners see. A fixed PBW-normalized decrement = a larger true (VT/PFVC) strain
# cut in the misdosed, by DELTA x (PBW/PFVC). Floored so already-low days are untouched.
shift_fun <- function(data, t) { i <- sub("^A_", "", t); ex <- data[[paste0("exposed_", i)]]
  out <- data[[t]]; reduced <- log(pmax(exp(out) - DELTA, FLOOR))
  out[ex == 1L] <- reduced[ex == 1L]; out }
# realized intervention intensity (audit the feasible-MTP bite): mean cut + fraction floor-bound
.exp_days <- long %>% filter(exposed == 1L) %>% mutate(raw = exp(A))
cat(sprintf("    realized cut: median VT/PBW %.2f -> %.2f mL/kg; mean reduction %.1f%%; %.1f%% of exposed days hit the floor\n",
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

ite <- (fs$estimate@x + fs$estimate@eif) - (fo$estimate@x + fo$estimate@eif)
stopifnot(length(ite) == nrow(wide))
cat(sprintf("    ATE (VT/PBW -%.2f mL/kg additive titration): risk_natural %.1f%% (CANARY ~27-32 real), RD %+.2f pp [%.2f, %.2f]\n",
            DELTA, 100 * ate$ref, 100 * ate$estimate, 100 * ate$conf.low, 100 * ate$conf.high))
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
cate_one <- function(mod, lab) {
  dd  <- tibble(lx = log(wide[[mod]]), ite = itew) %>% filter(is.finite(lx))
  kn  <- attr(ns(dd$lx, 3), "knots"); bd <- attr(ns(dd$lx, 3), "Boundary.knots")  # FIX basis
  fit_spl <- function(d) lm(ite ~ ns(lx, knots = kn, Boundary.knots = bd), data = d)
  spl <- fit_spl(dd)
  tr  <- summary(lm(ite ~ lx, data = dd))$coefficients["lx", ]   # per-log-unit linear slope
  g   <- tibble(lx = quantile(dd$lx, seq(0.1, 0.9, 0.1)))
  pr  <- predict(spl, g, se.fit = TRUE)
  tbl <- g %>% transmute(modifier = lab, pctile = seq(10, 90, 10), value = round(exp(lx), 2),
                         cate_pp = round(100 * pr$fit, 2),
                         lo = round(100 * (pr$fit - 1.96 * pr$se.fit), 2),
                         hi = round(100 * (pr$fit + 1.96 * pr$se.fit), 2))
  qg  <- quantile(dd$lx, c(0.1, 0.9))
  grad_pt <- diff(as.numeric(predict(spl, tibble(lx = qg))))
  bs  <- vapply(seq_len(BGRAD), function(b) {
    fb <- tryCatch(fit_spl(dd[sample.int(nrow(dd), replace = TRUE), ]), error = function(e) NULL)
    if (is.null(fb)) NA_real_ else diff(as.numeric(predict(fb, tibble(lx = qg))))
  }, numeric(1))
  bs <- bs[is.finite(bs)]
  slope <- tibble(modifier = lab,
    statistic   = c("per_log_unit_slope", "grad_p90_minus_p10"),
    estimate_pp = round(100 * c(tr["Estimate"], grad_pt), 2),
    lo = round(100 * c(tr["Estimate"] - 1.96 * tr["Std. Error"], quantile(bs, 0.025)), 2),
    hi = round(100 * c(tr["Estimate"] + 1.96 * tr["Std. Error"], quantile(bs, 0.975)), 2),
    p  = c(signif(tr["Pr(>|t|)"], 2), NA_real_), n_boot = c(NA_integer_, length(bs)))
  fg  <- tibble(lx = seq(quantile(dd$lx, 0.02), quantile(dd$lx, 0.98), length.out = 60))
  fp  <- predict(spl, fg, se.fit = TRUE)
  fig <- tibble(modifier = lab, x = exp(fg$lx), cate = 100 * fp$fit,
                lo = 100 * (fp$fit - 1.96 * fp$se.fit), hi = 100 * (fp$fit + 1.96 * fp$se.fit))
  list(tbl = tbl, trend = tr, fig = fig, slope = slope)
}
ct_disc <- cate_one("discord", "PBW/PFVC discordance")
ct_pfvc <- cate_one("pfvc",    "PFVC (predicted size, L)")
cate_tbl <- bind_rows(ct_disc$tbl, ct_pfvc$tbl)
write_csv(cate_tbl, file.path(final_dir, paste0("vtpbw_titration_cate_", site_name, ".csv")))
slope_tbl <- bind_rows(ct_disc$slope, ct_pfvc$slope)
write_csv(slope_tbl, file.path(final_dir, paste0("vtpbw_titration_slope_", site_name, ".csv")))

# --- figure: one panel per modifier, ATE reference line ---------------------------------------
ate_pp <- 100 * ate$estimate   # official TMLE ATE as the figure reference (not the winsorized mean)
mk_panel <- function(ct, xlab) {
  ggplot(ct$fig, aes(x, cate)) +
    geom_hline(yintercept = ate_pp, linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, fill = "#0072B2") +
    geom_line(colour = "#0072B2", linewidth = 1) +
    labs(x = xlab, y = sprintf("CATE: 28-d mortality RD of VT/PBW -%.2f mL/kg (pp)", DELTA), title = ct$fig$modifier[1]) +
    theme_minimal(base_size = 11)
}
fig <- mk_panel(ct_disc, "PBW/PFVC discordance - higher = PBW oversizes (PFVC says smaller)") +
  mk_panel(ct_pfvc, "PFVC (L) - lower = smaller predicted lung") +
  patchwork::plot_annotation(
    title = paste0("CATE of an ADDITIVE VT/PBW titration (-", DELTA, " mL/kg), by PFVC-derived modifiers - ", site_name,
                   if (is_synthetic) " (SYNTHETIC)" else ""),
    subtitle = "Dashed = ATE; negative = benefit. Down-slope (left) = a fixed bedside VT/PBW cut helps most where PBW over-doses the lung (normalizer-dependent).")
ggsave(file.path(final_dir, paste0("vtpbw_titration_cate_", site_name, ".pdf")), fig, width = 12, height = 5)

cat("\n=== 11_vtpbw continuous CATE of an ADDITIVE VT/PBW titration, by PFVC-derived modifiers ===\n")
cat("--- effect-modification trends (read the SHAPE, not the level: confounding-by-severity inflates the ATE) ---\n")
print(as.data.frame(slope_tbl), row.names = FALSE)
cat("    (discordance grad_p90_minus_p10 MORE NEGATIVE => the misdosed gain more from a fixed bedside VT/PBW cut;\n")
cat("     per_log slope CI = analytic regression; gradient CI = patient bootstrap -- both conditional on the DR pseudo-outcomes)\n")
print(as.data.frame(cate_tbl), row.names = FALSE)
message("Wrote vtpbw_titration_{cate,slope}_", site_name, ".csv + cate .pdf to ", final_dir)
