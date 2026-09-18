# =============================================================================
# Supplement.Z: Does the discordance-HTE survive demographics? (residualized linchpin)
# =============================================================================
# THE objection: PFVC = f(height, age, sex, race) and PBW = g(height, sex), so PBW/PFVC discordance
# is a near-deterministic function of demographics(+height) -- maybe the 37_tte_discordance_benefit discordance-HTE is
# just demographic effect-modification relabelled. Because discordance and age/sex/race are nearly
# COLLINEAR, you cannot answer this by throwing demographics into the outcome model (the spline
# coefficients shred and the standardization goes off-support -- the v1 base-vs-modifier test did
# exactly that and was uninterpretable). The clean test is to ORTHOGONALIZE first:
#
# PRIMARY -- residualized-discordance slope. resid = residual of
#   log(PBW/PFVC) ~ ns(age,3) + sex + race  (the HEIGHT-driven part of discordance, orthogonal to
#   demographics BY CONSTRUCTION), then run the EXACT 37_tte_discordance_benefit continuous-CATE estimator (SOFA-only
#   standardization, no demographic terms, no demographic-marginal grid) on resid. Reported beside
#   the RAW-discordance slope (which reproduces 37_tte_discordance_benefit) on the same resamples:
#     raw slope clearly negative, residualized slope still negative & excludes 0
#        => discordance modifies the benefit BEYOND age/sex/race (height-driven lung size = PFVC's
#           added physiologic content -- exactly what we argue for; we do NOT residualize against
#           height, that IS the signal).
#     residualized slope collapses to 0 => the gradient was demographics.
#   The demographic R^2 of discordance is reported so the reader sees how much residual signal
#   remains (low residual variance => this is POWER-limited at 2 sites, resolved by the 8-site run).
#
# SUPPORT -- within-stratum slopes: the 37_tte_discordance_benefit slope INSIDE strata of sex, age tertile, race (where
#   demographics are ~fixed, so within-stratum discordance variation is the height-driven residual).
#   Negative across the well-powered strata = the same conclusion, transparently. Honest n +
#   discordant-support flags; a CAP-stabilized bootstrap drops degenerate rank-deficient resamples
#   (the v1 +0.5/+0.9 upper bounds) so the CIs are not driven by spline extrapolation blow-ups.
#
# Standalone supplement -- NOT in 32_tte_run_all. Synthetic mortality is simulated => plumbing.
# Env: PBWPFVC_DISC_CATE_DF (ns df, 3), PBWPFVC_WITHIN_MIN_N (min stratum n), PBWPFVC_BOOT_CAP.
# =============================================================================
library(here); library(splines); library(parallel); library(patchwork)
source(here::here("code", "31_tte_engine.R"))

KD    <- as.integer(Sys.getenv("PBWPFVC_DISC_CATE_DF", "3"))
MIN_N <- as.integer(Sys.getenv("PBWPFVC_WITHIN_MIN_N", "300"))   # min patients to fit a stratum slope
BCAP  <- as.numeric(Sys.getenv("PBWPFVC_BOOT_CAP", "0.5"))       # |slope| > BCAP (RD pp/100) = degenerate resample
OKABE <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
arm_f <- function() factor(c("permissive", "strain_limiting"), c("permissive", "strain_limiting"))

# --- patient-level: discordance, demographics, severity; residualize discordance ---------------
pin <- base %>% distinct(hospitalization_id, pbw, pfvc, age10, sex_category, race_category, age_grp, sofa_total) %>%
  mutate(ldisc = log(pbw / pfvc))
pin$lpfvc <- log(pin$pfvc)
res_fit  <- lm(ldisc ~ ns(age10, 3) + sex_category + race_category, data = pin)
pin$rdisc <- as.numeric(residuals(res_fit))                       # height-driven, orthogonal to demographics
demo_r2  <- summary(res_fit)$r.squared
# The SAME orthogonalization applied to log PFVC. This is the better-powered version of the test:
# demographics explain ~99% of log-discordance but far less of log PFVC, because PFVC depends
# strongly on height while the ratio barely does. So residualized PFVC keeps several times the
# variance residualized discordance does, and what survives is the height direction the
# negative-control cohorts showed to be confined to ventilated patients.
res_fit_p  <- lm(lpfvc ~ ns(age10, 3) + sex_category + race_category, data = pin)
pin$rpfvc  <- as.numeric(residuals(res_fit_p))
demo_r2_p  <- summary(res_fit_p)$r.squared
cat(sprintf("tte_within_demographic_hte: demographics (age/sex/race) explain R^2 = %.3f of log-discordance and %.3f of log-PFVC; residual = the height-driven part\n",
            demo_r2, demo_r2_p))
cat(sprintf("      residual SD as a share of raw: discordance %.3f, PFVC %.3f -- the ratio of these is how much better powered the PFVC test is\n",
            sqrt(1 - demo_r2), sqrt(1 - demo_r2_p)))

long_s     <- long_all %>% left_join(pin %>% select(hospitalization_id, sofa_total), by = "hospitalization_id")
sofa_cells <- long_s %>% distinct(hospitalization_id, sofa_total) %>% count(sofa_total, name = "wt")

# add a fixed ns basis (knots from full data) for modifier column `m`; basis passed explicitly so
# the closure serializes cleanly to the bootstrap workers (no captured-environment surprises).
add_basis <- function(d, B, zc) { mm <- predict(B, d$m); for (j in seq_along(zc)) d[[zc[j]]] <- mm[, j]; d }

# build the clean 37_tte_discordance_benefit-style estimator for one modifier (SOFA-only standardization).
prep_modifier <- function(mod_tbl, prefix) {
  mv <- mod_tbl$m
  B  <- ns(mv, df = KD); zc <- paste0(prefix, seq_len(KD)); zt <- paste(zc, collapse = "+")
  long_m  <- add_basis(long_s %>% left_join(mod_tbl, by = "hospitalization_id"), B, zc)
  FORM <- as.formula(paste0("died ~ arm*(", zt, ") + arm*ns(day,4) + (", zt, ")*ns(day,4) + sofa_total"))
  DGRID <- seq(quantile(mv, .025), quantile(mv, .975), length.out = 40)
  qd <- quantile(mv, c(.10, .90)); EVAL <- sort(unique(c(DGRID, qd)))
  grid_full <- add_basis(tidyr::crossing(m = EVAL, sofa_cells, day = 1:HORIZON, arm = arm_f()), B, zc)
  grid_anc  <- add_basis(tidyr::crossing(m = qd,   sofa_cells, day = 1:HORIZON, arm = arm_f()), B, zc)
  list(long = long_m, FORM = FORM, grid_full = grid_full, grid_anc = grid_anc,
       qd = qd, DGRID = DGRID, B = B, zc = zc)
}
# standardized CATE(m) = avg daily hazard over the SOFA cells per arm -> CIF -> RD, by modifier value
cate_from_fit <- function(fit, grid) {
  g <- grid; g$haz <- predict(fit, g, type = "response")
  g %>% group_by(m, arm, day) %>% summarise(h = weighted.mean(haz, wt), .groups = "drop") %>%
    group_by(m, arm) %>% arrange(day) %>% summarise(cif = 1 - prod(1 - h), .groups = "drop") %>%
    group_by(m) %>% summarise(rd = cif[arm == "strain_limiting"] - cif[arm == "permissive"], .groups = "drop") %>%
    arrange(m)
}
slope_of <- function(fit, P) {                     # RD(p90) - RD(p10), at the 2 anchors (fast)
  cc <- cate_from_fit(fit, P$grid_anc)
  cc$rd[which.min(abs(cc$m - P$qd[2]))] - cc$rd[which.min(abs(cc$m - P$qd[1]))]
}

P_raw  <- prep_modifier(pin %>% transmute(hospitalization_id, m = ldisc), "zr")
P_res  <- prep_modifier(pin %>% transmute(hospitalization_id, m = rdisc), "zd")
P_praw <- prep_modifier(pin %>% transmute(hospitalization_id, m = lpfvc), "zp")
P_pres <- prep_modifier(pin %>% transmute(hospitalization_id, m = rpfvc), "zq")
# The p90-minus-p10 gradient is RANGE-DEPENDENT, so it cannot be compared between a raw modifier
# and its residual: residualizing removes most of the variance, so the same per-unit effect yields
# a much smaller gradient. Every modifier therefore also reports its p10-p90 span and the implied
# PER-UNIT slope (gradient / span), which is the quantity that IS comparable across axes.
mod_span <- function(P) unname(diff(P$qd))

# within-stratum estimator: 37_tte_discordance_benefit slope on the RAW modifier (within a stratum demographics ~fixed,
# so raw within-stratum variation IS the height-driven residual), at the SHARED p10/p90 anchors.
FORM_strat <- P_raw$FORM
strat_slope <- function(dat) {
  if (n_distinct(dat$hospitalization_id) < MIN_N) return(NA_real_)
  sof <- dat %>% distinct(hospitalization_id, sofa_total) %>% count(sofa_total, name = "wt")
  fit <- tryCatch(suppressWarnings(glm(FORM_strat, dat, family = binomial, weights = ipcw)),
                  error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  grid <- add_basis(tidyr::crossing(m = P_raw$qd, sof, day = 1:HORIZON, arm = arm_f()), P_raw$B, P_raw$zc)
  cc   <- tryCatch(cate_from_fit(fit, grid), error = function(e) NULL)
  if (is.null(cc) || nrow(cc) < 2) return(NA_real_) else cc$rd[2] - cc$rd[1]   # qd sorted: p10,p90
}
strat_defs <- bind_rows(
  tibble(var = "sex",     col = "sex_category",  lvl = sort(unique(as.character(pin$sex_category)))),
  tibble(var = "age_grp", col = "age_grp",       lvl = sort(unique(as.character(pin$age_grp)))),
  tibble(var = "race",    col = "race_category", lvl = sort(unique(as.character(pin$race_category))))) %>%
  rowwise() %>%
  mutate(n_pts      = sum(as.character(pin[[col]]) == lvl),
         n_disc_p90 = sum(as.character(pin[[col]]) == lvl & pin$ldisc >= P_raw$qd[2])) %>% ungroup()
strat_label <- paste0(strat_defs$var, "=", strat_defs$lvl)

# --- point estimates ---------------------------------------------------------------------------
fit_of <- function(P) suppressWarnings(glm(P$FORM, P$long, family = binomial, weights = ipcw))
fit_raw <- fit_of(P_raw); fit_res <- fit_of(P_res); fit_praw <- fit_of(P_praw); fit_pres <- fit_of(P_pres)
slope_raw  <- slope_of(fit_raw,  P_raw);  slope_res  <- slope_of(fit_res,  P_res)
slope_praw <- slope_of(fit_praw, P_praw); slope_pres <- slope_of(fit_pres, P_pres)
curve_of <- function(fit, P, lab) cate_from_fit(fit, P$grid_full) %>% filter(m %in% P$DGRID) %>% mutate(modifier = lab)
curve_raw  <- curve_of(fit_raw,  P_raw,  "raw")
curve_res  <- curve_of(fit_res,  P_res,  "residualized")
curve_praw <- curve_of(fit_praw, P_praw, "pfvc_raw")
curve_pres <- curve_of(fit_pres, P_pres, "pfvc_residualized")
strat_pt  <- vapply(seq_len(nrow(strat_defs)), function(k)
  strat_slope(P_raw$long %>% filter(as.character(.data[[strat_defs$col[k]]]) == strat_defs$lvl[k])), numeric(1))

# --- paired, cap-stabilized bootstrap: one resample -> raw, residualized, all strata -----------
boot_one <- function() {
  samp <- tibble(hospitalization_id = sample(ids, replace = TRUE))
  one <- function(P) {
    l <- P$long %>% inner_join(samp, by = "hospitalization_id", relationship = "many-to-many")
    f <- tryCatch(suppressWarnings(glm(P$FORM, l, family = binomial, weights = ipcw)), error = function(e) NULL)
    if (is.null(f)) NA_real_ else tryCatch(slope_of(f, P), error = function(e) NA_real_)
  }
  lr <- P_raw$long %>% inner_join(samp, by = "hospitalization_id", relationship = "many-to-many")
  ss <- vapply(seq_len(nrow(strat_defs)), function(k)
    strat_slope(lr %>% filter(as.character(.data[[strat_defs$col[k]]]) == strat_defs$lvl[k])), numeric(1))
  c(raw = one(P_raw), resid = one(P_res), pfvc_raw = one(P_praw), pfvc_resid = one(P_pres),
    setNames(ss, strat_label))
}
n_cores_used <- min(N_CORES, N_BOOT)
message("tte_within_demographic_hte residualized bootstrap (", N_BOOT, " reps across ", n_cores_used,
        " core(s); raw+residualized+", nrow(strat_defs), " strata per resample) ...")
boot_t0 <- Sys.time()
chunks  <- split(seq_len(N_BOOT), cut(seq_len(N_BOOT), min(20L, N_BOOT), labels = FALSE))
bl      <- vector("list", N_BOOT); done <- 0L
if (n_cores_used > 1) {
  cl <- makeCluster(n_cores_used, type = "PSOCK")
  clusterEvalQ(cl, { library(tidyverse); library(splines) })
  clusterExport(cl, envir = .GlobalEnv, varlist = c(
    "ids", "P_raw", "P_res", "P_praw", "P_pres", "cate_from_fit", "slope_of", "strat_slope", "add_basis",
    "FORM_strat", "KD", "HORIZON", "MIN_N", "arm_f", "strat_defs", "strat_label", "boot_one"))
  clusterSetRNGStream(cl, 20260626)
  tryCatch(for (ch in chunks) {
    bl[ch] <- parLapply(cl, ch, function(bb) boot_one())
    done <- done + length(ch)
    message(sprintf("  boot %d/%d (%2d%%) | %4.0fs", done, N_BOOT,
                    as.integer(round(100 * done / N_BOOT)),
                    as.numeric(difftime(Sys.time(), boot_t0, units = "secs"))))
  }, finally = stopCluster(cl))
} else for (ch in chunks) {
  for (b in ch) bl[[b]] <- boot_one()
  done <- done + length(ch); message(sprintf("  boot %d/%d", done, N_BOOT))
}
bts <- do.call(rbind, bl)
# cap-stabilized percentile CI: drop non-finite + |slope|>BCAP (degenerate rank-deficient resamples)
boot_ci <- function(col) {
  v  <- bts[, col]; vv <- v[is.finite(v) & abs(v) <= BCAP]
  c(lo = unname(quantile(vv, .025)), hi = unname(quantile(vv, .975)),
    n_valid = length(vv), n_drop = sum(is.finite(v)) - length(vv))
}

# --- outputs -----------------------------------------------------------------------------------
gr <- boot_ci("raw"); gd <- boot_ci("resid"); gpr <- boot_ci("pfvc_raw"); gpd <- boot_ci("pfvc_resid")
global_tbl <- tibble(
  modifier = c("discordance, raw (= 11.X reproduction)", "discordance, residualized vs age/sex/race",
               "PFVC, raw", "PFVC, residualized vs age/sex/race"),
  base_variable = c("log PBW/PFVC", "log PBW/PFVC", "log PFVC", "log PFVC"),
  residualized = c(FALSE, TRUE, FALSE, TRUE),
  slope = c(slope_raw, slope_res, slope_praw, slope_pres),
  lo = c(gr["lo"], gd["lo"], gpr["lo"], gpd["lo"]),
  hi = c(gr["hi"], gd["hi"], gpr["hi"], gpd["hi"]),
  boot_valid = c(gr["n_valid"], gd["n_valid"], gpr["n_valid"], gpd["n_valid"]),
  boot_dropped = c(gr["n_drop"], gd["n_drop"], gpr["n_drop"], gpd["n_drop"]),
  # p10-p90 span of the modifier on its own (log or residual-log) scale, and the gradient divided
  # by it: the gradient is range-dependent and NOT comparable across these four rows, the per-unit
  # slope is. A residualized row with a small gradient but a large per-unit slope is UNDERPOWERED,
  # not null -- read the two together.
  span_p10_p90 = c(mod_span(P_raw), mod_span(P_res), mod_span(P_praw), mod_span(P_pres)),
  demographic_r2 = c(demo_r2, demo_r2, demo_r2_p, demo_r2_p),
  disc_p10 = exp(P_raw$qd[1]), disc_p90 = exp(P_raw$qd[2])) %>%
  mutate(slope_per_unit = slope / span_p10_p90, lo_per_unit = lo / span_p10_p90,
         hi_per_unit = hi / span_p10_p90, site = site_name)
write_csv(global_tbl, file.path(final_dir, paste0("tte_ccw_within_demo_slope_", site_name, ".csv")))

strata_tbl <- strat_defs %>%
  transmute(stratum = strat_label, n_patients = n_pts, n_disc_p90, slope = strat_pt,
            lo = vapply(strat_label, function(s) boot_ci(s)["lo"], numeric(1)),
            hi = vapply(strat_label, function(s) boot_ci(s)["hi"], numeric(1)),
            boot_dropped = vapply(strat_label, function(s) boot_ci(s)["n_drop"], numeric(1)),
            # a stratum whose point estimate falls OUTSIDE its own bootstrap interval, or whose
            # bootstrap discarded a large share of resamples, is a degenerate fit, not a finding;
            # flagged here and dropped from the figure below.
            flag = dplyr::case_when(n_pts < MIN_N ~ "thin: n < min",
                                    n_disc_p90 < 10 ~ "thin: <10 at p90 discordance",
                                    !is.na(slope) & !is.na(lo) & !is.na(hi) &
                                      (slope < lo | slope > hi) ~ "degenerate: point outside its own CI",
                                    boot_dropped > 0.2 * N_BOOT ~ "unstable: >20% of resamples dropped",
                                    TRUE ~ "ok"))
write_csv(strata_tbl, file.path(final_dir, paste0("tte_ccw_within_demo_strata_", site_name, ".csv")))
write_csv(bind_rows(curve_raw, curve_res, curve_praw, curve_pres) %>%
            transmute(modifier, modifier_value = m, rd, site = site_name),
          file.path(final_dir, paste0("tte_ccw_within_demo_curve_", site_name, ".csv")))

# --- figure: CATE curves for both modifiers, the comparable per-unit slopes, and the forest ----
curve_panel <- function(cv, xlab, ttl, colour, exp_x = FALSE) {
  ggplot(cv, aes(if (exp_x) exp(m) else m, 100 * rd)) + geom_hline(yintercept = 0, colour = "grey80") +
    geom_line(colour = colour, linewidth = 1) +
    geom_rug(data = cv, aes(x = if (exp_x) exp(m) else m), inherit.aes = FALSE, alpha = 0.15,
             length = unit(0.02, "npc")) +
    labs(x = xlab, y = "CATE: strain-limiting RD (pp)", title = ttl) + theme_minimal(base_size = 11)
}
p_raw  <- curve_panel(curve_raw,  "PBW/PFVC discordance (raw)", "Discordance, raw (= 11.X)", OKABE[1], TRUE)
p_res  <- curve_panel(curve_res,  "residualized discordance (height-driven)",
                      sprintf("Discordance, residualized (demographics R2=%.2f)", demo_r2), OKABE[5])
p_praw <- curve_panel(curve_praw, "PFVC (L, raw)", "PFVC, raw", OKABE[4], TRUE)
p_pres <- curve_panel(curve_pres, "residualized log PFVC (height-driven)",
                      sprintf("PFVC, residualized (demographics R2=%.2f)", demo_r2_p), OKABE[6])

# the four global modifiers on the ONLY comparable scale (per log unit of the modifier)
perunit_df <- global_tbl %>%
  transmute(modifier, slope = slope_per_unit, lo = lo_per_unit, hi = hi_per_unit) %>%
  mutate(modifier = factor(modifier, rev(modifier)))
p_perunit <- ggplot(perunit_df, aes(100 * slope, modifier)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
  geom_pointrange(aes(xmin = 100 * lo, xmax = 100 * hi), colour = OKABE[6]) +
  labs(x = "slope per log unit of the modifier (pp)", y = NULL,
       title = "Per-unit slope (comparable across axes)",
       subtitle = "The p90-p10 gradient is range-dependent; residualizing shrinks the range, not necessarily the effect") +
  theme_minimal(base_size = 11)

forest_df <- bind_rows(
  tibble(stratum = "ALL | discordance, raw",          slope = slope_raw, lo = gr["lo"],  hi = gr["hi"],  grp = "global"),
  tibble(stratum = "ALL | discordance, residualized", slope = slope_res, lo = gd["lo"],  hi = gd["hi"],  grp = "global"),
  strat_tbl_ok <- strata_tbl %>% filter(flag == "ok") %>% transmute(stratum, slope, lo, hi, grp = "within-stratum")) %>%
  mutate(stratum = factor(stratum, rev(stratum)))
n_flagged <- sum(strata_tbl$flag != "ok")
p_forest <- ggplot(forest_df, aes(100 * slope, stratum, colour = grp)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
  geom_pointrange(aes(xmin = 100 * lo, xmax = 100 * hi)) +
  scale_colour_manual(values = c(global = OKABE[6], `within-stratum` = OKABE[3]), guide = "none") +
  labs(x = "discordance gradient RD(p90)-RD(p10) (pp); negative = misdosed benefit more", y = NULL,
       title = "Gradient, within demographic strata",
       subtitle = if (n_flagged > 0) paste0(n_flagged, " stratum/strata suppressed (thin or degenerate; see the CSV)") else NULL) +
  theme_minimal(base_size = 11)
fig <- (p_raw | p_res) / (p_praw | p_pres) / (p_perunit | p_forest) + plot_annotation(
  title = paste0("tte_within_demographic_hte: does the HTE survive demographics? discordance and PFVC, residualized - ", site_name,
                 if (is_synthetic) " (SYNTHETIC)" else ""),
  subtitle = "Within a stratum demographics are ~fixed, so within-stratum variation is the height-driven residual.")
ggsave(file.path(final_dir, paste0("tte_ccw_within_demo_", site_name, ".pdf")), fig, width = 12, height = 14)

cat("\n=== tte_within_demographic_hte residualized within-demographic discordance-HTE (negative slope = misdosed benefit more) ===\n")
cat(sprintf("--- PRIMARY: raw vs residualized slope (demographics explain R^2=%.3f of discordance) ---\n", demo_r2))
print(as.data.frame(global_tbl %>% transmute(modifier, demo_r2 = round(demographic_r2, 3),
        span = round(span_p10_p90, 3), gradient_pp = round(100 * slope, 2),
        gradient_ci = sprintf("[%.2f, %.2f]", 100 * lo, 100 * hi),
        per_unit_pp = round(100 * slope_per_unit, 1),
        per_unit_ci = sprintf("[%.0f, %.0f]", 100 * lo_per_unit, 100 * hi_per_unit),
        boot_dropped)), row.names = FALSE)
cat("    READ THE PER-UNIT COLUMN when comparing a raw modifier with its residual: the gradient is\n")
cat("    RD(p90)-RD(p10) and residualizing shrinks the p10-p90 SPAN, so a small residualized gradient\n")
cat("    beside a large per-unit slope means UNDERPOWERED, not null. PFVC keeps far more variance than\n")
cat("    the discordance ratio after the same adjustment, so its residualized test is the better-powered one.\n")
cat("\n--- SUPPORT: discordance slope within demographic strata (pp; 'ok' = passes n + p90-support) ---\n")
print(as.data.frame(strata_tbl %>% mutate(slope_pp = round(100 * slope, 2),
        ci = sprintf("[%.2f, %.2f]", 100 * lo, 100 * hi)) %>%
        select(stratum, n_patients, n_disc_p90, slope_pp, ci, boot_dropped, flag)), row.names = FALSE)
message("Wrote tte_ccw_within_demo_{slope,strata,curve}_", site_name, ".csv + .pdf to ", final_dir)
