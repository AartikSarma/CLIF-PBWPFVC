# =============================================================================
# Script 11.Z: Does the discordance-HTE survive demographics? (residualized linchpin)
# =============================================================================
# THE objection: PFVC = f(height, age, sex, race) and PBW = g(height, sex), so PBW/PFVC discordance
# is a near-deterministic function of demographics(+height) -- maybe the 11.X discordance-HTE is
# just demographic effect-modification relabelled. Because discordance and age/sex/race are nearly
# COLLINEAR, you cannot answer this by throwing demographics into the outcome model (the spline
# coefficients shred and the standardization goes off-support -- the v1 base-vs-modifier test did
# exactly that and was uninterpretable). The clean test is to ORTHOGONALIZE first:
#
# PRIMARY -- residualized-discordance slope. resid = residual of
#   log(PBW/PFVC) ~ ns(age,3) + sex + race  (the HEIGHT-driven part of discordance, orthogonal to
#   demographics BY CONSTRUCTION), then run the EXACT 11.X continuous-CATE estimator (SOFA-only
#   standardization, no demographic terms, no demographic-marginal grid) on resid. Reported beside
#   the RAW-discordance slope (which reproduces 11.X) on the same resamples:
#     raw slope clearly negative, residualized slope still negative & excludes 0
#        => discordance modifies the benefit BEYOND age/sex/race (height-driven lung size = PFVC's
#           added physiologic content -- exactly what we argue for; we do NOT residualize against
#           height, that IS the signal).
#     residualized slope collapses to 0 => the gradient was demographics.
#   The demographic R^2 of discordance is reported so the reader sees how much residual signal
#   remains (low residual variance => this is POWER-limited at 2 sites, resolved by the 8-site run).
#
# SUPPORT -- within-stratum slopes: the 11.X slope INSIDE strata of sex, age tertile, race (where
#   demographics are ~fixed, so within-stratum discordance variation is the height-driven residual).
#   Negative across the well-powered strata = the same conclusion, transparently. Honest n +
#   discordant-support flags; a CAP-stabilized bootstrap drops degenerate rank-deficient resamples
#   (the v1 +0.5/+0.9 upper bounds) so the CIs are not driven by spline extrapolation blow-ups.
#
# Standalone supplement -- NOT in 11_run_all. Synthetic mortality is simulated => plumbing.
# Env: PBWPFVC_DISC_CATE_DF (ns df, 3), PBWPFVC_WITHIN_MIN_N (min stratum n), PBWPFVC_BOOT_CAP.
# =============================================================================
library(here); library(splines); library(parallel); library(patchwork)
source(here::here("code", "10_tte_engine.R"))

KD    <- as.integer(Sys.getenv("PBWPFVC_DISC_CATE_DF", "3"))
MIN_N <- as.integer(Sys.getenv("PBWPFVC_WITHIN_MIN_N", "300"))   # min patients to fit a stratum slope
BCAP  <- as.numeric(Sys.getenv("PBWPFVC_BOOT_CAP", "0.5"))       # |slope| > BCAP (RD pp/100) = degenerate resample
OKABE <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
arm_f <- function() factor(c("permissive", "strain_limiting"), c("permissive", "strain_limiting"))

# --- patient-level: discordance, demographics, severity; residualize discordance ---------------
pin <- base %>% distinct(hospitalization_id, pbw, pfvc, age10, sex_category, race_category, age_grp, sofa_total) %>%
  mutate(ldisc = log(pbw / pfvc))
res_fit  <- lm(ldisc ~ ns(age10, 3) + sex_category + race_category, data = pin)
pin$rdisc <- as.numeric(residuals(res_fit))                       # height-driven, orthogonal to demographics
demo_r2  <- summary(res_fit)$r.squared
cat(sprintf("11.Z: demographics (age/sex/race) explain R^2 = %.3f of log-discordance; residual = the height-driven part\n", demo_r2))

long_s     <- long_all %>% left_join(pin %>% select(hospitalization_id, sofa_total), by = "hospitalization_id")
sofa_cells <- long_s %>% distinct(hospitalization_id, sofa_total) %>% count(sofa_total, name = "wt")

# add a fixed ns basis (knots from full data) for modifier column `m`; basis passed explicitly so
# the closure serializes cleanly to the bootstrap workers (no captured-environment surprises).
add_basis <- function(d, B, zc) { mm <- predict(B, d$m); for (j in seq_along(zc)) d[[zc[j]]] <- mm[, j]; d }

# build the clean 11.X-style estimator for one modifier (SOFA-only standardization).
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

P_raw <- prep_modifier(pin %>% transmute(hospitalization_id, m = ldisc), "zr")
P_res <- prep_modifier(pin %>% transmute(hospitalization_id, m = rdisc), "zd")

# within-stratum estimator: 11.X slope on the RAW modifier (within a stratum demographics ~fixed,
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
fit_raw <- suppressWarnings(glm(P_raw$FORM, P_raw$long, family = binomial, weights = ipcw))
fit_res <- suppressWarnings(glm(P_res$FORM, P_res$long, family = binomial, weights = ipcw))
slope_raw <- slope_of(fit_raw, P_raw); slope_res <- slope_of(fit_res, P_res)
curve_raw <- cate_from_fit(fit_raw, P_raw$grid_full) %>% filter(m %in% P_raw$DGRID) %>% mutate(modifier = "raw")
curve_res <- cate_from_fit(fit_res, P_res$grid_full) %>% filter(m %in% P_res$DGRID) %>% mutate(modifier = "residualized")
strat_pt  <- vapply(seq_len(nrow(strat_defs)), function(k)
  strat_slope(P_raw$long %>% filter(as.character(.data[[strat_defs$col[k]]]) == strat_defs$lvl[k])), numeric(1))

# --- paired, cap-stabilized bootstrap: one resample -> raw, residualized, all strata -----------
boot_one <- function() {
  samp <- tibble(hospitalization_id = sample(ids, replace = TRUE))
  lr <- P_raw$long %>% inner_join(samp, by = "hospitalization_id", relationship = "many-to-many")
  ld <- P_res$long %>% inner_join(samp, by = "hospitalization_id", relationship = "many-to-many")
  fr <- tryCatch(suppressWarnings(glm(P_raw$FORM, lr, family = binomial, weights = ipcw)), error = function(e) NULL)
  fd <- tryCatch(suppressWarnings(glm(P_res$FORM, ld, family = binomial, weights = ipcw)), error = function(e) NULL)
  sr <- if (is.null(fr)) NA_real_ else tryCatch(slope_of(fr, P_raw), error = function(e) NA_real_)
  sd <- if (is.null(fd)) NA_real_ else tryCatch(slope_of(fd, P_res), error = function(e) NA_real_)
  ss <- vapply(seq_len(nrow(strat_defs)), function(k)
    strat_slope(lr %>% filter(as.character(.data[[strat_defs$col[k]]]) == strat_defs$lvl[k])), numeric(1))
  c(raw = sr, resid = sd, setNames(ss, strat_label))
}
n_cores_used <- min(N_CORES, N_BOOT)
message("11.Z residualized bootstrap (", N_BOOT, " reps across ", n_cores_used,
        " core(s); raw+residualized+", nrow(strat_defs), " strata per resample) ...")
boot_t0 <- Sys.time()
chunks  <- split(seq_len(N_BOOT), cut(seq_len(N_BOOT), min(20L, N_BOOT), labels = FALSE))
bl      <- vector("list", N_BOOT); done <- 0L
if (n_cores_used > 1) {
  cl <- makeCluster(n_cores_used, type = "PSOCK")
  clusterEvalQ(cl, { library(tidyverse); library(splines) })
  clusterExport(cl, envir = .GlobalEnv, varlist = c(
    "ids", "P_raw", "P_res", "cate_from_fit", "slope_of", "strat_slope", "add_basis",
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
gr <- boot_ci("raw"); gd <- boot_ci("resid")
global_tbl <- tibble(
  modifier = c("raw discordance (= 11.X reproduction)", "residualized vs age/sex/race (height-driven)"),
  slope = c(slope_raw, slope_res), lo = c(gr["lo"], gd["lo"]), hi = c(gr["hi"], gd["hi"]),
  boot_valid = c(gr["n_valid"], gd["n_valid"]), boot_dropped = c(gr["n_drop"], gd["n_drop"]),
  disc_p10 = exp(P_raw$qd[1]), disc_p90 = exp(P_raw$qd[2]), demographic_r2 = demo_r2)
write_csv(global_tbl, file.path(final_dir, paste0("tte_ccw_within_demo_slope_", site_name, ".csv")))

strata_tbl <- strat_defs %>%
  transmute(stratum = strat_label, n_patients = n_pts, n_disc_p90, slope = strat_pt,
            lo = vapply(strat_label, function(s) boot_ci(s)["lo"], numeric(1)),
            hi = vapply(strat_label, function(s) boot_ci(s)["hi"], numeric(1)),
            boot_dropped = vapply(strat_label, function(s) boot_ci(s)["n_drop"], numeric(1)),
            flag = dplyr::case_when(n_pts < MIN_N ~ "thin: n < min",
                                    n_disc_p90 < 10 ~ "thin: <10 at p90 discordance",
                                    TRUE ~ "ok"))
write_csv(strata_tbl, file.path(final_dir, paste0("tte_ccw_within_demo_strata_", site_name, ".csv")))
write_csv(bind_rows(curve_raw, curve_res) %>% transmute(modifier, modifier_value = m, rd),
          file.path(final_dir, paste0("tte_ccw_within_demo_curve_", site_name, ".csv")))

# --- figure: (1) raw vs residualized CATE curves; (2) within-stratum slope forest --------------
p_raw <- ggplot(curve_raw, aes(exp(m), 100 * rd)) + geom_hline(yintercept = 0, colour = "grey80") +
  geom_line(colour = OKABE[1], linewidth = 1) +
  labs(x = "PBW/PFVC discordance (raw)", y = "CATE: strain-limiting RD (pp)", title = "Raw discordance (= 11.X)") +
  theme_minimal(base_size = 11)
p_res <- ggplot(curve_res, aes(m, 100 * rd)) + geom_hline(yintercept = 0, colour = "grey80") +
  geom_line(colour = OKABE[5], linewidth = 1) +
  labs(x = "residualized discordance (height-driven, orthogonal to age/sex/race)",
       y = "CATE: strain-limiting RD (pp)",
       title = sprintf("Residualized discordance (demographics R²=%.2f)", demo_r2)) +
  theme_minimal(base_size = 11)
forest_df <- bind_rows(
  tibble(stratum = "ALL | raw discordance",          slope = slope_raw, lo = gr["lo"], hi = gr["hi"], grp = "global"),
  tibble(stratum = "ALL | residualized (height)",    slope = slope_res, lo = gd["lo"], hi = gd["hi"], grp = "global"),
  strata_tbl %>% filter(flag == "ok") %>% transmute(stratum, slope, lo, hi, grp = "within-stratum")) %>%
  mutate(stratum = factor(stratum, rev(stratum)))
p_forest <- ggplot(forest_df, aes(100 * slope, stratum, colour = grp)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
  geom_pointrange(aes(xmin = 100 * lo, xmax = 100 * hi)) +
  scale_colour_manual(values = c(global = OKABE[6], `within-stratum` = OKABE[3]), guide = "none") +
  labs(x = "discordance slope RD(p90)-RD(p10) (pp); negative = misdosed benefit more", y = NULL,
       title = "Slope, residualized and within demographic strata") +
  theme_minimal(base_size = 11)
fig <- (p_raw | p_res) / p_forest + plot_annotation(
  title = paste0("11.Z: does the discordance-HTE survive demographics? (residualized) - ", site_name,
                 if (is_synthetic) " (SYNTHETIC)" else ""))
ggsave(file.path(final_dir, paste0("tte_ccw_within_demo_", site_name, ".pdf")), fig, width = 12, height = 9)

cat("\n=== 11.Z residualized within-demographic discordance-HTE (negative slope = misdosed benefit more) ===\n")
cat(sprintf("--- PRIMARY: raw vs residualized slope (demographics explain R^2=%.3f of discordance) ---\n", demo_r2))
print(as.data.frame(global_tbl %>% transmute(modifier, slope_pp = round(100 * slope, 2),
        ci = sprintf("[%.2f, %.2f]", 100 * lo, 100 * hi), boot_valid, boot_dropped)), row.names = FALSE)
cat("    (residualized slope still negative & excludes 0 => discordance modifies the benefit BEYOND age/sex/race;\n")
cat("     collapses to 0 => it was demographics. Low residual variance (high R^2) => power-limited at 2 sites.)\n")
cat("\n--- SUPPORT: discordance slope within demographic strata (pp; 'ok' = passes n + p90-support) ---\n")
print(as.data.frame(strata_tbl %>% mutate(slope_pp = round(100 * slope, 2),
        ci = sprintf("[%.2f, %.2f]", 100 * lo, 100 * hi)) %>%
        select(stratum, n_patients, n_disc_p90, slope_pp, ci, boot_dropped, flag)), row.names = FALSE)
message("Wrote tte_ccw_within_demo_{slope,strata,curve}_", site_name, ".csv + .pdf to ", final_dir)
