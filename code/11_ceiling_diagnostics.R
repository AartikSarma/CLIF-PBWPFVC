# =============================================================================
# 11_ceiling_diagnostics: the PRIMARY interval (weight-refit bootstrap) + the placebo contrast
# =============================================================================
# The design pieces (11_ceiling_headtohead.R / 11_ceiling_cap.R) use the fixed-weight cluster
# bootstrap every 11.* leaf uses: it resamples clones with the weights held fixed, so for two arms
# that share nearly all their clones the CI on their contrast conditions on one weight fit and is
# far too narrow. (The cap-on-top design made this visible: a Concordant RD of about -0.7 pp although
# the cap binds beyond the PBW ceiling on only ~1% of Concordant days -- the difference lived in the
# weights.) This piece therefore supplies, for EVERY design in PBWPFVC_CEIL_DIAG_DESIGNS (default
# "ceiling,cap"):
#
#   B. WEIGHT-REFIT bootstrap. Resample patients, REBUILD both arms (deviation models, weights,
#      trim) inside every replicate, refit the marginal and pooled SOFA-adjusted MSMs, and take
#      percentile CIs for the overall RD, the tertile RDs and the gradient. This is the honest CI
#      for a paired-arm contrast and is the PRIMARY interval for the manuscript; it is written as
#      tte_<fam>_<design>_overall_refit_* (one row per statistic) beside the design piece's
#      fixed-weight CI and the width ratio between them.
#   A. PLACEBO contrast. Concordant patients whose two clones are day-for-day identical in exposure
#      and deviation and differ ONLY in weights. Any nonzero RD there is pure weight artifact and
#      calibrates the floor the Mid and Discordant RDs must be read against. Reported with the
#      fixed-weight bootstrap CI (so the reader sees how narrow that CI is for a contrast that
#      should be exactly zero) and alongside the same contrast for Concordant patients WITH a
#      differing day. Most informative for the cap design (~95% of Concordant patients qualify).
#
# Refit reps: PBWPFVC_CEIL_REFIT_BOOT (default N_BOOT). Each replicate rebuilds both arms, so this
# is the expensive piece; at a large site lower the rep count via the env var rather than skip it.
# =============================================================================
library(here)
source(here::here("code", "11_ceiling_common.R"))

DIAG_DESIGNS <- strsplit(Sys.getenv("PBWPFVC_CEIL_DIAG_DESIGNS", "ceiling,cap"), ",")[[1]]
stopifnot(all(DIAG_DESIGNS %in% c("cap", "ceiling")))
N_REFIT <- suppressWarnings(as.integer(Sys.getenv("PBWPFVC_CEIL_REFIT_BOOT", unset = NA)))
if (is.na(N_REFIT)) N_REFIT <- N_BOOT
sofa    <- base %>% select(hospitalization_id, sofa_total)
ids_all <- unique(panel_x$hospitalization_id)

run_diagnostics <- function(DIAG) {
  tagd   <- function(x) paste0(PFX, "diag_", x, "_", DIAG, "_", site_name)
  prefix <- design_prefix(DIAG)
  des  <- build_ceiling_design(DIAG, keep_pday = FALSE)
  long <- des$long; ids <- unique(long$hospitalization_id)

  # ---------------------------------------------------------------------------
  # A. Placebo contrast: Concordant patients whose two clones differ only in the weights
  # ---------------------------------------------------------------------------
  # "differ" = the PFVC-informed arm's ceiling binds on a post-grace day where the PBW ceiling does
  # not. cap: cap_only; ceiling: any disagreement (either direction makes the clones differ).
  differ_ids <- if (DIAG == "cap") post %>% filter(cap_only) else post %>% filter(above_pfvc != above_pbw)
  differ_ids <- unique(differ_ids$hospitalization_id)
  conc_ids   <- base %>% filter(disc_grp == "Concordant") %>% pull(hospitalization_id)
  plac_ids   <- setdiff(intersect(conc_ids, ids), differ_ids)
  with_ids   <- intersect(intersect(conc_ids, ids), differ_ids)
  rd_sub <- function(id_set) {
    d <- long %>% filter(hospitalization_id %in% id_set)
    if (n_distinct(d$hospitalization_id) < 100) return(c(risk_sl = NA_real_, risk_pm = NA_real_, rd = NA_real_))
    rd_from(d)
  }
  w_sub <- function(b, id_set) mean(b$idsum$ipcw_term[b$idsum$hospitalization_id %in% id_set])
  pt_plac <- rd_sub(plac_ids); pt_with <- rd_sub(with_ids)
  boot_fixed <- if (length(plac_ids) >= 100) vapply(seq_len(N_BOOT), function(b) {
    samp <- tibble(hospitalization_id = sample(plac_ids, replace = TRUE))
    d <- long %>% inner_join(samp, by = "hospitalization_id", relationship = "many-to-many")
    tryCatch(unname(rd_from(d)["rd"]), error = function(e) NA_real_)
  }, numeric(1)) else NA_real_
  placebo <- tibble(
    family = FAM, design = DIAG,
    subset = c("concordant_no_differing_day (placebo)", "concordant_with_differing_day"),
    n_patients = c(length(plac_ids), length(with_ids)),
    risk_pfvc_arm = c(pt_plac["risk_sl"], pt_with["risk_sl"]), risk_pbw_arm = c(pt_plac["risk_pm"], pt_with["risk_pm"]),
    rd = c(pt_plac["rd"], pt_with["rd"]),
    rd_lo_fixed_weight_boot = c(quantile(boot_fixed, .025, na.rm = TRUE), NA_real_),
    rd_hi_fixed_weight_boot = c(quantile(boot_fixed, .975, na.rm = TRUE), NA_real_),
    mean_stab_w_pfvc_arm = c(w_sub(des$bl, plac_ids), w_sub(des$bl, with_ids)),
    mean_stab_w_pbw_arm  = c(w_sub(des$bh, plac_ids), w_sub(des$bh, with_ids)))
  write_csv(placebo, file.path(final_dir, paste0(tagd("placebo"), ".csv")))

  # ---------------------------------------------------------------------------
  # B. Weight-refit cluster bootstrap: rebuild both arms inside every replicate
  # ---------------------------------------------------------------------------
  refit_one <- function() {
    samp  <- tibble(hospitalization_id = sample(ids_all, replace = TRUE)) %>% mutate(k = row_number())
    pnl_b <- samp %>% inner_join(panel_x, by = "hospitalization_id", relationship = "many-to-many") %>%
      mutate(hospitalization_id = paste(hospitalization_id, k, sep = "_")) %>% select(-k)   # draws are distinct clusters
    out <- setNames(rep(NA_real_, 5), c("overall", DISC_LEVELS, "gradient"))
    des_b <- tryCatch(build_ceiling_design(DIAG, pnl = pnl_b), error = function(e) NULL)
    if (is.null(des_b)) return(out)
    lb <- des_b$long
    out["overall"] <- tryCatch(unname(rd_from(lb)["rd"]), error = function(e) NA_real_)
    ls <- lb %>% left_join(pnl_b %>% distinct(hospitalization_id, sofa_total), by = "hospitalization_id")
    pr <- ls %>% distinct(hospitalization_id, disc_grp, sofa_total)
    fb <- suppressWarnings(tryCatch(glm(FORM, data = ls, family = binomial, weights = ipcw), error = function(e) NULL))
    if (!is.null(fb)) {
      r <- vapply(DISC_LEVELS, function(g) unname(std_rd(fb, pr, g)["rd"]), numeric(1))
      out[DISC_LEVELS] <- r; out["gradient"] <- unname(r["Discordant"] - r["Concordant"])
    }
    out
  }
  n_cores_used <- min(N_CORES, N_REFIT)
  message("11_ceiling_diagnostics [", FAM, "]: weight-REFIT bootstrap for design '", DIAG, "' (", N_REFIT,
          " reps across ", n_cores_used, " core(s); each rep rebuilds both arms) ...")
  boot_t0 <- Sys.time()
  chunks  <- split(seq_len(N_REFIT), cut(seq_len(N_REFIT), min(20L, N_REFIT), labels = FALSE))
  rf_list <- vector("list", N_REFIT); done <- 0L
  if (n_cores_used > 1) {
    cl <- makeCluster(n_cores_used, type = "PSOCK")
    clusterEvalQ(cl, { library(tidyverse); library(splines); library(survival) })
    clusterExport(cl, envir = .GlobalEnv, varlist = c(
      "ids_all", "panel_x", "build_ceiling_design", "build_design", "arm_build", "make_long",
      "trunc_w", "rd_from", "ci_curve", "std_rd", "arm_f", "FORM", "DISC_LEVELS", "sf_for", "SHARED_LAGS",
      "tau_pbw", "tau_pfvc", "tau_cap", "HORIZON", "WT_TRUNC", "GRACE", "DAYW_CAP", "TRIM_ALPHA", "DEESC_FRAC"))
    clusterExport(cl, envir = environment(), varlist = c("DIAG", "refit_one"))
    clusterSetRNGStream(cl, if (DIAG == "cap") 20260906L else 20260907L)
    tryCatch(
      for (ch in chunks) {
        rf_list[ch] <- parLapply(cl, ch, function(bb) refit_one())
        done <- done + length(ch)
        message(sprintf("  refit bootstrap %d/%d (%2d%%) | elapsed %4.0fs", done, N_REFIT,
                        as.integer(round(100 * done / N_REFIT)),
                        as.numeric(difftime(Sys.time(), boot_t0, units = "secs"))))
      }, finally = stopCluster(cl))
  } else {
    for (ch in chunks) {
      for (b in ch) rf_list[[b]] <- refit_one()
      done <- done + length(ch); message(sprintf("  refit bootstrap %d/%d", done, N_REFIT))
    }
  }
  rf <- do.call(rbind, rf_list)
  n_fail <- sum(is.na(rf[, "overall"]))
  if (n_fail > 0) message(sprintf("  %d/%d refit replicates failed (design build or MSM) -> dropped", n_fail, N_REFIT))

  # point estimates from the ORIGINAL design (identical to the design piece's), CIs from the refit
  long_s <- long %>% left_join(sofa, by = "hospitalization_id")
  prof   <- long_s %>% distinct(hospitalization_id, disc_grp, sofa_total)
  fit0   <- suppressWarnings(glm(FORM, data = long_s, family = binomial, weights = ipcw))
  pt_t   <- vapply(DISC_LEVELS, function(g) unname(std_rd(fit0, prof, g)["rd"]), numeric(1))
  pt_all <- c(overall = unname(rd_from(long)["rd"]), pt_t, gradient = unname(pt_t["Discordant"] - pt_t["Concordant"]))
  # the design piece's fixed-weight CIs, for side-by-side (if that piece has been run)
  fixed <- local({
    f_ov <- file.path(final_dir, paste0(prefix, "overall_", site_name, ".csv"))
    f_ht <- file.path(final_dir, paste0(prefix, "disc_hte_", site_name, ".csv"))
    f_gr <- file.path(final_dir, paste0(prefix, "disc_gradient_", site_name, ".csv"))
    if (!all(file.exists(c(f_ov, f_ht, f_gr)))) return(tibble(statistic = character(), lo_fixed = numeric(), hi_fixed = numeric()))
    ov <- read_csv(f_ov, show_col_types = FALSE); ht <- read_csv(f_ht, show_col_types = FALSE); gr <- read_csv(f_gr, show_col_types = FALSE)
    bind_rows(tibble(statistic = "overall", lo_fixed = ov$rd_lo, hi_fixed = ov$rd_hi),
              ht %>% transmute(statistic = as.character(disc_grp), lo_fixed = rd_lo, hi_fixed = rd_hi),
              tibble(statistic = "gradient", lo_fixed = gr$lo, hi_fixed = gr$hi))
  })
  refit <- tibble(family = FAM, design = DIAG, statistic = names(pt_all), estimate = unname(pt_all),
                  lo_refit = apply(rf, 2, quantile, .025, na.rm = TRUE)[names(pt_all)],
                  hi_refit = apply(rf, 2, quantile, .975, na.rm = TRUE)[names(pt_all)],
                  sd_refit = apply(rf, 2, sd, na.rm = TRUE)[names(pt_all)],
                  n_reps_used = N_REFIT - n_fail, n_patients = length(ids),
                  ci_type = "weight_refit_bootstrap (PRIMARY)") %>%
    left_join(fixed, by = "statistic") %>%
    mutate(ci_width_ratio_refit_over_fixed = (hi_refit - lo_refit) / (hi_fixed - lo_fixed))
  # the primary-interval table the pooling script reads
  write_csv(refit, file.path(final_dir, paste0(prefix, "overall_refit_", site_name, ".csv")))

  cat(sprintf("\n=== 11_ceiling_diagnostics [%s] (design '%s') ===\n", FAM, DIAG))
  cat("--- A. Placebo: Concordant clones identical except for the weights (RD should be 0) ---\n")
  print(as.data.frame(placebo %>% select(-family, -design) %>% mutate(across(where(is.numeric), ~ round(., 4)))), row.names = FALSE)
  cat("--- B. PRIMARY weight-refit bootstrap vs the design piece's fixed-weight CI (pp) ---\n")
  print(as.data.frame(refit %>% transmute(statistic, estimate = round(100 * estimate, 2),
          refit_ci = sprintf("[%.2f, %.2f]", 100 * lo_refit, 100 * hi_refit),
          fixed_ci = ifelse(is.na(lo_fixed), NA, sprintf("[%.2f, %.2f]", 100 * lo_fixed, 100 * hi_fixed)),
          width_ratio = round(ci_width_ratio_refit_over_fixed, 2))), row.names = FALSE)
  cat("    (placebo RD != 0 with a fixed-weight CI excluding 0, but a refit CI covering 0 => weight-estimation noise the\n")
  cat("     design bootstrap cannot see; placebo RD != 0 with a refit CI ALSO excluding 0 => systematic weight bias)\n")
  message("Wrote ", prefix, "overall_refit_", site_name, ".csv (PRIMARY interval) and ", tagd("placebo"), ".csv to ", final_dir)
  invisible(list(placebo = placebo, refit = refit))
}

diag_results <- lapply(DIAG_DESIGNS, run_diagnostics)
names(diag_results) <- DIAG_DESIGNS
