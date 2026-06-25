# =============================================================================
# Script 11.A: PRIMARY -- cluster bootstrap + overall RD + E-value + subgroup CIs
#              + cumulative-incidence figure (§10e bootstrap, §10h figure)
# =============================================================================
# Sources the shared common file (10_tte_common.R), which builds `des`, `long_all`,
# `lib_all`, `point`, `lib_pt`, `sg_point`, `ids`, the engine functions, and all
# knobs. This script runs ONLY the parts that follow the design build: the cluster
# bootstrap, the overall-RD/E-value/subgroup writes, and the cumulative-incidence PDF.
# =============================================================================
library(here)
source(here::here("code", "10_tte_engine.R"))

# =============================================================================
# 10e. PRIMARY design + cluster bootstrap (overall RD + liberation + subgroup CIs)
# =============================================================================
# one cluster-bootstrap replicate -> named vector: overall, lib, and each subgroup
sg_keys <- sg_point$key
boot_template <- setNames(rep(0, 2 + length(sg_keys)), c("overall", "lib", sg_keys))
boot_one <- function() {
  samp <- tibble(hospitalization_id = sample(ids, replace = TRUE))
  bl   <- long_all %>% inner_join(samp, by = "hospitalization_id", relationship = "many-to-many")
  blib <- lib_all  %>% inner_join(samp, by = "hospitalization_id", relationship = "many-to-many")
  out <- boot_template
  out["overall"] <- tryCatch(unname(rd_from(bl)["rd"]), error = function(e) NA_real_)
  out["lib"]     <- tryCatch(lib_diff(blib),            error = function(e) NA_real_)
  sg <- tryCatch(sg_rd(bl), error = function(e) NULL)
  if (!is.null(sg) && nrow(sg)) out[sg$key[sg$key %in% sg_keys]] <-
    sg$rd[sg$key %in% sg_keys]
  out
}
n_cores_used <- min(N_CORES, N_BOOT)
message("Cluster bootstrap (", N_BOOT, " reps across ", n_cores_used, " core(s); ",
        "overall + liberation + ", length(sg_keys), " subgroups) ...")
boot_t0 <- Sys.time()
report_boot <- function(done) {
  el <- as.numeric(difftime(Sys.time(), boot_t0, units = "secs"))
  message(sprintf("  bootstrap %d/%d (%2d%%) | elapsed %4.0fs | eta %4.0fs",
                  done, N_BOOT, as.integer(round(100 * done / N_BOOT)),
                  el, if (done < N_BOOT) el / done * (N_BOOT - done) else 0))
}
# Chunked so a progress line prints after each chunk (~20 ticks) -- parLapply has
# no incremental callback, so we dispatch the reps in chunks. We use parLapply, NOT
# parLapplyLB: non-load-balanced dispatch gives each worker a fixed contiguous block of
# reps, so with clusterSetRNGStream the bootstrap is reproducible at a given worker count
# (LB's dynamic, timing-dependent dispatch yields non-reproducible CIs). PSOCK (separate
# processes) avoids the fork + multithreaded-BLAS instability that can crash
# mclapply on macOS (cf. 06); finally{} tears the workers down on any exit.
chunks   <- split(seq_len(N_BOOT), cut(seq_len(N_BOOT), min(20L, N_BOOT), labels = FALSE))
bts_list <- vector("list", N_BOOT); done <- 0L
if (n_cores_used > 1) {
  cl <- makeCluster(n_cores_used, type = "PSOCK")
  clusterEvalQ(cl, { library(tidyverse); library(splines); library(survival) })
  # export only what boot_one() needs (NOT the big raw waterfall/vitals/meds)
  clusterExport(cl, envir = .GlobalEnv, varlist = c(
    "ids", "long_all", "lib_all", "boot_template", "sg_keys", "sub_vars",
    "HORIZON", "WT_TRUNC", "boot_one", "rd_from", "ci_curve", "cif_lib",
    "lib_diff", "sg_rd", "trunc_w"))
  clusterSetRNGStream(cl, 20260617)
  tryCatch(
    for (ch in chunks) {
      bts_list[ch] <- parLapply(cl, ch, function(bb) boot_one())
      done <- done + length(ch); report_boot(done)
    },
    finally = stopCluster(cl))
} else {
  for (ch in chunks) {
    for (b in ch) bts_list[[b]] <- boot_one()
    done <- done + length(ch); report_boot(done)
  }
}
bts <- do.call(rbind, bts_list)
ci  <- function(col) quantile(bts[, col], c(.025, .975), na.rm = TRUE)

overall <- tibble(
  risk_strain_limiting = unname(point["risk_sl"]), risk_permissive = unname(point["risk_pm"]),
  rd = unname(point["rd"]), rd_lo = ci("overall")[1], rd_hi = ci("overall")[2],
  lib_diff = lib_pt, lib_lo = ci("lib")[1], lib_hi = ci("lib")[2], n_patients = length(ids),
  # [T9] share of clones censored by the common-support trim at the PRIMARY TRIM_ALPHA.
  # Reported in the headline because the trim restricts the estimand to the overlap
  # region, so the trimmed fraction is part of describing WHICH population this RD is for
  # (and is the cross-site-comparability flag for pooling -- it differs MIMIC vs UCSF).
  frac_trimmed_strain     = mean(is.finite(des$bl$idsum$trim_day)),
  frac_trimmed_permissive = mean(is.finite(des$bh$idsum$trim_day)))
write_csv(overall, file.path(final_dir, paste0("tte_ccw_overall_", site_name, ".csv")))
message(sprintf("28-day MORTALITY RD: %.3f [%.3f, %.3f]  (%.3f vs %.3f)",
        overall$rd, overall$rd_lo, overall$rd_hi, overall$risk_strain_limiting, overall$risk_permissive))
message(sprintf("28-day LIBERATION CIF diff: %.3f [%.3f, %.3f]", overall$lib_diff, overall$lib_lo, overall$lib_hi))
message(sprintf("Common-support trim removed strain %.1f%% / permissive %.1f%% of clones (TRIM_ALPHA=%g)",
        100 * overall$frac_trimmed_strain, 100 * overall$frac_trimmed_permissive, TRIM_ALPHA))

# [T2b] E-value for the primary RD: the minimum association (risk-ratio scale) an
# unmeasured time-varying confounder of the ADHERENCE-censoring would need with BOTH
# deviation and mortality to explain the effect away. CAVEAT: the textbook E-value is
# for a point exposure; here it bounds residual confounding of the informative
# censoring, so read it as an approximate robustness index, not an exact bound. The CI
# bound nearest the null is mapped through the permissive risk; a null-crossing CI -> 1.
evalue <- function(rr) { rr <- if (rr >= 1) rr else 1 / rr; rr + sqrt(rr * (rr - 1)) }
rr_point   <- overall$risk_strain_limiting / overall$risk_permissive
rr_ci_lo   <- (overall$risk_permissive + overall$rd_lo) / overall$risk_permissive
rr_ci_hi   <- (overall$risk_permissive + overall$rd_hi) / overall$risk_permissive
rr_ci_near <- if (rr_point >= 1) min(rr_ci_lo, rr_ci_hi) else max(rr_ci_lo, rr_ci_hi)
rr_ci_near <- if ((rr_point >= 1) != (rr_ci_near >= 1)) 1 else rr_ci_near  # CI crosses null
eval_tbl <- tibble(rr_point = rr_point, evalue_point = evalue(rr_point),
                   rr_ci_bound = rr_ci_near, evalue_ci = evalue(rr_ci_near))
write_csv(eval_tbl, file.path(final_dir, paste0("tte_ccw_evalue_", site_name, ".csv")))
message(sprintf("E-value (approx; censoring-confounding): point %.2f, CI bound %.2f",
        eval_tbl$evalue_point, eval_tbl$evalue_ci))

# subgroup table WITH bootstrap CIs ([T6])
sub <- sg_point %>% rowwise() %>%
  mutate(rd_lo = ci(key)[1], rd_hi = ci(key)[2]) %>% ungroup() %>%
  select(subgroup, level, rd, rd_lo, rd_hi, n)
write_csv(sub, file.path(final_dir, paste0("tte_ccw_subgroup_", site_name, ".csv")))

# =============================================================================
# 10h. Figure: MSM per-protocol cumulative mortality by arm
# =============================================================================
cc <- ci_curve(long_all) %>%
  pivot_longer(c(strain_limiting, permissive), names_to = "arm", values_to = "cuminc")
p <- ggplot(cc, aes(day, 100 * cuminc, colour = arm)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c(strain_limiting = okabe[3], permissive = okabe[1]),
                      labels = c(strain_limiting = "strain-limiting (<=11%)",
                                 permissive = "permissive (<=16%)"), name = NULL) +
  labs(x = "Days from index ventilation", y = "Per-protocol cumulative mortality (%)",
       title = "Longitudinal TTE (CCW + IPC-weighted MSM): strain-limiting vs permissive",
       subtitle = paste0(site_name, if (is_synthetic) " (SYNTHETIC - plumbing only)" else "",
         " - 28-d mortality, common-support trimmed (P(adhere) >= ", TRIM_ALPHA,
         "); + de-escalation MTP / trim / weight-cap / ceiling-grace / rule sensitivities")) +
  theme_minimal(base_size = 10) + theme(legend.position = "top")
ggsave(file.path(final_dir, paste0("tte_ccw_cuminc_", site_name, ".pdf")), p, width = 8, height = 5)
message("Wrote CCW primary tables + 1 figure to ", final_dir,
        "  [T1-T10 done; PRIMARY = time-only-numerator weights (V-balancing, marginal ",
        "MSM valid) trimmed at TRIM_ALPHA=", TRIM_ALPHA, "]")
