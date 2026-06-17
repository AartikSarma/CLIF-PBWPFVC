# =============================================================================
# Script 10: Longitudinal target trial emulation of a strain-limiting strategy
#            (clone-censor-weight)  ***DRAFT SCAFFOLD -- harden before reporting***
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# Emulates the trial the cross-sectional analysis could not (positivity probe:
# decision-level overlap AUC ~0.66-0.70 vs cross-sectional c=0.996). Two
# sustained ceiling strategies on the size-relative dose:
#   * LOW ceiling  (strain-limiting): keep VT/PFVC <= C_LOW (11% predicted FVC)
#   * HIGH ceiling (permissive ~ usual care): keep VT/PFVC <= C_HIGH
# Per-protocol effect on 60-day mortality, by clone-censor-weight: clone each
# patient into both arms at t0 (index IMV), censor a clone at first deviation
# (exceeds its ceiling after a grace period), and inverse-probability-of-
# censoring weight to correct the informative (adherence) censoring. PFVC is the
# IMPLEMENTATION lever; the estimand is the strain-limiting STRATEGY effect.
#
# WHY THIS DESIGN: identified WITHOUT demographic positivity (the wall that
# killed the cross-sectional VT/PFVC contrast) and WITHOUT the height exclusion
# restriction (the IV's vulnerability) -- it triangulates with both under
# different assumptions. Bridges to the separate Bayesian reanalysis of
# completed-trial IPD (e.g. ARMA), reanalyzed by VT/PFVC.
#
# *** SCAFFOLD STATUS -- decisions to LOCK / pieces to HARDEN before any claim ***
#   [T1 DONE] MORTALITY follows to 60 d -- extubation is NOT a censoring event,
#        because post-extubation deaths are observed via death_dttm; only deviation
#        + admin censor. (Censoring at extubation was informative -> downward bias.)
#        LIBERATION is reported separately as a competing-risk secondary
#        (Aalen-Johansen, death competing, IPC-weighted).
#   [T2 DONE] Time-varying IPC-weighted pooled-logistic MSM (discrete-day hazard,
#        marginal cumulative incidence) + a proper cluster bootstrap (resamples
#        patients with replacement). NOTE: the bootstrap holds the IPCW model
#        fixed -> a modest under-estimate of uncertainty; a fully proper bootstrap
#        would refit the weight model per replicate.
#   [T6 DONE] Subgroup effect modification over age / within-sex height / sex /
#        race (point estimates; subgroup CIs are a fast follow).
#   [T7 DONE] Per-arm positivity/weight diagnostics (frac deviated, weight p99/max,
#        ESS fraction) written out; inspect before trusting an arm's estimate.
#   [T3 DONE] Weight-cap (3/5/10/Inf) + ceiling-grace grid sensitivities written
#        (tte_ccw_sens_weightcap_*, tte_ccw_sens_ceiling_grace_*).
#   [T4 DONE] Deviation rule reported both ways: "simple" (any post-grace
#        exceedance) and "corrected" (exceedance not brought back under by next
#        assessment) -- tte_ccw_sens_rule_*.
#   [T5 CLOSED by decision] Lactate deliberately excluded: noisy biomarker
#        (timing/clearance/indication confounding) that would cost more sample
#        than it adds. Confounder set = resp settings + S/F + MAP + vasopressor.
#   [T5b DONE] pH sensitivity: arterial pH (venous imputed as venous + 0.05) added
#        as a LAGGED time-varying confounder on the pH-covered subset. Respiratory
#        acidosis (permissive hypercapnia) is THE feedback that drives deviation
#        from a low-VT arm, so it is the most decision-relevant gas. Reported as a
#        SENSITIVITY, not a core-panel member, because ABG/VBG sampling is
#        indication-driven (missing-not-at-random). tte_ccw_sens_ph_*: full-cohort
#        vs pH-subset(no pH) vs pH-subset(+ lagged pH).
# =============================================================================

# Pin BLAS to one thread per process so the PSOCK bootstrap workers don't each
# run all-core multithreaded BLAS (N workers x cores threads -> oversubscription).
# Must precede any BLAS use; inherited by the workers spawned later. (cf. script 06.)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")

library(tidyverse)
library(arrow)
library(here)
library(splines)
library(survival)
library(parallel)

source("utils/config.R")
site_name  <- config$site_name
output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
RNGkind("L'Ecuyer-CMRG"); set.seed(20260617)

# --- design knobs (provisional; see [T3]) ------------------------------------
C_LOW        <- 11      # strain-limiting ceiling, VT/PFVC %
C_HIGH       <- 16      # permissive ceiling (~ usual care)
GRACE        <- 2L      # days allowed above ceiling before deviation
DAYW_CAP     <- 5       # day-weight truncation (IPCW clamped to [1/cap, cap]); see weight-cap sensitivity
HORIZON      <- 60L     # days, primary outcome
MAX_VENT_DAY <- 27L     # ventilation/adherence window
WT_TRUNC     <- c(0.01, 0.99)
is_synthetic <- identical(site_name, "synthetic_clif")
# cluster bootstrap reps; override on real data via PBWPFVC_NBOOT (each rep refits
# the MSM on a resampled long panel, so this is the dominant cost at real scale)
N_BOOT       <- if (is_synthetic) 100L else
  suppressWarnings(as.integer(Sys.getenv("PBWPFVC_NBOOT", "500")))
if (is.na(N_BOOT)) N_BOOT <- 500L
# bootstrap worker count; override with PBWPFVC_CORES (default leaves 1 core free)
N_CORES      <- suppressWarnings(as.integer(Sys.getenv("PBWPFVC_CORES", unset = NA)))
if (is.na(N_CORES)) N_CORES <- max(1L, detectCores() - 1L)
message("Site: ", site_name, " | ceilings ", C_LOW, "/", C_HIGH, " | grace ", GRACE, "d")

# =============================================================================
# 10a. Baseline (PFVC, demographics, t0, death day) -- mortality handling
#      identical to scripts 06-08 incl. synthetic-only workaround.
# =============================================================================
cs <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"))
rtrunc_lnorm <- function(n_needed, meanlog, sdlog, lo, hi) {
  acc <- numeric(0)
  while (length(acc) < n_needed) {
    cand <- rlnorm(max(n_needed * 2L, 1000L), meanlog, sdlog)
    cand <- cand[cand > lo & cand <= hi]; acc <- c(acc, cand) }
  acc[seq_len(n_needed)]
}
if (is_synthetic) {
  message("*** SYNTHETIC SITE: simulated survival (plumbing only; see script 06). ***")
  set.seed(20260615); n <- nrow(cs); died60 <- rbinom(n, 1L, 0.35)
  tte <- rep(NA_real_, n); tte[died60 == 1L] <- rtrunc_lnorm(sum(died60), log(9), 0.95, 0.04, HORIZON)
  cs <- cs %>% mutate(death_day = if_else(died60 == 1L, floor(tte), NA_real_))
} else {
  cs <- cs %>% mutate(idx = as.numeric(difftime(death_dttm, recorded_dttm, units = "days")),
                      death_day = if_else(!is.na(idx) & idx >= 0 & idx <= HORIZON, floor(idx), NA_real_))
}
age_breaks <- quantile(cs$age_at_admission, c(1/3, 2/3), na.rm = TRUE)
base <- cs %>%
  filter(!is.na(pfvc), pfvc > 0, !is.na(age_at_admission), !is.na(sex_category),
         !is.na(race_category), !is.na(sofa_total), !is.na(height_cm)) %>%
  group_by(sex_category) %>% mutate(height_z = as.numeric(scale(height_cm))) %>% ungroup() %>%
  transmute(hospitalization_id, t0 = recorded_dttm, pfvc, death_day,
            age10 = age_at_admission / 10, sex_category, race_category, sofa_total,
            age_grp = cut(age_at_admission, c(-Inf, age_breaks, Inf),
                          labels = c("Young", "Middle", "Old")),
            height_grp = cut(height_z, c(-Inf, quantile(height_z, c(1/3, 2/3), na.rm = TRUE), Inf),
                             labels = c("Short", "Middle", "Tall")))

# =============================================================================
# 10b. Daily exposure + time-varying confounder panel (reuses the probe build)
# =============================================================================
wf <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet")) %>%
  select(hospitalization_id, recorded_dttm, tidal_volume_set, fio2_set, peep_set, resp_rate_set) %>%
  filter(!is.na(tidal_volume_set), tidal_volume_set > 0) %>%
  inner_join(base %>% select(hospitalization_id, t0, pfvc), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  mutate(vtpfvc = tidal_volume_set / pfvc * 0.1)
daily <- wf %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(vtpfvc = median(vtpfvc, na.rm = TRUE), fio2 = median(fio2_set, na.rm = TRUE),
            peep = median(peep_set, na.rm = TRUE), rr = median(resp_rate_set, na.rm = TRUE),
            .groups = "drop")
vit <- read_parquet(file.path(output_dir, "cohort_vitals_clean.parquet")) %>%
  filter(vital_category %in% c("spo2", "map")) %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(recorded_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  group_by(hospitalization_id, vent_day, vital_category) %>%
  summarise(v = median(vital_value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = vital_category, values_from = v)
med <- read_parquet(file.path(output_dir, "cohort_meds.parquet")) %>%
  filter(med_group == "vasoactives") %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(admin_dttm, t0, units = "days")))) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  distinct(hospitalization_id, vent_day) %>% mutate(on_pressor = 1L)

panel <- daily %>%
  left_join(vit, by = c("hospitalization_id", "vent_day")) %>%
  left_join(med, by = c("hospitalization_id", "vent_day")) %>%
  left_join(base, by = "hospitalization_id") %>%
  mutate(on_pressor = coalesce(on_pressor, 0L),
         sf = spo2 / if_else(fio2 > 1.5, fio2 / 100, fio2)) %>%
  filter(is.finite(vtpfvc), is.finite(fio2), is.finite(peep), is.finite(rr),
         is.finite(sf), is.finite(map))
# extubation proxy = last observed vent-day + 1 ([T1])
extub <- panel %>% group_by(hospitalization_id) %>%
  summarise(extub_day = max(vent_day) + 1L, .groups = "drop")
panel <- panel %>% left_join(extub, by = "hospitalization_id")
message("Panel: ", nrow(panel), " patient-days, ", n_distinct(panel$hospitalization_id), " patients")

# Daily arterial-equivalent pH for the [T5b] sensitivity: pool arterial + venous
# gases, imputing arterial from venous as venous + 0.05 (venous pH runs ~0.03-0.05
# below arterial); take the daily median. Kept SEPARATE from the core panel because
# gas sampling is indication-driven (missing-not-at-random) -> sensitivity only.
ph_daily <- read_parquet(file.path(output_dir, "cohort_labs_clean.parquet")) %>%
  filter(lab_category %in% c("ph_arterial", "ph_venous"), !is.na(lab_value_numeric)) %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(vent_day = floor(as.numeric(difftime(lab_result_dttm, t0, units = "days"))),
         ph_art_eq = if_else(lab_category == "ph_venous",
                             lab_value_numeric + 0.05, lab_value_numeric)) %>%
  filter(vent_day >= 0, vent_day <= MAX_VENT_DAY) %>%
  group_by(hospitalization_id, vent_day) %>%
  summarise(ph = median(ph_art_eq, na.rm = TRUE), .groups = "drop")
message("pH panel: ", nrow(ph_daily), " patient-days with a gas, ",
        n_distinct(ph_daily$hospitalization_id), " patients")

# =============================================================================
# 10c. Parameterized arm builder: IPCW day-weights (lagged-confounder model)
# =============================================================================
# Deviation hazard from LAGGED confounders only (never the concurrent vtpfvc that
# DEFINES deviation). Parameterized by ceiling, grace, day-weight cap, and rule so
# the sensitivities below are one-liners.
#   rule = "simple":    any exceedance after grace is a deviation.
#   rule = "corrected": an exceedance is a deviation ONLY if NOT brought back under
#                       the ceiling by the next assessment ([T4] -- transient
#                       excursions that the clinician corrects are not deviations).
trunc_w <- function(w) { q <- quantile(w, WT_TRUNC, na.rm = TRUE); pmin(pmax(w, q[1]), q[2]) }
ess_frac <- function(w) { w <- trunc_w(w); (sum(w)^2 / sum(w^2)) / length(w) }

# pnl: the panel to build on (default = global `panel`; the pH sensitivity passes a
#      pH-augmented, pH-covered-patient subset). use_ph: add lagged pH (l_ph) to the
#      IPCW denominator and require it non-missing for eligibility ([T5b]).
arm_build <- function(ceiling, grace = GRACE, cap = DAYW_CAP, rule = "simple",
                      pnl = panel, use_ph = FALSE) {
  p <- pnl %>% group_by(hospitalization_id) %>% arrange(vent_day) %>%
    mutate(above = vent_day > grace & vtpfvc > ceiling, lead_vt = lead(vtpfvc),
           viol = if (rule == "corrected") above & (is.na(lead_vt) | lead_vt > ceiling) else above,
           prior_dev = lag(cumsum(viol), default = 0) > 0,
           l_vtpfvc = lag(vtpfvc), l_fio2 = lag(fio2), l_peep = lag(peep),
           l_rr = lag(rr), l_sf = lag(sf), l_map = lag(map), l_pressor = lag(on_pressor)) %>%
    ungroup()
  if (use_ph) p <- p %>% group_by(hospitalization_id) %>% arrange(vent_day) %>%
    mutate(l_ph = lag(ph)) %>% ungroup()
  fr <- p %>% filter(!prior_dev, vent_day > grace, !is.na(l_vtpfvc), !is.na(l_sf), !is.na(l_map))
  if (use_ph) fr <- fr %>% filter(!is.na(l_ph))
  num <- glm(viol ~ ns(vent_day, 3) + age10 + sex_category + race_category + sofa_total,
             data = fr, family = binomial)
  den_rhs <- paste("ns(vent_day, 3) + l_vtpfvc + l_fio2 + l_peep + l_rr + l_sf + l_map +",
                   "l_pressor + age10 + sex_category + race_category + sofa_total",
                   if (use_ph) "+ l_ph" else "")
  den <- glm(as.formula(paste("viol ~", den_rhs)), data = fr, family = binomial)
  p <- p %>% mutate(
    elig  = !prior_dev & vent_day > grace & !is.na(l_vtpfvc) & !is.na(l_sf) & !is.na(l_map) &
            (if (use_ph) !is.na(l_ph) else TRUE),
    p_num = predict(num, newdata = ., type = "response"),
    p_den = predict(den, newdata = ., type = "response"),
    day_w = if_else(elig, pmin(pmax((1 - p_num) / (1 - p_den), 1 / cap), cap), 1)) %>%
    group_by(hospitalization_id) %>% arrange(vent_day) %>%
    mutate(cumw = cumprod(day_w)) %>% ungroup()
  idsum <- p %>% group_by(hospitalization_id) %>% arrange(vent_day) %>%
    summarise(dev_day = { w <- which(viol); if (length(w)) min(vent_day[w]) else Inf },
              last_vent = max(vent_day), death_day = first(death_day), ipcw_term = last(cumw),
              age10 = first(age10), sex_category = first(sex_category),
              race_category = first(race_category), age_grp = first(age_grp),
              height_grp = first(height_grp), .groups = "drop")
  wday <- p %>% select(hospitalization_id, vent_day, cumw) %>%
    group_by(hospitalization_id) %>%
    complete(vent_day = full_seq(c(0, max(vent_day)), 1)) %>%
    fill(cumw, .direction = "down") %>% mutate(cumw = coalesce(cumw, 1)) %>% ungroup()
  list(idsum = idsum, wday = wday)
}

# =============================================================================
# 10d. Long panel + IPC-weighted MSM + competing-risk liberation, by design
# =============================================================================
make_long <- function(b, arm_lab) {
  b$idsum %>%
    mutate(mort_cens = pmin(dev_day, HORIZON),
           event_id  = is.finite(death_day) & death_day <= mort_cens,
           fu_end    = pmax(pmin(ifelse(event_id, death_day, mort_cens), HORIZON), 1)) %>%
    transmute(hospitalization_id, fu_end, event_id, last_vent, age_grp, height_grp,
              sex_category, race_category, arm = arm_lab) %>%
    cross_join(tibble(day = 1:HORIZON)) %>% filter(day <= fu_end) %>%
    mutate(died = as.integer(event_id & day == fu_end), vent_day = pmin(day, last_vent)) %>%
    left_join(b$wday, by = c("hospitalization_id", "vent_day")) %>%
    mutate(ipcw = coalesce(cumw, 1)) %>% select(-cumw)
}
build_design <- function(c_low, c_high, grace = GRACE, cap = DAYW_CAP, rule = "simple",
                         pnl = panel, use_ph = FALSE) {
  bl <- arm_build(c_low, grace, cap, rule, pnl, use_ph)
  bh <- arm_build(c_high, grace, cap, rule, pnl, use_ph)
  long <- bind_rows(make_long(bl, "strain_limiting"), make_long(bh, "permissive")) %>%
    mutate(arm = factor(arm, levels = c("permissive", "strain_limiting")))
  lib <- bind_rows(bl$idsum %>% mutate(arm = "strain_limiting"),
                   bh$idsum %>% mutate(arm = "permissive")) %>%
    mutate(dd = ifelse(is.finite(death_day), death_day, Inf), extub_day = last_vent + 1L,
           cr_end = pmin(dd, extub_day, dev_day, HORIZON),
           cr_status = case_when(dd == cr_end ~ 1L, extub_day == cr_end & extub_day < dd ~ 2L, TRUE ~ 0L),
           cr_time = pmax(cr_end, 1))
  list(long = long, lib = lib, bl = bl, bh = bh)
}
ci_curve <- function(dat) {
  fit <- suppressWarnings(glm(died ~ arm * ns(day, 4), data = dat, family = binomial, weights = ipcw))
  haz <- function(a) predict(fit, tibble(arm = factor(a, c("permissive", "strain_limiting")),
                                         day = 1:HORIZON), type = "response")
  tibble(day = 1:HORIZON, strain_limiting = 1 - cumprod(1 - haz("strain_limiting")),
         permissive = 1 - cumprod(1 - haz("permissive")))
}
rd_from <- function(dat) {
  cc <- ci_curve(dat); rs <- cc$strain_limiting[HORIZON]; rp <- cc$permissive[HORIZON]
  c(risk_sl = unname(rs), risk_pm = unname(rp), rd = unname(rs - rp))
}
cif_lib <- function(df, day = 28L) {
  if (sum(df$cr_status == 2L) < 5) return(NA_real_)
  df <- df %>% mutate(w = trunc_w(ipcw_term),
                      st = factor(cr_status, c(0L, 1L, 2L), c("censored", "death", "liberation")))
  fit <- survfit(Surv(cr_time, st) ~ 1, data = df, weights = w)
  if (!"liberation" %in% colnames(fit$pstate)) return(NA_real_)
  fit$pstate[max(which(fit$time <= day)), "liberation"]
}
lib_diff <- function(df) cif_lib(df %>% filter(arm == "strain_limiting")) -
                         cif_lib(df %>% filter(arm == "permissive"))

# per-(var, level) subgroup RD, used for point estimates AND inside the bootstrap
sub_vars <- c(age_grp = "Age tertile", height_grp = "Height tertile (w/in sex)",
              sex_category = "Sex", race_category = "Race")
sg_rd <- function(long) {
  imap_dfr(sub_vars, function(lab, v) {
    map_dfr(as.character(na.omit(unique(long[[v]]))), function(g) {
      d <- long %>% filter(as.character(.data[[v]]) == g)
      n <- n_distinct(d$hospitalization_id); if (n < 100) return(NULL)
      tibble(subgroup = lab, level = g, key = paste(lab, g, sep = "||"),
             rd = unname(rd_from(d)["rd"]), n = n)
    })
  })
}

# =============================================================================
# 10e. PRIMARY design + cluster bootstrap (overall RD + liberation + subgroup CIs)
# =============================================================================
des <- build_design(C_LOW, C_HIGH, GRACE, DAYW_CAP, "simple")
long_all <- des$long; lib_all <- des$lib
point   <- rd_from(long_all); lib_pt <- lib_diff(lib_all)
sg_point <- sg_rd(long_all)
ids <- unique(long_all$hospitalization_id)

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
# no incremental callback, so we dispatch the reps in chunks. PSOCK (separate
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
      bts_list[ch] <- parLapplyLB(cl, ch, function(bb) boot_one())
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
  lib_diff = lib_pt, lib_lo = ci("lib")[1], lib_hi = ci("lib")[2], n_patients = length(ids))
write_csv(overall, file.path(final_dir, paste0("tte_ccw_overall_", site_name, ".csv")))
message(sprintf("60-day MORTALITY RD: %.3f [%.3f, %.3f]  (%.3f vs %.3f)",
        overall$rd, overall$rd_lo, overall$rd_hi, overall$risk_strain_limiting, overall$risk_permissive))
message(sprintf("28-day LIBERATION CIF diff: %.3f [%.3f, %.3f]", overall$lib_diff, overall$lib_lo, overall$lib_hi))

# subgroup table WITH bootstrap CIs ([T6])
sub <- sg_point %>% rowwise() %>%
  mutate(rd_lo = ci(key)[1], rd_hi = ci(key)[2]) %>% ungroup() %>%
  select(subgroup, level, rd, rd_lo, rd_hi, n)
write_csv(sub, file.path(final_dir, paste0("tte_ccw_subgroup_", site_name, ".csv")))

# =============================================================================
# 10f. Sensitivities: weight cap + ceiling/grace grid ([T3]) + deviation rule ([T4])
# =============================================================================
cap_sens <- map_dfr(c(3, 5, 10, 1e6), function(cp) {
  d <- build_design(C_LOW, C_HIGH, GRACE, cp, "simple")
  tibble(weight_cap = if (cp >= 1e6) Inf else cp, rd = unname(rd_from(d$long)["rd"]),
         ess_strain_limiting = ess_frac(d$bl$idsum$ipcw_term),
         ess_permissive = ess_frac(d$bh$idsum$ipcw_term))
})
write_csv(cap_sens, file.path(final_dir, paste0("tte_ccw_sens_weightcap_", site_name, ".csv")))

t3_grid <- expand.grid(c_low = c(10, 11, 12), c_high = c(14, 16), grace = c(1L, 2L, 3L)) %>%
  filter(c_low < c_high)
t3_sens <- pmap_dfr(t3_grid, function(c_low, c_high, grace)
  tibble(c_low, c_high, grace,
         rd = unname(rd_from(build_design(c_low, c_high, grace, DAYW_CAP, "simple")$long)["rd"])))
write_csv(t3_sens, file.path(final_dir, paste0("tte_ccw_sens_ceiling_grace_", site_name, ".csv")))

t4_sens <- map_dfr(c("simple", "corrected"), function(rl) {
  d <- build_design(C_LOW, C_HIGH, GRACE, DAYW_CAP, rl)
  tibble(deviation_rule = rl, rd = unname(rd_from(d$long)["rd"]),
         frac_deviated_strain_limiting = mean(is.finite(d$bl$idsum$dev_day)))
})
write_csv(t4_sens, file.path(final_dir, paste0("tte_ccw_sens_rule_", site_name, ".csv")))

cat("\n=== weight-cap sensitivity ===\n");       print(as.data.frame(cap_sens %>% mutate(across(where(is.numeric), ~round(.,3)))), row.names = FALSE)
cat("\n=== deviation-rule sensitivity ([T4]) ===\n"); print(as.data.frame(t4_sens %>% mutate(across(where(is.numeric), ~round(.,3)))), row.names = FALSE)
cat("=== ceiling/grace grid ([T3]): RD range ", round(min(t3_sens$rd),3), " to ", round(max(t3_sens$rd),3), " ===\n")

# =============================================================================
# 10i. pH sensitivity ([T5b]): add lagged arterial pH (venous imputed +0.05) to the
#      IPCW denominator, on the pH-covered subset. Three rows isolate the question:
#        full_cohort_primary    -- the headline RD (no pH, all patients)
#        ph_subset_no_ph_adj    -- same RD re-estimated on pH-covered patients only
#                                  (shows whether that subset is itself selected)
#        ph_subset_with_ph_adj  -- pH-covered patients WITH lagged pH in the weight
#                                  model (the actual acidosis-adjusted estimate)
#      Stability across the last two = the strategy effect is not driven by
#      unmeasured respiratory acidosis (the permissive-hypercapnia feedback).
# =============================================================================
ph_ids   <- unique(ph_daily$hospitalization_id)
panel_ph <- panel %>% filter(hospitalization_id %in% ph_ids) %>%
  left_join(ph_daily, by = c("hospitalization_id", "vent_day")) %>%
  group_by(hospitalization_id) %>% arrange(vent_day) %>%
  fill(ph, .direction = "downup") %>% ungroup()   # LOCF + backfill within patient
n_ph <- n_distinct(panel_ph$hospitalization_id)
if (n_ph >= 100) {
  d_ph_base <- build_design(C_LOW, C_HIGH, GRACE, DAYW_CAP, "simple", panel_ph, use_ph = FALSE)
  d_ph_adj  <- build_design(C_LOW, C_HIGH, GRACE, DAYW_CAP, "simple", panel_ph, use_ph = TRUE)
  ph_sens <- tibble(
    spec = c("full_cohort_primary", "ph_subset_no_ph_adj", "ph_subset_with_ph_adj"),
    rd = c(unname(point["rd"]), unname(rd_from(d_ph_base$long)["rd"]),
           unname(rd_from(d_ph_adj$long)["rd"])),
    n_patients = c(length(ids), n_ph, n_ph),
    frac_of_cohort = round(c(1, n_ph / length(ids), n_ph / length(ids)), 3))
  write_csv(ph_sens, file.path(final_dir, paste0("tte_ccw_sens_ph_", site_name, ".csv")))
  cat("\n=== pH sensitivity ([T5b]: arterial; venous imputed +0.05) ===\n")
  print(as.data.frame(ph_sens %>% mutate(rd = round(rd, 3))), row.names = FALSE)
} else {
  message("pH sensitivity ([T5b]) skipped: ", n_ph,
          " patients with a gas (<100; e.g. synthetic CLIF has no pH labs).")
}

# =============================================================================
# 10g. Diagnostics: deviation, weights, per-arm positivity ([T7])
# =============================================================================
diag <- bind_rows(des$bl$idsum %>% mutate(arm = "strain_limiting"),
                  des$bh$idsum %>% mutate(arm = "permissive")) %>%
  group_by(arm) %>%
  summarise(n_patients = n(), frac_deviated = mean(is.finite(dev_day)),
            wt_p99 = quantile(trunc_w(ipcw_term), 0.99), wt_max = max(trunc_w(ipcw_term)),
            ess_frac = ess_frac(ipcw_term), .groups = "drop")
write_csv(diag, file.path(final_dir, paste0("tte_ccw_diagnostics_", site_name, ".csv")))
cat("\n=== CCW diagnostics (per-arm positivity / weights) ===\n")
print(as.data.frame(diag %>% mutate(across(where(is.numeric), ~ round(., 3)))), row.names = FALSE)

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
         " - 60-d mortality; subgroup CIs + weight-cap/ceiling-grace/rule sensitivities written")) +
  theme_minimal(base_size = 10) + theme(legend.position = "top")
ggsave(file.path(final_dir, paste0("tte_ccw_cuminc_", site_name, ".pdf")), p, width = 8, height = 5)
message("Wrote CCW tables + 1 figure to ", final_dir,
        "  [T1-T4,T6,T7 done; T5b pH sensitivity added; T5 lactate excluded by decision]")
