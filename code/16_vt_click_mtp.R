# =============================================================================
# Script 16: one-click tidal-volume MTP at the index -- PFVC in one arm, PBW in the other
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Standalone; NOT in 00_run_pipeline.R. Reads the analytic cohort (analysis_cross_sectional,
# script 03: VT/PBW 6-8 mL/kg, SF < 315). Writes aggregates only.
#
# THE TARGET TRIAL. Time zero is the analytic index. At that moment the clinician has set a
# tidal volume. The trial randomizes what happens to that setting:
#   natural     the setting as observed
#   click_all   one click lower (CLICK_ML, 50 mL) for every patient in whom the lower value
#               stays inside the cohort's support (VT - 50 >= 6 mL/kg PBW); others unchanged
#   click_pfvc  one click lower only if VT/PFVC > TAU_PFVC (11% of predicted FVC, the ARMA
#               low-VT arm's p75): the PFVC-anchored arm
#   click_pbw   one click lower only if VT/PBW > TAU_PBW, with TAU_PBW set so the SAME share
#               of patients is shifted as under click_pfvc: the PBW-anchored arm, bite-matched
#   click_strain  a STRAIN-matched click for everyone: VT lower by CLICK_ML x PFVC / mean(PFVC)
#               mL, so every shifted patient loses the same VT/PFVC (CLICK_ML x 0.1 / mean PFVC
#               points) and the average cut in mL equals the fixed click's. The fallback's
#               companion (see CATE below).
# Every later setting follows the natural course given the (shifted) index setting, so the
# estimand is the effect of the initial setting together with its persistence. Outcome:
# all-cause death within 28 days of the index (out-of-hospital deaths included).
#
# WHY A POINT INTERVENTION, NOT A DAILY ONE. Set VT rarely changes from one day to the next, so
# a daily shift asks, given yesterday's VT, for a value that is almost never observed: the
# longitudinal density ratio collapses (11_vtpbw_titration.R, the sustained-dosing companion,
# fits its treatment density with a glm for exactly this reason, which stabilizes the variance by
# smoothing the positivity problem away). At the index the clinician's choice varies across
# patients with the same covariates, so a one-click shift stays in observed support.
#
# PRIMARY CONTRAST: click_pfvc - click_pbw. The two arms apply the identical shift to the same
# number of patients and differ only in WHO is shifted. Positivity is the same requirement in
# both arms: a click-lower setting must be observed among similar patients. No adherence process,
# no tail trimming. WHAT IT IDENTIFIES: whether a fixed VT cut helps more in the patients PFVC
# flags than in those PBW flags. That is the targeting question the paper asks, and it is
# identified. It is NOT lung mechanics: the selection difference between the arms is ~99%
# age, sex and race (11.Z), so a PFVC advantage says who to treat, not why. The mechanism under
# test is dose: a 50-mL click removes more strain from a smaller lung (1.7 points of VT/PFVC at
# 3 L, 1.0 at 5 L), and the PFVC arm selects smaller lungs, so its clicks deliver more strain
# reduction per patient. That is what PFVC anchoring is for.
#
# FALLBACK, PRE-SPECIFIED. The head-to-head is the headline only if, for BOTH targeted policies
# and in EVERY PFVC tertile, at most POS_MAX_SHARE_GT10 (5%) of density ratios exceed 10 and the
# effective sample size is at least POS_MIN_ESS_FRAC (50%) of the tertile, AND the two arms'
# selected sets overlap by at most OVERLAP_MAX_JACCARD (0.9; beyond that the arms are the same
# policy). Otherwise the headline is click_all versus natural with its CATE by PFVC. The script
# evaluates the rule and writes the decision (click_decision_{site}.csv); it is not a judgment call.
#
# CATE BY PFVC (always reported). DR-learner: per-patient pseudo-outcome difference between a
# click fit and the natural fit (estimate + influence function), smoothed on log PFVC and on log
# PBW/PFVC with a natural spline; p90 - p10 gradient with a patient bootstrap; and the relative
# (risk-ratio) curve, because low PFVC marks high baseline risk, so a sloped risk difference with a
# flat risk ratio is baseline risk, not effect modification.
# Read the two shifts together. Under click_all the dose itself varies with PFVC (a fixed 50 mL is
# a bigger strain cut in a smaller lung), so a CATE sloping toward small PFVC is EXPECTED with no
# effect modification at all; the risk-ratio check does not catch this. click_strain holds the
# strain cut constant. If the PFVC slope flattens under click_strain, the click_all slope was dose
# heterogeneity, which is itself evidence that strain is the dose; if it persists, it is effect
# modification. The strain-matched click lands between the 50-mL heaps, so check its positivity
# row. This block follows the titration script's CATE code in compact form (copied rather than moved to utils/, per the project's
# preference for inline analysis code).
#
# COVARIATES W: ns(age, 4), sex, race, ns(height, 4) x sex, SOFA, log SF ratio, ventilator mode.
# No BMI, per the project rule that BMI enters only pressure outcomes (it carries height, the
# exposure's identifying variation). --with_bmi adds it as a sensitivity. PBW, PFVC and the policy
# thresholds (functions of height, sex, age, race) ride along as columns because lmtp passes the
# shift function only model columns.
#
# COHORT. PRIMARY: index modes where the clinician sets VT (assist control-volume control,
# pressure-regulated volume control, SIMV). NOTE: script 03's VCV class holds only AC-VC and puts
# PRVC with pressure control; PRVC has a set VT, so it belongs here. That tidal_volume_set is the
# SET value in PRVC and SIMV rows is a CLIF convention this script cannot check; the mode table it
# writes shows what each site has. SENSITIVITIES: --modes all; --on_grid_only (VT on a 50-mL
# multiple, where a click lands on the next setting clinicians actually use).
#
# ESTIMATOR. lmtp_tmle with one time point (mtp = TRUE, binomial), cross-fit (5 folds). Outcome
# library {glm, gam, ranger, mean}; treatment-density library {glm, ranger, mean}: the ranger stays
# because set VT is heaped on 50-mL steps and a glm density would smooth the heaps away.
# Synthetic CLIF mortality is unreliable, so on the synthetic site survival is SIMULATED exactly as
# the TTE engine does (10_panel_common.R): plumbing only, never a result.
#
# Outputs (aggregates only; cells under MIN_CELL patients are suppressed):
#   final/click_policies_{site}.csv    per policy: share shifted, mean cut (mL, VT/PFVC points),
#                                      risk, RD vs natural; and the head-to-head contrast
#   final/click_selection_{site}.csv   share shifted by PFVC and discordance tertile, per policy
#   final/click_positivity_{site}.csv  density ratios per policy x PFVC tertile: share > 10, ESS
#   final/click_decision_{site}.csv    the fallback rule, evaluated
#   final/click_cate_{site}.csv        CATE curves (risk difference and risk ratio) by modifier
#   final/click_cate_slope_{site}.csv  p90 - p10 gradients
#   final/click_modes_{site}.csv       index ventilator modes in the analytic cohort
#   final/click_mtp_{site}.pdf
#   Non-default options add a suffix before the site name, so runs never overwrite each other.
# Usage: Rscript code/16_vt_click_mtp.R [--site_name NAME] [--output_root DIR]
#          [--modes vt_set|all] [--on_grid_only] [--with_bmi] [--tau_pfvc 11] [--click_ml 50]
#   Env: PBWPFVC_CLICK_FOLDS (5), PBWPFVC_CLICK_BOOT (500), PBWPFVC_CORES.
# =============================================================================
rm(list = ls())

# --- Command-line arguments (parsed the way code/00_run_pipeline.R parses its own) ------
parse_script_args <- function(args) {
  valued <- c("site_name", "output_root", "modes", "tau_pfvc", "click_ml")
  flags  <- c("on_grid_only", "with_bmi")
  usage <- paste("Usage: Rscript code/16_vt_click_mtp.R [--site_name NAME] [--output_root DIR]",
                 "[--modes vt_set|all] [--on_grid_only] [--with_bmi] [--tau_pfvc PCT] [--click_ml ML]")
  parsed <- list(); i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("--help", "-h")) { message(usage); quit(save = "no", status = 0) }
    if (a %in% paste0("--", flags)) { parsed[[sub("^--", "", a)]] <- TRUE; i <- i + 1L; next }
    if (grepl("^--[a-z_]+=", a)) {                       # --flag=value
      key <- sub("^--([a-z_]+)=.*$", "\\1", a); val <- sub("^--[a-z_]+=", "", a); i <- i + 1L
    } else if (grepl("^--[a-z_]+$", a)) {                # --flag value
      key <- sub("^--", "", a)
      if (!key %in% valued) stop("Unknown option --", key, "\n", usage)
      if (i == length(args) || grepl("^--", args[[i + 1L]]))
        stop("Missing value for --", key, "\n", usage)
      val <- args[[i + 1L]]; i <- i + 2L
    } else stop("Unrecognized argument: ", a, "\n", usage)
    if (!key %in% valued) stop("Unknown option --", key, "\n", usage)
    if (!nzchar(val)) stop("Empty value for --", key, "\n", usage)
    parsed[[key]] <- val
  }
  parsed
}
cli_args <- parse_script_args(commandArgs(trailingOnly = TRUE))

MODE_SET     <- if (is.null(cli_args$modes)) "vt_set" else cli_args$modes
if (!MODE_SET %in% c("vt_set", "all")) stop("--modes must be vt_set or all; got '", MODE_SET, "'")
ON_GRID_ONLY <- isTRUE(cli_args$on_grid_only)
USE_BMI      <- isTRUE(cli_args$with_bmi)
TAU_PFVC     <- if (is.null(cli_args$tau_pfvc)) 11 else as.numeric(cli_args$tau_pfvc)
CLICK_ML     <- if (is.null(cli_args$click_ml)) 50 else as.numeric(cli_args$click_ml)
if (!is.finite(TAU_PFVC) || TAU_PFVC <= 0) stop("--tau_pfvc must be a positive number (% of predicted FVC)")
if (!is.finite(CLICK_ML) || CLICK_ML <= 0) stop("--click_ml must be a positive number of mL")
run_suffix <- paste0(if (MODE_SET != "vt_set") "_allmodes" else "", if (ON_GRID_ONLY) "_ongrid" else "",
                     if (USE_BMI) "_bmi" else "", if (TAU_PFVC != 11) sprintf("_tau%g", TAU_PFVC) else "",
                     if (CLICK_ML != 50) sprintf("_click%g", CLICK_ML) else "")

# Run from the repository root whatever the caller's working directory. --output_root is
# resolved against the CALLER's directory first, before the setwd.
if (!is.null(cli_args$output_root))
  cli_args$output_root <- normalizePath(path.expand(cli_args$output_root), mustWork = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(script_file) == 1) setwd(normalizePath(file.path(dirname(script_file), "..")))
if (!is.null(cli_args$site_name)) Sys.setenv(PBWPFVC_SITE_NAME = cli_args$site_name)

suppressPackageStartupMessages({
  library(tidyverse); library(arrow); library(here); library(splines)
  library(lmtp); library(patchwork)
})
source("utils/config.R")
site_name    <- config$site_name
is_synthetic <- identical(site_name, "synthetic_clif")
output_root  <- if (is.null(cli_args$output_root)) here("output") else cli_args$output_root
output_dir   <- file.path(output_root, paste0(site_name, "_output"), "intermediate")
final_dir    <- file.path(output_root, paste0(site_name, "_output"), "final")
cohort_path  <- file.path(output_dir, "analysis_cross_sectional.parquet")
if (!file.exists(cohort_path))
  stop("Analytic cohort not found: ", cohort_path,
       "\nRun scripts 01-03 for site '", site_name, "' first, or check --site_name / --output_root.")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- function(stub, ext = "csv")
  file.path(final_dir, paste0(stub, run_suffix, "_", site_name, ".", ext))
message(sprintf("Site: %s | modes %s | on-grid only %s | BMI in W %s | tau_pfvc %g%% | click %g mL",
                site_name, MODE_SET, ON_GRID_ONLY, USE_BMI, TAU_PFVC, CLICK_ML))

HORIZON_DAYS <- 28
VTPBW_FLOOR  <- 6        # the analytic cohort's lower band edge: a click may not leave the support
MIN_CELL     <- 10
FOLDS        <- as.integer(Sys.getenv("PBWPFVC_CLICK_FOLDS", "5"))
N_BOOT       <- as.integer(Sys.getenv("PBWPFVC_CLICK_BOOT", "500"))
POS_MAX_SHARE_GT10  <- 0.05
POS_MIN_ESS_FRAC    <- 0.50
OVERLAP_MAX_JACCARD <- 0.90
VT_SET_MODES <- c("assist control-volume control", "pressure-regulated volume control", "simv")
NTHREAD <- suppressWarnings(as.integer(Sys.getenv("PBWPFVC_CORES", unset = NA)))
if (is.na(NTHREAD)) NTHREAD <- max(1L, parallel::detectCores())
SL.ranger.mc <- function(...) SuperLearner::SL.ranger(..., num.threads = NTHREAD)
LEARNERS_OUTCOME <- c("SL.glm", "SL.gam", "SL.ranger.mc", "SL.mean")
LEARNERS_TRT     <- c("SL.glm", "SL.ranger.mc", "SL.mean")
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
POLICY_LEVELS <- c("natural", "click_all", "click_pfvc", "click_pbw", "click_strain")
POLICY_COLOURS <- setNames(okabe[c(8, 3, 5, 6, 7)], POLICY_LEVELS)
RNGkind("L'Ecuyer-CMRG"); set.seed(20260918)

# =============================================================================
# 1. Cohort
# =============================================================================
cross_sectional <- read_parquet(cohort_path)
mode_table <- cross_sectional %>% mutate(mode = tolower(mode_category)) %>% count(mode, name = "n_patients") %>%
  mutate(share = n_patients / sum(n_patients), vt_set_mode = mode %in% VT_SET_MODES, site = site_name) %>%
  filter(n_patients >= MIN_CELL) %>% arrange(desc(n_patients))
write_csv(mode_table, out_file("click_modes"))

if (is_synthetic) {
  # 10_panel_common.R's synthetic survival, same seed and distribution: plumbing only
  message("*** SYNTHETIC SITE: simulated survival (plumbing only; synthetic CLIF mortality is unreliable). ***")
  rtrunc_lnorm <- function(n_needed, meanlog, sdlog, lo, hi) {
    accepted <- numeric(0)
    while (length(accepted) < n_needed) {
      candidates <- rlnorm(max(n_needed * 2L, 1000L), meanlog, sdlog)
      accepted <- c(accepted, candidates[candidates > lo & candidates <= hi])
    }
    accepted[seq_len(n_needed)]
  }
  set.seed(20260615); n_rows <- nrow(cross_sectional); died <- rbinom(n_rows, 1L, 0.35)
  days_to_death <- rep(NA_real_, n_rows)
  days_to_death[died == 1L] <- rtrunc_lnorm(sum(died), log(9), 0.95, 0.04, 60)
  cross_sectional$days_index_to_death <- days_to_death
  set.seed(20260918)
} else {
  cross_sectional <- cross_sectional %>%
    mutate(days_index_to_death = as.numeric(difftime(death_dttm, recorded_dttm, units = "days")))
}

cohort <- cross_sectional %>%
  mutate(mode = tolower(mode_category),
         died_28 = as.integer(!is.na(days_index_to_death) & days_index_to_death >= 0 &
                                days_index_to_death <= HORIZON_DAYS)) %>%
  filter(is.finite(tidal_volume_set), tidal_volume_set > 0, is.finite(pbw), pbw > 0,
         is.finite(pfvc), pfvc > 0, is.finite(age_at_admission), !is.na(sex_category),
         !is.na(race_category), is.finite(height_cm), is.finite(sofa_total), is.finite(sf_ratio),
         sf_ratio > 0, !is.na(mode))
n_before_filters <- nrow(cohort)
if (USE_BMI) cohort <- cohort %>% filter(is.finite(bmi))
if (MODE_SET == "vt_set") cohort <- cohort %>% filter(mode %in% VT_SET_MODES)
if (ON_GRID_ONLY) cohort <- cohort %>% filter(abs(tidal_volume_set - round(tidal_volume_set / 50) * 50) < 0.5)
message(sprintf("Cohort: %d of %d analytic patients with complete covariates (%s)", nrow(cohort), n_before_filters,
                paste(c(if (USE_BMI) "BMI required", if (MODE_SET == "vt_set") "VT-set modes",
                        if (ON_GRID_ONLY) "VT on 50-mL grid"), collapse = ", ")))
if (nrow(cohort) < 500) stop("Fewer than 500 patients after the cohort filters: not estimable. ",
                             "Check the mode table (", basename(out_file("click_modes")), ").")

# =============================================================================
# 2. Policies: thresholds, bite matching, and who each arm shifts
# =============================================================================
cohort <- cohort %>%
  mutate(vt = tidal_volume_set, vtpbw_index = vt / pbw, vtpfvc_index = vt / pfvc * 0.1,
         discordance = pbw / pfvc,
         vt_floor = VTPBW_FLOOR * pbw,             # lowest VT inside the cohort's support (mL)
         vt_cut_pfvc = 10 * TAU_PFVC * pfvc,       # VT above this <=> VT/PFVC > TAU_PFVC (%)
         feasible = vt - CLICK_ML >= vt_floor)
share_pfvc_arm <- mean(cohort$feasible & cohort$vt > cohort$vt_cut_pfvc)
if (share_pfvc_arm <= 0) stop("click_pfvc shifts nobody at tau_pfvc = ", TAU_PFVC, "%: lower --tau_pfvc.")
# bite matching: the VT/PBW threshold that shifts the same share of patients
TAU_PBW <- unname(quantile(cohort$vtpbw_index[cohort$feasible], 1 - share_pfvc_arm / mean(cohort$feasible), type = 1))
# the strain-matched click: the same VT/PFVC cut for every patient, equal on average in mL
STRAIN_CUT_POINTS <- CLICK_ML * 0.1 / mean(cohort$pfvc)
cohort <- cohort %>% mutate(vt_cut_pbw = TAU_PBW * pbw, strain_click_ml = CLICK_ML * pfvc / mean(pfvc))

# per-patient click size (mL) and who is selected, for each policy
click_size <- function(policy, d) if (policy == "click_strain") d$strain_click_ml else rep(CLICK_ML, nrow(d))
select_for <- function(policy, d) switch(policy,
  natural = rep(FALSE, nrow(d)), click_all = rep(TRUE, nrow(d)), click_strain = rep(TRUE, nrow(d)),
  click_pfvc = d$vt > d$vt_cut_pfvc, click_pbw = d$vt > d$vt_cut_pbw)
shifted_by <- function(policy, d) select_for(policy, d) & (d$vt - click_size(policy, d) >= d$vt_floor)
# the shift lmtp applies (it sees only model columns, which carry the thresholds and click sizes)
make_shift <- function(policy) function(data, trt) {
  d <- data.frame(vt = data[[trt]], vt_floor = data$vt_floor, vt_cut_pfvc = data$vt_cut_pfvc,
                  vt_cut_pbw = data$vt_cut_pbw, strain_click_ml = data$strain_click_ml)
  ifelse(shifted_by(policy, d), d$vt - click_size(policy, d), d$vt)
}

tertile <- function(x, labels) cut(x, quantile(x, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = labels)
cohort <- cohort %>%
  mutate(pfvc_tertile = tertile(pfvc, c("Small PFVC", "Mid PFVC", "Large PFVC")),
         discordance_tertile = tertile(discordance, c("Concordant", "Mid", "Discordant")))
for (policy in POLICY_LEVELS[-1]) cohort[[paste0("shifted_", policy)]] <- shifted_by(policy, cohort)
jaccard_arms <- with(cohort, sum(shifted_click_pfvc & shifted_click_pbw) / sum(shifted_click_pfvc | shifted_click_pbw))
message(sprintf("Thresholds: VT/PFVC > %g%% (PFVC arm) vs VT/PBW > %.2f mL/kg (PBW arm, bite-matched); each shifts %.1f%%; arms' overlap (Jaccard) %.2f",
                TAU_PFVC, TAU_PBW, 100 * share_pfvc_arm, jaccard_arms))

selection_table <- map_dfr(POLICY_LEVELS[-1], function(policy) {
  flag <- cohort[[paste0("shifted_", policy)]]
  bind_rows(
    cohort %>% mutate(shifted = flag) %>% group_by(modifier = "PFVC tertile", group = as.character(pfvc_tertile)) %>%
      summarise(n_patients = n(), share_shifted = mean(shifted), .groups = "drop"),
    cohort %>% mutate(shifted = flag) %>% group_by(modifier = "Discordance tertile", group = as.character(discordance_tertile)) %>%
      summarise(n_patients = n(), share_shifted = mean(shifted), .groups = "drop"),
    cohort %>% summarise(modifier = "All", group = "All", n_patients = n(), share_shifted = mean(flag))) %>%
    mutate(policy = policy)
}) %>% filter(n_patients >= MIN_CELL) %>% mutate(site = site_name)
write_csv(selection_table, out_file("click_selection"))

# =============================================================================
# 3. Fit the policies (lmtp, one time point)
# =============================================================================
age_basis    <- ns(cohort$age_at_admission / 10, df = 4)
height_basis <- ns(cohort$height_cm, df = 4)
model_data <- cohort %>%
  transmute(vt, Y = died_28, male = as.integer(tolower(sex_category) == "male"),
            race_category = factor(race_category), mode = factor(mode),
            sofa_total, log_sf = log(sf_ratio), bmi = if (USE_BMI) bmi else NA_real_,
            vt_floor, vt_cut_pfvc, vt_cut_pbw, strain_click_ml) %>%
  as.data.frame()
for (j in 1:4) {
  model_data[[paste0("age_ns", j)]] <- age_basis[, j]
  model_data[[paste0("height_ns", j)]] <- height_basis[, j]
  model_data[[paste0("height_ns", j, "_male")]] <- height_basis[, j] * model_data$male
}
if (!USE_BMI) model_data$bmi <- NULL
if (nlevels(model_data$mode) < 2) model_data$mode <- NULL
baseline_covariates <- setdiff(names(model_data), c("vt", "Y"))
message("W: ", paste(baseline_covariates, collapse = ", "))

fit_policy <- function(policy) {
  started <- Sys.time()
  fit <- lmtp_tmle(model_data, trt = "vt", outcome = "Y", baseline = baseline_covariates,
                   shift = if (policy == "natural") NULL else make_shift(policy), mtp = TRUE,
                   outcome_type = "binomial", learners_outcome = LEARNERS_OUTCOME,
                   learners_trt = LEARNERS_TRT, folds = FOLDS)
  message(sprintf("  %-10s fitted in %.1f min", policy, as.numeric(difftime(Sys.time(), started, units = "mins"))))
  fit
}
message("Fitting ", length(POLICY_LEVELS), " policies (", FOLDS, "-fold cross-fitting) ...")
fits <- setNames(map(POLICY_LEVELS, fit_policy), POLICY_LEVELS)

# =============================================================================
# 4. Positivity, the fallback rule, and the policy estimates
# =============================================================================
effective_sample_fraction <- function(ratio) (sum(ratio)^2 / sum(ratio^2)) / length(ratio)
positivity_table <- map_dfr(POLICY_LEVELS[-1], function(policy) {
  ratio <- as.numeric(fits[[policy]]$density_ratios[, 1])
  tibble(ratio, pfvc_tertile = as.character(cohort$pfvc_tertile)) %>%
    bind_rows(mutate(., pfvc_tertile = "All")) %>%
    group_by(pfvc_tertile) %>%
    summarise(n_patients = n(), share_ratio_gt10 = mean(ratio > 10), max_ratio = max(ratio),
              ess_fraction = effective_sample_fraction(ratio), .groups = "drop") %>%
    mutate(policy = policy)
}) %>% filter(n_patients >= MIN_CELL) %>% mutate(site = site_name)
write_csv(positivity_table, out_file("click_positivity"))

targeted_positivity <- positivity_table %>% filter(policy %in% c("click_pfvc", "click_pbw"), pfvc_tertile != "All")
positivity_ok <- all(targeted_positivity$share_ratio_gt10 <= POS_MAX_SHARE_GT10) &&
  all(targeted_positivity$ess_fraction >= POS_MIN_ESS_FRAC)
arms_distinct <- jaccard_arms <= OVERLAP_MAX_JACCARD
headline <- if (positivity_ok && arms_distinct) "head_to_head: click_pfvc - click_pbw" else
  "fallback: click_all - natural, with CATE by PFVC"
decision <- tibble(
  criterion = c("targeted arms: max share of density ratios > 10 across PFVC tertiles",
                "targeted arms: min effective-sample fraction across PFVC tertiles",
                "Jaccard overlap of the two arms' shifted patients", "headline"),
  value = c(sprintf("%.3f", max(targeted_positivity$share_ratio_gt10)),
            sprintf("%.3f", min(targeted_positivity$ess_fraction)), sprintf("%.3f", jaccard_arms), headline),
  threshold = c(sprintf("<= %.2f", POS_MAX_SHARE_GT10), sprintf(">= %.2f", POS_MIN_ESS_FRAC),
                sprintf("<= %.2f", OVERLAP_MAX_JACCARD), ""),
  passed = c(max(targeted_positivity$share_ratio_gt10) <= POS_MAX_SHARE_GT10,
             min(targeted_positivity$ess_fraction) >= POS_MIN_ESS_FRAC, arms_distinct,
             positivity_ok && arms_distinct),
  site = site_name)
write_csv(decision, out_file("click_decision"))

contrast_row <- function(fit, reference_fit, label) {
  est <- lmtp_contrast(fit, ref = reference_fit, type = "additive")$estimates
  tibble(contrast = label, risk_policy = est$shift, risk_reference = est$ref, rd = est$estimate,
         rd_lo = est$conf.low, rd_hi = est$conf.high, rd_se = est$std.error, p = est$p.value)
}
policy_summary <- map_dfr(POLICY_LEVELS[-1], function(policy) {
  flag <- cohort[[paste0("shifted_", policy)]]
  cut_ml <- click_size(policy, cohort)[flag]
  contrast_row(fits[[policy]], fits$natural, paste(policy, "- natural")) %>%
    mutate(policy = policy, share_shifted = mean(flag), n_shifted = sum(flag),
           mean_cut_ml_shifted = mean(cut_ml),
           mean_cut_vtpfvc_points_shifted = mean(cut_ml / cohort$pfvc[flag] * 0.1),
           mean_cut_vtpbw_mlkg_shifted = mean(cut_ml / cohort$pbw[flag]))
}) %>%
  bind_rows(contrast_row(fits$click_pfvc, fits$click_pbw, "click_pfvc - click_pbw (head-to-head)") %>%
              mutate(policy = "head_to_head")) %>%
  mutate(tau_pfvc = TAU_PFVC, tau_pbw = TAU_PBW, click_ml = CLICK_ML, strain_cut_points = STRAIN_CUT_POINTS,
         jaccard_arms = jaccard_arms,
         n_patients = nrow(cohort), headline = headline, site = site_name)
write_csv(policy_summary, out_file("click_policies"))
message("\nPolicy estimates (28-day mortality; negative RD = fewer deaths):")
print(as.data.frame(policy_summary %>% transmute(contrast, share_shifted = round(share_shifted, 3),
  rd_pp = round(100 * rd, 2), lo = round(100 * rd_lo, 2), hi = round(100 * rd_hi, 2))), row.names = FALSE)
message("Headline by the pre-specified rule: ", headline)

# =============================================================================
# 5. CATE by PFVC and by discordance (DR-learner), for the fixed and the strain-matched click
# =============================================================================
pseudo_outcome <- function(fit) fit$estimate@x + fit$estimate@eif
pseudo_natural <- pseudo_outcome(fits$natural)
RR_FLOOR <- 0.05   # the risk ratio is reported only where the fitted natural-course risk exceeds this

cate_by <- function(shift_policy, modifier_values, label) {
  pseudo_click <- pseudo_outcome(fits[[shift_policy]])
  d <- tibble(log_mod = log(modifier_values), effect = pseudo_click - pseudo_natural,
              click = pseudo_click, natural = pseudo_natural)
  basis <- ns(d$log_mod, df = 3)
  knots <- attr(basis, "knots"); boundary <- attr(basis, "Boundary.knots")
  spline_fit <- function(data, response)
    lm(as.formula(paste(response, "~ ns(log_mod, knots = knots, Boundary.knots = boundary)")), data = data)
  grid <- tibble(log_mod = seq(quantile(d$log_mod, 0.02), quantile(d$log_mod, 0.98), length.out = 60))
  ends <- tibble(log_mod = quantile(d$log_mod, c(0.1, 0.9)))
  curves_from <- function(data) {
    effect_fit <- spline_fit(data, "effect")
    risk_ratio <- predict(spline_fit(data, "click"), grid) / predict(spline_fit(data, "natural"), grid)
    risk_ratio[predict(spline_fit(data, "natural"), grid) <= RR_FLOOR] <- NA_real_
    ends_rr <- predict(spline_fit(data, "click"), ends) / predict(spline_fit(data, "natural"), ends)
    list(rd = predict(effect_fit, grid), rr = risk_ratio,
         gradient_rd = diff(predict(effect_fit, ends)), gradient_rr = ends_rr[2] / ends_rr[1])
  }
  point <- curves_from(d)
  boot <- map(seq_len(N_BOOT), ~ curves_from(d[sample.int(nrow(d), replace = TRUE), ]))
  band <- function(key, prob) apply(sapply(boot, `[[`, key), 1, quantile, prob, na.rm = TRUE)
  curve <- tibble(shift = shift_policy, modifier = label, value = exp(grid$log_mod),
                  rd = point$rd, rd_lo = band("rd", 0.025), rd_hi = band("rd", 0.975),
                  rr = point$rr, rr_lo = band("rr", 0.025), rr_hi = band("rr", 0.975))
  gradient_rd <- map_dbl(boot, "gradient_rd"); gradient_rr <- map_dbl(boot, "gradient_rr")
  slope <- tibble(shift = shift_policy, modifier = label,
    statistic = c("RD(p90) - RD(p10)", "RR(p90) / RR(p10)"),
    p10 = exp(ends$log_mod[1]), p90 = exp(ends$log_mod[2]),
    estimate = c(point$gradient_rd, point$gradient_rr),
    lo = c(quantile(gradient_rd, 0.025), quantile(gradient_rr, 0.025, na.rm = TRUE)),
    hi = c(quantile(gradient_rd, 0.975), quantile(gradient_rr, 0.975, na.rm = TRUE)))
  list(curve = curve, slope = slope)
}
message("CATE by PFVC and by discordance (", N_BOOT, " bootstrap replicates) ...")
cate_results <- list()
for (shift_policy in c("click_all", "click_strain")) {
  cate_results[[length(cate_results) + 1]] <- cate_by(shift_policy, cohort$pfvc, "PFVC (L)")
  cate_results[[length(cate_results) + 1]] <- cate_by(shift_policy, cohort$discordance, "PBW/PFVC discordance")
}
cate_curves <- map_dfr(cate_results, "curve") %>% mutate(site = site_name)
cate_slopes <- map_dfr(cate_results, "slope") %>% mutate(site = site_name)
write_csv(cate_curves, out_file("click_cate"))
write_csv(cate_slopes, out_file("click_cate_slope"))
message("CATE gradients (each click - natural; compare click_all with click_strain for dose heterogeneity):")
print(as.data.frame(cate_slopes %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)

# =============================================================================
# 6. Figure
# =============================================================================
theme_set(theme_minimal(base_size = 10))
panel_selection <- selection_table %>% filter(modifier == "PFVC tertile") %>%
  mutate(group = factor(group, c("Small PFVC", "Mid PFVC", "Large PFVC")),
         policy = factor(policy, POLICY_LEVELS)) %>%
  ggplot(aes(group, share_shifted, fill = policy)) + geom_col(position = "dodge") +
  scale_fill_manual(values = POLICY_COLOURS, name = NULL) +
  labs(title = "A. Who each arm shifts", subtitle = "share given one click lower, by PFVC tertile",
       x = NULL, y = "share shifted")

panel_effects <- policy_summary %>%
  mutate(contrast = fct_rev(factor(contrast, contrast))) %>%
  ggplot(aes(100 * rd, contrast)) + geom_vline(xintercept = 0, colour = "grey70") +
  geom_errorbar(aes(xmin = 100 * rd_lo, xmax = 100 * rd_hi), width = 0.2, orientation = "y", colour = okabe[5]) +
  geom_point(colour = okabe[5], size = 2) +
  labs(title = "B. 28-day mortality risk difference", subtitle = paste("headline:", headline),
       x = "risk difference (percentage points)", y = NULL)

panel_positivity <- positivity_table %>% filter(pfvc_tertile != "All") %>%
  mutate(pfvc_tertile = factor(pfvc_tertile, c("Small PFVC", "Mid PFVC", "Large PFVC")),
         policy = factor(policy, POLICY_LEVELS)) %>%
  ggplot(aes(pfvc_tertile, ess_fraction, fill = policy)) + geom_col(position = "dodge") +
  geom_hline(yintercept = POS_MIN_ESS_FRAC, linetype = 2) +
  scale_fill_manual(values = POLICY_COLOURS, name = NULL) +
  labs(title = "C. Positivity: effective-sample fraction", subtitle = "dashed = the fallback rule's floor",
       x = NULL, y = "ESS / n")

cate_panel <- function(curves, y, lo, hi, reference, ylab, title, log_y = FALSE) {
  plot <- ggplot(curves, aes(value, .data[[y]], colour = shift, fill = shift)) +
    geom_hline(yintercept = reference, colour = "grey70") +
    geom_ribbon(aes(ymin = .data[[lo]], ymax = .data[[hi]]), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.9) + facet_wrap(~ modifier, scales = "free_x") +
    scale_colour_manual(values = POLICY_COLOURS, name = NULL) +
    scale_fill_manual(values = POLICY_COLOURS, name = NULL) +
    labs(title = title, x = NULL, y = ylab)
  if (log_y) plot + scale_y_log10() else plot
}
panel_cate_rd <- cate_panel(cate_curves %>% mutate(across(c(rd, rd_lo, rd_hi), ~ 100 * .x)),
                            "rd", "rd_lo", "rd_hi", 0, "risk difference (pp)",
                            sprintf("D. CATE of one click lower for everyone: fixed %g mL vs strain-matched (%.2f points of VT/PFVC)",
                                    CLICK_ML, STRAIN_CUT_POINTS))
panel_cate_rr <- cate_panel(cate_curves, "rr", "rr_lo", "rr_hi", 1, "risk ratio",
                            "E. Same, relative scale (flat here + sloped above = baseline risk)", log_y = TRUE) +
  labs(caption = paste("A PFVC slope under the fixed click that flattens under the strain-matched click is dose",
                       "heterogeneity (strain is the dose); a slope that persists is effect modification."))

figure <- (panel_selection + panel_positivity) / panel_effects / panel_cate_rd / panel_cate_rr +
  plot_annotation(title = sprintf("One-click tidal-volume MTP at the index (%s%s)", site_name,
                                  if (is_synthetic) ", SYNTHETIC: simulated survival" else ""),
                  subtitle = sprintf("n = %d; PFVC arm VT/PFVC > %g%%, PBW arm VT/PBW > %.2f mL/kg (bite-matched); click %g mL",
                                     nrow(cohort), TAU_PFVC, TAU_PBW, CLICK_ML))
ggsave(out_file("click_mtp", "pdf"), figure, width = 11, height = 15)
message("16_vt_click_mtp complete -> ", final_dir)
