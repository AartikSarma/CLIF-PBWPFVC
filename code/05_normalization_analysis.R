# =============================================================================
# Script 05: PBW vs PFVC normalization of physiologic injury metrics
# Discordance + prognostic utility
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Pipeline script (run after 03). Reads the script 03 cross-sectional dataset and
# writes per-site tables to final/cross_sectional/ with a <site> suffix;
# pooled_estimates.R pools them across sites.
#
# For a multiplicative normalizer the discordance between the PBW- and the
# PFVC-normalized metric is exactly the size ratio PBW/PFVC, whatever the metric:
#    log(Ers x PFVC) - log(Ers x PBW) = -log(PBW/PFVC)
#    log(MP / PFVC)  - log(MP / PBW)  = +log(PBW/PFVC)
# The script reports two things about the normalizations of elastance (Goligher's
# Ers x PBW) and mechanical power (Gattinoni's MP/PBW):
#
#   PART 1 -- Discordance: the distribution of PBW/PFVC, and the share of patients
#     who change injury tertile when the normalizer is switched, overall and by age,
#     sex and race (Claim 3, supplement).
#   PART 2 -- Prognostic fit of each normalization for in-hospital death: the mechanic
#     alone, PBW-locked, PFVC-locked, and mechanic and size as free log terms, each
#     unadjusted and adjusted for age, sex and race (AIC, in-sample C, OR per log unit
#     and per SD).
#
# The rest of what this script once wrote (age-interaction ladders, encompassing
# tests, discrimination, value over driving pressure, prediction disagreement and
# calibration) was cut on 2026-09-24: no claim of the manuscript rests on it, and
# supplement/xsec_mortality_prediction.R and xsec_dp_vtpfvc_additive.R answer the
# prediction and driving-pressure questions (docs/output_manifest.md).
#
# QC: rows with dp <= 0 (plateau < PEEP; nonphysiologic measurement error) are
# dropped and counted.
#
# Input  : analysis_cross_sectional.parquet (script 03)
# Outputs: final/cross_sectional/norm_{discordance_summary, discordance_reclassification,
#          prognostic_fit, prognostic_coefs}_<site>.csv
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(broom)

source("utils/config.R")
site_name <- config$site_name

output_dir <- config$output_dir
final_dir  <- final_dir_for("cross_sectional")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

cross_sectional <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"))

# In-sample C-statistic (Mann-Whitney AUC). In-sample, so optimistic in absolute
# terms but fine for RELATIVE comparison of specs on the same data.
auc_fn <- function(y, p) {
  ok <- !is.na(y) & !is.na(p); y <- y[ok]; p <- p[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(p)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# -----------------------------------------------------------------------------
# QC + modelling frames
# -----------------------------------------------------------------------------
n_dp_bad <- cross_sectional %>%
  filter(!is.na(dp), dp <= 0) %>% nrow()
message("QC: ", n_dp_bad, " rows with dp <= 0 (plateau < PEEP) dropped.")

base <- cross_sectional %>%
  filter(!is.na(dp), dp > 0, !is.na(pfvc), pfvc > 0, !is.na(pbw), pbw > 0,
         !is.na(vtpbw), !is.na(bmi), !is.na(sofa_total), !is.na(sf_ratio),
         !is.na(deceased)) %>%
  mutate(
    sex_category  = factor(sex_category,  levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
    age10   = age_at_admission / 10,
    sf10    = sf_ratio / 10,
    pbwpfvc = pbw / pfvc                       # the size-estimate discordance factor
  )

ers_data <- base %>% filter(!is.na(ers), ers > 0) %>%
  mutate(ers_pbw = ers * pbw, ers_pfvc = ers * pfvc)
mp_data  <- base %>% filter(!is.na(mechanical_power), mechanical_power > 0,
                            !is.na(mp_pbw), !is.na(mp_pfvc))

message("Frames: base ", nrow(base), " | Ers ", nrow(ers_data),
        " | MP ", nrow(mp_data))

# =============================================================================
# PART 1 -- Physiologic discordance (the volume-estimate component = pbwpfvc)
# =============================================================================
# Discordance summary. PBW/PFVC is in kg/L, so its spread is reported as quantiles;
# the extremes stay at the site (no min or max: each is one patient's value).
disc_overall <- base %>%
  summarise(metric = "pbwpfvc (PBW/PFVC, kg/L)", n = n(),
            median = median(pbwpfvc), q25 = quantile(pbwpfvc, .25),
            q75 = quantile(pbwpfvc, .75), p10 = quantile(pbwpfvc, .10),
            p90 = quantile(pbwpfvc, .90))
write_csv(disc_overall,
          file.path(final_dir, paste0("norm_discordance_summary_", site_name, ".csv")))

message("\nPART 1 -- size-estimate discordance pbwpfvc: median ",
        round(disc_overall$median, 2), " kg/L (IQR ", round(disc_overall$q25, 2), "-",
        round(disc_overall$q75, 2), ")")

# Reclassification across injury tertiles when switching PBW -> PFVC normalization.
reclassify <- function(data, pbw_var, pfvc_var, label) {
  d2 <- data %>% filter(!is.na(.data[[pbw_var]]), !is.na(.data[[pfvc_var]]))
  t_pbw  <- dplyr::ntile(d2[[pbw_var]], 3)
  t_pfvc <- dplyr::ntile(d2[[pfvc_var]], 3)
  d2 <- d2 %>% mutate(reclassified = t_pbw != t_pfvc,
                      age_grp = cut(age_at_admission, c(0, 50, 65, 200),
                                    labels = c("<50", "50-65", ">65")))
  overall <- tibble(metric = label, group_type = "Overall", group = "All",
                    n = nrow(d2), pct_reclassified = mean(d2$reclassified) * 100)
  by_grp <- bind_rows(
    d2 %>% group_by(group = as.character(age_grp)) %>%
      summarise(n = n(), pct_reclassified = mean(reclassified) * 100, .groups = "drop") %>%
      mutate(metric = label, group_type = "Age"),
    d2 %>% group_by(group = as.character(sex_category)) %>%
      summarise(n = n(), pct_reclassified = mean(reclassified) * 100, .groups = "drop") %>%
      mutate(metric = label, group_type = "Sex"),
    d2 %>% group_by(group = as.character(race_category)) %>%
      summarise(n = n(), pct_reclassified = mean(reclassified) * 100, .groups = "drop") %>%
      mutate(metric = label, group_type = "Race")
  )
  bind_rows(overall, by_grp)
}
recl_tbl <- bind_rows(
  reclassify(ers_data, "ers_pbw", "ers_pfvc", "Normalized elastance (Goligher)"),
  reclassify(mp_data,  "mp_pbw",  "mp_pfvc",  "Mechanical power (Gattinoni)")
) %>%
  mutate(site = site_name, .before = 1)
write_csv(recl_tbl,
          file.path(final_dir, paste0("norm_discordance_reclassification_", site_name, ".csv")))

message("Reclassification across injury tertiles (PBW -> PFVC):")
recl_tbl %>% filter(group_type == "Overall") %>%
  pwalk(function(metric, pct_reclassified, ...)
    message("  ", metric, ": ", round(pct_reclassified, 1), "% reclassified"))

# =============================================================================
# PART 2 -- Prognostic utility of the normalization (in-hospital mortality)
# =============================================================================
base_cov <- "vtpbw + sofa_total + sf10 + bmi"   # vtpbw adjusted (dose); see memory
demo_cov <- "age10 + sex_category + race_category"

# 2x2 design (form x size) so the gain from "separate" can be split into model FORM
# (free vs locked size coefficient) and PHYSIOLOGY (PBW vs PFVC), plus the mechanic
# alone. The separate-PBW specs are the control: if Separate(.+PFVC) beats
# Separate(.+PBW), the improvement is physiologic, not just from freeing the form.
prog_specs <- tribble(
  ~family,           ~spec,                    ~exposure,
  "Elastance",       "Mechanic only",          "log(ers)",
  "Elastance",       "PBW-locked (Goligher)",  "log(ers_pbw)",
  "Elastance",       "PFVC-locked",            "log(ers_pfvc)",
  "Elastance",       "Separate (Ers + PBW)",   "log(ers) + log(pbw)",
  "Elastance",       "Separate (Ers + PFVC)",  "log(ers) + log(pfvc)",
  "Mechanical power","Mechanic only",          "log(mechanical_power)",
  "Mechanical power","PBW-locked (Gattinoni)", "log(mp_pbw)",
  "Mechanical power","PFVC-locked",            "log(mp_pfvc)",
  "Mechanical power","Separate (MP + PBW)",    "log(mechanical_power) + log(pbw)",
  "Mechanical power","Separate (MP + PFVC)",   "log(mechanical_power) + log(pfvc)"
)
family_data <- list("Elastance" = ers_data, "Mechanical power" = mp_data)

fit_prog <- function(family, spec, exposure, adjusted) {
  data <- family_data[[family]]
  cov  <- if (adjusted) paste(base_cov, "+", demo_cov) else base_cov
  m <- glm(as.formula(paste("deceased ~", exposure, "+", cov)),
           data = data, family = binomial)
  list(meta = tibble(family = family, spec = spec,
                     adjusted = if (adjusted) "adjusted" else "unadjusted",
                     aic = AIC(m), auc = auc_fn(data$deceased, fitted(m)),
                     n = stats::nobs(m)),
       model = m)
}

prog_fits <- pmap(crossing(prog_specs, adjusted = c(FALSE, TRUE)),
                  function(family, spec, exposure, adjusted)
                    fit_prog(family, spec, exposure, adjusted))
prog_tbl <- bind_rows(map(prog_fits, "meta")) %>%
  group_by(family, adjusted) %>%
  mutate(delta_aic = aic - min(aic),
         spec = factor(spec, levels = unique(prog_specs$spec))) %>%
  ungroup() %>% mutate(site = site_name, .before = 1)
write_csv(prog_tbl,
          file.path(final_dir, paste0("norm_prognostic_fit_", site_name, ".csv")))

# Effect sizes behind the fit: the mortality OR (point estimate + 95% CI) for each
# normalization's exposure term(s), at BOTH adjustment levels. The AIC/AUC above say
# which spec fits best; these give the actual association with uncertainty so the
# point estimates and CIs can be pooled and plotted across cohorts. ORs are reported
# per log-unit (the model coefficient) and per 1 SD of the log-exposure (comparable
# across normalizations, since log(ers_pbw) and log(ers_pfvc) have different spreads).
z975 <- qnorm(0.975)                              # Wald CIs (stable/fast at these N)
prog_coefs <- map_dfr(prog_fits, function(f) {
  d <- family_data[[f$meta$family]]
  broom::tidy(f$model) %>%
    filter(str_detect(term, "^log\\(")) %>%       # exposure terms only (not covariates)
    rowwise() %>%
    mutate(sd_log = sd(log(d[[gsub("^log\\((.*)\\)$", "\\1", term)]]), na.rm = TRUE)) %>%
    ungroup() %>%
    transmute(family = f$meta$family, spec = as.character(f$meta$spec),
              adjusted = f$meta$adjusted, term,
              or_per_log    = exp(estimate),
              or_per_log_lo = exp(estimate - z975 * std.error),
              or_per_log_hi = exp(estimate + z975 * std.error),
              or_per_sd     = exp(estimate * sd_log),
              or_per_sd_lo  = exp((estimate - z975 * std.error) * sd_log),
              or_per_sd_hi  = exp((estimate + z975 * std.error) * sd_log),
              std_error = std.error, n = f$meta$n)
}) %>% mutate(site = site_name, .before = 1)
write_csv(prog_coefs,
          file.path(final_dir, paste0("norm_prognostic_coefs_", site_name, ".csv")))

message("\nPART 2 -- prognostic fit (lower AIC = better; dAIC vs best within family x adjustment):")
prog_tbl %>% arrange(family, adjusted, delta_aic) %>%
  pwalk(function(family, spec, adjusted, delta_aic, auc, ...)
    message(sprintf("  [%-15s | %-10s] %-24s dAIC = %5.1f  C = %.3f",
                    family, adjusted, spec, delta_aic, auc)))
