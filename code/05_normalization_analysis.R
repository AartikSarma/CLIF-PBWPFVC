# =============================================================================
# Script 05: PBW vs PFVC normalization of physiologic injury metrics
# Discordance + prognostic utility
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Pipeline script (run after 03). Reads the script 03 cross-sectional dataset,
# writes per-site outputs to output/<site>_output/final/ with a <site> suffix;
# pooled_estimates.R discovers and pools them across cohorts (norm_* files).
#
# Prognostic superiority of PFVC for VT dosing is already established (Sarma LRM
# 2025 + the analyses here). The next question is whether PFVC should also replace
# PBW in the PHYSIOLOGY metrics whose landmark papers use PBW (or no size
# reference): Goligher's normalized elastance (Ers x PBW), Gattinoni's MP/PBW, and
# Amato's driving pressure. This script asks two things for the size-normalized
# metrics (Ers, MP):
#
#   PART 1 -- Physiologic discordance: how different is the PBW- vs PFVC-normalized
#     metric, and for whom? For any multiplicative normalizer the discordance is
#     EXACTLY the size ratio pbwpfvc = PBW/PFVC (independent of the metric):
#        log(Ers x PFVC) - log(Ers x PBW) = -log(pbwpfvc)
#        log(MP / PFVC)  - log(MP / PBW)  = +log(pbwpfvc)
#     so the discordance is the published demographic size bias. We quantify its
#     distribution, its demographic patterning, and the resulting reclassification
#     across injury tertiles.
#
#     CAVEAT (two error sources): a pressure-based metric is distorted by BOTH (a)
#     errors in the pre-injury lung-volume estimate -- which switching PBW->PFVC
#     corrects (this is the pbwpfvc discordance) -- and (b) changes in elastic
#     recoil, especially with age, which affect the measured Ers/DP itself and are
#     NOT fixed by any size normalization. The pbwpfvc discordance here is the
#     VOLUME-estimate component; the recoil component is a separate axis (shown vs
#     age) addressed only by the adjusted models and the age-referenced analyses.
#
#   PART 2 -- Prognostic utility of the normalization (in-hospital mortality). Per
#     metric, compare: mechanic alone; PBW-locked composite (the published form);
#     PFVC-locked composite; and SEPARATE (mechanic + PFVC as free log terms).
#     Earlier model fits favored keeping the mechanic and PFVC SEPARATE over locking
#     them into a single specific-metric product, so the locked-vs-separate LRT is a
#     headline. Compared by AIC and C-statistic, EACH in demographic-unadjusted and
#     adjusted form (age/sex/race are PFVC's parents; the unadjusted model is policy-
#     relevant, the adjusted is the identification check). Also the demographic-
#     invariance test per locked normalization (the universal one is unmodified).
#
# QC: rows with dp <= 0 (plateau < PEEP; nonphysiologic measurement error) are
# dropped and counted. Kept in this exploratory script for now.
#
# Inputs : output/<site>/intermediate/analysis_cross_sectional.parquet (script 03)
# Outputs: final/norm_discordance_*<site>.csv  and  final/norm_discordance_*<site>.pdf
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(broom)
library(splines)
library(patchwork)

source("utils/config.R")
site_name <- config$site_name

output_dir <- config$output_dir
final_dir  <- config$final_dir
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

cross_sectional <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"))

zscore <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

# In-sample C-statistic (Mann-Whitney AUC). In-sample, so optimistic in absolute
# terms but fine for RELATIVE comparison of specs on the same data.
auc_fn <- function(y, p) {
  ok <- !is.na(y) & !is.na(p); y <- y[ok]; p <- p[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(p)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Per-log-unit OR of an exposure by age (delta method), for a model with the
# exposure and its `:age10` interaction. age10 here is age/10 (uncentered), so the
# curve is evaluated across the grid; age_years = age10*10. Resolves reversed names.
ratio_by_age <- function(model, main_term, int_term, age_grid) {
  b <- coef(model); V <- vcov(model)
  if (!(int_term %in% names(b))) {
    alt <- paste(rev(strsplit(int_term, ":")[[1]]), collapse = ":")
    if (alt %in% names(b)) int_term <- alt else stop("interaction not found: ", int_term)
  }
  lp <- b[[main_term]] + b[[int_term]] * age_grid
  se <- sqrt(V[main_term, main_term] + age_grid^2 * V[int_term, int_term] +
               2 * age_grid * V[main_term, int_term])
  tibble(age_years = age_grid * 10, ratio = exp(lp),
         conf_low = exp(lp - 1.96 * se), conf_high = exp(lp + 1.96 * se))
}

# Max per-term VIF of a model's design matrix (the full pairwise model approaches
# the factorial, so flag if collinearity is getting large).
vif_max <- function(model) {
  X <- model.matrix(model)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  max(vapply(seq_len(ncol(X)), function(j) {
    r2 <- suppressWarnings(summary(lm(X[, j] ~ X[, -j, drop = FALSE]))$r.squared)
    if (!is.finite(r2) || r2 >= 1) Inf else 1 / (1 - r2)
  }, numeric(1)))
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
# Discordance summary overall and by demographic group.
disc_overall <- base %>%
  summarise(metric = "pbwpfvc (PBW/PFVC)", n = n(),
            median = median(pbwpfvc), q25 = quantile(pbwpfvc, .25),
            q75 = quantile(pbwpfvc, .75), min = min(pbwpfvc), max = max(pbwpfvc),
            pct_PBW_over_15 = mean(pbwpfvc > 1.15) * 100,
            pct_PBW_under_85 = mean(pbwpfvc < 0.85) * 100)
write_csv(disc_overall,
          file.path(final_dir, paste0("norm_discordance_summary_", site_name, ".csv")))

# Systematic discordance: log(pbwpfvc) ~ demographics (the published size bias).
disc_model <- lm(log(pbwpfvc) ~ age10 + sex_category + race_category, data = base)
disc_coef <- broom::tidy(disc_model, conf.int = TRUE) %>%
  mutate(site = site_name, .before = 1)
write_csv(disc_coef,
          file.path(final_dir, paste0("norm_discordance_demographics_", site_name, ".csv")))

message("\nPART 1 -- size-estimate discordance pbwpfvc: median ",
        round(disc_overall$median, 3), " (IQR ", round(disc_overall$q25, 3), "-",
        round(disc_overall$q75, 3), "); ", round(disc_overall$pct_PBW_over_15, 1),
        "% with PBW >15% over PFVC.")

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
) %>% mutate(site = site_name, .before = 1)
write_csv(recl_tbl,
          file.path(final_dir, paste0("norm_discordance_reclassification_", site_name, ".csv")))

message("Reclassification across injury tertiles (PBW -> PFVC):")
recl_tbl %>% filter(group_type == "Overall") %>%
  pwalk(function(metric, pct_reclassified, ...)
    message("  ", metric, ": ", round(pct_reclassified, 1), "% reclassified"))

# Figure: discordance vs age (volume-estimate error grows as PFVC declines with
# age) and reclassification by demographic group.
p_disc_age <- ggplot(base, aes(age_at_admission, pbwpfvc)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
  geom_point(alpha = 0.12, color = "#0072B2") +
  geom_smooth(method = "loess", formula = y ~ x, color = "#D55E00", fill = "#D55E00") +
  labs(title = "Size-estimate discordance vs age",
       subtitle = paste0(site_name,
         " - PBW/PFVC > 1 means PBW overestimates lung size; the volume-estimate ",
         "error grows with age (recoil error is a separate, uncorrected axis)."),
       x = "Age (years)", y = "PBW / PFVC") +
  theme_minimal(base_size = 11)
p_recl <- recl_tbl %>% filter(group_type != "Overall") %>%
  ggplot(aes(pct_reclassified, group, fill = metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_grid(group_type ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("Normalized elastance (Goligher)" = "#009E73",
                               "Mechanical power (Gattinoni)" = "#E69F00"), name = NULL) +
  labs(title = "Injury-tertile reclassification when switching PBW -> PFVC",
       subtitle = paste0(site_name, " - % of patients changing injury tertile, by group."),
       x = "% reclassified", y = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(final_dir, paste0("norm_discordance_part1_", site_name, ".pdf")),
       p_disc_age / p_recl + patchwork::plot_layout(heights = c(1, 1.2)),
       width = 9, height = 9)

# =============================================================================
# PART 1b -- Direct age-interaction tests of the discordance
# =============================================================================
# (i) Normalization x age (paired). For a patient, log(metric_PBW) -
# log(metric_PFVC) = log(pbwpfvc) exactly, so regressing log(pbwpfvc) on age IS the
# normalization x age interaction. A nonzero slope means the PBW-vs-PFVC discordance
# scales with age -- older patients are most mis-sized by PBW. A spline checks
# whether the age dependence is nonlinear.
m_disc_lin <- lm(log(pbwpfvc) ~ age10, data = base)
m_disc_spl <- lm(log(pbwpfvc) ~ ns(age10, 3), data = base)
disc_spl_lrt <- anova(m_disc_lin, m_disc_spl)
disc_age_tbl <- broom::tidy(m_disc_lin, conf.int = TRUE) %>%
  filter(term == "age10") %>%
  transmute(site = site_name,
            test = "Discordance ~ age (normalization x age, paired)",
            per_decade_pct = (exp(estimate) - 1) * 100,
            lo_pct = (exp(conf.low) - 1) * 100, hi_pct = (exp(conf.high) - 1) * 100,
            p_value = p.value, nonlinearity_p = disc_spl_lrt$`Pr(>F)`[2])

# (ii) Heterogeneity: does the age-discordance slope differ by sex / race?
m_het0 <- lm(log(pbwpfvc) ~ age10 + sex_category + race_category, data = base)
m_het  <- lm(log(pbwpfvc) ~ age10 * sex_category + age10 * race_category, data = base)
het_lrt <- anova(m_het0, m_het)
het_tbl <- broom::tidy(m_het, conf.int = TRUE) %>%
  filter(grepl("age10:", term)) %>%
  transmute(site = site_name, term,
            per_decade_pct = (exp(estimate) - 1) * 100,
            lo_pct = (exp(conf.low) - 1) * 100, hi_pct = (exp(conf.high) - 1) * 100,
            p_value = p.value) %>%
  mutate(joint_heterogeneity_p = het_lrt$`Pr(>F)`[2])

# (iii) Reclassification x age (per metric): does the probability of changing injury
# tertile when switching PBW -> PFVC rise with age?
recl_flag <- function(data, pbw_var, pfvc_var) {
  mutate(data, reclassified = as.integer(
    dplyr::ntile(.data[[pbw_var]], 3) != dplyr::ntile(.data[[pfvc_var]], 3)))
}
recl_age_fit <- list(
  "Normalized elastance (Goligher)" = glm(
    reclassified ~ age10 + sex_category + race_category,
    data = recl_flag(ers_data, "ers_pbw", "ers_pfvc"), family = binomial),
  "Mechanical power (Gattinoni)" = glm(
    reclassified ~ age10 + sex_category + race_category,
    data = recl_flag(mp_data, "mp_pbw", "mp_pfvc"), family = binomial))
recl_age_tbl <- imap_dfr(recl_age_fit, function(m, lab)
  broom::tidy(m, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "age10") %>%
    transmute(site = site_name, metric = lab, reclass_OR_per_decade = estimate,
              lo = conf.low, hi = conf.high, p_value = p.value))

write_csv(disc_age_tbl,
          file.path(final_dir, paste0("norm_discordance_age_", site_name, ".csv")))
write_csv(het_tbl,
          file.path(final_dir, paste0("norm_discordance_age_heterogeneity_", site_name, ".csv")))
write_csv(recl_age_tbl,
          file.path(final_dir, paste0("norm_reclassification_age_", site_name, ".csv")))

message("\nPART 1b -- direct age-interaction tests:")
message(sprintf("  Discordance per decade: %+.1f%% (p = %s); nonlinearity p = %s",
                disc_age_tbl$per_decade_pct, signif(disc_age_tbl$p_value, 2),
                signif(disc_age_tbl$nonlinearity_p, 2)))
message("  Age-slope heterogeneity by sex/race: joint p = ",
        signif(het_lrt$`Pr(>F)`[2], 3))
recl_age_tbl %>% pwalk(function(metric, reclass_OR_per_decade, p_value, ...)
  message(sprintf("  Reclassification OR/decade [%s]: %.2f (p = %s)",
                  metric, reclass_OR_per_decade, signif(p_value, 2))))

# Figure 1b: discordance-by-age slope by sex (heterogeneity) + predicted
# reclassification probability vs age per metric.
age_grid <- tibble(age10 = seq(quantile(base$age10, 0.05, na.rm = TRUE),
                               quantile(base$age10, 0.95, na.rm = TRUE), length.out = 40))
het_pred <- crossing(age_grid, sex_category = levels(base$sex_category)) %>%
  mutate(sex_category  = factor(sex_category, levels = levels(base$sex_category)),
         race_category = factor("WHITE", levels = levels(base$race_category)))
het_pred$fit       <- exp(predict(m_het, newdata = het_pred))
het_pred$age_years <- het_pred$age10 * 10
p_het <- ggplot(het_pred, aes(age_years, fit, color = sex_category)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c("Male" = "#0072B2", "Female" = "#D55E00"), name = NULL) +
  labs(subtitle = "Discordance (PBW/PFVC) by age and sex",
       x = "Age (years)", y = "PBW / PFVC (fitted)") +
  theme_minimal(base_size = 11)

recl_pred <- imap_dfr(recl_age_fit, function(m, lab) {
  nd <- age_grid %>% mutate(
    sex_category = factor("Male", levels = levels(base$sex_category)),
    race_category = factor("WHITE", levels = levels(base$race_category)))
  pr <- predict(m, newdata = nd, type = "link", se.fit = TRUE)
  tibble(metric = lab, age_years = nd$age10 * 10, p = plogis(pr$fit),
         lo = plogis(pr$fit - 1.96 * pr$se.fit), hi = plogis(pr$fit + 1.96 * pr$se.fit))
})
p_recl_age <- ggplot(recl_pred, aes(age_years, p, color = metric, fill = metric)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c("Normalized elastance (Goligher)" = "#009E73",
                                "Mechanical power (Gattinoni)" = "#E69F00"), name = NULL) +
  scale_fill_manual(values = c("Normalized elastance (Goligher)" = "#009E73",
                               "Mechanical power (Gattinoni)" = "#E69F00"), guide = "none") +
  labs(subtitle = "P(reclassified across injury tertile) vs age",
       x = "Age (years)", y = "Predicted P(reclassified)") +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(final_dir, paste0("norm_discordance_part1b_", site_name, ".pdf")),
       p_het | p_recl_age, width = 11, height = 5)

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

# Form vs physiology decomposition. "Separate" may fit better for reasons of model
# FORM (a free size coefficient) rather than because PFVC is the better size
# reference. Disentangle within each family x adjustment:
#   FORM effect (within a size reference) = locked vs separate -- nested LRT
#     (locked = separate with the two log-coefficients constrained equal).
#   PHYSIOLOGY effect (within a form) = PBW vs PFVC -- non-nested, same df, so
#     compared by AIC. The decisive control is Separate(.+PBW) vs Separate(.+PFVC):
#     if PFVC wins THERE, the gain is physiologic, not just from freeing the form.
get_model <- function(fam, sp, adj) {
  i <- which(map_lgl(prog_fits, ~ .x$meta$family == fam & .x$meta$spec == sp &
                       .x$meta$adjusted == (if (adj) "adjusted" else "unadjusted")))
  prog_fits[[i]]$model
}
decomp_tbl <- crossing(family = c("Elastance", "Mechanical power"),
                       adjusted = c(FALSE, TRUE)) %>%
  pmap_dfr(function(family, adjusted) {
    nm <- if (family == "Elastance")
      c(pbw_lock = "PBW-locked (Goligher)", pfvc_lock = "PFVC-locked",
        pbw_sep = "Separate (Ers + PBW)", pfvc_sep = "Separate (Ers + PFVC)")
    else
      c(pbw_lock = "PBW-locked (Gattinoni)", pfvc_lock = "PFVC-locked",
        pbw_sep = "Separate (MP + PBW)", pfvc_sep = "Separate (MP + PFVC)")
    m <- lapply(nm, function(s) get_model(family, s, adjusted))
    form_pbw  <- anova(m$pbw_lock,  m$pbw_sep,  test = "LRT")
    form_pfvc <- anova(m$pfvc_lock, m$pfvc_sep, test = "LRT")
    tibble(site = site_name, family = family,
           adjusted = if (adjusted) "adjusted" else "unadjusted",
           form_pbw_lrt_p  = form_pbw$`Pr(>Chi)`[2],   # locked vs separate, PBW
           form_pfvc_lrt_p = form_pfvc$`Pr(>Chi)`[2],  # locked vs separate, PFVC
           phys_locked_dAIC   = AIC(m$pbw_lock) - AIC(m$pfvc_lock),  # >0 = PFVC better
           phys_separate_dAIC = AIC(m$pbw_sep)  - AIC(m$pfvc_sep),   # >0 = PFVC better (key)
           aic_pbw_sep = AIC(m$pbw_sep), aic_pfvc_sep = AIC(m$pfvc_sep))
  })
write_csv(decomp_tbl,
          file.path(final_dir, paste0("norm_form_vs_physiology_", site_name, ".csv")))
message("\nForm vs physiology (is the gain PFVC, or just the separate form?):")
decomp_tbl %>% pwalk(function(family, adjusted, form_pbw_lrt_p, form_pfvc_lrt_p,
                              phys_locked_dAIC, phys_separate_dAIC, ...)
  message(sprintf(paste0("  [%-15s | %-10s] form LRT p PBW/PFVC = %s / %s | ",
                         "physiology dAIC locked %+.1f, separate %+.1f (>0 = PFVC better)"),
                  family, adjusted, signif(form_pbw_lrt_p, 2), signif(form_pfvc_lrt_p, 2),
                  phys_locked_dAIC, phys_separate_dAIC)))

# Demographic invariance: is each locked normalization's mortality effect modified
# by age/sex/race? The normalization that is INVARIANT is the universal dosing/
# assessment target. LRT of the exposure x (age+sex+race) interaction block.
invar_specs <- tribble(~family, ~exposure, ~label,
  "Elastance", "log(ers_pbw)",  "Ers x PBW (Goligher)",
  "Elastance", "log(ers_pfvc)", "Ers x PFVC",
  "Mechanical power", "log(mp_pbw)",  "MP/PBW (Gattinoni)",
  "Mechanical power", "log(mp_pfvc)", "MP/PFVC")
invar_tbl <- pmap_dfr(invar_specs, function(family, exposure, label) {
  data <- family_data[[family]]
  m0 <- glm(as.formula(paste("deceased ~", exposure, "+", base_cov, "+", demo_cov)),
            data = data, family = binomial)
  m1 <- glm(as.formula(paste0("deceased ~ ", exposure, " * (", demo_cov, ") + ", base_cov)),
            data = data, family = binomial)
  lrt <- anova(m0, m1, test = "LRT")
  tibble(site = site_name, family = family, normalization = label,
         demo_interaction_df = lrt$Df[2], demo_interaction_p = lrt$`Pr(>Chi)`[2])
})
write_csv(invar_tbl,
          file.path(final_dir, paste0("norm_demographic_invariance_", site_name, ".csv")))
message("\nDemographic invariance (interaction LRT; NS = universal normalization):")
invar_tbl %>% pwalk(function(normalization, demo_interaction_p, ...)
  message("  ", normalization, ": demo-interaction p = ", signif(demo_interaction_p, 3),
          if (!is.na(demo_interaction_p) && demo_interaction_p > 0.05)
            "  -> demographically invariant" else "  -> demographically modified"))

# Figure: prognostic fit (dAIC) per spec, faceted family x adjustment; plus C-stat.
p_aic <- ggplot(prog_tbl, aes(delta_aic, spec, color = adjusted)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  facet_wrap(~ family, scales = "free_y") +
  scale_color_manual(values = c("unadjusted" = "#0072B2", "adjusted" = "#D55E00"),
                     name = NULL) +
  labs(title = "Prognostic fit of normalization choices (lower dAIC = better)",
       subtitle = paste0(site_name,
         " - dAIC vs best spec within family x adjustment, in-hospital mortality. ",
         "'Separate' frees the mechanic and PFVC; 'locked' is the published composite."),
       x = "delta AIC (vs best in group)", y = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
p_auc <- ggplot(prog_tbl, aes(auc, spec, color = adjusted)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  facet_wrap(~ family, scales = "free_y") +
  scale_color_manual(values = c("unadjusted" = "#0072B2", "adjusted" = "#D55E00"),
                     guide = "none") +
  labs(subtitle = "C-statistic (in-sample; relative comparison only)",
       x = "C-statistic", y = NULL) +
  theme_minimal(base_size = 11)
ggsave(file.path(final_dir, paste0("norm_prognostic_fit_", site_name, ".pdf")),
       p_aic / p_auc + patchwork::plot_layout(heights = c(1, 1)),
       width = 10, height = 9)

# =============================================================================
# PART 2b -- Age interaction with the mechanic (DP / Ers / MP)
# =============================================================================
# Respiratory mechanics change with age (recoil declines), so the mortality effect
# of a given DP / Ers / MP may itself be age-dependent. For each mechanic, fit a
# nested ladder (VT/PBW adjusted throughout; unadjusted + demographic-adjusted):
#   a: log(X)                            (mechanic only)
#   b: log(X) + age10                    (age main effect)
#   c: log(X) * age10                    (mechanic x age)
#   d: log(X) * age10 + log(pfvc)        (+ PFVC kept separate)
#   e: log(X) * age10 + log(pfvc)*age10  (+ PFVC x age)
# and test: age_interaction (c vs b), PFVC adds beyond age-modified mechanic
# (d vs c), and PFVC x age (e vs d). DP (Amato) is included here even though it has
# no PBW/PFVC-normalized form -- the PFVC question for DP is whether available
# volume modifies it, which these specs capture.
# Ladder rungs. Beyond the age-interaction rungs, two test a pressure x lung-size
# interaction (does the harm of the pressure depend on the volume it is distributed
# over): sep_int adds it age-naively, d_int adds it on top of the age-modified
# model. d_int is the joint two-modifier model (pressure modified by BOTH age/recoil
# and lung volume) -- d_int vs d is the decisive test that volume modifies the
# pressure effect beyond recoil/age.
# Rungs use CENTERED predictors (mech_c = log(X)-mean, pfvc_c = log(pfvc)-mean,
# age_c = age10-6 [age 60]); fit_ladder builds them per subset. Centering removes
# the product-with-component VIF artifact in the interaction rungs (LRTs and AIC are
# centering-invariant) and makes "hold the other modifier at its mean" = set it to 0.
# full = all three 2-way interactions (no 3-way), so recoil (mech:age), PFVC-age
# (pfvc:age) and volume (mech:pfvc) modifications coexist and the volume term is not
# confounded by an omitted PFVC:age effect.
ladder_rhs <- c(
  a       = "mech_c",
  b       = "mech_c + age_c",
  c       = "mech_c * age_c",
  sep     = "mech_c + pfvc_c",
  sep_int = "mech_c + pfvc_c + mech_c:pfvc_c",
  d       = "mech_c * age_c + pfvc_c",
  d_int   = "mech_c * age_c + pfvc_c + mech_c:pfvc_c",
  e       = "mech_c * age_c + pfvc_c * age_c",
  full    = "mech_c * age_c + pfvc_c * age_c + mech_c:pfvc_c")

exps <- tribble(~family, ~mx, ~data_name,
  "Driving pressure", "dp",               "base",
  "Elastance",        "ers",              "ers_data",
  "Mechanical power", "mechanical_power", "mp_data")
data_map <- list(base = base, ers_data = ers_data, mp_data = mp_data)

fit_ladder <- function(family, mx, data, adjusted) {
  # Centered predictors (see ladder_rhs). age_c uses age10 - 6 so age_c is THE age
  # variable in both the interactions and the demographic adjustment (avoids a
  # collinear age10/age_c pair); logs centered at their means.
  d2 <- data %>% mutate(
    mech_c = log(.data[[mx]]) - mean(log(.data[[mx]]), na.rm = TRUE),
    pfvc_c = log(pfvc) - mean(log(pfvc), na.rm = TRUE),
    age_c  = age10 - 6)
  cov <- if (adjusted)
    paste(base_cov, "+ age_c + sex_category + race_category") else base_cov
  ms  <- lapply(ladder_rhs, function(r)
    glm(as.formula(paste("deceased ~", r, "+", cov)), data = d2, family = binomial))
  lrtp <- function(s, b) anova(ms[[s]], ms[[b]], test = "LRT")$`Pr(>Chi)`[2]
  meta <- tibble(family = family, mx = mx,
                 adjusted = if (adjusted) "adjusted" else "unadjusted",
                 aic_mech = AIC(ms$a), aic_mech_x_age = AIC(ms$c),
                 aic_age_pfvc = AIC(ms$d), aic_joint = AIC(ms$d_int),
                 aic_full = AIC(ms$full),
                 age_interaction_p = lrtp("b", "c"),     # mechanic x age
                 pfvc_adds_p = lrtp("c", "d"),           # PFVC main adds beyond age-mod
                 pfvc_x_age_p = lrtp("d", "e"),          # PFVC x age
                 size_int_alone_p = lrtp("sep", "sep_int"),  # pressure x size, age-naive
                 size_int_joint_p = lrtp("d", "d_int"),  # pressure x size beyond mechanic x age
                 size_int_full_p = lrtp("e", "full"),    # pressure x size beyond BOTH age interactions
                 pfvc_age_in_full_p = lrtp("d_int", "full"),  # PFVC x age beyond pressure x size
                 max_vif_full = vif_max(ms$full), n = stats::nobs(ms$a))
  list(meta = meta, joint_model = ms$full,
       mean_log_pfvc = mean(log(data$pfvc), na.rm = TRUE))
}
ladder_in   <- crossing(exps, adjusted = c(FALSE, TRUE))
ladder_fits <- pmap(ladder_in, function(family, mx, data_name, adjusted)
  fit_ladder(family, mx, data_map[[data_name]], adjusted))
ladder_tbl  <- bind_rows(map(ladder_fits, "meta")) %>% mutate(site = site_name, .before = 1)
write_csv(ladder_tbl,
          file.path(final_dir, paste0("norm_age_interaction_ladder_", site_name, ".csv")))

message("\nPART 2b -- age + size interaction with the mechanic (LRT p; 'full' = all 2-way):")
ladder_tbl %>% pwalk(function(family, adjusted, age_interaction_p, pfvc_x_age_p,
                              size_int_joint_p, size_int_full_p, max_vif_full, ...)
  message(sprintf(paste0("  [%-15s | %-10s] mech x age = %s | PFVC x age = %s | ",
                         "press x size beyond mech-age = %s | beyond BOTH age ints = %s | max VIF = %.1f"),
                  family, adjusted, signif(age_interaction_p, 2), signif(pfvc_x_age_p, 2),
                  signif(size_int_joint_p, 2), signif(size_int_full_p, 2), max_vif_full)))

# Per-log-unit mechanic OR as a linear combination of coefficients across a grid,
# for the joint model (log(X)*age10 + log(pfvc) + log(X):log(pfvc)): the slope wrt
# log(X) is b_logX + b_{X:age}*age10 + b_{X:size}*log(pfvc). wfun(g) returns the
# (age, size) weights at grid point g; delta-method CI over the three coefficients.
or_combo <- function(model, base, t_age, t_size, grid, wfun) {
  b <- coef(model); V <- vcov(model)
  rsv <- function(t) { if (t %in% names(b)) return(t)
    alt <- paste(rev(strsplit(t, ":")[[1]]), collapse = ":")
    if (alt %in% names(b)) alt else stop("term not found: ", t) }
  nm <- c(rsv(base), rsv(t_age), rsv(t_size))
  rows <- lapply(grid, function(g) {
    w <- c(1, wfun(g))
    est <- sum(w * b[nm]); v <- as.numeric(crossprod(w, V[nm, nm] %*% w))
    c(est = est, se = sqrt(v))
  })
  est <- vapply(rows, `[[`, numeric(1), "est"); se <- vapply(rows, `[[`, numeric(1), "se")
  tibble(grid = grid, ratio = exp(est),
         conf_low = exp(est - 1.96 * se), conf_high = exp(est + 1.96 * se))
}
exp_levels <- c("Driving pressure", "Elastance", "Mechanical power")
# Centered terms: base mech_c, interactions mech_c:age_c and mech_c:pfvc_c. The
# other modifier is centered, so holding it "at its mean / age 60" = setting it to 0.
age_grid2 <- seq(-3.0, 2.5, length.out = 40)   # age_c (age 30-85; age = age_c*10 + 60)

# (i) Mechanic OR BY AGE, at mean lung size (pfvc_c = 0), from the joint model.
ladder_age <- imap_dfr(ladder_fits, function(f, i) {
  if (!ladder_in$adjusted[i]) return(NULL)
  or_combo(f$joint_model, "mech_c", "mech_c:age_c", "mech_c:pfvc_c",
           age_grid2, function(ac) c(ac, 0)) %>%
    transmute(exposure = f$meta$family, x = grid * 10 + 60, ratio, conf_low, conf_high)
}) %>% mutate(exposure = factor(exposure, levels = exp_levels))

# (ii) Mechanic OR BY LUNG SIZE, at age 60 (age_c = 0), from the joint model.
ladder_pfvc <- imap_dfr(ladder_fits, function(f, i) {
  if (!ladder_in$adjusted[i]) return(NULL)
  data <- data_map[[ladder_in$data_name[i]]]; mlp <- f$mean_log_pfvc
  pv <- seq(quantile(data$pfvc, 0.05, na.rm = TRUE),
            quantile(data$pfvc, 0.95, na.rm = TRUE), length.out = 40)
  or_combo(f$joint_model, "mech_c", "mech_c:age_c", "mech_c:pfvc_c",
           log(pv) - mlp, function(pc) c(0, pc)) %>%
    transmute(exposure = f$meta$family, x = exp(grid + mlp), ratio, conf_low, conf_high)
}) %>% mutate(exposure = factor(exposure, levels = exp_levels))

write_csv(bind_rows(mutate(ladder_age, axis = "Age (years)"),
                    mutate(ladder_pfvc, axis = "Predicted FVC (L)")),
          file.path(final_dir, paste0("norm_age_interaction_slopes_", site_name, ".csv")))

mod_curve <- function(df, xlab, line_col) {
  ggplot(df, aes(x, ratio)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.15, fill = line_col) +
    geom_line(linewidth = 0.9, color = line_col) +
    facet_wrap(~ exposure, scales = "free_y") +
    scale_y_log10() +
    labs(x = xlab, y = "Mortality OR per log-unit (log scale)") +
    theme_minimal(base_size = 11)
}
p_age  <- mod_curve(ladder_age, "Age (years)", "#0072B2") +
  labs(subtitle = "Recoil axis: mechanic OR by age (lung size at mean)")
p_pfvc <- mod_curve(ladder_pfvc, "Predicted FVC (L)", "#E69F00") +
  labs(subtitle = "Volume axis: mechanic OR by lung size (age 60)")
p_ladder <- (p_age / p_pfvc) +
  patchwork::plot_annotation(
    title = "Does the mechanic's mortality effect change with recoil (age) and lung volume?",
    subtitle = paste0(site_name,
      " - per-log-unit mortality OR from the joint model (X*age + PFVC + X:PFVC, ",
      "VT/PBW + demographic adjusted). Non-flat = the mechanic is modified by that axis."))
ggsave(file.path(final_dir, paste0("norm_age_interaction_slopes_", site_name, ".pdf")),
       p_ladder, width = 11, height = 8)

# =============================================================================
# PART 2c -- Encompassing test: does PFVC add information that PBW misses?
# =============================================================================
# Separate models (mechanic + PBW vs mechanic + PFVC) are NON-nested -- AIC says
# which fits better, not that one carries information the other misses. The clean
# test is incremental (nested): put both size references in one model and test each
# direction. PFVC strictly dominates if it adds beyond PBW (LRT significant) while
# PBW adds nothing beyond PFVC (NS). Equivalence: since log(pfvc) = log(pbw) -
# log(pbwpfvc), "PFVC beyond PBW" is the same test as adding the PBW/PFVC discordance
# term -- i.e. does the PBW sizing error carry prognostic information. PBW and PFVC
# are collinear, so rely on the block LRT + dAIC, not individual coefficients
# (max VIF of the both-in model is reported). VT/PBW adjusted; +/- demographics.
encompassing <- function(family, mx, data, adjusted) {
  cov <- if (adjusted)
    paste(base_cov, "+ age10 + sex_category + race_category") else base_cov
  f <- function(rhs) glm(as.formula(paste("deceased ~", rhs, "+", cov)),
                         data = data, family = binomial)
  m_pbw  <- f(sprintf("log(%s) + log(pbw)", mx))
  m_pfvc <- f(sprintf("log(%s) + log(pfvc)", mx))
  m_both <- f(sprintf("log(%s) + log(pbw) + log(pfvc)", mx))
  lrt_pfvc <- anova(m_pbw,  m_both, test = "LRT")   # PFVC beyond PBW
  lrt_pbw  <- anova(m_pfvc, m_both, test = "LRT")   # PBW  beyond PFVC
  tibble(site = site_name, family = family,
         adjusted = if (adjusted) "adjusted" else "unadjusted",
         pfvc_beyond_pbw_p = lrt_pfvc$`Pr(>Chi)`[2],
         pbw_beyond_pfvc_p = lrt_pbw$`Pr(>Chi)`[2],
         daic_pfvc_beyond_pbw = AIC(m_pbw)  - AIC(m_both),  # >0 = adding PFVC improves
         daic_pbw_beyond_pfvc = AIC(m_pfvc) - AIC(m_both),  # >0 = adding PBW improves
         max_vif_both = vif_max(m_both), n = stats::nobs(m_both))
}
enc_tbl <- pmap_dfr(crossing(exps, adjusted = c(FALSE, TRUE)),
                    function(family, mx, data_name, adjusted)
                      encompassing(family, mx, data_map[[data_name]], adjusted))
write_csv(enc_tbl, file.path(final_dir, paste0("norm_encompassing_", site_name, ".csv")))

message("\nPART 2c -- encompassing test (does PFVC add info PBW misses?):")
enc_tbl %>% pwalk(function(family, adjusted, pfvc_beyond_pbw_p, pbw_beyond_pfvc_p,
                           daic_pfvc_beyond_pbw, daic_pbw_beyond_pfvc, max_vif_both, ...)
  message(sprintf(paste0("  [%-15s | %-10s] PFVC beyond PBW: p=%s dAIC=%+.1f | ",
                         "PBW beyond PFVC: p=%s dAIC=%+.1f | max VIF=%.1f"),
                  family, adjusted, signif(pfvc_beyond_pbw_p, 2), daic_pfvc_beyond_pbw,
                  signif(pbw_beyond_pfvc_p, 2), daic_pbw_beyond_pfvc, max_vif_both)))

# Figure: dAIC from adding the second size reference, each direction. The asymmetry
# -- PFVC beyond PBW positive/large, PBW beyond PFVC near zero -- is the result.
enc_long <- enc_tbl %>%
  transmute(family = factor(family, levels = exp_levels), adjusted,
            `PFVC beyond PBW` = daic_pfvc_beyond_pbw,
            `PBW beyond PFVC` = daic_pbw_beyond_pfvc) %>%
  pivot_longer(c(`PFVC beyond PBW`, `PBW beyond PFVC`),
               names_to = "direction", values_to = "dAIC")
p_enc <- ggplot(enc_long, aes(dAIC, family, fill = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_wrap(~ adjusted) +
  scale_fill_manual(values = c("PFVC beyond PBW" = "#009E73",
                               "PBW beyond PFVC" = "#D55E00"), name = NULL) +
  labs(title = "Does PFVC add prognostic information that PBW misses?",
       subtitle = paste0(site_name,
         " - improvement in AIC from adding the second size reference to mechanic + ",
         "the first (>0 = improves fit). PFVC adds beyond PBW; PBW does not beyond PFVC."),
       x = "delta AIC from adding the size term (>0 = improves fit)", y = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(final_dir, paste0("norm_encompassing_", site_name, ".pdf")),
       p_enc, width = 10, height = 4.5)

# =============================================================================
# PART 2d -- Discrimination by normalization (optimism-corrected univariate C)
# =============================================================================
# Which size scale lets each mechanic best discriminate mortality? For DP, MP and
# Ers, compare raw vs PBW- vs PFVC-normalized on the optimism-corrected C-statistic
# (AUC) for in-hospital mortality. Logistic/in-hospital analogue of the Cox version
# previously in script 06 (which has been retired); the ordering -- PFVC-normalized
# discriminating best -- is the result. DP-normalized columns (DP/PBW, DP/PFVC) are
# derived here as dp / pbw and dp / pfvc.
B_DISC <- if (identical(site_name, "synthetic_clif")) 100L else 300L
auc_optimism <- function(d, var, B) {           # apparent AUC minus mean bootstrap optimism
  f   <- reformulate(var, "deceased")
  app <- auc_fn(d$deceased, fitted(glm(f, data = d, family = binomial)))
  opt <- numeric(0)
  for (b in seq_len(B)) {
    bd <- d[sample.int(nrow(d), replace = TRUE), , drop = FALSE]
    m  <- tryCatch(glm(f, data = bd, family = binomial), error = function(e) NULL)
    if (is.null(m)) next
    opt <- c(opt, auc_fn(bd$deceased, fitted(m)) -
                  auc_fn(d$deceased, predict(m, newdata = d, type = "response")))
  }
  o <- mean(opt, na.rm = TRUE)        # robust to occasional NA-AUC resamples
  if (!is.finite(o)) o <- 0
  tibble(c_apparent = app, optimism = o, c_corrected = app - o)
}
dp_disc <- base %>% mutate(dp_pbw = dp / pbw, dp_pfvc = dp / pfvc)
disc_reg <- tribble(
  ~family,            ~norm,  ~var,               ~metric,
  "Driving pressure", "raw",  "dp",               "DP",
  "Driving pressure", "pbw",  "dp_pbw",           "DP / PBW",
  "Driving pressure", "pfvc", "dp_pfvc",          "DP / PFVC",
  "Mechanical power", "raw",  "mechanical_power", "MP",
  "Mechanical power", "pbw",  "mp_pbw",           "MP / PBW",
  "Mechanical power", "pfvc", "mp_pfvc",          "MP / PFVC",
  "Elastance",        "raw",  "ers",              "Ers",
  "Elastance",        "pbw",  "ers_pbw",          "Ers x PBW",
  "Elastance",        "pfvc", "ers_pfvc",         "Ers x PFVC")
disc_data <- list("Driving pressure" = dp_disc, "Mechanical power" = mp_data, "Elastance" = ers_data)
set.seed(20260617)
disc_tbl <- pmap_dfr(disc_reg, function(family, norm, var, metric) {
  d <- disc_data[[family]] %>% filter(is.finite(.data[[var]]))
  bind_cols(tibble(site = site_name, family = family, norm = norm, metric = metric),
            auc_optimism(d, var, B_DISC), tibble(n = nrow(d)))
})
write_csv(disc_tbl, file.path(final_dir, paste0("norm_discrimination_", site_name, ".csv")))

message("\nPART 2d -- optimism-corrected discrimination (C for in-hospital mortality):")
disc_tbl %>% arrange(family, norm) %>%
  pwalk(function(family, metric, c_corrected, ...)
    message(sprintf("  [%-16s] %-12s C = %.3f", family, metric, c_corrected)))
best_norm <- disc_tbl %>% group_by(family) %>% slice_max(c_corrected, n = 1) %>% ungroup()
message("  best-discriminating normalization per family: ",
        paste(sprintf("%s=%s", best_norm$family, best_norm$norm), collapse = "; "))

# =============================================================================
# PART 2e -- Does a PFVC-normalized metric add prognostic value OVER driving pressure?
# =============================================================================
# Driving pressure (Amato) is the bedside-validated injury metric, and it is ALREADY partly
# size-aware (dP = VT / compliance, and compliance scales with lung size -- so a small lung
# shows a high dP at a "normal" VT/PBW). The clinical skeptic's question is therefore the
# pressure analogue of "is PFVC just age": is PFVC just driving pressure? For each PFVC-
# normalized metric we fit dP alone vs dP + the metric (both keep the standard covariates,
# incl. VT/PBW dose), and report the discrimination gain, the metric's mortality OR per SD
# NET of dP, and the nested LRT -- with and without demographic adjustment.
#   metric adds over dP  => size-relative strain/power captures injury dP misses.
#   metric adds nothing  => dP already captures it; PFVC's edge is then PRACTICAL -- it is
#                           computable from routine settings + height when no plateau exists
#                           (dP needs a measured plateau: intermittent, sedation-dependent, MNAR).
# NB: normalized ELASTANCE is deliberately NOT tested incrementally over dP -- dP IS
# mechanically VT x Ers, so the two carry the same information and the test is collinear
# by construction (its PFVC-vs-PBW question is answered in the PART 2 head-to-head instead).
# Raw MP and MP/PBW included alongside MP/PFVC: under the identifiability argument the
# incremental-over-dP signal is the POWER itself (demographic-independent: rate, PEEP,
# pressures), so adjusted, raw/PBW/PFVC should add over dP by ~the same -- the normalizer
# (a demographic function) is absorbed by adjustment, just like in the multiplicative shift.
# The unadjusted MP/PFVC vs MP/PBW gap, if any, is the discordance (demographic) encoding.
incr_specs <- tribble(
  ~metric_label, ~metric_var,        ~data_name,
  "VT/PFVC",     "vtpfvc",           "base",
  "MP (raw)",    "mechanical_power", "mp_data",
  "MP/PBW",      "mp_pbw",           "mp_data",
  "MP/PFVC",     "mp_pfvc",          "mp_data")
data_map_incr <- list(base = base, mp_data = mp_data, ers_data = ers_data)
incr_fn <- function(metric_label, metric_var, data_name, adjusted) {
  d <- data_map_incr[[data_name]] %>%
    filter(is.finite(dp), dp > 0, is.finite(.data[[metric_var]]), .data[[metric_var]] > 0)
  cov <- if (adjusted) paste(base_cov, "+", demo_cov) else base_cov
  mx  <- paste0("log(", metric_var, ")")
  m0  <- glm(as.formula(paste("deceased ~ log(dp) +", cov)), data = d, family = binomial)
  m1  <- glm(as.formula(paste("deceased ~ log(dp) +", mx, "+", cov)), data = d, family = binomial)
  co  <- summary(m1)$coefficients
  trm <- rownames(co)[grepl(metric_var, rownames(co), fixed = TRUE)][1]
  sdl <- sd(log(d[[metric_var]]), na.rm = TRUE)
  est <- co[trm, "Estimate"]; se <- co[trm, "Std. Error"]
  tibble(site = site_name, metric = metric_label,
         adjusted = if (adjusted) "adjusted" else "unadjusted",
         c_dp_alone = auc_fn(d$deceased, fitted(m0)),
         c_dp_plus_metric = auc_fn(d$deceased, fitted(m1)),
         c_gain = auc_fn(d$deceased, fitted(m1)) - auc_fn(d$deceased, fitted(m0)),
         metric_or_per_sd = exp(est * sdl),
         or_lo = exp((est - z975 * se) * sdl), or_hi = exp((est + z975 * se) * sdl),
         lrt_p = anova(m0, m1, test = "LRT")$`Pr(>Chi)`[2], n = nrow(d))
}
incr_tbl <- pmap_dfr(crossing(incr_specs, adjusted = c(FALSE, TRUE)),
  function(metric_label, metric_var, data_name, adjusted)
    incr_fn(metric_label, metric_var, data_name, adjusted))
write_csv(incr_tbl, file.path(final_dir, paste0("norm_dp_incremental_", site_name, ".csv")))
message("\nPART 2e -- incremental value of a PFVC metric OVER driving pressure (C + OR/SD net of dP):")
incr_tbl %>% arrange(metric, adjusted) %>%
  pwalk(function(metric, adjusted, c_dp_alone, c_gain, metric_or_per_sd, lrt_p, ...)
    message(sprintf("  [%-10s | %-10s] C(dP)=%.3f  +metric => +%.3f  OR/SD=%.2f  LRT p=%s",
                    metric, adjusted, c_dp_alone, c_gain, metric_or_per_sd, signif(lrt_p, 2))))

# =============================================================================
# PART 2f -- Partial-adjustment ladder + E-values: identifiability vs attenuation
# =============================================================================
# PFVC is a deterministic function of (height, age, sex, race), so adjusting for all of its
# demographic parents removes its variation and its independent effect becomes UNIDENTIFIABLE
# -- which is NOT the same as null. This walks the adjustment up one parent at a time (age,
# +sex, +race) and reads the SIGNATURE:
#   ATTENUATION   -- OR drifts to 1 with a TIGHT CI  => the effect was demographic confounding.
#   DESTABILIZATION-- OR stays away from 1 but the SE/CI EXPLODE (se_inflation rises) => the
#                     effect is UNIDENTIFIABLE (collinearity removed the variation), not absent.
# IMPORTANT: destabilization signals unidentifiability, which is AGNOSTIC to whether the true
# effect is real or null -- both a real-but-collinear effect and a null-but-collinear quantity
# blow the SE up. PBW/PFVC (the discordance) is included as the NEGATIVE CONTROL: it is a pure
# measurement-mismatch ratio of demographic functions with NO independent lung effect, yet it
# should destabilize just like VT/PFVC -- proving the signature alone cannot certify "real".
# What distinguishes VT/PFVC from PBW/PFVC is physiology + the instrument/trial, not this table.
# VT/PFVC vs MP/PFVC/raw MP (which keep demographic-independent variation and stay identified)
# is the contrast. NOTE: height (a parent) is only partly captured here via VT/PBW and BMI.
eval_or <- function(or) {                          # VanderWeele E-value, common-outcome OR (~ sqrt(OR) -> RR)
  o <- if (is.finite(or) && or >= 1) or else if (is.finite(or) && or > 0) 1 / or else return(NA_real_)
  rr <- sqrt(o); if (!is.finite(rr) || rr < 1) return(1); rr + sqrt(rr * (rr - 1))
}
ladder_steps <- list(base = "", "+age" = "+ age10", "+age+sex" = "+ age10 + sex_category",
                     "+age+sex+race" = "+ age10 + sex_category + race_category")
ladder_fn <- function(metric_label, metric_var, data_name) {
  d <- data_map_incr[[data_name]] %>% filter(is.finite(.data[[metric_var]]), .data[[metric_var]] > 0)
  sdl <- sd(log(d[[metric_var]]), na.rm = TRUE); mx <- paste0("log(", metric_var, ")")
  imap_dfr(ladder_steps, function(extra, step) {
    m  <- glm(as.formula(paste("deceased ~", mx, "+", base_cov, extra)), data = d, family = binomial)
    co <- summary(m)$coefficients[mx, ]; est <- co["Estimate"]; se <- co["Std. Error"]
    or <- exp(est * sdl); lo <- exp((est - z975 * se) * sdl); hi <- exp((est + z975 * se) * sdl)
    tibble(metric = metric_label, adjust = step, or_per_sd = or, lo = lo, hi = hi, se_log = se,
           evalue_point = eval_or(or),
           evalue_ci = if (lo > 1) eval_or(lo) else if (hi < 1) eval_or(hi) else 1)
  }) %>% mutate(se_inflation = se_log / se_log[adjust == "base"])
}
ladder_specs <- bind_rows(incr_specs,   # + PBW/PFVC discordance as the null-physiology negative control
  tibble(metric_label = "PBW/PFVC", metric_var = "pbwpfvc", data_name = "base"))
ladder_tbl <- pmap_dfr(ladder_specs, ladder_fn) %>% mutate(site = site_name, .before = 1)
write_csv(ladder_tbl, file.path(final_dir, paste0("norm_partial_adjust_evalue_", site_name, ".csv")))
message("\nPART 2f -- partial-adjustment ladder + E-values (SE-inflation = collinearity; CI-explode w/o OR->1 = unidentifiable, not null):")
ladder_tbl %>% pwalk(function(metric, adjust, or_per_sd, lo, hi, se_inflation, evalue_point, evalue_ci, ...)
  message(sprintf("  [%-9s | %-13s] OR/SD=%.2f [%.2f,%.2f]  SEx%.1f  E(pt)=%.2f  E(CI)=%.2f",
                  metric, adjust, or_per_sd, lo, hi, se_inflation, evalue_point, evalue_ci)))

message("\nExploratory normalization discordance + prognostic utility analysis complete.")

# =============================================================================
# PART 3 -- Do the two normalizations make DIFFERENT statements about the lung,
#           and is either of them a universal scale?
# =============================================================================
# PART 2 asked which normalization predicts death better and found that, once age,
# sex and race are adjusted, they are indistinguishable. That is a question about
# prognosis. This part asks the MEASUREMENT question the paper's second pillar is
# built on: Ers x PBW and Ers x PFVC (and MP/PBW vs MP/PFVC) are two scales for the
# same physiology, they differ by exactly the PBW/PFVC discordance, and a clinician
# reading one reaches a different conclusion about the same patient than one reading
# the other. Two reads:
#
#   3a. HOW MUCH they disagree, and WHERE. The two scales have different units, so
#       they are compared on the only common footing that needs no cohort-specific
#       rescaling of one into the other: the patient's PERCENTILE within the cohort
#       on each scale. The disagreement is the percentile shift, it is reported
#       against discordance and by demographic group, and the reclassification is
#       the share of patients who cross the top-tertile "high normalized elastance"
#       line when the normalizer is swapped. This is the consequence of the bias for
#       measurement, independent of any outcome.
#
#   3b. WHICH scale is universal, if either. A normalized measurement claims that a
#       given value means the same thing in every patient. Test it: fit mortality on
#       the normalized measure and severity ONLY (no demographics -- adjusting for
#       them would absorb the discordance and force the answer), then compare the
#       observed with the predicted mortality WITHIN discordance strata. A normalizer
#       that sizes the lung correctly is calibrated across those strata; one that
#       mis-sizes systematically under- or over-predicts where mis-sizing is worst.
#       Unlike the adjusted AIC/AUC comparisons in PART 2, this can separate them.
norm_families <- tribble(
  ~family,             ~frame,      ~pbw_var,  ~pfvc_var, ~direction,
  "Elastance",         "ers_data",  "ers_pbw", "ers_pfvc", "higher = stiffer",
  "Mechanical power",  "mp_data",   "mp_pbw",  "mp_pfvc",  "higher = more power")

disc_cuts <- function(x) cut(x, c(-Inf, quantile(x, c(1/3, 2/3)), Inf),
                            labels = c("Concordant", "Mid", "Discordant"))
MIN_CELL <- 10   # small-cell suppression for every stratum reported below (CLAUDE.md)

pred_disagree <- list(); pred_recl <- list(); pred_calib <- list()
for (i in seq_len(nrow(norm_families))) {
  fam <- norm_families$family[i]
  d   <- get(norm_families$frame[i]) %>%
    mutate(x_pbw = .data[[norm_families$pbw_var[i]]], x_pfvc = .data[[norm_families$pfvc_var[i]]]) %>%
    filter(is.finite(x_pbw), x_pbw > 0, is.finite(x_pfvc), x_pfvc > 0) %>%
    mutate(disc_grp = disc_cuts(pbwpfvc),
           pct_pbw  = 100 * percent_rank(x_pbw),
           pct_pfvc = 100 * percent_rank(x_pfvc),
           shift    = pct_pfvc - pct_pbw,               # + = the PFVC scale ranks this patient higher
           hi_pbw   = x_pbw  > quantile(x_pbw,  2/3),   # "high normalized elastance/power"
           hi_pfvc  = x_pfvc > quantile(x_pfvc, 2/3))
  if (nrow(d) < 200) next

  # --- 3a. disagreement, overall / by discordance decile / by demographic group ----
  summ <- function(g) summarise(g, n = n(), median_shift = median(shift),
                                median_abs_shift = median(abs(shift)),
                                p90_abs_shift = quantile(abs(shift), .9),
                                frac_shift_over_20 = mean(abs(shift) > 20), .groups = "drop")
  pred_disagree[[fam]] <- bind_rows(
    d %>% mutate(stratum_type = "Overall", stratum = "All") %>% group_by(stratum_type, stratum) %>% summ(),
    d %>% mutate(stratum_type = "Discordance decile",
                 stratum = as.character(ntile(pbwpfvc, 10))) %>% group_by(stratum_type, stratum) %>% summ(),
    d %>% mutate(stratum_type = "Sex",  stratum = as.character(sex_category)) %>% group_by(stratum_type, stratum) %>% summ(),
    d %>% mutate(stratum_type = "Race", stratum = as.character(race_category)) %>% group_by(stratum_type, stratum) %>% summ(),
    d %>% mutate(stratum_type = "Age group",
                 stratum = as.character(cut(age_at_admission, c(-Inf, 50, 65, Inf),
                                            labels = c("<50", "50-64", ">=65")))) %>%
      group_by(stratum_type, stratum) %>% summ()) %>%
    filter(n >= MIN_CELL) %>% mutate(family = fam, .before = 1)
  # how strongly the shift tracks discordance (it should: the scales differ BY discordance)
  pred_disagree[[fam]]$shift_vs_disc_r2 <- summary(lm(shift ~ log(pbwpfvc), data = d))$r.squared

  # --- 3a. reclassification across the top-tertile decision line -------------------
  pred_recl[[fam]] <- bind_rows(
    d %>% mutate(stratum_type = "Overall", stratum = "All"),
    d %>% mutate(stratum_type = "Discordance tertile", stratum = as.character(disc_grp)),
    d %>% mutate(stratum_type = "Sex",  stratum = as.character(sex_category)),
    d %>% mutate(stratum_type = "Race", stratum = as.character(race_category))) %>%
    group_by(stratum_type, stratum) %>%
    summarise(n = n(), frac_high_pbw = mean(hi_pbw), frac_high_pfvc = mean(hi_pfvc),
              frac_reclassified = mean(hi_pbw != hi_pfvc),
              frac_pfvc_only = mean(hi_pfvc & !hi_pbw), frac_pbw_only = mean(hi_pbw & !hi_pfvc),
              .groups = "drop") %>%
    filter(n >= MIN_CELL) %>% mutate(family = fam, .before = 1)

  # --- 3b. is either scale universal? calibration across discordance strata ---------
  # Severity-adjusted, demographics-FREE: the normalized measure is asked to carry the
  # size information on its own, which is what "normalized" claims.
  calib_one <- function(var, lab) {
    f <- glm(as.formula(paste("deceased ~ log(", var, ") + vtpbw + sofa_total + sf10 + bmi")),
             data = d, family = binomial)
    dd <- d %>% mutate(p = predict(f, type = "response"))
    # per stratum: observed vs expected, and the calibration intercept (logit offset
    # needed to correct the stratum; 0 = calibrated, >0 = the model UNDER-predicts risk)
    dd %>% group_by(disc_grp) %>%
      group_modify(~ {
        # a stratum with almost no events cannot be calibrated: the offset model runs
        # off to +-20 with an infinite CI, so report the counts and leave the intercept NA
        n_ev <- sum(.x$deceased == 1)
        base_out <- tibble(n = nrow(.x), events = n_ev, observed = mean(.x$deceased),
                           expected = mean(.x$p),
                           obs_over_exp = if (mean(.x$p) > 0) mean(.x$deceased) / mean(.x$p) else NA_real_)
        if (n_ev < MIN_CELL || n_ev == nrow(.x))
          return(bind_cols(base_out, tibble(calib_intercept = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_)))
        off <- qlogis(pmin(pmax(.x$p, 1e-6), 1 - 1e-6))
        m   <- glm(deceased ~ 1, offset = off, data = .x, family = binomial)
        ci  <- suppressMessages(confint.default(m))
        bind_cols(base_out, tibble(calib_intercept = unname(coef(m)[1]), ci_lo = ci[1, 1], ci_hi = ci[1, 2]))
      }) %>% ungroup() %>% mutate(normalizer = lab, .before = 1)
  }
  pred_calib[[fam]] <- bind_rows(calib_one("x_pbw", "PBW-normalized"),
                                 calib_one("x_pfvc", "PFVC-normalized")) %>%
    mutate(family = fam, .before = 1)
}
disagree_tbl <- bind_rows(pred_disagree) %>% mutate(site = site_name)
recl_tbl     <- bind_rows(pred_recl)     %>% mutate(site = site_name)
calib_tbl    <- bind_rows(pred_calib)    %>% mutate(site = site_name)
write_csv(disagree_tbl, file.path(final_dir, paste0("norm_prediction_disagreement_", site_name, ".csv")))
write_csv(recl_tbl,     file.path(final_dir, paste0("norm_prediction_reclassification_", site_name, ".csv")))
write_csv(calib_tbl,    file.path(final_dir, paste0("norm_prediction_calibration_", site_name, ".csv")))

message("\n=== PART 3a: the two normalizations rank the same patient differently ===")
print(as.data.frame(disagree_tbl %>% filter(stratum_type %in% c("Overall", "Discordance decile")) %>%
        transmute(family, stratum_type, stratum, n, median_shift = round(median_shift, 1),
                  median_abs = round(median_abs_shift, 1), pct_over_20 = round(100 * frac_shift_over_20))),
      row.names = FALSE)
message("=== PART 3a: who crosses the 'high' line when the normalizer is swapped ===")
print(as.data.frame(recl_tbl %>% transmute(family, stratum_type, stratum, n,
        pct_reclassified = round(100 * frac_reclassified), pct_pfvc_only = round(100 * frac_pfvc_only),
        pct_pbw_only = round(100 * frac_pbw_only))), row.names = FALSE)
message("=== PART 3b: is either scale universal? observed/expected within discordance strata ===")
message("    (severity-adjusted, NO demographics; O/E far from 1 = that normalizer mis-sizes there)")
print(as.data.frame(calib_tbl %>% transmute(family, normalizer, disc_grp, n, events,
        obs = round(100 * observed, 1), exp = round(100 * expected, 1),
        o_over_e = round(obs_over_exp, 3),
        calib_intercept = ifelse(is.na(calib_intercept), "(too few events)",
                                 sprintf("%+.2f [%.2f, %.2f]", calib_intercept, ci_lo, ci_hi)))), row.names = FALSE)
