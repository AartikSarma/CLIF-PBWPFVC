# =============================================================================
# Script 05: PBW vs PFVC normalization of physiologic injury metrics
# Discordance + prognostic utility
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Pipeline script (run after 03; before 06). Reads the script 03 cross-sectional
# dataset, writes per-site outputs to output/<site>_output/final/ with a <site>
# suffix; script 06 discovers and pools them across cohorts (norm_* files).
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

output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
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

message("\nExploratory normalization discordance + prognostic utility analysis complete.")
