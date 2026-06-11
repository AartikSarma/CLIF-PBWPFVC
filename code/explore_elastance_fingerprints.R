# =============================================================================
# Exploratory: statistical fingerprints of the two errors in normalized elastance
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Exploratory, standalone analysis (NOT part of the 00 pipeline runner). Heavily
# caveated -- see below.
#
# "Elastance normalized to predicted size" (e.g. Ers x PBW) tries to recover the
# intrinsic specific lung elastance E_spec,L. For the respiratory system:
#       E_rs = E_cw + E_L = E_cw + E_spec,L / V0
# (E_spec,L = intrinsic recoil, V0 = aerated resting lung volume, E_cw = chest
# wall). "E_rs x size" assumes (1) E_spec,L constant [recoil error] and (2) the
# surrogate is proportional to V0 [size-surrogate error], and ignores E_cw. We
# cannot measure E_spec,L or V0 here, so the two errors cannot be separated
# absolutely -- but they leave different fingerprints.
#
# KEY FRAMING (corrected): "flatter normalized elastance vs demographics" is NOT a
# clean criterion, because AGE is a common cause of BOTH errors -- it lowers V0
# (size) AND lowers E_spec,L (recoil loss with aging), which push E_rs in opposite
# directions. Writing g_X = d(log X)/d(age):
#       g_Ers = g_Espec - g_V0,    g_(Ers x S) = g_Espec - g_V0 + g_S
#   * PBW is age-independent (g_PBW = 0): Ers x PBW flat vs age means g_Espec ~=
#     g_V0, i.e. the recoil decline and the size decline CANCEL -- flatness is the
#     artifact, not a virtue.
#   * PFVC declines with age (g_PFVC ~= g_V0): Ers x PFVC ~ g_Espec, so it should
#     DECLINE with age, recovering the expected recoil loss.
# So age structure is confounded and must be read against the recoil prior, not
# minimized. The clean separation uses axes that load on only one error:
#   * SEX, RACE, HEIGHT move predicted SIZE but not intrinsic recoil -> residual
#     structure here is a near-pure size-surrogate signal (error 2); flatter after
#     PFVC IS meaningful on these axes.
#   * OXYGENATION / SOFA, and predicted FEV1/FVC, load on RECOIL/derangement, not
#     predicted size -> near-pure recoil signal (error 1). Predicted FEV1/FVC is a
#     demographic estimate of intrinsic recoil (it falls with age as elastic recoil
#     is lost): lower predicted FEV1/FVC = less recoil = lower specific elastance.
#
# CAVEATS: surrogate-on-surrogate (predicted, not measured, volumes; no esophageal
# pressure, so chest wall -- which STIFFENS with age, partly offsetting lung recoil
# loss in airway-measured E_rs -- is folded in; SF ratio is a crude aerated-fraction
# proxy). Predicted FVC and predicted FEV1/FVC are both deterministic functions of
# the same demographics, so they are collinear: their relative loading shows whether
# the demographic structure in E_rs looks more "size-shaped" or "recoil-shaped", but
# does not add information beyond age/sex/height/race. This is a variance/association
# partition, NOT a per-patient deconfound.
#
# Inputs : output/<site>/intermediate/analysis_cross_sectional.parquet (script 03)
# Outputs: final/elastance_fingerprints_summary_<site>.csv  (poolable statistics)
#          final/elastance_fingerprints_size_<site>.pdf      (size/demographic + age panels)
#          final/elastance_fingerprints_recoil_<site>.pdf    (recoil + variance partition)
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(patchwork)
library(broom)
library(rspiro)

source("utils/config.R")
site_name <- config$site_name

output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

cross_sectional <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"))

MIN_N <- 10  # do not report any model fit on fewer than this many subjects

# Mechanics subset: elastance (and its normalizations) require a measured plateau.
# Predicted FEV1/FVC (a demographic estimate of intrinsic recoil) is computed the
# same GLI-2012 way as PFVC. sex_numeric / race_numeric are recomputed here so the
# script does not depend on them surviving the parquet round-trip.
mech <- cross_sectional %>%
  filter(is.finite(ers), is.finite(ers_pbw), is.finite(ers_pfvc)) %>%
  mutate(
    sex_category  = factor(sex_category,  levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
    sex_numeric   = if_else(sex_category == "Male", 1L, 2L),
    race_numeric  = case_when(race_category == "WHITE" ~ 1L,
                              race_category == "BLACK" ~ 2L, TRUE ~ 5L),
    age10    = age_at_admission / 10,
    height10 = height_cm / 10,
    sf10     = sf_ratio / 10,
    # Predicted recoil surrogate: GLI-2012 predicted FEV1/FVC ratio.
    pred_fev1fvc = pred_GLI(age = age_at_admission, height = height_cm / 100,
                            gender = sex_numeric, ethnicity = race_numeric,
                            param = "FEV1FVC")
  )

message("Subjects with a measured respiratory-system elastance: ", nrow(mech))
if (nrow(mech) < MIN_N) {
  stop("Only ", nrow(mech), " subjects have a measured elastance (need >= ", MIN_N,
       "). Too few to explore the elastance fingerprints.")
}

# Helpers.
r2 <- function(formula, data) summary(lm(as.formula(formula), data = data))$r.squared
slope <- function(formula, data, term) coef(lm(as.formula(formula), data = data))[[term]]

results <- list()  # collected scalar findings -> tidy summary CSV

# =============================================================================
# (A) Size axis: which predicted-size surrogate explains elastance variation?
# =============================================================================
size_data <- mech %>% filter(is.finite(pbw), is.finite(pfvc))
r2_size_pbw  <- r2("ers ~ pbw",  size_data)
r2_size_pfvc <- r2("ers ~ pfvc", size_data)
results[["A_size_axis"]] <- tibble(
  fingerprint = "A: size axis (Ers ~ predicted size)",
  metric = c("R2_Ers_on_PBW", "R2_Ers_on_PFVC", "n"),
  value  = c(r2_size_pbw, r2_size_pfvc, nrow(size_data))
)

# =============================================================================
# (B) Error 2 fingerprint: residual SIZE-demographic structure (sex/race/height)
# =============================================================================
# Sex, race, and height move predicted lung SIZE but not intrinsic recoil, so they
# are a near-pure size-surrogate axis. AGE IS EXCLUDED here (it is confounded with
# recoil; handled separately in B2). Lower residual R^2 for Ers x PFVC means PFVC
# removes more of the size structure -- the size-surrogate error is PBW-specific.
size_demo_rhs <- "sex_category + race_category + height10"
demo_data <- mech %>% filter(is.finite(height10))
r2_sd_pbw  <- r2(paste("ers_pbw ~",  size_demo_rhs), demo_data)
r2_sd_pfvc <- r2(paste("ers_pfvc ~", size_demo_rhs), demo_data)
results[["B_residual_size_demographics"]] <- tibble(
  fingerprint = "B: residual sex/race/height structure (error 2)",
  metric = c("R2_sizedemo_ErsxPBW", "R2_sizedemo_ErsxPFVC",
             "R2_reduction_PBW_minus_PFVC", "n"),
  value  = c(r2_sd_pbw, r2_sd_pfvc, r2_sd_pbw - r2_sd_pfvc, nrow(demo_data))
)

# =============================================================================
# (B2) Age trend, read against the recoil prior (NOT minimized)
# =============================================================================
# Age is a common cause of both errors. The recoil prior says correctly-normalized
# specific elastance should DECLINE with age. A flat Ers x PBW means the size and
# recoil declines cancel (artifact); a declining Ers x PFVC recovers the recoil
# loss. Slopes are per +10 yr of age.
age_data <- mech %>% filter(is.finite(age10))
sl_ers   <- slope("ers ~ age10",      age_data, "age10")
sl_pbw   <- slope("ers_pbw ~ age10",  age_data, "age10")
sl_pfvc  <- slope("ers_pfvc ~ age10", age_data, "age10")
results[["B2_age_slopes"]] <- tibble(
  fingerprint = "B2: age slope per +10yr (read vs recoil prior)",
  metric = c("slope_Ers_per10yr", "slope_ErsxPBW_per10yr",
             "slope_ErsxPFVC_per10yr", "n"),
  value  = c(sl_ers, sl_pbw, sl_pfvc, nrow(age_data))
)

# =============================================================================
# (C) Error 1 fingerprint: residual specific elastance vs lung injury
# =============================================================================
# After PFVC normalization, does specific elastance still track injury severity
# (oxygenation, SOFA)? Injury markers are NON-mechanical (plateau/DP are excluded
# as algebraically tied to Ers). A non-zero association is the recoil signal.
injury_rhs  <- "sf10 + sofa_total"
injury_data <- mech %>% filter(is.finite(sf10), is.finite(sofa_total))
r2_injury_pfvc <- r2(paste("ers_pfvc ~", injury_rhs), injury_data)
injury_coefs <- broom::tidy(lm(as.formula(paste("ers_pfvc ~", injury_rhs)),
                               data = injury_data)) %>% filter(term != "(Intercept)")
results[["C_injury_association"]] <- tibble(
  fingerprint = "C: residual specific elastance vs injury (error 1)",
  metric = c("R2_injury_ErsxPFVC", paste0("beta_", injury_coefs$term),
             paste0("p_", injury_coefs$term), "n"),
  value  = c(r2_injury_pfvc, injury_coefs$estimate, injury_coefs$p.value,
             nrow(injury_data))
)

# =============================================================================
# (C2) Error 1 fingerprint: specific elastance vs PREDICTED recoil (FEV1/FVC)
# =============================================================================
# Predicted FEV1/FVC is a demographic estimate of intrinsic recoil (lower = less
# recoil). If specific elastance is recoil-driven, Ers x PFVC should rise with
# predicted FEV1/FVC (lower predicted recoil -> lower specific elastance). Also fit
# Ers on predicted SIZE (PFVC) and predicted RECOIL (FEV1/FVC) jointly to see which
# physiologic projection the demographic structure aligns with -- NOTE both are
# deterministic demographic functions and hence collinear (correlation reported).
recoil_data <- mech %>% filter(is.finite(pred_fev1fvc), is.finite(pfvc))
recoil_fit  <- lm(ers_pfvc ~ pred_fev1fvc, data = recoil_data)
recoil_coef <- broom::tidy(recoil_fit) %>% filter(term == "pred_fev1fvc")
cor_size_recoil <- cor(recoil_data$pfvc, recoil_data$pred_fev1fvc)
joint_fit <- lm(ers ~ scale(pfvc) + scale(pred_fev1fvc), data = recoil_data)
joint_coefs <- broom::tidy(joint_fit) %>% filter(term != "(Intercept)")
results[["C2_predicted_recoil"]] <- tibble(
  fingerprint = "C2: specific elastance vs predicted recoil (FEV1/FVC)",
  metric = c("beta_ErsxPFVC_on_predFEV1FVC", "p_ErsxPFVC_on_predFEV1FVC",
             "cor_PFVC_predFEV1FVC",
             "std_beta_Ers_on_PFVC", "std_beta_Ers_on_predFEV1FVC", "n"),
  value  = c(recoil_coef$estimate, recoil_coef$p.value, cor_size_recoil,
             joint_coefs$estimate[joint_coefs$term == "scale(pfvc)"],
             joint_coefs$estimate[joint_coefs$term == "scale(pred_fev1fvc)"],
             nrow(recoil_data))
)

# =============================================================================
# (D) Variance partition: Ers attributed to a size axis vs an injury axis
# =============================================================================
# Commonality analysis on a single complete-case set. Age is intentionally NOT a
# standalone axis (confounded); the size axis (PFVC / PBW) already carries its
# demographic size effect, and the injury axis (oxygenation, SOFA) is the clean
# error-1 probe. "Unique injury" variance is the cleanest recoil signal.
vp_data <- mech %>% filter(is.finite(pbw), is.finite(pfvc),
                           is.finite(sf10), is.finite(sofa_total))
partition <- function(size_var, data) {
  r2_size   <- r2(paste("ers ~", size_var), data)
  r2_injury <- r2(paste("ers ~", injury_rhs), data)
  r2_full   <- r2(paste("ers ~", size_var, "+", injury_rhs), data)
  tibble(size_surrogate = size_var,
         unique_size   = r2_full - r2_injury,
         unique_injury = r2_full - r2_size,
         shared        = r2_size + r2_injury - r2_full,
         total_r2      = r2_full)
}
vp <- bind_rows(partition("pbw", vp_data), partition("pfvc", vp_data)) %>%
  mutate(size_surrogate = recode(size_surrogate, pbw = "PBW", pfvc = "PFVC"))
results[["D_variance_partition"]] <- vp %>%
  pivot_longer(c(unique_size, unique_injury, shared, total_r2),
               names_to = "metric", values_to = "value") %>%
  transmute(fingerprint = "D: Ers variance partition (size vs injury)",
            metric = paste0(metric, "_", size_surrogate), value) %>%
  bind_rows(tibble(fingerprint = "D: Ers variance partition (size vs injury)",
                   metric = "n", value = nrow(vp_data)))

# =============================================================================
# Write the poolable summary table
# =============================================================================
summary_tbl <- bind_rows(results) %>% mutate(site = site_name, .before = 1)
write_csv(summary_tbl,
          file.path(final_dir, paste0("elastance_fingerprints_summary_", site_name, ".csv")))

message("Fingerprint summary:")
message("  (A) Ers explained by size: R2(PBW)=", round(r2_size_pbw, 3),
        ", R2(PFVC)=", round(r2_size_pfvc, 3))
message("  (B) Residual sex/race/height R2: ErsxPBW=", round(r2_sd_pbw, 3),
        ", ErsxPFVC=", round(r2_sd_pfvc, 3), " (error 2)")
message("  (B2) Age slope per +10yr: Ers=", round(sl_ers, 2),
        ", ErsxPBW=", round(sl_pbw, 1), " (flat = artifact), ErsxPFVC=",
        round(sl_pfvc, 1), " (decline = recoil recovered)")
message("  (C) Residual ErsxPFVC vs injury R2=", round(r2_injury_pfvc, 3))
message("  (C2) ErsxPFVC vs predicted FEV1/FVC: beta=", round(recoil_coef$estimate, 1),
        ", p=", signif(recoil_coef$p.value, 2),
        " (cor[PFVC,predFEV1FVC]=", round(cor_size_recoil, 2), ")")
message("  (D) Unique injury variance of Ers (PFVC size axis)=",
        round(vp$unique_injury[vp$size_surrogate == "PFVC"], 3))

# =============================================================================
# Figures
# =============================================================================
# Figure 1: size-demographic (error 2) on height + age trend (read vs recoil prior)
f1a <- ggplot(demo_data, aes(height_cm, ers_pfvc)) +
  geom_point(color = "#0072B2", alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#E69F00") +
  labs(x = "Height (cm)", y = "Ers x PFVC",
       title = "Ers x PFVC vs height (size axis)") + theme_minimal(base_size = 10)
f1b <- ggplot(demo_data, aes(height_cm, ers_pbw)) +
  geom_point(color = "#D55E00", alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#E69F00") +
  labs(x = "Height (cm)", y = "Ers x PBW",
       title = "Ers x PBW vs height (size axis)") + theme_minimal(base_size = 10)
f1c <- ggplot(age_data, aes(age_at_admission, ers_pfvc)) +
  geom_point(color = "#0072B2", alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#E69F00") +
  labs(x = "Age (years)", y = "Ers x PFVC",
       title = "Ers x PFVC vs age (decline = recoil recovered)") +
  theme_minimal(base_size = 10)
f1d <- ggplot(age_data, aes(age_at_admission, ers_pbw)) +
  geom_point(color = "#D55E00", alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#E69F00") +
  labs(x = "Age (years)", y = "Ers x PBW",
       title = "Ers x PBW vs age (flat = errors cancel)") +
  theme_minimal(base_size = 10)
fig_size <- (f1a | f1b) / (f1c | f1d) +
  plot_annotation(
    title = "Size-surrogate (error 2, top) and the age/recoil trend (bottom)",
    subtitle = paste0(site_name,
      " — on size axes (height) flatter Ers x PFVC = better surrogate; vs age, a ",
      "declining Ers x PFVC recovers recoil loss while a flat Ers x PBW is the artifact"))
ggsave(file.path(final_dir, paste0("elastance_fingerprints_size_", site_name, ".pdf")),
       fig_size, width = 10, height = 8)

# Figure 2: recoil signals (injury + predicted FEV1/FVC) and variance partition
f2a <- ggplot(injury_data, aes(sf_ratio, ers_pfvc)) +
  geom_point(color = "#0072B2", alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#E69F00") +
  labs(x = "SF ratio (lower = worse injury)", y = "Ers x PFVC",
       title = paste0("vs injury (R2 = ", round(r2_injury_pfvc, 3), ")")) +
  theme_minimal(base_size = 10)
f2b <- ggplot(recoil_data, aes(pred_fev1fvc, ers_pfvc)) +
  geom_point(color = "#009E73", alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#E69F00") +
  labs(x = "Predicted FEV1/FVC (lower = less recoil)", y = "Ers x PFVC",
       title = "vs predicted recoil") + theme_minimal(base_size = 10)
vp_long <- vp %>%
  pivot_longer(c(unique_size, unique_injury, shared),
               names_to = "component", values_to = "variance") %>%
  mutate(component = recode(component, unique_size = "Unique size",
                            unique_injury = "Unique injury", shared = "Shared"),
         component = factor(component, levels = c("Unique size", "Shared", "Unique injury")))
f2c <- ggplot(vp_long, aes(size_surrogate, variance, fill = component)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("Unique size" = "#0072B2", "Shared" = "#999999",
                               "Unique injury" = "#D55E00"), name = NULL) +
  labs(x = "Size axis", y = "Share of Ers variance (R2)",
       title = "Variance partition: size vs injury") + theme_minimal(base_size = 10)
fig_recoil <- (f2a | f2b | f2c) +
  plot_annotation(
    title = "Recoil signals (error 1): injury and predicted recoil, plus partition",
    subtitle = paste0(site_name,
      " — Ers x PFVC tracking injury and rising with predicted FEV1/FVC, and ",
      "'unique injury' variance, indicate non-constant recoil independent of size"))
ggsave(file.path(final_dir, paste0("elastance_fingerprints_recoil_", site_name, ".pdf")),
       fig_recoil, width = 13, height = 4.5)

message("Elastance-fingerprints figures written.")
message("Exploratory elastance-fingerprints analysis complete.")
