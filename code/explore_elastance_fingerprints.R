# =============================================================================
# Exploratory: statistical fingerprints of the two errors in normalized elastance
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Exploratory, standalone analysis (NOT part of the 00 pipeline runner). Heavily
# caveated -- see below.
#
# Driving pressure / "elastance normalized to predicted size" (e.g. Ers x PBW) is
# meant to recover the intrinsic specific lung elastance E_spec,L. For the
# respiratory system:
#       E_rs = E_cw + E_L = E_cw + E_spec,L / V0
# where E_spec,L is intrinsic recoil, V0 is the aerated (resting) lung volume, and
# E_cw is chest-wall elastance. "E_rs x size" assumes (1) E_spec,L is constant
# across patients [recoil error] and (2) the size surrogate is proportional to V0
# [size-surrogate error], and ignores E_cw. We cannot measure E_spec,L or V0 here
# (that needs esophageal pressure and a measured aerated volume), so the two errors
# cannot be separated absolutely. They do, however, leave different statistical
# fingerprints, which this script probes:
#
#   Error 2 (size surrogate) is DEMOGRAPHICALLY structured -- the part of normalized
#     elastance that (a) changes when PBW is swapped for PFVC and (b) still tracks
#     age/sex/race/height after normalization.
#   Error 1 (recoil non-constancy) tracks LUNG INJURY, not demographics -- residual
#     specific elastance that varies with oxygenation / global severity after the
#     best available size normalization.
#
# Approach: (A) which predicted-size surrogate explains elastance variation better;
# (B) residual demographic structure of Ers x PBW vs Ers x PFVC; (C) does residual
# specific elastance track injury (SF ratio, SOFA) after PFVC normalization;
# (D) a variance partition of Ers into a size axis vs an injury axis. Injury markers
# are deliberately NON-mechanical (oxygenation, SOFA) -- plateau/driving pressure
# are excluded because they are algebraically tied to Ers.
#
# CAVEATS: surrogate-on-surrogate (predicted, not measured, volumes; no esophageal
# pressure so chest wall is folded in; SF ratio is only a crude aerated-fraction
# proxy). This is a variance/association partition, NOT a per-patient deconfound.
#
# Inputs : output/<site>/intermediate/analysis_cross_sectional.parquet (script 03)
# Outputs: final/elastance_fingerprints_summary_<site>.csv  (poolable statistics)
#          final/elastance_fingerprints_size_<site>.pdf      (size + demographic panels)
#          final/elastance_fingerprints_injury_<site>.pdf    (recoil + variance partition)
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(patchwork)
library(broom)

source("utils/config.R")
site_name <- config$site_name

output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

cross_sectional <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"))

MIN_N <- 10  # do not report any model fit on fewer than this many subjects

# Mechanics subset: elastance (and its normalizations) require a measured plateau,
# so restrict to subjects with a finite Ers. Standardize the per-10-unit covariates.
mech <- cross_sectional %>%
  filter(is.finite(ers), is.finite(ers_pbw), is.finite(ers_pfvc)) %>%
  mutate(
    sex_category  = factor(sex_category,  levels = c("Male", "Female")),
    race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
    age10    = age_at_admission / 10,
    height10 = height_cm / 10,
    sf10     = sf_ratio / 10
  )

message("Subjects with a measured respiratory-system elastance: ", nrow(mech))
if (nrow(mech) < MIN_N) {
  stop("Only ", nrow(mech), " subjects have a measured elastance (need >= ", MIN_N,
       "). Too few to explore the elastance fingerprints.")
}

# Helper: R^2 of an lm fit on the given (complete-case) data.
r2 <- function(formula, data) summary(lm(as.formula(formula), data = data))$r.squared

results <- list()  # collected scalar findings -> tidy summary CSV

# =============================================================================
# (A) Size axis: which predicted-size surrogate explains elastance variation?
# =============================================================================
# E_rs ~ E_spec,L / V0, so elastance should fall as predicted lung size rises. The
# surrogate that explains more of the Ers variance is the better "size axis".
size_data <- mech %>% filter(is.finite(pbw), is.finite(pfvc))
r2_size_pbw  <- r2("ers ~ pbw",  size_data)
r2_size_pfvc <- r2("ers ~ pfvc", size_data)
results[["A_size_axis"]] <- tibble(
  fingerprint = "A: size axis (Ers ~ predicted size)",
  metric = c("R2_Ers_on_PBW", "R2_Ers_on_PFVC", "n"),
  value  = c(r2_size_pbw, r2_size_pfvc, nrow(size_data))
)

# =============================================================================
# (B) Error 2 fingerprint: residual demographic structure after normalization
# =============================================================================
# If the size surrogate were perfect and recoil constant, Ers x size would carry no
# demographic structure. The residual demographic R^2 measures the size-surrogate
# error; a LOWER value for Ers x PFVC means PFVC removes more of the size effect.
demo_rhs <- "age10 + sex_category + race_category + height10"
demo_data <- mech %>% filter(is.finite(height10))
r2_demo_pbw  <- r2(paste("ers_pbw ~",  demo_rhs), demo_data)
r2_demo_pfvc <- r2(paste("ers_pfvc ~", demo_rhs), demo_data)
results[["B_residual_demographics"]] <- tibble(
  fingerprint = "B: residual demographic structure (error 2)",
  metric = c("R2_demo_ErsxPBW", "R2_demo_ErsxPFVC",
             "R2_demo_reduction_PBW_minus_PFVC", "n"),
  value  = c(r2_demo_pbw, r2_demo_pfvc, r2_demo_pbw - r2_demo_pfvc, nrow(demo_data))
)

# =============================================================================
# (C) Error 1 fingerprint: residual specific elastance vs lung injury
# =============================================================================
# After the best available size normalization (PFVC), does the residual specific
# elastance still track injury severity (oxygenation, SOFA)? A non-zero association
# is the recoil-non-constancy signal -- present even with a perfect size surrogate.
injury_rhs  <- "sf10 + sofa_total"
injury_data <- mech %>% filter(is.finite(sf10), is.finite(sofa_total))
r2_injury_pfvc <- r2(paste("ers_pfvc ~", injury_rhs), injury_data)
injury_fit <- lm(as.formula(paste("ers_pfvc ~", injury_rhs)), data = injury_data)
injury_coefs <- broom::tidy(injury_fit) %>% filter(term != "(Intercept)")
results[["C_injury_association"]] <- tibble(
  fingerprint = "C: residual specific elastance vs injury (error 1)",
  metric = c("R2_injury_ErsxPFVC", paste0("beta_", injury_coefs$term),
             paste0("p_", injury_coefs$term), "n"),
  value  = c(r2_injury_pfvc, injury_coefs$estimate, injury_coefs$p.value,
             nrow(injury_data))
)

# =============================================================================
# (D) Variance partition: Ers attributed to a size axis vs an injury axis
# =============================================================================
# Commonality analysis on a single complete-case set (so the nested R^2 are
# comparable): split the variance in Ers into the part uniquely explained by
# predicted size, the part uniquely explained by injury severity, and the shared
# part. Run with PBW and with PFVC as the size axis to show how the partition
# shifts. Optionally adjusted for BMI (a partial chest-wall proxy) below.
vp_data <- mech %>% filter(is.finite(pbw), is.finite(pfvc),
                           is.finite(sf10), is.finite(sofa_total))

partition <- function(size_var, data) {
  r2_size   <- r2(paste("ers ~", size_var), data)
  r2_injury <- r2(paste("ers ~", injury_rhs), data)
  r2_full   <- r2(paste("ers ~", size_var, "+", injury_rhs), data)
  tibble(
    size_surrogate = size_var,
    unique_size    = r2_full - r2_injury,
    unique_injury  = r2_full - r2_size,
    shared         = r2_size + r2_injury - r2_full,
    total_r2       = r2_full
  )
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
message("  (B) Residual demographic R2: ErsxPBW=", round(r2_demo_pbw, 3),
        ", ErsxPFVC=", round(r2_demo_pfvc, 3),
        " (PFVC removes ", round((r2_demo_pbw - r2_demo_pfvc) * 100, 1), "% more)")
message("  (C) Residual ErsxPFVC vs injury R2=", round(r2_injury_pfvc, 3))
message("  (D) Unique injury variance of Ers (PFVC size axis)=",
        round(vp$unique_injury[vp$size_surrogate == "PFVC"], 3))

# =============================================================================
# Figures
# =============================================================================
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

# Figure 1: size relationship + residual demographic structure (error 2)
f1a <- ggplot(size_data, aes(pfvc, ers)) +
  geom_point(color = "#0072B2", alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#E69F00") +
  labs(x = "PFVC (L)", y = "Ers (cmH2O/L)",
       title = paste0("Ers vs PFVC (R2 = ", round(r2_size_pfvc, 3), ")")) +
  theme_minimal(base_size = 10)
f1b <- ggplot(size_data, aes(pbw, ers)) +
  geom_point(color = "#0072B2", alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#E69F00") +
  labs(x = "PBW (kg)", y = "Ers (cmH2O/L)",
       title = paste0("Ers vs PBW (R2 = ", round(r2_size_pbw, 3), ")")) +
  theme_minimal(base_size = 10)
f1c <- ggplot(demo_data, aes(age_at_admission, ers_pfvc)) +
  geom_point(color = "#009E73", alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#E69F00") +
  labs(x = "Age (years)", y = "Ers x PFVC",
       title = "Residual structure: Ers x PFVC vs age") +
  theme_minimal(base_size = 10)
f1d <- ggplot(demo_data, aes(age_at_admission, ers_pbw)) +
  geom_point(color = "#D55E00", alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#E69F00") +
  labs(x = "Age (years)", y = "Ers x PBW",
       title = "Residual structure: Ers x PBW vs age") +
  theme_minimal(base_size = 10)
fig_size <- (f1a | f1b) / (f1c | f1d) +
  plot_annotation(
    title = "Size axis and residual demographic structure (error 2)",
    subtitle = paste0(site_name,
      " — flatter Ers x PFVC vs demographics = PFVC the better size surrogate"))
ggsave(file.path(final_dir, paste0("elastance_fingerprints_size_", site_name, ".pdf")),
       fig_size, width = 10, height = 8)

# Figure 2: recoil-vs-injury + variance partition
f2a <- ggplot(injury_data, aes(sf_ratio, ers_pfvc)) +
  geom_point(color = "#0072B2", alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#E69F00") +
  labs(x = "SF ratio (oxygenation; lower = worse injury)", y = "Ers x PFVC",
       title = paste0("Residual specific elastance vs injury (R2 = ",
                      round(r2_injury_pfvc, 3), ")")) +
  theme_minimal(base_size = 10)
vp_long <- vp %>%
  pivot_longer(c(unique_size, unique_injury, shared),
               names_to = "component", values_to = "variance") %>%
  mutate(component = recode(component,
                            unique_size = "Unique size", unique_injury = "Unique injury",
                            shared = "Shared"),
         component = factor(component, levels = c("Unique size", "Shared", "Unique injury")))
f2b <- ggplot(vp_long, aes(size_surrogate, variance, fill = component)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("Unique size" = "#0072B2", "Shared" = "#999999",
                               "Unique injury" = "#D55E00"), name = NULL) +
  labs(x = "Size axis", y = "Share of Ers variance (R2)",
       title = "Variance partition of Ers: size vs injury") +
  theme_minimal(base_size = 10)
fig_injury <- (f2a | f2b) +
  plot_annotation(
    title = "Recoil-vs-injury signal (error 1) and variance partition",
    subtitle = paste0(site_name,
      " — Ers x PFVC still tracking injury, and 'unique injury' variance, "
      , "indicate non-constant recoil independent of size"))
ggsave(file.path(final_dir, paste0("elastance_fingerprints_injury_", site_name, ".pdf")),
       fig_injury, width = 10, height = 5)

message("Elastance-fingerprints figures written.")
message("Exploratory elastance-fingerprints analysis complete.")
