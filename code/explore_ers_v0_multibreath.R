# =============================================================================
# Exploratory: multi-breath Ers/V0 (slope of driving pressure vs tidal volume)
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Exploratory, standalone analysis (NOT part of the 00 pipeline runner).
#
# In the Chiumello / Gattinoni / Protti stress-strain framework, varying tidal
# volume at a fixed PEEP traces out the driving-pressure-vs-tidal-volume line
# whose slope is the respiratory-system elastance (Ers = dDP/dVT). This script
# finds subjects who happened to be ventilated at >= 2 distinct tidal volumes
# (spread >= 5%, to avoid unstable slopes from near-identical VTs) -- each with a
# recorded plateau pressure, under a volume-targeted mode (AC-VC or PRVC; pure
# pressure control is excluded) -- within the first 6 hours of ventilation, and
# for each such subject:
#   * collects the (VT, driving pressure) pairs,
#   * fits the per-subject slope of DP vs VT (the "Ers/V0" estimate, = Ers) and
#     checks it against the single-point DP/VT elastance,
#   * plots DP vs VT one panel per subject, and
#   * plots the per-subject slope against age, sex, race, and height.
# The per-subject slope -- not a pooled regression across subjects -- is the unit
# of analysis: subjects differ systematically by predicted lung size (the PBW vs
# PFVC finding), so a single line across subjects would be misleading.
#
# Caveats (exploratory): the slope is a clean elastance estimate only if PEEP is
# unchanged across the pairs (otherwise recruitment/derecruitment confounds it),
# so PEEP constancy is recorded and flagged. Plateau pressure is an observed value
# (never forward-filled), so only timepoints with a genuinely measured plateau
# contribute.
#
# Inputs : output/<site>/intermediate/resp_support_waterfall_clean.parquet (script 02)
#          output/<site>/intermediate/cohort_hospitalization_ids.rds        (script 01)
# Outputs: final/ers_v0_per_subject_<site>.csv   (per-subject; LOCAL ONLY, patient-level)
#          final/ers_v0_summary_<site>.csv        (aggregate distribution; poolable)
#          final/ers_v0_dp_vs_vt_<site>.pdf       (DP-vs-VT panels, one per subject)
#          final/ers_v0_vs_demographics_<site>.pdf (per-subject slope vs age/sex/race/height)
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(lubridate)
library(patchwork)

source("utils/config.R")
site_name <- config$site_name

output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

resp_waterfall <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet"))
cohort_ids     <- readRDS(file.path(output_dir, "cohort_hospitalization_ids.rds"))

# Subject-level demographics for the slope-vs-covariate plots.
demographics <- read_parquet(file.path(output_dir, "cohort_demographics.parquet")) %>%
  select(hospitalization_id, age_at_admission, sex_category, race_category)
heights <- read_parquet(file.path(output_dir, "cohort_heights.parquet")) %>%
  select(hospitalization_id, height_cm)

WINDOW_HOURS  <- 6     # "first six hours of ventilation" (anchored at first IMV timepoint)
MIN_VT_CHANGE <- 0.05  # require >= 5% spread in tidal volume between timepoints; smaller
                       # changes give unstable dDP/dVT slopes (often data-entry artifacts)

# =============================================================================
# Assemble qualifying (VT, driving pressure) measurements in the first 6 hours
# =============================================================================
# Keep IMV timepoints for cohort subjects that have a measured plateau pressure,
# a set PEEP, and a set tidal volume, then restrict to the first 6 h of each
# subject's ventilation and compute driving pressure (plateau - PEEP).
#
# Volume-targeted modes only: a "set tidal volume" is physiologically meaningful
# only when the ventilator targets a volume. Both assist control-volume control
# (AC-VC) and pressure-regulated volume control (PRVC) contain "volume control" in
# the CLIF mode_category, whereas pure "pressure control" does not -- a
# (forward-filled) set tidal volume recorded under pressure control is spurious.
# mode_category is lowercased by the waterfall.
imv_pressures <- resp_waterfall %>%
  filter(hospitalization_id %in% cohort_ids,
         tolower(device_category) == "imv",
         str_detect(coalesce(mode_category, ""), "volume control"),
         !is.na(plateau_pressure_obs), !is.na(peep_set),
         !is.na(tidal_volume_set), tidal_volume_set > 0) %>%
  select(hospitalization_id, recorded_dttm, mode_category,
         tidal_volume_set, peep_set, plateau_pressure_obs) %>%
  group_by(hospitalization_id) %>%
  mutate(imv_start_dttm = min(recorded_dttm)) %>%
  ungroup() %>%
  filter(recorded_dttm <= imv_start_dttm + lubridate::hours(WINDOW_HOURS)) %>%
  mutate(dp = plateau_pressure_obs - peep_set) %>%
  filter(dp > 0)

# Subjects qualify if they have >= 2 DISTINCT tidal volumes (each with a plateau)
# in the window -- i.e. at least two points on the DP-vs-VT line -- AND the spread
# between the smallest and largest VT is >= MIN_VT_CHANGE. The VT-spread filter
# drops near-identical VTs whose tiny denominator (dVT) produces improbably large
# slopes that are usually data-entry errors.
subject_counts <- imv_pressures %>%
  group_by(hospitalization_id) %>%
  summarise(
    n_vt_levels = n_distinct(tidal_volume_set),
    vt_min      = min(tidal_volume_set),
    vt_max      = max(tidal_volume_set),
    .groups = "drop"
  ) %>%
  mutate(pct_vt_change = (vt_max - vt_min) / vt_min)

qualifying_ids <- subject_counts %>%
  filter(n_vt_levels >= 2, pct_vt_change >= MIN_VT_CHANGE) %>%
  pull(hospitalization_id)
multibreath <- imv_pressures %>% filter(hospitalization_id %in% qualifying_ids)

message("Subjects with a measured plateau at >= 2 distinct tidal volumes ",
        "(>= ", round(MIN_VT_CHANGE * 100), "% VT spread) in the first ",
        WINDOW_HOURS, " h: ", length(qualifying_ids),
        " of ", length(cohort_ids), " cohort subjects (",
        sum(subject_counts$n_vt_levels >= 2 & subject_counts$pct_vt_change < MIN_VT_CHANGE),
        " dropped for < ", round(MIN_VT_CHANGE * 100), "% VT change).")

if (length(qualifying_ids) == 0) {
  stop("No subjects had a measured plateau pressure at two or more tidal volumes ",
       ">= ", round(MIN_VT_CHANGE * 100), "% apart within the first ", WINDOW_HOURS,
       " h of ventilation. This exploratory analysis requires within-patient ",
       "tidal-volume variation with concurrent plateau measurements.")
}

# =============================================================================
# Per-subject Ers/V0 = slope of driving pressure vs tidal volume
# =============================================================================
# Slope via the least-squares estimate cov(VT, DP) / var(VT) (equals the lm slope;
# for exactly two points it is the simple rise/run). Units: cmH2O/mL; also
# reported as cmH2O/L (x1000) to match the pipeline's elastance scale.
ers_v0 <- multibreath %>%
  group_by(hospitalization_id) %>%
  summarise(
    n_points       = n(),
    n_vt_levels    = n_distinct(tidal_volume_set),
    vt_min_ml      = min(tidal_volume_set),
    vt_max_ml      = max(tidal_volume_set),
    dp_min_cmh2o   = min(dp),
    dp_max_cmh2o   = max(dp),
    peep_constant  = n_distinct(peep_set) == 1,
    ers_v0_cmh2o_per_ml = cov(tidal_volume_set, dp) / var(tidal_volume_set),
    r_squared      = if (n_distinct(dp) > 1) cor(tidal_volume_set, dp)^2 else NA_real_,
    # Single-point elastance check: Ers/V0 = dDP/dVT should be close to the
    # point-wise DP/VT (mean across the subject's measurements). cmH2O/L.
    ers_single_point_cmh2o_per_l = mean(dp / tidal_volume_set) * 1000,
    .groups = "drop"
  ) %>%
  mutate(ers_v0_cmh2o_per_l = ers_v0_cmh2o_per_ml * 1000) %>%
  left_join(demographics, by = "hospitalization_id") %>%
  left_join(heights, by = "hospitalization_id") %>%
  arrange(desc(peep_constant), hospitalization_id)

# Validation: the multi-point slope should track the single-point DP/VT elastance.
valid <- ers_v0 %>% filter(is.finite(ers_v0_cmh2o_per_l),
                           is.finite(ers_single_point_cmh2o_per_l))
if (nrow(valid) >= 2) {
  message("Validation — per-subject slope (Ers/V0) vs single-point Ers (mean DP/VT): ",
          "Pearson r = ",
          round(cor(valid$ers_v0_cmh2o_per_l, valid$ers_single_point_cmh2o_per_l), 3),
          " across ", nrow(valid), " subjects.")
}

# Per-subject table is patient-level: written for LOCAL exploration only (the
# output/ tree is gitignored). Do not share without aggregation.
write_csv(ers_v0, file.path(final_dir, paste0("ers_v0_per_subject_", site_name, ".csv")))

# Aggregate distribution (poolable). Restricted to PEEP-constant subjects, for
# whom the slope is an unconfounded elastance estimate. Suppressed if n < 10.
ers_v0_clean <- ers_v0 %>% filter(peep_constant, is.finite(ers_v0_cmh2o_per_l))
MIN_CELL <- 10
ers_v0_summary <- tibble(
  site               = site_name,
  n_subjects_total   = nrow(ers_v0),
  n_subjects_peep_constant = nrow(ers_v0_clean),
  median_ers_v0_cmh2o_per_l = if (nrow(ers_v0_clean) >= MIN_CELL) median(ers_v0_clean$ers_v0_cmh2o_per_l) else NA_real_,
  q1_ers_v0_cmh2o_per_l     = if (nrow(ers_v0_clean) >= MIN_CELL) quantile(ers_v0_clean$ers_v0_cmh2o_per_l, 0.25) else NA_real_,
  q3_ers_v0_cmh2o_per_l     = if (nrow(ers_v0_clean) >= MIN_CELL) quantile(ers_v0_clean$ers_v0_cmh2o_per_l, 0.75) else NA_real_,
  suppressed         = nrow(ers_v0_clean) < MIN_CELL
)
write_csv(ers_v0_summary, file.path(final_dir, paste0("ers_v0_summary_", site_name, ".csv")))

message("Median Ers/V0 (PEEP-constant subjects, n=", nrow(ers_v0_clean), "): ",
        if (nrow(ers_v0_clean) >= MIN_CELL)
          paste0(round(median(ers_v0_clean$ers_v0_cmh2o_per_l), 1), " cmH2O/L")
        else "suppressed (n < 10)")

# =============================================================================
# Plot: driving pressure vs tidal volume, one panel per subject
# =============================================================================
# Okabe-Ito colors (per project convention): blue points, orange fitted slope.
# PEEP-varying subjects are marked in the panel strip so the confounded slopes are
# obvious. A short anonymous panel id is used instead of the hospitalization_id.
panel_ids <- ers_v0 %>%
  transmute(hospitalization_id,
            panel = paste0("S", row_number(),
                           if_else(peep_constant, "", " (PEEP varies)"),
                           "  Ers/V0=", round(ers_v0_cmh2o_per_l, 1), " cmH2O/L"))

plot_data <- multibreath %>% left_join(panel_ids, by = "hospitalization_id")

n_subj <- length(qualifying_ids)
n_col  <- min(4, max(1, ceiling(sqrt(n_subj))))

dp_vt_plot <- ggplot(plot_data, aes(x = tidal_volume_set, y = dp)) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x,
              color = "#E69F00", linewidth = 0.7) +
  geom_point(color = "#0072B2", size = 1.8) +
  facet_wrap(~ panel, scales = "free", ncol = n_col) +
  labs(
    title = "Driving pressure vs tidal volume (multi-breath Ers/V0)",
    subtitle = paste0(site_name, " — slope of each line estimates Ers (= Ers/V0); ",
                      "panels with 'PEEP varies' are confounded by recruitment"),
    x = "Set tidal volume (mL)",
    y = "Driving pressure, Pplat - PEEP (cmH2O)"
  ) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(size = 7))

ggsave(file.path(final_dir, paste0("ers_v0_dp_vs_vt_", site_name, ".pdf")),
       dp_vt_plot,
       width  = 2.6 * n_col + 1,
       height = 2.4 * ceiling(n_subj / n_col) + 1,
       limitsize = FALSE)

message("DP-vs-VT panel plot written for ", n_subj, " subjects.")

# =============================================================================
# Per-subject Ers/V0 vs demographics
# =============================================================================
# The unit of analysis is the per-subject slope, NOT pooled DP-VT points: a single
# regression line across subjects would conflate within- and between-subject
# variation, and subjects differ systematically by predicted lung size (the PBW vs
# PFVC finding). Restricted to constant-PEEP subjects, for whom the slope is an
# unconfounded elastance estimate.
demo_slopes <- ers_v0 %>% filter(peep_constant, is.finite(ers_v0_cmh2o_per_l))

if (nrow(demo_slopes) == 0) {
  message("No constant-PEEP subjects with a finite slope; skipping demographic plots.")
} else {
  okabe_ito_cat <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
  y_lab <- expression(E[rs] * "/V0 (slope of DP vs VT, cmH"[2] * "O/L)")

  p_age <- ggplot(demo_slopes, aes(age_at_admission, ers_v0_cmh2o_per_l)) +
    geom_point(color = "#0072B2", alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                color = "#E69F00", fill = "#E69F00", alpha = 0.15) +
    labs(x = "Age (years)", y = y_lab) + theme_minimal(base_size = 11)

  p_height <- ggplot(demo_slopes, aes(height_cm, ers_v0_cmh2o_per_l)) +
    geom_point(color = "#0072B2", alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                color = "#E69F00", fill = "#E69F00", alpha = 0.15) +
    labs(x = "Height (cm)", y = y_lab) + theme_minimal(base_size = 11)

  p_sex <- ggplot(demo_slopes, aes(sex_category, ers_v0_cmh2o_per_l, color = sex_category)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.12, alpha = 0.6, size = 2) +
    scale_color_manual(values = okabe_ito_cat, guide = "none") +
    labs(x = "Sex", y = y_lab) + theme_minimal(base_size = 11)

  p_race <- ggplot(demo_slopes, aes(race_category, ers_v0_cmh2o_per_l, color = race_category)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.12, alpha = 0.6, size = 2) +
    scale_color_manual(values = okabe_ito_cat, guide = "none") +
    labs(x = "Race", y = y_lab) + theme_minimal(base_size = 11)

  demo_fig <- (p_age | p_height) / (p_sex | p_race) +
    plot_annotation(
      title = "Per-subject Ers/V0 (multi-breath DP-vs-VT slope) by demographics",
      subtitle = paste0(site_name, " — one point per subject (constant-PEEP, n = ",
                        nrow(demo_slopes), "); slope = dDP/dVT")
    )

  ggsave(file.path(final_dir, paste0("ers_v0_vs_demographics_", site_name, ".pdf")),
         demo_fig, width = 11, height = 9)
  message("Ers/V0-vs-demographics figure written (n = ", nrow(demo_slopes),
          " constant-PEEP subjects).")
}

message("Exploratory Ers/V0 analysis complete.")
