# =============================================================================
# Exploratory: multi-breath Ers/V0 (slope of driving pressure vs tidal volume)
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Exploratory, standalone analysis (NOT part of the 00 pipeline runner).
#
# In the Chiumello / Gattinoni / Protti stress-strain framework, varying tidal
# volume at a fixed PEEP traces out the driving-pressure-vs-tidal-volume line
# whose slope is the respiratory-system elastance (Ers = dDP/dVT). This script
# finds subjects who happened to be ventilated at >= 2 distinct tidal volumes --
# each with a recorded plateau pressure -- within the first 6 hours of ventilation,
# and for each such subject:
#   * collects the (VT, driving pressure) pairs,
#   * fits the slope of DP vs VT (the "Ers/V0" estimate, = Ers), and
#   * plots DP vs VT with the fitted line, one panel per subject.
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
#          final/ers_v0_dp_vs_vt_<site>.pdf       (DP-vs-VT panels)
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(lubridate)

source("utils/config.R")
site_name <- config$site_name

output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

resp_waterfall <- read_parquet(file.path(output_dir, "resp_support_waterfall_clean.parquet"))
cohort_ids     <- readRDS(file.path(output_dir, "cohort_hospitalization_ids.rds"))

WINDOW_HOURS <- 6  # "first six hours of ventilation" (anchored at first IMV timepoint)

# =============================================================================
# Assemble qualifying (VT, driving pressure) measurements in the first 6 hours
# =============================================================================
# Keep IMV timepoints for cohort subjects that have a measured plateau pressure,
# a set PEEP, and a set tidal volume, then restrict to the first 6 h of each
# subject's ventilation and compute driving pressure (plateau - PEEP).
imv_pressures <- resp_waterfall %>%
  filter(hospitalization_id %in% cohort_ids,
         tolower(device_category) == "imv",
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
# in the window -- i.e. at least two points on the DP-vs-VT line.
subject_counts <- imv_pressures %>%
  group_by(hospitalization_id) %>%
  summarise(n_vt_levels = n_distinct(tidal_volume_set), .groups = "drop")

qualifying_ids <- subject_counts %>% filter(n_vt_levels >= 2) %>% pull(hospitalization_id)
multibreath <- imv_pressures %>% filter(hospitalization_id %in% qualifying_ids)

message("Subjects with a measured plateau at >= 2 distinct tidal volumes in the ",
        "first ", WINDOW_HOURS, " h: ", length(qualifying_ids),
        " of ", length(cohort_ids), " cohort subjects.")

if (length(qualifying_ids) == 0) {
  stop("No subjects had a measured plateau pressure at two or more distinct ",
       "tidal volumes within the first ", WINDOW_HOURS, " h of ventilation. ",
       "This exploratory analysis requires within-patient tidal-volume variation ",
       "with concurrent plateau measurements.")
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
    .groups = "drop"
  ) %>%
  mutate(ers_v0_cmh2o_per_l = ers_v0_cmh2o_per_ml * 1000) %>%
  arrange(desc(peep_constant), hospitalization_id)

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
message("Exploratory Ers/V0 analysis complete.")
