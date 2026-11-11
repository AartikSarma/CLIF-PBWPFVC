#!/usr/bin/env Rscript
# =============================================================================
# enrollment_feasibility.R -- grant preliminary data: hypoxemic ventilated adults per month
#                             at one site over the last 24 months, split by whether
#                             6 mL/kg PBW exceeds 11% of PFVC
# =============================================================================
# Population: patients in script 03's broad PFVC cohort (analysis_broad_pfvc), that is,
# adults with an ICU stay and a set tidal volume > 0 (script 01) whose height, age, sex
# and race give a PBW and a PFVC (height 150-210 cm), who were hypoxemic on the
# ventilator. No VT/PBW gate is applied: a trial that changes the tidal volume would
# approach patients whatever volume they were given. Script 01 keeps one
# hospitalization per patient (the last that qualifies), so the counts are patients,
# not encounters, and slightly undercount a site's ventilated admissions.
#
# Clock: the month of the first invasive-ventilation row (device imv or a set tidal
# volume > 0, script 01's rule) in script 01's waterfall, since a ventilation trial
# enrolls at intubation. The figure shows the last LOOKBACK_MONTHS complete months; the
# extract's final month is partial and is dropped.
#
# Hypoxemia: any invasive-ventilation timepoint with SF < 315 (Rice 2007's equivalent of
# PF 300; script 03's rule without its VT/PBW band). The SF values are script 03's
# pre-gate IMV rows: SpO2 80-97 paired with FiO2 recorded in the 4 h before it, carried
# for up to 4 h. A patient with no SF on the ventilator counts as not hypoxemic, so the
# counts are a lower bound.
#
# PBW/PFVC split. 6 mL/kg PBW delivers more than 11% of PFVC (about the 75th percentile
# of VT/PFVC in the ARMA low tidal volume arm; the constants of script 03) exactly when
# PBW/PFVC > 11 / 0.6 = 18.3 kg/L. The discordance threshold above which a PFVC-limited
# ceiling showed growing benefit in MIMIC and UCSF lies at about 17-18.
#
# Heights below 150 cm fall outside script 03's height range and so outside this count,
# though they would be among the smallest lungs. The script logs how many such patients
# the cohort holds so the grant can say so.
#
# The figure has no title or subtitle (the grant supplies a caption); the numbers a
# caption needs are in the summary CSV and the log.
#
# Writes to final/supplement/ (aggregates only; no small-cell masking, so check for
# cells < 10 before sharing):
#   enrollment_feasibility_monthly_{site}.csv   month x PBW/PFVC group counts
#   enrollment_feasibility_summary_{site}.csv   totals and monthly means
#   enrollment_feasibility_{site}.png / .pdf    the figure
#
# PBWPFVC_HYPOX_WINDOW_H=24 (say): count hypoxemia only within that many hours of first
# intubation; unset, the whole ventilation course counts.
#
# Usage: uvr run code/tools/enrollment_feasibility.R
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(arrow); library(lubridate); library(ragg)
})
source("utils/config.R")

site_name    <- config$site_name
intermediate <- config$output_dir
out_dir      <- final_dir_for("supplement")

VT_PER_KG_ARMA   <- 6    # mL/kg PBW, the ARMA low-VT arm's prescribed dose (as in script 03)
VT_PCT_PFVC_ARMA <- 11   # % of PFVC, about the 75th percentile delivered in that arm (script 03)
HIGH_RATIO_CUT   <- VT_PCT_PFVC_ARMA / (VT_PER_KG_ARMA / 10)   # 18.33 kg/L
MIN_HEIGHT_CM    <- 150  # script 03's lower height limit
LOOKBACK_MONTHS  <- 24
SF_HYPOXEMIA_THRESHOLD <- 315

hypoxemia_window_h <- as.numeric(Sys.getenv("PBWPFVC_HYPOX_WINDOW_H", "Inf"))
if (is.na(hypoxemia_window_h) || hypoxemia_window_h <= 0)
  stop("PBWPFVC_HYPOX_WINDOW_H must be a positive number of hours")

group_labels <- c(
  high = sprintf("6 mL/kg PBW > %d%% of PFVC (PBW/PFVC ≥ %.1f kg/L)", VT_PCT_PFVC_ARMA, HIGH_RATIO_CUT),
  not  = sprintf("6 mL/kg PBW ≤ %d%% of PFVC", VT_PCT_PFVC_ARMA)
)
# Okabe-Ito vermillion for the target group, grey for the rest; the target group is the
# last factor level, so it sits on the baseline of each stacked bar
group_colours <- setNames(c("#D55E00", "#999999"), unname(group_labels))

# ---- patients, their PBW/PFVC and their intubation time ----------------------------
broad_pfvc <- read_parquet(file.path(intermediate, "analysis_broad_pfvc.parquet"),
                           col_select = c("hospitalization_id", "pbw", "pfvc"))

# Script 01's invasive-ventilation signal (device imv, or any set tidal volume > 0), read
# from its waterfall before script 02's outlier cleaning, which can blank a patient's only
# set tidal volumes.
first_intubation <- read_parquet(file.path(intermediate, "resp_support_waterfall.parquet"),
                                 col_select = c("hospitalization_id", "recorded_dttm",
                                                "device_category", "tidal_volume_set")) %>%
  filter(replace_na(tolower(device_category) == "imv", FALSE) |
           replace_na(as.numeric(tidal_volume_set) > 0, FALSE)) %>%
  summarise(first_imv_dttm = min(recorded_dttm), .by = hospitalization_id)

patients <- broad_pfvc %>%
  left_join(first_intubation, by = "hospitalization_id") %>%
  mutate(pbwpfvc = pbw / pfvc)

# Script 01 admits a patient on any raw respiratory_support row with a set tidal volume,
# but the waterfall it then builds blanks set tidal volumes on pressure-support and
# tracheostomy rows, keeps one row per timestamp, and leaves a device other than imv in
# place. A patient whose only tidal volumes sat on such rows has no ventilation row
# left and so no intubation time. Such patients were never on volume-set invasive
# ventilation in the waterfall's terms; they are dropped and counted.
no_intubation_time <- sum(is.na(patients$first_imv_dttm))
message("Dropping ", no_intubation_time, " of ", nrow(patients), " patients with no invasive-ventilation ",
        "row left in script 01's waterfall (their set tidal volumes were on rows the waterfall cleans)")
patients <- patients %>% filter(!is.na(first_imv_dttm))

# ---- hypoxemia on the ventilator --------------------------------------------------
lowest_sf <- read_parquet(file.path(intermediate, "analysis_all_eligible_timepoints.parquet"),
                          col_select = c("hospitalization_id", "recorded_dttm", "sf_ratio")) %>%
  filter(!is.na(sf_ratio)) %>%
  inner_join(patients %>% select(hospitalization_id, first_imv_dttm), by = "hospitalization_id") %>%
  filter(difftime(recorded_dttm, first_imv_dttm, units = "hours") <= hypoxemia_window_h) %>%
  summarise(min_sf = min(sf_ratio), .by = hospitalization_id)

patients <- patients %>% left_join(lowest_sf, by = "hospitalization_id")
message(sum(is.na(patients$min_sf)), " of ", nrow(patients), " patients have no SF on the ventilator",
        if (is.finite(hypoxemia_window_h)) paste0(" within ", hypoxemia_window_h, " h of intubation") else "",
        " and count as not hypoxemic")

# ---- the last LOOKBACK_MONTHS complete months -----------------------------------
# The extract almost never ends on a month boundary, so its last month is partial and
# would read as a fall in admissions; the window ends the month before it.
partial_month <- floor_date(as_date(max(patients$first_imv_dttm)), unit = "month")
window_months <- seq(partial_month %m-% months(LOOKBACK_MONTHS), by = "month", length.out = LOOKBACK_MONTHS)
if (min(window_months) < floor_date(as_date(min(patients$first_imv_dttm)), unit = "month"))
  stop("The extract starts after ", min(window_months), "; it does not cover ", LOOKBACK_MONTHS, " complete months")
message("Window: ", min(window_months), " to ", max(window_months), " (dropping the partial month ", partial_month, ")")

hypoxemic_patients <- patients %>%
  mutate(month = floor_date(as_date(first_imv_dttm), unit = "month")) %>%
  filter(month %in% window_months, !is.na(min_sf), min_sf < SF_HYPOXEMIA_THRESHOLD) %>%
  mutate(group = if_else(pbwpfvc >= HIGH_RATIO_CUT, group_labels[["high"]], group_labels[["not"]]),
         group = factor(group, levels = rev(group_labels)))

# ---- patients the height range leaves out ------------------------------------------
demographics <- read_parquet(file.path(intermediate, "cohort_demographics.parquet"),
                             col_select = "hospitalization_id")
heights <- read_parquet(file.path(intermediate, "cohort_heights_clean.parquet")) %>%
  semi_join(demographics, by = "hospitalization_id")
below_height_range <- sum(heights$height_cm < MIN_HEIGHT_CM, na.rm = TRUE)
message(below_height_range, " cohort patients (all years) have a height below ", MIN_HEIGHT_CM,
        " cm and are not counted (outside script 03's height range)")

# ---- monthly counts ---------------------------------------------------------------
# Every month in the window appears, with zeros, so an empty month is visible.
monthly_counts <- hypoxemic_patients %>%
  count(month, group, name = "n_patients") %>%
  complete(month = window_months, group, fill = list(n_patients = 0L)) %>%
  arrange(month, group)

month_totals <- monthly_counts %>%
  summarise(n_total = sum(n_patients),
            n_high  = sum(n_patients[group == group_labels[["high"]]]), .by = month)

feasibility_summary <- tibble(
  site                           = site_name,
  first_month                    = min(window_months),
  last_month                     = max(window_months),
  hypoxemia_window_h             = hypoxemia_window_h,
  n_patients                     = sum(month_totals$n_total),
  n_high_ratio                   = sum(month_totals$n_high),
  pct_high_ratio                 = 100 * n_high_ratio / n_patients,
  high_ratio_cut_kg_per_l        = HIGH_RATIO_CUT,
  mean_total_per_month           = mean(month_totals$n_total),
  mean_high_per_month            = mean(month_totals$n_high),
  min_high_per_month             = min(month_totals$n_high),
  max_high_per_month             = max(month_totals$n_high),
  high_per_year                  = mean(month_totals$n_high) * 12,
  months_high_below_10           = sum(month_totals$n_high < 10),
  n_below_height_range_all_years = below_height_range
)

write_csv(monthly_counts %>% mutate(site = site_name, month = as.character(month)),
          file.path(out_dir, paste0("enrollment_feasibility_monthly_", site_name, ".csv")))
write_csv(feasibility_summary,
          file.path(out_dir, paste0("enrollment_feasibility_summary_", site_name, ".csv")))

message(sprintf("%d hypoxemic patients, %s to %s: %.1f a month (range %d-%d, about %.0f a year) with 6 mL/kg PBW > %d%% of PFVC (%.0f%%)",
                feasibility_summary$n_patients, format(min(window_months), "%b %Y"), format(max(window_months), "%b %Y"),
                feasibility_summary$mean_high_per_month, feasibility_summary$min_high_per_month,
                feasibility_summary$max_high_per_month, feasibility_summary$high_per_year,
                VT_PCT_PFVC_ARMA, feasibility_summary$pct_high_ratio))
if (feasibility_summary$months_high_below_10 > 0)
  message("NOTE: ", feasibility_summary$months_high_below_10, " of ", LOOKBACK_MONTHS,
          " months have fewer than 10 patients with 6 mL/kg > 11% PFVC; check before sharing")

# ---- figure -----------------------------------------------------------------------
feasibility_figure <- ggplot(monthly_counts, aes(month, n_patients, fill = group)) +
  geom_col(width = 25) +
  scale_fill_manual(values = group_colours, breaks = unname(group_labels), name = NULL) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y", expand = expansion(mult = 0.01)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Month", y = "Hypoxemic ventilated adults") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", legend.justification = "left",
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())

figure_stem <- file.path(out_dir, paste0("enrollment_feasibility_", site_name))
ggsave(paste0(figure_stem, ".png"), feasibility_figure, width = 8, height = 4.5, dpi = 300, device = agg_png)
ggsave(paste0(figure_stem, ".pdf"), feasibility_figure, width = 8, height = 4.5, device = cairo_pdf)
message("Wrote ", figure_stem, ".png/.pdf and the monthly and summary CSVs to ", out_dir)
