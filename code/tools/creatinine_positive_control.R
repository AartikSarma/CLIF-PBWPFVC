# =============================================================================
# Creatinine positive control: does creatinine rise where kidney injury is expected?
# =============================================================================
# The joint models found no creatinine divergence by predicted lung size. That is
# only informative if creatinine behaves like a kidney-injury marker in this panel
# at all. Two gaps could flatten it: patients with end-stage renal disease are not
# excluded (their creatinine follows dialysis sessions), and only CRRT censors the
# trajectory, so intermittent haemodialysis pulls creatinine down unseen. This
# script tabulates the observed change in creatinine from the day-0 value over the
# first week, by stratum:
#
#   all patients
#   ESRD (ICD-10 N18.5, N18.6, Z99.2; ICD-9 585.5, 585.6, V45.11), from hospital_diagnosis
#   any haemodialysis procedure in the hospitalization (CPT 90935, 90937, 90945,
#     90947; ICD-10-PCS 5A1D70Z, 5A1D80Z, 5A1D90Z), from patient_procedures
#   day-0 vasopressor dose (NE-equivalents, ug/kg/min: none, under 0.1, 0.1 or more),
#     excluding ESRD
#   vasopressor trend from day 0 to day 2 (rising or falling by 0.05 or more, stable,
#     none), excluding ESRD
#
# Since 2026-09-21 the panel censors ESRD patients at day 0 and treats any dialysis
# procedure as the start of RRT, so on a current panel the ESRD row is empty and the
# dialysis row shows the trajectory before dialysis only.
#
# Expected if creatinine works: it rises more with higher and rising pressor
# requirements. If it falls or stays flat even there, the marker (or the panel) is
# not measuring kidney injury, and the null divergence says nothing.
#
# Reads the 7-day panel (21_biotrauma_panel.R with PBWPFVC_JM_GRID=daily
# PBWPFVC_JM_HORIZON=7) and two CLIF tables. Writes AGGREGATES ONLY, every cell of
# fewer than 10 patients suppressed, to final/injury/:
#   creatinine_check_7d_{site}.csv          change from day 0 by stratum and day
#   creatinine_check_counts_7d_{site}.csv   patients, CRRT starts and deaths by stratum
#
# Usage: PBWPFVC_JM_GRID=daily PBWPFVC_JM_HORIZON=7 Rscript code/tools/creatinine_positive_control.R
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(here) })
rm(list = ls())
source("utils/config.R")
site_name <- config$site_name
source(here("code", "20_biotrauma_grid.R"))   # h_suffix
if (h_suffix != "7d") stop("run with PBWPFVC_JM_GRID=daily PBWPFVC_JM_HORIZON=7 (this is the 7-day panel's check)")
output_dir <- config$output_dir
final_dir  <- final_dir_for("injury")
tables_path <- path.expand(config$tables_path)
MIN_CELL <- 10L

long <- read_parquet(file.path(output_dir, paste0("jm_long_", h_suffix, ".parquet")))
surv <- read_parquet(file.path(output_dir, paste0("jm_surv_", h_suffix, ".parquet")))

# ---- ESRD and haemodialysis flags from the CLIF tables (NA when a table is absent)
read_table <- function(tbl, cols) {
  f <- file.path(tables_path, paste0("clif_", tbl, ".", config$file_type))
  if (!file.exists(f)) { message("*** no ", tbl, " table at ", tables_path, ": that flag is left out"); return(NULL) }
  switch(config$file_type,
         parquet = open_dataset(f) %>% select(all_of(cols)) %>% collect(),
         csv     = read_csv(f, show_col_types = FALSE, col_select = all_of(cols)),
         fst     = fst::read_fst(f, columns = cols))
}
ESRD_CODES <- c("N185", "N186", "Z992", "5855", "5856", "V4511")
HD_CODES   <- c("90935", "90937", "90945", "90947", "5A1D70Z", "5A1D80Z", "5A1D90Z")
dx <- read_table("hospital_diagnosis", c("hospitalization_id", "diagnosis_code"))
px <- read_table("patient_procedures", c("hospitalization_id", "procedure_code"))
esrd_ids <- if (is.null(dx)) NULL else
  dx %>% filter(toupper(gsub(".", "", diagnosis_code, fixed = TRUE)) %in% ESRD_CODES) %>% distinct(hospitalization_id) %>% pull()
hd_ids <- if (is.null(px)) NULL else
  px %>% filter(toupper(procedure_code) %in% HD_CODES) %>% distinct(hospitalization_id) %>% pull()

# ---- patient strata
dose2 <- long %>% filter(period == 2L) %>% select(hospitalization_id, ne_equiv_2 = ne_equiv_peak)
pts <- surv %>%
  filter(!is.na(creatinine_0), creatinine_0 > 0) %>%          # a creatinine trajectory exists (not on CRRT at the index)
  left_join(dose2, by = "hospitalization_id") %>%
  transmute(hospitalization_id, creatinine_0,
            esrd = if (is.null(esrd_ids)) NA else hospitalization_id %in% esrd_ids,
            hd_procedure = if (is.null(hd_ids)) NA else hospitalization_id %in% hd_ids,
            pressor_day0 = case_when(is.na(ne_equiv_0) ~ NA_character_, ne_equiv_0 == 0 ~ "none",
                                     ne_equiv_0 < 0.1 ~ "under 0.1", TRUE ~ "0.1 or more"),
            pressor_trend = case_when(is.na(ne_equiv_2) | is.na(ne_equiv_0) ~ "no day-2 value",
                                      ne_equiv_0 == 0 & ne_equiv_2 == 0 ~ "none",
                                      ne_equiv_2 - ne_equiv_0 >= 0.05 ~ "rising",
                                      ne_equiv_2 - ne_equiv_0 <= -0.05 ~ "falling",
                                      TRUE ~ "stable"),
            crrt_7d  = !is.na(rrt_day) & rrt_day >= 0 & rrt_day <= JM_HORIZON,
            death_7d = !is.na(death_day) & death_day <= JM_HORIZON)
not_esrd <- if (all(is.na(pts$esrd))) rep(TRUE, nrow(pts)) else !pts$esrd
strata <- bind_rows(
  pts %>% mutate(stratifier = "all", stratum = "all"),
  pts %>% filter(!is.na(esrd)) %>% mutate(stratifier = "ESRD", stratum = if_else(esrd, "ESRD", "no ESRD")),
  pts %>% filter(!is.na(hd_procedure)) %>% mutate(stratifier = "haemodialysis procedure",
                                                   stratum = if_else(hd_procedure, "any", "none")),
  pts[not_esrd, ] %>% filter(!is.na(pressor_day0)) %>% mutate(stratifier = "day-0 pressor (no ESRD)", stratum = pressor_day0),
  pts[not_esrd, ] %>% mutate(stratifier = "pressor trend day 0-2 (no ESRD)", stratum = pressor_trend))

# ---- change from day 0 by day (creatinine is already censored at CRRT start in the panel)
daily <- long %>% filter(period >= 1L, period <= JM_HORIZON, !is.na(creatinine), creatinine > 0) %>%
  select(hospitalization_id, day = period, creatinine)
change <- strata %>% select(hospitalization_id, stratifier, stratum, creatinine_0) %>%
  inner_join(daily, by = "hospitalization_id", relationship = "many-to-many") %>%
  mutate(log_change = log(creatinine) - log(creatinine_0),
         kdigo_ratio = creatinine >= 1.5 * creatinine_0, kdigo_rise = creatinine - creatinine_0 >= 0.3) %>%
  group_by(stratifier, stratum, day) %>%
  summarise(n_patients = n(), mean_log_change = mean(log_change), median_log_change = median(log_change),
            pct_1.5x_baseline = 100 * mean(kdigo_ratio), pct_rise_0.3 = 100 * mean(kdigo_rise), .groups = "drop") %>%
  mutate(across(c(mean_log_change, median_log_change, pct_1.5x_baseline, pct_rise_0.3),
                ~ if_else(n_patients < MIN_CELL, NA_real_, round(.x, 3))),
         n_patients = if_else(n_patients < MIN_CELL, NA_integer_, n_patients)) %>%
  # secondary suppression: a count hidden in one stratum could be recovered by subtracting
  # the others from the total, so when any cell of a stratifier (on that day) is under
  # the minimum, that count is hidden for the whole stratifier
  group_by(stratifier, day) %>% mutate(n_patients = if (anyNA(n_patients)) NA_integer_ else n_patients) %>%
  ungroup() %>% mutate(site = site_name)
counts <- strata %>% group_by(stratifier, stratum) %>%
  summarise(n_patients = n(), n_crrt_7d = sum(crrt_7d), n_death_7d = sum(death_7d), .groups = "drop") %>%
  group_by(stratifier) %>%
  mutate(across(c(n_patients, n_crrt_7d, n_death_7d), ~ if (any(.x < MIN_CELL)) NA_integer_ else .x)) %>%
  ungroup() %>% mutate(site = site_name)

write_csv(change, file.path(final_dir, paste0("creatinine_check_", h_suffix, "_", site_name, ".csv")))
write_csv(counts, file.path(final_dir, paste0("creatinine_check_counts_", h_suffix, "_", site_name, ".csv")))
options(width = 200)
message("--- patients, CRRT starts and deaths within 7 days, by stratum (cells under ", MIN_CELL, " suppressed)")
print(as.data.frame(counts %>% select(-site)), row.names = FALSE)
message("--- mean change in log creatinine from day 0, by day (0.1 is about a 10% rise)")
print(as.data.frame(change %>% select(stratifier, stratum, day, mean_log_change) %>%
                      pivot_wider(names_from = day, names_prefix = "day_", values_from = mean_log_change)), row.names = FALSE)
message("--- percent at 1.5 x the day-0 creatinine or more (KDIGO ratio criterion), by day")
print(as.data.frame(change %>% select(stratifier, stratum, day, pct_1.5x_baseline) %>%
                      pivot_wider(names_from = day, names_prefix = "day_", values_from = pct_1.5x_baseline)), row.names = FALSE)
message("creatinine positive control -> ", final_dir)
