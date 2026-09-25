# =============================================================================
# Creatinine positive control: does creatinine rise where kidney injury is expected?
# =============================================================================
# The panel censors ESRD at day 0 and ends creatinine at any RRT, including
# intermittent dialysis; this check confirms creatinine rises where kidney injury is
# expected. A creatinine result from the joint models is informative only if
# creatinine behaves like a kidney-injury marker in the panel. The script tabulates
# the observed change in creatinine from the day-0 value over the first week, by
# stratum:
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
# Because the panel censors ESRD at day 0 and starts RRT at any dialysis procedure,
# the ESRD row is empty and the dialysis row shows the trajectory before dialysis
# only.
#
# Expected if creatinine works: it rises more with higher and rising pressor
# requirements. If it falls or stays flat even there, the marker (or the panel) is
# not measuring kidney injury, and a creatinine result from the joint models says
# nothing.
#
# Reads the 7-day panel (21_biotrauma_panel.R with PBWPFVC_JM_GRID=daily
# PBWPFVC_JM_HORIZON=7) and two CLIF tables. Writes aggregates only, unmasked,
# to final/injury/:
#   creatinine_check_7d_{site}.csv          change from day 0 by stratum and day
#   creatinine_check_counts_7d_{site}.csv   patients, RRT starts (any modality) and deaths by stratum
#
# Usage: PBWPFVC_JM_GRID=daily PBWPFVC_JM_HORIZON=7 uvr run code/tools/creatinine_positive_control.R
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
# Stratum cutpoints for the vasopressor dose, NE-equivalents in ug/kg/min
PRESSOR_DAY0_HIGH  <- 0.1    # day-0 dose: "under 0.1" against "0.1 or more"
PRESSOR_TREND_STEP <- 0.05   # day 0 to day 2: a change this large is "rising" or "falling"
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
  filter(!is.na(creatinine_0), creatinine_0 > 0) %>%          # a creatinine trajectory exists (not on RRT at the index)
  left_join(dose2, by = "hospitalization_id") %>%
  transmute(hospitalization_id, creatinine_0,
            esrd = if (is.null(esrd_ids)) NA else hospitalization_id %in% esrd_ids,
            hd_procedure = if (is.null(hd_ids)) NA else hospitalization_id %in% hd_ids,
            pressor_day0 = case_when(is.na(ne_equiv_0) ~ NA_character_, ne_equiv_0 == 0 ~ "none",
                                     ne_equiv_0 < PRESSOR_DAY0_HIGH ~ paste("under", PRESSOR_DAY0_HIGH),
                                     TRUE ~ paste(PRESSOR_DAY0_HIGH, "or more")),
            pressor_trend = case_when(is.na(ne_equiv_2) | is.na(ne_equiv_0) ~ "no day-2 value",
                                      ne_equiv_0 == 0 & ne_equiv_2 == 0 ~ "none",
                                      ne_equiv_2 - ne_equiv_0 >= PRESSOR_TREND_STEP ~ "rising",
                                      ne_equiv_2 - ne_equiv_0 <= -PRESSOR_TREND_STEP ~ "falling",
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

# ---- change from day 0 by day (creatinine is already censored at RRT start in the panel)
daily <- long %>% filter(period >= 1L, period <= JM_HORIZON, !is.na(creatinine), creatinine > 0) %>%
  select(hospitalization_id, day = period, creatinine)
change <- strata %>% select(hospitalization_id, stratifier, stratum, creatinine_0) %>%
  inner_join(daily, by = "hospitalization_id", relationship = "many-to-many") %>%
  mutate(log_change = log(creatinine) - log(creatinine_0),
         kdigo_ratio = creatinine >= 1.5 * creatinine_0, kdigo_rise = creatinine - creatinine_0 >= 0.3) %>%
  group_by(stratifier, stratum, day) %>%
  summarise(n_patients = n(), mean_log_change = mean(log_change), median_log_change = median(log_change),
            pct_1.5x_baseline = 100 * mean(kdigo_ratio), pct_rise_0.3 = 100 * mean(kdigo_rise), .groups = "drop") %>%
  mutate(across(c(mean_log_change, median_log_change, pct_1.5x_baseline, pct_rise_0.3), ~ round(.x, 3)),
         site = site_name)
counts <- strata %>% group_by(stratifier, stratum) %>%
  summarise(n_patients = n(), n_rrt_7d = sum(crrt_7d), n_death_7d = sum(death_7d), .groups = "drop") %>%
  mutate(site = site_name)

write_csv(change, file.path(final_dir, paste0("creatinine_check_", h_suffix, "_", site_name, ".csv")))
write_csv(counts, file.path(final_dir, paste0("creatinine_check_counts_", h_suffix, "_", site_name, ".csv")))
options(width = 200)
message("--- patients, RRT starts (any modality) and deaths within 7 days, by stratum")
print(as.data.frame(counts %>% select(-site)), row.names = FALSE)
message("--- mean change in log creatinine from day 0, by day (0.1 is about a 10% rise)")
print(as.data.frame(change %>% select(stratifier, stratum, day, mean_log_change) %>%
                      pivot_wider(names_from = day, names_prefix = "day_", values_from = mean_log_change)), row.names = FALSE)
message("--- percent at 1.5 x the day-0 creatinine or more (KDIGO ratio criterion), by day")
print(as.data.frame(change %>% select(stratifier, stratum, day, pct_1.5x_baseline) %>%
                      pivot_wider(names_from = day, names_prefix = "day_", values_from = pct_1.5x_baseline)), row.names = FALSE)
message("creatinine positive control -> ", final_dir)
