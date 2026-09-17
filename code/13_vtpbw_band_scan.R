# =============================================================================
# Script 13 (band scan): who enters the cohort if the VT/PBW gate is widened?
# =============================================================================
# The analytic cohort takes patients with a complete, hypoxemic ventilator
# timepoint at VT/PBW 6-8 mL/kg (the lung-protective band). Many patients sit
# just below 6, and the small-lunged ones are dosed there on purpose, so the
# 6-8 gate thins exactly the discordant tertile the strain question needs.
# This script re-applies the gate over a range of bands to the pre-gate
# per-timepoint table script 03 saves (analysis_all_eligible_timepoints) and
# reports, per lower bound (the 8 mL/kg ceiling stays): patients qualifying, patients added over 6-8, deaths,
# the index VT/PBW and VT/PFVC (percent of predicted FVC), the share above
# the 11% ceiling, and the PBW/PFVC tertile mix. Same index rule as script 03
# (first qualifying timepoint; SF < 315), only the band changes.
#
# Outputs (aggregates only):
#   final/vtpbw_band_scan_{site}.csv   one row per band
#   final/vtpbw_band_scan_{site}.pdf   patients, discordant share and strain by band,
#                                      and the distribution of each patient's
#                                      lowest complete-data VT/PBW around 6
# Usage: Rscript code/13_vtpbw_band_scan.R   (PBWPFVC_BAND_LOWER="4,4.5,5,5.5,6";
#        PBWPFVC_BAND_UPPER="8")
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(here); library(patchwork) })
rm(list = ls())
source("utils/config.R")
site_name  <- config$site_name
output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
LOWERS <- as.numeric(strsplit(Sys.getenv("PBWPFVC_BAND_LOWER", "4,4.5,5,5.5,6"), ",")[[1]])
UPPERS <- as.numeric(strsplit(Sys.getenv("PBWPFVC_BAND_UPPER", "8"), ",")[[1]])   # 8: the lung-protective ceiling stays; above it is another population
SF_HYPOXEMIA_THRESHOLD <- 315
STRAIN_CEILING <- 11   # VT/PFVC, % of predicted FVC
okabe <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7")

ae <- read_parquet(file.path(output_dir, "analysis_all_eligible_timepoints.parquet")) %>%
  filter(has_all_data) %>%
  mutate(vtpfvc_pct = vtpfvc)   # script 03's vtpfvc is VT (mL) / PFVC (L) x 0.1 = percent of predicted FVC
message("Complete-data ventilator timepoints: ", nrow(ae), " on ", n_distinct(ae$hospitalization_id), " patients")

# the reference tertiles of PBW/PFVC, cut on every complete-data patient so the
# band comparison uses one fixed definition of "discordant"
pt <- ae %>% group_by(hospitalization_id) %>%
  summarise(pbwpfvc = first(pbwpfvc), deceased = first(deceased),
            min_vtpbw = min(vtpbw), median_vtpbw = median(vtpbw), .groups = "drop")
disc_breaks <- quantile(pt$pbwpfvc, c(1/3, 2/3), na.rm = TRUE)
disc_of <- function(x) cut(x, c(-Inf, disc_breaks, Inf), labels = c("Concordant", "Mid", "Discordant"))

# the index under a band: the first qualifying timepoint (script 03's tier 2)
index_under <- function(lo, hi) {
  ae %>% filter(vtpbw >= lo, vtpbw <= hi, sf_ratio < SF_HYPOXEMIA_THRESHOLD) %>%
    group_by(hospitalization_id) %>% slice_min(recorded_dttm, n = 1, with_ties = FALSE) %>% ungroup() %>%
    mutate(disc_grp = disc_of(pbwpfvc), band = sprintf("%g-%g", lo, hi), lo = lo, hi = hi)
}
bands <- expand_grid(lo = LOWERS, hi = UPPERS) %>% filter(lo < hi)
idx <- pmap_dfr(bands, index_under)
ref_ids <- idx %>% filter(lo == 6, hi == 8) %>% pull(hospitalization_id)
scan <- idx %>% group_by(band, lo, hi) %>%
  summarise(n_patients = n(),
            n_added_over_6_8 = sum(!hospitalization_id %in% ref_ids),
            n_deaths = sum(deceased == 1, na.rm = TRUE),
            index_vtpbw_median = median(vtpbw), index_vtpbw_q10 = quantile(vtpbw, 0.1),
            index_vtpfvc_pct_median = median(vtpfvc_pct, na.rm = TRUE),
            frac_above_ceiling = mean(vtpfvc_pct > STRAIN_CEILING, na.rm = TRUE),
            frac_discordant = mean(disc_grp == "Discordant"),
            n_discordant = sum(disc_grp == "Discordant"),
            frac_dp_observed = mean(!is.na(dp)),
            .groups = "drop") %>%
  arrange(hi, lo) %>% mutate(site = site_name)
write_csv(scan, file.path(final_dir, paste0("vtpbw_band_scan_", site_name, ".csv")))
message("\nCohort under each VT/PBW band (index = first qualifying timepoint, SF < 315):")
print(as.data.frame(scan %>% select(band, n_patients, n_added_over_6_8, n_deaths, index_vtpbw_median, index_vtpfvc_pct_median,
                                    frac_above_ceiling, frac_discordant, n_discordant) %>%
                      mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)

# who sits just below 6: each patient's LOWEST complete-data VT/PBW
below <- pt %>% mutate(bin = cut(min_vtpbw, c(-Inf, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, Inf))) %>%
  count(bin, disc = disc_of(pbwpfvc))
message("\nPatients by their lowest complete-data VT/PBW (all complete-data patients, not only the hypoxemic):")
print(as.data.frame(pt %>% mutate(bin = cut(min_vtpbw, c(-Inf, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, Inf))) %>% count(bin)), row.names = FALSE)

# figure
sc <- scan %>% mutate(upper = factor(paste("upper", hi)), lower = lo)
p1 <- ggplot(sc, aes(lower, n_patients, colour = upper)) + geom_line() + geom_point() +
  scale_colour_manual(values = okabe[1:2], name = NULL) +
  labs(title = "Patients qualifying", x = "lower VT/PBW bound (mL/kg)", y = "patients")
p2 <- ggplot(sc, aes(lower, frac_discordant, colour = upper)) + geom_line() + geom_point() +
  scale_colour_manual(values = okabe[1:2], name = NULL) +
  labs(title = "Share in the discordant (small-lung) tertile", x = "lower VT/PBW bound (mL/kg)", y = "fraction")
p3 <- ggplot(sc, aes(lower, frac_above_ceiling, colour = upper)) + geom_line() + geom_point() +
  scale_colour_manual(values = okabe[1:2], name = NULL) +
  labs(title = "Share above the 11% VT/PFVC ceiling at the index", x = "lower VT/PBW bound (mL/kg)", y = "fraction")
p4 <- ggplot(pt %>% filter(is.finite(min_vtpbw), min_vtpbw >= 2, min_vtpbw <= 12), aes(min_vtpbw, fill = disc_of(pbwpfvc))) +
  geom_histogram(binwidth = 0.25, boundary = 0, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = c(6, 8), linetype = 2) +
  scale_fill_manual(values = okabe[c(3, 1, 4)], name = "PBW/PFVC tertile", na.translate = FALSE) +
  labs(title = "Each patient's lowest complete-data VT/PBW", x = "VT/PBW (mL/kg)", y = "patients")
p <- (p1 + p2) / (p3 + p4) + plot_layout(guides = "collect") +
  plot_annotation(title = sprintf("Widening the VT/PBW inclusion band (%s)", site_name),
                  subtitle = "index = first complete, hypoxemic timepoint in the band; tertiles cut on all complete-data patients") &
  theme_minimal(base_size = 10)
ggsave(file.path(final_dir, paste0("vtpbw_band_scan_", site_name, ".pdf")), p, width = 11, height = 8)
message("13_vtpbw_band_scan complete -> ", final_dir)
