# =============================================================================
# Supplement (diagnostic): driving pressure against VT/PFVC, by VT/PBW bin
# =============================================================================
# A patient ventilated at 4 mL/kg PBW is not a patient at 7 mL/kg with a smaller
# number: the low dose is usually a response to something (a stiff lung, a high
# plateau, ECMO, a small patient the clinician has already down-sized). Before
# widening the VT/PBW inclusion band this asks what that something looks like:
# at the same VT/PFVC (percent of predicted FVC), do the low-dose bins carry a
# higher driving pressure (stiffer lungs) or the same (smaller lungs, dosed to
# size)? Timepoints come from the pre-gate table script 03 saves (all complete-
# data ventilator timepoints), restricted to those with a recorded plateau, and
# are reduced to one row per patient x VT/PBW bin (medians) so a long stay does
# not dominate. Everything drawn is a binned summary, not patient points.
#
# Outputs (aggregates only):
#   final/dp_by_vtpbw_{site}.csv  per VT/PBW bin: patients, timepoints, quantiles
#                                 of driving pressure, compliance, VT/PFVC, PEEP,
#                                 plateau; and the binned means of DP by VT/PFVC
#   final/dp_by_vtpbw_{site}.pdf  DP vs VT/PFVC by VT/PBW bin (binned means with
#                                 95% intervals and a linear fit), and DP,
#                                 compliance, VT/PFVC and PEEP by bin
# Usage: Rscript code/supplement/injury_dp_by_vtpbw.R
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(here); library(patchwork) })
rm(list = ls())
source("utils/config.R")
site_name  <- config$site_name
output_dir <- config$output_dir
final_dir  <- config$final_dir
okabe <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9")
BIN_BREAKS <- c(-Inf, 5, 6, 7, 8, Inf)
BIN_LABELS <- c("< 5", "5-6", "6-7", "7-8", "> 8")

ae <- read_parquet(file.path(output_dir, "analysis_all_eligible_timepoints.parquet")) %>%
  filter(has_all_data, !is.na(dp), dp > 0, !is.na(vtpbw), !is.na(vtpfvc)) %>%
  mutate(vtpbw_bin = cut(vtpbw, BIN_BREAKS, labels = BIN_LABELS),
         crs = tidal_volume_set / dp,                 # mL / cmH2O
         vtpfvc_pct = vtpfvc)                          # percent of predicted FVC
message("Timepoints with a recorded plateau: ", nrow(ae), " on ", n_distinct(ae$hospitalization_id), " patients")

# one row per patient x bin: the patient's median values while dosed in that bin
pb <- ae %>% group_by(hospitalization_id, vtpbw_bin) %>%
  summarise(n_tp = n(), dp = median(dp), crs = median(crs), vtpfvc_pct = median(vtpfvc_pct),
            peep = median(peep_set, na.rm = TRUE), pplat = median(plateau_pressure_obs, na.rm = TRUE),
            pbwpfvc = first(pbwpfvc), .groups = "drop")

q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
by_bin <- pb %>% group_by(vtpbw_bin) %>%
  summarise(n_patients = n(), n_timepoints = sum(n_tp),
            dp_q25 = q(dp, .25), dp_median = median(dp), dp_q75 = q(dp, .75),
            crs_q25 = q(crs, .25), crs_median = median(crs), crs_q75 = q(crs, .75),
            vtpfvc_q25 = q(vtpfvc_pct, .25), vtpfvc_median = median(vtpfvc_pct), vtpfvc_q75 = q(vtpfvc_pct, .75),
            peep_median = median(peep, na.rm = TRUE), pplat_median = median(pplat, na.rm = TRUE),
            pbwpfvc_median = median(pbwpfvc, na.rm = TRUE),
            frac_dp_over_15 = mean(dp > 15), .groups = "drop") %>%
  mutate(site = site_name)
message("\nBy VT/PBW bin (one row per patient x bin; medians):")
print(as.data.frame(by_bin %>% select(vtpbw_bin, n_patients, dp_median, crs_median, vtpfvc_median, peep_median, pplat_median, frac_dp_over_15) %>%
                      mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)

# DP against VT/PFVC within each VT/PBW bin: binned means (1-point bins of VT/PFVC) with 95% intervals,
# and the within-bin linear slope of DP on VT/PFVC (cmH2O per point)
binned <- pb %>% mutate(vtpfvc_bin = floor(vtpfvc_pct)) %>%
  group_by(vtpbw_bin, vtpfvc_bin) %>%
  summarise(n = n(), dp_mean = mean(dp), dp_se = sd(dp) / sqrt(n()), .groups = "drop") %>%
  filter(n >= 10)
slopes <- pb %>% group_by(vtpbw_bin) %>% filter(n() >= 30) %>%
  group_modify(~ { f <- lm(dp ~ vtpfvc_pct, data = .x); co <- summary(f)$coefficients
    tibble(n = nrow(.x), slope = co[2, 1], slope_se = co[2, 2], intercept = co[1, 1],
           dp_at_9 = co[1, 1] + 9 * co[2, 1], dp_at_11 = co[1, 1] + 11 * co[2, 1]) }) %>% ungroup()
message("\nDP on VT/PFVC within each bin (cmH2O per point of predicted FVC), and fitted DP at 9% and 11%:")
print(as.data.frame(slopes %>% mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
write_csv(bind_rows(by_bin %>% mutate(table = "by_bin"),
                    binned %>% mutate(table = "dp_by_vtpfvc_bin", site = site_name),
                    slopes %>% mutate(table = "slopes", site = site_name)),
          file.path(final_dir, paste0("dp_by_vtpbw_", site_name, ".csv")))

# figure
p1 <- ggplot(binned, aes(vtpfvc_bin + 0.5, dp_mean, colour = vtpbw_bin)) +
  geom_pointrange(aes(ymin = dp_mean - 1.96 * dp_se, ymax = dp_mean + 1.96 * dp_se), size = 0.3) +
  geom_smooth(data = pb, aes(vtpfvc_pct, dp, colour = vtpbw_bin), method = "lm", se = FALSE, linewidth = 0.7) +
  geom_vline(xintercept = 11, linetype = 2, colour = "grey50") +
  scale_colour_manual(values = okabe, name = "VT/PBW (mL/kg)") +
  coord_cartesian(xlim = c(4, 20)) +
  labs(title = "Driving pressure against VT/PFVC, by VT/PBW bin",
       subtitle = "binned means with 95% intervals (bins with >= 10 patients) and the within-bin linear fit; dashed = 11% ceiling",
       x = "VT/PFVC (% of predicted FVC)", y = "driving pressure (cmH2O)")
sumplot <- function(var, ylab, title) {
  d <- pb %>% filter(is.finite(.data[[var]]))
  ggplot(d, aes(vtpbw_bin, .data[[var]], colour = vtpbw_bin)) +
    stat_summary(fun = median, fun.min = function(x) quantile(x, .25), fun.max = function(x) quantile(x, .75),
                 geom = "pointrange") +
    stat_summary(fun = median, geom = "text", aes(label = after_stat(sprintf("%.1f", y))), vjust = -0.8, size = 3, show.legend = FALSE) +
    scale_colour_manual(values = okabe, guide = "none") +
    labs(title = title, x = "VT/PBW (mL/kg)", y = ylab)
}
p2 <- sumplot("dp", "cmH2O", "Driving pressure (median, IQR)")
p3 <- sumplot("crs", "mL / cmH2O", "Respiratory-system compliance")
p4 <- sumplot("vtpfvc_pct", "% of predicted FVC", "VT/PFVC")
p5 <- sumplot("peep", "cmH2O", "PEEP")
counts_lab <- paste(sprintf("%s: %d pts", by_bin$vtpbw_bin, by_bin$n_patients), collapse = "   ")
p <- p1 / (p2 + p3 + p4 + p5 + plot_layout(nrow = 1)) + plot_layout(heights = c(2, 1), guides = "collect") +
  plot_annotation(title = sprintf("What a low VT/PBW means (%s)", site_name),
                  subtitle = paste0("complete-data ventilator timepoints with a recorded plateau, one row per patient x bin;  ", counts_lab)) &
  theme_minimal(base_size = 10)
ggsave(file.path(final_dir, paste0("dp_by_vtpbw_", site_name, ".pdf")), p, width = 12, height = 9)
message("injury_dp_by_vtpbw complete -> ", final_dir)
