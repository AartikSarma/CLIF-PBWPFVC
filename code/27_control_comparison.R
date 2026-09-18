# =============================================================================
# Script 27 (control comparison): the divergence by lung size, arm by arm
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# The biotrauma claim is that a smaller predicted FVC marks a lung more vulnerable
# to the ventilator, so the organ markers of patients with small lungs DIVERGE
# from those with large lungs over the first week (the log_pfvc_sd:vent_day term
# of the "pfvc" joint model). This script lines that term up across the arms that
# test it, reading only the aggregate tables 23_biotrauma_report.R wrote:
#
#   ventilated                 the analytic cohort
#   ventilated, SF band        the same, within a baseline SF class (PBWPFVC_JM_SF_BAND)
#   ventilated, anchor floor   the same, above the severity floor (PBWPFVC_JM_SEV_MIN)
#   no support                 the negative control: no positive pressure, no strain
#   no support, anchor floor   the control restricted to the ventilated severity range
#   noninvasive                NOT a control: uncontrolled tidal volumes, so a point on
#                              the strain gradient between the other two
#
# A control arm is informative only if its marker MOVES. The movement panel (mean
# change from baseline by day, from jm_movement_*) is therefore read before the
# divergence: no movement, no possible divergence, and the arm cannot adjudicate.
#
# Arms are discovered from the files present. Each cohort lives in its own site
# folder ({site}, {site}_nosupport, {site}_niv); each restriction carries its own
# tag. Run with the base config (no PBWPFVC_COHORT / PBWPFVC_SITE_NAME override).
#
# Outputs, in the base site's final/ folder:
#   jm_control_comparison_{form}_{tag}_{site}.csv
#   jm_control_movement_{form}_{tag}_{site}.csv
#   jm_control_comparison_{form}_{tag}_{site}.pdf
#
# Usage: PBWPFVC_JM_GRID=daily PBWPFVC_JM_HORIZON=7 Rscript code/27_control_comparison.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(here)
})
rm(list = ls())
source("utils/config.R")
if (config$cohort != "imv")
  stop("27_control_comparison.R reads every cohort's tables itself: unset PBWPFVC_COHORT")
base_site <- config$site_name
source(here("code", "20_biotrauma_grid.R"))   # h_suffix
MOD_FORM  <- Sys.getenv("PBWPFVC_JM_MODIFIER", "pfvc")
stopifnot(MOD_FORM %in% c("pfvc", "channels"))
RHAT_GATE <- 1.1
final_dir <- config$final_dir                     # the site's final/; the controls sit in final/controls/
okabe <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#000000", "#F0E442")

cohort_folders <- tibble(cohort = c("imv", "niv", "nosupport"),
                         cohort_label = c("Ventilated", "Noninvasive", "No support"),
                         site = c(base_site, paste0(base_site, "_niv"), paste0(base_site, "_nosupport")))

# ---- discover the arms: one per (cohort folder, restriction tag) with an estimates table
file_stub <- function(site) paste0(MOD_FORM, "_", h_suffix, "_", site, ".csv")
arms <- pmap_dfr(cohort_folders, function(cohort, cohort_label, site) {
  folder <- if (cohort == "imv") final_dir else file.path(final_dir, "controls")
  found <- list.files(folder, pattern = paste0("^jm_estimates_.*", file_stub(site), "$"))
  restriction <- sub(paste0(file_stub(site), "$"), "", sub("^jm_estimates_", "", found))
  # only the cohort-restriction tags; the rrtcause_ and offset_ variants are other analyses
  is_arm <- grepl("^(sev[0-9.]+_|sev(_[a-z_]+?[0-9.]+)+_)?(sf[0-9.]+to[0-9.]+_)?$", restriction, perl = TRUE)
  tibble(cohort, cohort_label, site, folder, restriction = restriction[is_arm])
})
if (!nrow(arms)) stop("no jm_estimates_*", MOD_FORM, "_", h_suffix, "_* tables found for ", base_site)
arms <- arms %>%
  mutate(sev_only = sub("sf[0-9.]+to[0-9.]+_$", "", restriction),      # the SF band is labelled separately
         sev_part = if_else(grepl("^sev_", restriction),
                            # per-marker floors, "sev_bilirubin1_platelets2_" -> "bilirubin >= 1, platelets >= 2"
                            sev_only %>% str_remove("^sev_") %>%
                              str_replace_all("([a-z_]+?)([0-9.]+)(_|$)", "\\1 >= \\2, ") %>% str_remove(", $"),
                            paste0("anchor >= ", str_match(restriction, "^sev([0-9.]+)_")[, 2])),
         sev_part = if_else(grepl("^sev", restriction), sev_part, NA_character_),
         sf_part  = str_match(restriction, "sf([0-9.]+)to([0-9.]+)_")[, 2:3, drop = FALSE] %>%
           apply(1, function(limits) if (anyNA(limits)) NA_character_ else
             if (as.numeric(limits[1]) == 0) paste0("SF <= ", limits[2]) else paste0("SF ", limits[1], "-", limits[2])),
         arm = paste0(cohort_label,
                      if_else(is.na(sf_part), "", paste0(", ", sf_part)),
                      if_else(is.na(sev_part), "", paste0(", matched (", sev_part, ")"))))
message("=== 27_control_comparison (", MOD_FORM, ", ", h_suffix, ", ", base_site, "): ", nrow(arms), " arms ===")
message(paste0("  ", arms$arm, collapse = "\n"))

read_arm <- function(table_name, arm_row) {
  path <- file.path(arm_row$folder, paste0("jm_", table_name, "_", arm_row$restriction, file_stub(arm_row$site)))
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>% mutate(arm = arm_row$arm, cohort = arm_row$cohort)
}
each_arm <- function(table_name) map_dfr(seq_len(nrow(arms)), function(arm_i) read_arm(table_name, arms[arm_i, ]))

estimates <- each_arm("estimates")
manifest  <- each_arm("manifest")
movement  <- each_arm("movement")

# ---- the comparison table: divergence and level per SD of log PFVC, with what qualifies them
size_terms <- estimates %>%
  filter(model == "main", block == "longitudinal", term %in% c("log_pfvc_sd", "log_pfvc_sd:vent_day")) %>%
  mutate(quantity = if_else(term == "log_pfvc_sd", "level", "divergence")) %>%
  select(arm, cohort, marker, adjustment, quantity, estimate, lo, hi, rhat) %>%
  pivot_wider(names_from = quantity, values_from = c(estimate, lo, hi, rhat), names_glue = "{quantity}_{.value}")
# Movement is summarised over ALL days, not the last one: platelets fall and recover
# inside a week, and a marker that moved and came back would read as still on day 7.
last_movement <- if (nrow(movement)) movement %>%
  filter(model == "main") %>%
  group_by(arm, marker, adjustment) %>% arrange(day, .by_group = TRUE) %>%
  summarise(movement_peak_mean_change = mean_change[which.max(abs(mean_change))],
            movement_peak_day = day[which.max(abs(mean_change))],
            movement_last_day = last(day), movement_last_n = last(n_patients),
            movement_last_mean_change = last(mean_change),
            movement_mean_abs_change = mean(mean_abs_change), movement_mean_sd_change = mean(sd_change),
            .groups = "drop") else
  tibble(arm = character(), marker = character(), adjustment = character())
comparison <- size_terms %>%
  left_join(manifest %>% filter(model == "main") %>%
              select(arm, marker, adjustment, n_patients, n_deaths, n_competing = n_extubations, status,
                     hazard_rhat, any_of(c("sev_floor", "sev_anchor", "sf_band", "n_iter"))),
            by = c("arm", "marker", "adjustment")) %>%
  left_join(last_movement, by = c("arm", "marker", "adjustment")) %>%
  mutate(divergence_converged = divergence_rhat <= RHAT_GATE,
         arm = factor(arm, levels = arms$arm), form = MOD_FORM, panel = h_suffix, site = base_site) %>%
  arrange(marker, adjustment, arm)
stopifnot(all(comparison$n_patients >= 10L, na.rm = TRUE))
out_stub <- paste0(MOD_FORM, "_", h_suffix, "_", base_site)
write_csv(comparison, file.path(final_dir, paste0("jm_control_comparison_", out_stub, ".csv")))
if (nrow(movement)) write_csv(movement %>% filter(model == "main"),
                              file.path(final_dir, paste0("jm_control_movement_", out_stub, ".csv")))

message("--- divergence per day per SD of log PFVC (log marker units), adjusted")
print(as.data.frame(comparison %>% filter(adjustment == "adjusted") %>%
                      transmute(marker, arm, n_patients, n_deaths, divergence = signif(divergence_estimate, 3),
                                lo = signif(divergence_lo, 3), hi = signif(divergence_hi, 3),
                                rhat = round(divergence_rhat, 2), hazard_rhat = round(hazard_rhat, 2),
                                peak_change = signif(movement_peak_mean_change, 3), peak_day = movement_peak_day,
                                mean_abs_change = signif(movement_mean_abs_change, 3))), row.names = FALSE)

# ---- figure: the divergence by arm, above how much the marker moves in each arm
arm_colours <- setNames(okabe[seq_len(nrow(arms))], arms$arm)
if (nrow(arms) > length(okabe)) stop("more arms than Okabe-Ito colours; restrict the arms before plotting")
forest <- ggplot(comparison, aes(divergence_estimate, fct_rev(arm), colour = arm, shape = adjustment, linetype = divergence_converged)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = divergence_lo, xmax = divergence_hi), position = position_dodge(width = 0.6)) +
  facet_wrap(~ marker, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = arm_colours, guide = "none") +
  scale_shape_manual(values = c(adjusted = 16, unadjusted = 1), name = NULL) +
  scale_linetype_manual(values = c(`TRUE` = 1, `FALSE` = 3), labels = c(`TRUE` = "R-hat <= 1.1", `FALSE` = "R-hat > 1.1"), name = NULL) +
  labs(title = "Divergence by lung size, arm by arm",
       subtitle = "Change in the log marker per day, per SD of log predicted FVC (95% credible interval)",
       x = "log marker per day per SD", y = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
panels <- list(forest)
if (nrow(movement)) {
  movement_plot <- movement %>%
    filter(model == "main", adjustment == "adjusted", !binary) %>%
    mutate(arm = factor(arm, levels = arms$arm)) %>%
    ggplot(aes(day, mean_change, colour = arm)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
    facet_wrap(~ marker, scales = "free_y", nrow = 1) +
    scale_colour_manual(values = arm_colours, name = NULL) +
    scale_x_continuous(breaks = seq(1, 28)) +
    labs(title = "Does the marker move in this arm?",
         subtitle = "Observed mean change from baseline among patients still observed (log units; days with under 10 patients dropped)",
         x = "Day", y = "Mean change in log marker") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  panels <- c(panels, list(movement_plot))
}
ggsave(file.path(final_dir, paste0("jm_control_comparison_", out_stub, ".pdf")),
       wrap_plots(panels, ncol = 1), width = 12, height = 4.5 * length(panels))
message("27_control_comparison complete: ", nrow(comparison), " rows -> ", final_dir)
