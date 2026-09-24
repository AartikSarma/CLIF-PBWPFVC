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
#   no support                 the negative control: no positive pressure, no strain.
#                              Every patient is kept, and the divergence is read at
#                              the ventilated cohort's severity (PBWPFVC_JM_SEV_CENTER,
#                              the "sevstd_" tables); the severity x divergence term
#                              says whether sicker controls diverge faster.
#   no support, SF <= 315      the same control restricted to patients hypoxemic on
#                              the index day ("sevstd_sf0to315_" tables), so that it
#                              differs from the ventilated cohort in ventilation and
#                              not in hypoxemia; its DiD is written separately
# These are the arms 29_run_figure4.R fits. Tables from designs it no longer runs
# (severity floors, the unstandardised control, the noninvasive cohort, which is not a
# control) are ignored even when present on disk, so an old run cannot add an arm.
#
# A control arm is informative only if its marker MOVES. The movement panel (mean
# change from baseline by day, from jm_movement_*) is therefore read before the
# divergence: no movement, no possible divergence, and the arm cannot adjudicate.
#
# Arms are discovered from the files present: the ventilated cohort's tables in
# final/injury/, the control cohorts' in final/controls/ (file names tagged
# {site}_nosupport, {site}_niv); each restriction carries its own tag. Run with
# PBWPFVC_COHORT unset.
#
# Outputs, in final/injury/:
#   jm_control_comparison_{form}_{tag}_{site}.csv
#   jm_control_comparison_{form}_{tag}_{site}.pdf
#   jm_control_did_{form}_{h}_{site}.csv            ventilated minus the control
#   jm_hypoxemic_control_did_{form}_{h}_{site}.csv  ventilated minus the hypoxemic control
#
# Usage: PBWPFVC_JM_GRID=daily PBWPFVC_JM_HORIZON=7 uvr run code/27_control_comparison.R
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
final_dir <- final_dir_for("injury")              # the ventilated tables; the controls sit in final/controls/
okabe <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#000000", "#F0E442")

# the arm tags of the current design, per cohort (regular expressions on the restriction tag)
cohort_folders <- tibble(cohort = c("imv", "nosupport"),
                         cohort_label = c("Ventilated", "No support, at ventilated severity"),
                         site = c(base_site, paste0(base_site, "_nosupport")),
                         arm_pattern = c("^(sf[0-9.]+to[0-9.]+_)?$", "^sevstd_(sf0to315_)?$"))

# ---- discover the arms: one per (cohort folder, restriction tag) with an estimates table
file_stub <- function(site) paste0(MOD_FORM, "_", h_suffix, "_", site, ".csv")
arms <- pmap_dfr(cohort_folders, function(cohort, cohort_label, site, arm_pattern) {
  folder <- if (cohort == "imv") final_dir else file.path(config$final_root, "controls")
  found <- list.files(folder, pattern = paste0("^jm_estimates_.*", file_stub(site), "$"))
  restriction <- sub(paste0(file_stub(site), "$"), "", sub("^jm_estimates_", "", found))
  # the rrtcause_ and offset_ variants are other analyses, not arms
  tibble(cohort, cohort_label, site, folder, restriction = restriction[grepl(arm_pattern, restriction)])
})
if (!nrow(arms)) stop("no jm_estimates_*", MOD_FORM, "_", h_suffix, "_* tables found for ", base_site)
arms <- arms %>%
  mutate(sf_part  = str_match(restriction, "sf([0-9.]+)to([0-9.]+)_")[, 2:3, drop = FALSE] %>%
           apply(1, function(limits) if (anyNA(limits)) NA_character_ else
             if (as.numeric(limits[1]) == 0) paste0("SF <= ", limits[2]) else paste0("SF ", limits[1], "-", limits[2])),
         arm = paste0(cohort_label, if_else(is.na(sf_part), "", paste0(", ", sf_part))))
message("=== 27_control_comparison (", MOD_FORM, ", ", h_suffix, ", ", base_site, "): ", nrow(arms), " arms ===")
message(paste0("  ", arms$arm, collapse = "\n"))

# Each arm's table, with creatinine taken from its dialysis-as-third-cause fit
# (the rrtcause_ twin, as 29_run_figure4.R fits it) whenever that twin exists.
read_arm <- function(table_name, arm_row) {
  path <- file.path(arm_row$folder, paste0("jm_", table_name, "_", arm_row$restriction, file_stub(arm_row$site)))
  twin <- file.path(arm_row$folder, paste0("jm_", table_name, "_rrtcause_", arm_row$restriction, file_stub(arm_row$site)))
  as_text <- function(d) d %>% mutate(across(everything(), as.character))
  main <- if (file.exists(path)) read_csv(path, show_col_types = FALSE) else NULL
  rrt  <- if (file.exists(twin)) read_csv(twin, show_col_types = FALSE) %>% filter(marker == "creatinine") else NULL
  if (!is.null(main) && !is.null(rrt) && nrow(rrt)) main <- main %>% filter(marker != "creatinine")
  out <- bind_rows(if (!is.null(main)) as_text(main), if (!is.null(rrt)) as_text(rrt))
  if (!nrow(out)) return(NULL)
  has_rrt <- !is.null(rrt) && nrow(rrt) > 0
  out %>% type_convert(guess_integer = TRUE, na = c("", "NA")) %>%
    mutate(arm = arm_row$arm, cohort = arm_row$cohort,
           creatinine_model = if_else(marker == "creatinine" & has_rrt, "dialysis as a third cause", NA_character_))
}
each_arm <- function(table_name) map_dfr(seq_len(nrow(arms)), function(arm_i) read_arm(table_name, arms[arm_i, ]))

estimates <- each_arm("estimates")
manifest  <- each_arm("manifest")
movement  <- each_arm("movement")

# ---- the comparison table: divergence and level per SD of log PFVC, with what qualifies them
# sev_modification = the control's log_pfvc_sd:vent_day:sev_anchor_c: the change in the
# divergence per point of the severity anchor (NA in arms without it)
components <- function(term) vapply(strsplit(term, ":"), function(p) paste(sort(p), collapse = ":"), character(1))
term_quantity <- c(log_pfvc_sd = "level", "log_pfvc_sd:vent_day" = "divergence",
                   "log_pfvc_sd:sev_anchor_c:vent_day" = "sev_modification")
size_terms <- estimates %>%
  filter(model == "main", block == "longitudinal") %>%
  mutate(quantity = unname(term_quantity[components(term)])) %>%
  filter(!is.na(quantity)) %>%
  select(arm, cohort, marker, adjustment, any_of("creatinine_model"), quantity, estimate, sd, lo, hi, rhat) %>%
  pivot_wider(names_from = quantity, values_from = c(estimate, sd, lo, hi, rhat), names_glue = "{quantity}_{.value}")
# Movement is summarised over ALL days, not the last one: platelets fall and recover
# inside a week, and a marker that moved and came back would read as still on day 7.
last_movement <- if (nrow(movement)) movement %>%
  # a yes/no marker (any vasopressor) has no change from baseline: its rows are empty
  filter(model == "main", is.finite(mean_change)) %>%
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
                     hazard_rhat, any_of(c("sev_center", "sev_anchor", "sf_band", "n_iter"))),
            by = c("arm", "marker", "adjustment")) %>%
  left_join(last_movement, by = c("arm", "marker", "adjustment")) %>%
  mutate(divergence_converged = divergence_rhat <= RHAT_GATE,
         arm = factor(arm, levels = arms$arm), form = MOD_FORM, panel = h_suffix, site = base_site) %>%
  arrange(marker, adjustment, arm)
out_stub <- paste0(MOD_FORM, "_", h_suffix, "_", base_site)
write_csv(comparison, file.path(final_dir, paste0("jm_control_comparison_", out_stub, ".csv")))

# ---- the difference-in-differences: ventilated divergence minus control divergence
# First difference: at a fixed VT/PBW the PBW formula, which omits age and race,
# assigns the strain, so the divergence by predicted lung size is not chosen by
# indication. Second difference: that lung-size variation also exists without a
# ventilator, where it can act only through the patient (demographics, organ reserve);
# the no-support control, read at the ventilated severity, estimates that path. The
# difference is the divergence the ventilator adds, per SD of log PFVC per day.
# Identifying assumption, the parallel-trends analogue: absent ventilation, lung size
# would shape the marker's trajectory equally in both cohorts at equal severity.
# The cohorts are different patients, so their posteriors are independent and the
# difference's SD is the root sum of squares (normal approximation to the posterior).
# Units: each cohort's rate is per SD of log PFVC in ITS OWN panel (22 writes the SDs
# to jm_scale_*), so the control's rate is put on the ventilated cohort's unit before
# the subtraction: per ventilated SD = per control SD x (SD ventilated / SD control).
# Without both SDs there is no difference-in-differences, never one that assumes them equal.
did_path <- file.path(final_dir, paste0("jm_control_did_", out_stub, ".csv"))
scale_vent <- file.path(final_dir, paste0("jm_scale_", h_suffix, "_", base_site, ".csv"))
scale_ctrl <- file.path(config$final_root, "controls", paste0("jm_scale_", h_suffix, "_", base_site, "_nosupport.csv"))
to_vent_sd <- if (file.exists(scale_vent) && file.exists(scale_ctrl))
  read_csv(scale_vent, show_col_types = FALSE)$sd_log_pfvc / read_csv(scale_ctrl, show_col_types = FALSE)$sd_log_pfvc else NA_real_
if (is.na(to_vent_sd)) {
  message("--- no difference-in-differences: the SD of log PFVC is missing for ",
          paste(c(if (!file.exists(scale_vent)) "the ventilated cohort", if (!file.exists(scale_ctrl)) "the control"), collapse = " and "),
          " (rerun 22_biotrauma_fit.R, e.g. its anchor step, for that cohort)")
  unlink(did_path)
}
# The same difference against any control arm: the primary one, and the hypoxemic
# control (index-day SF <= 315, the ventilated cohort's own gate; PBWPFVC_JM_SF_BAND=0,315
# on the control fit, 2026-09-23). The ventilated cohort is hypoxemic by construction and
# the whole control mostly is not, so the primary difference also contrasts hypoxemia;
# against the hypoxemic control the arms differ in ventilation alone. Both use the
# control's whole-panel SD of log PFVC (a restricted fit keeps it, 22_biotrauma_fit.R).
# No difference exists without both sides (the SD rescaling missing, or no fit on
# one side): that is an empty table, which the callers report, not an error.
did_against <- function(control_arm) {
  both_arms <- comparison %>%
    filter(!is.na(to_vent_sd), arm %in% c("Ventilated", control_arm)) %>%
    mutate(side = if_else(cohort == "imv", "ventilated", "control"))
  if (!all(c("ventilated", "control") %in% both_arms$side)) return(tibble())
  both_arms %>%
  # the control on the ventilated cohort's unit (the ventilated rows are multiplied by 1)
  mutate(unit_factor = if_else(side == "control", to_vent_sd, 1),
         divergence_estimate = divergence_estimate * unit_factor, divergence_sd = divergence_sd * unit_factor) %>%
  select(marker, adjustment, side, divergence_estimate, divergence_sd, divergence_rhat, n_patients, any_of("creatinine_model")) %>%
  pivot_wider(names_from = side, values_from = c(divergence_estimate, divergence_sd, divergence_rhat, n_patients, any_of("creatinine_model"))) %>%
  filter(!is.na(divergence_estimate_ventilated), !is.na(divergence_estimate_control)) %>%
  mutate(did_estimate = divergence_estimate_ventilated - divergence_estimate_control,
         did_sd = sqrt(divergence_sd_ventilated^2 + divergence_sd_control^2),
         did_lo = did_estimate - 1.96 * did_sd, did_hi = did_estimate + 1.96 * did_sd,
         p_did_gt0 = pnorm(did_estimate / did_sd),
         both_converged = divergence_rhat_ventilated <= RHAT_GATE & divergence_rhat_control <= RHAT_GATE,
         control_to_ventilated_sd = to_vent_sd, control_arm = control_arm,
         unit = "log marker per day per SD of log PFVC in the ventilated cohort", form = MOD_FORM, panel = h_suffix, site = base_site)
}
did <- did_against("No support, at ventilated severity")
HYPOXEMIC_CONTROL_ARM <- "No support, at ventilated severity, SF <= 315"
hypoxemic_did_path <- file.path(final_dir, paste0("jm_hypoxemic_control_did_", out_stub, ".csv"))
if (HYPOXEMIC_CONTROL_ARM %in% arms$arm) {
  hypoxemic_did <- did_against(HYPOXEMIC_CONTROL_ARM)
  if (nrow(hypoxemic_did)) {
    write_csv(hypoxemic_did, hypoxemic_did_path)
    message("--- difference-in-differences against the hypoxemic control (index-day SF <= 315)")
    print(as.data.frame(hypoxemic_did %>% transmute(marker, adjustment, ventilated = signif(divergence_estimate_ventilated, 3),
                                                    control = signif(divergence_estimate_control, 3), did = signif(did_estimate, 3),
                                                    lo = signif(did_lo, 3), hi = signif(did_hi, 3), both_converged)), row.names = FALSE)
  } else unlink(hypoxemic_did_path)
} else unlink(hypoxemic_did_path)
if (nrow(did)) {
  write_csv(did, did_path)
  message("--- difference-in-differences: ventilated minus no-support divergence (log marker per day per SD of log PFVC)")
  print(as.data.frame(did %>% transmute(marker, adjustment, ventilated = signif(divergence_estimate_ventilated, 3),
                                        control = signif(divergence_estimate_control, 3), did = signif(did_estimate, 3),
                                        lo = signif(did_lo, 3), hi = signif(did_hi, 3), p_did_gt0 = signif(p_did_gt0, 3),
                                        both_converged)), row.names = FALSE)
} else {
  message("--- no marker has both a ventilated and a severity-standardised control fit: no difference-in-differences")
  unlink(did_path)
}

message("--- divergence per day per SD of log PFVC (log marker units), adjusted")
print(as.data.frame(comparison %>% filter(adjustment == "adjusted") %>%
                      transmute(marker, arm, n_patients, n_deaths, divergence = signif(divergence_estimate, 3),
                                lo = signif(divergence_lo, 3), hi = signif(divergence_hi, 3),
                                rhat = round(divergence_rhat, 2), hazard_rhat = round(hazard_rhat, 2),
                                sev_mod = if ("sev_modification_estimate" %in% names(comparison)) signif(sev_modification_estimate, 3) else NA_real_,
                                sev_mod_lo = if ("sev_modification_lo" %in% names(comparison)) signif(sev_modification_lo, 3) else NA_real_,
                                sev_mod_hi = if ("sev_modification_hi" %in% names(comparison)) signif(sev_modification_hi, 3) else NA_real_,
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
         subtitle = "Observed mean change from baseline among patients still observed (log units)",
         x = "Day", y = "Mean change in log marker") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  panels <- c(panels, list(movement_plot))
}
ggsave(file.path(final_dir, paste0("jm_control_comparison_", out_stub, ".pdf")),
       wrap_plots(panels, ncol = 1), width = 12, height = 4.5 * length(panels))
message("27_control_comparison complete: ", nrow(comparison), " rows -> ", final_dir)
