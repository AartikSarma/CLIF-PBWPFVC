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
#   ventilated                 the analytic cohort (figure 4)
#   ventilated, at ICU admission  the ventilated patients on IMV at ICU admission
#                              ("day0_" tables): the ventilated side of every
#                              difference-in-differences, so that both arms are
#                              assigned their status at the same moment, ICU admission
#   ventilated, SF band        the same, within an index SF class (PBWPFVC_JM_SF_BAND)
#   no support                 the negative control: no positive pressure, no strain.
#                              Every patient is kept, and the divergence is read at
#                              the ventilated cohort's severity (PBWPFVC_JM_SEV_CENTER,
#                              the "sevstd_" tables); the severity x divergence term
#                              says whether sicker controls diverge faster. Follow-up
#                              ends at escalation to invasive ventilation, NIPPV or
#                              another advanced support (the competing event), and no
#                              patient is in both arms: 03 removes from the control
#                              every patient of the ventilated ICU-admission arm.
#   no support, SF < 315       the same control restricted to patients hypoxemic at
#                              the index ("sevstd_sf0to315_" tables), so that it
#                              differs from the ventilated cohort in ventilation and
#                              not in hypoxemia; its DiD is written separately
# Every arm runs on one clock, days from its index (the first qualifying ventilator
# row; ICU admission in the control), with delayed entry at the first trajectory day.
# These are the arms 29_run_figure4.R fits. Only these arm tags are read (arm_pattern
# below); any other tables in the folders are ignored.
#
# Both differences-in-differences are written and pooled: the hypoxemic one (SF < 315)
# isolates ventilation, the all-patients one is the larger sample.
#
# Figure 4 = PBWPFVC_JM_MODIFIER=pfvc, daily grid, 7 days, as set by 29_run_figure4.R.
#
# Convergence: the figure's estimates are gated on the lung-size terms (size_gate, the
# level and divergence with R-hat <= 1.1, 20_biotrauma_grid.R); both_converged is that
# gate in both arms. The hazard-link convergence is reported (hazard_rhat) and read with
# the longitudinal-only comparison: the same difference-in-differences from each arm's
# longitudinal submodel fitted alone (did_lme_*, from jm_lme_check_*).
#
# Vasopressors are a two-part (hurdle) outcome: on/off (any_pressor, every patient-day,
# a logistic mixed model, its rate in log-odds per day) and the dose on the days a
# pressor runs (pressor_dose). Being on a pressor is itself an outcome, so the dose part
# is conditional on it and is never compared with a control: the difference-in-
# differences uses the on/off part. Rescaling a control's rate by the ratio of the SDs
# of log PFVC is the same for the log-odds of a yes/no marker as for a log marker.
#
# A control arm is informative only if its marker MOVES. The movement panel (mean
# change from baseline by day, from jm_movement_*) is therefore read before the
# divergence: no movement, no possible divergence, and the arm cannot adjudicate.
#
# Arms are discovered from the files present: the ventilated cohort's tables in
# final/injury/, the control cohort's in final/controls/ (file names tagged
# {site}_nosupport); each restriction carries its own tag. Run with
# PBWPFVC_COHORT unset.
#
# Outputs, in final/injury/:
#   jm_control_comparison_{form}_{h}_{site}.csv    each arm's divergence, per SD of log
#                                                  PFVC in that arm's own cohort
#   jm_control_comparison_{form}_{h}_{site}.pdf
#   jm_control_did_{form}_{h}_{site}.csv            ventilated at ICU admission minus the control,
#                                                  from the joint models (did_*) and from the
#                                                  longitudinal submodels alone (did_lme_*)
#   jm_hypoxemic_control_did_{form}_{h}_{site}.csv  the same against the hypoxemic control
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
final_dir <- final_dir_for("injury")              # the ventilated tables; the controls sit in final/controls/
okabe <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#000000", "#F0E442")

# the arm tags of the current design, per cohort (regular expressions on the restriction tag)
cohort_folders <- tibble(cohort = c("imv", "nosupport"),
                         cohort_label = c("Ventilated", "No support, at ventilated severity"),
                         site = c(base_site, paste0(base_site, "_nosupport")),
                         arm_pattern = c("^(day0_)?(sf[0-9.]+to[0-9.]+_)?$", "^sevstd_(sf0to315_)?$"))

# ---- discover the arms: one per (cohort folder, restriction tag) with an estimates table
file_stub <- function(site) paste0(MOD_FORM, "_", h_suffix, "_", site, ".csv")
arms <- pmap_dfr(cohort_folders, function(cohort, cohort_label, site, arm_pattern) {
  folder <- if (cohort == "imv") final_dir else file.path(config$final_root, "controls")
  found <- list.files(folder, pattern = paste0("^jm_estimates_.*", file_stub(site), "$"))
  restriction <- sub(paste0(file_stub(site), "$"), "", sub("^jm_estimates_", "", found))
  # the rrtcause_, offset_ and nolag_ variants are other analyses, not arms
  tibble(cohort, cohort_label, site, folder, restriction = restriction[grepl(arm_pattern, restriction)])
})
if (!nrow(arms)) stop("no jm_estimates_*", MOD_FORM, "_", h_suffix, "_* tables found for ", base_site)
arms <- arms %>%
  mutate(sf_part  = str_match(restriction, "sf([0-9.]+)to([0-9.]+)_")[, 2:3, drop = FALSE] %>%
           apply(1, function(limits) if (anyNA(limits)) NA_character_ else
             if (as.numeric(limits[1]) == 0) paste0("SF < ", limits[2]) else paste0("SF ", limits[1], "-", limits[2])),
         at_icu_admission = grepl("^day0_", restriction),
         arm = paste0(cohort_label, if_else(at_icu_admission, ", at ICU admission", ""),
                      if_else(is.na(sf_part), "", paste0(", ", sf_part))))
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
lme_check <- each_arm("lme_check")
# the gate of every fit in every arm, from its own estimates (20_biotrauma_grid.R)
convergence <- fit_convergence(estimates, MOD_FORM, by = c("arm", "marker", "model", "adjustment"))

# the unit of a marker's rate: log-odds for the yes/no marker, the log marker otherwise
marker_scale <- function(marker) if_else(marker == "any_pressor", "log-odds of any vasopressor", "log marker")
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
movement_summary <- if (nrow(movement)) movement %>%
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
                     any_of(c("sev_center", "sev_anchor", "sf_band", "n_iter"))),
            by = c("arm", "marker", "adjustment")) %>%
  left_join(movement_summary, by = c("arm", "marker", "adjustment")) %>%
  left_join(convergence %>% filter(model == "main") %>% select(arm, marker, adjustment, size_terms_rhat, size_gate, hazard_rhat),
            by = c("arm", "marker", "adjustment")) %>%
  mutate(arm = factor(arm, levels = arms$arm),
         unit = paste(marker_scale(marker), "per day per SD of log PFVC in this arm's cohort"),
         form = MOD_FORM, panel = h_suffix, site = base_site) %>%
  arrange(marker, adjustment, arm)
out_stub <- paste0(MOD_FORM, "_", h_suffix, "_", base_site)
write_csv(comparison, file.path(final_dir, paste0("jm_control_comparison_", out_stub, ".csv")))

# ---- the difference-in-differences: ventilated divergence (patients on IMV at ICU
#      admission) minus control divergence
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
# The ventilated side is the ICU-admission arm, not the whole ventilated cohort: both
# arms are then defined by their support at ICU admission, and no patient is in both.
# Units: each cohort's rate is per SD of log PFVC in ITS OWN panel (22 writes the SDs
# to jm_scale_*), so the control's rate is put on the ventilated cohort's unit before
# the subtraction: per ventilated SD = per control SD x (SD ventilated / SD control).
# The ICU-admission arm's rate is already per SD of the whole ventilated panel.
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
# control (index SF < 315, the ventilated cohort's own gate; PBWPFVC_JM_SF_BAND=0,315
# on the control fit). The ventilated cohort is hypoxemic by construction and
# the whole control mostly is not, so the primary difference also contrasts hypoxemia;
# against the hypoxemic control the arms differ in ventilation alone. Both use the
# control's whole-panel SD of log PFVC (a restricted fit keeps it, 22_biotrauma_fit.R).
# No difference exists without both sides (the SD rescaling missing, or no fit on
# one side): that is an empty table, which the callers report, not an error.
VENTILATED_DID_ARM <- "Ventilated, at ICU admission"
# the dose part is conditional on being on a pressor (an outcome) and is not compared
# with a control, even where an older run left its rows in an arm's tables
CONDITIONAL_MARKERS <- "pressor_dose"
# The same difference from the longitudinal submodels fitted alone (each arm's
# jm_lme_check_*: maximum likelihood, no correction for patients leaving the panel), so
# the difference-in-differences can be read without the hazard links. The LME's
# standard errors take the place of the posterior SDs; the control is put on the
# ventilated unit in the same way.
lme_did_against <- function(control_arm) {
  none <- tibble(marker = character(), adjustment = character())
  sides <- if (is.null(lme_check) || !nrow(lme_check)) tibble() else lme_check %>%
    filter(model == "main", exposure == "log_pfvc_sd", term == "divergence per day",
           arm %in% c(VENTILATED_DID_ARM, control_arm), !marker %in% CONDITIONAL_MARKERS)
  missing <- setdiff(c(VENTILATED_DID_ARM, control_arm), unique(sides$arm))
  if (length(missing)) {
    message("--- no longitudinal-only difference-in-differences against ", control_arm, ": no jm_lme_check_* divergence for ",
            paste(missing, collapse = " and "))
    return(none)
  }
  sides %>%
    mutate(side = if_else(arm == VENTILATED_DID_ARM, "ventilated", "control"),
           unit_factor = if_else(side == "control", to_vent_sd, 1),
           lme_estimate = lme_estimate * unit_factor, lme_se = lme_se * unit_factor) %>%
    select(marker, adjustment, side, lme_estimate, lme_se) %>%
    pivot_wider(names_from = side, values_from = c(lme_estimate, lme_se)) %>%
    filter(!is.na(lme_estimate_ventilated), !is.na(lme_estimate_control)) %>%
    transmute(marker, adjustment,
              divergence_lme_ventilated = lme_estimate_ventilated, divergence_lme_control = lme_estimate_control,
              did_lme_estimate = lme_estimate_ventilated - lme_estimate_control,
              did_lme_sd = sqrt(lme_se_ventilated^2 + lme_se_control^2),
              did_lme_lo = did_lme_estimate - 1.96 * did_lme_sd, did_lme_hi = did_lme_estimate + 1.96 * did_lme_sd)
}
did_against <- function(control_arm) {
  both_arms <- comparison %>%
    filter(!is.na(to_vent_sd), arm %in% c(VENTILATED_DID_ARM, control_arm), !marker %in% CONDITIONAL_MARKERS) %>%
    mutate(side = if_else(cohort == "imv", "ventilated", "control"))
  if (!all(c("ventilated", "control") %in% both_arms$side)) return(tibble())
  both_arms %>%
  # the control on the ventilated cohort's unit (the ventilated rows are multiplied by 1)
  mutate(unit_factor = if_else(side == "control", to_vent_sd, 1),
         divergence_estimate = divergence_estimate * unit_factor, divergence_sd = divergence_sd * unit_factor) %>%
  select(marker, adjustment, side, divergence_estimate, divergence_sd, divergence_rhat, size_terms_rhat, size_gate,
         n_patients, any_of("creatinine_model")) %>%
  pivot_wider(names_from = side, values_from = c(divergence_estimate, divergence_sd, divergence_rhat, size_terms_rhat, size_gate,
                                                 n_patients, any_of("creatinine_model"))) %>%
  filter(!is.na(divergence_estimate_ventilated), !is.na(divergence_estimate_control)) %>%
  rename(ventilated_size_rhat = size_terms_rhat_ventilated, control_size_rhat = size_terms_rhat_control) %>%
  mutate(did_estimate = divergence_estimate_ventilated - divergence_estimate_control,
         did_sd = sqrt(divergence_sd_ventilated^2 + divergence_sd_control^2),
         did_lo = did_estimate - 1.96 * did_sd, did_hi = did_estimate + 1.96 * did_sd,
         p_did_gt0 = pnorm(did_estimate / did_sd),
         both_converged = size_gate_ventilated & size_gate_control,
         control_to_ventilated_sd = to_vent_sd, ventilated_arm = VENTILATED_DID_ARM, control_arm = control_arm,
         unit = paste(marker_scale(marker), "per day per SD of log PFVC in the ventilated cohort"), form = MOD_FORM, panel = h_suffix, site = base_site) %>%
  select(-size_gate_ventilated, -size_gate_control) %>%
  left_join(lme_did_against(control_arm), by = c("marker", "adjustment"))
}
did <- did_against("No support, at ventilated severity")
HYPOXEMIC_CONTROL_ARM <- "No support, at ventilated severity, SF < 315"
hypoxemic_did_path <- file.path(final_dir, paste0("jm_hypoxemic_control_did_", out_stub, ".csv"))
if (HYPOXEMIC_CONTROL_ARM %in% arms$arm) {
  hypoxemic_did <- did_against(HYPOXEMIC_CONTROL_ARM)
  if (nrow(hypoxemic_did)) {
    write_csv(hypoxemic_did, hypoxemic_did_path)
    message("--- difference-in-differences, ventilated at ICU admission against the hypoxemic control (index SF < 315)")
    print(as.data.frame(hypoxemic_did %>% transmute(marker, adjustment, ventilated = signif(divergence_estimate_ventilated, 3),
                                                    control = signif(divergence_estimate_control, 3), did = signif(did_estimate, 3),
                                                    lo = signif(did_lo, 3), hi = signif(did_hi, 3),
                                                    did_lme = signif(did_lme_estimate, 3), both_converged)), row.names = FALSE)
  } else unlink(hypoxemic_did_path)
} else unlink(hypoxemic_did_path)
if (nrow(did)) {
  write_csv(did, did_path)
  message("--- difference-in-differences: ventilated at ICU admission minus no-support divergence (log marker, or log-odds for any_pressor, per day per SD of log PFVC)")
  print(as.data.frame(did %>% transmute(marker, adjustment, ventilated = signif(divergence_estimate_ventilated, 3),
                                        control = signif(divergence_estimate_control, 3), did = signif(did_estimate, 3),
                                        lo = signif(did_lo, 3), hi = signif(did_hi, 3), p_did_gt0 = signif(p_did_gt0, 3),
                                        did_lme = signif(did_lme_estimate, 3), did_lme_lo = signif(did_lme_lo, 3),
                                        did_lme_hi = signif(did_lme_hi, 3), both_converged)), row.names = FALSE)
} else {
  message("--- no marker has both a ventilated ICU-admission fit (day0_) and a severity-standardised control fit: no difference-in-differences")
  unlink(did_path)
}

message("--- divergence per day per SD of log PFVC (log marker units; log-odds for any_pressor), adjusted")
print(as.data.frame(comparison %>% filter(adjustment == "adjusted") %>%
                      transmute(marker, arm, n_patients, n_deaths, divergence = signif(divergence_estimate, 3),
                                lo = signif(divergence_lo, 3), hi = signif(divergence_hi, 3),
                                size_rhat = round(size_terms_rhat, 2), hazard_rhat = round(hazard_rhat, 2),
                                sev_mod = if ("sev_modification_estimate" %in% names(comparison)) signif(sev_modification_estimate, 3) else NA_real_,
                                sev_mod_lo = if ("sev_modification_lo" %in% names(comparison)) signif(sev_modification_lo, 3) else NA_real_,
                                sev_mod_hi = if ("sev_modification_hi" %in% names(comparison)) signif(sev_modification_hi, 3) else NA_real_,
                                peak_change = signif(movement_peak_mean_change, 3), peak_day = movement_peak_day,
                                mean_abs_change = signif(movement_mean_abs_change, 3))), row.names = FALSE)

# ---- figure: the divergence by arm, above how much the marker moves in each arm
arm_colours <- setNames(okabe[seq_len(nrow(arms))], arms$arm)
if (nrow(arms) > length(okabe)) stop("more arms than Okabe-Ito colours; restrict the arms before plotting")
forest <- ggplot(comparison, aes(divergence_estimate, fct_rev(arm), colour = arm, shape = adjustment, linetype = size_gate)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = divergence_lo, xmax = divergence_hi), position = position_dodge(width = 0.6)) +
  facet_wrap(~ marker, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = arm_colours, guide = "none") +
  scale_shape_manual(values = c(adjusted = 16, unadjusted = 1), name = NULL) +
  scale_linetype_manual(values = c(`TRUE` = 1, `FALSE` = 3), labels = c(`TRUE` = "size terms R-hat <= 1.1", `FALSE` = "size terms R-hat > 1.1"), name = NULL) +
  labs(title = "Divergence by lung size, arm by arm",
       subtitle = paste("Change in the log marker (log-odds for any vasopressor) per day, per SD of log predicted FVC (95% credible interval);",
                        "each arm on its own cohort's SD; the difference-in-differences table rescales the control"),
       x = "log marker (log-odds for any vasopressor) per day per SD", y = NULL) +
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
