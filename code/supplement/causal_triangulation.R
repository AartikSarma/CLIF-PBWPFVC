# =============================================================================
# Supplement (causal): does a marker's slope follow PBW/PFVC only while the
# patient is ventilated? Complementary designs, side by side
# =============================================================================
#
# At a fixed VT/PBW, VT/PFVC moves only with lambda = log(PBW/PFVC), a fixed function
# of age, sex, race and height. Adjusting for age, sex and race leaves lambda almost
# no variation (a positivity problem); not adjusting leaves every other effect of
# those demographics in the estimate. This script fits the designs of
# code/supplement/causal_common.R to the site's own platelet and creatinine
# trajectories. Each design estimates the same quantity from a different source of
# variation and fails under a different assumption, so agreement across them is the
# evidence, not any one row:
#
#   unadjusted, adjusted, fingerprint, control_did   the designs the paper already has,
#       refitted here on one model so they are comparable with the new ones
#   within_onoff    each ventilated patient's days before intubation against their
#                   days on the ventilator: demographics are fixed within a patient
#   within_triple   that change minus the controls' change around their ICU admission
#   off_reversal    on the ventilator against the days after extubation
#   free_<input>    the ladder: one of age, sex, race or height set free at a time,
#                   lambda identified by the other three; and each as a
#                   ventilated-minus-control contrast (_did)
#
# The estimate, in every row: the change in the marker's slope (log units per day)
# per log unit of PBW/PFVC that belongs to ventilation. Biotrauma predicts a negative
# value for platelets and a positive one for creatinine; injury_estimate carries the
# sign so that positive is always more injury. partial_r2 is the share of the
# identifying term's variation that the design's other terms take away (0.99 in the
# adjusted model is the positivity problem made visible; in a difference design part of
# it is the comparison phase or arm, which is the design and not a loss), and se_ratio
# is the design's standard error over the unadjusted model's: what its assumption costs.
#
# What the script does not do: it is the longitudinal model alone, conditional on a
# patient still being measured, with no joint model for death, liberation or
# escalation (22_biotrauma_fit.R has that for the on-ventilator divergence). The
# pre-intubation window holds only the patients who were in hospital and had the
# marker drawn twice before intubation, a sicker-on-the-ward subset; n_patients says
# how many. causal_testbed.R plants a known effect on the same patients and
# measurement times and reports which designs recover it at this site's sample.
#
# Writes to final/supplement/ (final/ungated/supplement/ with PBWPFVC_VTPBW_GATE=0):
#   causal_triangulation_estimates_<site>.csv   one row per marker x design
#   causal_triangulation_phases_<site>.csv      patients and measurements per arm x phase
#   causal_triangulation_<site>.pdf             the estimates, one panel per marker
#
# Runtime: seconds. Needs scripts 01-03 for the ventilated cohort and for
# PBWPFVC_COHORT=nosupport. PBWPFVC_DID_VTPBW=0 drops the ventilated arm's VT/PBW term
# (tables go to a novtpbw/ subfolder), as in the other ventilated-versus-control
# supplements.
#
# Usage: uvr run code/supplement/causal_triangulation.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(here)
})

rm(list = ls())
source("utils/config.R")
site_name <- config$site_name
final_dir <- final_dir_for("supplement")
source(here::here("code", "20_biotrauma_grid.R"))            # height_fingerprint(), sorted_term()
source(here::here("code", "supplement", "causal_common.R"))
if (!CAUSAL_VTPBW) { final_dir <- file.path(final_dir, "novtpbw"); dir.create(final_dir, showWarnings = FALSE) }

patients <- load_causal_patients(config)
designs  <- causal_designs()
message("Patients: ", sum(patients$vent == 1), " ventilated, ", sum(patients$vent == 0), " no support")

estimates <- list(); phases <- list()
for (m in seq_len(nrow(CAUSAL_MARKERS))) {
  marker <- CAUSAL_MARKERS$marker[m]
  panel <- load_causal_panel(config, patients, CAUSAL_MARKERS$lab_category[m])
  phases[[marker]] <- panel %>% inner_join(patients %>% select(id, arm), by = "id") %>%
    group_by(arm, phase) %>% summarise(n_patients = n_distinct(id), n_measurements = n(), .groups = "drop") %>%
    mutate(marker = marker, .before = 1)
  estimates[[marker]] <- fit_all_designs(designs, patients, panel) %>%
    mutate(marker = marker, injury_sign = CAUSAL_MARKERS$injury_sign[m]) %>%
    left_join(designs %>% select(design, formula, target, fails_if), by = "design")
}
estimates <- bind_rows(estimates) %>%
  mutate(site = site_name, cohort = if (config$vtpbw_gate == "0") "ungated" else "gated",
         conf_low = estimate - qnorm(0.975) * se, conf_high = estimate + qnorm(0.975) * se,
         p_value = 2 * pnorm(-abs(estimate / se)),
         injury_estimate = injury_sign * estimate,
         status = if_else(is.na(estimate), paste0("not estimable: fewer than ", MIN_PATIENTS, " patients in an arm with the phases"), "ok")) %>%
  # what a design pays for its assumption: its standard error over the unadjusted model's
  group_by(marker) %>% mutate(se_ratio = se / se[design == "unadjusted"]) %>% ungroup() %>%
  select(site, cohort, marker, design, estimate, se, conf_low, conf_high, p_value, injury_estimate,
         n_patients, n_obs, partial_r2, se_ratio, status, target, formula, fails_if)
phases <- bind_rows(phases) %>% mutate(site = site_name, .before = 1)

write_csv(estimates, file.path(final_dir, paste0("causal_triangulation_estimates_", site_name, ".csv")))
write_csv(phases, file.path(final_dir, paste0("causal_triangulation_phases_", site_name, ".csv")))

plot_data <- estimates %>% filter(status == "ok") %>%
  mutate(design = factor(design, levels = rev(designs$design)))
if (nrow(plot_data)) {
  p <- ggplot(plot_data, aes(estimate, design)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = conf_low, xmax = conf_high), size = 0.3) +
    facet_wrap(~ marker, scales = "free_x") +
    labs(x = "Ventilation-specific change in slope per log unit of PBW/PFVC (log marker per day), 95% CI",
         y = NULL, title = paste0("Complementary designs, ", site_name),
         subtitle = "Biotrauma predicts negative for platelets, positive for creatinine") +
    theme_bw(base_size = 10)
  ggsave(file.path(final_dir, paste0("causal_triangulation_", site_name, ".pdf")), p, width = 9, height = 5)
}

message("\nPatients with the phases each design needs:")
print(as.data.frame(phases), row.names = FALSE)
message("\nEstimates:")
print(as.data.frame(estimates %>% transmute(marker, design, estimate = signif(estimate, 3), se = signif(se, 3),
                                            p_value = signif(p_value, 2), n_patients, partial_r2 = round(partial_r2, 3))),
      row.names = FALSE)
message("\ncausal_triangulation complete: ", final_dir)
