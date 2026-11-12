# =============================================================================
# Supplement (causal): death as the far end of the injury, not as a nuisance
# An illness-death model of organ injury and death along PBW/PFVC
# =============================================================================
#
# causal_triangulation.R reads a marker's slope while the patient is alive and measured.
# Death ends that record, and where strain kills it does so partly through the same
# injury the markers measure: the patients whose markers worsen fastest leave first.
# Treating death as censoring pulls the slope toward zero; treating it as an unrelated
# competing event misstates what it is. This script gives death its own place. Each
# patient moves between three states,
#
#     not injured  --(1)-->  injured  --(3)-->  dead
#          \------------------(2)---------------^
#
# and the lambda = log(PBW/PFVC) gradient is estimated on each transition's hazard:
#   injury                 (1) the first time the marker crosses its threshold on the
#                              ventilator (the control: before escalation), within
#                              PBWPFVC_CAUSAL_ON_DAYS of t0
#   death without injury   (2) death before any injury, to PBWPFVC_CAUSAL_HORIZON days
#   death after injury     (3) death once injured, on the same clock (days from t0),
#                              entering at the time of injury
#   death, any             for reference: death whatever came before
#
# What the hypothesis predicts, and what the rival does. If strain injures and injury
# kills, lambda raises the hazard of injury, and raises death mainly by way of the
# injured state: the gradient on death WITHOUT injury is smaller than the gradient on
# injury ("contrast" rows: injury minus death without injury, with a standard error
# that allows for the two coming from the same patients). If lambda marks frailty,
# older patients die more by every route, injured or not, and the gradient on death
# without injury is as large as any. No sharp prediction is made for (3): strain may
# or may not matter once a patient is injured.
#
# Injury, per marker (the first threshold of each is the default; settings below):
#   creatinine   KDIGO stage 2 or worse: creatinine at least 2x baseline, or renal
#                replacement started. PBWPFVC_CAUSAL_AKI_STAGE=1 uses stage 1 (1.5x
#                baseline, or a rise of 0.3 mg/dL within 48 hours).
#   platelets    below 100, or a fall of at least half from baseline
#                (PBWPFVC_CAUSAL_PLT_BELOW, PBWPFVC_CAUSAL_PLT_FALL)
#   pressor      any vasopressor after t0; PBWPFVC_CAUSAL_PRESSOR_MIN=0.1 asks for a
#                norepinephrine-equivalent rate of at least 0.1 mcg/kg/min
# Baseline is the best value from hospital admission to PBWPFVC_CAUSAL_GRACE_H hours
# after t0 (lowest creatinine, highest platelet count), so that the sample drawn at
# intubation counts as baseline and not as an outcome. Only NEW injury is counted: a
# patient already past the threshold at the end of that window (for a pressor: one
# running in the 24 hours before t0) is left out of that marker's analysis, as is a
# patient with no value to take a baseline from. Once injured a patient stays injured.
# The oxygen saturation index is not used: it exists only on the ventilator.
#
# Designs. The single-phase designs of causal_common.R, as Cox models stratified by
# cohort with the same covariates: unadjusted, adjusted, fingerprint, control_did (the
# ventilated minus the control gradient on the same transition) and the ladder
# (free_<input>, free_<input>_did). The before-and-after designs have no counterpart: a
# patient is injured, and dies, once. Each is fitted with the index VT/PBW held and
# dropped (column vtpbw), as in causal_triangulation.R.
#
# The estimate is a log hazard ratio per log unit of PBW/PFVC; hr_per_0.1 is the hazard
# ratio for 0.1 log units (about 10% in the ratio). Positive is more injury or death.
# A transition with fewer than MIN_EVENTS events in an arm is not estimated.
#
# What this does not do. A control's follow-up ends at escalation of support, so the
# control arm's deaths are those before any escalation. Injury can only be seen when
# the marker is measured, so a patient measured more often is seen to be injured
# sooner. Death after discharge is counted where the site records it. And a threshold
# discards the continuous information the slope designs use: this is the complement
# of causal_triangulation.R, not a replacement.
#
# Writes to final/supplement/ (final/ungated/supplement/ with PBWPFVC_VTPBW_GATE=0):
#   causal_illness_death_estimates_<site>.csv   marker x vtpbw x design x transition
#   causal_illness_death_counts_<site>.csv      patients and events per marker x arm
#   causal_illness_death_<site>.pdf
#
# Runtime: about a minute. Needs scripts 01-03 for the ventilated cohort and for
# PBWPFVC_COHORT=nosupport. PBWPFVC_CAUSAL_MARKERS picks markers (default: all).
#
# Usage: uvr run code/supplement/causal_illness_death.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(here)
  library(survival)
})

rm(list = ls())
source("utils/config.R")
site_name <- config$site_name
final_dir <- final_dir_for("supplement")
source(here::here("code", "20_biotrauma_grid.R"))            # height_fingerprint(), pfvc_channels(), sorted_term()
source(here::here("code", "supplement", "causal_common.R"))  # also the illness-death functions

patients <- load_causal_patients(config)
message("Patients: ", sum(patients$vent == 1), " ventilated, ", sum(patients$vent == 0), " no support")
VTPBW <- c(held = TRUE, dropped = FALSE)

estimates <- list(); counts <- list()
MARKERS <- strsplit(Sys.getenv("PBWPFVC_CAUSAL_MARKERS", paste(names(INJURY_MARKERS), collapse = ",")), ",")[[1]]
MARKERS <- intersect(names(INJURY_MARKERS), MARKERS)
if (!length(MARKERS)) stop("PBWPFVC_CAUSAL_MARKERS names none of: ", paste(names(INJURY_MARKERS), collapse = ", "))
for (marker in MARKERS) {
  injury <- load_injury_times(config, patients, marker)
  states <- illness_death_states(patients, injury)
  definition <- INJURY_MARKERS[[marker]]
  counts[[marker]] <- illness_death_counts(states) %>% mutate(marker = marker, definition = definition, .before = 1)
  for (v in names(VTPBW))
    estimates[[paste(marker, v)]] <- fit_illness_death(causal_designs(VTPBW[[v]]), states) %>%
      mutate(marker = marker, vtpbw = v, .before = 1)
  message("  ", marker, " done")
}
cohort <- if (config$vtpbw_gate == "0") "ungated" else "gated"
estimates <- bind_rows(estimates) %>%
  mutate(site = site_name, cohort = cohort, .before = 1) %>%
  mutate(conf_low = estimate - qnorm(0.975) * se, conf_high = estimate + qnorm(0.975) * se,
         p_value = 2 * pnorm(-abs(estimate / se)),
         hr_per_0.1 = if_else(transition == "contrast: injury minus death without injury", NA_real_, exp(0.1 * estimate)),
         status = if_else(is.na(estimate), paste0("not estimable: fewer than ", MIN_EVENTS, " events in an arm, or the model could not be fitted"), "ok")) %>%
  select(site, cohort, marker, vtpbw, design, transition, estimate, se, conf_low, conf_high, p_value, hr_per_0.1,
         n_patients, n_events, status)
counts <- bind_rows(counts) %>% mutate(site = site_name, cohort = cohort, .before = 1)

write_csv(estimates, file.path(final_dir, paste0("causal_illness_death_estimates_", site_name, ".csv")))
write_csv(counts, file.path(final_dir, paste0("causal_illness_death_counts_", site_name, ".csv")))

plot_data <- estimates %>% filter(status == "ok", vtpbw == "held", design %in% c("unadjusted", "control_did", "free_age_did", "adjusted")) %>%
  mutate(transition = factor(transition, levels = rev(c(ILLNESS_DEATH_TRANSITIONS, "contrast: injury minus death without injury"))))
if (nrow(plot_data)) {
  p <- ggplot(plot_data, aes(estimate, transition)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = conf_low, xmax = conf_high), size = 0.25) + facet_grid(design ~ marker, scales = "free_x") +
    labs(x = "Log hazard ratio per log unit of PBW/PFVC, 95% CI (VT/PBW held)", y = NULL,
         title = paste0("Injury and death along PBW/PFVC, ", site_name),
         subtitle = "Strain predicts a gradient on injury larger than on death without injury; frailty predicts one on every route") +
    theme_bw(base_size = 9)
  ggsave(file.path(final_dir, paste0("causal_illness_death_", site_name, ".pdf")), p, width = 11, height = 7)
}

message("\nPatients and events:")
print(as.data.frame(counts %>% select(-site, -cohort, -definition)), row.names = FALSE)
message("\nEstimates (VT/PBW held):")
print(as.data.frame(estimates %>% filter(vtpbw == "held", design %in% c("unadjusted", "control_did", "free_age_did")) %>%
                      transmute(marker, design, transition, estimate = signif(estimate, 3), se = signif(se, 3),
                                p_value = signif(p_value, 2), n_events)), row.names = FALSE)
message("\ncausal_illness_death complete: ", final_dir)
