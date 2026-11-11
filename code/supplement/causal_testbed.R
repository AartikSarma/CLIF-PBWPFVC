# =============================================================================
# Supplement (causal): a planted-truth test of the complementary designs
# =============================================================================
#
# causal_triangulation.R lays several designs side by side, each failing under a
# different assumption. This script asks what each of them returns when the truth is
# known. It keeps the site's own patients, cohorts, phases and measurement times, and
# replaces only the marker values with a simulated trajectory (plant_outcome() in
# causal_common.R) under each of these scenarios:
#
#   null             no strain effect, no demographic effect
#   strain           the hypothesis: the slope on the ventilator follows
#                    log VT/PFVC (= lambda + log VT/PBW); nothing else
#   age_everywhere   NO strain effect. Older patients' markers move faster in every
#                    phase and in both cohorts
#   age_sicker       NO strain effect. The same, but weaker in the controls
#                    (arm_scale): the age gradient scales with how sick a cohort is
#   age_acute        NO strain effect. The age gradient appears at t0 (intubation, or
#                    the control's ICU admission) in both cohorts: it follows the
#                    acute illness, not the ventilator
#   age_ventilation  NO strain effect. The age gradient exists only on the ventilator,
#                    by a path that is not strain (sedation, say)
#   sex_ventilation  NO strain effect. The same through sex
#   strain_plus_age  the strain effect and age_everywhere together
#
# Every confounded scenario is sized so that the unadjusted model sees exactly the
# gradient a strain effect of the planted size would give: the naive analysis cannot
# tell them from the hypothesis. A design passes a scenario when it rejects at the
# nominal rate under "NO strain effect" and recovers beta under strain. No design
# passes all of them; the table shows which combination of results each truth leaves
# behind, which is what the observed estimates are then read against. Two things to
# read off it before trusting a row of causal_triangulation: a rung of the ladder
# removes its own input's confounding and AMPLIFIES every other input's (the lambda it
# keeps is a smaller share of the same confounded demographics), so the rungs are read
# together, never one alone; and a design marked underpowered cannot speak against the
# hypothesis at this site.
#
# The planted strain effect is PBWPFVC_CAUSAL_EFFECT (default 0.15 log units per day
# per log unit of VT/PFVC, in the direction of injury: about 2% a day per SD of
# lambda). Noise comes from the site's observed marker (residual SD about each
# patient-phase line, between-patient SD of level and slope), so power and the
# minimum detectable effect are this site's, at its sample and its sampling density.
#
# Writes to final/supplement/ (final/ungated/supplement/ with PBWPFVC_VTPBW_GATE=0):
#   causal_testbed_scenarios_<site>.csv   the planted parameters of each scenario
#   causal_testbed_operating_<site>.csv   per marker x scenario x design: truth, mean
#                                         estimate, bias, empirical SD, mean SE,
#                                         coverage, rejection rate, MDE at 80% power
#   causal_testbed_<site>.pdf             rejection rate and bias, scenario x design
# No patient-level row is written.
#
# Runtime: about a minute per 100 replicates. PBWPFVC_CAUSAL_REPS sets the number
# (default 200). Needs scripts 01-03 for the ventilated cohort and for
# PBWPFVC_COHORT=nosupport.
#
# Usage: uvr run code/supplement/causal_testbed.R
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

N_REPS <- as.integer(Sys.getenv("PBWPFVC_CAUSAL_REPS", "200"))
EFFECT <- as.numeric(Sys.getenv("PBWPFVC_CAUSAL_EFFECT", "0.15"))
ARM_SCALE <- 0.3            # age_sicker: the control's share of the ventilated age gradient
TESTBED_SEED <- 20261111
ALPHA <- 0.05

patients <- load_causal_patients(config)
designs  <- causal_designs()

# the scenarios for a strain effect of size beta (signed: negative for platelets)
scenarios_for <- function(beta) {
  d_age <- confounding_equal_to(beta, patients, "age")
  d_sex <- confounding_equal_to(beta, patients, "sex")
  s <- function(scenario, beta = 0, index = "age", d_all = 0, d_acute = 0, d_vent = 0, arm_scale = 1)
    tibble(scenario, beta, index, d_all, d_acute, d_vent, arm_scale)
  bind_rows(
    s("null"),
    s("strain", beta = beta),
    s("age_everywhere", d_all = d_age),
    s("age_sicker", d_all = d_age, arm_scale = ARM_SCALE),
    s("age_acute", d_acute = d_age),
    s("age_ventilation", d_vent = d_age),
    s("sex_ventilation", index = "sex", d_vent = d_sex),
    s("strain_plus_age", beta = beta, d_all = d_age))
}

operating <- list(); scenario_log <- list()
for (m in seq_len(nrow(CAUSAL_MARKERS))) {
  marker <- CAUSAL_MARKERS$marker[m]
  panel  <- load_causal_panel(config, patients, CAUSAL_MARKERS$lab_category[m])
  noise  <- observed_noise(panel)
  scenarios <- scenarios_for(CAUSAL_MARKERS$injury_sign[m] * EFFECT)
  scenario_log[[marker]] <- scenarios %>%
    mutate(marker = marker, sd_resid = noise$sd_resid, sd_slope = noise$sd_slope, sd_level = noise$sd_level, .before = 1)
  preps <- lapply(seq_len(nrow(designs)), function(i) prepare_design(designs[i, ], patients, panel))
  rows  <- lapply(preps, function(p) if (p$ok) match(p$row_key, paste(panel$id, panel$phase, panel$t)) else integer(0))
  message(marker, ": ", sum(vapply(preps, function(p) p$ok, logical(1))), " of ", length(preps),
          " designs estimable; residual SD ", signif(noise$sd_resid, 3), ", slope SD ", signif(noise$sd_slope, 3))
  for (s in seq_len(nrow(scenarios))) {
    truth <- as.list(scenarios[s, ])
    fits <- vector("list", N_REPS)
    for (r in seq_len(N_REPS)) {
      sim <- plant_outcome(panel, patients, truth, noise, rng_seed = TESTBED_SEED + 1000L * s + r)
      fits[[r]] <- bind_rows(lapply(seq_along(preps), function(i) fit_design(preps[[i]], sim$y_sim[rows[[i]]])))
    }
    operating[[length(operating) + 1]] <- bind_rows(fits) %>% filter(!is.na(estimate)) %>%
      group_by(design) %>%
      summarise(n_reps = n(), mean_estimate = mean(estimate), empirical_sd = sd(estimate), mean_se = mean(se),
                coverage = mean(abs(estimate - truth$beta) <= qnorm(1 - ALPHA / 2) * se),
                rejection_rate = mean(abs(estimate / se) > qnorm(1 - ALPHA / 2)), .groups = "drop") %>%
      mutate(marker = marker, scenario = truth$scenario, truth = truth$beta, bias = mean_estimate - truth,
             # the effect a design would detect with 80% power at this site's sample
             mde_80 = (qnorm(1 - ALPHA / 2) + qnorm(0.80)) * empirical_sd)
    message("  ", truth$scenario, " done")
  }
}
operating <- bind_rows(operating) %>%
  left_join(designs %>% select(design, fails_if), by = "design") %>%
  mutate(site = site_name, cohort = if (config$vtpbw_gate == "0") "ungated" else "gated",
         planted_effect = EFFECT,
         # a design's verdict on a scenario. "biased": the mean estimate is off by more
         # than a quarter of the planted effect and by more than three Monte Carlo
         # standard errors. "underpowered": not biased, but this site's sample would
         # detect the planted effect less than 80% of the time, so a null from this
         # design is not evidence against strain.
         mc_se = empirical_sd / sqrt(n_reps),
         verdict = case_when(abs(bias) > 0.25 * EFFECT & abs(bias) > 3 * mc_se ~ "biased",
                             mde_80 > EFFECT ~ "underpowered",
                             TRUE ~ "ok")) %>%
  select(site, cohort, marker, scenario, design, truth, mean_estimate, bias, empirical_sd, mean_se,
         coverage, rejection_rate, mde_80, mc_se, verdict, n_reps, planted_effect, fails_if)
scenario_log <- bind_rows(scenario_log) %>% mutate(site = site_name, planted_effect = EFFECT, n_reps = N_REPS, .before = 1)

write_csv(operating, file.path(final_dir, paste0("causal_testbed_operating_", site_name, ".csv")))
write_csv(scenario_log, file.path(final_dir, paste0("causal_testbed_scenarios_", site_name, ".csv")))

plot_data <- operating %>%
  mutate(design = factor(design, levels = rev(designs$design)),
         scenario = factor(scenario, levels = unique(scenario_log$scenario)),
         label = sprintf("%.2f\n(%+.2f)", rejection_rate, bias / EFFECT))
p <- ggplot(plot_data, aes(scenario, design, fill = verdict)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = label), size = 2.3, lineheight = 0.85) +
  facet_wrap(~ marker) +
  scale_fill_manual(values = c(ok = "#cfe8d5", underpowered = "#f3e9c0", biased = "#f1c9c5"), name = NULL) +
  labs(x = NULL, y = NULL,
       title = paste0("Planted-truth test of the designs, ", site_name, " (", N_REPS, " replicates)"),
       subtitle = paste0("Rejection rate at 5%, and (bias as a share of the planted effect, ", EFFECT,
                         " log units/day per log unit). Only strain and strain_plus_age hold a true effect.")) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom", panel.grid = element_blank())
ggsave(file.path(final_dir, paste0("causal_testbed_", site_name, ".pdf")), p, width = 13, height = 7)

message("\nVerdicts (rows: design; columns: scenario), by marker:")
for (mk in unique(operating$marker)) {
  message("\n", mk, ": rejection rate (bias / planted effect)")
  tab <- operating %>% filter(marker == mk) %>%
    transmute(design = factor(design, levels = designs$design), scenario = factor(scenario, levels = unique(scenario_log$scenario)),
              cell = sprintf("%.2f (%+.2f)", rejection_rate, bias / EFFECT)) %>%
    arrange(design) %>% pivot_wider(names_from = scenario, values_from = cell)
  print(as.data.frame(tab), row.names = FALSE)
}
message("\ncausal_testbed complete: ", final_dir)
