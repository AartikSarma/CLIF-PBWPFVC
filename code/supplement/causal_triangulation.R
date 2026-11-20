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
# The estimate, in every row: the change in the outcome's slope (log units per day)
# per log unit of PBW/PFVC that belongs to ventilation. injury_estimate carries the
# sign so that positive is always more injury. partial_r2 is the share of the
# identifying term's variation that the design's other terms take away (0.99 in the
# adjusted model is the positivity problem made visible; in a difference design part of
# it is the comparison phase or arm, which is the design and not a loss), and se_ratio
# is the design's standard error over the unadjusted model's: what its assumption costs.
#
# Outcomes (CAUSAL_MARKERS in causal_common.R):
#   platelets, creatinine   log lab value, least squares
#   pressor_dose   the norepinephrine-equivalent dose with its zero bins kept: a model of
#                  the log of the EXPECTED dose (Poisson quasi-likelihood with a
#                  patient-phase intercept). One number for "how much pressor", with no
#                  added constant and no conditioning on being on one.
#   pressor_any    whether any pressor runs in a bin (linear probability; the estimate
#                  is a change in probability per day, not in log units)
#   osi            log oxygen saturation index. It needs a mean airway pressure, so it
#                  exists on the ventilator only: the control and before/after designs
#                  are not estimable for it, and the decision rule cannot be applied.
#                  Mean airway pressure also rises with tidal volume arithmetically, so
#                  read its level beside its components (21_biotrauma_panel.R).
#
# Every design is fitted twice (column vtpbw): with the index VT/PBW "held" (the
# hypothesis as stated, at a fixed VT/PBW) and "dropped". Where tidal volume is chosen
# by how sick the lung is, holding VT/PBW lets severity into the lambda gradient (the
# testbed's vt_by_indication scenario); dropping it avoids that and estimates the
# gradient as practice delivers it. The decision rule is applied both ways and the
# data support the hypothesis only if it holds both ways. With VT/PBW held the table
# also carries vt_gradient, the slope's gradient on log VT/PBW itself: under strain,
# with dose = log VT - log lung size, it equals the lambda gradient (p_vt_equal).
#
# With VT/PBW dropped, lambda alone also assumes the index tidal volume was set at the
# same multiple of PBW in everyone (rho = 1): only then do log PBW and log PFVC act
# equal and opposite. A reported check, not part of the decision rule, frees log PBW
# beside lambda in every dropped design (fit_pbw_separate() in causal_common.R):
#   pbw_separate   estimate is the gradient along PFVC at a fixed PBW, whatever rho is;
#                  ratio_estimate is the lambda-only fit beside it; pbw_extra is the
#                  coefficient on log PBW, b (rho - 1) under strain (so close to
#                  estimate x (rho - 1), given as pbw_extra_strain), p_ratio_adequate
#                  tests it against 0 and p_strain_consistent against that prediction.
#                  Where rho is not 1 read p_strain_consistent: a non-zero pbw_extra
#                  is then what strain predicts. It
#                  also takes any direct effect of height or sex, so read it beside
#                  rho, not alone. Not simulated: unit tests only. At a fixed PBW lambda moves
#                  mainly with PFVC's age and race terms, so these are less precise and
#                  lean on age. In the free_height and free_sex rows log PBW is nearly
#                  a function of the freed input, and its coefficient may be lost.
#   vt_pbw         rho itself: the slope of log VT on log PBW in ventilated patients.
#                  The default cohort is gated on VT/PBW 6-8 mL/kg, which pulls rho
#                  toward 1; PBWPFVC_VTPBW_GATE=0 shows how tidal volume was set.
#
# Aeration (fit_aeration(), compliance_scaling() in causal_common.R). Compliance over
# predicted lung size estimates the fraction of the lung still aerated, and the
# driving pressure splits into log VT/PBW + lambda - log(Crs/PFVC). The held designs are
# refitted with log(Crs/PFVC) and log BMI: strain gives the three terms equal weight
# (p_strain_restriction), and the VT/PBW gradient is then adjusted for tidal volume
# lowered for stiff lungs. Specific elastance is similar on average across patients,
# not constant (about two-fold from low to high), so the adjustment is incomplete; Crs
# includes the chest wall; and it exists only where a plateau pressure was charted
# (base_* columns: the same patients without the terms). compliance_scaling asks which
# predicted volume compliance follows on average.
#
# Two questions the main table does not answer have their own tables:
#   denominators   Is PFVC the volume that sets strain? The designs again with lambda =
#                  log(PBW / volume) for PFVC at age 25, TLC, alveolar volume and FRC
#                  (script 03), all on the patients in whom every volume is defined
#                  (age <= 80). The volumes differ mostly in how they change with age
#                  (FVC falls by a third from 40 to 80, TLC hardly, FRC rises), so
#                  this comparison and the age question are the same question. The
#                  "encompassing" rows put a candidate's lambda and what PFVC adds to
#                  it (ldiff) in one model: under the candidate the coefficient on
#                  ldiff is 0 (p_candidate_enough), under PFVC it equals the
#                  candidate's (p_pfvc_enough).
#   timecourse     Does the effect accrue steadily, early, or late? The on-ventilator
#                  phase cut at PBWPFVC_CAUSAL_KNOTS (days 2 and 4): the gradient within
#                  each segment and the cumulative effect by the end of each
#                  (fit_time_course() in causal_common.R). No shape is assumed, so fast
#                  and slow markers can be read on the same scale.
#
# What the script does not do: it is the longitudinal model alone, conditional on a
# patient still being under observation, with no joint model for death, liberation or
# escalation (22_biotrauma_fit.R has that for the on-ventilator divergence). The
# pre-intubation window holds only the patients who were in hospital before
# intubation, a sicker-on-the-ward subset; n_patients says how many. causal_testbed.R
# plants a known effect and scores the designs and the rule for platelets and
# creatinine; the vasopressor and OSI models have unit tests but no simulation yet.
#
# Writes to final/supplement/ (final/ungated/supplement/ with PBWPFVC_VTPBW_GATE=0):
#   causal_triangulation_estimates_<site>.csv      marker x vtpbw x design
#   causal_triangulation_rule_<site>.csv           the decision rule per marker
#   causal_triangulation_denominators_<site>.csv   marker x vtpbw x denominator x design
#   causal_triangulation_timecourse_<site>.csv     marker x vtpbw x design x window
#   causal_triangulation_phases_<site>.csv         patients and measurements per arm x phase
#   causal_triangulation_pbw_separate_<site>.csv   marker x design: VT/PBW dropped, PBW freed
#   causal_triangulation_vt_pbw_<site>.csv         rho, overall and by sex
#   causal_triangulation_aeration_<site>.csv       the held designs with log(Crs/PFVC)
#                                                  and log BMI (fit_aeration())
#   causal_triangulation_compliance_scaling_<site>.csv  log Crs on log PBW and each
#                                                  predicted volume (compliance_scaling())
#   causal_triangulation_age_mediation_<site>.csv  lambda split into its age part and
#                                                  the rest (fit_age_mediation())
#   causal_triangulation_<site>.pdf
#
# Runtime: about a minute. Needs scripts 01-03 for the ventilated cohort and for
# PBWPFVC_COHORT=nosupport. PBWPFVC_CAUSAL_MARKERS picks outcomes (default: all).
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

MARKERS <- strsplit(Sys.getenv("PBWPFVC_CAUSAL_MARKERS", paste(CAUSAL_MARKERS$marker, collapse = ",")), ",")[[1]]
stopifnot(all(MARKERS %in% CAUSAL_MARKERS$marker))
patients <- load_causal_patients(config)
VTPBW <- c(held = TRUE, dropped = FALSE)
designs <- lapply(VTPBW, causal_designs)
TIMECOURSE_DESIGNS <- c("unadjusted", "control_did", "free_age_did", "within_onoff", "off_reversal")
vt_pbw_fit <- vt_follows_pbw(patients)
RHO <- vt_pbw_fit$rho[vt_pbw_fit$group == "all"]
message("Patients: ", sum(patients$vent == 1), " ventilated, ", sum(patients$vent == 0), " no support")

# a candidate's lambda and what PFVC adds to it, in one model
encompassing_designs <- function(vtpbw) {
  dose <- if (vtpbw) " + lv_on" else ""
  tibble(design = c("encompassing", "encompassing_did"),
         arms = list("ventilated", c("ventilated", "no support")), phases = list("on", "on"),
         formula = c(paste0("~ ldisc + ldiff + sofa_z + sf_z", dose),
                     paste0("~ vent * (sofa_z + sf_z) + ldisc + ldiff + vent:ldisc + vent:ldiff", dose)),
         target = c("ldisc", "vent:ldisc"), extra = c("ldiff", "vent:ldiff"), sign = 1, fails_if = "")
}
fit_encompassing <- function(design, pats, panel, family) {
  prep <- prepare_design(design, pats, panel, raw = family == "poisson")
  fit <- if (prep$ok) fit_coefs(prep, align_outcome(prep, panel), family) else NULL
  out <- tibble(design = design$design, estimate = NA_real_, se = NA_real_, pfvc_extra = NA_real_, pfvc_extra_se = NA_real_,
                p_candidate_enough = NA_real_, p_pfvc_enough = NA_real_, n_patients = prep$n_patients, n_obs = prep$n_obs)
  e <- which(names(fit$beta) == design$extra)
  if (is.null(fit) || length(e) != 1) return(out)
  k <- prep$k; b <- fit$beta; V <- fit$V
  out$estimate <- b[[k]]; out$se <- sqrt(V[k, k]); out$pfvc_extra <- b[[e]]; out$pfvc_extra_se <- sqrt(V[e, e])
  out$p_candidate_enough <- 2 * pnorm(-abs(b[[e]]) / sqrt(V[e, e]))
  out$p_pfvc_enough <- 2 * pnorm(-abs(b[[e]] - b[[k]]) / sqrt(V[e, e] + V[k, k] - 2 * V[e, k]))
  out
}

estimates <- list(); phases <- list(); rule <- list(); denominators <- list(); timecourse <- list(); pbw_separate <- list()
age_mediation <- list(); aeration <- list()
for (marker in MARKERS) {
  m <- CAUSAL_MARKERS[CAUSAL_MARKERS$marker == marker, ]
  panel <- load_causal_outcome(config, patients, marker)
  phases[[marker]] <- panel %>% inner_join(patients %>% select(id, arm), by = "id") %>%
    group_by(arm, phase) %>% summarise(n_patients = n_distinct(id), n_measurements = n(), mean_y = mean(y), .groups = "drop") %>%
    mutate(marker = marker, .before = 1)
  fits <- lapply(designs, fit_all_designs, patients = patients, panel = panel, family = m$family)
  estimates[[marker]] <- bind_rows(fits, .id = "vtpbw") %>%
    left_join(designs$held %>% select(design, target, fails_if), by = "design") %>%
    mutate(marker = marker, family = m$family, injury_sign = m$injury_sign)
  both <- triangulation_supports_both(fits$held, fits$dropped, m$injury_sign)
  # VT/PBW dropped, with log PBW freed beside lambda (a reported check, not in the rule)
  pbw_separate[[marker]] <- fit_pbw_separate(designs$dropped, patients, panel, m$family, rho = RHO) %>%
    left_join(fits$dropped %>% select(design, ratio_estimate = estimate, ratio_se = se), by = "design") %>%
    mutate(marker = marker, injury_sign = m$injury_sign, .before = 1)
  age_mediation[[marker]] <- bind_rows(lapply(designs, fit_age_mediation, patients = patients, panel = panel, family = m$family),
                                       .id = "vtpbw") %>% mutate(marker = marker, .before = 1)
  aeration[[marker]] <- fit_aeration(designs$held, patients, panel, m$family) %>% mutate(marker = marker, .before = 1)
  rule[[marker]] <- tibble(marker = marker, supports = unname(both[["rule_both"]]), supports_vtpbw_held = unname(both[["rule_held"]]),
                           supports_vtpbw_dropped = unname(both[["rule_dropped"]]),
                           can_be_applied = !anyNA(fits$held$estimate[fits$held$design %in% c("control_did", paste0("free_", c("age", "sex", "race", "height"), "_did"))]))
  for (v in names(VTPBW)) {
    for (den in names(CAUSAL_DENOMINATORS)) {
      pats <- use_denominator(patients, den)
      pan <- panel[panel$id %in% pats$id, , drop = FALSE]
      rows <- fit_all_designs(designs[[v]] %>% filter(design != "fingerprint"), pats, pan, family = m$family)
      if (den != "pfvc") {
        enc <- encompassing_designs(VTPBW[[v]])
        rows <- bind_rows(rows, map_dfr(seq_len(nrow(enc)), ~ fit_encompassing(enc[.x, ], pats, pan, m$family)))
      }
      denominators[[paste(marker, v, den)]] <- rows %>% mutate(marker = marker, vtpbw = v, denominator = den, .before = 1)
    }
    timecourse[[paste(marker, v)]] <- map_dfr(TIMECOURSE_DESIGNS, function(dn)
      fit_time_course(designs[[v]] %>% filter(design == dn), patients, panel, m$family)) %>%
      mutate(marker = marker, vtpbw = v, injury_sign = m$injury_sign, .before = 1)
  }
  message("  ", marker, " done")
}
cohort <- if (config$vtpbw_gate == "0") "ungated" else "gated"
with_interval <- function(d) d %>%
  mutate(site = site_name, cohort = cohort, .before = 1) %>%
  mutate(conf_low = estimate - qnorm(0.975) * se, conf_high = estimate + qnorm(0.975) * se,
         p_value = 2 * pnorm(-abs(estimate / se)))
NOT_ESTIMABLE <- paste0("not estimable: fewer than ", MIN_PATIENTS, " patients in an arm with the phases, or the model could not be fitted")
estimates <- bind_rows(estimates) %>% with_interval() %>%
  mutate(injury_estimate = injury_sign * estimate, status = if_else(is.na(estimate), NOT_ESTIMABLE, "ok")) %>%
  # what a design pays for its assumption: its standard error over the unadjusted model's
  group_by(marker, vtpbw) %>% mutate(se_ratio = se / se[design == "unadjusted"]) %>% ungroup() %>%
  select(site, cohort, marker, family, vtpbw, design, estimate, se, conf_low, conf_high, p_value, injury_estimate,
         vt_gradient, vt_gradient_se, p_vt_equal, n_patients, n_obs, partial_r2, se_ratio, status, target, fails_if)
rule <- bind_rows(rule) %>% mutate(site = site_name, cohort = cohort, .before = 1) %>%
  mutate(across(starts_with("supports"), ~ if_else(can_be_applied, .x, NA)))
denominators <- bind_rows(denominators) %>% with_interval() %>%
  select(site, cohort, marker, vtpbw, denominator, design, estimate, se, conf_low, conf_high, p_value, resid_sd,
         pfvc_extra, pfvc_extra_se, p_candidate_enough, p_pfvc_enough, n_patients, n_obs)
timecourse <- bind_rows(timecourse) %>% with_interval() %>%
  select(site, cohort, marker, vtpbw, design, quantity, window, from, to, estimate, se, conf_low, conf_high, p_value, n_patients, n_obs)
phases <- bind_rows(phases) %>% mutate(site = site_name, .before = 1)
pbw_separate <- bind_rows(pbw_separate) %>% with_interval() %>%
  select(site, cohort, marker, design, estimate, se, conf_low, conf_high, p_value, pbw_extra, pbw_extra_se,
         p_ratio_adequate, rho, pbw_extra_strain, p_strain_consistent, ratio_estimate, ratio_se, n_patients, n_obs, injury_sign)
age_mediation <- bind_rows(age_mediation) %>% mutate(site = site_name, cohort = cohort, .before = 1)
aeration <- bind_rows(aeration) %>% mutate(site = site_name, cohort = cohort, .before = 1)
scaling <- compliance_scaling(patients) %>% mutate(site = site_name, cohort = cohort, .before = 1)
vt_pbw <- vt_pbw_fit %>% mutate(site = site_name, cohort = cohort, .before = 1)

out_file <- function(name, ext = "csv") file.path(final_dir, paste0("causal_triangulation_", name, site_name, ".", ext))
write_csv(estimates, out_file("estimates_"))
write_csv(rule, out_file("rule_"))
write_csv(denominators, out_file("denominators_"))
write_csv(timecourse, out_file("timecourse_"))
write_csv(phases, out_file("phases_"))
write_csv(pbw_separate, out_file("pbw_separate_"))
write_csv(vt_pbw, out_file("vt_pbw_"))
write_csv(age_mediation, out_file("age_mediation_"))
write_csv(aeration, out_file("aeration_"))
write_csv(scaling, out_file("compliance_scaling_"))

pdf(out_file("", "pdf"), width = 11, height = 7.5)
plot_data <- estimates %>% filter(status == "ok") %>% mutate(design = factor(design, levels = rev(designs$held$design)))
if (nrow(plot_data))
  print(ggplot(plot_data, aes(estimate, design, colour = vtpbw)) +
          geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
          geom_pointrange(aes(xmin = conf_low, xmax = conf_high), size = 0.25, position = position_dodge(width = 0.6)) +
          facet_wrap(~ marker, scales = "free_x", nrow = 1) + scale_colour_brewer(palette = "Dark2", name = "Index VT/PBW") +
          labs(x = "Ventilation-specific change in slope per log unit of PBW/PFVC, 95% CI", y = NULL,
               title = paste0("Complementary designs, ", site_name),
               subtitle = "Injury is negative for platelets and positive for the others. pressor_any is a change in probability per day.") +
          theme_bw(base_size = 9) + theme(legend.position = "bottom"))
plot_data <- denominators %>% filter(!is.na(estimate), design %in% c("unadjusted", "control_did", "free_age_did", "within_onoff"), vtpbw == "held")
if (nrow(plot_data))
  print(ggplot(plot_data, aes(estimate, denominator)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
          geom_pointrange(aes(xmin = conf_low, xmax = conf_high), size = 0.25) + facet_grid(design ~ marker, scales = "free_x") +
          labs(x = "Gradient per log unit of PBW / predicted volume, 95% CI (VT/PBW held)", y = NULL,
               title = "The strain denominator", subtitle = "Patients in whom every volume is defined (age <= 80)") + theme_bw(base_size = 9))
plot_data <- timecourse %>% filter(!is.na(estimate), vtpbw == "held")
if (nrow(plot_data))
  print(ggplot(plot_data, aes(to, estimate, colour = quantity)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
          geom_pointrange(aes(ymin = conf_low, ymax = conf_high), size = 0.25, position = position_dodge(width = 0.4)) +
          geom_line(position = position_dodge(width = 0.4)) + facet_grid(marker ~ design, scales = "free_y") +
          scale_colour_brewer(palette = "Dark2", name = NULL) +
          labs(x = "End of segment (days from intubation)", y = "Segment: gradient per day. Cumulative: gradient of the total change",
               title = "The time course (VT/PBW held)") + theme_bw(base_size = 9) + theme(legend.position = "bottom"))
invisible(dev.off())

message("\nPatients with the phases each design needs:")
print(as.data.frame(phases), row.names = FALSE)
message("\nEstimates:")
print(as.data.frame(estimates %>% transmute(marker, vtpbw, design, estimate = signif(estimate, 3), se = signif(se, 3),
                                            p_value = signif(p_value, 2), n_patients, partial_r2 = round(partial_r2, 3))),
      row.names = FALSE)
message("\nThe gradient on log VT/PBW itself, VT/PBW held (under strain it equals the lambda gradient):")
print(as.data.frame(estimates %>% filter(vtpbw == "held", !is.na(vt_gradient)) %>%
                      transmute(marker, design, estimate = signif(estimate, 3), vt_gradient = signif(vt_gradient, 3),
                                vt_gradient_se = signif(vt_gradient_se, 3), p_vt_equal = signif(p_vt_equal, 2))), row.names = FALSE)
message("\nOther strain denominators (patients aged 80 or under in whom every volume is defined):")
print(as.data.frame(denominators %>% filter(design %in% c("unadjusted", "adjusted", "control_did", "free_age_did", "free_sex_did",
                                                          "free_race_did", "free_height_did", "within_onoff", "off_reversal",
                                                          "encompassing", "encompassing_did")) %>%
                      transmute(marker, vtpbw, denominator, design, estimate = signif(estimate, 3), se = signif(se, 3),
                                p_value = signif(p_value, 2), pfvc_extra = signif(pfvc_extra, 3),
                                p_candidate_enough = signif(p_candidate_enough, 2), n_patients)), row.names = FALSE)
message("\nVT/PBW dropped, PBW and PFVC as separate terms (estimate: at a fixed PBW; ratio_estimate: lambda alone):")
print(as.data.frame(pbw_separate %>% transmute(marker, design, estimate = signif(estimate, 3), se = signif(se, 3),
                                               ratio_estimate = signif(ratio_estimate, 3), pbw_extra = signif(pbw_extra, 3),
                                               p_ratio_adequate = signif(p_ratio_adequate, 2),
                                               pbw_extra_strain = signif(pbw_extra_strain, 3),
                                               p_strain_consistent = signif(p_strain_consistent, 2), n_patients)), row.names = FALSE)
message("\nDoes PFVC carry the age gradient? b_rest: lambda's non-age part; b_age: its age part (share_mediated = b_rest / b_age):")
print(as.data.frame(age_mediation %>% filter(design %in% c("unadjusted", "control_did", "within_onoff", "off_reversal")) %>%
                      transmute(marker, vtpbw, design, b_rest = signif(b_rest, 3), b_rest_se = signif(b_rest_se, 3),
                                b_age = signif(b_age, 3), b_age_se = signif(b_age_se, 3), p_full_mediation = signif(p_full_mediation, 2),
                                share_mediated = signif(share_mediated, 2), share_se = signif(share_mediated_se, 2),
                                age_per_decade = signif(age_gradient_per_decade, 2), via_pfvc_per_decade = signif(through_pfvc_per_decade, 2),
                                n_patients)), row.names = FALSE)
message("\nAeration, VT/PBW held (log(Crs/PFVC) and log BMI added; base_*: the same patients without them):")
print(as.data.frame(aeration %>% transmute(marker, design, estimate = signif(estimate, 3), se = signif(se, 3),
                                           vt_gradient = signif(vt_gradient, 3), vt_se = signif(vt_gradient_se, 3),
                                           aeration = signif(aeration_gradient, 3), aeration_se = signif(aeration_gradient_se, 3),
                                           p_strain_restriction = signif(p_strain_restriction, 2),
                                           base_estimate = signif(base_estimate, 3), base_vt_gradient = signif(base_vt_gradient, 3),
                                           n_patients, n_with_compliance)), row.names = FALSE)
message("\nWhat compliance scales with: log Crs on log PBW, the log predicted volume and log BMI:")
print(as.data.frame(scaling %>% transmute(subset, volume, n_patients, volume_coef = signif(volume_coef, 3), volume_se = signif(volume_se, 2),
                                          pbw_coef = signif(pbw_coef, 3), pbw_se = signif(pbw_se, 2), bmi_coef = signif(bmi_coef, 3),
                                          residual_sd = signif(residual_sd, 3), aic = round(aic, 1))), row.names = FALSE)
message("\nHow closely the index tidal volume followed PBW (rho = 1: the same multiple of PBW in everyone):")
print(as.data.frame(vt_pbw %>% select(group, n_patients, rho, se, conf_low, conf_high, p_rho_one) %>%
                      mutate(across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)
message("\nDecision rule:")
print(as.data.frame(rule), row.names = FALSE)
message("\ncausal_triangulation complete: ", final_dir)
