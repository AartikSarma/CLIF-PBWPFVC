# =============================================================================
# Supplement (cross-sectional): does the VT/PBW 6-8 gate make predicted lung size
# look sicker than it is?
# =============================================================================
# Clinicians set VT/PBW, and two things push it. Illness pushes it down: SOFA is
# highest at 6 mL/kg and falls to about 7.3. Predicted lung size pushes it up: the
# PBW/PFVC ratio (high where PBW oversizes the lung: short, older, female patients)
# rises across the band, most likely because round-number volumes land at a high
# mL/kg on a small PBW. VT/PBW is then a common effect of the ratio and of illness,
# and selecting on it or adjusting for it links the two. At a given VT/PBW, a
# high-ratio patient who is not at the high VT/PBW that dosing habit gives them was
# more often turned down for illness, so conditioning makes a high ratio (a low
# PFVC) look sicker. That is the direction of the paper's finding, and the paper's
# cohort does both: it enters only patients with a VT/PBW of 6-8 at some complete,
# hypoxemic timepoint, and several of its models adjust for VT/PBW.
#
# The check measures the association that the gate and the adjustment create, with
# no outcome model. The reference is the ungated index (script 03, 3e3): the paper's
# two-tier index rule applied to every complete-data hypoxemic IMV timepoint,
# whatever its VT/PBW, with SOFA scored over the 24 h from that index. SOFA is
# regressed on the size term under six designs:
#
#   none                   ungated index, no VT/PBW term: the reference
#   conditioned            ungated index, + ns(vtpbw, 4)
#   index 6-8              ungated index, patients whose VT/PBW there is 6-8
#   index 6-8 + VT/PBW     the same, + vtpbw (linear, as 04 enters it)
#   paper cohort           the paper's own rows: index at a 6-8 timepoint, SOFA
#                          from that index (analysis_cross_sectional)
#   paper cohort + VT/PBW  the same, + vtpbw
#
# The gate acts on the timepoint more than on the patient. Almost every ventilated
# patient has some complete hypoxemic timepoint at 6-8 mL/kg and so enters the
# paper's cohort (the _gate_ table counts those who do not); what the gate decides is
# WHICH timepoint becomes the index. The "index 6-8" designs isolate the selection on
# VT/PBW at one fixed index. The "paper cohort" designs are what the paper fitted;
# they differ from the reference both by that selection and by indexing some
# patients later, so they bound the gate's total effect rather than isolate it.
#
# If the gate and the adjustment are harmless on this path, every design gives the
# reference's slope. A slope that moves away from the reference, positive for the
# ratio or negative for log PFVC, is the bias the gate adds. Every design's change
# from the reference has a patient bootstrap interval: the designs share patients,
# so their standard errors cannot be differenced.
#
# SOFA is split into its respiratory component and the rest. If the respiratory
# component carries the association, the selection is titration to oxygenation or
# mechanics, which predicted lung size itself moves (a small lung has a higher
# driving pressure at a given VT/PBW), and adjusting for severity will not remove
# it. If the non-respiratory rest carries it, the selection is general illness,
# which the models' SOFA term absorbs.
#
# Exposures: the PBW/PFVC ratio and log PFVC (04's primary size term), each per SD of
# the ungated sample, so all designs share one scale (the SD is written beside the
# estimate). Every model is fitted unadjusted and adjusted for ns(age10, 4), sex and
# race. Adjusted, log PFVC is identified largely by height; the adjusted ratio keeps
# little variance once age is a spline (supplement/xsec_age_form_check.R), so read
# it for direction only.
#
# Practice. VT/PBW habits differ widely between providers and between sites: the
# share of set volumes at or below 6 mL/kg runs from almost none to most. Practice
# that follows neither the ratio nor illness is noise in VT/PBW's collider role and
# weakens the bias within a site. The practice mix therefore sets how much the gate
# bends the association at each site, so read this check site by site before pooling
# it. The two arms of the collider (VT/PBW on the size term and on SOFA) and the
# share of VT/PBW they leave unexplained are reported for that reason.
#
# Writes to final/supplement/:
#   vtpbw_gate_collider_estimates_{site}.csv  the size term's SOFA slope under each
#                                             design, with its bootstrap change from
#                                             the reference
#   vtpbw_gate_collider_arms_{site}.csv       VT/PBW on the size term and SOFA
#                                             (both arms, and R-squared)
#   vtpbw_gate_collider_bins_{site}.csv       at the ungated index, by 0.5 mL/kg
#                                             VT/PBW bin [lower, upper): counts, how
#                                             many are in the paper's cohort, and the
#                                             mean and SD of the ratio, log PFVC and
#                                             SOFA
#   vtpbw_gate_collider_gate_{site}.csv       the gate's counts: VT/PBW below, in and
#                                             above the band at the ungated index, by
#                                             membership of the paper's cohort
#   vtpbw_gate_collider_{site}.pdf            the slopes by design
#
# Needs scripts 01-03 (the ventilated cohort; 03 writes analysis_ungated_index).
# Usage: uvr run code/supplement/xsec_vtpbw_gate_collider.R
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(splines) })
rm(list = ls())
source("utils/config.R")
if (config$cohort != "imv") stop("The VT/PBW gate exists only in the ventilated cohort; unset PBWPFVC_COHORT")
site_name  <- config$site_name
output_dir <- config$output_dir
final_dir  <- final_dir_for("supplement")

VTPBW_BAND      <- c(6, 8)     # mL/kg, the paper's lung-protective gate (03, 3e)
VTPBW_BIN_WIDTH <- 0.5         # mL/kg, the descriptive bins
BOOTSTRAP_REPS  <- 500
set.seed(20260930)

# --- the two frames, one row per patient: the ungated index and the paper's own index
prepare <- function(frame) {
  frame %>%
    mutate(sex_category  = factor(sex_category,  levels = c("Male", "Female")),
           race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
           age10 = age_at_admission / 10,
           log_pfvc = log(pfvc),
           sofa_nonresp = sofa_total - sofa_resp)
}
FRAME_COLUMNS <- c("hospitalization_id", "vtpbw", "pbwpfvc", "pfvc", "age_at_admission", "sex_category",
                   "race_category", "sofa_total", "sofa_resp")
ungated <- read_parquet(file.path(output_dir, "analysis_ungated_index.parquet")) %>%
  select(all_of(FRAME_COLUMNS), in_paper_cohort) %>%
  prepare() %>%
  mutate(index_in_band = vtpbw >= VTPBW_BAND[1] & vtpbw <= VTPBW_BAND[2])
paper <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet"), col_select = all_of(FRAME_COLUMNS)) %>%
  prepare()
# every paper patient has a complete hypoxemic timepoint, so has an ungated index
if (!all(paper$hospitalization_id %in% ungated$hospitalization_id))
  stop(sum(!paper$hospitalization_id %in% ungated$hospitalization_id),
       " patients in the paper's cohort have no ungated index; 03's 3e3 and 3e disagree")
message("=== xsec_vtpbw_gate_collider: ", site_name, ", ", nrow(ungated), " patients at the ungated index, ",
        sum(ungated$index_in_band), " with VT/PBW ", VTPBW_BAND[1], "-", VTPBW_BAND[2], " there; ",
        nrow(paper), " in the paper's cohort")

EXPOSURES <- c("PBW/PFVC" = "pbwpfvc", "log PFVC" = "log_pfvc")
OUTCOMES  <- c("SOFA" = "sofa_total", "Respiratory SOFA" = "sofa_resp", "Non-respiratory SOFA" = "sofa_nonresp")
ADJUSTMENTS <- c(unadjusted = "", adjusted = "+ ns(age10, 4) + sex_category + race_category")
# both frames on the ungated sample's SD, so every design shares one scale
exposure_sd <- map_dbl(EXPOSURES, ~ sd(ungated[[.x]]))
standardise <- function(frame) {
  for (exposure in EXPOSURES) frame[[paste0(exposure, "_sd")]] <- frame[[exposure]] / exposure_sd[[match(exposure, EXPOSURES)]]
  frame
}
ungated <- standardise(ungated)
paper   <- standardise(paper)

# --- the designs: which frame and patients each keeps, and its VT/PBW term
DESIGNS <- tribble(
  ~design,                  ~keep,            ~vtpbw_term,
  "none",                   "ungated",        "",
  "conditioned",            "ungated",        "+ ns(vtpbw, 4)",
  "index 6-8",              "index_in_band",  "",
  "index 6-8 + VT/PBW",     "index_in_band",  "+ vtpbw",
  "paper cohort",           "paper",          "",
  "paper cohort + VT/PBW",  "paper",          "+ vtpbw")
MODELS <- expand_grid(DESIGNS, exposure = names(EXPOSURES), outcome = names(OUTCOMES),
                      adjustment = names(ADJUSTMENTS)) %>%
  mutate(size_term = paste0(EXPOSURES[exposure], "_sd"),
         formula = paste(OUTCOMES[outcome], "~", size_term, vtpbw_term, ADJUSTMENTS[adjustment]))

# the size term's slope in every model, for one set of patients (a bootstrap draw
# repeats a patient in both frames)
fit_all <- function(patient_ids) {
  ungated_rows <- ungated[match(patient_ids, ungated$hospitalization_id), ]
  paper_rows   <- paper[na.omit(match(patient_ids, paper$hospitalization_id)), ]
  samples <- list(ungated = ungated_rows, index_in_band = filter(ungated_rows, index_in_band), paper = paper_rows)
  pmap_dfr(MODELS, function(keep, formula, size_term, ...) {
    coefficients <- summary(lm(as.formula(formula), data = samples[[keep]]))$coefficients
    tibble(estimate = coefficients[size_term, 1], se = coefficients[size_term, 2], n = nrow(samples[[keep]]))
  }) %>% bind_cols(MODELS %>% select(design, exposure, outcome, adjustment), .)
}
estimates <- fit_all(ungated$hospitalization_id)

# --- each design's change from the reference, by patient bootstrap
change_from_reference <- function(fits) {
  fits %>%
    group_by(exposure, outcome, adjustment) %>%
    mutate(change = estimate - estimate[design == "none"]) %>%
    ungroup() %>%
    select(design, exposure, outcome, adjustment, change)
}
message("Bootstrap: ", BOOTSTRAP_REPS, " resamples of ", nrow(ungated), " patients")
bootstrap_changes <- map_dfr(seq_len(BOOTSTRAP_REPS), function(rep) {
  change_from_reference(fit_all(sample(ungated$hospitalization_id, replace = TRUE))) %>% mutate(rep = rep)
})
change_summary <- bootstrap_changes %>%
  group_by(design, exposure, outcome, adjustment) %>%
  summarise(change_se = sd(change), change_lo = quantile(change, 0.025), change_hi = quantile(change, 0.975),
            .groups = "drop")
estimates <- estimates %>%
  left_join(change_from_reference(estimates), by = c("design", "exposure", "outcome", "adjustment")) %>%
  left_join(change_summary, by = c("design", "exposure", "outcome", "adjustment")) %>%
  mutate(across(starts_with("change"), ~ if_else(design == "none", NA_real_, .x)),
         exposure_sd = exposure_sd[exposure], site = site_name, .before = 1) %>%
  mutate(design = factor(design, levels = DESIGNS$design)) %>%
  arrange(exposure, outcome, adjustment, design)
write_csv(estimates, file.path(final_dir, paste0("vtpbw_gate_collider_estimates_", site_name, ".csv")))

# --- the collider's two arms, in every patient: VT/PBW on the size term and on SOFA
arms <- expand_grid(exposure = names(EXPOSURES), outcome = names(OUTCOMES), adjustment = names(ADJUSTMENTS)) %>%
  pmap_dfr(function(exposure, outcome, adjustment) {
    size_term <- paste0(EXPOSURES[[exposure]], "_sd")
    sofa_term <- OUTCOMES[[outcome]]
    model <- lm(as.formula(paste("vtpbw ~", size_term, "+", sofa_term, ADJUSTMENTS[[adjustment]])), data = ungated)
    coefficients <- summary(model)$coefficients
    tibble(site = site_name, exposure, outcome, adjustment,
           arm = c("size term (per SD)", "SOFA (per point)"),
           estimate = coefficients[c(size_term, sofa_term), 1],
           se = coefficients[c(size_term, sofa_term), 2],
           r_squared = summary(model)$r.squared, n = stats::nobs(model))
  })
write_csv(arms, file.path(final_dir, paste0("vtpbw_gate_collider_arms_", site_name, ".csv")))

# --- descriptive: the ratio, log PFVC and SOFA across VT/PBW, in bins
bin_breaks <- seq(floor(min(ungated$vtpbw) / VTPBW_BIN_WIDTH) * VTPBW_BIN_WIDTH,
                  ceiling(max(ungated$vtpbw) / VTPBW_BIN_WIDTH) * VTPBW_BIN_WIDTH + VTPBW_BIN_WIDTH,
                  by = VTPBW_BIN_WIDTH)
bins <- ungated %>%
  mutate(vtpbw_bin_lower = bin_breaks[findInterval(vtpbw, bin_breaks)]) %>%
  group_by(vtpbw_bin_lower) %>%
  summarise(n = n(), n_paper_cohort = sum(in_paper_cohort),
            across(c(pbwpfvc, log_pfvc, sofa_total, sofa_resp, sofa_nonresp),
                   list(mean = mean, sd = sd), .names = "{.col}_{.fn}"),
            .groups = "drop") %>%
  mutate(vtpbw_bin_upper = vtpbw_bin_lower + VTPBW_BIN_WIDTH, .after = vtpbw_bin_lower) %>%
  mutate(site = site_name, .before = 1)
write_csv(bins, file.path(final_dir, paste0("vtpbw_gate_collider_bins_", site_name, ".csv")))

# --- the gate's counts
gate <- ungated %>%
  mutate(index_vtpbw = case_when(vtpbw < VTPBW_BAND[1] ~ "below 6", vtpbw > VTPBW_BAND[2] ~ "above 8",
                                 TRUE ~ "6-8")) %>%
  count(index_vtpbw, in_paper_cohort, name = "n") %>%
  mutate(site = site_name, .before = 1)
write_csv(gate, file.path(final_dir, paste0("vtpbw_gate_collider_gate_", site_name, ".csv")))

# --- figure: the size term's SOFA slope under each design
ADJUSTMENT_COLOURS <- c(unadjusted = "#E69F00", adjusted = "#0072B2")    # Okabe-Ito
collider_plot <- estimates %>%
  mutate(design = fct_rev(design)) %>%
  ggplot(aes(estimate, design, colour = adjustment)) +
  geom_vline(data = estimates %>% filter(design == "none"),
             aes(xintercept = estimate, colour = adjustment), linetype = "dashed", linewidth = 0.3) +
  geom_pointrange(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se),
                  position = position_dodge(width = 0.5), size = 0.25) +
  facet_grid(exposure ~ outcome, scales = "free_x") +
  scale_colour_manual(values = ADJUSTMENT_COLOURS, name = NULL) +
  labs(x = "SOFA points per SD of the size term (dashed: every patient, no VT/PBW term)", y = NULL,
       title = paste0("Does the VT/PBW gate bend the size-SOFA association? (", site_name, ")")) +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(final_dir, paste0("vtpbw_gate_collider_", site_name, ".pdf")), collider_plot, width = 9, height = 5)

message("Wrote vtpbw_gate_collider_* to ", final_dir)
print(estimates %>% filter(outcome == "SOFA") %>%
        select(exposure, adjustment, design, estimate, change, change_lo, change_hi, n), n = Inf)
