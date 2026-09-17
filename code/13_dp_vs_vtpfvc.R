# =============================================================================
# Script 13 (mechanics): driving pressure on VT/PFVC with a spline in VT/PBW
# =============================================================================
# Driving pressure is tidal volume times elastance. If lung size sets elastance,
# elastance scales with 1 / PFVC and driving pressure is proportional to VT/PFVC:
# at a fixed VT/PBW (spline), the remaining variation in VT/PFVC is the PBW/PFVC
# discordance, and a positive coefficient says the lung that is smaller than
# PBW predicts is stiffer for the same dose per kilogram. If PBW were the right
# normalizer, VT/PFVC would add nothing once the dose spline is in. The mirror
# model (VT/PBW linear, spline in VT/PFVC) closes the argument: whichever term
# survives its rival's spline carries the mechanics. A head-to-head of the two
# normalizers with measured pressures, on the index timepoint.
#
#   dp ~ vtpfvc (percent of predicted FVC) + ns(vtpbw, 3) + bmi
#        [+ log SF + non-respiratory SOFA]  [+ ns(age, 4) + sex + race]
#   dp ~ vtpbw + ns(vtpfvc, 3) + bmi [+ ...]                           (mirror)
#
# BMI is in every model (driving pressure includes the chest wall). Fitted in
# three strata of the index VT/PBW: 6-8 (the analytic cohort), 5-8 (the widened
# band), and < 5 (the titrated, stiff-lung population, reported separately, not
# pooled). Index = the patient's first complete, hypoxemic timepoint with a
# recorded plateau in the stratum's band. Titration caveat: clinicians lower
# VT/PBW when driving pressure is high, so conditioning on the dose biases the
# VT/PFVC coefficient toward zero; the index timepoint and the 5-8 band limit it.
#
# Outputs (aggregates only):
#   final/dp_vs_vtpfvc_{site}.csv  per stratum x adjustment set x model: the
#                                  exposure coefficient (cmH2O per point of
#                                  predicted FVC, or per mL/kg), CI, p, n, R2
#   final/dp_vs_vtpfvc_{site}.pdf  forest of both exposures across strata
# Usage: Rscript code/13_dp_vs_vtpfvc.R
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(here); library(splines); library(patchwork) })
rm(list = ls())
source("utils/config.R")
site_name  <- config$site_name
output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
okabe <- c("#0072B2", "#E69F00", "#009E73", "#D55E00")
SF_HYPOXEMIA_THRESHOLD <- 315
STRATA <- list("6-8 (analytic)" = c(6, 8), "5-8 (widened)" = c(5, 8), "< 5 (titrated)" = c(0, 5))

ae <- read_parquet(file.path(output_dir, "analysis_all_eligible_timepoints.parquet")) %>%
  filter(has_all_data, !is.na(dp), dp > 0, sf_ratio < SF_HYPOXEMIA_THRESHOLD)
weights <- read_parquet(file.path(output_dir, "cohort_weights.parquet"))
ae <- ae %>% left_join(weights, by = "hospitalization_id") %>%
  mutate(bmi = if_else(!is.na(weight_kg) & height_cm > 0, weight_kg / (height_cm / 100)^2, NA_real_),
         np_sofa = if ("sofa_resp" %in% names(.)) sofa_total - sofa_resp else sofa_total,
         log_sf = log(sf_ratio), age10 = age_at_admission / 10,
         vtpfvc_pct = vtpfvc)   # percent of predicted FVC (script 03: VT mL / PFVC L x 0.1)
message("Complete, hypoxemic timepoints with a recorded plateau: ", nrow(ae), " on ", n_distinct(ae$hospitalization_id), " patients")

index_in <- function(band) ae %>% filter(vtpbw >= band[1], vtpbw < band[2] | (band[2] == 8 & vtpbw <= 8)) %>%
  group_by(hospitalization_id) %>% slice_min(recorded_dttm, n = 1, with_ties = FALSE) %>% ungroup() %>%
  filter(!is.na(bmi), !is.na(np_sofa), is.finite(log_sf), !is.na(age10), !is.na(sex_category), !is.na(race_category))

SETS <- list(none = NULL, severity = c("log_sf", "np_sofa"),
             demographics = c("ns(age10, 4)", "sex_category", "race_category"),
             both = c("log_sf", "np_sofa", "ns(age10, 4)", "sex_category", "race_category"))
fit_one <- function(d, exposure, spline_of, extra, stratum, set) {
  f <- as.formula(paste("dp ~", exposure, "+ ns(", spline_of, ", 3) + bmi", if (length(extra)) paste("+", paste(extra, collapse = " + ")) else ""))
  m <- lm(f, data = d); co <- summary(m)$coefficients[exposure, ]; ci <- confint(m)[exposure, ]
  tibble(stratum = stratum, adjustment = set, exposure = exposure,
         unit = if (exposure == "vtpfvc_pct") "cmH2O per point of predicted FVC (VT/PBW spline held)" else "cmH2O per mL/kg PBW (VT/PFVC spline held)",
         estimate = unname(co[1]), se = unname(co[2]), lo = unname(ci[1]), hi = unname(ci[2]), p = unname(co[4]),
         n = nrow(d), r2 = summary(m)$r.squared,
         # the exposure's own contrast across its interquartile range, for scale
         iqr = unname(diff(quantile(d[[exposure]], c(.25, .75)))), estimate_per_iqr = unname(co[1]) * unname(diff(quantile(d[[exposure]], c(.25, .75)))))
}
res <- imap_dfr(STRATA, function(band, stratum) {
  d <- index_in(band)
  if (nrow(d) < 50) { message(stratum, ": ", nrow(d), " patients, skipped"); return(NULL) }
  message(stratum, ": ", nrow(d), " patients; VT/PBW median ", signif(median(d$vtpbw), 3), ", VT/PFVC median ", signif(median(d$vtpfvc_pct), 3), "%")
  imap_dfr(SETS, function(extra, set) bind_rows(
    fit_one(d, "vtpfvc_pct", "vtpbw", extra, stratum, set),
    fit_one(d, "vtpbw", "vtpfvc_pct", extra, stratum, set)))
}) %>% mutate(site = site_name)
write_csv(res, file.path(final_dir, paste0("dp_vs_vtpfvc_", site_name, ".csv")))
message("\nDriving pressure on each normalizer's exposure with a spline in the other (cmH2O per unit; positive = stiffer for the same dose):")
print(as.data.frame(res %>% select(stratum, adjustment, exposure, estimate, lo, hi, p, estimate_per_iqr, n) %>%
                      mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)

# forest: the two exposures side by side, one panel per stratum, adjustment sets down the axis
fd <- res %>% mutate(adjustment = factor(adjustment, names(SETS)), stratum = factor(stratum, names(STRATA)),
                     exposure = factor(if_else(exposure == "vtpfvc_pct", "VT/PFVC (per point of predicted FVC)", "VT/PBW (per mL/kg)")))
p <- ggplot(fd, aes(estimate, adjustment, colour = exposure)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5)) +
  facet_grid(stratum ~ exposure, scales = "free_x") +
  scale_colour_manual(values = okabe[c(4, 1)], guide = "none") +
  labs(title = sprintf("Driving pressure on VT/PFVC (spline in VT/PBW) and on VT/PBW (spline in VT/PFVC), %s", site_name),
       subtitle = "index timepoint with a recorded plateau; BMI in every model; positive = stiffer lung for the same dose by the other normalizer",
       x = "coefficient (cmH2O per unit, 95% CI)", y = "adjustment set") +
  theme_minimal(base_size = 10) + theme(strip.text.y = element_text(angle = 0))
ggsave(file.path(final_dir, paste0("dp_vs_vtpfvc_", site_name, ".pdf")), p, width = 10, height = 2 + 1.6 * n_distinct(fd$stratum))
message("13_dp_vs_vtpfvc complete -> ", final_dir)
