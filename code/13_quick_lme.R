# =============================================================================
# Script 13 (quick): the longitudinal submodel alone, no joint model
# =============================================================================
# Fits the PFVC-level question as a plain linear mixed model on the panel
# tables 13_biotrauma_panel.R wrote, in seconds:
#
#   log marker ~ time + log PFVC (per SD) + log PFVC x time
#                + dose (patient mean, within-patient change)
#                + baseline marker + lagged SF and pressor + non-respiratory
#                  SOFA + BMI  [+ ns(age, 4) + sex + race]
#   random intercept and slope per patient (pdDiag)
#
# and prints the marker difference per SD of log PFVC at 24, 48 and 72 hours
# (level + divergence x time). This is the joint model's longitudinal part
# without the survival linkage, so death before the horizon is NOT accounted
# for: it is the fast look, and the joint-model contrast is the read.
#
# Usage: PBWPFVC_INJ_MARKER=creatinine Rscript code/13_quick_lme.R
#        (PBWPFVC_JM_GRID / _HORIZON_H select the panel, as for the fit;
#         PBWPFVC_QUICK_EXPO=ldisc_sd swaps in log PBW/PFVC)
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(here); library(splines); library(nlme) })
rm(list = ls())
source("utils/config.R")
site_name  <- config$site_name
output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
source(here("code", "13_biotrauma_grid.R"))

MARKER <- Sys.getenv("PBWPFVC_INJ_MARKER", "creatinine")
EXPO   <- Sys.getenv("PBWPFVC_QUICK_EXPO", "log_pfvc_sd")
y_col  <- c(creatinine = "creatinine", ne_equiv = "ne_equiv_peak", platelets = "platelets",
            bilirubin = "bilirubin", sf = "sf", dp = "dp")[[MARKER]]
y0_col <- c(creatinine = "creatinine_0", ne_equiv = "ne_equiv_0", platelets = "platelet_0",
            bilirubin = "bilirubin_0", sf = "sf_0", dp = "dp_0")[[MARKER]]
offset <- if (MARKER == "ne_equiv") 0.01 else 0
own_lag <- c(sf = "l_log_sf", ne_equiv = "l_pressor")[MARKER]

long <- read_parquet(file.path(output_dir, paste0("jm_long_", h_suffix, ".parquet")))
surv <- read_parquet(file.path(output_dir, paste0("jm_surv_", h_suffix, ".parquet")))
d <- long %>%
  filter(period >= 1L, !is.na(.data[[y_col]]), !is.na(l_vtpbw_within), !is.na(l_sf), !is.na(l_pressor)) %>%
  inner_join(surv %>% select(hospitalization_id, np_sofa, bmi, age10, sex_category, race_category,
                             vtpbw_pt_mean, log_pfvc_sd, ldisc_sd, all_of(y0_col)), by = "hospitalization_id") %>%
  filter(!is.na(.data[[y0_col]]), !is.na(np_sofa), !is.na(bmi)) %>%
  mutate(log_y = log(.data[[y_col]] + offset), log_y0 = log(.data[[y0_col]] + offset), l_log_sf = log(l_sf)) %>%
  group_by(hospitalization_id) %>% filter(n() >= 2L) %>% ungroup() %>%
  mutate(id = factor(hospitalization_id))
message(MARKER, " on the ", h_suffix, " grid: ", nrow(d), " rows, ", n_distinct(d$id), " patients")

lags <- setdiff(c("l_log_sf", "l_pressor"), own_lag)
rhs  <- function(adjusted) paste(c("vent_day", EXPO, paste0(EXPO, ":vent_day"), "l_vtpbw_within", "vtpbw_pt_mean",
                                   "log_y0", lags, "np_sofa", "bmi",
                                   if (adjusted) c("ns(age10, 4)", "sex_category", "race_category")), collapse = " + ")
fits <- map(c(adjusted = TRUE, unadjusted = FALSE), function(adj)
  lme(as.formula(paste("log_y ~", rhs(adj))), random = list(id = pdDiag(~ vent_day)), data = d,
      control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)))

hours <- c(24, 48, 72); hours <- hours[hours <= JM_HORIZON * 24]
out <- imap_dfr(fits, function(f, adj) {
  b <- fixef(f); V <- vcov(f)
  tn <- intersect(c(paste0(EXPO, ":vent_day"), paste0("vent_day:", EXPO)), names(b))
  map_dfr(hours, function(hh) {
    w <- setNames(rep(0, length(b)), names(b)); w[EXPO] <- 1; w[tn] <- hh / 24
    est <- sum(w * b); se <- sqrt(as.numeric(t(w) %*% V %*% w))
    tibble(marker = MARKER, exposure = EXPO, adjustment = adj, horizon_h = hh,
           estimate = est, lo = est - 1.96 * se, hi = est + 1.96 * se,
           level = unname(b[EXPO]), divergence_per_day = unname(if (length(tn)) b[tn] else 0),
           dose_within = unname(b["l_vtpbw_within"]), n_patients = n_distinct(d$id), n_rows = nrow(d))
  })
})
message("\nMarker difference per SD of ", EXPO, " (log units; a lower PFVC is the negative of this). ",
        "No correction for death before the horizon: compare with jm_level_contrast_*.")
print(as.data.frame(out %>% mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
write_csv(out %>% mutate(grid = JM_GRID, site = site_name),
          file.path(final_dir, paste0("quick_lme_", MARKER, "_", h_suffix, "_", site_name, ".csv")))
