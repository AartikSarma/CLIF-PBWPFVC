# =============================================================================
# DEFERRED: Federated conditional-bias export (per site)
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Standalone, NOT part of the 00 pipeline runner. Produces the per-site aggregate
# that drives the POOLED cross-cohort conditional-bias plots
# (code/cbias_pooled_plots.R).
#
# *** DO NOT RUN / SHARE YET ***
# After stratifying by demographics, some (stratum x percentile) cells fall below
# n = 10. The resulting aggregate must therefore be passed through the consortium's
# deterministic additive-masking pipeline before it can leave a site. That step is
# held until ALL sites have confirmed they will run the code, so this export is
# deferred. Run it (and apply masking) only once that confirmation is in place.
#
# For each (plot, dependent var, independent var, demographic stratum, percentile)
# the export records the count and the sum / sum-of-squares of the dependent
# variable. Percentiles are ntiles of the independent variable over ALL subjects,
# so a stratum's per-percentile counts show where that group sits in the
# distribution (the bias signal, and the input to the pooled density band). The
# pooled script forms N-weighted percentile means, aggregates ten percentiles per
# decile for the circle/CI points, and weights the loess from these.
#
# Inputs : output/<site>/intermediate/analysis_cross_sectional.parquet (script 03)
# Outputs: final/cbias_export_<site>.csv  (must be masked before sharing)
# =============================================================================

library(tidyverse)
library(arrow)
library(here)

source("utils/config.R")
site_name <- config$site_name

output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

# Rebuild the diagnostic frame exactly as script 04 does for the per-site plots.
diag_data <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet")) %>%
  mutate(age_group = cut(age_at_admission,
                         breaks = c(18, 40, 60, 80, Inf),
                         labels = c("18-39", "40-59", "60-79", "80+"),
                         right = FALSE)) %>%
  filter(!is.na(age_group))

# Same dependent x independent configurations as the per-site plots (script 04, 4h).
cbias_configs <- list(
  list(id = "mortality",     dep = "deceased", dep_label = "Mortality",
       indep = c("vtpfvc", "vtpbw", "pbwpfvc")),
  list(id = "pfvc_vs_pbw",   dep = "pfvc",     dep_label = "PFVC (L)",
       indep = "pbw"),
  list(id = "elastance",     dep = "ers",      dep_label = "Elastance (cmH2O/L)",
       indep = c("vtpfvc", "vtpbw")),
  list(id = "compliance",    dep = "crs",      dep_label = "Compliance (mL/cmH2O)",
       indep = "pbwpfvc"),
  list(id = "vfd28",         dep = "vfd_28",   dep_label = "28-day VFDs",
       indep = "pbwpfvc"),
  list(id = "pbwpfvc_ratio", dep = "pbwpfvc",  dep_label = "PBW/PFVC Ratio",
       indep = c("pfvc", "pbw"))
)
cbias_grouping    <- c("race_category", "sex_category", "age_group")
CBIAS_PERCENTILES <- 100

# Per-(stratum, percentile) aggregate for one dependent x independent x grouping
# cell. NO suppression here — the small cells are handled by the downstream
# additive-masking pipeline (see header).
cbias_percentile_export <- function(data, dep, dep_label, indep, grp, plot_id) {
  d <- data %>%
    filter(!is.na(.data[[dep]]), !is.na(.data[[indep]]), !is.na(.data[[grp]])) %>%
    mutate(percentile = ntile(.data[[indep]], CBIAS_PERCENTILES),
           stratum    = as.character(.data[[grp]]),
           depv       = as.numeric(.data[[dep]]))
  if (nrow(d) == 0) return(tibble())
  is_binary <- all(d$depv %in% c(0, 1))
  d %>%
    group_by(stratum, percentile) %>%
    summarise(n = n(), sum_dep = sum(depv), sumsq_dep = sum(depv^2), .groups = "drop") %>%
    transmute(plot_id, dependent = dep, dep_label, independent = indep, grouping = grp,
              stratum, percentile, n, sum_dep, sumsq_dep, is_binary)
}

cbias_export <- map_dfr(cbias_configs, function(cfg) {
  map_dfr(cfg$indep, function(iv) {
    map_dfr(cbias_grouping, function(gv) {
      cbias_percentile_export(diag_data, cfg$dep, cfg$dep_label, iv, gv, cfg$id)
    })
  })
}) %>% mutate(site = site_name, .before = 1)

n_small <- sum(cbias_export$n < 10)
write_csv(cbias_export, file.path(final_dir, paste0("cbias_export_", site_name, ".csv")))
message("Conditional-bias federated export: ", nrow(cbias_export),
        " (plot x independent x stratum x percentile) rows; ", n_small,
        " cells with n < 10 — APPLY THE ADDITIVE-MASKING PIPELINE BEFORE SHARING.")
