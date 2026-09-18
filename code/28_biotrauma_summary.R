# =============================================================================
# Script 28 (summary): one table from the overnight biotrauma run
# =============================================================================
# Reads only the site's final/ CSVs and writes overnight_summary_{site}.csv, one
# row per estimate in long format (source file, family, marker, form, panel
# horizon, contrast horizon, model / adjustment, exposure or term, estimate,
# interval, p, status), and prints the fit status table from every manifest so
# the morning read is one file and one screen.
#
# Usage: Rscript code/28_biotrauma_summary.R   (PBWPFVC_FIG_DIR points at another final/)
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(here) })
rm(list = ls())
source("utils/config.R")
site_name <- config$site_name
final_dir <- Sys.getenv("PBWPFVC_FIG_DIR", here("output", paste0(site_name, "_output"), "final"))

files <- list.files(final_dir, "\\.csv$", full.names = TRUE)
fam_of <- function(f) str_match(basename(f), "^(jm_manifest|jm_level_contrast|injury_at_horizon_counts|injury_at_horizon|injury_channels|injury_nested|injury_dose_channels|injury_sf_channels|injury_negctrl|quick_lme|quick_channels|quick_nested|quick_dose_channels|quick_sf_channels)_")[, 2]
fam <- fam_of(files)
keep <- !is.na(fam) & fam != "injury_at_horizon_counts" & !grepl("^quick_lme_[a-z_]+_\\d+h_", basename(files))   # stale horizon-tagged quick files
files <- files[keep]; fam <- fam[keep]
if (!length(files)) stop("no biotrauma tables in ", final_dir)
rd <- function(f) read_csv(f, show_col_types = FALSE, guess_max = 1e5)
jm_tag <- function(f) {
  m <- str_match(basename(f), "^jm_\\w+?_(?:(offset)_)?(?:(pfvc|channels|disc_level|vtpfvc|saturated|none)_)?(\\d+[hd])_")
  list(form = replace_na(m[, 3], "disc"), panel = m[, 4], baseline = replace_na(m[, 2], "free"))
}
num <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)
chr <- function(x) if (is.null(x)) NA_character_ else as.character(x)

rows <- list(); manifests <- list()
for (i in seq_along(files)) {
  f <- files[i]; d <- rd(f); fm <- fam[i]
  if (!nrow(d)) next
  base <- tibble(source_file = basename(f), family = fm)
  r <- switch(fm,
    jm_manifest = { t <- jm_tag(f); manifests[[length(manifests) + 1]] <- d %>% mutate(form = t$form, panel_h = t$panel); NULL },
    jm_level_contrast = { t <- jm_tag(f)
      tibble(marker = d$marker, form = t$form, panel_h = t$panel, horizon_h = num(d$horizon_h), model = chr(d$adjustment),
             exposure_or_term = d$exposure, estimate = d$estimate, lo = d$lo, hi = d$hi, p = num(d$p_equal), p_gt0 = num(d$p_gt0),
             status = if_else(d$rhat_gate, "converged", "rhat_fail")) },
    injury_at_horizon = tibble(marker = d$marker, form = "comparator", panel_h = NA, horizon_h = num(d$horizon_h),
                               model = paste(d$outcome_type, d$adjustment, sep = " / "), exposure_or_term = d$exposure,
                               estimate = d$estimate, lo = d$lo, hi = d$hi, p = num(d$p), p_gt0 = NA, status = "ok"),
    injury_channels = tibble(marker = d$marker, form = "comparator channels", panel_h = NA, horizon_h = num(d$horizon_h),
                             model = paste(d$outcome_type, d$exposure, sep = " / "), exposure_or_term = d$channel,
                             estimate = d$estimate, lo = d$lo, hi = d$hi, p = num(d$p_equal), p_gt0 = NA, status = "ok"),
    injury_nested = tibble(marker = d$marker, form = "comparator nested", panel_h = NA, horizon_h = num(d$horizon_h),
                           model = paste(d$small, d$big, sep = " -> "), exposure_or_term = d$test,
                           estimate = d$d_aic, lo = NA, hi = NA, p = num(d$p), p_gt0 = NA, status = "ok"),
    injury_dose_channels = , injury_sf_channels = { if (!"estimate" %in% names(d)) NULL else
      tibble(marker = d$marker, form = paste("comparator", sub("^injury_", "", fm)), panel_h = NA, horizon_h = num(d$horizon_h),
             model = d$model, exposure_or_term = d$term, estimate = d$estimate, lo = d$lo, hi = d$hi, p = num(d$p), p_gt0 = NA, status = "ok") },
    injury_negctrl = tibble(marker = d$marker, form = "comparator negative control", panel_h = NA, horizon_h = 0,
                            model = d$model, exposure_or_term = d$term, estimate = d$estimate, lo = d$lo, hi = d$hi,
                            p = num(d$p), p_gt0 = NA, status = if_else(d$separation_flag, "separation", "ok")),
    quick_lme = tibble(marker = d$marker, form = "quick lme", panel_h = NA, horizon_h = num(d$model_horizon_h),
                       model = chr(d$adjustment), exposure_or_term = d$exposure, estimate = d$estimate, lo = d$lo, hi = d$hi,
                       p = NA, p_gt0 = NA, status = "ok"),
    quick_channels = tibble(marker = d$marker, form = "quick channels", panel_h = NA, horizon_h = num(d$model_horizon_h),
                            model = d$exposure, exposure_or_term = d$channel, estimate = d$estimate, lo = d$lo, hi = d$hi,
                            p = num(d$p_equal), p_gt0 = NA, status = "ok"),
    quick_nested = tibble(marker = d$marker, form = "quick nested", panel_h = NA, horizon_h = num(d$model_horizon_h),
                          model = paste(d$exposure, d$small, d$big, sep = " / "), exposure_or_term = d$test,
                          estimate = d$d_aic, lo = NA, hi = NA, p = num(d$p), p_gt0 = NA, status = "ok"),
    quick_dose_channels = , quick_sf_channels = { if (!"estimate" %in% names(d)) NULL else
      tibble(marker = d$marker, form = paste("quick", sub("^quick_", "", fm)), panel_h = NA, horizon_h = num(d$model_horizon_h),
             model = paste(d$exposure, d$model, sep = " / "), exposure_or_term = d$term, estimate = d$estimate, lo = d$lo, hi = d$hi,
             p = num(d$p), p_gt0 = NA, status = "ok") })
  if (!is.null(r)) rows[[length(rows) + 1]] <- bind_cols(base[rep(1, nrow(r)), ], r)
}
summary_tbl <- bind_rows(rows) %>% mutate(site = site_name)
write_csv(summary_tbl, file.path(final_dir, paste0("overnight_summary_", site_name, ".csv")))
message("overnight_summary: ", nrow(summary_tbl), " rows from ", length(files), " tables -> ", final_dir)

if (length(manifests)) {
  st <- bind_rows(manifests) %>%
    select(any_of(c("panel_h", "form", "marker", "adjustment", "status", "n_patients", "n_deaths",
                    "longitudinal_rhat", "association_rhat", "hazard_rhat", "max_rhat", "worst_terms", "reason"))) %>%
    arrange(panel_h, form, marker)
  message("\n---- joint-model fits")
  print(as.data.frame(st %>% select(-any_of(c("worst_terms", "reason"))) %>%
                        mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
  bad <- st %>% filter(status %in% c("rhat_fail", "failed"))
  if (nrow(bad)) {
    message("\n---- worst terms / reasons of the fits that did not converge")
    for (i in seq_len(nrow(bad))) message(sprintf("  %s %s %s %s: %s", bad$panel_h[i], bad$form[i], bad$marker[i], bad$adjustment[i],
                                                   substr(coalesce(bad$worst_terms[i], bad$reason[i], ""), 1, 110)))
  }
  message("\n", sum(st$status == "converged"), " converged, ", sum(st$status == "rhat_fail"), " R-hat fail, ",
          sum(st$status == "failed"), " failed, ", sum(st$status == "skipped"), " skipped")
  if ("longitudinal_rhat" %in% names(st))
    message("by block (R-hat <= 1.1): longitudinal ", sum(st$longitudinal_rhat <= 1.1, na.rm = TRUE), " of ", sum(is.finite(st$longitudinal_rhat)),
            ", association ", sum(st$association_rhat <= 1.1, na.rm = TRUE), " of ", sum(is.finite(st$association_rhat)),
            ", hazard ", sum(st$hazard_rhat <= 1.1, na.rm = TRUE), " of ", sum(is.finite(st$hazard_rhat)),
            "  (the trajectory figures need the longitudinal block; the hazard ratios need the hazard block)")
}
message("\n---- PFVC-level contrasts at 48 h (joint model, main, adjusted / channels; log units)")
print(as.data.frame(summary_tbl %>% filter(family == "jm_level_contrast", horizon_h == 48, model %in% c("adjusted", "channels")) %>%
                      select(marker, form, model, exposure_or_term, estimate, lo, hi, p_gt0, status) %>%
                      mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
