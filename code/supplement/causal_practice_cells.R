# =============================================================================
# Supplement (causal): is there enough exogenous variation in how practice scales tidal
# volume with PBW? A diagnostic for a unit x period design. No outcome is read.
# =============================================================================
#
# The design this prepares (not fitted here). Each practice cell c (an ICU in a period)
# sets tidal volume by its own rule, log VT = a_c + rho_c log PBW: rho_c near 1 dials
# mL/kg by PBW, near 0 gives everyone much the same volume. Under strain the marker's
# slope is b a_c + b rho_c log PBW - b log PFVC + ..., so in
#     slope ~ cell + beta1 log PBW + beta2 (rho_c x log PBW) + beta3 log PFVC + ...
# the cell term absorbs everything a cell does on average (case mix, quality, the level
# a_c), beta1 the direct effects of height and sex, beta3 the age confounding of PFVC,
# and beta2 is identified only by cells differing in how strongly tidal volume follows
# PBW. Strain predicts beta2 = b = -beta3. The instrument is rho_c x (log PBW about the
# cell's mean), with rho_c estimated without the patient's own data and shrunken toward
# the mean over cells.
#
# Checked on synthetic data with between-hospital differences in rho planted (SD 0.24):
# tau came back 0.25 and the first-stage coefficient 0.97. The synthetic data as they
# come have no such differences (tau = 0), and the script then says so.
#
# Whether that can work rests on three things this script measures, before any outcome
# is looked at:
#   1. Do cells differ in rho at all, beyond noise? Per-cell rho with its standard
#      error, the between-cell SD tau and I^2 (DerSimonian-Laird), and the F test that
#      the slope of log VT on log PBW is the same in every cell.
#   2. Is the instrument strong? The first-stage F of the leave-one-out, shrunken
#      rho_c x (log PBW about the cell mean) for log VT, given cell, log PBW and log PFVC
#      (cluster-robust on the cell); its coefficient should be near 1. Below about 10 the design is weak; 20 or more is comfortable.
#   3. Is it plausibly exogenous? Balance: the same first-stage model with a baseline
#      characteristic in place of log VT (SOFA, S/F, age, BMI, female). A design that
#      works has instrument coefficients near 0 on all of them. And whether a cell's rho
#      goes with its case mix (cell-level correlations).
# A period effect on rho alone (practice changing over time everywhere) is reported
# too: it is variation, but variation shared with everything else that changed then.
#
# Cells. The ICU where the patient was at the index (ADT row covering index_dttm), or,
# for an index outside an ICU (ED, operating room), the first ICU entered within 24
# hours; patients with neither are left out (unit_source says how each was placed).
# The unit is the site's own location_name within hospital_id (PBWPFVC_CELL_UNIT=
# location_type uses the mCIDE ICU type instead, coarser and comparable across sites).
# Periods are calendar years of the index (PBWPFVC_CELL_YEARS=2 groups two years).
# Cells with fewer than PBWPFVC_CELL_MIN patients (default 30) are left out of the
# estimates and counted.
#
# Cohort. Run it UNGATED (PBWPFVC_VTPBW_GATE=0): the default cohort keeps only an index
# at VT/PBW 6-8 mL/kg, which truncates exactly the variation in question. The gated
# run is reported as such in the cohort column.
#
# Version. Reads clif_adt from config$tables_path: hospital_id, location_name and
# location_category (CLIF 2.1 and 3.0), location_type where present (needed only for
# PBWPFVC_CELL_UNIT=location_type). location_category is compared as a normalised key
# ("icu"); the source values are kept in the output.
#
# Writes to final/supplement/ (final/ungated/supplement/ with PBWPFVC_VTPBW_GATE=0):
#   causal_practice_cells_<site>.csv     one row per cell: n, VT/PBW level, rho (raw and
#                                        shrunken), case mix
#   causal_practice_summary_<site>.csv   heterogeneity, first-stage F, balance
#   causal_practice_placement_<site>.csv how patients were placed in cells, and lost
#   causal_practice_<site>.pdf           rho by cell and period
# Cell-level aggregates only: no patient-level row. Small cells are not masked here
# (see utils/config.R): a site masks final/ before sharing it.
#
# Usage (bash/zsh):  PBWPFVC_VTPBW_GATE=0 uvr run code/supplement/causal_practice_cells.R
#       (PowerShell): $env:PBWPFVC_VTPBW_GATE = "0"; uvr run code/supplement/causal_practice_cells.R
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
cohort <- if (config$vtpbw_gate == "0") "ungated" else "gated"

CELL_UNIT  <- Sys.getenv("PBWPFVC_CELL_UNIT", "location_name")
CELL_YEARS <- as.integer(Sys.getenv("PBWPFVC_CELL_YEARS", "1"))
CELL_MIN   <- as.integer(Sys.getenv("PBWPFVC_CELL_MIN", "30"))
ICU_AFTER_H <- 24
stopifnot(CELL_UNIT %in% c("location_name", "location_type"), CELL_YEARS >= 1, CELL_MIN >= 10)

# ---- The ventilated cohort at its index ----------------------------------------------
cs_file <- file.path(config$output_dir, paste0("analysis_cross_sectional", config$cs_suffix, ".parquet"))
if (!file.exists(cs_file)) stop("missing ", cs_file, ": run scripts 01-03 first")
vent <- read_parquet(cs_file) %>%
  select(hospitalization_id, index_dttm, vtpbw, pbw, pfvc, age_at_admission, sex_category, sofa_total, sf_ratio, any_of("bmi")) %>%
  filter(!is.na(index_dttm), vtpbw > 0, pbw > 0, pfvc > 0)
if (!"bmi" %in% names(vent)) vent$bmi <- NA_real_
message("Ventilated patients with an index tidal volume: ", nrow(vent), " (", cohort, ")")

# ---- The ICU at the index -----------------------------------------------------------
adt_file <- file.path(config$tables_path, paste0("clif_adt.", config$file_type))
if (!file.exists(adt_file)) stop("missing ", adt_file, " (config$tables_path)")
adt <- switch(config$file_type,
              parquet = read_parquet(adt_file), csv = readr::read_csv(adt_file, show_col_types = FALSE),
              fst = fst::read_fst(adt_file)) %>%
  filter(hospitalization_id %in% vent$hospitalization_id)
need <- c("hospitalization_id", "in_dttm", "out_dttm", "location_category", "location_name")
if (!all(need %in% names(adt))) stop("clif_adt lacks ", paste(setdiff(need, names(adt)), collapse = ", "))
if (CELL_UNIT == "location_type" && !"location_type" %in% names(adt))
  stop("PBWPFVC_CELL_UNIT=location_type needs clif_adt$location_type, which this site's ADT table does not have")
if (!"hospital_id" %in% names(adt)) { adt$hospital_id <- "site"; message("clif_adt has no hospital_id: one hospital assumed") }
icu <- adt %>%
  filter(trimws(tolower(location_category)) == "icu", !is.na(in_dttm)) %>%
  mutate(unit = .data[[CELL_UNIT]]) %>%
  select(hospitalization_id, hospital_id, unit, in_dttm, out_dttm)

at_index <- vent %>% select(hospitalization_id, index_dttm) %>%
  inner_join(icu, by = "hospitalization_id", relationship = "many-to-many") %>%
  filter(in_dttm <= index_dttm, is.na(out_dttm) | out_dttm > index_dttm) %>%
  arrange(hospitalization_id, desc(in_dttm)) %>% distinct(hospitalization_id, .keep_all = TRUE) %>%
  mutate(unit_source = "ICU at the index")
after <- vent %>% select(hospitalization_id, index_dttm) %>%
  filter(!hospitalization_id %in% at_index$hospitalization_id) %>%
  inner_join(icu, by = "hospitalization_id", relationship = "many-to-many") %>%
  filter(in_dttm > index_dttm, in_dttm <= index_dttm + ICU_AFTER_H * 3600) %>%
  arrange(hospitalization_id, in_dttm) %>% distinct(hospitalization_id, .keep_all = TRUE) %>%
  mutate(unit_source = paste0("first ICU within ", ICU_AFTER_H, " h after the index"))
placed <- bind_rows(at_index, after) %>% select(hospitalization_id, hospital_id, unit, unit_source)

d <- vent %>% left_join(placed, by = "hospitalization_id") %>%
  mutate(unit_source = coalesce(unit_source, "no ICU at or within 24 h of the index"),
         year = as.integer(format(index_dttm, "%Y", tz = "UTC")),
         period = paste0(year - (year - min(year, na.rm = TRUE)) %% CELL_YEARS,
                         if (CELL_YEARS > 1) paste0("-", year - (year - min(year, na.rm = TRUE)) %% CELL_YEARS + CELL_YEARS - 1) else ""),
         cell = paste(hospital_id, unit, period, sep = " | "),
         lvt = log(vtpbw * pbw), lvtpbw = log(vtpbw), lpbw = log(pbw), lpfvc = log(pfvc), lambda = log(pbw / pfvc),
         female = as.numeric(sex_category == "Female"), age10 = age_at_admission / 10,
         lbmi = if_else(bmi > 0, log(bmi), NA_real_))
cell_n <- d %>% filter(!is.na(unit)) %>% count(cell)
d <- d %>% left_join(cell_n %>% rename(cell_n = n), by = "cell")
placement <- d %>% mutate(status = case_when(is.na(unit) ~ unit_source,
                                             cell_n < CELL_MIN ~ paste0("in a cell of fewer than ", CELL_MIN),
                                             TRUE ~ paste0("analysed (", unit_source, ")"))) %>%
  count(status, name = "n_patients") %>% mutate(share = n_patients / sum(n_patients))
a <- d %>% filter(!is.na(unit), cell_n >= CELL_MIN)
if (n_distinct(a$cell) < 2) stop("fewer than two cells with ", CELL_MIN, " or more patients: nothing to compare")
message("Cells with ", CELL_MIN, " or more patients: ", n_distinct(a$cell), ", holding ", nrow(a), " patients")

# ---- 1. Per-cell rho and its heterogeneity --------------------------------------------
cells <- a %>% group_by(hospital_id, unit, period, cell) %>%
  group_modify(function(x, k) {
    f <- lm(lvt ~ lpbw, data = x); s <- summary(f)$coefficients
    tibble(n_patients = nrow(x), mean_vtpbw = exp(mean(x$lvtpbw)), sd_log_vtpbw = sd(x$lvtpbw),
           rho = s["lpbw", 1], rho_se = s["lpbw", 2], sd_log_pbw = sd(x$lpbw),
           mean_age = mean(x$age_at_admission), share_female = mean(x$female), mean_sofa = mean(x$sofa_total, na.rm = TRUE),
           mean_sf = mean(x$sf_ratio, na.rm = TRUE), mean_bmi = mean(x$bmi, na.rm = TRUE), mean_lambda = mean(x$lambda))
  }) %>% ungroup()
dl <- function(est, se) {                        # DerSimonian-Laird
  w <- 1 / se^2; m <- sum(w * est) / sum(w); Q <- sum(w * (est - m)^2); k <- length(est)
  tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  list(mean = m, tau = sqrt(tau2), Q = Q, df = k - 1, p = pchisq(Q, k - 1, lower.tail = FALSE), I2 = max(0, (Q - (k - 1)) / Q))
}
het <- dl(cells$rho, cells$rho_se)
cells <- cells %>% mutate(shrink = het$tau^2 / (het$tau^2 + rho_se^2), rho_shrunk = het$mean + shrink * (rho - het$mean))
# the same question at the patient level: one slope for every cell, or one per cell?
f_common <- lm(lvt ~ factor(cell) + lpbw, data = a)
f_cell   <- lm(lvt ~ factor(cell) * lpbw, data = a)
slopes_F <- anova(f_common, f_cell)

# ---- 2. The instrument: leave-one-out, shrunken rho_c x log PBW --------------------------
loo <- a %>% group_by(cell) %>%
  mutate(n = n(), sx = sum(lpbw), sy = sum(lvt), sxx = sum(lpbw^2), sxy = sum(lpbw * lvt),
         rho_loo = ((sxy - lpbw * lvt) - (sx - lpbw) * (sy - lvt) / (n - 1)) / ((sxx - lpbw^2) - (sx - lpbw)^2 / (n - 1))) %>%
  ungroup() %>% left_join(cells %>% select(cell, shrink), by = "cell") %>%
  # the instrument is the cell's slope applied to the patient's PBW about the cell's mean
  # PBW: on uncentred log PBW (about 4.1), the leave-one-out noise in rho times 4.1
  # would swamp the signal, rho_c times the spread of log PBW within the cell (about 0.17)
  group_by(cell) %>% mutate(lpbw_c = lpbw - mean(lpbw)) %>% ungroup() %>%
  mutate(rho_iv = het$mean + shrink * (rho_loo - het$mean), z = rho_iv * lpbw_c, z_raw = rho_loo * lpbw_c)

# coefficient on z given cell, log PBW and log PFVC, with the standard error clustered on the cell
cluster_fit <- function(y, data, instrument = "z") {
  ok <- is.finite(data[[y]]) & is.finite(data[[instrument]])
  dd <- data[ok, , drop = FALSE]
  dd$z <- dd[[instrument]]
  X <- model.matrix(~ factor(cell) + lpbw + lpfvc + z, data = dd)
  qx <- qr(X); X <- X[, qx$pivot[seq_len(qx$rank)], drop = FALSE]
  A <- solve(crossprod(X)); b <- A %*% crossprod(X, dd[[y]]); e <- dd[[y]] - X %*% b
  S <- rowsum(X * as.numeric(e), dd$cell); G <- nrow(S)
  V <- A %*% crossprod(S) %*% A * G / (G - 1)
  j <- which(colnames(X) == "z")
  # no between-cell variation in rho: the instrument is a multiple of log PBW and drops out
  if (length(j) != 1) return(tibble(outcome = y, coef = NA_real_, se = NA_real_, t = NA_real_, n_patients = nrow(dd), n_cells = G))
  tibble(outcome = y, coef = b[j, 1], se = sqrt(V[j, j]), t = b[j, 1] / sqrt(V[j, j]), n_patients = nrow(dd), n_cells = G)
}
first <- cluster_fit("lvt", loo) %>% mutate(F = t^2)
first_raw <- cluster_fit("lvt", loo, "z_raw") %>% mutate(F = t^2)    # without shrinkage: noisier, for comparison
# how much the instrument moves beside log PBW and log PFVC (its residual spread sets the precision)
z_res <- resid(lm(z ~ factor(cell) + lpbw + lpfvc, data = loo, na.action = na.exclude))
lam_res <- resid(lm(lambda ~ factor(cell) + age10 + female, data = loo))

# ---- 3. Balance -------------------------------------------------------------------------
loo <- loo %>% mutate(sofa = sofa_total, sf = sf_ratio)
balance <- map_dfr(c("sofa", "sf", "age10", "lbmi", "female"), cluster_fit, data = loo) %>%
  mutate(p_value = 2 * pnorm(-abs(t)), outcome = recode(outcome, sofa = "SOFA", sf = "S/F", age10 = "age (decades)",
                                                        lbmi = "log BMI", female = "female"))
case_mix <- map_dfr(c("mean_age", "share_female", "mean_sofa", "mean_sf", "mean_bmi", "mean_lambda"), function(v) {
  ok <- is.finite(cells[[v]])
  if (sum(ok) < 4 || sd(cells$rho_shrunk[ok]) == 0 || sd(cells[[v]][ok]) == 0)
    return(tibble(cell_characteristic = v, correlation_with_rho = NA_real_, p_value = NA_real_, n_cells = sum(ok)))
  ct <- cor.test(cells$rho_shrunk[ok], cells[[v]][ok])
  tibble(cell_characteristic = v, correlation_with_rho = unname(ct$estimate), p_value = ct$p.value, n_cells = sum(ok))
})
# a period effect: rho over calendar time, pooled across units
by_period <- a %>% group_by(period) %>%
  group_modify(~ { s <- summary(lm(lvt ~ factor(cell) + lpbw, data = .x))$coefficients
                   tibble(n_patients = nrow(.x), n_cells = n_distinct(.x$cell), rho = s["lpbw", 1], rho_se = s["lpbw", 2]) }) %>% ungroup()

summary_tbl <- bind_rows(
  tibble(quantity = "cells analysed", value = nrow(cells)),
  tibble(quantity = "patients analysed", value = nrow(a)),
  tibble(quantity = "rho, precision-weighted mean over cells", value = het$mean),
  tibble(quantity = "rho, between-cell SD (tau)", value = het$tau),
  tibble(quantity = "rho, I^2", value = het$I2),
  tibble(quantity = "rho, Q test p (all cells equal)", value = het$p),
  tibble(quantity = "rho shrunken, 10th centile over cells", value = unname(quantile(cells$rho_shrunk, 0.1))),
  tibble(quantity = "rho shrunken, 90th centile over cells", value = unname(quantile(cells$rho_shrunk, 0.9))),
  tibble(quantity = "slopes differ by cell: F", value = slopes_F$F[2]),
  tibble(quantity = "slopes differ by cell: p", value = slopes_F$`Pr(>F)`[2]),
  tibble(quantity = "first stage: coefficient of rho_c x log PBW on log VT", value = first$coef),
  tibble(quantity = "first stage: F (cluster-robust, 1 df)", value = first$F),
  tibble(quantity = "first stage without shrinkage: F", value = first_raw$F),
  tibble(quantity = "instrument residual SD (beside cell, log PBW, log PFVC)", value = sd(z_res, na.rm = TRUE)),
  tibble(quantity = "lambda residual SD (beside cell, age, sex), for comparison", value = sd(lam_res))) %>%
  mutate(site = site_name, cohort = cohort, cell_unit = CELL_UNIT, cell_years = CELL_YEARS, cell_min = CELL_MIN, .before = 1)

out_file <- function(name, ext = "csv") file.path(final_dir, paste0("causal_practice_", name, site_name, ".", ext))
tag <- function(x) x %>% mutate(site = site_name, cohort = cohort, .before = 1)
write_csv(tag(cells), out_file("cells_"))
write_csv(bind_rows(summary_tbl,
                    tag(balance %>% transmute(quantity = paste0("balance: instrument coefficient on ", outcome), value = coef,
                                              se, p_value, n_patients, n_cells)),
                    tag(case_mix %>% transmute(quantity = paste0("cell-level correlation of rho with ", cell_characteristic),
                                               value = correlation_with_rho, p_value, n_cells)),
                    tag(by_period %>% transmute(quantity = paste0("rho in period ", period), value = rho, se = rho_se,
                                                n_patients, n_cells))), out_file("summary_"))
write_csv(tag(placement), out_file("placement_"))

pdf(out_file("", "pdf"), width = 10, height = 7)
multi <- cells %>% group_by(hospital_id, unit) %>% filter(n() > 1) %>% ungroup()   # lines only where a unit has several periods
print(ggplot(cells, aes(period, rho, colour = unit, group = unit)) +
        geom_hline(yintercept = c(0, 1), linetype = 2, colour = "grey60") +
        geom_pointrange(aes(ymin = rho - 1.96 * rho_se, ymax = rho + 1.96 * rho_se), position = position_dodge(width = 0.5), size = 0.25) +
        geom_line(data = multi, position = position_dodge(width = 0.5), alpha = 0.5) + facet_wrap(~ hospital_id) +
        labs(x = "Period", y = "rho: slope of log VT on log PBW, 95% CI",
             title = paste0("How tidal volume followed PBW, by ICU and period: ", site_name, " (", cohort, ")"),
             subtitle = "1: mL/kg by PBW; 0: the same volume for everyone") +
        theme_bw(base_size = 9) + theme(legend.position = "bottom"))
print(ggplot(cells %>% arrange(rho_shrunk) %>% mutate(k = row_number()), aes(k, rho_shrunk)) +
        geom_hline(yintercept = het$mean, linetype = 2) + geom_point(aes(size = n_patients), alpha = 0.6) +
        geom_errorbar(aes(ymin = rho - 1.96 * rho_se, ymax = rho + 1.96 * rho_se), alpha = 0.25, width = 0) +
        labs(x = "Cell, ranked", y = "rho, shrunken (points) and raw 95% CI (bars)",
             title = "Between-cell spread of rho", subtitle = sprintf("tau = %.2f, I^2 = %.0f%%, first-stage F = %.1f", het$tau, 100 * het$I2, first$F)) +
        theme_bw(base_size = 9))
invisible(dev.off())

message("\nHow patients were placed in cells:")
print(as.data.frame(placement %>% mutate(share = round(share, 3))), row.names = FALSE)
if (het$tau == 0) message("\nNo between-cell variation in rho beyond noise (tau = 0): the shrunken instrument is not defined.")
message("\nIs there variation, is it strong, is it plausibly exogenous?")
print(as.data.frame(summary_tbl %>% transmute(quantity, value = signif(value, 3))), row.names = FALSE)
message("\nBalance (instrument coefficient on baseline characteristics; a usable design has these near 0):")
print(as.data.frame(balance %>% transmute(outcome, coef = signif(coef, 3), se = signif(se, 3), p_value = signif(p_value, 2))), row.names = FALSE)
message("\nDo cells with a higher rho have a different case mix?")
print(as.data.frame(case_mix %>% mutate(across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)
message("\nrho by period (pooled across units):")
print(as.data.frame(by_period %>% mutate(across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)
message("\ncausal_practice_cells complete: ", final_dir)
