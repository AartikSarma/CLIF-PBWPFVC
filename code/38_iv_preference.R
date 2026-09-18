# =============================================================================
# Script 38: practice-variation instrument for the strain-limiting strategy
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# The target-trial emulation (the TTE (30_tte_common)) estimates the 60-day mortality effect of
# a strain-limiting strategy (VT/PFVC held under C_LOW) from measured confounders
# and is exposed to confounding by indication and by quality of care. This is
# the second, independent design: an instrument that shifts who receives the
# strategy for reasons unrelated to the patient. The instrument is practice
# variation (Brookhart's preference instrument): the strain-limiting rate among
# the OTHER cohort patients admitted to the same ICU in the same calendar period.
# A patient admitted to a unit in a season when that unit was holding strain low
# is more likely to be held low, whatever their own severity.
#
#   exposure   A = 1 if the patient's VT/PFVC stayed <= C_LOW on every panel day
#                  after the grace day through K_DAYS (the early strategy); the
#                  continuous companion D = mean VT/PFVC over those days
#   instrument Z = the strain-limiting rate of the patient's ICU in the ADJACENT
#                  half-years (the periods before and after the patient's own), a
#                  rate the patient is not part of; companions: the leave-one-out
#                  rate of the patient's own ICU x period cell (mechanically biased
#                  toward zero and below in small cells once fixed effects absorb
#                  the between-cell variation, kept for the record), the ICU's
#                  leave-one-out rate over all time (for sites whose dates are
#                  shifted per patient, as MIMIC's are, where calendar cells do not
#                  group contemporaries), and the site x period rate (single-unit
#                  sites). Cells need MIN_CELL patients.
#   outcome    Y = 60-day all-cause mortality
#   covariates V = ns(age, 4) + sex + race + SOFA + log index SF + log PBW/PFVC,
#                  with ICU and period fixed effects, so the instrument is the
#                  ICU x period interaction: how far a unit's practice in that
#                  season sat from its own average and from the other units'
#
# The DOSE instrument (the read at UCSF, 2026-09-17). Whether a patient is
# strain-limited is decided by their lung size, not by their clinicians: within
# 6-8 mL/kg PBW dosing, VT/PFVC <= 11% is the concordant-lung patient, and the
# strain-limiting rate of a unit in neighbouring seasons predicted nothing about
# the next patient (first stage 0, F 0.002). What clinicians do vary, between
# units and seasons, is the dose itself. So the second block instruments the
# continuous strain (mean VT/PFVC, D) with the unit's dose preference: the mean
# VT/PBW among the OTHER patients (adjacent periods; the unit over all time; the
# site x period), and reports the mortality change per point of VT/PFVC and the
# implied effect of the strain-limiting policy = coefficient x the mean shift
# needed to bring every patient under the ceiling (iv_height_policy's construction).
#
# Estimators: the naive adjusted risk difference (linear probability, HC1), the
# 2SLS risk difference per switch to strain-limiting (the complier effect), the
# reduced form, and the first stage with its F statistic. Instrument checks: the
# Brookhart balance table (standardized differences of the covariates across the
# exposure and across the instrument), and a falsification regression of every
# pre-index covariate on the instrument (the instrument must not predict who
# arrives sicker). Exclusion (the unit's season affects mortality only through
# strain) is an assumption; unit and period fixed effects narrow it to the
# interaction.
#
# Outputs (aggregates only; cells with fewer than MIN_CELL patients are dropped):
#   final/iv_preference_{site}.csv          the estimates
#   final/iv_preference_balance_{site}.csv  the Brookhart table and falsification
#   final/iv_preference_cells_{site}.csv    ICU x period cells: n, strain-limiting rate
#   final/iv_preference_{site}.pdf          forest, first stage, cell rates over time
# Usage: Rscript code/38_iv_preference.R   (PBWPFVC_TTE_CLOW=11, PBWPFVC_IV_DAYS=3,
#        PBWPFVC_IV_PERIOD_MONTHS=6, PBWPFVC_IV_MIN_CELL=20)
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(tidyverse); library(arrow); library(here); library(splines)
  library(lubridate); library(sandwich); library(lmtest); library(patchwork)
})
rm(list = ls())
source("utils/config.R")
site_name  <- config$site_name
output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
okabe <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#000000")   # the full Okabe-Ito set: up to eight instruments

C_LOW     <- as.numeric(Sys.getenv("PBWPFVC_TTE_CLOW", "11"))        # the strain ceiling of the TTE (30_tte_common)
GRACE     <- 1L
K_DAYS    <- as.integer(Sys.getenv("PBWPFVC_IV_DAYS", "3"))         # the early strategy: days GRACE+1 .. K_DAYS
PERIOD_M  <- as.integer(Sys.getenv("PBWPFVC_IV_PERIOD_MONTHS", "6"))
MIN_CELL  <- as.integer(Sys.getenv("PBWPFVC_IV_MIN_CELL", "20"))    # other patients a cell needs to yield an instrument
Y_DAYS    <- 60

# the shared daily panel (10_panel_common.R): base (index t0, size, severity, death)
# and the daily VT/PFVC panel
HORIZON      <- 28L
MAX_VENT_DAY <- 27L
is_synthetic <- grepl("^synthetic_clif", site_name)
PANEL_NORM   <- "pfvc"
source(here("code", "10_panel_common.R"))
tables_path <- path.expand(config$tables_path); file_type <- config$file_type
open_clif <- function(tbl) {
  fpath <- file.path(tables_path, paste0("clif_", tbl, ".", file_type))
  if (file_type == "parquet") return(arrow::open_dataset(fpath))
  if (file_type == "csv")     return(readr::read_csv(fpath, show_col_types = FALSE))
  if (file_type == "fst")     return(fst::read_fst(fpath))
  stop("Unsupported file_type: ", file_type)
}

# ---- exposure: the early strategy from the daily panel
strat <- panel %>% filter(vent_day > GRACE, vent_day <= K_DAYS, is.finite(vtpfvc)) %>%
  group_by(hospitalization_id) %>%
  summarise(n_days = n(), max_vtpfvc = max(vtpfvc), mean_vtpfvc = mean(vtpfvc),
            mean_vtpbw = mean(vt_ml / pbw), .groups = "drop") %>%
  mutate(A = as.integer(max_vtpfvc <= C_LOW))

# ---- outcome: 60-day all-cause mortality from the index (synthetic: the simulated survival)
cs <- read_parquet(file.path(output_dir, "analysis_cross_sectional.parquet")) %>%
  select(hospitalization_id, death_dttm)
died <- base %>% select(hospitalization_id, t0, death_time_days) %>% left_join(cs, by = "hospitalization_id") %>%
  mutate(Y = if (is_synthetic) as.integer(!is.na(death_time_days)) else
           as.integer(!is.na(death_dttm) & as.numeric(difftime(death_dttm, t0, units = "days")) <= Y_DAYS &
                        as.numeric(difftime(death_dttm, t0, units = "days")) >= 0)) %>%
  select(hospitalization_id, Y)

# ---- the ICU of the index: the ICU stay covering t0, else the first ICU stay after it, else the first
adt <- open_clif("adt") %>% filter(hospitalization_id %in% base$hospitalization_id) %>% collect() %>%
  filter(tolower(location_category) == "icu") %>%
  transmute(hospitalization_id, icu = coalesce(as.character(location_name), "ICU"), in_dttm, out_dttm) %>%
  inner_join(base %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(covers = !is.na(in_dttm) & in_dttm <= t0 & (is.na(out_dttm) | out_dttm > t0),
         after  = !is.na(in_dttm) & in_dttm > t0) %>%
  group_by(hospitalization_id) %>% arrange(desc(covers), !after, in_dttm, .by_group = TRUE) %>%
  slice(1) %>% ungroup() %>% select(hospitalization_id, icu)

# ---- the analysis frame
sf0 <- sf_daily %>% filter(vent_day == 0L) %>% select(hospitalization_id, sf_0 = sf)
d <- base %>%
  inner_join(strat, by = "hospitalization_id") %>%
  left_join(died, by = "hospitalization_id") %>%
  left_join(adt, by = "hospitalization_id") %>%
  left_join(sf0, by = "hospitalization_id") %>%
  mutate(icu = coalesce(icu, "ICU"),
         period = paste0(year(t0), "-", sprintf("%02d", (month(t0) - 1) %/% PERIOD_M + 1)),
         cell = paste(icu, period, sep = " | "),
         log_sf_0 = log(sf_0), ldisc_c = log(pbw / pfvc_gli) - median(log(pbw / pfvc_gli), na.rm = TRUE)) %>%
  filter(!is.na(Y), !is.na(sofa_total), is.finite(log_sf_0), !is.na(age10), !is.na(sex_category), !is.na(race_category))
message("Cohort with an early-strategy assessment: ", nrow(d), " patients; strain-limiting ", sum(d$A), " (",
        round(100 * mean(d$A)), "%); deaths ", sum(d$Y), "; ICUs ", n_distinct(d$icu), "; periods ", n_distinct(d$period))

# ---- the instruments
# is the calendar real? per-patient date shifting (MIMIC) spreads admissions over a
# century and makes calendar cells meaningless; flag it and prefer the unit instrument
calendar_ok <- diff(range(year(d$t0), na.rm = TRUE)) <= 30 && n_distinct(d$period) <= 60
if (!calendar_ok) message("NOTE: admission dates span ", diff(range(year(d$t0), na.rm = TRUE)), " years over ",
                          n_distinct(d$period), " periods: per-patient date shifting; period-based instruments are not meaningful here")
d <- d %>% mutate(period_idx = as.integer(factor(period, levels = sort(unique(period)))))
cell_sums <- d %>% group_by(icu, period_idx) %>% summarise(s = sum(A), n = n(), .groups = "drop")
adj <- cell_sums %>% select(icu, period_idx, s, n) %>%
  mutate(period_idx = period_idx - 1L) %>% rename(s_next = s, n_next = n) %>%   # the next period, keyed to this one
  full_join(cell_sums %>% mutate(period_idx = period_idx + 1L) %>% rename(s_prev = s, n_prev = n), by = c("icu", "period_idx")) %>%
  mutate(s_adj = coalesce(s_prev, 0L) + coalesce(s_next, 0L), n_adj = coalesce(n_prev, 0L) + coalesce(n_next, 0L)) %>%
  select(icu, period_idx, s_adj, n_adj)
d <- d %>% left_join(adj, by = c("icu", "period_idx")) %>%
  mutate(z_adj = if_else(!is.na(n_adj) & n_adj >= MIN_CELL, s_adj / n_adj, NA_real_)) %>%      # PRIMARY: the unit's neighbouring seasons
  group_by(cell) %>%
  mutate(n_cell = n(), z_icu = if_else(n_cell - 1L >= MIN_CELL, (sum(A) - A) / (n_cell - 1L), NA_real_)) %>%   # own cell, leave-one-out
  ungroup() %>%
  group_by(icu) %>%
  mutate(n_icu = n(), z_unit = if_else(n_icu - 1L >= MIN_CELL, (sum(A) - A) / (n_icu - 1L), NA_real_)) %>%     # the unit over all time
  ungroup() %>%
  group_by(period) %>%
  mutate(n_period = n(), z_site = if_else(n_period - 1L >= MIN_CELL, (sum(A) - A) / (n_period - 1L), NA_real_)) %>%   # site x period
  ungroup()
# the dose instruments: the unit's mean VT/PBW (mL/kg) among the other patients
dose_sums <- d %>% group_by(icu, period_idx) %>% summarise(sd_ = sum(mean_vtpbw), n = n(), .groups = "drop")
dose_adj <- dose_sums %>% mutate(period_idx = period_idx - 1L) %>% rename(s_next = sd_, n_next = n) %>%
  full_join(dose_sums %>% mutate(period_idx = period_idx + 1L) %>% rename(s_prev = sd_, n_prev = n), by = c("icu", "period_idx")) %>%
  transmute(icu, period_idx, sd_adj = coalesce(s_prev, 0) + coalesce(s_next, 0), nd_adj = coalesce(n_prev, 0L) + coalesce(n_next, 0L))
d <- d %>% left_join(dose_adj, by = c("icu", "period_idx")) %>%
  mutate(zd_adj = if_else(!is.na(nd_adj) & nd_adj >= MIN_CELL, sd_adj / nd_adj, NA_real_)) %>%
  group_by(icu)    %>% mutate(nn = n(), zd_unit = if_else(nn - 1L >= MIN_CELL, (sum(mean_vtpbw) - mean_vtpbw) / (nn - 1L), NA_real_)) %>% ungroup() %>%
  group_by(period) %>% mutate(nn = n(), zd_site = if_else(nn - 1L >= MIN_CELL, (sum(mean_vtpbw) - mean_vtpbw) / (nn - 1L), NA_real_)) %>% ungroup() %>%
  select(-nn)
cells <- d %>% group_by(icu, period) %>% summarise(n = n(), rate = mean(A), deaths = sum(Y), .groups = "drop") %>%
  filter(n >= 10) %>% mutate(site = site_name)
write_csv(cells, file.path(final_dir, paste0("iv_preference_cells_", site_name, ".csv")))
message("Cells (ICU x ", PERIOD_M, "-month period) with >= 10 patients: ", nrow(cells), "; strain-limiting rate range ",
        signif(min(cells$rate), 2), " to ", signif(max(cells$rate), 2))

# ---- estimators
V_RHS <- "ns(age10, 4) + sex_category + race_category + sofa_total + log_sf_0 + ldisc_c"
hc <- function(m, term) { V <- sandwich::vcovHC(m, type = "HC1"); b <- coef(m)[term]; se <- sqrt(V[term, term])
  tibble(estimate = unname(b), se = unname(se), lo = unname(b - 1.96 * se), hi = unname(b + 1.96 * se),
         p = unname(2 * pnorm(-abs(b / se)))) }
# two-stage least squares by hand (first-stage fitted values, second-stage residuals from the
# actual exposure), HC1 on the correct residuals, so no package beyond sandwich is needed
tsls <- function(dd, y, x, z, rhs) {
  fs  <- lm(as.formula(paste(x, "~", z, "+", rhs)), data = dd)
  dd$x_hat <- fitted(fs)
  ss  <- lm(as.formula(paste(y, "~ x_hat +", rhs)), data = dd)
  X   <- model.matrix(ss); Xa <- X; Xa[, "x_hat"] <- dd[[x]]   # actual exposure in place of the fitted
  res <- dd[[y]] - as.numeric(Xa %*% coef(ss))
  bread <- solve(crossprod(X)); meat <- crossprod(X * res)
  Vh <- bread %*% meat %*% bread * nrow(X) / (nrow(X) - ncol(X))
  b <- coef(ss)["x_hat"]; se <- sqrt(Vh["x_hat", "x_hat"])
  Fz <- (coef(fs)[z] / sqrt(sandwich::vcovHC(fs, type = "HC1")[z, z]))^2
  list(est = tibble(estimate = unname(b), se = unname(se), lo = unname(b - 1.96 * se), hi = unname(b + 1.96 * se),
                    p = unname(2 * pnorm(-abs(b / se)))), first_stage = hc(fs, z), F = unname(Fz))
}
# fixed effects that exist in the frame: ICU when more than one unit, period only
# for the ICU x period instrument (a period effect would absorb the site-wide one)
fe_terms <- function(dd, with_period, unit_fe = TRUE) paste(c(if (unit_fe && n_distinct(dd$icu) > 1) "+ factor(icu)",
                                             if (with_period && n_distinct(dd$period) > 1) "+ factor(period)"), collapse = " ")
run_iv <- function(z, label, with_period, unit_fe = TRUE) {
  dd <- d %>% filter(is.finite(.data[[z]]))
  if (nrow(dd) < 100 || sd(dd[[z]]) == 0) {
    message("instrument '", label, "' unusable here: ", nrow(dd), " patients in cells with >= ", MIN_CELL,
            " others (a small site or too many units for the period; lower PBWPFVC_IV_MIN_CELL or lengthen PBWPFVC_IV_PERIOD_MONTHS)")
    return(NULL)
  }
  rhs <- paste(V_RHS, fe_terms(dd, with_period, unit_fe))
  naive <- hc(lm(as.formula(paste("Y ~ A +", rhs)), data = dd), "A")
  rf    <- hc(lm(as.formula(paste("Y ~", z, "+", rhs)), data = dd), z)
  iv_a  <- tsls(dd, "Y", "A", z, rhs)
  iv_d  <- tsls(dd, "Y", "mean_vtpfvc", z, rhs)
  bind_rows(
    naive %>% mutate(estimator = "naive adjusted RD", exposure = "strain-limiting (A)"),
    iv_a$est %>% mutate(estimator = "2SLS RD per switch (complier)", exposure = "strain-limiting (A)", first_stage_F = iv_a$F,
                        first_stage_coef = iv_a$first_stage$estimate),
    rf %>% mutate(estimator = "reduced form (Y on Z)", exposure = "instrument"),
    iv_a$first_stage %>% mutate(estimator = "first stage (A on Z)", exposure = "instrument", first_stage_F = iv_a$F),
    iv_d$est %>% mutate(estimator = "2SLS RD per point of VT/PFVC", exposure = "mean VT/PFVC (D)", first_stage_F = iv_d$F)) %>%
    mutate(instrument = label, n = nrow(dd), n_deaths = sum(dd$Y), n_strain_limited = sum(dd$A),
           n_icus = n_distinct(dd$icu), n_periods = n_distinct(dd$period), .before = 1)
}
# the dose block: D (mean VT/PFVC, per point) instrumented by the unit's dose preference;
# the policy effect = coefficient x mean shift to the ceiling among those above it
run_iv_dose <- function(z, label, with_period, unit_fe = TRUE) {
  dd <- d %>% filter(is.finite(.data[[z]]), is.finite(mean_vtpfvc))
  if (nrow(dd) < 100 || sd(dd[[z]]) == 0) { message("dose instrument '", label, "' unusable here (", nrow(dd), " patients)"); return(NULL) }
  rhs <- paste(V_RHS, fe_terms(dd, with_period, unit_fe))
  naive <- hc(lm(as.formula(paste("Y ~ mean_vtpfvc +", rhs)), data = dd), "mean_vtpfvc")
  iv    <- tsls(dd, "Y", "mean_vtpfvc", z, rhs)
  shift <- mean(pmax(dd$mean_vtpfvc - C_LOW, 0))          # mean points above the ceiling (0 for those under it)
  bind_rows(
    naive %>% mutate(estimator = "naive adjusted RD per point of VT/PFVC", exposure = "mean VT/PFVC (D)"),
    iv$est %>% mutate(estimator = "2SLS RD per point of VT/PFVC (dose instrument)", exposure = "mean VT/PFVC (D)",
                      first_stage_F = iv$F, first_stage_coef = iv$first_stage$estimate),
    iv$first_stage %>% mutate(estimator = "first stage (D on unit dose preference)", exposure = "dose instrument", first_stage_F = iv$F),
    # the policy LOWERS strain by `shift`, so its risk difference is minus the per-point effect times the shift
    iv$est %>% mutate(estimate = -estimate * shift, se = se * shift, lo = -hi * shift, hi = -lo * shift,
                      estimator = sprintf("2SLS policy RD: everyone to <= %g%% (mean shift %.2f points)", C_LOW, shift),
                      exposure = "strain-limiting policy", first_stage_F = iv$F)) %>%
    mutate(instrument = label, n = nrow(dd), n_deaths = sum(dd$Y), n_strain_limited = sum(dd$A),
           n_icus = n_distinct(dd$icu), n_periods = n_distinct(dd$period), mean_shift_points = shift, .before = 1)
}
res <- bind_rows(
  run_iv_dose("zd_adj",  "unit dose preference, adjacent periods (primary)", with_period = TRUE),
  run_iv_dose("zd_unit", "unit dose preference, all periods (no unit fixed effect)", with_period = calendar_ok, unit_fe = FALSE),
  run_iv_dose("zd_site", "site dose preference by period", with_period = FALSE),
  run_iv("z_adj",  "ICU rate in the adjacent periods (primary)", with_period = TRUE),
  run_iv("z_icu",  "ICU x period leave-one-out rate (mechanically biased in small cells)", with_period = TRUE),
  run_iv("z_unit", "ICU leave-one-out rate, all periods (no unit fixed effect)", with_period = calendar_ok, unit_fe = FALSE),
  run_iv("z_site", "site x period leave-one-out rate", with_period = FALSE)) %>%
  mutate(c_low = C_LOW, k_days = K_DAYS, period_months = PERIOD_M, calendar_reliable = calendar_ok, site = site_name)
write_csv(res, file.path(final_dir, paste0("iv_preference_", site_name, ".csv")))
message("\nDose instruments (60-day mortality per point of VT/PFVC; the policy row = per-point effect x the mean shift to the ceiling):")
print(as.data.frame(res %>% filter(grepl("dose", instrument)) %>% select(instrument, estimator, estimate, lo, hi, p, first_stage_F, n) %>%
                      mutate(instrument = substr(instrument, 1, 30), across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
message("\nStrategy instruments (RD per switch to strain-limiting):")
print(as.data.frame(res %>% filter(!grepl("dose", instrument)) %>% select(instrument, estimator, estimate, lo, hi, p, first_stage_F, n) %>%
                      mutate(instrument = substr(instrument, 1, 30), across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)

# ---- instrument checks: Brookhart balance and falsification
smd <- function(x, g) { m1 <- mean(x[g == 1], na.rm = TRUE); m0 <- mean(x[g == 0], na.rm = TRUE)
  s <- sqrt((var(x[g == 1], na.rm = TRUE) + var(x[g == 0], na.rm = TRUE)) / 2); (m1 - m0) / s }
# on the primary instrument when it has the patients, else the unit one, else the site one
z_use <- if (calendar_ok && sum(is.finite(d$z_adj)) >= 100) "z_adj" else if (sum(is.finite(d$z_unit)) >= 100) "z_unit" else "z_site"
dz <- d %>% filter(is.finite(.data[[z_use]])) %>% mutate(z = .data[[z_use]], z_hi = as.integer(z > median(z)),
                                                        female = as.integer(sex_category == "Female"))
bal_fe <- fe_terms(dz, with_period = z_use == "z_adj" || (z_use == "z_unit" && calendar_ok), unit_fe = z_use != "z_unit")
covs <- c(age10 = "age (decades)", female = "female", sofa_total = "SOFA", log_sf_0 = "log index SF",
          ldisc_c = "log PBW/PFVC", bmi = "BMI")
balance <- imap_dfr(covs, function(lab, v) {
  f <- lm(as.formula(paste(v, "~ z", bal_fe)), data = dz)
  co <- coeftest(f, vcov = sandwich::vcovHC(f, type = "HC1"))["z", ]
  tibble(covariate = lab, smd_across_exposure = smd(dz[[v]], dz$A), smd_across_instrument = smd(dz[[v]], dz$z_hi),
         falsification_coef = unname(co[1]), falsification_p = unname(co[4]))
}) %>% mutate(instrument = z_use, n = nrow(dz), site = site_name)
write_csv(balance, file.path(final_dir, paste0("iv_preference_balance_", site_name, ".csv")))
message("\nBrookhart balance (standardized differences) and falsification (covariate on the instrument ", z_use, ", within the fixed effects):")
print(as.data.frame(balance %>% mutate(across(where(is.numeric), ~ signif(., 2)))), row.names = FALSE)

# ---- figure
fr <- res %>% filter(estimator %in% c("naive adjusted RD", "2SLS RD per switch (complier)") | grepl("policy RD", estimator)) %>%
  mutate(estimator = if_else(grepl("policy RD", estimator), "2SLS policy RD (dose instrument)", estimator),
         estimator = factor(estimator, c("naive adjusted RD", "2SLS RD per switch (complier)", "2SLS policy RD (dose instrument)")),
         instrument = str_wrap(instrument, 28))
p1 <- ggplot(fr, aes(100 * estimate, estimator, colour = instrument)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = 100 * lo, xmax = 100 * hi), position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = okabe, name = NULL) +
  coord_cartesian(xlim = c(-60, 60)) +
  labs(title = "60-day mortality risk difference, strain-limiting vs not", x = "percentage points (95% CI; axis clipped at +/- 60)", y = NULL) +
  theme(legend.position = "bottom")
p2 <- ggplot(cells, aes(rate, deaths / n, size = n)) +
  geom_point(alpha = 0.6, colour = okabe[1]) + geom_smooth(method = "lm", se = TRUE, colour = okabe[4], linewidth = 0.7) +
  scale_size_continuous(name = "patients") +
  labs(title = "Cell strain-limiting rate against cell mortality (the reduced form, unadjusted)",
       x = "strain-limiting rate in the ICU x period cell", y = "60-day mortality in the cell")
p3 <- ggplot(cells, aes(period, rate, group = icu, colour = icu)) +
  geom_line(alpha = 0.7) + geom_point(aes(size = n), alpha = 0.7) +
  scale_colour_viridis_d(guide = "none") + scale_size_continuous(guide = "none") +
  labs(title = "Strain-limiting rate by ICU over time (the instrument's variation)", x = "period", y = "rate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p <- (p1 + p2) / p3 + plot_annotation(title = sprintf("Practice-variation instrument (%s)", site_name),
                                       subtitle = sprintf("A = VT/PFVC <= %g%% on days %d-%d; Z = leave-one-out rate of the ICU x %d-month cell; V = age, sex, race, SOFA, index SF, PBW/PFVC; ICU and period fixed effects",
                                                          C_LOW, GRACE + 1, K_DAYS, PERIOD_M)) &
  theme_minimal(base_size = 10)
ggsave(file.path(final_dir, paste0("iv_preference_", site_name, ".pdf")), p, width = 12, height = 9)
message("38_iv_preference complete -> ", final_dir)
