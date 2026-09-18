# =============================================================================
# Script 26 (quick): the longitudinal submodel alone, no joint model
# =============================================================================
# Fits the PFVC-level question as a plain linear mixed model on the panel
# tables 21_biotrauma_panel.R wrote, in seconds:
#
#   log marker ~ time + size exposure (per SD) + size exposure x time
#                + dose (patient mean, within-patient change)
#                + baseline marker + lagged SF and pressor + non-respiratory
#                  SOFA + BMI  [+ ns(age, 4) + sex + race]
#   random intercept and slope per patient (pdDiag)
#
# for BOTH size exposures (log PFVC per SD, and log PBW/PFVC per SD = VT/PFVC at
# a given VT/PBW), fitted SEPARATELY for the 24, 48 and 72-hour windows: each
# model uses only the rows up to its horizon and only patients whose baseline
# was observed before it, and reports the marker difference per SD at that
# horizon (level + divergence x time). This is the joint model's longitudinal
# part without the survival linkage, so death before the horizon is NOT
# accounted for: it is the fast look, and the joint-model contrast is the read.
#
# Also written, per window:
#   quick_channels_{marker}_{site}.csv       the channel decomposition of each
#       exposure (20_biotrauma_grid.R): the contrast at the horizon identified
#       through height, age, sex and race separately, with a test that they agree
#   quick_nested_{marker}_{site}.csv         nested ladder (ML fits): demographics
#       only, size only, both, pieces free, pieces + free demographics; LR and AIC
#   quick_dose_channels_{marker}_{site}.csv  dose x piece (between-patient mean
#       VT/PBW, centred): does each piece's contrast scale with the dose?
#   quick_sf_channels_{marker}_{site}.csv    severity x piece (baseline SF, the
#       panel's index-day worst SF, centred); skipped for sf
#
# Needs the 72-hour panel: PBWPFVC_JM_HORIZON_H=72 Rscript code/21_biotrauma_panel.R
# (written beside the 48-hour files, not over them).
#
# Usage: PBWPFVC_INJ_MARKER=creatinine Rscript code/26_quick_lme.R
#        (PBWPFVC_QUICK_HORIZONS_H=24,48,72; PBWPFVC_JM_HORIZON_H picks the panel, default 72)
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(here); library(splines); library(nlme) })
rm(list = ls())
source("utils/config.R")
site_name  <- config$site_name
output_dir <- config$output_dir
final_dir  <- final_dir_for("injury")
if (!nzchar(Sys.getenv("PBWPFVC_JM_HORIZON_H"))) Sys.setenv(PBWPFVC_JM_HORIZON_H = "72")
source(here("code", "20_biotrauma_grid.R"))
QUICK_HOURS <- as.numeric(strsplit(Sys.getenv("PBWPFVC_QUICK_HORIZONS_H", "24,48,72"), ",")[[1]])
QUICK_HOURS <- QUICK_HOURS[QUICK_HOURS <= JM_HORIZON * 24]

MARKER <- Sys.getenv("PBWPFVC_INJ_MARKER", "creatinine")
# the never-intubated control has no ventilator dose: dose terms and the dose x piece block are dropped
HAS_DOSE <- config$cohort == "imv"
if (!HAS_DOSE && MARKER == "dp") stop("driving pressure does not exist outside the ventilated cohort")
EXPOS  <- c(log_pfvc_sd = "log_pfvc", ldisc_sd = "ldisc")   # exposure column -> channel base
# the size contrasts: the two above plus VT/PFVC in percent of predicted FVC (patient
# mean over the window, centred, per point) at a given VT/PBW, the reader's form
SIZE_EXPOS <- c(names(EXPOS), if (HAS_DOSE) "vtpfvc_c")
y_col  <- c(creatinine = "creatinine", ne_equiv = "ne_equiv_peak", platelets = "platelets",
            bilirubin = "bilirubin", sf = "sf", dp = "dp", oi = "oi", osi = "osi")[[MARKER]]
y0_col <- c(creatinine = "creatinine_0", ne_equiv = "ne_equiv_0", platelets = "platelet_0",
            bilirubin = "bilirubin_0", sf = "sf_0", dp = "dp_0", oi = "oi_0", osi = "osi_0")[[MARKER]]
offset <- if (MARKER == "ne_equiv") 0.01 else 0
# SF is a component of OSI and of OI's ratio, so the lagged SF covariate is
# dropped for the oxygenation indices as it is for SF itself: adjusting an
# outcome for a piece of itself attenuates the exposure it is there to isolate.
own_lag <- c(sf = "l_log_sf", ne_equiv = "l_pressor", oi = "l_log_sf", osi = "l_log_sf")[MARKER]
SF_IS_OUTCOME <- MARKER %in% c("sf", "oi", "osi")
expo_label <- c(log_pfvc = "log PFVC", ldisc = "log PBW/PFVC (VT/PFVC at a given VT/PBW)")

panel_path <- file.path(output_dir, paste0("jm_long_", h_suffix, ".parquet"))
if (!file.exists(panel_path))
  stop("no ", h_suffix, " panel: run  PBWPFVC_JM_HORIZON_H=", as.integer(JM_HORIZON * 24), " Rscript code/21_biotrauma_panel.R")
long <- read_parquet(panel_path)
surv <- read_parquet(file.path(output_dir, paste0("jm_surv_", h_suffix, ".parquet")))
y0_day <- paste0(y0_col, "_day")
d_all <- long %>%
  filter(period >= 1L, !is.na(.data[[y_col]]), if (HAS_DOSE) !is.na(l_vtpbw_within) else TRUE, !is.na(l_sf), !is.na(l_pressor)) %>%
  inner_join(surv %>% select(hospitalization_id, np_sofa, bmi, age10, sex_category, race_category, height_cm, pfvc_gli,
                             sf_0, vtpbw_pt_mean, log_pfvc_sd, ldisc_sd, vtpfvc_c, all_of(c(y0_col, y0_day))), by = "hospitalization_id") %>%
  filter(!is.na(.data[[y0_col]]), !is.na(np_sofa), if (MARKER == "dp") !is.na(bmi) else TRUE) %>%
  mutate(log_y = log(.data[[y_col]] + offset), log_y0 = log(.data[[y0_col]] + offset), l_log_sf = log(l_sf))

lags <- setdiff(c("l_log_sf", "l_pressor"), own_lag)
ctrl <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)
# the common right-hand side: time, the size terms (level + divergence), the dose,
# the baseline, the lags, non-respiratory SOFA; BMI only for the pressure-derived marker
rhs_of <- function(size_terms, extra = NULL, dose_mean = "vtpbw_pt_mean")
  paste(c("vent_day", size_terms, if (length(size_terms)) paste0(size_terms, ":vent_day"), extra,
          if (HAS_DOSE) c("l_vtpbw_within", dose_mean),
          "log_y0", lags, "np_sofa", if (MARKER == "dp") "bmi"), collapse = " + ")
DEMO <- "ns(age10, 4) + sex_category + race_category"
window_data <- function(dat, hh) dat %>%
  filter(vent_day <= hh / 24, .data[[y0_day]] * STEP < hh / 24) %>%
  group_by(hospitalization_id) %>% filter(n() >= 2L) %>% ungroup() %>%
  mutate(id = factor(hospitalization_id))
fit_lme <- function(rhs, d, method = "REML")
  lme(as.formula(paste("log_y ~", rhs)), random = list(id = pdDiag(~ vent_day)), data = d, control = ctrl, method = method)
# a fixed-effect term by its components, whatever order R put them in
term_by <- function(b, comps) names(b)[sapply(strsplit(names(b), ":"), setequal, comps)]
# the contrast at the horizon for each of `terms`: level + divergence x (hh / 24)
contrast <- function(f, terms, hh) {
  b <- fixef(f); V <- as.matrix(vcov(f))
  W <- sapply(terms, function(tm) {
    tn <- term_by(b, c(strsplit(tm, ":")[[1]], "vent_day"))
    w <- setNames(rep(0, length(b)), names(b)); w[term_by(b, strsplit(tm, ":")[[1]])] <- 1; w[tn] <- hh / 24; w })
  list(est = as.numeric(t(W) %*% b), V = t(W) %*% V %*% W)
}

# =============================================================================
# 1. the size contrast per SD, per exposure and window, adjusted and unadjusted
# =============================================================================
out <- map_dfr(SIZE_EXPOS, function(EXPO) map_dfr(QUICK_HOURS, function(hh) {
  d <- window_data(d_all, hh) %>% filter(!is.na(.data[[EXPO]]))
  if (EXPO == names(EXPOS)[1]) message(MARKER, ", ", hh, "-hour window: ", nrow(d), " rows, ", n_distinct(d$id), " patients")
  if (n_distinct(d$id) < 50) return(NULL)
  imap_dfr(c(adjusted = TRUE, unadjusted = FALSE), function(adj, adj_lab) {
    f <- fit_lme(rhs_of(EXPO, if (adj) DEMO), d)
    cf <- contrast(f, EXPO, hh)
    tibble(marker = MARKER, exposure = EXPO, adjustment = adj_lab, model_horizon_h = hh,
           estimate = cf$est, lo = cf$est - 1.96 * sqrt(cf$V[1, 1]), hi = cf$est + 1.96 * sqrt(cf$V[1, 1]),
           level = unname(fixef(f)[EXPO]), divergence_per_day = unname(fixef(f)[term_by(fixef(f), c(EXPO, "vent_day"))]),
           dose_within = unname(fixef(f)["l_vtpbw_within"]), n_patients = n_distinct(d$id), n_rows = nrow(d))
  })
}))
if (nrow(out) == 0) stop("no window had 50 patients with two or more rows")

# =============================================================================
# 2. channel decomposition per exposure and window, with the supports
# =============================================================================
pt_all <- d_all %>% distinct(hospitalization_id, .keep_all = TRUE)
chan_rows <- list(); nested_rows <- list(); dose_rows <- list(); sf_rows <- list()
for (EXPO in names(EXPOS)) {
  EXPO_BASE <- EXPOS[[EXPO]]
  d_ch <- d_all %>%
    left_join(bind_cols(pt_all %>% select(hospitalization_id), pfvc_channels(pt_all, EXPO_BASE),
                        pt_all %>% transmute(vtpbw_c = vtpbw_pt_mean - median(vtpbw_pt_mean),
                                             log_sf_0_c = log(sf_0) - median(log(sf_0), na.rm = TRUE))),
              by = "hospitalization_id")
  for (hh in QUICK_HOURS) {
    d <- window_data(d_ch, hh)
    if (n_distinct(d$id) < 50) next
    pt <- d %>% distinct(hospitalization_id, .keep_all = TRUE)
    sd_expo <- sd(pt$ch_sum + pt$ch_remainder)
    # --- the pieces, each with level + divergence, beside the one-beta model on their sum
    free <- fit_lme(rhs_of(CHANNELS), d); one <- fit_lme(rhs_of("ch_sum"), d)
    cf <- contrast(free, CHANNELS, hh); c1 <- contrast(one, "ch_sum", hh)
    chan_rows[[length(chan_rows) + 1]] <-
      bind_rows(tibble(channel = c("height", "age", "sex", "race"), estimate = cf$est, se = sqrt(diag(cf$V))),
                tibble(channel = "all (one beta)", estimate = c1$est, se = sqrt(diag(c1$V)))) %>%
      mutate(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se, per_sd = estimate * sd_expo,
             per_100ml_at_median = if (EXPO_BASE == "log_pfvc") estimate * 0.1 / median(pt$pfvc_gli) else NA_real_,
             marker = MARKER, exposure = expo_label[[EXPO_BASE]], exposure_col = EXPO,
             model_horizon_h = hh, p_equal = channels_equal_p(cf$est, cf$V),
             n_patients = n_distinct(d$id), remainder_sd = sd(pt$ch_remainder))
    # --- nested ladder (ML fits so the likelihoods compare); every pair exactly nested
    fits <- list(
      demo_only     = fit_lme(rhs_of(NULL, DEMO), d, "ML"),
      one_beta      = fit_lme(rhs_of("ch_sum"), d, "ML"),
      one_beta_demo = fit_lme(rhs_of("ch_sum", DEMO), d, "ML"),
      pieces_free   = fit_lme(rhs_of(CHANNELS), d, "ML"),
      pieces_demo   = fit_lme(rhs_of(c("ch_height", "ch_age"), DEMO), d, "ML"),   # ch_sex/ch_race are spanned by sex/race
      expo_demo     = fit_lme(rhs_of(EXPO, DEMO), d, "ML"))
    lr_row <- function(small, big, test, nested = TRUE) {
      ll <- function(nm) logLik(fits[[nm]])
      lr <- if (nested) as.numeric(2 * (ll(big) - ll(small))) else NA_real_
      df <- if (nested) attr(ll(big), "df") - attr(ll(small), "df") else NA_real_
      tibble(test = test, small = small, big = big, lr = lr, df = df,
             p = if (nested) pchisq(lr, df, lower.tail = FALSE) else NA_real_,
             d_aic = AIC(fits[[big]]) - AIC(fits[[small]]), aic_small = AIC(fits[[small]]), aic_big = AIC(fits[[big]]))
    }
    nested_rows[[length(nested_rows) + 1]] <- bind_rows(
      lr_row("demo_only", "one_beta_demo", "a_vs_c: size adds to demographics"),
      lr_row("one_beta", "one_beta_demo", "b_vs_c: demographics add to size"),
      lr_row("one_beta", "pieces_free", "b_vs_d: the pieces disagree"),
      lr_row("pieces_free", "pieces_demo", "d_vs_top: free demographic shape beyond the pieces"),
      lr_row("one_beta_demo", "pieces_demo", "c_vs_top: pieces beyond one beta + demographics"),
      lr_row("one_beta_demo", "pieces_free", "d_vs_c: not nested, AIC only", nested = FALSE),
      imap_dfr(fits, ~ tibble(test = "model", small = NA_character_, big = .y, lr = NA_real_, df = NA_real_, p = NA_real_,
                              d_aic = NA_real_, aic_small = NA_real_, aic_big = AIC(.x)))) %>%
      mutate(marker = MARKER, exposure = expo_label[[EXPO_BASE]], exposure_col = EXPO, model_horizon_h = hh,
             n_patients = n_distinct(d$id), n_rows = nrow(d), remainder_sd = sd(pt$ch_remainder))
    # --- dose x piece and severity x piece: the pieces' contrasts as a function of the
    #     between-patient dose (mean VT/PBW, centred) and of the baseline SF (centred)
    interactions <- function(modifier, label, dose_mean = "vtpbw_pt_mean") {
      ints <- paste0(CHANNELS, ":", modifier)
      free <- fit_lme(rhs_of(c(CHANNELS, ints), dose_mean = dose_mean), d)
      one  <- fit_lme(rhs_of(c(EXPO, paste0(EXPO, ":", modifier)), DEMO, dose_mean = dose_mean), d)
      cf <- contrast(free, ints, hh); c1 <- contrast(one, paste0(EXPO, ":", modifier), hh)
      hm <- cf$est[1] - cf$est[2]; hm_se <- sqrt(cf$V[1, 1] + cf$V[2, 2] - 2 * cf$V[1, 2])
      bind_rows(
        tibble(model = "pieces_free", term = ints, estimate = cf$est, se = sqrt(diag(cf$V)), p_int_equal = channels_equal_p(cf$est, cf$V)),
        tibble(model = "pieces_free", term = "height_minus_age", estimate = hm, se = hm_se, p_int_equal = NA_real_),
        tibble(model = "single_index_adjusted", term = paste0(EXPO, ":", modifier), estimate = c1$est, se = sqrt(diag(c1$V)), p_int_equal = NA_real_)) %>%
        mutate(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se, p = 2 * pnorm(-abs(estimate / se)),
               marker = MARKER, exposure = expo_label[[EXPO_BASE]], exposure_col = EXPO, model_horizon_h = hh,
               modifier = label, n_patients = n_distinct(d$id))
    }
    if (HAS_DOSE) dose_rows[[length(dose_rows) + 1]] <- interactions("vtpbw_c", "patient-mean VT/PBW, mL/kg, centred", dose_mean = "vtpbw_c")
    if (!SF_IS_OUTCOME) sf_rows[[length(sf_rows) + 1]] <- interactions("log_sf_0_c", "log index-day worst SF, centred")
  }
}
chan   <- bind_rows(chan_rows);   nested <- bind_rows(nested_rows)
dose_ch <- bind_rows(dose_rows);  sf_ch  <- bind_rows(sf_rows)
if (!nrow(dose_ch)) dose_ch <- tibble(marker = MARKER, note = "skipped: no ventilator dose outside the ventilated cohort")
if (!nrow(sf_ch)) sf_ch <- tibble(marker = MARKER, note = "skipped: SF is this marker's own baseline")

message("\nMarker difference per SD of each size exposure at each window's horizon, one model per window ",
        "(log units; a lower PFVC is the negative of this). No correction for death before the horizon.")
print(as.data.frame(out %>% mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
message("\nChannel decomposition (log PFVC): the size effect per log unit, identified through each GLI input ",
        "(equal = lung size is the operative quantity; p_equal)")
print(as.data.frame(chan %>% filter(exposure_col == "log_pfvc_sd") %>%
                      select(model_horizon_h, channel, estimate, lo, hi, per_sd, per_100ml_at_median, p_equal, n_patients) %>%
                      mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
message("\nNested ladder (log PFVC, ML; d_aic < 0 favours the bigger model)")
print(as.data.frame(nested %>% filter(exposure_col == "log_pfvc_sd", test != "model") %>%
                      select(model_horizon_h, test, lr, df, p, d_aic) %>% mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
if ("estimate" %in% names(dose_ch)) {
  message("\nDose x piece (log PFVC): interaction contrast at the horizon per log unit of the piece per mL/kg")
  print(as.data.frame(dose_ch %>% filter(exposure_col == "log_pfvc_sd") %>%
                        select(model_horizon_h, model, term, estimate, lo, hi, p, p_int_equal) %>%
                        mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
}
tag <- function(x) x %>% mutate(grid = JM_GRID, site = site_name)
write_csv(tag(out),     file.path(final_dir, paste0("quick_lme_", MARKER, "_", site_name, ".csv")))
write_csv(tag(chan),    file.path(final_dir, paste0("quick_channels_", MARKER, "_", site_name, ".csv")))
write_csv(tag(nested),  file.path(final_dir, paste0("quick_nested_", MARKER, "_", site_name, ".csv")))
write_csv(tag(dose_ch), file.path(final_dir, paste0("quick_dose_channels_", MARKER, "_", site_name, ".csv")))
write_csv(tag(sf_ch),   file.path(final_dir, paste0("quick_sf_channels_", MARKER, "_", site_name, ".csv")))
message("26_quick_lme complete -> ", final_dir)
