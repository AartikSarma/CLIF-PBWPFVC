# =============================================================================
# Supplement (cross-sectional): does measured compliance scale like predicted FVC,
# input by input?
# =============================================================================
# Respiratory-system compliance (Crs = VT / driving pressure) measures how much
# volume the lung and chest wall take per unit pressure, and in a healthy lung it
# scales with lung size. If PFVC tracks lung size, log Crs rises one-for-one with
# log PFVC: an exponent of 1. Across sites, adjusted for age, sex and race, Crs
# rises about 9 mL/cmH2O per litre of PFVC, close to proportional (pooled figure,
# 2026-09-23). That estimate is identified by height alone. This script asks the
# same question input by input, and against PBW.
#
# 1. Channels. GLI-2012 log PFVC is, to a small remainder, a sum of a height piece,
#    an age curve, a sex shift and a race shift (pfvc_channels() in
#    20_biotrauma_grid.R, each in log-PFVC units). log Crs is regressed on the four
#    pieces at once, so each coefficient is the Crs exponent through that input.
#    Stated before the data (2026-09-23):
#      height, sex, race   the lung-size inputs: exponents near 1, and equal
#      age                 no directional prediction. GLI FVC falls with age, and
#                          on one reading of strain the usable tidal range (FVC
#                          limited, as residual volume rises) falls with it; on
#                          another the aerated end-expiratory volume at ZEEP is
#                          preserved. Bedside Crs also mixes a looser lung with a
#                          stiffer chest wall as age rises. Reported as it comes.
#    Tests: height = sex = race (Wald, 2 df); each piece against 1; all four equal
#    (3 df, for reference).
#
# 2. Head-to-head (Claim 3). log Crs on log PFVC against log Crs on log PBW, with
#    identical covariates and no demographics (the exposures are the demographic
#    scalings being compared): AIC and each exponent. Height moves PBW and PFVC
#    nearly in proportion, so the two models differ where the formulas part: short
#    women, where Devine's height elasticity (2.8 to 3.4) departs most from GLI's
#    (about 2.3). The head-to-head is repeated in women shorter than the median woman.
#    log PBW and log PFVC are never entered together. A third exposure, PFVC at age
#    25 (GLI's height, sex and race scaling without its age decline; script 03's
#    pfvc_age25), and a second version of every model with ns(age, 4) added, followed
#    MIMIC's first run (2026-09-23), where Crs tracked PFVC's height piece but not its
#    age piece and PBW won the head-to-head.
#
# 2b. Specific elastance (Ers x predicted size) as a second outcome throughout. Stress
#    = specific elastance x strain, so the two are interchangeable only where specific
#    elastance is constant. Its exponent through a GLI piece is 1 minus that piece's
#    Crs exponent, tested against 0: a non-zero value says the formula's predicted
#    size difference and the pressure per unit of relative distension part company
#    through that input. Age is where a loss of elastic recoil would show.
#
# 3. Height elasticity by sex: d log Crs / d log height within each sex, beside
#    GLI's height elasticity (2.41 men, 2.26 women) and Devine's at the sex's median
#    height (a line with an intercept, so its elasticity falls as height rises).
#
# 4. The figure's model on the log scale: log Crs on log PFVC with ns(age, 4), sex and
#    race, the height-identified exponent.
#
# Covariates in every model: log SF, SOFA, PEEP (compliance depends on the volume
# PEEP holds), BMI (Crs includes the chest wall; BMI enters because the outcome is
# a pressure-derived measure). Sample: the plateau-measured index timepoint
# (dp > 0; pressures are never forward-filled). Sensitivity: driving pressure >= 5
# cmH2O, because a driving pressure near 1 gives compliances up to about 530 that
# pass QC today. A model that warns stops the script.
#
# Inputs : intermediate/analysis_cross_sectional.parquet (script 03)
# Outputs: final/supplement/
#   crs_channels_estimates_{site}.csv   every exponent: model, term, estimate, SE,
#                                       CI, p, the value predicted before the data;
#                                       the specific-elastance rows carry their own
#                                       null (0) and the spread of log(Ers x size)
#   crs_channels_tests_{site}.csv       Wald tests of the channels, AIC head-to-head
#   crs_channels_{site}.pdf             channel exponents, head-to-head, height
#                                       elasticity by sex
# Usage: Rscript code/supplement/xsec_crs_channels.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(splines)
  library(patchwork)
  library(here)
})

options(width = 220)
source("utils/config.R")
site_name <- config$site_name
final_dir <- final_dir_for("supplement")
source(here("code", "20_biotrauma_grid.R"))   # pfvc_channels(), CHANNELS, channels_equal_p()

DP_FLOOR_SENSITIVITY <- 5   # cmH2O
GLI_HEIGHT_ELASTICITY <- c(Male = 2.41, Female = 2.26)
OKABE_ITO <- c(height = "#0072B2", age = "#E69F00", sex = "#009E73", race = "#CC79A7",
               PFVC = "#0072B2", PBW = "#D55E00")

# a model that warns stops the script (no silent fallback)
fit_strict <- function(expr) withCallingHandlers(expr, warning = function(w)
  stop("model warning, stopping: ", conditionMessage(w), call. = FALSE))

# =============================================================================
# Data: the plateau-measured index timepoint
# =============================================================================
mechanics_all <- read_parquet(file.path(config$output_dir, "analysis_cross_sectional.parquet")) %>%
  filter(!is.na(crs), crs > 0, !is.na(dp), dp > 0) %>%
  filter(!is.na(pfvc), pfvc > 0, !is.na(pfvc_age25), pfvc_age25 > 0, !is.na(pbw), pbw > 0, !is.na(height_cm), !is.na(age_at_admission),
         !is.na(sex_category), !is.na(race_category), !is.na(sf_ratio), sf_ratio > 0,
         !is.na(sofa_total), !is.na(peep_set), !is.na(bmi)) %>%
  mutate(sex_category  = factor(sex_category, levels = c("Male", "Female")),
         race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")),
         age10 = age_at_admission / 10,
         log_crs = log(crs), log_pfvc = log(pfvc), log_pfvc25 = log(pfvc_age25), log_pbw = log(pbw),
         log_height = log(height_cm), log_sf = log(sf_ratio),
         # specific elastance: Ers x the predicted size, the pressure it takes to
         # double the lung's volume. log Espec = log Ers + log PFVC = log PFVC - log Crs
         # (plus the constant that carries Crs from mL to L), so its exponent through a
         # GLI piece is 1 minus that piece's Crs exponent (a piece moves log PFVC by one
         # unit and log Crs by its exponent), and the two models are one model read two
         # ways, up to the small GLI remainder. Both are fitted, because the test
         # differs: a Crs
         # exponent of 1 is the same statement as an Espec exponent of 0, and the
         # second is the one a reader weighs against "specific elastance is constant".
         log_espec = log(ers * pfvc), log_espec25 = log(ers * pfvc_age25))
if (nrow(mechanics_all) < 100) stop("fewer than 100 patients with a measured plateau and every covariate")
SAMPLES <- list(`all plateau-measured` = mechanics_all,
                `driving pressure >= 5` = mechanics_all %>% filter(dp >= DP_FLOOR_SENSITIVITY))
message(sprintf("=== xsec_crs_channels, %s: %d plateau-measured patients (%d with driving pressure >= %d) ===",
                site_name, nrow(mechanics_all), nrow(SAMPLES[[2]]), DP_FLOOR_SENSITIVITY))

COVARIATES <- "log_sf + sofa_total + peep_set + bmi"
DEMOGRAPHICS <- "ns(age10, 4) + sex_category + race_category"
coef_rows <- function(fit, terms, model, sample_label, predicted = NA_real_) {
  b <- coef(fit)[terms]; se <- sqrt(diag(vcov(fit))[terms])
  tibble(sample = sample_label, model = model, term = terms, estimate = unname(b), se = unname(se),
         lo = estimate - 1.96 * se, hi = estimate + 1.96 * se, p = 2 * pnorm(-abs(estimate / se)),
         predicted = predicted, n_patients = nobs(fit))
}
wald <- function(fit, contrast_matrix, rhs, test_label, sample_label) {
  b <- coef(fit)[colnames(contrast_matrix)]; V <- vcov(fit)[colnames(contrast_matrix), colnames(contrast_matrix)]
  d <- contrast_matrix %*% b - rhs
  stat <- as.numeric(t(d) %*% solve(contrast_matrix %*% V %*% t(contrast_matrix)) %*% d)
  tibble(sample = sample_label, test = test_label, statistic = stat, df = nrow(contrast_matrix),
         p = pchisq(stat, nrow(contrast_matrix), lower.tail = FALSE), n_patients = nobs(fit))
}

analyse_sample <- function(dat, sample_label) {
  estimates <- list(); tests <- list()

  # ---- 1. channels: log Crs, and log specific elastance, on the four GLI pieces of
  #        log PFVC. Compliance answers "does the formula's predicted size difference
  #        show up at the bedside?" (exponent 1 = proportional); specific elastance
  #        answers "does the pressure per unit of relative distension differ by group?"
  #        (exponent 0 = constant specific elastance, so stress tracks strain through
  #        that input). Stress = specific elastance x strain, so an input with a
  #        non-zero Espec exponent breaks the proportionality between them.
  pieces <- bind_cols(dat, pfvc_channels(dat, "log_pfvc"))
  size_channels <- c("ch_height", "ch_sex", "ch_race")
  OUTCOMES <- tribble(
    ~outcome,            ~model,            ~null_value, ~scale,
    "log_crs",           "channels",        1,           "exponent of log Crs (1 = proportional to predicted FVC)",
    "log_espec",         "espec channels",  0,           "exponent of log specific elastance (0 = constant specific elastance)")
  channel_fits <- list()
  for (row_i in seq_len(nrow(OUTCOMES))) {
    outcome <- OUTCOMES$outcome[row_i]; model_label <- OUTCOMES$model[row_i]; null_value <- OUTCOMES$null_value[row_i]
    fit <- fit_strict(lm(as.formula(paste(outcome, "~", paste(CHANNELS, collapse = " + "), "+", COVARIATES)), data = pieces))
    channel_fits[[model_label]] <- fit
    estimates[[model_label]] <- coef_rows(fit, CHANNELS, model_label, sample_label,
                                          predicted = if (null_value == 1) c(1, NA, 1, 1) else c(0, NA, 0, 0)) %>%
      mutate(term = sub("^ch_", "", term), channel_sd = sapply(pieces[CHANNELS], sd), scale = OUTCOMES$scale[row_i])
    tests[[model_label]] <- bind_rows(
      wald(fit, rbind(c(1, -1, 0), c(1, 0, -1)) %>% `colnames<-`(size_channels), c(0, 0),
           paste0(model_label, ": height = sex = race (the lung-size inputs agree)"), sample_label),
      map_dfr(CHANNELS, ~ wald(fit, matrix(1, 1, 1, dimnames = list(NULL, .x)), null_value,
                               sprintf("%s: %s exponent = %g", model_label, sub("^ch_", "", .x), null_value), sample_label)),
      tibble(sample = sample_label, test = paste0(model_label, ": all four equal (height, age, sex, race)"),
             statistic = NA_real_, df = 3, p = channels_equal_p(coef(fit)[CHANNELS], vcov(fit)[CHANNELS, CHANNELS]),
             n_patients = nobs(fit)))
  }
  channel_fit <- channel_fits[["channels"]]

  # ---- 2. head-to-head: log PFVC against log PBW, identical covariates, no demographics
  # three exposures, each on its own, against PBW: PFVC; PFVC at age 25 (GLI's height,
  # sex and race scaling with its age decline removed, script 03's pfvc_age25); PBW.
  # Fitted without demographics, and again with ns(age, 4) in every model, which asks
  # whether the formulas differ once age is held fixed (MIMIC, 2026-09-23: Crs follows
  # PFVC's height piece but not its age piece, so PFVC lost to PBW overall).
  H2H_EXPOSURES <- c(PFVC = "log_pfvc", `PFVC at age 25` = "log_pfvc25", PBW = "log_pbw")
  # the same head-to-head on specific elastance: which scaling makes Ers x size most
  # nearly constant across patients (the Chiumello reading of a correct normaliser)
  ESPEC_EXPOSURES <- c(PFVC = "log_espec", `PFVC at age 25` = "log_espec25")
  head_to_head <- function(sub_dat, subgroup, age_adjusted) {
    label <- paste0("head-to-head, ", subgroup, if (age_adjusted) ", age spline in every model" else "")
    fits <- map(H2H_EXPOSURES, function(exposure)
      fit_strict(lm(as.formula(paste("log_crs ~", exposure, "+", if (age_adjusted) "ns(age10, 4) +", COVARIATES)), data = sub_dat)))
    list(estimates = imap_dfr(fits, ~ coef_rows(.x, H2H_EXPOSURES[[.y]], label, sample_label, 1) %>% mutate(exposure = .y)),
         tests = imap_dfr(fits[names(fits) != "PBW"], ~ tibble(
           sample = sample_label, test = paste0(label, ": AIC, ", .y, " minus PBW (< 0 favours ", .y, ")"),
           statistic = AIC(.x) - AIC(fits$PBW), df = NA_real_, p = NA_real_, n_patients = nobs(.x))))
  }
  female_median_height <- median(dat$height_cm[dat$sex_category == "Female"])
  short_women <- dat %>% filter(sex_category == "Female", height_cm < female_median_height)
  short_women_label <- sprintf("women shorter than %.0f cm", female_median_height)
  h2h <- list(head_to_head(dat, "everyone", FALSE), head_to_head(dat, "everyone", TRUE),
              if (nrow(short_women) >= 50) head_to_head(short_women, short_women_label, FALSE),
              if (nrow(short_women) >= 50) head_to_head(short_women, short_women_label, TRUE))
  estimates$head_to_head <- map_dfr(compact(h2h), "estimates")
  tests$head_to_head <- map_dfr(compact(h2h), "tests")

  # ---- 2b. is specific elastance more nearly constant under one scaling than the other?
  #      The spread of log(Ers x size) across patients, and its residual spread after the
  #      covariates: a smaller spread means the normaliser leaves less unexplained
  #      variation in the pressure per unit of relative distension.
  estimates$espec_spread <- imap_dfr(ESPEC_EXPOSURES, function(column, label) {
    fit <- fit_strict(lm(as.formula(paste(column, "~", COVARIATES)), data = dat))
    tibble(sample = sample_label, model = "specific elastance spread", term = label,
           estimate = sd(dat[[column]]), se = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_,
           predicted = NA_real_, n_patients = nrow(dat), residual_sd = sd(resid(fit)),
           scale = "SD of log specific elastance (smaller = the scaling leaves less unexplained)")
  })

  # ---- 3. height elasticity by sex, beside GLI's and Devine's
  estimates$height_by_sex <- map_dfr(c("Male", "Female"), function(sx) {
    sub_dat <- dat %>% filter(sex_category == sx)
    if (nrow(sub_dat) < 50) return(NULL)
    fit <- fit_strict(lm(as.formula(paste("log_crs ~ log_height + ns(age10, 4) + race_category +", COVARIATES)), data = sub_dat))
    median_height <- median(sub_dat$height_cm)
    # Devine: PBW = a + 2.3 (h/2.54 - 60), so d log PBW / d log h = (2.3 h / 2.54) / PBW
    devine_at_median <- (2.3 * median_height / 2.54) / (if_else(sx == "Male", 50, 45.5) + 2.3 * (median_height / 2.54 - 60))
    coef_rows(fit, "log_height", paste("height elasticity,", sx), sample_label) %>%
      mutate(gli_elasticity = GLI_HEIGHT_ELASTICITY[[sx]], devine_elasticity_at_median = devine_at_median,
             median_height_cm = median_height)
  })

  # ---- 4. the figure's model on the log scale: the height-identified exponent
  demo_fit <- fit_strict(lm(as.formula(paste("log_crs ~ log_pfvc +", DEMOGRAPHICS, "+", COVARIATES)), data = dat))
  estimates$demographic_adjusted <- coef_rows(demo_fit, "log_pfvc", "log PFVC + ns(age, 4), sex, race", sample_label, 1)

  list(estimates = bind_rows(estimates), tests = bind_rows(tests))
}
results <- imap(SAMPLES, analyse_sample)
estimates <- map_dfr(results, "estimates") %>% mutate(site = site_name)
tests <- map_dfr(results, "tests") %>% mutate(site = site_name)

message("\nExponents of Crs (per log unit; 1 = proportional scaling):")
print(as.data.frame(estimates %>% select(sample, model, term, estimate, lo, hi, p, predicted, n_patients) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
message("\nTests:")
print(as.data.frame(tests %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
write_csv(mask_small_counts(estimates), file.path(final_dir, paste0("crs_channels_estimates_", site_name, ".csv")))
write_csv(mask_small_counts(tests), file.path(final_dir, paste0("crs_channels_tests_", site_name, ".csv")))

# =============================================================================
# Figure
# =============================================================================
primary <- names(SAMPLES)[1]
channel_panel <- estimates %>% filter(model == "channels") %>%
  mutate(term = factor(term, levels = rev(c("height", "sex", "race", "age"))),
         sample = factor(sample, levels = names(SAMPLES))) %>%
  ggplot(aes(estimate, term, colour = term, shape = sample)) +
  geom_vline(xintercept = c(0, 1), linetype = c(2, 1), colour = c("grey60", "grey30")) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = OKABE_ITO, guide = "none") + scale_shape_manual(values = c(16, 1), name = NULL) +
  labs(title = "A. Crs exponent through each input of log PFVC",
       subtitle = "1 = Crs proportional to predicted FVC; predicted near 1 for height, sex, race; age not predicted",
       x = "d log Crs / d log PFVC (piece)", y = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
h2h_panel <- estimates %>% filter(startsWith(model, "head-to-head"), sample == primary) %>%
  mutate(model = sub("head-to-head, ", "", sub(", age spline in every model", ",\nage spline", model))) %>%
  ggplot(aes(estimate, model, colour = exposure)) +
  geom_vline(xintercept = 1, colour = "grey30") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c(PFVC = "#0072B2", `PFVC at age 25` = "#009E73", PBW = "#D55E00"), name = NULL) +
  labs(title = "B. Head-to-head: Crs exponent on log PFVC, PFVC at age 25 and PBW",
       subtitle = "identical covariates; with and without an age spline", x = "exponent", y = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
elasticity_panel <- estimates %>% filter(startsWith(model, "height elasticity"), sample == primary) %>%
  mutate(sex = sub("height elasticity, ", "", model)) %>%
  ggplot(aes(estimate, sex)) +
  geom_pointrange(aes(xmin = lo, xmax = hi), colour = "grey20") +
  geom_point(aes(x = gli_elasticity, colour = "GLI FVC"), shape = 17, size = 3) +
  geom_point(aes(x = devine_elasticity_at_median, colour = "Devine PBW at the median height"), shape = 15, size = 3) +
  scale_colour_manual(values = c(`GLI FVC` = "#0072B2", `Devine PBW at the median height` = "#D55E00"), name = NULL) +
  labs(title = "C. Height elasticity of Crs within sex",
       subtitle = "d log Crs / d log height (age, race, covariates adjusted)", x = "elasticity", y = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(final_dir, paste0("crs_channels_", site_name, ".pdf")),
       channel_panel / (h2h_panel | elasticity_panel), width = 12, height = 8)
message("xsec_crs_channels complete -> ", final_dir)
