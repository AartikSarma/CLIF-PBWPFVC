# =============================================================================
# Supplement (cross-sectional): which dose or mechanics measure, alone, predicts
# death best?
# =============================================================================
# Twelve measures at the index timepoint, each in its own model with no covariates:
#   dose        VT/PBW, VT/PFVC, VT/PFVC at age 25
#   elastance   Ers x PBW, Ers x PFVC, Ers x PFVC at age 25 (specific elastance by
#               each size scaling)
#   power       MP, MP/Crs, MP/PBW, MP/PFVC, MP/PFVC at age 25
#   pressure    driving pressure
# PFVC at age 25 is GLI's prediction at a common reference age (script 03's
# pfvc_age25): the height, sex and race scaling without GLI's age decline.
#
# This is a question about prediction, not cause. With no covariates, each measure
# carries its own demographic content: PFVC-normalised measures inherit GLI's age
# term, which predicts death through everything age does. The PFVC-at-age-25
# versions separate that age content from the structural scaling, so reading PFVC
# against PFVC at age 25 shows how much of a normalisation's predictive edge is age.
#
# Methods
#   outcomes   in-hospital death, and death by day 60 (binary), logistic
#   form       each measure as a 3-df natural spline of its log (no measure is
#              penalised for a shape the others do not have); a linear-in-log form
#              is reported beside it
#   sample     the patients with every measure (plateau-measured, mechanical power
#              defined), so all twelve are compared on the same people; each
#              measure is also fitted on every patient who has it, as a secondary
#              read (not comparable across measures)
#   adjustment each measure alone, and given VT/PBW (beside a 3-df spline of log
#              VT/PBW, the dose the clinician set: does the measure add anything once
#              the dose is known?). Given VT/PBW, VT/PFVC carries the PBW/PFVC ratio.
#              Third, given VT/PBW, sex and race: the structural ratio (VT/PFVC at
#              age 25) is then mostly height's sex-specific hump.
#   sign       the linear-in-log form reports each measure's log-odds per log unit
#   metrics    10-fold cross-validated AUC, Brier score and log loss (out-of-fold
#              predictions, the same folds for every measure), in-sample AIC, and
#              the AUC difference from the adjustment's base model (VT/PBW; VT/PBW
#              with sex and race in the third) with a paired bootstrap interval
#              (patients resampled, the out-of-fold predictions kept); within each
#              family, the PBW, PFVC and PFVC-at-age-25 scalings against each other
# A model that warns stops the script.
#
# Inputs : intermediate/analysis_cross_sectional.parquet (script 03)
# Outputs: final/supplement/
#   mortality_prediction_{site}.csv   one row per measure x outcome x form x
#                                     adjustment x sample
#   mortality_prediction_pairwise_{site}.csv  the size scalings within each family
#                                     against each other (common sample)
#   mortality_prediction_{site}.pdf   cross-validated AUC with the difference from
#                                     VT/PBW alone, alone and given VT/PBW
# Usage: Rscript code/supplement/xsec_mortality_prediction.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(splines)
})

options(width = 220)
source("utils/config.R")
site_name <- config$site_name
final_dir <- final_dir_for("supplement")

N_FOLDS <- 10
N_BOOT <- 1000
MIN_DEATHS <- 10L
set.seed(20260923)
MEASURES <- tribble(
  ~column,          ~label,                     ~family,
  "vtpbw",          "VT/PBW",                   "dose",
  "vtpfvc",         "VT/PFVC",                  "dose",
  "vtpfvc_age25",   "VT/PFVC at age 25",        "dose",
  "ers_pbw",        "Ers x PBW",                "elastance",
  "ers_pfvc",       "Ers x PFVC",               "elastance",
  "ers_pfvc_age25", "Ers x PFVC at age 25",     "elastance",
  "mechanical_power", "MP",                     "power",
  "mp_crs",         "MP/Crs",                   "power",
  "mp_pbw",         "MP/PBW",                   "power",
  "mp_pfvc",        "MP/PFVC",                  "power",
  "mp_pfvc_age25",  "MP/PFVC at age 25",        "power",
  "dp",             "Driving pressure",         "pressure")
OUTCOMES <- c(deceased = "in-hospital death", mortality_event_60 = "death by day 60")
FORMS <- c(spline = "ns(log_x, 3)", linear = "log_x")
FAMILY_COLOURS <- c(dose = "#0072B2", elastance = "#E69F00", power = "#009E73", pressure = "#CC79A7")

# a model that warns stops the script (no silent fallback)
fit_strict <- function(expr) withCallingHandlers(expr, warning = function(w)
  stop("model warning, stopping: ", conditionMessage(w), call. = FALSE))

# =============================================================================
# Data
# =============================================================================
cross_sectional <- read_parquet(file.path(config$output_dir, "analysis_cross_sectional.parquet")) %>%
  select(hospitalization_id, all_of(MEASURES$column), deceased, mortality_event_60, sex_category, race_category) %>%
  mutate(sex_category = factor(sex_category, levels = c("Male", "Female")),
         race_category = factor(race_category, levels = c("WHITE", "BLACK", "OTHER")))

# SYNTHETIC SITE ONLY: synthetic CLIF mortality is unreliable, so death is simulated
# independently of every measure (35%), as the other supplement scripts do. The run
# exercises the machinery; every AUC should sit near 0.5. Never runs at a real site.
if (grepl("^synthetic_clif", site_name)) {
  message("*** SYNTHETIC SITE: simulated mortality (plumbing only; synthetic CLIF mortality is unreliable). ***")
  simulated_death <- rbinom(nrow(cross_sectional), 1L, 0.35)
  cross_sectional <- cross_sectional %>% mutate(deceased = simulated_death, mortality_event_60 = simulated_death)
}

# a measure is usable where it is positive and finite (every one is logged)
usable <- cross_sectional %>%
  mutate(across(all_of(MEASURES$column), ~ if_else(is.finite(.x) & .x > 0, .x, NA_real_)))
common_ids <- usable %>% filter(if_all(all_of(MEASURES$column), ~ !is.na(.x)), !is.na(sex_category), !is.na(race_category)) %>%
  pull(hospitalization_id)
message(sprintf("=== xsec_mortality_prediction, %s: %d patients; %d with every measure (the common sample) ===",
                site_name, nrow(usable), length(common_ids)))
print(as.data.frame(tibble(measure = MEASURES$label,
                           n_patients = map_int(MEASURES$column, ~ sum(!is.na(usable[[.x]])))), row.names = FALSE))

# =============================================================================
# Metrics
# =============================================================================
auc_of <- function(y, p) {   # rank (Mann-Whitney) AUC
  r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
log_loss_of <- function(y, p) { p <- pmin(pmax(p, 1e-12), 1 - 1e-12); -mean(y * log(p) + (1 - y) * log(1 - p)) }

# out-of-fold predictions for one measure, outcome and form, on the rows of `dat`
# (folds are assigned once per sample, so every measure shares them)
cross_validated <- function(dat, folds, rhs) {
  predictions <- numeric(nrow(dat))
  for (k in seq_len(N_FOLDS)) {
    training <- dat[folds != k, ]; held_out <- dat[folds == k, ]
    # the spline basis is built on the training fold and carried to the held-out one
    fit <- fit_strict(glm(as.formula(paste("y ~", rhs)), family = binomial, data = training))
    predictions[folds == k] <- predict(fit, newdata = held_out, type = "response")
  }
  predictions
}

# Two adjustments. "alone": the measure by itself. "given VT/PBW": the measure beside a
# 3-df spline of log VT/PBW, the dose the clinician set, which is the clinical question
# (does this measure add anything once the dose is known?). VT/PBW itself is the base
# model of the second: its row there is the VT/PBW-alone model, and every other row's
# difference from it is the measure's gain over the dose. Given VT/PBW, log VT/PFVC =
# log VT/PBW + log PBW/PFVC, so VT/PFVC carries the PBW/PFVC ratio and its age-25
# version the ratio's structural part.
# A third adjustment adds sex and race to the second (2026-09-23): given VT/PBW, the
# structural ratio (VT/PFVC at age 25) is sex, race and height's small sex-specific
# hump, and sex and race predict death through routes other than the ventilator. Its
# base model is VT/PBW with sex and race, and every row's difference is from that base.
ADJUSTMENTS <- c(alone = "", given_vtpbw = " + ns(log_vtpbw, 3)",
                 given_vtpbw_sex_race = " + ns(log_vtpbw, 3) + sex_category + race_category")
# what the VT/PBW row carries in each adjustment: it is the base model
BASE_EXTRA <- c(alone = "", given_vtpbw = "", given_vtpbw_sex_race = " + sex_category + race_category")
ADJUSTMENT_LABELS <- c(alone = "alone", given_vtpbw = "given VT/PBW", given_vtpbw_sex_race = "given VT/PBW, sex and race")
# the three size scalings within each family, compared pairwise on the common sample
FAMILY_TRIPLETS <- list(dose = c("vtpbw", "vtpfvc", "vtpfvc_age25"),
                        elastance = c("ers_pbw", "ers_pfvc", "ers_pfvc_age25"),
                        power = c("mp_pbw", "mp_pfvc", "mp_pfvc_age25"))

evaluate_sample <- function(sample_label, ids_for) {
  pairwise <- list()
  rows <- map_dfr(names(OUTCOMES), function(outcome) {
    map_dfr(names(FORMS), function(form) {
      map_dfr(names(ADJUSTMENTS), function(adjustment) {
        per_measure <- map(MEASURES$column, function(measure) {
          dat <- usable %>% filter(hospitalization_id %in% ids_for(measure), !is.na(.data[[outcome]]), !is.na(vtpbw)) %>%
            transmute(hospitalization_id, y = .data[[outcome]], log_x = log(.data[[measure]]), log_vtpbw = log(vtpbw),
                      sex_category, race_category) %>%
            filter(!is.na(sex_category), !is.na(race_category))
          if (sum(dat$y == 1) < MIN_DEATHS || sum(dat$y == 0) < MIN_DEATHS) return(NULL)
          # VT/PBW given VT/PBW is VT/PBW alone: the base model
          # VT/PBW in each adjustment is that adjustment's base model
          rhs <- paste0(FORMS[[form]], if (measure != "vtpbw") ADJUSTMENTS[[adjustment]] else BASE_EXTRA[[adjustment]])
          folds <- sample(rep_len(seq_len(N_FOLDS), nrow(dat)))
          predicted <- cross_validated(dat, folds, rhs)
          full_fit <- fit_strict(glm(as.formula(paste("y ~", rhs)), family = binomial, data = dat))
          # the sign: log-odds per log unit of the measure, in the linear-in-log form
          slope <- if (form == "linear") coef(full_fit)[["log_x"]] else NA_real_
          slope_se <- if (form == "linear") sqrt(vcov(full_fit)["log_x", "log_x"]) else NA_real_
          list(ids = dat$hospitalization_id, y = dat$y, predicted = predicted, aic = AIC(full_fit),
               slope = slope, slope_se = slope_se)
        }) %>% set_names(MEASURES$column)
        present <- compact(per_measure)
        rows_here <- imap_dfr(present, function(m, measure) tibble(
          measure = measure, n_patients = length(m$y), n_deaths = sum(m$y == 1),
          auc = auc_of(m$y, m$predicted), brier = mean((m$predicted - m$y)^2),
          log_loss = log_loss_of(m$y, m$predicted), aic = m$aic,
          log_or_per_log_unit = m$slope, log_or_se = m$slope_se))
        # paired bootstrap: on the common sample every measure holds the same patients in
        # the same order, so one resample of rows serves all twelve
        if (sample_label == "common" && !is.null(present$vtpbw)) {
          n <- length(present$vtpbw$y); y <- present$vtpbw$y
          stopifnot(all(map_lgl(present, ~ identical(.x$ids, present$vtpbw$ids))))
          boot <- replicate(N_BOOT, {
            i <- sample.int(n, n, replace = TRUE)
            map_dbl(present, ~ auc_of(y[i], .x$predicted[i]))
          })
          versus_vtpbw <- boot - rep(boot["vtpbw", ], each = nrow(boot))
          rows_here <- rows_here %>% mutate(
            auc_lo = apply(boot, 1, quantile, 0.025)[measure], auc_hi = apply(boot, 1, quantile, 0.975)[measure],
            delta_auc_vs_vtpbw = auc - auc[measure == "vtpbw"],
            delta_auc_lo = apply(versus_vtpbw, 1, quantile, 0.025)[measure],
            delta_auc_hi = apply(versus_vtpbw, 1, quantile, 0.975)[measure],
            delta_aic_vs_vtpbw = aic - aic[measure == "vtpbw"])
          # pairwise within each family: PFVC against PBW, PFVC at age 25 against PBW,
          # PFVC against PFVC at age 25 (the age term's share)
          pairwise[[length(pairwise) + 1]] <<- imap_dfr(FAMILY_TRIPLETS, function(triplet, family) {
            if (!all(triplet %in% rownames(boot))) return(NULL)
            pairs <- list(c(triplet[2], triplet[1]), c(triplet[3], triplet[1]), c(triplet[2], triplet[3]))
            map_dfr(pairs, function(pair) {
              difference <- boot[pair[1], ] - boot[pair[2], ]
              observed <- rows_here$auc[rows_here$measure == pair[1]] - rows_here$auc[rows_here$measure == pair[2]]
              tibble(family = family, measure = pair[1], minus = pair[2], delta_auc = observed,
                     lo = quantile(difference, 0.025), hi = quantile(difference, 0.975))
            })
          }) %>% mutate(outcome = OUTCOMES[[outcome]], form = form, adjustment = ADJUSTMENT_LABELS[[adjustment]])
        }
        rows_here %>% mutate(outcome = OUTCOMES[[outcome]], form = form,
                             adjustment = ADJUSTMENT_LABELS[[adjustment]], sample = sample_label)
      })
    })
  })
  list(rows = rows, pairwise = bind_rows(pairwise))
}
common_results <- evaluate_sample("common", function(measure) common_ids)
own_results <- evaluate_sample("each measure's own", function(measure) usable$hospitalization_id[!is.na(usable[[measure]])])
results <- bind_rows(common_results$rows, own_results$rows) %>%
  left_join(MEASURES, by = c("measure" = "column")) %>%
  mutate(sample = if_else(sample == "common", "common (every measure present)", "each measure's own (not comparable)"),
         site = site_name) %>%
  relocate(sample, outcome, form, adjustment, family, label)
measure_label <- setNames(MEASURES$label, MEASURES$column)
pairwise <- common_results$pairwise %>%
  mutate(contrast = paste(measure_label[measure], "minus", measure_label[minus]), site = site_name) %>%
  relocate(outcome, form, adjustment, family, contrast)

for (adjustment_label in ADJUSTMENT_LABELS) {
  message("\nCross-validated AUC on the common sample, spline form, ", adjustment_label,
          " (difference from the base model, with a paired bootstrap interval):")
  print(as.data.frame(results %>% filter(startsWith(sample, "common"), form == "spline", adjustment == adjustment_label) %>%
                        arrange(outcome, desc(auc)) %>%
                        select(outcome, label, auc, auc_lo, auc_hi, delta_auc_vs_vtpbw, delta_auc_lo, delta_auc_hi,
                               delta_aic_vs_vtpbw, brier, n_patients, n_deaths) %>%
                        mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
}
message("\nWithin each family, the size scalings against each other (spline form, AUC difference, paired bootstrap):")
print(as.data.frame(pairwise %>% filter(form == "spline") %>%
                      select(outcome, adjustment, contrast, delta_auc, lo, hi) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
message("\nThe sign of the dose measures (linear in log, log-odds per log unit; positive = higher value, more death):")
print(as.data.frame(results %>% filter(startsWith(sample, "common"), form == "linear", family == "dose") %>%
                      transmute(outcome, adjustment, label, log_or_per_log_unit, lo = log_or_per_log_unit - 1.96 * log_or_se,
                                hi = log_or_per_log_unit + 1.96 * log_or_se) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
write_csv(mask_small_counts(results), file.path(final_dir, paste0("mortality_prediction_", site_name, ".csv")))
write_csv(pairwise, file.path(final_dir, paste0("mortality_prediction_pairwise_", site_name, ".csv")))

# =============================================================================
# Figure: cross-validated AUC and its difference from VT/PBW, common sample
# =============================================================================
figure_rows <- results %>% filter(startsWith(sample, "common"), form == "spline", outcome == OUTCOMES[["deceased"]]) %>%
  mutate(label = factor(label, levels = rev(MEASURES$label)), adjustment = factor(adjustment, levels = ADJUSTMENT_LABELS))
auc_panel <- ggplot(figure_rows, aes(auc, label, colour = family)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(xmin = auc_lo, xmax = auc_hi)) +
  facet_wrap(~ adjustment) +
  scale_colour_manual(values = FAMILY_COLOURS, name = NULL) +
  labs(title = "Cross-validated AUC for in-hospital death: each measure alone, given VT/PBW, and given VT/PBW, sex and race",
       subtitle = sprintf("3-df spline of the log measure (and of log VT/PBW); %d-fold cross-validation; the same patients for every measure", N_FOLDS),
       x = "AUC (95% bootstrap interval)", y = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
delta_panel <- ggplot(figure_rows, aes(delta_auc_vs_vtpbw, label, colour = family)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(xmin = delta_auc_lo, xmax = delta_auc_hi)) +
  facet_wrap(~ adjustment) +
  scale_colour_manual(values = FAMILY_COLOURS, guide = "none") +
  labs(title = "Difference from the base model (VT/PBW; with sex and race in the third panel)",
       x = "AUC minus the base model's AUC (paired bootstrap)", y = NULL) +
  theme_minimal(base_size = 10)
ggsave(file.path(final_dir, paste0("mortality_prediction_", site_name, ".pdf")),
       patchwork::wrap_plots(auc_panel, delta_panel, ncol = 1), width = 14, height = 10)
message("xsec_mortality_prediction complete -> ", final_dir)
