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
#   metrics    10-fold cross-validated AUC, Brier score and log loss (out-of-fold
#              predictions, the same folds for every measure), in-sample AIC, and
#              the AUC difference from VT/PBW with a paired bootstrap interval
#              (patients resampled, the out-of-fold predictions kept)
# A model that warns stops the script.
#
# Inputs : intermediate/analysis_cross_sectional.parquet (script 03)
# Outputs: final/supplement/
#   mortality_prediction_{site}.csv   one row per measure x outcome x form x sample
#   mortality_prediction_{site}.pdf   cross-validated AUC with the difference from
#                                     VT/PBW, on the common sample
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
  select(hospitalization_id, all_of(MEASURES$column), deceased, mortality_event_60)

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
common_ids <- usable %>% filter(if_all(all_of(MEASURES$column), ~ !is.na(.x))) %>% pull(hospitalization_id)
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

evaluate_sample <- function(sample_label, ids_for) {
  map_dfr(names(OUTCOMES), function(outcome) {
    map_dfr(names(FORMS), function(form) {
      per_measure <- map(MEASURES$column, function(measure) {
        dat <- usable %>% filter(hospitalization_id %in% ids_for(measure), !is.na(.data[[outcome]])) %>%
          transmute(hospitalization_id, y = .data[[outcome]], log_x = log(.data[[measure]]))
        if (sum(dat$y == 1) < MIN_DEATHS || sum(dat$y == 0) < MIN_DEATHS) return(NULL)
        folds <- sample(rep_len(seq_len(N_FOLDS), nrow(dat)))
        predicted <- cross_validated(dat, folds, FORMS[[form]])
        full_fit <- fit_strict(glm(as.formula(paste("y ~", FORMS[[form]])), family = binomial, data = dat))
        list(ids = dat$hospitalization_id, y = dat$y, predicted = predicted, aic = AIC(full_fit))
      }) %>% set_names(MEASURES$column)
      rows <- imap_dfr(compact(per_measure), function(m, measure) tibble(
        measure = measure, n_patients = length(m$y), n_deaths = sum(m$y == 1),
        auc = auc_of(m$y, m$predicted), brier = mean((m$predicted - m$y)^2),
        log_loss = log_loss_of(m$y, m$predicted), aic = m$aic))
      # paired bootstrap: on the common sample every measure holds the same patients in
      # the same order, so one resample of rows serves all twelve
      if (sample_label == "common" && !is.null(per_measure$vtpbw)) {
        n <- length(per_measure$vtpbw$y); y <- per_measure$vtpbw$y
        stopifnot(all(map_lgl(compact(per_measure), ~ identical(.x$ids, per_measure$vtpbw$ids))))
        boot <- replicate(N_BOOT, {
          i <- sample.int(n, n, replace = TRUE)
          map_dbl(compact(per_measure), ~ auc_of(y[i], .x$predicted[i]))
        })
        rows <- rows %>% mutate(
          auc_lo = apply(boot, 1, quantile, 0.025)[measure], auc_hi = apply(boot, 1, quantile, 0.975)[measure],
          delta_auc_vs_vtpbw = auc - auc[measure == "vtpbw"],
          delta_auc_lo = apply(boot - rep(boot["vtpbw", ], each = nrow(boot)), 1, quantile, 0.025)[measure],
          delta_auc_hi = apply(boot - rep(boot["vtpbw", ], each = nrow(boot)), 1, quantile, 0.975)[measure],
          delta_aic_vs_vtpbw = aic - aic[measure == "vtpbw"])
      }
      rows %>% mutate(outcome = OUTCOMES[[outcome]], form = form, sample = sample_label)
    })
  })
}
results <- bind_rows(
  evaluate_sample("common", function(measure) common_ids),
  evaluate_sample("each measure's own", function(measure) usable$hospitalization_id[!is.na(usable[[measure]])])) %>%
  left_join(MEASURES, by = c("measure" = "column")) %>%
  mutate(sample = if_else(sample == "common", "common (every measure present)", "each measure's own (not comparable)"),
         site = site_name) %>%
  relocate(sample, outcome, form, family, label)

message("\nCross-validated AUC on the common sample, spline form (difference from VT/PBW with a paired bootstrap interval):")
print(as.data.frame(results %>% filter(startsWith(sample, "common"), form == "spline") %>%
                      arrange(outcome, desc(auc)) %>%
                      select(outcome, label, auc, auc_lo, auc_hi, delta_auc_vs_vtpbw, delta_auc_lo, delta_auc_hi,
                             delta_aic_vs_vtpbw, brier, n_patients, n_deaths) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
write_csv(mask_small_counts(results), file.path(final_dir, paste0("mortality_prediction_", site_name, ".csv")))

# =============================================================================
# Figure: cross-validated AUC and its difference from VT/PBW, common sample
# =============================================================================
figure_rows <- results %>% filter(startsWith(sample, "common"), form == "spline") %>%
  mutate(label = factor(label, levels = rev(MEASURES$label)))
auc_panel <- ggplot(figure_rows, aes(auc, label, colour = family)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(xmin = auc_lo, xmax = auc_hi)) +
  facet_wrap(~ outcome) +
  scale_colour_manual(values = FAMILY_COLOURS, name = NULL) +
  labs(title = "Cross-validated AUC of each measure alone",
       subtitle = sprintf("no covariates; 3-df spline of the log measure; %d-fold cross-validation; the same patients for every measure", N_FOLDS),
       x = "AUC (95% bootstrap interval)", y = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
delta_panel <- ggplot(figure_rows, aes(delta_auc_vs_vtpbw, label, colour = family)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(xmin = delta_auc_lo, xmax = delta_auc_hi)) +
  facet_wrap(~ outcome) +
  scale_colour_manual(values = FAMILY_COLOURS, guide = "none") +
  labs(title = "Difference from VT/PBW", x = "AUC minus the AUC of VT/PBW (paired bootstrap)", y = NULL) +
  theme_minimal(base_size = 10)
ggsave(file.path(final_dir, paste0("mortality_prediction_", site_name, ".pdf")),
       patchwork::wrap_plots(auc_panel, delta_panel, ncol = 1), width = 11, height = 10)
message("xsec_mortality_prediction complete -> ", final_dir)
