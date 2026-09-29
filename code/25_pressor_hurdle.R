# =============================================================================
# Script 25 (pressor hurdle): the vasopressor dose with its zero days kept
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
#
# Figure 4 reads vasopressors as two separate fits: any_pressor, a logistic model
# of whether a pressor runs, and pressor_dose, a Gaussian model of the log dose on
# the days one does. The second conditions on being on a pressor, which is a
# post-baseline state the exposure may itself affect, so the on-pressor days of
# small- and large-lung patients need not be comparable populations.
#
# Folding the zero days in as log(dose + c) does not fix it: with about half the
# patient-days at zero the spike at log(c) carries most of the variance, the slope
# moves with the analyst's c, and one coefficient then mixes whether a pressor runs
# with how much runs. (That model exists as the ne_equiv_peak marker in 22 and is
# not in the figure.)
#
# This script fits the parameterization that keeps the zeros without inventing a
# constant: a HURDLE LOG-NORMAL model (GLMMadaptive), one fit holding
#
#   zero part      the log-odds that a pressor runs at all that day
#   positive part  the log dose on the days one does, log-normal
#
# with their own random effects, and from them the quantity the zero days are
# wanted for:
#
#   marginal       the effect on the EXPECTED DAILY DOSE, zero days included.
#                  For a hurdle log-normal, E[Y] = P(Y > 0) x exp(mu + sigma^2/2),
#                  so d log E[Y] / dx = (1 - p) * beta_on + beta_positive, read at
#                  the sample's mean p. Both terms and their covariance come from
#                  the one fit, so the standard error is a delta-method one and not
#                  an assumption that the parts are independent.
#
# Two things to know before this is pooled or read across sites. The marginal takes p
# at the cohort's own mean, so (1 - p) weights the zero part differently where
# pressors are common than where they are rare: two sites' marginals are the same
# quantity only if their p is similar, and mean_p_on is in the table so that can be
# checked. And GLMMadaptive's hurdle families model the probability of a ZERO, so the
# zero part's coefficients are negated here to read as the odds that a pressor runs.
# That convention was verified against these data rather than taken from the manual:
# sicker patients are on pressors far more often (1.5% to 21.6% of days across SOFA
# tertiles) and the zero part's SOFA coefficient is correspondingly negative.
#
# Longitudinal only, by design and not by limitation: on this panel the joint
# model's death and extubation correction moves every lab divergence by under a
# fifth of a standard error (jm_lme_check_*), and no hurdle family is available
# inside JMbayes2. What this cannot do is correct for patients leaving as they are
# weaned, which is the same limit the dose part has today.
#
# Writes final/pressor_hurdle_{tag}.csv: each part's level and divergence per SD of
# log PFVC, adjusted and unadjusted, in the units the pooling converts.
#
# Usage:  Rscript code/25_pressor_hurdle.R
#         PBWPFVC_COHORT=nosupport PBWPFVC_JM_CLOCK=icu Rscript code/25_pressor_hurdle.R
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(arrow); library(splines); library(GLMMadaptive)
})
rm(list = ls())
source("utils/config.R")

site_name  <- config$site_name
output_dir <- config$output_dir
final_dir  <- final_dir_for("injury")
source(here("code", "20_biotrauma_grid.R"))   # JM_GRID, JM_HORIZON, h_suffix, panel_path, clock_tag

out_tag <- paste0(clock_tag, h_suffix, "_", site_name)
message("=== 25_pressor_hurdle: ", site_name, ", cohort ", config$cohort, ", clock ", JM_CLOCK)

long <- read_parquet(panel_path("long"))
surv <- read_parquet(panel_path("surv"))

# The same frame the joint models use for the pressor outcome, with the zero days
# kept: every panel day of every patient, the dose being 0 when none runs.
d <- long %>%
  select(hospitalization_id, vent_day, ne_equiv_peak, l_pressor) %>%
  inner_join(surv %>% select(hospitalization_id, np_sofa, sf_0, age10, sex_category,
                             race_category, log_pfvc_sd, ne_equiv_0),
             by = "hospitalization_id") %>%
  # the previous day's pressor state is missing on the first day of each trajectory,
  # as it is for every marker in 22: those rows leave here too, so the two agree
  filter(!is.na(ne_equiv_peak), !is.na(np_sofa), !is.na(sf_0), !is.na(log_pfvc_sd),
         !is.na(l_pressor)) %>%
  mutate(log_sf_0 = as.numeric(scale(log(sf_0))), np_sofa = as.numeric(scale(np_sofa)),
         # the day-0 dose enters as two terms, as it does in 22: a pressor running
         # at day 0, and the log dose when one was
         on_y0  = as.numeric(ne_equiv_0 > 0),
         log_y0 = if_else(ne_equiv_0 > 0, log(ne_equiv_0), 0),
         y = ne_equiv_peak) %>%
  group_by(hospitalization_id) %>% filter(n() >= 2) %>% ungroup()
p_on <- mean(d$y > 0)
message("Frame: ", nrow(d), " patient-days for ", n_distinct(d$hospitalization_id),
        " patients; a pressor runs on ", sum(d$y > 0), " of them (", round(100 * p_on), "%)")
if (n_distinct(d$hospitalization_id) < 50 || sum(d$y > 0) < 50) {
  message("too few pressor days to fit a hurdle model; nothing written"); quit(save = "no", status = 0)
}

DEMO    <- "ns(age10, 4) + sex_category + race_category"
DEMO_ZI <- DEMO   # the same age spline in both parts, as everywhere else
# The previous day's pressor state belongs in the dose part, where it is the lag the
# joint models use. It cannot go in the zero part: a pressor running yesterday almost
# determines one running today, the starting logistic fit separates, and the hurdle
# optimiser then diverges. The zero part keeps the baseline state instead.
BASE    <- "np_sofa + log_sf_0 + on_y0 + log_y0 + l_pressor"
BASE_ZI <- "np_sofa + log_sf_0 + on_y0"
# A cohort where nobody who was off pressors at the baseline ever starts one separates
# the zero part completely: on_y0 predicts it perfectly, its coefficient runs off, and
# the fit fails. The term is dropped in that case and the table records it, so a
# specification chosen by the data is visible rather than assumed.
zi_separates <- with(d, {
  tab <- table(factor(on_y0, 0:1), factor(y > 0, c(FALSE, TRUE)))
  any(rowSums(tab) > 0 & (tab[, 1] == 0 | tab[, 2] == 0))
})
if (zi_separates) {
  message("The baseline pressor state separates the zero part completely ",
          "(no patient crosses it); it is dropped from that part")
  BASE_ZI <- "np_sofa + log_sf_0"
}
# The positive part is fitted on the pressor days alone, where a baseline term can be
# constant -- in a cohort whose pressors all start on day 0, every positive day has
# on_y0 = 1 and the term is the intercept again. Any such term is dropped, and named.
constant_on_positive <- d %>% filter(y > 0) %>%
  summarise(across(c(on_y0, log_y0, l_pressor), ~ n_distinct(.x) < 2)) %>%
  unlist() %>% (\(x) names(x)[x])
if (length(constant_on_positive)) {
  message("Constant on the pressor days, so dropped from the dose part: ",
          paste(constant_on_positive, collapse = ", "))
  BASE <- paste(setdiff(c("np_sofa", "log_sf_0", "on_y0", "log_y0", "l_pressor"),
                        constant_on_positive), collapse = " + ")
}
EXPOSURE <- "log_pfvc_sd"
SIZE <- paste(EXPOSURE, "+", EXPOSURE, ": vent_day")
ZI_RANDOM <- identical(Sys.getenv("PBWPFVC_HURDLE_ZI_RANDOM", "0"), "1")
switchers <- d %>% group_by(hospitalization_id) %>%
  summarise(on = sum(y > 0), days = n(), .groups = "drop") %>%
  summarise(switch = sum(on > 0 & on < days), never = sum(on == 0), always = sum(on == days))
message("Pressor pattern: ", switchers$switch, " patients switch on or off, ",
        switchers$never, " never on, ", switchers$always, " always on",
        "; zero part random intercept: ", if (ZI_RANDOM) "on" else "off")

rows <- list()
for (adjusted in c(TRUE, FALSE)) {
  # The dose part keeps the day spline; the zero part takes time linearly, which costs
  # nothing since the term it reports is the exposure by day interaction, linear in
  # either case. Both parts use the age spline, the one age form used everywhere.
  # The zero part's age was briefly linear, a choice made on the synthetic panel, whose
  # pressors run on 9% of days and never start after day 0. MIMIC (8,735 pressor days
  # of 17,019) identified the model with that linear term; the spline there is
  # untested until the next run, and a failure stops the script rather than falling
  # back to the line.
  rhs    <- paste("ns(vent_day, 3) +", SIZE, "+", BASE,    if (adjusted) paste("+", DEMO) else "")
  rhs_zi <- paste("vent_day +", SIZE, "+", BASE_ZI,        if (adjusted) paste("+", DEMO_ZI) else "")
  adj_label <- if (adjusted) "adjusted" else "unadjusted"
  # Whether the zero part carries its own patient random intercept (ZI_RANDOM).
  # It is off by default: most patients are never on a pressor at all, and a
  # patient-level intercept in a binary part whose clusters are almost all constant
  # is not identified -- the intercept runs off and the optimiser reports a large
  # coefficient. Where the switching is dense enough to support it, set
  # PBWPFVC_HURDLE_ZI_RANDOM=1 and the table records which was used.
  hurdle_fit <- function(...)
    mixed_model(fixed = as.formula(paste("y ~", rhs)), random = ~ vent_day | hospitalization_id,
                zi_fixed = as.formula(paste("~", rhs_zi)),
                zi_random = if (ZI_RANDOM) ~ 1 | hospitalization_id else NULL,
                data = d, family = hurdle.lognormal(), ...)
  fit <- tryCatch(hurdle_fit(),
    error = function(e) {
      message("  ", adj_label, ": diverged from the EM start; refitting with iter_EM = 0")
      tryCatch(hurdle_fit(control = list(iter_EM = 0)),
               error = function(e2) { message("  ", adj_label, ": hurdle fit failed: ",
                                              conditionMessage(e2)); NULL })
    })
  if (is.null(fit)) next

  # GLMMadaptive's hurdle families model the probability of a ZERO in the zi part,
  # so its coefficients are flipped to read as the odds that a pressor runs.
  fe_pos <- fixef(fit); fe_zi <- fixef(fit, sub_model = "zero_part")
  # A fit can converge and still leave a singular Hessian, which means the design is
  # not identified on this cohort: said so, rather than reported without intervals.
  V <- tryCatch(vcov(fit), error = function(e) {
    message("  ", adj_label, ": the fit converged but its covariance is singular, so the ",
            "design is not identified on this cohort; nothing reported for it")
    NULL })
  if (is.null(V)) next
  # The two parts share one covariance matrix, whose rows run: the dose part's
  # coefficients, then the random-effects variances (D_11, ...), then the zero part's
  # with a zi_ prefix. The variances sit BETWEEN the two blocks, so a row cannot be
  # found by counting along from the start -- each is matched by name.
  # R names an interaction in the order its variables first appear in the formula,
  # so the same term is log_pfvc_sd:vent_day in the dose part, where a spline comes
  # first, and vent_day:log_pfvc_sd in the zero part, where vent_day does. Each
  # part's term is therefore found by its components, never by the other part's
  # spelling. A term that cannot be found is an error, not a skipped row: the first
  # MIMIC run lost its divergence rows to exactly this, with only a log line to say so.
  same_term <- function(a, b) identical(sort(strsplit(a, ":", fixed = TRUE)[[1]]),
                                        sort(strsplit(b, ":", fixed = TRUE)[[1]]))
  find_name <- function(names_in_part, wanted) {
    hit <- names_in_part[vapply(names_in_part, same_term, logical(1), b = wanted)]
    if (length(hit)) hit[1] else NA_character_
  }
  TERMS <- c("level" = EXPOSURE, "divergence per day" = paste0(EXPOSURE, ":vent_day"))

  for (which_term in names(TERMS)) {
    nm_pos <- find_name(names(fe_pos), TERMS[[which_term]])
    nm_zi  <- find_name(names(fe_zi),  TERMS[[which_term]])
    i_pos  <- if (is.na(nm_pos)) NA_integer_ else match(nm_pos, rownames(V))
    i_zi   <- if (is.na(nm_zi))  NA_integer_ else match(paste0("zi_", nm_zi), rownames(V))
    if (anyNA(c(i_pos, i_zi)))
      stop(adj_label, ": the ", which_term, " term (", TERMS[[which_term]], ") was not found in ",
           if (is.na(i_pos)) "the dose part" else "the zero part", "; its names are: ",
           paste(if (is.na(i_pos)) names(fe_pos) else names(fe_zi), collapse = ", "))
    b_on  <- -unname(fe_zi[[nm_zi]])         # log-odds that a pressor RUNS
    b_pos <-  unname(fe_pos[[nm_pos]])
    se_on <- sqrt(V[i_zi, i_zi]); se_pos <- sqrt(V[i_pos, i_pos])
    # d log E[Y] / dx = (1 - p) * b_on + b_pos, with the covariance of the two parts
    # from the one fit; the sign flip makes the zero part's covariance negative
    g <- c(-(1 - p_on), 1)                   # gradient on (zi coefficient, positive coefficient)
    Vsub <- V[c(i_zi, i_pos), c(i_zi, i_pos)]
    marg <- (1 - p_on) * b_on + b_pos
    se_marg <- sqrt(drop(t(g) %*% Vsub %*% g))
    # the loop variable is not called `part`: tibble() evaluates its arguments in
    # order and each sees the columns already made, so a column named part would
    # shadow it from the second argument on
    for (piece in list(list("zero part (odds a pressor runs)", b_on, se_on),
                       list("positive part (log dose when on)", b_pos, se_pos),
                       list("marginal (expected daily dose, zeros included)", marg, se_marg)))
      rows[[length(rows) + 1L]] <- tibble(
        part = piece[[1]], term = which_term, adjustment = adj_label,
        estimate = piece[[2]], se = piece[[3]],
        lo = piece[[2]] - 1.96 * piece[[3]], hi = piece[[2]] + 1.96 * piece[[3]],
        p = 2 * pnorm(-abs(piece[[2]] / piece[[3]])),
        n_patients = n_distinct(d$hospitalization_id), n_obs = nrow(d),
        n_pressor_days = sum(d$y > 0), mean_p_on = p_on,
        zero_part_baseline = !zi_separates, zero_part_random_intercept = ZI_RANDOM)
  }
}

if (!length(rows)) { message("no hurdle fit produced estimates; nothing written"); quit(save = "no", status = 0) }
out <- bind_rows(rows) %>%
  mutate(exposure = EXPOSURE,
         unit = "per SD of log PFVC (log-odds for the zero part, log dose otherwise)",
         cohort = config$cohort, clock = JM_CLOCK, panel = h_suffix, site = site_name)
write_csv(out, file.path(final_dir, paste0("pressor_hurdle_", out_tag, ".csv")))
message("pressor hurdle -> ", file.path(final_dir, paste0("pressor_hurdle_", out_tag, ".csv")))
print(as.data.frame(out %>% filter(adjustment == "adjusted") %>%
                      transmute(part, term, estimate = signif(estimate, 3),
                                lo = signif(lo, 3), hi = signif(hi, 3), p = signif(p, 3))), row.names = FALSE)
message("=== 25_pressor_hurdle complete")
