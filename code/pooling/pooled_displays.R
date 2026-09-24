# =============================================================================
# Pooled displays: figures 1C, 2 and 3C, and the ratio's height channel
# =============================================================================
# Run centrally, like pooled_estimates.R and pooled_biotrauma.R. Two of the displays
# need no patient data at all; the other two pool tables every site returns.
#
#   Figure 1C  the counterfactual dosing grid (formulas only): the tidal volume
#              6 mL/kg PBW prescribes against the volume that delivers 11% of
#              predicted FVC, by age, sex, race and height
#   Claim 2    the ratio's height channel (formulas only): log PBW/PFVC across height
#              by sex at a fixed age. GLI FVC is a power law in height and Devine PBW a
#              line with an intercept, so the log ratio is an inverted U whose peak
#              differs by sex
#   Figure 2   delivered strain inside the 6-8 mL/kg band: the pooled VT/PFVC
#              distribution by sex, race, age and height (dose_vtpfvc_histograms_*,
#              04 section 4l), and the variance of log VT/PFVC split into the
#              clinician's dose and the label's mis-sizing, pooled exactly from each
#              site's moments (dose_variance_decomposition_*)
#   Figure 3C  the worked example (crs_channels_estimates_*, supplement/xsec_crs_channels.R
#              section 5): a Black woman of 160 cm against a white man of 180 cm, both 60,
#              by Devine PBW, GLI PFVC at 60 and at 25, and measured compliance
#
# The formulas mirror script 03: Devine PBW = 50 (men) or 45.5 (women) + 2.3 kg per
# inch over 60; race-specific GLI-2012 FVC (rspiro::pred_GLI, White = 1, Black = 2).
# Colours: viridis for age (continuous), Okabe-Ito for sex, race and sites.
#
# Usage: PBWPFVC_RESULTS_ROOT=/path/to/results Rscript code/pooling/pooled_displays.R
#   (one subfolder per site holding its final/ tables; output to <root>/All sites/)
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(here); library(patchwork); library(rspiro) })

root <- Sys.getenv("PBWPFVC_RESULTS_ROOT", here("results"))
if (!dir.exists(root)) stop("results root not found: ", root)
out_dir <- file.path(root, "All sites"); dir.create(out_dir, showWarnings = FALSE)
sites <- setdiff(list.dirs(root, recursive = FALSE, full.names = FALSE), "All sites")
sites <- sites[!startsWith(sites, ".")]
source(here("utils", "site_anonymization.R"))
site_aliases <- build_site_aliases(file.path(root, sites))
print_site_alias_key(site_aliases)

OKABE_ITO <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
SEX_COLOURS  <- c(Male = "#0072B2", Female = "#D55E00")
VTPFVC_TARGET <- 11   # % of predicted FVC: ARMA's low-VT arm near its 75th percentile
theme_set(theme_minimal(base_size = 10))

# every table a site returns, whichever block subfolder it sits in
read_sites <- function(prefix) {
  map_dfr(sites, function(s) {
    files <- unlist(map(c("", "final", "cross_sectional", "supplement", "final/cross_sectional", "final/supplement"),
                        ~ Sys.glob(file.path(root, s, .x, paste0(prefix, "*.csv")))))
    if (!length(files)) return(NULL)
    read_csv(files[1], show_col_types = FALSE) %>% mutate(site = anonymize_site(s, site_aliases$aliases))
  })
}
devine_pbw <- function(height_cm, sex) if_else(sex == "Male", 50, 45.5) + 2.3 * (height_cm / 2.54 - 60)
gli_fvc <- function(age, height_cm, sex, race) pred_GLI(age = age, height = height_cm / 100,
                                                       gender = if_else(sex == "Male", 1, 2),
                                                       ethnicity = if_else(race == "Black", 2, 1), param = "FVC")

# =============================================================================
# Figure 1C: the dosing grid (no patient data)
# =============================================================================
dosing_grid <- expand_grid(age = seq(40, 90, by = 10), height_cm = seq(150, 200, by = 5),
                           sex = c("Female", "Male"), race = c("White", "Black")) %>%
  mutate(pbw = devine_pbw(height_cm, sex), pfvc_l = gli_fvc(age, height_cm, sex, race),
         vt_pbw_ml = 6 * pbw, vt_strain_ml = VTPFVC_TARGET / 100 * pfvc_l * 1000,
         vtpfvc_at_6_ml_kg = vt_pbw_ml / (pfvc_l * 1000) * 100,
         excess_ml = vt_pbw_ml - vt_strain_ml)
write_csv(dosing_grid, file.path(out_dir, "figure1c_dosing_grid.csv"))
fig_1c <- ggplot(dosing_grid, aes(height_cm, vtpfvc_at_6_ml_kg, colour = age, group = age)) +
  geom_hline(yintercept = VTPFVC_TARGET, linetype = 2, colour = "grey40") +
  geom_line(linewidth = 0.8) +
  facet_grid(race ~ sex) +
  scale_colour_viridis_c(name = "Age (years)") +
  labs(title = "C. The strain 6 mL/kg PBW delivers, by age, sex, race and height",
       subtitle = sprintf("VT/PFVC (%% of predicted FVC) at 6 mL/kg PBW; dashed line: %d%%. Devine PBW, race-specific GLI-2012 FVC", VTPFVC_TARGET),
       x = "Height (cm)", y = "VT/PFVC at 6 mL/kg PBW (%)")
ggsave(file.path(out_dir, "figure1c_dosing_grid.pdf"), fig_1c, width = 8, height = 6)

# =============================================================================
# Claim 2: the ratio's height channel (no patient data)
# =============================================================================
HEIGHT_CHANNEL_AGE <- 60
height_channel <- expand_grid(height_cm = seq(150, 200, by = 0.5), sex = c("Female", "Male"), race = "White") %>%
  mutate(log_ratio = log(devine_pbw(height_cm, sex) / gli_fvc(HEIGHT_CHANNEL_AGE, height_cm, sex, race))) %>%
  group_by(sex) %>%
  mutate(log_ratio_vs_peak = log_ratio - max(log_ratio), peak_height_cm = height_cm[which.max(log_ratio)],
         range_150_200 = max(log_ratio) - min(log_ratio)) %>% ungroup()
write_csv(height_channel, file.path(out_dir, "claim2_height_channel.csv"))
message("Height channel at age ", HEIGHT_CHANNEL_AGE, ": peak ",
        paste(distinct(height_channel, sex, peak_height_cm) %>% transmute(txt = sprintf("%s %.1f cm", sex, peak_height_cm)) %>% pull(txt), collapse = ", "),
        "; range 150-200 cm ",
        paste(distinct(height_channel, sex, range_150_200) %>% transmute(txt = sprintf("%s %.3f", sex, range_150_200)) %>% pull(txt), collapse = ", "))
fig_height <- ggplot(height_channel, aes(height_cm, log_ratio_vs_peak, colour = sex)) +
  geom_line(linewidth = 0.9) +
  geom_vline(data = distinct(height_channel, sex, peak_height_cm), aes(xintercept = peak_height_cm, colour = sex), linetype = 3) +
  scale_colour_manual(values = SEX_COLOURS, name = NULL) +
  labs(title = "The ratio's height channel is small, curved and sex-reversed",
       subtitle = sprintf("log PBW/PFVC relative to its peak, white patients aged %d; dotted: the peak", HEIGHT_CHANNEL_AGE),
       x = "Height (cm)", y = "log PBW/PFVC minus its maximum")
ggsave(file.path(out_dir, "claim2_height_channel.pdf"), fig_height, width = 7, height = 4.5)

# =============================================================================
# Figure 2: delivered strain inside the band (pooled from the sites)
# =============================================================================
dose_hist <- read_sites("dose_vtpfvc_histograms_")
dose_moments <- read_sites("dose_variance_decomposition_")
if (nrow(dose_hist) && nrow(dose_moments)) {
  pooled_hist <- dose_hist %>% group_by(group_type, group_value, bin_left, bin_right) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    group_by(group_type, group_value) %>%
    mutate(n = sum(count), share = count / n, mid = (bin_left + bin_right) / 2) %>% ungroup()
  over_target <- pooled_hist %>% group_by(group_type, group_value, n) %>%
    summarise(pct_over_target = 100 * sum(count[bin_left >= VTPFVC_TARGET]) / first(n),
              median = mid[which(cumsum(count) >= first(n) / 2)[1]], .groups = "drop")
  write_csv(over_target, file.path(out_dir, "figure2a_vtpfvc_by_group.csv"))
  group_order <- c(sex = "Sex", race = "Race", age_bin = "Age (years)", height_bin = "Height (cm)")
  # one panel per grouping, each with its own legend: the groupings share no levels
  group_panel <- function(g) {
    d <- pooled_hist %>% filter(group_type == g)
    levels_in_order <- unique(d$group_value[order(suppressWarnings(as.numeric(str_extract(d$group_value, "[0-9]+"))), d$group_value)])
    ggplot(d %>% mutate(group_value = factor(group_value, levels_in_order)), aes(mid, share, colour = group_value)) +
      geom_vline(xintercept = VTPFVC_TARGET, linetype = 2, colour = "grey40") +
      geom_line(linewidth = 0.7) +
      scale_colour_manual(values = OKABE_ITO, name = NULL) +
      labs(subtitle = group_order[[g]], x = "VT/PFVC (%)", y = if (g == "sex") "share per 0.5-point bin" else NULL) +
      theme(legend.position = "bottom") + guides(colour = guide_legend(ncol = 2))
  }
  fig_2a <- wrap_plots(map(names(group_order), group_panel), nrow = 1) +
    plot_annotation(title = "A. Delivered VT/PFVC inside the 6-8 mL/kg band",
                    subtitle = sprintf("pooled across %d site(s); dashed line: %d%% of predicted FVC",
                                       n_distinct(dose_hist$site), VTPFVC_TARGET))

  # exact pooled variance and covariance from each site's moments: within-site sums of
  # squares plus the spread of the site means around the grand mean
  pool_moments <- function(d) {
    N <- sum(d$n)
    grand_vtpbw <- sum(d$n * d$mean_log_vtpbw) / N; grand_ratio <- sum(d$n * d$mean_log_ratio) / N
    var_vtpbw <- (sum((d$n - 1) * d$var_log_vtpbw) + sum(d$n * (d$mean_log_vtpbw - grand_vtpbw)^2)) / (N - 1)
    var_ratio <- (sum((d$n - 1) * d$var_log_ratio) + sum(d$n * (d$mean_log_ratio - grand_ratio)^2)) / (N - 1)
    cov_vr <- (sum((d$n - 1) * d$cov_log_vtpbw_ratio) +
                 sum(d$n * (d$mean_log_vtpbw - grand_vtpbw) * (d$mean_log_ratio - grand_ratio))) / (N - 1)
    total <- var_vtpbw + var_ratio + 2 * cov_vr
    tibble(n = N, var_log_vtpfvc = total, share_clinician = var_vtpbw / total,
           share_missizing = var_ratio / total, share_covariance = 2 * cov_vr / total)
  }
  decomposition <- bind_rows(
    dose_moments %>% select(site, group_type, group_value, n, var_log_vtpfvc, share_clinician, share_missizing, share_covariance),
    dose_moments %>% group_by(group_type, group_value) %>% group_modify(~ pool_moments(.x)) %>% ungroup() %>%
      mutate(site = "Pooled"))
  write_csv(decomposition, file.path(out_dir, "figure2b_variance_decomposition.csv"))
  fig_2b <- decomposition %>% filter(group_type == "overall") %>%
    pivot_longer(starts_with("share_"), names_to = "component", values_to = "share") %>%
    mutate(component = factor(recode(component, share_missizing = "Mis-sizing (PBW/PFVC)",
                                      share_clinician = "Clinician's dose (VT/PBW)", share_covariance = "Covariance"),
                              c("Mis-sizing (PBW/PFVC)", "Clinician's dose (VT/PBW)", "Covariance")),
           site = factor(site, c("Pooled", sort(setdiff(unique(site), "Pooled"))))) %>%
    ggplot(aes(share, fct_rev(site), fill = component)) +
    geom_col(width = 0.6) +
    geom_vline(xintercept = 0, colour = "grey30") +
    scale_fill_manual(values = c(`Mis-sizing (PBW/PFVC)` = "#D55E00", `Clinician's dose (VT/PBW)` = "#0072B2",
                                 Covariance = "#999999"), name = NULL) +
    scale_x_continuous(labels = scales::percent) +
    labs(title = "B. Where the variance in delivered strain comes from",
         subtitle = "Var(log VT/PFVC) = Var(log VT/PBW) + Var(log PBW/PFVC) + 2 Cov; inside the band",
         x = "share of the variance of log VT/PFVC", y = NULL) +
    theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "figure2_delivered_strain.pdf"), wrap_elements(fig_2a) / fig_2b + plot_layout(heights = c(1.4, 1)),
         width = 13, height = 8.5)
  message("Figure 2 -> ", file.path(out_dir, "figure2_delivered_strain.pdf"))
} else message("Figure 2 skipped: no dose_vtpfvc_histograms_ or dose_variance_decomposition_ tables (04 section 4l)")

# =============================================================================
# Figure 3C: the worked example (formula ratios; measured compliance pooled)
# =============================================================================
crs <- read_sites("crs_channels_estimates_")
if (nrow(crs) && any(crs$model == "worked example")) {
  worked <- crs %>% filter(model == "worked example", sample == "all plateau-measured")
  formula_rows <- worked %>% filter(is.na(se)) %>% distinct(term, estimate)
  contrast <- worked %>% filter(startsWith(term, "log Crs ratio"), !is.na(se))
  w <- 1 / contrast$se^2
  pooled_contrast <- tibble(estimate = sum(w * contrast$estimate) / sum(w), se = sqrt(1 / sum(w)))
  bars <- bind_rows(
    formula_rows %>% transmute(measure = recode(term, "log PBW ratio, A / B (Devine)" = "PBW (Devine)",
                                                "log PFVC ratio, A / B (GLI at the profiles' age)" = "PFVC (GLI, at 60)",
                                                "log PFVC ratio, A / B (GLI at age 25)" = "Age-standardised PFVC (GLI at 25)"),
                               ratio = exp(estimate), lo = NA_real_, hi = NA_real_, source = "formula"),
    tibble(measure = "Measured compliance (pooled)", ratio = exp(pooled_contrast$estimate),
           lo = exp(pooled_contrast$estimate - 1.96 * pooled_contrast$se),
           hi = exp(pooled_contrast$estimate + 1.96 * pooled_contrast$se), source = "patients"),
    contrast %>% transmute(measure = paste0("Measured compliance, ", site), ratio = exp(estimate),
                           lo = exp(lo), hi = exp(hi), source = "patients"))
  write_csv(bars, file.path(out_dir, "figure3c_worked_example.csv"))
  fig_3c <- ggplot(bars, aes(ratio, fct_rev(fct_inorder(measure)), colour = source)) +
    geom_vline(xintercept = 1, colour = "grey30") +
    geom_point(size = 2.8) + geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, orientation = "y", na.rm = TRUE) +
    scale_colour_manual(values = c(formula = "#999999", patients = "#0072B2"), name = NULL) +
    labs(title = "C. A Black woman of 160 cm against a white man of 180 cm, both 60",
         subtitle = "her size as a fraction of his: each formula, and measured compliance\n(model contrast at each site's median covariates)",
         x = "ratio, woman / man", y = NULL) +
    theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "figure3c_worked_example.pdf"), fig_3c, width = 8, height = 3.5 + 0.25 * nrow(contrast))
  message("Figure 3C -> ", file.path(out_dir, "figure3c_worked_example.pdf"))
} else message("Figure 3C skipped: no worked-example rows in crs_channels_estimates_ (xsec_crs_channels.R section 5)")

message("pooled_displays complete -> ", out_dir)
