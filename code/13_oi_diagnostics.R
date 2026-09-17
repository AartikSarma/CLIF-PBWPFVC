# =============================================================================
# Script 13 (OI diagnostics): is an oxygenation-index signal the lung, or the
# ventilator setting in its own numerator?
# =============================================================================
# The oxygenation indices are worth having because they credit the SUPPORT
# needed to reach an oxygenation, which the SF ratio ignores. That same property
# is their hazard here. Both are
#
#   index = 100 x mean airway pressure / ratio      (ratio = SF for OSI, P/F for OI)
#
# and mean airway pressure is a ventilator setting that the exposure moves
# arithmetically: a larger tidal volume at the same PEEP and compliance raises
# it. So an index that worsens with VT/PFVC may be reporting the numerator, with
# nothing happening in the lung.
#
# In logs the index separates exactly:
#
#   log index = log 100 + log(mean airway pressure) - log(ratio)
#
# so the same contrast fitted on all three outcomes, over identical rows with an
# identical design, decomposes exactly:
#
#   contrast(log index) = contrast(log mean airway pressure) - contrast(log ratio)
#
# Read it as: the first term is the pressure the exposure buys, the second is the
# oxygenation it buys. An index effect carried by the first term is arithmetic,
# and belongs in the paper as a mechanics result, not as biotrauma. One carried
# by the second term is oxygenation, and the index is telling you what SF would
# have told you. One carried by both, in the same direction, is the interesting
# case: more pressure AND worse gas exchange per unit of it.
#
# The fits are ordinary least squares with patient-clustered standard errors,
# not mixed models, precisely so the decomposition is an algebraic identity
# rather than an approximation (a mixed model estimates its variance components
# per outcome, which breaks the identity by a little and invites the reader to
# wonder by how much).
#
# Also written: coverage by day for each ingredient (the constraint on OI), and
# the agreement between the two indices where both exist.
#
# Needs a panel carrying map_aw / oi / osi (13_biotrauma_panel.R, any grid):
#   PBWPFVC_JM_HORIZON_H=72 Rscript code/13_biotrauma_panel.R
# Usage:
#   PBWPFVC_JM_HORIZON_H=72 Rscript code/13_oi_diagnostics.R
#
# Outputs to final/ (aggregates only, every reported cell >= MIN_CELL patients):
#   oi_decomposition_{site}.csv   the contrast on the index, the numerator and
#       the ratio, per index x adjustment x horizon, with the identity check
#   oi_coverage_{site}.csv        patients per day with each ingredient
#   oi_agreement_{site}.csv       OI against OSI where both exist, binned
#   oi_diagnostics_{site}.pdf     four pages
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(here)
                                 library(splines); library(patchwork) })
rm(list = ls())
source("utils/config.R")
site_name  <- config$site_name
output_dir <- here("output", paste0(site_name, "_output"), "intermediate")
final_dir  <- here("output", paste0(site_name, "_output"), "final")
dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)
if (!nzchar(Sys.getenv("PBWPFVC_JM_HORIZON_H"))) Sys.setenv(PBWPFVC_JM_HORIZON_H = "72")
source(here("code", "13_biotrauma_grid.R"))

MIN_CELL <- 10L
HOURS    <- as.numeric(strsplit(Sys.getenv("PBWPFVC_OI_HORIZONS_H", "24,48,72"), ",")[[1]])
HOURS    <- HOURS[HOURS <= JM_HORIZON * 24]
okabe    <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#000000")
theme_set(theme_minimal(base_size = 10))

panel_path <- file.path(output_dir, paste0("jm_long_", h_suffix, ".parquet"))
if (!file.exists(panel_path))
  stop("no ", h_suffix, " panel: run  PBWPFVC_JM_HORIZON_H=", as.integer(JM_HORIZON * 24),
       " Rscript code/13_biotrauma_panel.R")
long <- read_parquet(panel_path)
surv <- read_parquet(file.path(output_dir, paste0("jm_surv_", h_suffix, ".parquet")))
if (!all(c("map_aw", "oi", "osi") %in% names(long)))
  stop("this panel predates the oxygenation indices: rebuild it with the current 13_biotrauma_panel.R")

d_all <- long %>%
  filter(period >= 1L, !is.na(l_pressor)) %>%
  inner_join(surv %>% select(hospitalization_id, np_sofa, bmi, age10, sex_category, race_category,
                             vtpfvc_c, vtpbw_pt_mean, osi_0, oi_0),
             by = "hospitalization_id") %>%
  filter(!is.na(np_sofa), !is.na(vtpfvc_c)) %>%
  mutate(pf = if_else(is.finite(oi) & oi > 0, 100 * map_aw / oi, NA_real_))   # the P/F the OI implies
message("=== 13_oi_diagnostics (", h_suffix, " panel): ", nrow(d_all), " patient-periods, ",
        n_distinct(d_all$hospitalization_id), " patients ===")

# =============================================================================
# 1. coverage of each ingredient, by day
# =============================================================================
days <- sort(unique(round(d_all$vent_day / STEP) * STEP))
coverage <- d_all %>%
  mutate(day = round(vent_day / STEP) * STEP) %>%
  group_by(day) %>%
  summarise(n_patients = n_distinct(hospitalization_id),
            n_map_aw = n_distinct(hospitalization_id[is.finite(map_aw)]),
            n_sf     = n_distinct(hospitalization_id[is.finite(sf)]),
            n_pf     = n_distinct(hospitalization_id[is.finite(pf)]),
            n_osi    = n_distinct(hospitalization_id[is.finite(osi)]),
            n_oi     = n_distinct(hospitalization_id[is.finite(oi)]), .groups = "drop") %>%
  mutate(across(starts_with("n_") & !all_of("n_patients"),
                ~ if_else(n_patients >= MIN_CELL, ., NA_integer_)),
         site = site_name)
write_csv(coverage, file.path(final_dir, paste0("oi_coverage_", site_name, ".csv")))

# =============================================================================
# 2. the exact decomposition
# =============================================================================
# patient-clustered sandwich: the rows of one patient are one observation
cluster_vcov <- function(fit, cluster) {
  X <- model.matrix(fit); u <- residuals(fit)
  bread <- solve(crossprod(X))
  meat <- matrix(0, ncol(X), ncol(X))
  for (g in split(seq_along(u), cluster)) {
    s <- crossprod(X[g, , drop = FALSE], u[g])
    meat <- meat + tcrossprod(s)
  }
  G <- n_distinct(cluster); n <- nrow(X); k <- ncol(X)
  bread %*% meat %*% bread * (G / (G - 1)) * ((n - 1) / (n - k))
}
# the contrast at hour hh: level + divergence x days, per point of VT/PFVC
contrast_at <- function(fit, V, hh) {
  b <- coef(fit)
  tn <- intersect(c("vtpfvc_c:vent_day", "vent_day:vtpfvc_c"), names(b))
  k <- as.numeric(names(b) == "vtpfvc_c") + as.numeric(names(b) %in% tn) * (hh / 24)
  keep <- !is.na(b)
  est <- sum(b[keep] * k[keep])
  se  <- sqrt(as.numeric(t(k[keep]) %*% V[keep, keep, drop = FALSE] %*% k[keep]))
  c(estimate = est, se = se)
}

demo_rhs <- "ns(age10, 4) + sex_category + race_category"
base_rhs <- "vent_day + vtpfvc_c + vtpfvc_c:vent_day + vtpbw_pt_mean + l_vtpbw_within + np_sofa + l_pressor"

decompose <- function(index) {
  ratio_col   <- if (index == "osi") "sf" else "pf"
  ratio_label <- if (index == "osi") "SF ratio" else "P/F ratio"
  dd <- d_all %>% filter(is.finite(.data[[index]]), is.finite(map_aw), is.finite(.data[[ratio_col]]),
                         is.finite(l_vtpbw_within), .data[[index]] > 0, map_aw > 0, .data[[ratio_col]] > 0)
  if (n_distinct(dd$hospitalization_id) < 5 * MIN_CELL) {
    message("  ", index, ": only ", n_distinct(dd$hospitalization_id),
            " patients with all three ingredients; decomposition skipped")
    return(tibble())
  }
  message("  ", index, ": ", nrow(dd), " patient-periods, ", n_distinct(dd$hospitalization_id), " patients")
  map_dfr(c("adjusted", "unadjusted"), function(adj) {
    rhs <- paste(base_rhs, if (adj == "adjusted") paste("+", demo_rhs) else "")
    # identical rows and identical design for all three outcomes: that is what
    # makes the decomposition an identity rather than an approximation
    fits <- map(c(index = index, numerator = "map_aw", ratio = ratio_col), function(y) {
      f <- lm(as.formula(paste("log(", y, ") ~", rhs)), data = dd)
      list(fit = f, V = cluster_vcov(f, dd$hospitalization_id))
    })
    map_dfr(HOURS, function(hh) {
      rows <- imap_dfr(fits, function(o, nm) {
        cc <- contrast_at(o$fit, o$V, hh)
        tibble(part = nm, estimate = cc[["estimate"]], se = cc[["se"]])
      })
      idty <- rows$estimate[rows$part == "numerator"] - rows$estimate[rows$part == "ratio"]
      rows %>% mutate(index = .env$index, ratio_name = ratio_label,
                      adjustment = adj, horizon_h = hh,
                      lo = estimate - 1.96 * se, hi = estimate + 1.96 * se,
                      p = 2 * pnorm(-abs(estimate / se)),
                      identity_check = idty - rows$estimate[rows$part == "index"],
                      n_obs = nrow(dd), n_patients = n_distinct(dd$hospitalization_id))
    })
  })
}
decomp <- bind_rows(decompose("osi"), decompose("oi")) %>%
  mutate(unit = "log units per point of VT/PFVC at a given VT/PBW", site = site_name)
write_csv(decomp, file.path(final_dir, paste0("oi_decomposition_", site_name, ".csv")))

# =============================================================================
# 3. do the two indices agree where both exist?
# =============================================================================
both <- d_all %>% filter(is.finite(oi), is.finite(osi), oi > 0, osi > 0)
agreement <- if (n_distinct(both$hospitalization_id) < MIN_CELL) tibble() else
  both %>%
  mutate(bin = cut(osi, breaks = unique(quantile(osi, seq(0, 1, 0.1), na.rm = TRUE)), include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(n_obs = n(), n_patients = n_distinct(hospitalization_id),
            osi_median = median(osi), oi_median = median(oi),
            oi_q25 = quantile(oi, 0.25), oi_q75 = quantile(oi, 0.75), .groups = "drop") %>%
  filter(n_patients >= MIN_CELL) %>%
  mutate(spearman_all = suppressWarnings(cor(both$oi, both$osi, method = "spearman")),
         n_patients_both = n_distinct(both$hospitalization_id), site = site_name)
write_csv(agreement, file.path(final_dir, paste0("oi_agreement_", site_name, ".csv")))

# =============================================================================
# 4. the figures
# =============================================================================
p1 <- coverage %>%
  select(day, `any patient-day` = n_patients, `mean airway pressure` = n_map_aw,
         `SF ratio` = n_sf, `P/F ratio` = n_pf) %>%
  pivot_longer(-day) %>% filter(!is.na(value)) %>%
  ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1) + geom_point() +
  scale_colour_manual(values = okabe[c(8, 5, 2, 4)], name = NULL) +
  labs(title = "Coverage of each ingredient", x = "ventilator day", y = "patients with a value")

p2 <- coverage %>%
  mutate(OSI = n_osi / n_patients, OI = n_oi / n_patients) %>%
  select(day, OSI, OI) %>% pivot_longer(-day) %>% filter(!is.na(value)) %>%
  ggplot(aes(day, value, fill = name)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = okabe[c(3, 1)], name = NULL) +
  labs(title = "Which index is usable here", subtitle = "share of the day's patients with the index",
       x = "ventilator day", y = "coverage")

p3 <- if (!nrow(decomp)) plot_spacer() else {
  lab <- c(index = "the index", numerator = "mean airway pressure (numerator)", ratio = "the ratio (denominator)")
  decomp %>% filter(horizon_h == max(HOURS)) %>%
    mutate(part_lab = factor(lab[part], rev(lab)),
           panel = paste0(toupper(index), " = 100 x mean airway pressure / ", ratio_name)) %>%
    ggplot(aes(estimate, part_lab, colour = adjustment)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5)) +
    facet_wrap(~ panel, ncol = 1, scales = "free_x") +
    scale_colour_manual(values = okabe[c(1, 2)], name = NULL) +
    labs(title = paste0("Where an index effect comes from, at ", max(HOURS), " hours"),
         subtitle = "the index contrast is the numerator contrast minus the ratio contrast, exactly; a numerator-only effect is arithmetic",
         x = "log units per point of VT/PFVC (95% CI)", y = NULL)
}

p4 <- if (!nrow(decomp)) plot_spacer() else
  decomp %>% filter(part != "index") %>%
  mutate(part_lab = if_else(part == "numerator", "mean airway pressure", "ratio"),
         panel = toupper(index)) %>%
  ggplot(aes(horizon_h, estimate, colour = part_lab, fill = part_lab)) +
  geom_hline(yintercept = 0, colour = "grey55") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + geom_point() +
  facet_grid(panel ~ adjustment, scales = "free_y") +
  scale_colour_manual(values = okabe[c(5, 2)], name = NULL) +
  scale_fill_manual(values = okabe[c(5, 2)], name = NULL) +
  labs(title = "The two components over the horizons",
       x = "hours from the index", y = "log units per point of VT/PFVC")

p5 <- if (!nrow(agreement)) plot_spacer() else
  ggplot(agreement, aes(osi_median, oi_median)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_linerange(aes(ymin = oi_q25, ymax = oi_q75), colour = okabe[4]) +
  geom_point(size = 2.4, colour = okabe[4]) +
  labs(title = "Do the two indices agree where both exist?",
       subtitle = paste0("decile bins of OSI, median OI with IQR; Spearman ",
                         signif(agreement$spearman_all[1], 3), " over ",
                         agreement$n_patients_both[1], " patients. Dashed line is equality."),
       x = "oxygen saturation index (bin median)", y = "oxygenation index (median, IQR)")

pdf_path <- file.path(final_dir, paste0("oi_diagnostics_", site_name, ".pdf"))
pdf(pdf_path, width = 12, height = 8, onefile = TRUE)
print((p1 + p2) +
        plot_annotation(title = paste0(site_name, ": is the oxygenation index available here?"),
                        subtitle = paste0("groups under ", MIN_CELL, " patients suppressed")) &
        theme(legend.position = "top"))
print(p3 + plot_annotation(title = paste0(site_name, ": the lung or the ventilator?"),
                           subtitle = "ordinary least squares, patient-clustered intervals, identical rows and design across the three outcomes") &
        theme(legend.position = "top"))
print(p4 + plot_annotation(title = paste0(site_name, ": the components over time")) &
        theme(legend.position = "top"))
print(p5 + plot_annotation(title = paste0(site_name, ": agreement between the indices")))
invisible(dev.off())

# =============================================================================
message("\nCoverage by day:")
print(as.data.frame(coverage %>% select(day, n_patients, n_map_aw, n_sf, n_pf, n_osi, n_oi)), row.names = FALSE)
if (nrow(decomp)) {
  message("\nDecomposition (log units per point of VT/PFVC; index = numerator - ratio):")
  print(as.data.frame(decomp %>% filter(horizon_h == max(HOURS)) %>%
                        select(index, adjustment, part, estimate, lo, hi, p, identity_check, n_patients) %>%
                        mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)
  chk <- max(abs(decomp$identity_check), na.rm = TRUE)
  message(sprintf("Largest departure from the identity: %.2e (should be numerically zero)", chk))
}
message("\n13_oi_diagnostics complete -> ", final_dir)
