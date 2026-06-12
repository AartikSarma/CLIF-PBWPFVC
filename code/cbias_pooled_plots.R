# =============================================================================
# DEFERRED: Pooled conditional-bias plots (cross-cohort)
# PBW vs PFVC Replication Using CLIF Data
# =============================================================================
# Standalone, NOT part of the 00 pipeline runner. Reconstructs the styled
# conditional-bias plots across cohorts from each site's MASKED
# cbias_export_<site>.csv (code/cbias_federated_export.R).
#
# *** DEFERRED ***
# The per-site exports contain (stratum x percentile) cells with n < 10, so they
# may only be collected after passing through the consortium's deterministic
# additive-masking pipeline (held until all sites confirm). Run this script once
# the masked exports have been gathered into the local output/ tree.
#
# From the per-(plot, dependent, independent, stratum, percentile) counts and sums,
# pooled N-weighted, this reproduces the per-site plot style:
#   - faint percentile x's     (mean per percentile),
#   - a weighted loess line     (geom_smooth on the percentile means, weight = N),
#   - decile circles + 95% CI   (ten percentiles aggregated per point; Wilson for
#                                binary outcomes, normal for continuous), and
#   - a per-group density band  (federated stand-in for the per-subject beeswarm).
# Percentiles are within-site, so a given percentile blends slightly different
# absolute ranges across cohorts (accepted imprecision).
#
# Inputs : output/*_output/final/cbias_export_<site>.csv  (masked)
# Outputs: output/cross_cohort/conditional_bias/conditional_bias_pooled_all.pdf
#          output/cross_cohort/cbias_pooled_all_cohorts.csv
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(patchwork)

cross_dir <- here("output", "cross_cohort")
cbias_dir <- file.path(cross_dir, "conditional_bias")
dir.create(cbias_dir, recursive = TRUE, showWarnings = FALSE)

okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#000000")

cbias_files <- Sys.glob(here("output", "*_output", "final", "cbias_export_*.csv"))
if (length(cbias_files) == 0) {
  stop("No per-site cbias_export_*.csv found. Run code/cbias_federated_export.R ",
       "for each cohort (and apply masking) first.")
}

cbias_all <- map_dfr(cbias_files, read_csv, show_col_types = FALSE)
n_cbias_sites <- n_distinct(cbias_all$site)
z <- qnorm(0.975)

# Pool to per-(config, stratum, percentile): N-weighted sums.
perc_pooled <- cbias_all %>%
  group_by(plot_id, dependent, dep_label, independent, grouping, stratum,
           percentile, is_binary) %>%
  summarise(n = sum(n), s = sum(sum_dep), ss = sum(sumsq_dep), .groups = "drop") %>%
  mutate(mean = s / n, x_pos = percentile - 0.5)

# Deciles: aggregate ten percentiles per point, then mean + 95% CI.
dec_pooled <- perc_pooled %>%
  mutate(decile = pmin(ceiling(percentile / 10), 10)) %>%
  group_by(plot_id, dependent, dep_label, independent, grouping, stratum, decile,
           is_binary) %>%
  summarise(n = sum(n), s = sum(s), ss = sum(ss), .groups = "drop") %>%
  mutate(
    mean    = s / n,
    x_pos   = decile * 10 - 5,
    cont_se = sqrt(pmax(ss - s^2 / n, 0) / pmax(n - 1, 1) / n),
    wil_c   = (mean + z^2 / (2 * n)) / (1 + z^2 / n),
    wil_s   = z * sqrt(pmax(mean * (1 - mean) + z^2 / (4 * n), 0) / n) / (1 + z^2 / n),
    ci_lo   = if_else(is_binary, wil_c - wil_s, mean - z * cont_se),
    ci_hi   = if_else(is_binary, wil_c + wil_s, mean + z * cont_se)
  )

write_csv(dec_pooled %>% select(plot_id, dependent, dep_label, independent,
                                grouping, stratum, decile, n, mean, ci_lo, ci_hi),
          file.path(cross_dir, "cbias_pooled_all_cohorts.csv"))

grouping_labels    <- c(race_category = "Race", sex_category = "Sex",
                        age_group = "Age group")
independent_labels <- c(vtpfvc = "VT/PFVC", vtpbw = "VT/PBW", pbwpfvc = "PBW/PFVC",
                        pbw = "PBW", pfvc = "PFVC")
strata_levels <- sort(unique(perc_pooled$stratum))
cbias_pal <- setNames(grDevices::colorRampPalette(okabe_ito)(length(strata_levels)),
                      strata_levels)

# Mirrored horizontal density band per stratum from its per-percentile counts —
# the federated analog of the per-subject beeswarm (no individuals leave a site).
density_band <- function(per_df, strata_order) {
  per_df %>%
    group_by(independent, stratum) %>%
    group_modify(function(d, key) {
      if (nrow(d) < 2) return(tibble(x = numeric(0), off = numeric(0)))
      grid <- seq(0, 100, length.out = 120)
      dens <- approx(d$x_pos, d$n, xout = grid, rule = 2)$y
      dens[dens < 0] <- 0
      if (max(dens) > 0) dens <- dens / max(dens)
      half <- 0.43 * dens
      tibble(x = c(grid, rev(grid)), off = c(half, -rev(half)))
    }) %>%
    ungroup() %>%
    mutate(y = match(stratum, strata_order) + off,
           poly = paste(independent, stratum))
}

make_cbias_page <- function(pid, gv) {
  pp <- perc_pooled %>% filter(plot_id == pid, grouping == gv)
  dd <- dec_pooled  %>% filter(plot_id == pid, grouping == gv)
  if (nrow(pp) == 0) return(NULL)
  dep_label <- pp$dep_label[1]
  gv_label  <- recode(gv, !!!grouping_labels)
  relab <- function(d) d %>% mutate(independent = recode(independent, !!!independent_labels))
  pp <- relab(pp); dd <- relab(dd)
  strata_here <- sort(unique(pp$stratum))

  # y-range focused on the decile points/CIs; faint percentile x's squished in.
  yr  <- range(c(dd$ci_lo, dd$ci_hi), na.rm = TRUE)
  pad <- diff(yr) * 0.15
  y_range <- c(yr[1] - pad, yr[2] + pad)

  main <- ggplot() +
    geom_point(data = pp, aes(x_pos, scales::squish(mean, y_range), color = stratum),
               shape = 4, size = 1, alpha = 0.4) +
    geom_smooth(data = pp, aes(x_pos, mean, color = stratum, weight = n),
                method = "loess", formula = y ~ x, se = FALSE, span = 0.75,
                linewidth = 0.7) +
    geom_point(data = dd, aes(x_pos, mean, color = stratum), size = 2,
               position = position_dodge(width = 4)) +
    geom_errorbar(data = dd, aes(x_pos, ymin = ci_lo, ymax = ci_hi, color = stratum),
                  width = 2.5, position = position_dodge(width = 4)) +
    facet_wrap(~ independent, nrow = 1) +
    scale_color_manual(values = cbias_pal, name = gv_label) +
    coord_cartesian(xlim = c(0, 100), ylim = y_range) +
    labs(x = NULL, y = dep_label) +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(face = "bold"),
          panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.4))

  dens <- density_band(pp, strata_here)
  bottom <- ggplot(dens, aes(x, y, group = poly, fill = stratum)) +
    geom_polygon(alpha = 0.75, color = NA) +
    facet_wrap(~ independent, nrow = 1) +
    scale_fill_manual(values = cbias_pal, guide = "none") +
    scale_y_continuous(breaks = seq_along(strata_here), labels = strata_here) +
    coord_cartesian(xlim = c(0, 100)) +
    labs(x = "Percentile of independent variable", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.4))

  (main / bottom) +
    plot_layout(heights = c(3, 1)) +
    plot_annotation(
      title = paste0("Pooled conditional bias: ", dep_label, " by ", gv_label),
      subtitle = paste0("Pooled across ", n_cbias_sites,
        " cohort(s), N-weighted. x = percentile mean, line = weighted loess, ",
        "circles = decile mean +/- 95% CI; lower band = per-group density ",
        "(federated beeswarm). Percentiles are within-site."),
      theme = theme(plot.title = element_text(face = "bold")))
}

cbias_gvs <- intersect(c("race_category", "sex_category", "age_group"),
                       unique(perc_pooled$grouping))
pdf(file.path(cbias_dir, "conditional_bias_pooled_all.pdf"), width = 11, height = 8)
n_cbias_pages <- 0
for (pid in unique(perc_pooled$plot_id)) {
  for (gv in cbias_gvs) {
    pg <- make_cbias_page(pid, gv)
    if (is.null(pg)) next
    print(pg); n_cbias_pages <- n_cbias_pages + 1
  }
}
invisible(dev.off())
message("Pooled conditional-bias plots: ", n_cbias_pages,
        " pages (config x stratifier) across ", n_cbias_sites, " cohort(s): ",
        file.path(cbias_dir, "conditional_bias_pooled_all.pdf"))
