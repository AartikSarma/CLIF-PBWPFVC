# =============================================================================
# CONSORT-style inclusion flow diagram (ggplot2 only)
# =============================================================================
# Renders the cohort attrition funnel from the structured attrition tibble
# produced by utils/attrition_log.R, using only ggplot2 primitives so no
# DiagrammeR / consort package dependency is introduced.
#
# The main vertical chain shows n_remaining at each step; right-hand-side branch
# boxes show the number excluded and the reason at each transition.

library(tidyverse)

render_consort <- function(attrition_tbl, title = NULL) {
  at <- attrition_tbl %>% arrange(step_order)
  n <- nrow(at)
  if (n == 0) stop("render_consort: empty attrition table")

  # Evenly space main boxes down a left-of-center column.
  ys <- seq(1, 0, length.out = n)
  main_x <- 0.30
  excl_x <- 0.80
  # Half the vertical gap between boxes, used to start/stop connector arrows.
  half_gap <- if (n > 1) (ys[1] - ys[2]) / 2 * 0.55 else 0.1

  main <- tibble(
    x = main_x,
    y = ys,
    label = paste0(at$step_label, "\nn = ", format(at$n_remaining, big.mark = ","))
  )

  # Vertical arrows between consecutive main boxes.
  v <- tibble(
    x = main_x, xend = main_x,
    y = head(ys, -1) - half_gap,
    yend = tail(ys, -1) + half_gap
  )

  excl_idx <- which(!is.na(at$n_excluded) & at$n_excluded > 0)
  has_excl <- length(excl_idx) > 0

  if (has_excl) {
    excl <- tibble(
      x = excl_x,
      y = (ys[excl_idx] + ys[excl_idx - 1]) / 2,
      label = paste0(
        coalesce(at$exclusion_reason[excl_idx], "Excluded"),
        "\nn = ", format(at$n_excluded[excl_idx], big.mark = ",")
      )
    )
    # Horizontal arrows branching off the main chain to each exclusion box.
    h <- tibble(
      x = main_x + 0.02, xend = excl_x - 0.16,
      y = excl$y, yend = excl$y
    )
  }

  p <- ggplot() +
    geom_segment(data = v, aes(x = x, y = y, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
                 linewidth = 0.4, colour = "grey30")

  if (has_excl) {
    p <- p +
      geom_segment(data = h, aes(x = x, y = y, xend = xend, yend = yend),
                   arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
                   linewidth = 0.4, colour = "#D55E00") +
      geom_label(data = excl, aes(x = x, y = y, label = label),
                 fill = "#FBE7DD", colour = "#D55E00", linewidth = 0.3,
                 size = 2.8, lineheight = 0.95,
                 label.padding = unit(0.35, "lines"))
  }

  p <- p +
    geom_label(data = main, aes(x = x, y = y, label = label),
               fill = "grey95", colour = "black", linewidth = 0.3,
               size = 3.1, lineheight = 0.95,
               label.padding = unit(0.45, "lines")) +
    scale_x_continuous(limits = c(0.05, 1.05)) +
    scale_y_continuous(limits = c(-0.08, 1.08)) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))

  if (!is.null(title)) p <- p + ggtitle(title)
  p
}
