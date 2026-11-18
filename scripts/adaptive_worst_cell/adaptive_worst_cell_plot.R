# Redraws the adaptive_worst_cell figure from the saved scores, without
# refitting anything.
#
# Three panels, one per population mean (mu_V, mu_delt_p, mu_delt_c): signed
# bias (posterior mean - generating value) against the number of points added
# per infection, 0 being the starting fit, with one line per iteration.
#
# Reads:   throughput/adaptive_worst_cell_scores.csv
# Writes:  fig/adaptive_worst_cell.png
#
# Run from the repository root:
#   Rscript scripts/adaptive_worst_cell/adaptive_worst_cell_plot.R

library(ggplot2)
library(gridExtra)

if (!interactive()) pdf(NULL)

SCORE_FILE <- "throughput/adaptive_worst_cell_scores.csv"
FIG_FILE   <- "fig/adaptive_worst_cell.png"

PLOT_PARS  <- c("mu_V", "mu_delt_p", "mu_delt_c")

scores <- read.csv(SCORE_FILE)

# ---- Figure -----------------------------------------------------------------

SURFACE  <- "#fcfcfb"
LINE_INK <- "#2a78d6"
MUTED    <- "#6b6b68"

# Plotmath panel headings, as nonadaptive_experiment_plot.R.
PAR_LABELS <- list(mu_V      = quote(mu[V]),
                   mu_delt_p = quote(mu[Delta[p]]),
                   mu_delt_c = quote(mu[Delta[c]]))

# One panel: the signed bias of one parameter after each addition, one line
# per iteration. The lines share one colour because the iterations are
# exchangeable replicates. The y axis is the panel's own, since the three
# biases differ several-fold in size; the line at 0 marks an unbiased fit.
bias_panel <- function(scores, par) {
  x <- scores[scores$parameter == par, ]
  ggplot(x, aes(round, bias, group = rep)) +
    geom_hline(yintercept = 0, colour = MUTED, linewidth = 0.4) +
    geom_line(colour = LINE_INK, linewidth = 0.6, alpha = 0.6) +
    geom_point(colour = LINE_INK, size = 1.2, alpha = 0.6) +
    scale_x_continuous(breaks = sort(unique(scores$round))) +
    labs(x = "points added per infection", y = "signed bias",
         title = PAR_LABELS[[par]]) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.title       = element_text(size = 11, hjust = 0.5))
}

fig <- arrangeGrob(
  grobs = lapply(PLOT_PARS, function(p) bias_panel(scores, p)),
  nrow  = 1)

ggsave(FIG_FILE, fig, width = 11, height = 3.6, dpi = 600, bg = SURFACE)

cat(sprintf("wrote %s\n", FIG_FILE))
