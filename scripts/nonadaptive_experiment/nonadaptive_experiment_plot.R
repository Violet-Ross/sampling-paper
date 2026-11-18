# Redraws the non-adaptive experiment figure from saved results, without
# refitting anything. The plotting code is copied from
# nonadaptive_experiment.R, which cannot be sourced because it fits at top
# level; only the labels differ.
#
# Two figures are written. They share the right-hand draws panels, which are
# always the least and most BIASED cells; only the left-hand heatmaps differ.
#
# Inputs (written by nonadaptive_experiment.R):
#   throughput/nonadaptive_experiment_scores.csv       bias per fit and parameter
#   throughput/nonadaptive_experiment_diagnostics.csv  Rhat and n_eff per fit
#   throughput/nonadaptive_experiment_examples.rds     the two fits drawn below
#
# Outputs:
#   fig/nonadaptive_experiment_plot.png              heatmaps of mean signed bias
#   fig/nonadaptive_experiment_plot_supplementary.png  heatmaps of mean Rhat and
#                                                    mean n_eff, side by side
#
# Run from the repository root:
#   Rscript scripts/nonadaptive_experiment/nonadaptive_experiment_plot.R

library(ggplot2)
library(gridExtra)

if (!interactive()) pdf(NULL)

SCORE_FILE    <- "throughput/nonadaptive_experiment_scores.csv"
DIAG_FILE     <- "throughput/nonadaptive_experiment_diagnostics.csv"
EXAMPLES_FILE <- "throughput/nonadaptive_experiment_examples.rds"
FIG_FILE      <- "fig/nonadaptive_experiment_plot.png"
DIAG_FIG_FILE <- "fig/nonadaptive_experiment_plot_supplementary.png"

# ---- Labels -----------------------------------------------------------------

# Plotmath expressions. The script-l is the Unicode character U+2113, since
# plotmath has no symbol for it.
LAB_L_0      <- quote("\u2113"[0])
LAB_PRIOR_SD <- quote(sigma[prior])

# Heatmap panel headers, keyed by the parameter names in the scores file.
# The order here is the top-to-bottom order of the heatmaps.
PAR_LABELS <- list(mu_V      = quote(mu[V]),
                   mu_delt_p = quote(mu[Delta[p]]),
                   mu_delt_c = quote(mu[Delta[c]]))
PLOT_PARS  <- names(PAR_LABELS)

# ---- Curve ------------------------------------------------------------------

witch <- function(t, V, delt_p, delt_c) {
  ifelse(t > -delt_p & t <= 0, V * (1 + t / delt_p),
  ifelse(t >  0      & t <  delt_c, V * (1 - t / delt_c), 0))
}

# One row per (person, draw, x): that draw's witch hat on a fine grid over the
# person's record, in the record-local days the model was given.
draw_curves <- function(ex) {
  do.call(rbind, lapply(seq_along(ex$persons), function(k) {
    p    <- ex$persons[k]
    span <- max(ex$points$x[ex$points$person == p])
    xx   <- seq(0, span, by = 0.1)
    do.call(rbind, lapply(seq_len(nrow(ex$psi)), function(d) {
      data.frame(person = p, draw = d, x = xx,
                 y = witch(xx - ex$psi[d, k], ex$V[d, k],
                           ex$delt_p[d, k], ex$delt_c[d, k]))
    }))
  }))
}

# ---- Figure -----------------------------------------------------------------

SURFACE   <- "#fcfcfb"
MISSING   <- "#e8e8e6"   # cells where every fit failed
DRAW_INK  <- "#2a78d6"
DATA_INK  <- "#0b0b0b"

# Bias is signed, so the fill diverges about 0: hue carries the sign and an
# unbiased cell reads as blank.
BIAS_RAMP <- c("#8c510a", "#c88a2e", "#e8d3a0", "#f7f7f7",
               "#a8c9ee", "#3d7fc1", "#0d366b")

# Rhat and effective sample size are unsigned, so each gets one hue running
# light to dark: purple for Rhat, teal for effective sample size. Neither
# reuses the brown or blue of the bias ramp.
RHAT_RAMP <- c("#f7f7f7", "#dcd0ea", "#b199d1", "#7a52a8", "#3f1d6e")
NEFF_RAMP <- c("#f7f7f7", "#cfe8e2", "#8fcabd", "#3a9a88", "#0f5c4e")

# Label ink flips to white on the dark ends of the ramp.
ink_for <- function(pos) ifelse(!is.na(pos) & pos > 0.6, "#ffffff", "#0b0b0b")

# What differs between the three heatmap columns. Each metric gives:
#   title   column heading, as a function of the number of iterations
#   fmt     sprintf format of the value printed on each tile
#   scale   the fill scale for one panel, given that panel's finite values
#   dark    how far along the ramp towards its dark end each value sits, in
#           [0, 1]; used only to pick the label ink
# In every case the scale is taken from one parameter alone, so panels cannot
# be compared by colour; every tile carries its value for that reason.
METRICS <- list(

  # Symmetric about 0.
  bias = list(
    title = function(n) sprintf("Mean signed bias over %d iterations", n),
    fmt   = "%+.2f",
    scale = function(v) {
      m <- max(abs(v), 0); if (m == 0) m <- 1
      scale_fill_gradientn(colours = BIAS_RAMP, limits = c(-m, m),
                           na.value = MISSING)
    },
    dark  = function(v) {
      m <- max(abs(v), 0, na.rm = TRUE); if (m == 0) m <- 1
      abs(v) / m
    }),

  # Log colour scale from 1 upward: a few fits with a stuck chain have Rhat in
  # the tens or hundreds and pull their cell's mean far above the rest, which
  # on a linear scale would leave every other tile blank. Dark is bad.
  rhat = list(
    title = function(n) bquote("Mean " * hat(R) * " over " * .(n) * " iterations"),
    fmt   = "%.2f",
    scale = function(v) {
      lim <- c(min(1, v), max(v, 1.01))
      scale_fill_gradientn(colours = RHAT_RAMP, limits = lim,
                           transform = "log10",
                           breaks = scales::breaks_log(4),
                           na.value = MISSING)
    },
    dark  = function(v) {
      lim <- log10(c(min(1, v, na.rm = TRUE), max(v, 1.01, na.rm = TRUE)))
      (log10(v) - lim[1]) / diff(lim)
    }),

  # Log colour scale, since n_eff runs from single digits to thousands. The
  # ramp is reversed so that dark is bad here too: dark means FEW effective
  # draws.
  neff = list(
    title = function(n) sprintf("Mean effective sample size over %d iterations", n),
    fmt   = "%.0f",
    scale = function(v) {
      scale_fill_gradientn(colours = rev(NEFF_RAMP), limits = range(v),
                           transform = "log10",
                           breaks = scales::breaks_log(4),
                           na.value = MISSING)
    },
    dark  = function(v) {
      lim <- log10(range(v, na.rm = TRUE))
      if (diff(lim) == 0) return(rep(0, length(v)))
      1 - (log10(v) - lim[1]) / diff(lim)
    })
)

# One heatmap: the per-cell mean of one metric for one parameter. `cells` has
# columns l_0, prior_sd, parameter and value.
#
# `starred` holds the l_0 and prior_sd of the two cells drawn in the right-hand
# panels: the least and the most biased, ranked by mean absolute relative bias
# over the three population means. The same two cells are starred on every
# heatmap of every figure, in the tile's top-right corner. The star's colour
# is fixed by which cell it marks (column star_ink), not by the tile under it:
# black for the least biased cell, white for the most biased.
heat_panel <- function(cells, par, starred, metric) {
  x  <- cells[cells$parameter == par, ]
  ok <- is.finite(x$value)
  x$label <- ifelse(ok, sprintf(metric$fmt, x$value), "NA")
  x$ink   <- ink_for(metric$dark(x$value))
  stars   <- merge(x, starred)

  ggplot(x, aes(factor(l_0), factor(prior_sd), fill = value)) +
    geom_tile(colour = SURFACE, linewidth = 1.5) +
    geom_text(aes(label = label, colour = ink), size = 2.6) +
    geom_text(data = stars, aes(colour = star_ink), label = "\u2605", size = 3.2,
              position = position_nudge(x = 0.3, y = 0.28)) +
    scale_colour_identity() +
    metric$scale(x$value[ok]) +
    coord_equal() +
    labs(x = LAB_L_0, y = LAB_PRIOR_SD, fill = NULL,
         title = PAR_LABELS[[par]]) +
    theme_minimal(base_size = 10) +
    theme(panel.grid       = element_blank(),
          plot.title       = element_text(size = 11, hjust = 0.5),
          legend.key.width = grid::unit(0.25, "cm"))
}

plot_heat <- function(cells, n_reps, starred, metric) {
  arrangeGrob(
    grobs = lapply(PLOT_PARS,
                   function(p) heat_panel(cells, p, starred, metric)),
    ncol  = 1,
    top   = grid::textGrob(metric$title(n_reps),
                           gp = grid::gpar(fontsize = 11)))
}

plot_draws <- function(ex, heading) {
  curves <- draw_curves(ex)
  points <- ex$points
  lab    <- function(p) sprintf("person %d", p)
  curves$panel <- factor(lab(curves$person), levels = lab(ex$persons))
  points$panel <- factor(lab(points$person), levels = lab(ex$persons))

  ggplot() +
    geom_line(data = curves, aes(x, y, group = draw),
              colour = DRAW_INK, alpha = 0.15, linewidth = 0.4) +
    geom_point(data = points, aes(x, y),
               colour = DATA_INK, size = 1.4) +
    facet_wrap(~ panel, nrow = 1) +
    labs(x = "days since first assayed day of record", y = "viral load",
         title = heading) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.title       = element_text(size = 11),
          strip.text       = element_text(size = 8))
}

# ---- Run --------------------------------------------------------------------

scores   <- read.csv(SCORE_FILE)
diags    <- read.csv(DIAG_FILE)
examples <- readRDS(EXAMPLES_FILE)
n_reps   <- length(unique(scores$rep))

# Mean per cell and parameter of one column of `df`, returned as `value`.
# Merged onto the full grid so a cell whose fits all failed shows as missing,
# not as a gap.
grid <- expand.grid(prior_sd = sort(unique(scores$prior_sd)),
                    l_0      = sort(unique(scores$l_0)))
full <- merge(grid, data.frame(parameter = PLOT_PARS))

cell_mean <- function(df, column) {
  df$value <- df[[column]]
  merge(full,
        aggregate(value ~ l_0 + prior_sd + parameter, data = df, FUN = mean),
        all.x = TRUE)
}

# The cells the two right-hand panels were drawn from.
starred <- data.frame(l_0      = c(examples$low$l_0,      examples$high$l_0),
                      prior_sd = c(examples$low$prior_sd, examples$high$prior_sd),
                      star_ink = c("#0b0b0b", "#ffffff"))

# The draws panels are identical in all three figures, so they are built once.
draws_low  <- plot_draws(examples$low,
                         "Least biased iteration of the least biased cell")
draws_high <- plot_draws(examples$high,
                         "Most biased iteration of the most biased cell")

# `heat_cols` is a list of heatmap columns, each from plot_heat(). They sit
# side by side on the left, each as wide as the single column of the bias
# figure, and the page is widened to fit so the draws panels keep their size.
write_figure <- function(heat_cols, file) {
  n_heat <- length(heat_cols)
  fig <- arrangeGrob(
    arrangeGrob(grobs = heat_cols, nrow = 1),
    draws_low, draws_high,
    layout_matrix = cbind(c(1, 1), c(2, 3)),
    widths = c(n_heat, 2.6))

  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  ggsave(file, fig, width = 15 * (n_heat + 2.6) / 3.6, height = 8.5,
         dpi = 600, bg = SURFACE, device = ragg::agg_png)
  cat(sprintf("wrote %s\n", file))
}

write_figure(
  list(plot_heat(cell_mean(scores, "bias"), n_reps, starred, METRICS$bias)),
  FIG_FILE)

write_figure(
  list(plot_heat(cell_mean(diags, "Rhat"),  n_reps, starred, METRICS$rhat),
       plot_heat(cell_mean(diags, "n_eff"), n_reps, starred, METRICS$neff)),
  DIAG_FIG_FILE)
