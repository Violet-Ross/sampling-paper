# Non-adaptive sampling experiment: a sweep over the interval-sampling spacing
# l_0 and the prior standard deviation prior_sd, with N_REPS iterations per
# cell.
#
# Each iteration draws ONE synthetic sample matrix, and every l_0 x prior_sd
# cell is fitted to it:
#   1. infection discovery: pool each row, then assay the positive rows on days
#      0, l_0, 2*l_0, ... and keep the assays within W_DAYS of the first and
#      last positive;
#   2. fit the BHM in src/bhm.stan, with prior_sd as
#      the prior standard deviation and the generating values as prior centres;
#   3. record bias, Rhat and n_eff.
#
# Outputs:
#   throughput/nonadaptive_experiment_scores.csv       bias per fit and parameter
#   throughput/nonadaptive_experiment_diagnostics.csv  Rhat, n_eff, divergences
#   throughput/nonadaptive_experiment_examples.rds     the two fits drawn below
#   fig/nonadaptive_experiment.png
#
# The figure has three parts. Left: one heatmap per population mean of the
# signed bias averaged over iterations. Top right: HMC draws and assayed points
# for N_SHOW infections from the least biased iteration of the least biased
# cell. Bottom right: the same for the most biased iteration of the most
# biased cell. "Biased" for choosing cells and iterations means the mean
# absolute relative bias over the three population means.
#
# Run from the repository root:
#   Rscript scripts/nonadaptive_experiment/nonadaptive_experiment.R
# or on the cluster:
#   sbatch scripts/nonadaptive_experiment/nonadaptive_experiment.sh

library(rstan)
library(ggplot2)
library(gridExtra)

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# Under Rscript, assembling the figure would otherwise open a default device
# and leave an Rplots.pdf in the working directory.
if (!interactive()) pdf(NULL)

# ---- Settings ---------------------------------------------------------------

L_0_GRID      <- c(1, 2, 3, 4, 5)   # interval-sampling spacings, in days
PRIOR_SD_GRID <- c(1, 2, 3)         # prior standard deviations
N_REPS        <- 10                 # iterations per l_0 x prior_sd cell

W_DAYS      <- 10      # assays kept this many days either side of the positives

# Generating values. These are also the prior centres passed to Stan.
TRUE_V      <- 4       # peak viral load
TRUE_DELT_P <- 3.5       # proliferation time, days
TRUE_DELT_C <- 10       # clearance time, days
SIGMA_IND   <- 0.1     # person-to-person spread
SIGMA_EPS   <- 0.1     # observation noise

M           <- 500     # days (matrix columns)
N_PEOPLE    <- 500     # people (matrix rows)
PREV        <- 0.3     # fraction of people infected

POS_THRESH  <- 0       # an assay is positive above this value
PSI_SHAPE   <- 2       # Beta(a, a) prior on the peak time within its bracket

N_CHAINS    <- 4
N_ITER      <- 2000    # per chain; Stan warms up for half
SEED        <- 1

N_SHOW      <- 5       # example infections per draws panel
N_DRAWS     <- 100     # HMC draws per example infection, chains pooled

STAN_FILE     <- "src/bhm.stan"
SCORE_FILE    <- "throughput/nonadaptive_experiment_scores.csv"
DIAG_FILE     <- "throughput/nonadaptive_experiment_diagnostics.csv"
EXAMPLES_FILE <- "throughput/nonadaptive_experiment_examples.rds"
FIG_FILE      <- "fig/nonadaptive_experiment.png"

stopifnot(W_DAYS >= max(L_0_GRID), N_REPS >= 1,
          N_DRAWS <= N_CHAINS * (N_ITER %/% 2))
for (f in c(SCORE_FILE, DIAG_FILE, EXAMPLES_FILE, FIG_FILE))
  dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)

# ---- Curve ------------------------------------------------------------------

witch <- function(t, V, delt_p, delt_c) {
  ifelse(t > -delt_p & t <= 0, V * (1 + t / delt_p),
  ifelse(t >  0      & t <  delt_c, V * (1 - t / delt_c), 0))
}

rlnorm_at <- function(median, sigma) median * exp(rnorm(1, 0, sigma / median))

# ---- Generation -------------------------------------------------------------

# Rows are people, columns are days 0..M-1. Uninfected rows and off-cone days
# are exactly 0; noise is added on the cone only and truncated at 0.
generate_matrix <- function() {
  mat      <- matrix(0, nrow = N_PEOPLE, ncol = M)
  is_infx  <- rbinom(N_PEOPLE, 1, PREV)
  days     <- 0:(M - 1)
  infected <- which(is_infx == 1)

  truth <- data.frame(person = infected, V = NA_real_,
                      delt_p = NA_real_, delt_c = NA_real_)

  for (k in seq_along(infected)) {
    i      <- infected[k]
    V      <- rlnorm_at(TRUE_V,      SIGMA_IND)
    delt_p <- rlnorm_at(TRUE_DELT_P, SIGMA_IND)
    delt_c <- rlnorm_at(TRUE_DELT_C, SIGMA_IND)
    peak   <- runif(1, delt_p + 1, M - delt_c - 1)

    mu     <- witch(days - peak, V, delt_p, delt_c)

    inside <- mu > 0
    mat[i, inside] <- pmax(mu[inside] + rnorm(sum(inside), 0, SIGMA_EPS), 0)

    truth[k, c("V", "delt_p", "delt_c")] <- c(V, delt_p, delt_c)
  }

  list(mat = mat, truth = truth, n_infx = length(infected))
}

# ---- Infection discovery ----------------------------------------------------

pool_positive <- function(row) any(row > 0)

interval_days <- function(l_0) seq(0, M - 1, by = l_0)

sample_rows <- function(mat, l_0, w_days = W_DAYS) {
  stopifnot(w_days >= l_0)

  days      <- interval_days(l_0)
  cols      <- days + 1
  records   <- list()
  brackets  <- list()
  n_pool    <- 0
  n_missed  <- 0
  n_clamped <- 0

  for (i in seq_len(nrow(mat))) {
    row <- mat[i, ]
    if (!pool_positive(row)) next
    n_pool <- n_pool + 1

    vals <- as.numeric(row[cols])
    pos  <- which(vals > POS_THRESH)

    # Pool positive but no positive interval assay: the infection is missed.
    if (length(pos) == 0) {
      n_missed <- n_missed + 1
      next
    }

    lo   <- days[min(pos)] - w_days
    hi   <- days[max(pos)] + w_days
    keep <- which(days >= lo & days <= hi)

    records[[length(records) + 1]] <- data.frame(
      person = i,
      x      = days[keep] - days[keep[1]],
      y      = vals[keep]
    )

    # The bracketing negatives: the last assayed negative before the first
    # positive and the first assayed negative after the last positive. The
    # peak lies between them. Both are inside `keep` whenever they exist,
    # since w_days >= l_0. A person positive on the first or last assayed day
    # of the study has no bracketing negative on that side, so that side falls
    # back to the record edge; n_clamped counts those.
    org  <- days[keep[1]]
    lo_i <- min(pos) - 1L
    hi_i <- max(pos) + 1L
    n_clamped <- n_clamped + (lo_i < 1L) + (hi_i > length(days))

    brackets[[length(brackets) + 1]] <- data.frame(
      person = i,
      psi_lo = if (lo_i < 1L) 0 else days[lo_i] - org,
      psi_hi = if (hi_i > length(days)) days[keep[length(keep)]] - org
               else days[hi_i] - org
    )
  }

  list(records = records, brackets = brackets, n_pool = n_pool,
       n_missed = n_missed, n_clamped = n_clamped)
}

build_stan_data <- function(records, brackets, prior_sd) {
  df      <- do.call(rbind, records)
  persons <- unique(df$person)
  df$id   <- match(df$person, persons)

  # Brackets reordered onto `persons`, the order `id` indexes.
  b       <- do.call(rbind, brackets)
  b       <- b[match(persons, b$person), ]
  stopifnot(!anyNA(b$psi_lo), !anyNA(b$psi_hi), all(b$psi_hi > b$psi_lo))

  list(persons = persons, df = df,
       data = list(N = nrow(df), J = length(persons), id = df$id,
                   x = df$x, y = df$y,
                   psi_lo = b$psi_lo, psi_hi = b$psi_hi,
                   psi_shape = PSI_SHAPE,
                   prior_sd = prior_sd,
                   prior_mu_V      = TRUE_V,
                   prior_mu_delt_p = TRUE_DELT_P,
                   prior_mu_delt_c = TRUE_DELT_C))
}

# ---- Fit --------------------------------------------------------------------

# psi, V, delt_p and delt_c are kept because the draws panels need them; only
# the _raw offsets are dropped.
fit_model <- function(model, stan_data, seed) {
  sampling(model, data = stan_data,
           chains = N_CHAINS, iter = N_ITER, seed = seed,
           pars = c("V_raw", "delt_p_raw", "delt_c_raw"),
           include = FALSE)
}

# ---- Scoring ----------------------------------------------------------------

PLOT_PARS <- c("mu_V", "mu_delt_p", "mu_delt_c")

NOMINAL <- list(mu_V        = TRUE_V,
                mu_delt_p   = TRUE_DELT_P,
                mu_delt_c   = TRUE_DELT_C,
                mu_duration = TRUE_DELT_P + TRUE_DELT_C,
                sigma_eps   = SIGMA_EPS)

# One row per parameter. `bias` is the posterior mean minus the generating
# value; `rmse` is taken over the posterior draws.
score <- function(fit) {
  post <- rstan::extract(fit, pars = names(NOMINAL))
  do.call(rbind, lapply(names(NOMINAL), function(p) {
    d <- post[[p]]
    n <- NOMINAL[[p]]
    data.frame(parameter = p,
               truth     = n,
               estimate  = mean(d),
               bias      = mean(d - n),
               rmse      = sqrt(mean((d - n)^2)))
  }))
}

# The single number used to rank fits and cells: the mean over the three
# population means of |bias| / truth.
abs_rel_bias <- function(sc) {
  s <- sc[sc$parameter %in% PLOT_PARS, ]
  mean(abs(s$bias) / s$truth)
}

DIAG_PARS <- c("mu_V", "mu_delt_p", "mu_delt_c",
               "sigma_log_V", "sigma_log_delt_p", "sigma_log_delt_c",
               "sigma_eps", "mu_duration")

diagnostics <- function(fit) {
  s <- summary(fit, pars = DIAG_PARS)$summary
  data.frame(parameter = rownames(s),
             n_eff     = s[, "n_eff"],
             Rhat      = s[, "Rhat"],
             row.names = NULL)
}

n_divergent <- function(fit) {
  sum(sapply(get_sampler_params(fit, inc_warmup = FALSE),
             function(p) sum(p[, "divergent__"])))
}

# ---- Example draws ----------------------------------------------------------

# What a draws panel needs from one fit, small enough to keep for every fit:
# N_DRAWS posterior draws of each shown person's curve parameters, and the
# assayed (x, y) points the model was given for those people.
example_payload <- function(fit, prep, show, draw_idx) {
  post <- rstan::extract(fit, pars = c("psi", "V", "delt_p", "delt_c"))
  j    <- match(show, prep$persons)
  list(persons = show,
       psi     = post$psi[draw_idx, j, drop = FALSE],
       V       = post$V[draw_idx, j, drop = FALSE],
       delt_p  = post$delt_p[draw_idx, j, drop = FALSE],
       delt_c  = post$delt_c[draw_idx, j, drop = FALSE],
       points  = prep$df[prep$df$person %in% show, c("person", "x", "y")])
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
# unbiased cell reads as blank. Brown against blue stays separable under the
# common forms of colour blindness.
BIAS_RAMP <- c("#8c510a", "#c88a2e", "#e8d3a0", "#f7f7f7",
               "#a8c9ee", "#3d7fc1", "#0d366b")

# Label ink flips to white on the dark ends of the ramp.
ink_for <- function(pos) ifelse(!is.na(pos) & pos > 0.6, "#ffffff", "#0b0b0b")

# One heatmap: mean signed bias of one parameter over the grid. Limits are
# symmetric about 0 and taken from this parameter alone, so panels cannot be
# compared by colour; every tile carries its value for that reason.
# coord_equal() keeps the tiles square.
bias_panel <- function(cells, par) {
  x <- cells[cells$parameter == par, ]
  m <- max(abs(x$bias[is.finite(x$bias)]), 0)
  if (!is.finite(m) || m == 0) m <- 1
  x$label <- ifelse(is.finite(x$bias), sprintf("%+.2f", x$bias), "NA")
  x$ink   <- ink_for(abs(x$bias) / m)

  ggplot(x, aes(factor(l_0), factor(prior_sd), fill = bias)) +
    geom_tile(colour = SURFACE, linewidth = 1.5) +
    geom_text(aes(label = label, colour = ink), size = 2.6) +
    scale_colour_identity() +
    scale_fill_gradientn(colours = BIAS_RAMP, limits = c(-m, m),
                         na.value = MISSING) +
    coord_equal() +
    labs(x = "l_0", y = "prior_sd", fill = NULL, title = par) +
    theme_minimal(base_size = 10) +
    theme(panel.grid       = element_blank(),
          plot.title       = element_text(size = 10, hjust = 0.5),
          legend.key.width = grid::unit(0.25, "cm"))
}

plot_bias <- function(cells) {
  arrangeGrob(
    grobs = lapply(PLOT_PARS, function(p) bias_panel(cells, p)),
    ncol  = 1,
    top   = grid::textGrob(
      sprintf("Mean signed bias over %d iterations\n(posterior mean - truth)",
              N_REPS),
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
         title = heading,
         subtitle = sprintf("l_0 = %g, prior_sd = %g, iteration %d; mean absolute relative bias %.3f; %d HMC draws per infection",
                            ex$l_0, ex$prior_sd, ex$rep, ex$abs_rel_bias,
                            nrow(ex$psi))) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.title       = element_text(size = 11),
          plot.subtitle    = element_text(size = 8),
          strip.text       = element_text(size = 8))
}

# ---- Run --------------------------------------------------------------------

started <- Sys.time()
set.seed(SEED)

model    <- stan_model(STAN_FILE)
grid     <- expand.grid(prior_sd = PRIOR_SD_GRID, l_0 = L_0_GRID)
n_fits   <- N_REPS * nrow(grid)
draw_idx <- sort(sample(N_CHAINS * (N_ITER %/% 2), N_DRAWS))

score_rows <- list()
diag_rows  <- list()
examples   <- list()   # one small payload per successful fit
fit_no     <- 0

for (rep in seq_len(N_REPS)) {
  sim <- generate_matrix()

  for (l_0 in L_0_GRID) {
    # Discovery depends on the matrix and l_0 only, so it is done once and
    # shared by every prior_sd, as are the N_SHOW example infections.
    sampling_out <- sample_rows(sim$mat, l_0)
    fitted       <- vapply(sampling_out$records,
                           function(r) r$person[1], numeric(1))
    stopifnot(length(fitted) >= N_SHOW)
    show <- sort(fitted[sample.int(length(fitted), N_SHOW)])

    for (prior_sd in PRIOR_SD_GRID) {
      fit_no <- fit_no + 1
      prep   <- build_stan_data(sampling_out$records, sampling_out$brackets,
                                prior_sd)

      cat(sprintf("\n=== fit %d/%d | iteration %d | l_0 = %g | prior_sd = %g ===\n",
                  fit_no, n_fits, rep, l_0, prior_sd))
      cat(sprintf("infected %d | pool positive %d | missed %d | fitted %d | assays %d | psi bound clamped %d\n",
                  sim$n_infx, sampling_out$n_pool, sampling_out$n_missed,
                  prep$data$J, prep$data$N, sampling_out$n_clamped))

      fit <- tryCatch(fit_model(model, prep$data, seed = SEED + fit_no),
                      error = function(e) {
                        warning(sprintf("iteration %d, l_0 = %g, prior_sd = %g failed: %s",
                                        rep, l_0, prior_sd, conditionMessage(e)))
                        NULL
                      })
      if (is.null(fit)) next

      sc   <- score(fit)
      arb  <- abs_rel_bias(sc)
      ndiv <- n_divergent(fit)
      key  <- data.frame(rep = rep, l_0 = l_0, prior_sd = prior_sd)

      score_rows[[fit_no]] <- cbind(key, sc, abs_rel_bias = arb,
                                    n_fitted = prep$data$J,
                                    n_missed = sampling_out$n_missed)
      diag_rows[[fit_no]]  <- cbind(key, diagnostics(fit), n_divergent = ndiv)
      examples[[fit_no]]   <- c(list(rep = rep, l_0 = l_0, prior_sd = prior_sd,
                                     abs_rel_bias = arb),
                                example_payload(fit, prep, show, draw_idx))

      cat(sprintf("mean absolute relative bias %.4f | divergences %d\n",
                  arb, ndiv))

      rm(fit); invisible(gc())
    }
  }
}

scores <- do.call(rbind, score_rows)
diags  <- do.call(rbind, diag_rows)
if (is.null(scores)) stop("every fit failed; nothing to summarise")

write.csv(scores, SCORE_FILE, row.names = FALSE)
write.csv(diags,  DIAG_FILE,  row.names = FALSE)

# ---- Summaries --------------------------------------------------------------

# Mean signed bias per cell and parameter, for the heatmaps. Merged onto the
# full grid so a cell whose fits all failed shows as missing, not as a gap.
cell_bias <- aggregate(bias ~ l_0 + prior_sd + parameter, data = scores,
                       FUN = mean)
full      <- merge(grid, data.frame(parameter = PLOT_PARS))
cell_bias <- merge(full, cell_bias, all.x = TRUE)

# One row per fit, then the mean per cell, for choosing the examples.
per_fit  <- unique(scores[, c("rep", "l_0", "prior_sd", "abs_rel_bias")])
per_cell <- aggregate(abs_rel_bias ~ l_0 + prior_sd, data = per_fit, FUN = mean)

pick_example <- function(which_fn) {
  cell <- per_cell[which_fn(per_cell$abs_rel_bias), ]
  fits <- per_fit[per_fit$l_0 == cell$l_0 & per_fit$prior_sd == cell$prior_sd, ]
  best <- fits[which_fn(fits$abs_rel_bias), ]
  hit  <- Filter(function(e) !is.null(e) && e$rep == best$rep &&
                   e$l_0 == best$l_0 && e$prior_sd == best$prior_sd, examples)
  c(hit[[1]], list(cell_abs_rel_bias = cell$abs_rel_bias))
}

low  <- pick_example(which.min)
high <- pick_example(which.max)
saveRDS(list(low = low, high = high), EXAMPLES_FILE)

# ---- Figure -----------------------------------------------------------------

fig <- arrangeGrob(
  plot_bias(cell_bias),
  plot_draws(low,  "Least biased iteration of the least biased cell"),
  plot_draws(high, "Most biased iteration of the most biased cell"),
  layout_matrix = cbind(c(1, 1), c(2, 3)),
  widths = c(1, 2.6))

ggsave(FIG_FILE, fig, width = 15, height = 8.5, dpi = 150, bg = SURFACE)

cat(sprintf("\nleast biased cell: l_0 = %g, prior_sd = %g (cell mean %.4f), iteration %d (%.4f)\n",
            low$l_0, low$prior_sd, low$cell_abs_rel_bias, low$rep,
            low$abs_rel_bias))
cat(sprintf("most biased cell:  l_0 = %g, prior_sd = %g (cell mean %.4f), iteration %d (%.4f)\n",
            high$l_0, high$prior_sd, high$cell_abs_rel_bias, high$rep,
            high$abs_rel_bias))
cat(sprintf("%d of %d fits succeeded\n", nrow(per_fit), n_fits))
cat(sprintf("wrote %s, %s, %s and %s\n",
            SCORE_FILE, DIAG_FILE, EXAMPLES_FILE, FIG_FILE))
cat(sprintf("total time %.1f min\n",
            as.numeric(difftime(Sys.time(), started, units = "mins"))))
