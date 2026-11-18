# Adaptive sampling in the most biased cell of the non-adaptive experiment.
#
# The non-adaptive experiment's most biased cell is l_0 = 5, prior_sd = 1.
# This script runs N_REPS independent iterations in that cell. Each iteration
# draws a fresh synthetic matrix, fits the l_0 design (round 0), then runs the
# adaptive loop of vk-modeling's adaptive_grid/run_adaptive_condition.R for
# N_ROUNDS rounds. Each round scans the predictive entropy of every unassayed
# day in each infection's record, adds the matrix observation at the
# highest-entropy day (one new assay per infection), rebrackets and refits
# from scratch.
#
# After the starting fit and after every round, three diagnostics are
# recorded: signed bias against the generating values, Rhat and n_eff.
#
#
# The model is src/bhm.stan throughout. The sampling, bracketing and entropy
# functions are copied from run_adaptive_condition.R; generation, Stan-data
# assembly, fit, scoring and diagnostics from nonadaptive_experiment.R.
#
# Outputs, every row tagged with rep (the iteration) and round:
#   throughput/adaptive_worst_cell_scores.csv       bias per parameter
#   throughput/adaptive_worst_cell_diagnostics.csv  Rhat, n_eff, divergences
#   throughput/adaptive_worst_cell_additions.csv    the assays added
#   throughput/adaptive_worst_cell_rounds.csv       one row per fit
#   fig/adaptive_worst_cell.png
#
# Run from the repository root:
#   Rscript scripts/adaptive_worst_cell/adaptive_worst_cell.R
# or on the cluster:
#   sbatch scripts/adaptive_worst_cell/adaptive_worst_cell.sh

library(rstan)
library(ggplot2)
library(gridExtra)

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

if (!interactive()) pdf(NULL)

# ---- Settings ---------------------------------------------------------------

L_0         <- 5      # interval-sampling spacing of the starting design, days
PRIOR_SD    <- 1      # prior standard deviation

N_REPS      <- 10     # iterations, each on its own synthetic matrix
N_ROUNDS    <- 5      # add-and-refit cycles after each starting fit

W_DAYS      <- 10

# Generating values, as nonadaptive_experiment.R. These are also the prior
# centres passed to Stan.
TRUE_V      <- 4
TRUE_DELT_P <- 3.5
TRUE_DELT_C <- 10
SIGMA_IND   <- 0.1
SIGMA_EPS   <- 0.1

M           <- 500
N_PEOPLE    <- 500
PREV        <- 0.3

POS_THRESH  <- 0
PSI_SHAPE   <- 2

N_CHAINS    <- 4
N_ITER      <- 2000
SEED        <- 1

# Entropy scan, as run_adaptive_condition.R.
Y_SD_MULT   <- 6      # integral runs this many sigma_eps beyond the extreme means
SIGMA_FLOOR <- 1e-3   # draws with sigma_eps below this leave the predictive mixture

STAN_FILE   <- "src/bhm.stan"
SCORE_FILE  <- "throughput/adaptive_worst_cell_scores.csv"
DIAG_FILE   <- "throughput/adaptive_worst_cell_diagnostics.csv"
ADD_FILE    <- "throughput/adaptive_worst_cell_additions.csv"
ROUND_FILE  <- "throughput/adaptive_worst_cell_rounds.csv"
FIG_FILE    <- "fig/adaptive_worst_cell.png"

# Resuming. START_REP=k in the environment keeps the saved rows of iterations
# before k and reruns from iteration k. The earlier iterations' draws on R's
# random stream are replayed without fitting (each matrix, then one tie-break
# draw per saved addition), so the stream and the Stan seeds are where an
# uninterrupted run would have them. Each regenerated design is checked
# against the saved assay count, which fails if the replay has drifted.
START_REP   <- as.integer(Sys.getenv("START_REP", "1"))

stopifnot(W_DAYS >= L_0, N_REPS >= 1, N_ROUNDS >= 1,
          START_REP >= 1, START_REP <= N_REPS)
for (f in c(SCORE_FILE, FIG_FILE))
  dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)

# ---- Curve ------------------------------------------------------------------

witch <- function(t, V, delt_p, delt_c) {
  ifelse(t > -delt_p & t <= 0, V * (1 + t / delt_p),
  ifelse(t >  0      & t <  delt_c, V * (1 - t / delt_c), 0))
}

rlnorm_at <- function(median, sigma) median * exp(rnorm(1, 0, sigma / median))

# ---- Generation -------------------------------------------------------------

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

# ---- Sampling ---------------------------------------------------------------

pool_positive <- function(row) any(row > 0)

interval_days <- function(l_0) seq(0, M - 1, by = l_0)

# As nonadaptive_experiment.R's sample_rows, plus an `org` column on each
# record: the absolute matrix day of the record's first entry.
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

    if (length(pos) == 0) {
      n_missed <- n_missed + 1
      next
    }

    lo   <- days[min(pos)] - w_days
    hi   <- days[max(pos)] + w_days
    keep <- which(days >= lo & days <= hi)
    org  <- days[keep[1]]

    records[[length(records) + 1]] <- data.frame(
      person = i,
      x      = days[keep] - org,
      y      = vals[keep],
      org    = org
    )

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

# The psi bracket of one record from the record alone: the last negative test
# before the first positive (0, the record start, if there is none) and the
# first negative after the last positive (the record end if there is none).
# On the records sample_rows builds this reproduces sample_rows' own brackets,
# which the run section asserts; it is what lets the augmented records be
# bracketed the same way.
bracket_record <- function(rec) {
  rec <- rec[order(rec$x), ]
  pos <- which(rec$y > POS_THRESH)
  stopifnot(length(pos) > 0)
  n    <- nrow(rec)
  lo_i <- min(pos) - 1L
  hi_i <- max(pos) + 1L
  data.frame(person = rec$person[1],
             psi_lo = if (lo_i < 1L) 0 else rec$x[lo_i],
             psi_hi = if (hi_i > n)  rec$x[n] else rec$x[hi_i])
}

build_stan_data <- function(records, brackets, prior_sd) {
  df      <- do.call(rbind, records)
  persons <- unique(df$person)
  df$id   <- match(df$person, persons)

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

# ---- Entropy ----------------------------------------------------------------

# Mixture-of-normals predictive density at each y, averaged over draws.
# `mu` and `sigma` are one value per draw; `y` is the vector integrate
# supplies. Rows are draws, columns are y; `mu` and `sigma` recycle down the
# columns because their length equals nrow.
posterior_predictive <- function(mu, sigma, y) {
  Y    <- matrix(y, nrow = length(mu), ncol = length(y), byrow = TRUE)
  dens <- dnorm(Y, mean = mu, sd = sigma)
  colMeans(dens)
}

# Differential entropy of the predictive at one time t, over the whole real
# line. mu depends on t and the draws only, so it is computed once here, not
# once per integrand call, and the integration limits come from it: the
# mixture has no mass beyond Y_SD_MULT sigma of its extreme means. A
# quadrature failure returns NA with a warning rather than aborting.
H <- function(draws, t) {
  mu    <- witch(t - draws$psi, draws$V, draws$delt_p, draws$delt_c)
  reach <- Y_SD_MULT * max(draws$sigma_eps)
  lower <- min(mu) - reach
  upper <- max(mu) + reach

  f <- function(y) {
    p <- posterior_predictive(mu, draws$sigma_eps, y)
    ifelse(p > 0, -p * log(p), 0)
  }

  tryCatch(
    integrate(f, lower = lower, upper = upper)$value,
    error = function(e) {
      warning(sprintf("person %d, t = %s: integrate failed (%s); returning NA",
                      draws$person, format(t), conditionMessage(e)),
              call. = FALSE)
      NA_real_
    }
  )
}

# sample() on a length-one vector samples from 1:x; this never does.
resample <- function(v) v[sample.int(length(v), 1)]

# The row choose_addition returns when nothing is added.
no_addition <- function(p, org, span, n_cand, n_failed) {
  data.frame(person = p, org = org, span = span,
             n_candidates = n_cand, n_failed = n_failed,
             n_tied = NA_integer_, x_added = NA_real_,
             day_added = NA_real_, y_added = NA_real_,
             entropy_max = NA_real_)
}

# For one record: entropy at every unassayed integer day in [0, span], then
# the matrix observation at the highest-entropy day. Every post-warmup draw,
# chains pooled, enters the predictive mixture. Returns a one-row data frame;
# x_added is NA if the record has no unassayed day left or none could be
# evaluated.
choose_addition <- function(rec, j, post, mat) {
  p     <- rec$person[1]
  org   <- rec$org[1]
  span  <- max(rec$x)
  cand  <- setdiff(0:span, rec$x)

  if (length(cand) == 0) return(no_addition(p, org, span, 0L, 0L))

  draws <- list(person    = p,
                psi       = post$psi[, j],
                V         = post$V[, j],
                delt_p    = post$delt_p[, j],
                delt_c    = post$delt_c[, j],
                sigma_eps = post$sigma_eps)

  ent <- vapply(cand, function(t) H(draws, t), numeric(1))
  ok  <- !is.na(ent)

  if (!any(ok)) {
    warning(sprintf("person %d: entropy failed at every candidate day; no addition", p),
            call. = FALSE)
    return(no_addition(p, org, span, length(cand), sum(!ok)))
  }

  best  <- which(ok & ent == max(ent[ok]))
  x_new <- cand[resample(best)]
  day   <- org + x_new

  data.frame(person = p, org = org, span = span,
             n_candidates = length(cand), n_failed = sum(!ok),
             n_tied = length(best), x_added = x_new, day_added = day,
             y_added = mat[p, day + 1], entropy_max = max(ent[ok]))
}

# The record with the chosen observation appended, sorted by x.
augment_record <- function(rec, add) {
  if (is.na(add$x_added)) return(rec)
  stopifnot(!(add$x_added %in% rec$x))
  rec <- rbind(rec, data.frame(person = rec$person[1], x = add$x_added,
                               y = add$y_added, org = rec$org[1]))
  rec[order(rec$x), ]
}

# ---- Run --------------------------------------------------------------------

started <- Sys.time()
elapsed <- function() as.numeric(difftime(Sys.time(), started, units = "mins"))

set.seed(SEED)
model <- stan_model(STAN_FILE)

tag <- function(df, rep, r) cbind(rep = rep, round = r, df)
score_rows <- list(); diag_rows <- list(); add_rows <- list(); round_rows <- list()

saved_rounds <- NULL
if (START_REP > 1) {
  before       <- function(path) { x <- read.csv(path); x[x$rep < START_REP, ] }
  saved_rounds <- before(ROUND_FILE)
  stopifnot(nrow(saved_rounds) == (START_REP - 1) * (N_ROUNDS + 1))
  score_rows[[1]] <- before(SCORE_FILE)
  diag_rows[[1]]  <- before(DIAG_FILE)
  saved_adds      <- before(ADD_FILE)
  add_rows[[1]]   <- saved_adds
  round_rows[[1]] <- saved_rounds
  cat(sprintf("resuming at iteration %d; kept %d saved fits\n",
              START_REP, nrow(saved_rounds)))
}

write_all <- function() {
  write.csv(do.call(rbind, score_rows), SCORE_FILE, row.names = FALSE)
  write.csv(do.call(rbind, diag_rows),  DIAG_FILE,  row.names = FALSE)
  write.csv(do.call(rbind, round_rows), ROUND_FILE, row.names = FALSE)
  if (length(add_rows) > 0)
    write.csv(do.call(rbind, add_rows), ADD_FILE, row.names = FALSE)
}

# Score and diagnose one fit, and rewrite the CSVs so a run that dies part-way
# keeps the fits it finished.
record_round <- function(rep, r, fit, prep, n_dropped) {
  sc   <- score(fit)
  d    <- diagnostics(fit)
  ndiv <- n_divergent(fit)
  arb  <- abs_rel_bias(sc)
  k    <- length(score_rows) + 1

  score_rows[[k]] <<- tag(cbind(sc, abs_rel_bias = arb), rep, r)
  diag_rows[[k]]  <<- tag(cbind(d, n_divergent = ndiv), rep, r)
  round_rows[[k]] <<- tag(data.frame(n_obs = prep$data$N,
                                     n_fitted = prep$data$J,
                                     abs_rel_bias = arb,
                                     max_Rhat = max(d$Rhat),
                                     min_n_eff = min(d$n_eff),
                                     divergences = ndiv,
                                     n_draws_dropped = n_dropped,
                                     elapsed_min = elapsed()), rep, r)
  write_all()

  cat(sprintf("iteration %d | round %d | assays %d | abs rel bias %.4f | max Rhat %.2f | min n_eff %.0f | divergences %d | %.1f min\n",
              rep, r, prep$data$N, arb, max(d$Rhat), min(d$n_eff), ndiv,
              elapsed()))
  flush(stdout())
}

# Stan seeds run SEED + 1, SEED + 2, ... over the fits in the order made, so
# the whole run can be repeated. run_adaptive_condition.R leaves them unset.
fit_no <- 0

for (rep in seq_len(N_REPS)) {
  sim <- generate_matrix()

  # -- Round 0: the l_0 design ------------------------------------------------

  sampling_out <- sample_rows(sim$mat, L_0)
  records      <- sampling_out$records
  stopifnot(isTRUE(all.equal(do.call(rbind, lapply(records, bracket_record)),
                             do.call(rbind, sampling_out$brackets),
                             check.attributes = FALSE)))

  prep    <- build_stan_data(records, sampling_out$brackets, PRIOR_SD)
  persons <- prep$persons

  cat(sprintf("\n=== iteration %d/%d | infected %d | missed %d | fitted %d | assays %d | psi bound clamped %d ===\n",
              rep, N_REPS, sim$n_infx, sampling_out$n_missed,
              prep$data$J, prep$data$N, sampling_out$n_clamped))

  if (rep < START_REP) {
    stopifnot(prep$data$N == saved_rounds$n_obs[saved_rounds$rep == rep &
                                                  saved_rounds$round == 0])
    # choose_addition draws once, through resample(), for every addition
    # it makes, in the order the additions were saved.
    tied <- saved_adds$n_tied[saved_adds$rep == rep]
    for (n in tied[!is.na(tied)]) sample.int(n, 1)
    fit_no <- fit_no + N_ROUNDS + 1
    next
  }

  fit_no <- fit_no + 1
  fit    <- fit_model(model, prep$data, seed = SEED + fit_no)
  record_round(rep, 0, fit, prep, NA_integer_)

  # -- Rounds 1..N_ROUNDS -----------------------------------------------------

  for (r in seq_len(N_ROUNDS)) {
    post <- rstan::extract(fit, pars = c("psi", "V", "delt_p", "delt_c", "sigma_eps"))

    # Drop degenerate draws from the mixture: a draw with sigma_eps below
    # SIGMA_FLOOR is a point mass and has no finite differential entropy.
    ok        <- post$sigma_eps >= SIGMA_FLOOR
    n_dropped <- sum(!ok)
    if (n_dropped > 0) {
      post <- list(psi = post$psi[ok, , drop = FALSE], V = post$V[ok, , drop = FALSE],
                   delt_p = post$delt_p[ok, , drop = FALSE],
                   delt_c = post$delt_c[ok, , drop = FALSE],
                   sigma_eps = post$sigma_eps[ok])
    }
    stopifnot(nrow(post$psi) > 0)

    adds <- do.call(rbind, lapply(seq_along(records), function(k) {
      rec <- records[[k]]
      j   <- match(rec$person[1], persons)
      choose_addition(rec, j, post, sim$mat)
    }))
    add_rows[[length(add_rows) + 1]] <- tag(adds, rep, r)

    cat(sprintf("iteration %d | round %d | additions: %d made, %d skipped | %d positive values | %d with ties | %d candidate days failed | %d draws dropped\n",
                rep, r, sum(!is.na(adds$x_added)), sum(is.na(adds$x_added)),
                sum(adds$y_added > POS_THRESH, na.rm = TRUE),
                sum(adds$n_tied > 1, na.rm = TRUE), sum(adds$n_failed),
                n_dropped))

    records <- Map(augment_record, records, split(adds, seq_len(nrow(adds))))
    prep    <- build_stan_data(records, lapply(records, bracket_record), PRIOR_SD)
    stopifnot(identical(prep$persons, persons))

    rm(fit); invisible(gc())
    fit_no <- fit_no + 1
    fit    <- fit_model(model, prep$data, seed = SEED + fit_no)
    record_round(rep, r, fit, prep, n_dropped)
  }

  rm(fit); invisible(gc())
}

scores <- do.call(rbind, score_rows)
rounds <- do.call(rbind, round_rows)

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

# ---- Summary ----------------------------------------------------------------

cat("\nmean over iterations, by round:\n")
print(format(aggregate(cbind(abs_rel_bias, max_Rhat, min_n_eff, divergences) ~ round,
                       data = rounds, FUN = mean), digits = 4),
      row.names = FALSE)
cat(sprintf("\nwrote %s, %s, %s, %s and %s\n",
            SCORE_FILE, DIAG_FILE, ADD_FILE, ROUND_FILE, FIG_FILE))
cat(sprintf("total time %.1f min\n", elapsed()))
