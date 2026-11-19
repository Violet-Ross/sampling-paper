# Case Study 1 worker script -- run once per l_0 condition.
# Called by the SLURM array job; reads its l_0 from param_grid.csv using the
# SLURM_ARRAY_TASK_ID environment variable.
#
# For this l_0:
#   1. Generate ONE synthetic sample matrix in which each infected person is
#      vaccinated with probability P_VAX. Vaccinated and unvaccinated
#      infections differ only in the mean of the peak viral load.
#   2. Infection discovery: pool each row, then assay the positive rows on days
#      0, l_0, 2*l_0, ... and keep the assays within W_DAYS of the first and
#      last positive.
#   3. Sort the discovered infections into V (vaccinated) and U (unvaccinated)
#      and fit the BHM in src/bhm.stan to each group separately.
#   4. Subtract the mu_V posterior draws of the two fits entrywise (draw i of V
#      minus draw i of U) and save the resulting vector.
#   5. Save the signed error (posterior mean minus generating value) of
#      each fit against its own group's generating values.
#
# The matrix is drawn from the same seed in every task, so all l_0 conditions
# are discovered on the identical matrix.
#
# The generation, discovery and Stan-data code is copied from
# scripts/nonadaptive_experiment/nonadaptive_experiment.R, which cannot be
# sourced because it fits at top level.
#
# Run from the repository root, after setup_param_array.R:
#   SLURM_ARRAY_TASK_ID=1 Rscript scripts/question1/run_condition.R
# or on the cluster:
#   sbatch scripts/question1/run_experiment.sh

library(rstan)
library(readr)

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# ---- Settings ---------------------------------------------------------------

STAN_FILE <- "src/bhm.stan"

# Generating values. The two groups differ in the mean peak viral load only.
TRUE_V_VAX   <- 5      # peak viral load, vaccinated
TRUE_V_UNVAX <- 6      # peak viral load, unvaccinated
TRUE_DELT_P  <- 3.5    # proliferation time, days
TRUE_DELT_C  <- 10     # clearance time, days
SIGMA_IND    <- 0.1    # person-to-person spread
SIGMA_EPS    <- 0.1    # observation noise

M            <- 500    # days (matrix columns)
N_PEOPLE     <- 500    # people (matrix rows)
PREV         <- 0.3    # fraction of people infected
P_VAX        <- 0.5    # fraction of people vaccinated

# Assays kept this many days either side of the positives. sample_rows needs
# W_DAYS >= l_0, and the l_0 grid runs to 20.
W_DAYS       <- 20

POS_THRESH   <- 0      # an assay is positive above this value
PSI_SHAPE    <- 2      # Beta(a, a) prior on the peak time within its bracket

# Priors. Both fits get the SAME centre for mu_V, the midpoint of the two
# generating peaks, so the priors carry no difference between the groups. The
# two durations are shared by the groups and centred on their generating values.
PRIOR_SD     <- 2
PRIOR_MU_V   <- mean(c(TRUE_V_VAX, TRUE_V_UNVAX))

N_CHAINS     <- 4
N_ITER       <- 4000   # per chain; Stan warms up for half
SEED         <- 1

# ---- Read condition for this task -------------------------------------------

task_id    <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
param_grid <- read_csv("throughput/question1/param_grid.csv", show_col_types = FALSE)
stopifnot(!is.na(task_id), task_id >= 1, task_id <= nrow(param_grid))
l_0        <- param_grid$l_0[task_id]
stopifnot(W_DAYS >= l_0)

cat(sprintf("Task %d | l_0 = %d | V vaccinated = %.2f | V unvaccinated = %.2f | prior_mu_V = %.2f | prior_sd = %.2f\n",
            task_id, l_0, TRUE_V_VAX, TRUE_V_UNVAX, PRIOR_MU_V, PRIOR_SD))

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
  is_vax   <- rbinom(N_PEOPLE, 1, P_VAX)
  days     <- 0:(M - 1)
  infected <- which(is_infx == 1)

  truth <- data.frame(person = infected, vaccinated = is_vax[infected] == 1,
                      V = NA_real_, delt_p = NA_real_, delt_c = NA_real_)

  for (k in seq_along(infected)) {
    i      <- infected[k]
    V      <- rlnorm_at(if (truth$vaccinated[k]) TRUE_V_VAX else TRUE_V_UNVAX,
                        SIGMA_IND)
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
                   prior_mu_V      = PRIOR_MU_V,
                   prior_mu_delt_p = TRUE_DELT_P,
                   prior_mu_delt_c = TRUE_DELT_C))
}

# ---- Fit --------------------------------------------------------------------

fit_model <- function(model, stan_data, seed) {
  sampling(model, data = stan_data,
           chains = N_CHAINS, iter = N_ITER, seed = seed, refresh = 0,
           pars = c("V_raw", "delt_p_raw", "delt_c_raw"),
           include = FALSE)
}

# ---- Scoring ----------------------------------------------------------------

SCORE_PARS <- c("mu_V", "mu_delt_p", "mu_delt_c", "mu_duration")

# Signed error of each population mean, as in the other experiments: the
# posterior mean minus the generating value. `true_v` is the generating peak
# of the group the fit belongs to.
compute_bias <- function(fit, true_v) {
  post  <- rstan::extract(fit, pars = SCORE_PARS)
  truth <- c(mu_V        = true_v,
             mu_delt_p   = TRUE_DELT_P,
             mu_delt_c   = TRUE_DELT_C,
             mu_duration = TRUE_DELT_P + TRUE_DELT_C)
  out <- lapply(SCORE_PARS, function(p) mean(post[[p]] - truth[[p]]))
  names(out) <- paste0(SCORE_PARS, "_bias")
  as.data.frame(out)
}

DIAG_PARS <- c("mu_V", "mu_delt_p", "mu_delt_c",
               "sigma_log_V", "sigma_log_delt_p", "sigma_log_delt_c",
               "sigma_eps", "mu_duration")

max_rhat <- function(fit) max(summary(fit, pars = DIAG_PARS)$summary[, "Rhat"])

n_divergent <- function(fit) {
  sum(sapply(get_sampler_params(fit, inc_warmup = FALSE),
             function(p) sum(p[, "divergent__"])))
}

# ---- Generate and discover --------------------------------------------------

# The same seed in every task: one matrix shared by all l_0 conditions.
set.seed(SEED)
sim <- generate_matrix()

sampling_out <- sample_rows(sim$mat, l_0)
found        <- vapply(sampling_out$records, function(r) r$person[1], numeric(1))
found_vax    <- sim$truth$vaccinated[match(found, sim$truth$person)]
stopifnot(!anyNA(found_vax), any(found_vax), any(!found_vax))

cat(sprintf("infected %d (vaccinated %d) | pool positive %d | missed %d | discovered V %d, U %d | psi bound clamped %d\n",
            sim$n_infx, sum(sim$truth$vaccinated), sampling_out$n_pool,
            sampling_out$n_missed, sum(found_vax), sum(!found_vax),
            sampling_out$n_clamped))

model <- stan_model(STAN_FILE)

# One group: build its Stan data, fit, and score against its own truth.
run_arm <- function(label, keep, true_v, seed) {
  cat(sprintf("  Fitting %s...\n", label))
  prep <- build_stan_data(sampling_out$records[keep],
                          sampling_out$brackets[keep], PRIOR_SD)
  fit  <- fit_model(model, prep$data, seed)
  out  <- list(bias        = compute_bias(fit, true_v),
               peak_draws  = as.numeric(rstan::extract(fit, pars = "mu_V")$mu_V),
               n_fitted    = prep$data$J,
               n_assays    = prep$data$N,
               max_rhat    = max_rhat(fit),
               n_divergent = n_divergent(fit))
  cat(sprintf("  %s: fitted %d infections, %d assays | mu_V bias %+.3f | max Rhat %.3f | divergences %d\n",
              label, out$n_fitted, out$n_assays, out$bias$mu_V_bias,
              out$max_rhat, out$n_divergent))
  out
}

# ---- V and U arms -----------------------------------------------------------

arm_v <- run_arm("V", found_vax,  TRUE_V_VAX,   seed = SEED + task_id)
arm_u <- run_arm("U", !found_vax, TRUE_V_UNVAX, seed = SEED + 1000 + task_id)

# ---- Entrywise posterior difference (peak viral load) -----------------------

n_draws <- min(length(arm_v$peak_draws), length(arm_u$peak_draws))
stopifnot(length(arm_v$peak_draws) == length(arm_u$peak_draws))
peak_diff <- arm_v$peak_draws[1:n_draws] - arm_u$peak_draws[1:n_draws]

# ---- Write this task's results ----------------------------------------------

diff_df <- data.frame(
  task_id   = task_id,
  l_0       = l_0,
  draw      = seq_len(n_draws),
  peak_diff = peak_diff
)
diff_path <- sprintf("throughput/question1/results_diff/diff_task_%04d.csv", task_id)
write_csv(diff_df, diff_path)

arm_cols <- function(prefix, arm) {
  out <- cbind(arm$bias,
               n_fitted    = arm$n_fitted,
               n_assays    = arm$n_assays,
               max_rhat    = arm$max_rhat,
               n_divergent = arm$n_divergent)
  names(out) <- paste0(prefix, "_", names(out))
  out
}

summary_row <- cbind(
  data.frame(task_id      = task_id,
             l_0          = l_0,
             true_v_vax   = TRUE_V_VAX,
             true_v_unvax = TRUE_V_UNVAX,
             prior_sd     = PRIOR_SD,
             prior_mu_V   = PRIOR_MU_V),
  arm_cols("v", arm_v),
  arm_cols("u", arm_u),
  data.frame(n_missed       = sampling_out$n_missed,
             peak_diff_mean = mean(peak_diff),
             peak_diff_sd   = sd(peak_diff))
)
summary_path <- sprintf("throughput/question1/results_summary/summary_task_%04d.csv", task_id)
write_csv(summary_row, summary_path)

cat(sprintf("peak difference (V - U): mean %.3f, sd %.3f | P(diff > 0) = %.4f\n",
            mean(peak_diff), sd(peak_diff), mean(peak_diff > 0)))
cat(sprintf("Written to %s and %s\n", summary_path, diff_path))
