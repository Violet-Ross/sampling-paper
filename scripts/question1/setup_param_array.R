# Run this once, from the repository root, before submitting the job array.
# Writes param_grid.csv, which the array jobs index into, and compiles the
# Stan model.
#
# The experiment fixes prior_sd and the generating values as constants inside
# run_condition.R and only loops over l_0. One array task per l_0.

library(readr)
library(rstan)

rstan_options(auto_write = TRUE)

param_grid <- data.frame(l_0 = 1:20)

dir.create("throughput/question1/results_summary", showWarnings = FALSE, recursive = TRUE)
dir.create("throughput/question1/results_diff",    showWarnings = FALSE, recursive = TRUE)

write_csv(param_grid, "throughput/question1/param_grid.csv")
cat(sprintf("Wrote %d conditions to throughput/question1/param_grid.csv\n", nrow(param_grid)))

# Compiled here so the array tasks all load the cached src/bhm.rds. Without
# this they would each compile the model and write that file at the same time.
invisible(stan_model("src/bhm.stan"))
cat("Compiled src/bhm.stan\n")
