# Run this after all array jobs have finished.
# Combines the per-task signed-error summary CSVs and the per-task peak-viral-load
# posterior-difference CSVs into two combined files (20 tasks expected,
# one per l_0).

library(readr)
library(dplyr)

N_TASKS <- 20
all_ids <- 1:N_TASKS

# ── Signed-error summaries ─────────────────────────────────────────────────────────────

summary_files <- list.files("throughput/question1/results_summary", pattern = "summary_task_.*\\.csv",
                            full.names = TRUE)

cat(sprintf("Found %d summary files (expected %d)\n", length(summary_files), N_TASKS))

found_summary_ids   <- as.integer(gsub(".*task_(\\d+)\\.csv", "\\1", summary_files))
missing_summary_ids <- setdiff(all_ids, found_summary_ids)

if (length(missing_summary_ids) > 0) {
  cat(sprintf("WARNING: %d summary tasks missing: %s\n",
              length(missing_summary_ids), paste(missing_summary_ids, collapse = ", ")))
} else {
  cat("All summary tasks accounted for.\n")
}

summary_list <- bind_rows(lapply(summary_files, read_csv, show_col_types = FALSE))
write_csv(summary_list, "throughput/question1/bias_summary.csv")
cat(sprintf("Combined %d summary rows written to throughput/question1/bias_summary.csv\n",
            nrow(summary_list)))

# ── Peak viral load posterior differences ──────────────────────────────────────

diff_files <- list.files("throughput/question1/results_diff", pattern = "diff_task_.*\\.csv",
                         full.names = TRUE)

cat(sprintf("Found %d diff files (expected %d)\n", length(diff_files), N_TASKS))

found_diff_ids   <- as.integer(gsub(".*task_(\\d+)\\.csv", "\\1", diff_files))
missing_diff_ids <- setdiff(all_ids, found_diff_ids)

if (length(missing_diff_ids) > 0) {
  cat(sprintf("WARNING: %d diff tasks missing: %s\n",
              length(missing_diff_ids), paste(missing_diff_ids, collapse = ", ")))
} else {
  cat("All diff tasks accounted for.\n")
}

diff_list <- bind_rows(lapply(diff_files, read_csv, show_col_types = FALSE))
write_csv(diff_list, "throughput/question1/peak_diff.csv")
cat(sprintf("Combined %d diff rows written to throughput/question1/peak_diff.csv\n",
            nrow(diff_list)))