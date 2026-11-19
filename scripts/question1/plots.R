# Draws the Case Study 1 figure from the combined results, without refitting.
#
# Inputs (written by combine_results.R):
#   throughput/question1/peak_diff.csv     mu_V draws of V minus those of U
#   throughput/question1/bias_summary.csv  signed error of each fit (posterior
#                                          mean minus truth), one row per l_0
#
# Output:
#   fig/question1.png
#
# Run from the repository root:
#   Rscript scripts/question1/plots.R

library(tidyverse)
library(latex2exp)

diffs <- read_csv("throughput/question1/peak_diff.csv", show_col_types = FALSE)
bias  <- read_csv("throughput/question1/bias_summary.csv", show_col_types = FALSE)

# A difference is detected at an l_0 when less than 5% of the posterior
# differences fall on one side of 0.
detect_summary <- diffs %>%
  group_by(l_0) %>%
  summarise(found_diff = min(mean(peak_diff > 0), mean(peak_diff < 0)) < 0.05)

bias_long <- bias %>%
  select(l_0, v_mu_V_bias, u_mu_V_bias) %>%
  pivot_longer(cols = 2:3, names_to = "param", values_to = "bias")

# The error is signed, so the y range straddles 0; it always includes 0.
y_lo <- min(bias_long$bias, 0)
y_hi <- max(bias_long$bias, 0)

p <- ggplot(bias_long)

# Shade the l_0 range where the difference is detected, and label either side
# of its upper edge.
if (any(detect_summary$found_diff)) {
  shade_range <- detect_summary %>%
    filter(found_diff) %>%
    summarise(xmin = min(l_0), xmax = max(l_0))

  arrow_y <- y_lo + 0.95 * (y_hi - y_lo)
  label_y <- y_lo + 1.00 * (y_hi - y_lo)

  p <- p +
    annotate("rect",
             xmin = shade_range$xmin, xmax = shade_range$xmax,
             ymin = -Inf, ymax = Inf,
             fill = "grey70", alpha = 0.3) +
    annotate("segment",
             xend = shade_range$xmax - 2.2, x = shade_range$xmax - 0.5,
             y = arrow_y, yend = arrow_y,
             arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.6) +
    annotate("text",
             x = shade_range$xmax - 2.8, y = label_y,
             label = "difference\ndetected",
             hjust = 0, size = 3.2) +
    annotate("segment",
             xend = shade_range$xmax + 2.2, x = shade_range$xmax + 0.5,
             y = arrow_y, yend = arrow_y,
             arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.6) +
    annotate("text",
             x = shade_range$xmax + 0.5, y = label_y,
             label = "difference\nundetected",
             hjust = 0, size = 3.2)
}

p <- p +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(aes(x = l_0, y = bias, color = param)) +
  theme_classic() +
  scale_color_manual(
    values = c("v_mu_V_bias" = "#F8766D", "u_mu_V_bias" = "#00BFC4"),
    labels = c("v_mu_V_bias" = "vaccinated population", "u_mu_V_bias" = "unvaccinated population"),
    name = "population to which \nthe BHM is fit"
  ) +
  xlab(TeX("$\\ell_0$")) +
  ylab(TeX("signed error of $\\mu_V$"))

dir.create("fig", showWarnings = FALSE)
ggsave("fig/question1.png", p, dpi = 600)
