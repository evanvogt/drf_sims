##########
# title: figures for the thesis chapter - continuous outcomes, missing covariates
##########
# Labels, palette and summaries come from R/figures.R. This script carries
# only the paths and this study's filters.

library(here)
library(patchwork)
source(here("R", "figures.R"))

# paths
path <- here()
res_path <- file.path(dirname(path), "results", "missing", "continuous")
fig_path <- file.path(dirname(path), "results", "thesis_figures", "miss_cts")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

metrics <- readRDS(file.path(res_path, "cts_miss_metrics.RDS"))

# the reduced scenario set, renamed for the chapter
MISS_SCENARIO_LABELS <- c(`1` = "Null", `5` = "Simple", `3` = "Complex",
                          `4` = "Non-linear")

metrics <- metrics %>%
  filter((scenario %in% c(1, 3, 4, 5) & type == "both") | (scenario == 1)) %>%
  apply_labels(MISS_SCENARIO_LABELS) %>%
  mutate(n = factor(n), prop = factor(prop))

metrics_summary <- summarise_metrics(
  metrics,
  c("scenario", "n", "type", "prop", "mechanism", "method", "model"),
  cols = c(bias = "bias", mse = "mse", corr = "corr")
)

bias_sum_plot <- point_range_plot(metrics_summary, "bias", "Bias")
save_fig("miss_cts_bias.png", fig_path)

mse_sum_plot <- point_range_plot(metrics_summary, "mse", "MSE",
                                 facet_scales = "free_y")
save_fig("miss_cts_mse.png", fig_path)

b_new <- bias_sum_plot + theme(axis.text.x = element_blank())
(b_new / mse_sum_plot) + plot_layout(guides = "collect")
ggsave("miss_cts_bias_mse.png", path = fig_path, width = 21, height = 15, units = "cm")

# every panel on its own y scale, so the variation within a scenario shows -
# the shared-scale rows above compress it. Not saved.
point_range_plot(metrics_summary, "bias", "Bias", independent_y = TRUE)
