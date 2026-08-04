##########
# title: collect up results - binary outcome CIs
##########
# The grid and the results path come from the study config. Produces one row per
# parameter combination with a list-column of per-run results, which is the shape
# the metrics script unnests.

library(here)
source(here("confidence_intervals/binary/bin_ci_config.R"))

workers <- 2

all_results_df <- get_results(study, workers = workers)

saveRDS(all_results_df, file.path(study$res_path, "ci_bin_all.RDS"))
print("Collection complete!")
