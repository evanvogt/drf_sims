##########
# title: collect up results - correlated binary outcome CIs
##########
# The grid and the results path come from the study config. Produces one row per
# parameter combination with a list-column of per-run results, which is the shape
# the metrics script unnests.

library(here)
source(here("sample_size/correlated/confidence_intervals/binary/bin_corr_ci_config.R"))

workers <- 2

all_results_df <- get_results(study, workers = workers)

saveRDS(all_results_df, file.path(study$res_path, "bin_corr_ci_all.RDS"))
print("Collection complete!")
