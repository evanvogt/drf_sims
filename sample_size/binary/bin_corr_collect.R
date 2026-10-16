##########
# title: collect up results - correlated binary outcome
##########
# The grid and the results path come from the study config. Produces one row per
# parameter combination with a list-column of per-run results, which is the shape
# the metrics script unnests.

library(here)
source(here("sample_size/correlated/binary/bin_corr_config.R"))

workers <- 2

all_results_df <- get_results(study, workers = workers)

saveRDS(all_results_df, file.path(study$res_path, "bin_corr_all.RDS"))
print("Collection complete!")
