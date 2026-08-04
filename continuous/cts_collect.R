##########
# title: collect up results - continuous outcome
##########
# The grid and the results path come from the study config. Produces one row per
# parameter combination with a list-column of per-run results, which is the shape
# the metrics script unnests.

library(here)
source(here("continuous/cts_config.R"))

workers <- 2

all_results_df <- get_results(study, workers = workers)

saveRDS(all_results_df, file.path(study$res_path, "cts_all.RDS"))
print("Collection complete!")
