##########
# title: collect up results - interim-analysis validation study
##########
# The grid and the results path come from the study config. Produces one row
# per parameter combination with a list-column of per-run results, which is
# the shape cts_val_metrics.R unnests.

library(here)
source(here("validation/cts_val_config.R"))

workers <- 2

all_results_df <- get_results(study, workers = workers)

saveRDS(all_results_df, file.path(study$res_path, "cts_val_all.RDS"))
print("Collection complete!")
