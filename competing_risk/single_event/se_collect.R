##########
# title: collect up results - single-event survival
##########
# One row per parameter combination with a list-column of per-run results, the
# shape se_metrics.R unnests.

library(here)
source(here("competing_risk", "single_event", "se_config.R"))

workers <- 2

all_results_df <- get_results(study, workers = workers)

saveRDS(all_results_df, file.path(study$res_path, "se_all.RDS"))
print("Collection complete!")
