##########
# title: collect up results - competing-risk SNR check
##########
# One row per (n, hte_mult) with a list-column of per-run results, the shape
# snr_metrics.R unnests.

library(here)
source(here("competing_risk", "snr", "snr_config.R"))

workers <- 2

all_results_df <- get_results(study, workers = workers)

saveRDS(all_results_df, file.path(study$res_path, "snr_all.RDS"))
print("Collection complete!")
