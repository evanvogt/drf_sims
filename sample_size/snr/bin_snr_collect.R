##########
# title: collect up results - correlated binary SNR check
##########
# One row per (scenario, n, hte_mult) with a list-column of per-run results,
# the shape bin_snr_metrics.R unnests.

library(here)
source(here("sample_size", "correlated", "binary", "snr", "bin_snr_config.R"))

workers <- 2

all_results_df <- get_results(study, workers = workers)

saveRDS(all_results_df, file.path(study$res_path, "bin_snr_all.RDS"))
print("Collection complete!")
