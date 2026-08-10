##########
# title: collect up results - model evaluation study
##########
# The grid and the results path come from the study config. Produces one row
# per parameter combination with a list-column of per-run results, which is
# the shape me_metrics.R unnests.
#
# Per-run objects here are heavier than continuous's (9 model results + 2
# nuisance pipelines x 2 fold regimes of prediction data - all confirmed to
# be plain numeric vectors/data.frames, not raw fitted-model objects, so
# nothing here is unsafe to have saveRDS()'d). If this job runs low on
# memory at mem=20gb, that's the first place to look - see README.md.

library(here)
source(here("model_evaluation/me_config.R"))

workers <- 2

all_results_df <- get_results(study, workers = workers)

saveRDS(all_results_df, file.path(study$res_path, "me_all.RDS"))
print("Collection complete!")
