##########
# title: collect up results - model evaluation study
##########
# The grid and the results path come from the study config. Produces one row
# per parameter combination with a list-column of per-run results, which is
# the shape me_metrics.R unnests.
#
#   Rscript me_collect.R                # the main study
#   Rscript me_collect.R strategies     # the nuisance-arm pass
#   Rscript me_collect.R split          # the 80:20 arm
#
# Per-run objects here are heavier than continuous's (9 model results + 2
# nuisance pipelines x N arms of prediction data - all confirmed to be plain
# numeric vectors/data.frames, not raw fitted-model objects, so nothing here
# is unsafe to have saveRDS()'d). The strategies tree carries FOUR arms per
# pipeline rather than two, so its per-run objects are roughly twice the size
# of the main study's - if this job runs low on memory at mem=20gb, that is
# the first place to look, and the strategies collect is the one to raise
# first.

library(here)
source(here("model_evaluation/me_config.R"))

which_study <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(which_study)) which_study <- "main"
st <- me_study(which_study)

workers <- 2

all_results_df <- get_results(st, workers = workers)

saveRDS(all_results_df, file.path(st$res_path, paste0(st$prefix, "_all.RDS")))
print(paste0("Collection complete! (", st$name, ")"))
