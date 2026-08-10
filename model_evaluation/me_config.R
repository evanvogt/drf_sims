##########
# title: model evaluation study - the one definition of its parameter grid
##########
# Sourced by me_analysis.R, me_check.R, me_collect.R and me_metrics.R.
# The array index handed to me_analysis.R is a row number of `grid`, so `grid`
# must never be filtered or reordered after construction.
#
# Scenario set is deliberately 4, not all 10: this study's research question
# (do cheap proxy losses rank 9 candidate models the way true PEHE would?)
# doesn't need every CATE-structure re-litigated, the same reasoning
# crossfitting/cf_analysis.R uses for its own scenario = c(1, 4, 6, 9).

library(here)
source(here("R", "pipeline.R"))

study <- study_config(
  name     = "model_evaluation",
  res_path = file.path(dirname(here()), "results", "model_evaluation"),
  grid = expand.grid(
    scenario = c(1, 4, 6, 9),
    n = c(250, 500, 1000),
    run = c(1:30),
    stringsAsFactors = FALSE
  ),
  path_cols   = c("scenario", "n"),
  n_sims      = 30,
  failed_file = here("model_evaluation", "jobscripts", "failed_ids.txt")
)

# The 9 candidate CATE-learner configurations, as returned by
# run_all_candidate_models() (me_models.R) and saved as top-level keys of
# each run's results (me_analysis.R). Declared once here so me_models.R and
# me_metrics.R cannot disagree about the model list - the same reason the
# grid above is declared once rather than re-typed in every script.
CANDIDATE_MODELS <- c("rf1", "rf2", "rf3", "net1", "net2", "net3", "SL1", "SL2", "SL3")
