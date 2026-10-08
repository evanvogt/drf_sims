##########
# title: interim-analysis validation - binary outcome
##########
# Fits the three estimators (causal forest, DR random forest, DR SuperLearner)
# on a trial split into two chronological chunks (before and after an "interim
# analysis" at `interim_prop` of the way through), then checks whether
# subgroups, CATE variance and variable-importance ranking found on the first
# chunk hold up on the second. continuous/cts_val_analysis.R with a binary
# outcome: the splitting, fitting and comparing are shared - see
# validation/val_common.R.

library(dplyr)
library(furrr)
library(grf)
library(SuperLearner)
library(rpart)
library(here)

# functions
source(here("R", "utils.R"))
source(here("validation", "binary", "bin_val_dgms.R"))
source(here("validation", "binary", "bin_val_models.R"))
source(here("validation", "binary", "bin_val_config.R"))

# simulation parameters
#
# workers and grf_threads come off the jobscript's Rscript line, so the PBS
# resource request and the R-level parallelism cannot drift apart - as
# continuous/. A bare `Rscript bin_val_analysis.R <i>` still works as a local
# smoke test (5 workers, grf unthrottled).
args <- commandArgs(trailingOnly = TRUE)

i <- as.numeric(args[1])
workers <- if (length(args) >= 2) as.integer(args[2]) else 5L
grf_threads <- if (length(args) >= 3) as.integer(args[3]) else NULL

# The parameter grid lives in the study config, so this script and the
# check/collect scripts cannot disagree about what index i means.
param <- study$grid[i, ]
print(param)

scenario <- param$scenario
n <- param$n
interim_prop <- param$interim_prop
run <- param$run
rho <- param$rho

# set up simulation seed - the run alone, so every interim_prop of a run splits
# the same trial
setup_rng_stream(run)

# One trial of n, split at the interim analysis - see split_trial()
gen <- generate_binary_scenario_data(scenario, n, rho)
chunks <- split_trial(gen, interim_prop)
data1 <- chunks$data1
data2 <- chunks$data2

# Folds (chunk_folds()) and SuperLearner libraries for the DR SuperLearner,
# each sized to its own chunk
n_folds1 <- chunk_folds(nrow(data1))
n_folds2 <- chunk_folds(nrow(data2))
sl_lib1 <- sl_libraries(nrow(data1))
sl_lib2 <- sl_libraries(nrow(data2))

# multisession workers are new R processes and inherit this, so it keeps their
# OpenMP thread pools in step with grf's num.threads - as continuous/
if (!is.null(grf_threads)) Sys.setenv(OMP_NUM_THREADS = grf_threads)

metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

# Fit the three estimators on each chunk
results1 <- run_all_cate_methods(data = data1, n_folds = n_folds1,
                                 num.threads = grf_threads, sl_lib = sl_lib1)
results1$data <- data1
results1$truth <- chunks$truth1

results2 <- run_all_cate_methods(data = data2, n_folds = n_folds2,
                                 num.threads = grf_threads, sl_lib = sl_lib2)
results2$data <- data2
results2$truth <- chunks$truth2

# The four chunk comparisons: subgroups, variances, var_imps, top_var_tests.
# robust = TRUE: with a binary Y every interaction test is a linear probability
# model - linear in the risk, the scale the treatment effect is on - whose
# errors are heteroskedastic by construction, so its standard errors are HC3
# (coef_pval(), validation/val_common.R). The continuous arm uses HC3 too.
validations <- chunk_validations(results1, results2, data1, data2, robust = TRUE)

results <- list(results1 = results1, results2 = results2, validations = validations)

# under the path check/collect expect (combo_dir())
output_dir <- combo_dir(study, param)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print(paste0("All methods for rho ", rho, " scenario ", scenario, "_", n,
             " interim ", interim_prop, " run ", run, " completed successfully!"))
