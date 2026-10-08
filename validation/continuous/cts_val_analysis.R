##########
# title: interim-analysis validation - continuous outcome
##########
# Fits both estimators on a trial split into two chronological chunks (before
# and after an "interim analysis" at `interim_prop` of the way through), then
# checks whether subgroups, CATE variance and variable-importance ranking
# found on the first chunk hold up on the second.

library(dplyr)
library(furrr)
library(grf)
library(rpart)
library(here)

# functions
source(here("R", "utils.R"))
source(here("validation", "continuous", "cts_val_dgms.R"))
source(here("validation", "continuous", "cts_val_models.R"))
source(here("validation/continuous/cts_val_config.R"))

# simulation parameters
#
# workers and grf_threads come off the jobscript's Rscript line rather than
# being hardcoded here, so the PBS resource request and the R-level parallelism
# cannot drift apart - same arrangement as sample_size/correlated/continuous/cts_corr_analysis.R. The
# defaults reproduce what this script did before they were arguments, so a bare
# `Rscript cts_val_analysis.R <i>` (cts_val_testing.R check 7) still works.
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

# One trial of n, split at the interim analysis: the first n * interim_prop
# participants are chunk 1, the rest chunk 2. Rows are iid, so the first n1 are
# a valid "enrolled by the interim" cohort. This used to draw the two chunks as
# two separate datasets, and generate_scenario_data() calibrates bW to 80% power
# at the n it is given - so each chunk was its own trial with its own ATE, not
# two halves of one. Splitting one draw also nests the chunks across
# interim_prop: run r's chunk 1 at 0.25 is a subset of its chunk 1 at 0.30, so
# curves over interim_prop are paired within run.
#
# round() because n * interim_prop is not always an exact integer in floating
# point (1000 * 0.35 need not be 350), and an index sequence would truncate it.
gen <- generate_continuous_scenario_data(scenario, n, rho)

n1 <- round(n * interim_prop)
chunk1 <- seq_len(n1)
chunk2 <- (n1 + 1):n

data1 <- gen$dataset[chunk1, ]
data2 <- gen$dataset[chunk2, ]
rownames(data2) <- NULL

truth1 <- gen$truth[chunk1, , drop = FALSE]
truth2 <- gen$truth[chunk2, , drop = FALSE]
rownames(truth2) <- NULL

n_folds1 <- ifelse(nrow(data1) < 250, 5, 10)
n_folds2 <- ifelse(nrow(data2) < 250, 5, 10)

# multisession workers are new R processes and inherit this, so setting it here
# does control their OpenMP thread pools even though this process's own libraries
# have already initialised - matches sample_size/correlated/continuous/cts_corr_analysis.R
if (!is.null(grf_threads)) Sys.setenv(OMP_NUM_THREADS = grf_threads)

metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

# Fit both estimators on each chunk
results1 <- run_all_cate_methods(data = data1, n_folds = n_folds1,
                                 num.threads = grf_threads)
results1$data <- data1
results1$truth <- truth1

results2 <- run_all_cate_methods(data = data2, n_folds = n_folds2,
                                 num.threads = grf_threads)
results2$data <- data2
results2$truth <- truth2

##########
# subgroups based on top and bottom responders
##########

X2 <- data2[, -c(1, 2)]

# "timings" is only present when run_all_cate_methods() was called with
# verbose_timing = TRUE, which this script does not do, but it is excluded here
# anyway so every non-data element of `results1` is guaranteed to be a model.
models <- setdiff(names(results1), c("data", "truth", "timings"))
subgroups <- list()
for (model in models) {
  fit <- results1[[model]]
  tau1 <- fit$tau
  X1 <- data1[, -c(1, 2)]

  group <- cut(tau1,
               breaks = quantile(tau1, probs = c(0, 0.1, 0.9, 1)),
               labels = c("bottom10", "middle", "top10"),
               include.lowest = TRUE)
  df_train <- data.frame(group = group, X1)

  tree_group <- rpart(group ~ ., data = df_train, method = "class")
  group_pred <- predict(tree_group, newdata = X2, type = "class")

  data2[paste0(model, "_top10")] <- as.numeric(group_pred == "top10")
  data2[paste0(model, "_bottom10")] <- as.numeric(group_pred == "bottom10")

  # interaction_pval (cts_val_models.R) indexes the W:v coefficient by name.
  # This used to be two positional lookups, and the bottom one read
  # `pvals_bottom[1]` - the intercept - so bottom_pval was never a subgroup
  # test at all. Results from before this fix are not comparable.
  subgroups[[model]] <- c(
    top    = interaction_pval(data2$Y, data2$W, data2[[paste0(model, "_top10")]]),
    bottom = interaction_pval(data2$Y, data2$W, data2[[paste0(model, "_bottom10")]])
  )
}

##########
# Compare variance between early and later chunks
##########
variances <- list()
for (model in models) {
  fit1 <- results1[[model]]
  tau1 <- fit1$tau

  fit2 <- results2[[model]]
  tau2 <- fit2$tau

  vt1 <- var(tau1)
  vt2 <- var(tau2)

  variances[[model]] <- c(vt1 = unname(vt1), vt2 = unname(vt2))
}

##########
# Compare variable importance between early and late chunks
##########
# Two measures now: the TE-VIMs and surrogate TreeSHAP (both in
# cts_val_models.R). Both are larger-is-more-important, so rank() means the same
# thing for each - rank 1 is the least important covariate.
measure_fields <- c(tevim = "te_vims", shap = "shap_vims")

var_imps <- list()
for (model in models) {
  fit1 <- results1[[model]]
  fit2 <- results2[[model]]

  per_measure <- lapply(names(measure_fields), function(measure) {
    field <- measure_fields[[measure]]
    imp1 <- unlist(fit1[[field]][1, ])
    imp2 <- unlist(fit2[[field]][1, ])

    data.frame(variables = colnames(fit1[[field]]),
               measure = measure,
               vi1 = rank(imp1),
               vi2 = rank(imp2),
               stringsAsFactors = FALSE) %>%
      mutate(diff = vi2 - vi1)
  })

  var_imps[[model]] <- do.call(rbind, per_measure)
}

##########
# Carry the top-ranked covariate into the remaining participants
##########
# The point of ranking covariates is whether the winner means anything, so take
# each measure's chunk-1 most important covariate and interaction-test it in
# chunk 2. Two forms side by side: the continuous W x X_top interaction (works
# for X1/X2/X4 and the already-binary X01-X05 alike, no arbitrary cut point),
# and a median split, which is directly parallel to the top10/bottom10 tests
# above. x_top2 is chunk 2's own winner, kept so the report can ask how often
# the two chunks even agree on which covariate matters most.
#
# p_cts_adj is the continuous test again, adjusted for every other covariate
# (interaction_pval_adj, cts_val_models.R). The covariates are correlated
# (rho = 0.5), so a non-modifier correlated with X4 shows a real *marginal*
# W x X interaction in p_cts; p_cts_adj asks whether x_top modifies the effect
# given the rest.
top_var_tests <- list()
for (model in models) {
  vi <- var_imps[[model]]

  rows <- lapply(split(vi, vi$measure), function(v) {
    x_top <- v$variables[which.max(v$vi1)]
    xt <- data2[[x_top]]

    data.frame(measure = v$measure[1],
               x_top = x_top,
               x_top2 = v$variables[which.max(v$vi2)],
               p_cts = interaction_pval(data2$Y, data2$W, xt),
               p_cts_adj = interaction_pval_adj(data2$Y, data2$W, X2, x_top),
               p_split = interaction_pval(data2$Y, data2$W,
                                          as.numeric(xt > median(xt))),
               stringsAsFactors = FALSE)
  })

  top_var_tests[[model]] <- do.call(rbind, rows)
}

##########
# Compare HTE tests between chunks
##########

# TODO: both estimators now carry BLP_whole/independence_cate/independence_po
# (see R/cate_models.R, R/metrics.R::hte_test_metrics()) in a shape a chunk-vs-
# chunk comparison could use directly. Not implemented yet - see
# validation/README.md's Status section.

validations <- list(subgroups = subgroups, variances = variances,
                    var_imps = var_imps, top_var_tests = top_var_tests)

results <- list(results1 = results1, results2 = results2, validations = validations)

# under the path check/collect expect (combo_dir()) - built by hand here it would
# not follow path_cols, which now lead with rho
output_dir <- combo_dir(study, param)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print(paste0("All methods for rho ", rho, " scenario ", scenario, "_", n,
             " interim ", interim_prop, " run ", run, " completed successfully!"))
