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
i <- as.numeric(commandArgs(trailingOnly = TRUE))

workers <- 5

# The parameter grid lives in the study config, so this script and the
# check/collect scripts cannot disagree about what index i means.
param <- study$grid[i, ]
print(param)

scenario <- param$scenario
n <- param$n
interim_prop <- param$interim_prop
run <- param$run

# set up simulation seed
setup_rng_stream(run)

# Generate data from trial chunks before and after the interim analysis
gen1 <- generate_continuous_scenario_data(scenario, n * interim_prop)
gen2 <- generate_continuous_scenario_data(scenario, n * (1 - interim_prop))

data1 <- gen1$dataset
data2 <- gen2$dataset

n_folds1 <- ifelse(n * interim_prop < 250, 5, 10)
n_folds2 <- ifelse(n * (1 - interim_prop) < 250, 5, 10)

metaplan <- plan(multisession, workers = workers)
on.exit(plan(metaplan), add = TRUE)

# Fit both estimators on each chunk
results1 <- run_all_cate_methods(data = data1, n_folds = n_folds1)
results1$data <- data1
results1$truth <- gen1$truth

results2 <- run_all_cate_methods(data = data2, n_folds = n_folds2)
results2$data <- data2
results2$truth <- gen2$truth

##########
# subgroups based on top and bottom responders
##########

X2 <- data2[, -c(1, 2)]

models <- setdiff(names(results1), c("data", "truth"))
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

output_dir <- file.path(study$res_path, paste0("scenario_", scenario), n, interim_prop)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, paste0("res_sim_", run, ".RDS")))

print(paste0("All methods for scenario ", scenario, "_", n, " interim ", interim_prop,
            " run ", run, " completed successfully!"))
