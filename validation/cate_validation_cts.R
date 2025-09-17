###############
# validating CATEs within a trial - continuous
# 
###############

library(dplyr)
library(furrr)
library(grf)
library(rpart)


path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# functions
source("/rds/general/project/nihr_drf_simulations/live/scripts/drf_sims/validation/cts_dgm_validation.R")
source("/rds/general/project/nihr_drf_simulations/live/scripts/drf_sims/validation/CATE_validation_functions.R")
source("/rds/general/project/nihr_drf_simulations/live/scripts/drf_sims/utils.R")

# arguments to get scenario and simulation number
args <- commandArgs(trailingOnly = T)
scenario <- 3
n <- 1000
sim <- as.numeric(args[1])
interim_prop <- as.numeric(args[2])
workers <- 5

# Set up simulation seed
setup_rng_stream(sim)

# Generate data from trials before and after interim analysis
gen1 <- generate_continuous_scenario_data(scenario, n*interim_prop)
gen2 <- generate_continuous_scenario_data(scenario, n*(1-interim_prop))

# Analyse first chunk
data1 <- gen1$dataset

# All methods on first chunk
n_folds <- ifelse(n*interim_prop < 250, 5, 10)

results1 <- run_all_cate_methods(
  data = data1, 
  n_folds = n_folds, 
  workers = workers
)

results1$data <- data1
results1$truth <- gen1$truth

# analyse second chunk - naiive
data2 <- gen2$dataset

# All methods on first chunk
n_folds <- ifelse(n*(1-interim_prop) < 250, 5, 10)

results2 <- run_all_cate_methods(
  data = data2, 
  n_folds = n_folds, 
  workers = workers
)

results2$data <- data2
results2$truth <- gen2$truth

##########
# subgroups based on top and bottom responders
##########

X2 <- data2[,-c(1,2)]

models <- setdiff(names(results1), c("data", "truth"))
subgroups <- list()
for (model in models) {
  fit <- results1[[model]]
  tau1 <- fit$tau
  X1 <- data1[,-c(1,2)]
  
  group <- cut(tau1,
               breaks = quantile(tau1, probs = c(0, 0.1, 0.9, 1)),
               labels = c("bottom10", "middle", "top10"),
               include.lowest = TRUE)
  df_train <- data.frame(group = group, X1)
  
  tree_group <- rpart(group ~ ., data = df_train, method = "class")
  group_pred <- predict(tree_group, newdata = X2, type = "class")
  
  data2[paste0(model, "_top10")] <- as.numeric(group_pred == "top10")
  data2[paste0(model, "_bottom10")] <- as.numeric(group_pred == "bottom10")
  
  lm_top <- lm(Y ~ W * get(paste0(model, "_top10")), data = data2)
  lm_bottom <- lm(Y ~ W * get(paste0(model, "_bottom10")), data = data2)
  
  sum_top <- summary(lm_top)
  sum_bottom <- summary(lm_bottom)
  
  pvals_top <- sum_top$coefficients[, 4]
  pvals_bottom <- sum_bottom$coefficients[, 4]
  subgroups[[model]] <- c(top = unname(pvals_top[4]), bottom = unname(pvals_bottom[1]))
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
var_imps <- list()
for (model in models) {
  fit1 <- results1[[model]]
  tevim1 <- unlist(fit1$te_vims[1,])
  
  varnames <- colnames(fit1$te_vims[1,])
  
  fit2 <- results2[[model]]
  tevim2 <- unlist(fit2$te_vims[1,])
  
  vi_df <- data.frame(variables = varnames,
                      vi1 = rank(tevim1),
                      vi2 = rank(tevim2))
  vi_df <- vi_df %>%
    mutate(
      diff = vi2 - vi1
    )
  var_imps[[model]] <- vi_df
}

##########
# Compare HTE tests between chunks
##########

# will add later


validations <- list(subgroups = subgroups, variances = variances, var_imps = var_imps)

results <- list(results1 = results1, results2 = results2, validations = validations)

output_dir <- paste0("live/results/validation/scenario_", scenario, "/", n, "/", interim_prop, "/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, paste0(output_dir, "res_sim_", sim, ".RDS"))

print(paste0("All methods for scenario ", scenario, "_", n, " sim ", sim, " completed successfully!"))
