######################
# title: hierarchical lasso
# date started: 11/02/25
# date finished:
# author: Ellie Van Vogt
#####################
set.seed(1998)
# wondering if I should apply a crossfitting approach to all of these too??

#libraries
library(stats)
library(dplyr)
library(furrr)
library(glinternet)

# paths 
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# arguments and parameters
args <- commandArgs(trailingOnly = TRUE)
scenario <- as.character(args[1])
n <- as.numeric(args[2])

n_cores <- 10 #floor(future::availableCores() *0.9)
n_folds <- 10

oldplan <- plan(multicore, workers = n_cores)

# load in the data
datasets <- readRDS(paste0(c("live/data/continuous/", scenario, "_", n, ".RDS"), collapse = ""))
datasets <- lapply(datasets, `[[`, 1) # just want the data not the truth



# functions
get_levels <- function(x) {
  ifelse(max(x) - min(x) == 1, 2, 1)
}
#levels <- apply(X, 2, get_levels) %>% unname()

#main
hlasso <- function(data) {
  X <- as.matrix(data[,-1])  # treatment and covariates for interactions
  Y <- data$Y

  model.cv <- glinternet.cv(X, Y, numLevels = rep(1, ncol(X)), nFolds = n_folds, nLambda = 100, interactionCandidates = 1, family = "gaussian")
  
  lam <- model.cv$lambdaHat1Std
  ind <- which(model.cv$lambda == lam)

  coefs <- coef(model.cv, lambdaIndex = ind)
  
  # extract main effect terms:
  if (length(coefs[["mainEffects"]]$cont) == 0) {
    #print("model did not keep any main effects")
    main <- NULL
  } else {
    main_names <- unique(colnames(X)[coefs[["mainEffects"]]$cont])
    
    main <- unlist(coefs[["mainEffectsCoef"]]$cont)
    names(main) <- main_names
  }
  # extract the interactions terms
  if (length(coefs[["interactions"]]$contcont) == 0) {
    #print("model did not keep any interactions")
    interactions <- NULL
  } else {
    int_names <- unique(colnames(X)[coefs[["interactions"]]$contcont])
    int_names <- int_names[int_names != "W"]
    
    interactions <- unlist(coefs[["interactionsCoef"]]$contcont)
    names(interactions) <- int_names
  }
  list(main = main, interactions = interactions)
}

# Parallelise the function
t0 <- Sys.time()
results <- future_map(datasets, hlasso, .progress = T, .options = furrr_options(seed = T))
t1 <- Sys.time()
print(t1-t0)
plan(oldplan)

# Save results
saveRDS(results, paste0(c("live/results/continuous/", scenario, "/", n, "/H_lasso/", "HL_interactions.RDS"), collapse = ""))
