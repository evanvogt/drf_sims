###############
# title: performance metrics: BLP test metrics
# date started: 20/02/25
# date finished:
# author: Ellie Van Vogt
###############

rm(list = ls(all = TRUE))

library(dplyr)
library(tidyr)
library(purrr)

# paths
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# lists and categories
samplesizes <- c(250, 500, 1000)#, 5000)
models <- c("CF", "DR_oracle", "DR_RF", "T_RF")#, "Logistic", "H_lasso")
scens <- paste0("scenario_", seq_along(1:10))


# collated datasets
data_all <- readRDS("live/data/continuous/all_data.RDS")
BLPs_whole <- readRDS("live/results/continuous/all/BLP_wholes.RDS")

# maybe just look at the whole set sets for now?
beta_2_whole <- setNames(vector("list", length(scens)), scens)
beta_1_whole <- setNames(vector("list", length(scens)), scens)

for (scenario in scens) {
  beta_1_whole[[scenario]] <- setNames(vector("list", length(samplesizes)), as.character(samplesizes))
  beta_2_whole[[scenario]] <- setNames(vector("list", length(samplesizes)), as.character(samplesizes))
  
  for (n in samplesizes) {
    beta_1_whole[[scenario]][[as.character(n)]] <- setNames(vector("list", length(models)), models)
    beta_2_whole[[scenario]][[as.character(n)]] <- setNames(vector("list", length(models)), models)
    
    for (model in models) {
      # specify BLP whole list
      BLP_list <- BLPs_whole[[scenario]][[as.character(n)]][[model]]
      
      #check this model does generate taus that we can check out
      if (!is.null(BLP_list) & length(BLP_list)==1000) {
        beta1 <- rep(NA, 1000)
        beta2 <- rep(NA, 1000)
        
        for (i in 1:1000) {
          beta1[i] <- BLP_list[[i]][3,2]
          beta2[i] <- BLP_list[[i]][4,2]
        }
        beta_1_whole[[scenario]][[as.character(n)]][[model]] <- beta1
        beta_2_whole[[scenario]][[as.character(n)]][[model]] <- beta2
      }
    }
  }
}


# ok so now we've got the beta's we can get out the performance metrics
BLP_performance <- setNames(vector("list", length(scens)), scens)

# set a significance threshold
threshold <- 0.1
for (scenario in scens) {
  BLP_performance[[scenario]] <- setNames(vector("list", length(samplesizes)), as.character(samplesizes))

  for (n in samplesizes) {
    BLP_performance[[scenario]][[as.character(n)]] <- setNames(vector("list", length(models)), models)

    for (model in models) {
      # specify BLP whole list
      test_vals <- beta_2_whole[[scenario]][[as.character(n)]][[model]]
      
      #check this model does generate taus that we can check out
      if (!is.null(test_vals) & length(test_vals)==1000) {
        
        signif <- test_vals < threshold
        
        
        # not super sure about this bit
        if (scenario == "scenario_1") {
          truth <- rep(FALSE, 1000)
        } else {
          truth <- rep(TRUE, 1000)
        }
        
        TP <- sum(signif & truth)
        FP <- sum(signif & !truth)
        TN <- sum(!signif & !truth)
        FN <- sum(!signif & truth)
        
        metrics <- setNames(c(TP, FP, TN, FN), c("TP", "FP", "TN", "FN"))

        BLP_performance[[scenario]][[as.character(n)]][[model]] <- metrics
      }
    }
  }
}

# squish this into one data frame

performance_df <- map_dfr(names(BLP_performance), function(scenario) {
  map_dfr(names(BLP_performance[[scenario]]), function(n) {
    map_dfr(names(BLP_performance[[scenario]][[n]]), function(model) {
      entry <- BLP_performance[[scenario]][[n]][[model]]
      
      data.frame(
        scenario = scenario,
        n = as.numeric(n),
        model = model,
        TP = if (!is.null(entry)) as.numeric(entry["TP"]) else NA,
        FP = if (!is.null(entry)) as.numeric(entry["FP"]) else NA,
        TN = if (!is.null(entry)) as.numeric(entry["TN"]) else NA,
        FN = if (!is.null(entry)) as.numeric(entry["FN"]) else NA
      )
    })
  })
})

performance_df <- performance_df %>%
  mutate(
    sensitivity = TP/(TP+FN),
    specificity = TN/(TN+FP)
  )

saveRDS(performance_df, "live/results/continuous/all/BLP_tests_metrics.RDS")
