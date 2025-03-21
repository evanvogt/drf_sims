###############
# title: performance metrics: CATE correlations, MSE, and coverage
# date started: 20/02/25
# date finished:
# author: Ellie Van Vogt
###############

rm(list = ls(all = TRUE))

library(dplyr)
library(tidyr)

# paths
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# lists and categories
samplesizes <- c(250, 500, 1000, 5000)
models <-  c("CF", "DR_RF", "T_RF")#, "Logistic", "H_lasso")
scens <- paste0("scenario_", seq_along(1:10))


# collated datasets
data_all <- readRDS("live/data/all_data.RDS")
cate_preds <- readRDS("live/results/all/cate_preds.RDS")

cate_corr <- setNames(vector("list", length(scens)), scens)
cate_mse <- setNames(vector("list", length(scens)), scens)
cate_bias <- setNames(vector("list", length(scens)), scens)
cate_coverage <- setNames(vector("list", length(scens)), scens)

for (scenario in scens) {
  cate_corr[[scenario]] <- setNames(vector("list", length(samplesizes)), as.character(samplesizes))
  cate_mse[[scenario]] <- setNames(vector("list", length(samplesizes)), as.character(samplesizes))
  cate_coverage[[scenario]] <- setNames(vector("list", length(samplesizes)), as.character(samplesizes))
  
  for (n in samplesizes) {
    cate_corr[[scenario]][[as.character(n)]] <- setNames(vector("list", length(models)), models)
    cate_mse[[scenario]][[as.character(n)]] <- setNames(vector("list", length(models)), models)
    cate_coverage[[scenario]][[as.character(n)]] <- setNames(vector("list", length(models)), models)
    
    datasets <- data_all[[scenario]][[as.character(n)]]
    
    # list of all the true taus
    truth <- lapply(datasets, `[[`, 2)
    
    for (model in models) {
      cates <- cate_preds[[scenario]][[as.character(n)]][[model]]

      #check this model does generate taus that we can check out
      if (!is.null(cates) & length(cates)==1) {
        sim_cate <- sim_results[[scenario]][[as.character(n)]][[model]][[cates]]
        
        corrs <- rep(NA, 1000)
        mse <- rep(NA, 1000)
        coverage <- rep(NA, 1000)
        for (i in 1:1000) {
          corrs[i] <- corr(truth[[i]][["tau"]], sim_cate[[i]][["tau"]])
          mse[i] <- ((truth[[i]][["tau"]] - sim_cate[[i]][["tau"]])^2)/n
          coverage[i] <- sum(data.table::between(truth[[i]][["tau"]], sim_cate[[i]][["lb"]], sim_cate[[i]][["ub"]]))/1000
        }
        cate_corr[[scenario]][[as.character(n)]][[model]] <- corrs
        cate_mse[[scenario]][[as.character(n)]][[model]] <- mse
        cate_coverage[[scenario]][[as.character(n)]][[model]] <- coverage
      }
    }
  }
}

saveRDS(cate_corr, "live/results/all/cate_correlations.RDS")
saveRDS(cate_mse, "live/results/all/cate_mses.RDS")
saveRDS(cate_coverage, "live/results/all/cate_coverage.RDS")
