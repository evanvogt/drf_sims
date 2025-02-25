###############
# title: creating a file structure for sims
# date started: 07/01/2025
# date finished:
# author: Ellie Van Vogt
###############

rm(list = ls(all = TRUE))

# paths
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# lists and categories
samplesizes <- c(250, 500, 1000, 5000)
models <- c("CF", "DR_RF", "DR_oracle", "T_RF", "Logistic", "H_lasso")


# currently has 10 scenarios - this might change
for (i in 1:10){
  scenario <- paste0(c("scenario", i), collapse = "_")
  dir.create(paste0(c("live", "results", scenario), collapse = "/"))

  # folder for each sample size within each scenario
  for (s in samplesizes){
    dir.create(paste0(c("live", "results", scenario, s), collapse = "/"))

    # folders for each model type
    for (m in models){
      dir.create(paste0(c("live", "results", scenario, s, m), collapse = "/"))
    }
  }
}
