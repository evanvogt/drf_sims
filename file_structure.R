###############
# title: creating a file structure for sims
# date started: 07/01/2025
# date finished:
# author: Ellie Van Vogt
###############

rm(list = ls(all = TRUE))

# paths
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live"
setwd(path)

dir.create("results")
dir.create("scripts")
dir.create("data")
# lists and categories
samplesizes <- c(100, 250, 500, 1000)
models <- c("CF", "DR_RF", "DR_oracle", "DR_semi_oracle", "DR_superlearner", "T_RF", "Logistic", "H_lasso")
outcomes <- c("binary", "continuous", "competing_risk")


# currently has 10 scenarios - this might change
for (outcome in outcomes) {
  dir.create(paste0(c("results", outcome), collapse = "/"))
  dir.create(paste0(c("scripts", outcome), collapse = "/"))
  dir.create(paste0(c("data", outcome), collapse = "/"))
  
  for (i in 1:10){
    scenario <- paste0(c("scenario", i), collapse = "_")
    dir.create(paste0(c("results", outcome, scenario), collapse = "/"))
    
    # folder for each sample size within each scenario
    for (s in samplesizes){
      dir.create(paste0(c("results", outcome, scenario, s), collapse = "/"))
      
      # folders for each model type
      for (m in models){
        dir.create(paste0(c("results", outcome, scenario, s, m), collapse = "/"))
      }
    }
  }
}

# make folders for all of the job scripts?
for (o in outcomes) {
  for (m in models) {
    dir.create(paste0(c("scripts/", o, "/", m), collapse = ""))
    dir.create(paste0(c("scripts/", o, "/", m, "/jobscripts"), collapse = ""))
    for (s in samplesizes) {
      dir.create(paste0(c("scripts/", o, "/", m, "/jobscripts/logs_", s), collapse = ""))
    }
  }
}
