###############
# title: generating job scripts across all scenarios and sample sizes
# date started: 07/01/2025
# date finished:
# author: Ellie Van Vogt
###############
set.seed(1998)

# libraries ----
library(dplyr)
library(syrup)
library(hms)

# paths ----
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# parameters -----

samplesizes <- c(250, 500, 1000, 5000)
models <- c("CF", "DR_RF", "T_RF")# , "Logistic", "H_lasso") # not doing the last two because they are run by scenario anyways
scens <- paste0("scenario_", seq_along(1:10))
test_scenarios <- paste0("test_scenario_", seq_along(1:10))


resource_data <- read.csv("resource_estimates.csv")


# stuff for changing the time for one run to 1000 runs across 10 cores

for (i in 1:nrow(resource_data)) {
  script <- resource_data$script[i]
  scenario <- resource_data$scenario[i]
  n <- resource_data$n[i]
  time <- resource_data$time[i]
  cpu <- resource_data$cpu
  mem <- ceiling(resource_data$mem_gb[i])

  
  job_script <- paste0(script, "_", scenario, "_", method, ".sh")
  
  writeLines(c(
    paste0("#PBS -l walltime=", time),
    paste0("#PBS -l select=1:ncpus=", cpu, ":mem=", mem),
    paste0("#PBS -N ", script, "_", scenario, "_", n),
    "",
    "module purge",
    "module add tools/prod",
    "module add R/4.2.1-foss-2022a",
    "eval \"$(~/miniforge3/bin/conda shell.bash hook)\"",
    "conda activate drf-env",
    "",
    paste0("cd /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/", model),
    "",
    paste0("scenario=", "\"", scenario, "\""),
    paste0("n=", "\"", n, "\""),
    "",
    paste0("Rscript ", script, ".R", scenario, n)
  ), con = paste0("live/scripts/", model, "/", job_script))
}
