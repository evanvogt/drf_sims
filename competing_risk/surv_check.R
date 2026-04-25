##########
# title: check for failed simulations - competing risks
##########

# libraries
library(dplyr)
library(here)
library(purrr)

# path
res_path <- file.path(dirname(here()), "results", "competing_risk")
failed_file <- here("competing_risk" ,"jobscripts", "failed_ids.txt")

# parameters
n_sims <- 100

params <- expand.grid(scenario = c(1:7),
                      n = c(500),
                      censoring = c(TRUE, FALSE),
                      stringsAsFactors = F)

check_failed <- function(scenario, n, censoring) {
  folder <- file.path(res_path, paste0("scenario_", scenario), n, paste0("censor_", censoring))
  result_files <- list.files(folder, pattern = "^res_sim_\\d+\\.RDS$", full.names = TRUE)
  if (length(result_files) < n_sims) {
    complete_runs <- as.numeric(gsub(".*res_sim_(\\d+)\\.RDS$", "\\1", result_files))
    failed_runs <- setdiff(1:n_sims, complete_runs)
    return(data.frame(
      scenario = scenario,
      n = n,
      censoring = censoring,
      run = failed_runs
    ))
  }
  return(NULL)
}

failed <- pmap(params, check_failed) %>% bind_rows()


if (nrow(failed) > 0) {
  full_params <- expand.grid(scenario = c(1:7),
                             n = c(500),
                             censoring = c(TRUE, FALSE),
                             run = c(1:100),
                             stringsAsFactors = F)
  
  failed_idx <- which(interaction(full_params) %in% interaction(failed))
  
  cat(failed_idx, file = failed_file, sep = "\n")
  print(paste0("failed runs found (", nrow(failed), ") saved to jobscripts folder"))
} else {
  print("All simulations complete! Go ahead and collect up the results")
}
