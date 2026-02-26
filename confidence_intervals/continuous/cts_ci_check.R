##########
# title: check for failed simulations - cts CIs
##########

# libraries
library(dplyr)
library(here)
library(purrr)

# path
res_path <- file.path(dirname(here()), "results", "confidence_intervals", "continuous")
failed_file <- here("confidence_intervals", "continuous" ,"jobscripts", "failed_ids.txt")

# parameters
n_sims <- 100

params <- expand.grid(scenario = c(1:10),
                      n = c(500, 1000),
                      CI_sf = seq(0.05, 0.5, 0.05),
                      stringsAsFactors = F)

check_failed <- function(scenario, n, CI_sf) {
  folder <- file.path(res_path, paste0("scenario_", scenario), n, CI_sf)
  result_files <- list.files(folder, pattern = "^res_sim_\\d+\\.RDS$", full.names = TRUE)
  if (length(result_files) < n_sims) {
    complete_runs <- as.numeric(gsub(".*res_sim_(\\d+)\\.RDS$", "\\1", result_files))
    failed_runs <- setdiff(1:n_sims, complete_runs)
    return(data.frame(
      scenario = scenario,
      n = n,
      CI_sf = CI_sf,
      run = failed_runs
    ))
  }
  return(NULL)
}

failed <- pmap(params, check_failed) %>% bind_rows()


if (nrow(failed) > 0) {
  full_params <- expand.grid(scenario = c(1:10),
                             n = c(500, 1000),
                             CI_sf = seq(0.05, 0.5, 0.05),
                             run = c(1:100),
                             stringsAsFactors = F)
  
  failed_idx <- which(interaction(full_params) %in% interaction(failed))
  
  cat(failed_idx, file = failed_file, sep = "\n")
  print(paste0("failed runs found (", nrow(failed), ") saved to jobscripts folder"))
} else {
  print("All simulations complete! Go ahead and collect up the results")
}