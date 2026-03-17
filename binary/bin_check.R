##########
# title: check for failed simulations - bin outcome
##########

# libraries
library(dplyr)
library(here)
library(purrr)

# path
res_path <- file.path(dirname(here()), "results", "binary")
failed_file <- here("binary" ,"jobscripts", "failed_ids.txt")

# parameters
n_sims <- 100

params <- expand.grid(scenario = c(1, 3, 8, 9),
                      n = c(100, 250, 500, 1000),
                      stringsAsFactors = F)

check_failed <- function(scenario, n) {
  folder <- file.path(res_path, paste0("scenario_", scenario), n)
  result_files <- list.files(folder, pattern = "^res_sim_\\d+\\.RDS$", full.names = TRUE)
  if (length(result_files) < n_sims) {
    complete_runs <- as.numeric(gsub(".*res_sim_(\\d+)\\.RDS$", "\\1", result_files))
    failed_runs <- setdiff(1:n_sims, complete_runs)
    return(data.frame(
      scenario = scenario,
      n = n,
      run = failed_runs
    ))
  }
  return(NULL)
}

failed <- pmap(params, check_failed) %>% bind_rows()


if (nrow(failed) > 0) {
  full_params <- expand.grid(scenario = c(1, 3, 8, 9),
                             n = c(100, 250, 500, 1000),
                             run = c(1:100),
                             stringsAsFactors = F)
  
  failed_idx <- which(interaction(full_params) %in% interaction(failed))
  
  cat(failed_idx, file = failed_file, sep = "\n")
  print(paste0("failed runs found (", nrow(failed), ") saved to jobscripts folder"))
} else {
  print("All simulations complete! Go ahead and collect up the results")
}
