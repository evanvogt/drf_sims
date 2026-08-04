##########
# title: check for failed simulations - crossfitting comparison
##########

# libraries
library(dplyr)
library(here)
library(purrr)

# path
res_path <- file.path(dirname(here()), "results", "crossfitting")
failed_file <- here("crossfitting", "jobscripts", "failed_ids.txt")

# parameters - must match cf_analysis.R exactly
n_sims <- 500

full_params <- expand.grid(scenario = c(1, 4, 6, 9),
                           n = 500,
                           run = c(1:n_sims),
                           stringsAsFactors = F)

params <- distinct(full_params, scenario, n)

check_failed <- function(scenario, n) {
  folder <- file.path(res_path, paste0("scenario_", scenario))
  result_files <- list.files(folder, pattern = "^res_sim_\\d+\\.RDS$", full.names = TRUE)
  if (length(result_files) < n_sims) {
    complete_runs <- as.numeric(gsub(".*res_sim_(\\d+)\\.RDS$", "\\1", result_files))
    failed_runs <- setdiff(1:n_sims, complete_runs)
    return(data.frame(scenario = scenario, n = n, run = failed_runs))
  }
  return(NULL)
}

failed <- pmap(params, check_failed) %>% bind_rows()

if (nrow(failed) > 0) {
  failed_idx <- which(interaction(full_params) %in% interaction(failed))

  cat(failed_idx, file = failed_file, sep = "\n")
  print(paste0("failed runs found (", nrow(failed), ") saved to jobscripts folder"))
} else {
  print("All simulations complete! Go ahead and collect up the results")
}
