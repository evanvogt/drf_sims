##########
# title: check for failed simulations - cts + missing
##########

# libraries
library(dplyr)
library(here)
library(purrr)

# path
res_path <- file.path(dirname(here()), "results", "missing", "continuous")
failed_file <- here("missing", "continuous" ,"jobscripts", "failed_ids.txt")

# parameters
n_sims <- 100

params <- expand.grid(scenario = c(1:5),
                      n = 500,
                      type = c("prognostic", "predictive", "both"),
                      prop = 0.3,
                      mechanism = c("MAR", "AUX", "AUX-Y"),
                      method = c("complete_cases", "mean_imputation", "missforest", "regression", 
                                 "missing_indicator", "IPW", "multiple_imputation", "none"),
                      stringsAsFactors = F)
params <- params %>%
  filter(!(scenario == 1 & (type != "prognostic" | mechanism == "AUX-Y")))

check_failed <- function(scenario, n, type, prop, mechanism, method) {
  folder <- file.path(res_path, paste0("scenario_", scenario), n, type, prop, mechanism, method)
  result_files <- list.files(folder, pattern = "^res_sim_\\d+\\.RDS$", full.names = TRUE)
  if (length(result_files) < n_sims) {
    complete_runs <- as.numeric(gsub(".*res_sim_(\\d+)\\.RDS$", "\\1", result_files))
    failed_runs <- setdiff(1:n_sims, complete_runs)
    return(data.frame(
      scenario = scenario,
      n = n,
      type = type,
      prop = prop,
      mechanism = mechanism,
      method = method,
      run = failed_runs
    ))
  }
  return(NULL)
}

failed <- pmap(params, check_failed) %>% bind_rows()


if (nrow(failed) > 0) {
  full_params <- expand.grid(scenario = c(1:5),
                             n = 500,
                             type = c("prognostic", "predictive", "both"),
                             prop = 0.3,
                             mechanism = c("MAR", "AUX", "AUX-Y"),
                             method = c("complete_cases", "mean_imputation", "missforest", "regression", 
                                        "missing_indicator", "IPW", "multiple_imputation", "none"),
                             run = c(1:100),
                             stringsAsFactors = F)
  full_params <- full_params %>%
    filter(!(scenario == 1 & (type != "prognostic" | mechanism == "AUX-Y")))
  
  failed_idx <- which(interaction(full_params) %in% interaction(failed))
  
  cat(failed_idx, file = failed_file, sep = "\n")
  print(paste0("failed runs found (", nrow(failed), ") saved to jobscripts folder"))
} else {
  print("All simulations complete! Go ahead and collect up the results")
}
