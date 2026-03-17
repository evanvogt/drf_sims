##########
# title: check for failed simulations - bin + missing
##########

# libraries
library(dplyr)
library(here)
library(purrr)

# path
res_path <- file.path(dirname(here()), "results", "missing", "binary")
failed_file <- here("missing", "binary" ,"jobscripts", "failed_ids.txt")

# parameters
n_sims <- 100

params <- expand.grid(
  scenario = c(1, 2, 4, 5),
  n = c(500),
  type = c("both"),
  prop = c(0.3),
  mechanism = c("MAR", "MNAR", "MNAR-Y"), # only for 1 scenario?
  method = c("complete_cases", "mean_imputation", "missforest", "regression", 
             "missing_indicator", "IPW", "multiple_imputation", "none"),
  stringsAsFactors = F)

# redundant scenarios
params <- params %>%
  filter(!(scenario == 1 & mechanism == "MNAR-Y"))

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
  full_params <- expand.grid(
    scenario = c(1, 2, 4, 5),
    n = c(500),
    type = c("both"),
    prop = c(0.3),
    mechanism = c("MAR", "MNAR", "MNAR-Y"), # only for 1 scenario?
    method = c("complete_cases", "mean_imputation", "missforest", "regression", 
               "missing_indicator", "IPW", "multiple_imputation", "none"),
    run = c(1:100),
    stringsAsFactors = F)
  
  # redundant scenarios
  full_params <- full_params %>%
    filter(!(scenario == 1 & mechanism == "MNAR-Y"))  
  failed_idx <- which(interaction(full_params) %in% interaction(failed))
  
  cat(failed_idx, file = failed_file, sep = "\n")
  print(paste0("failed runs found (", nrow(failed), ") saved to jobscripts folder"))
} else {
  print("All simulations complete! Go ahead and collect up the results")
}
