##########
# title: check for failed simulations - cts + missing CI
##########

# libraries
library(dplyr)
library(here)
library(purrr)

# path
res_path <- file.path(dirname(here()), "results", "missing", "ci_example")
failed_file <- here("missing", "ci_example" ,"jobscripts", "failed_ids.txt")

# parameters
n_sims <- 100

params <- expand.grid(
  scenario = c(1:5),
  n = c(500),
  type = c("both"),
  prop = c(0.3),
  mechanism = c("MAR"),
  method = c("multiple_imputation"),
  stringsAsFactors = F
)

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
    scenario = c(1:5),
    n = c(500),
    type = c("both"),
    prop = c(0.3),
    mechanism = c("MAR"),
    method = c("multiple_imputation"),
    run = c(1:100),
    stringsAsFactors = F
  )
  full_params <- full_params %>%
    arrange(scenario, run)
  
  failed_idx <- which(interaction(full_params) %in% interaction(failed))
  
  cat(failed_idx, file = failed_file, sep = "\n")
  print(paste0("failed runs found (", nrow(failed), ") saved to jobscripts folder"))
} else {
  print("All simulations complete! Go ahead and collect up the results")
}
