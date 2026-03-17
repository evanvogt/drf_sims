##########
# title: collect up binary outcome results
##########
#libraries
library(dplyr)
library(tidyr)
library(here)
library(future)
library(furrr)

# paths
res_path <- file.path(dirname(here()), "results", "binary")
out_file <- file.path(res_path, "bin_all.RDS")

# parameters
n_sims <- 100
workers <- 2

params <- expand.grid(scenario = c(1, 3, 8, 9),
                      n = c(100, 250, 500, 1000),
                      stringsAsFactors = F)

get_results <- function(scenario, n) {
  folder <- file.path(res_path, paste0("scenario_", scenario), n)
  result_files <- list.files(folder, pattern = "^res_sim_\\d+\\.RDS$", full.names = TRUE)
  
  if (length(result_files) == 0) return(NULL)
  
  temp <- map(result_files, function(res_file) {
    sim_res <- readRDS(res_file)
    sim_num <- as.integer(gsub(".*res_sim_(\\d+)\\.RDS$", "\\1", res_file))
    list(run = sim_num, result = sim_res)
  })
  gc()
  tibble(
    scenario = scenario,
    n = n,
    results = list(temp)
  )
}

plan(multisession, workers = workers)
all_results <- future_pmap(params, get_results)
plan(sequential)

# remove nulls
all_results_df <- bind_rows(all_results[!sapply(all_results, is.null)])

# Save results
saveRDS(all_results_df, out_file)
print("Collection complete!")