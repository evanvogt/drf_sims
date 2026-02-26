##########
# title: collect all binary outcome confidence interval results
##########

# libraries
library(dplyr)
library(here)
library(future)
library(furrr)

# path
res_path <- file.path(dirname(here()), "results", "confidence_intervals", "binary")
out_file <- file.path(res_path, "ci_bin_all.RDS")

# parameters
n_sims <- 100
workers <- 4

params <- expand.grid(scenario = c(1:10),
                      n = c(500, 1000),
                      CI_sf = seq(0.05, 0.5, 0.05),
                      stringsAsFactors = F)

get_results <- function(scenario, n, CI_sf) {
  folder <- file.path(res_path, paste0("scenario_", scenario), n, CI_sf)
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
    CI_sf = CI_sf,
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