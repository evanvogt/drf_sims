##########
# Title: collect up the results - missing bin outcome
##########

# libraries
library(here)
library(dplyr)
library(future)
library(furrr)

# path
res_path <- file.path(dirname(here()), "results", "missing", "binary")
out_file <- file.path(res_path, "bin_miss_all.RDS")

# parameters
n_sims <- 100
workers <- 2

params <- expand.grid(
  scenario = c(1, 2, 4, 5),
  n = c(500),
  type = c("both"),
  prop = c(0.3),
  mechanism = c("MAR", "MNAR", "MNAR-Y"), # only for 1 scenario?
  method = c("complete_cases", "mean_imputation", "missforest", "regression", 
             "missing_indicator", "IPW", "multiple_imputation", "none"),
  stringsAsFactors = F)

params <- params %>%
  filter(!(scenario == 1 & mechanism == "MNAR-Y"))

get_results <- function(scenario, n, type, prop, mechanism, method) {
  folder <- file.path(res_path, paste0("scenario_", scenario), n, type, prop, mechanism, method)
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
    type = type,
    prop = prop,
    mechanism = mechanism,
    method = method,
    results = list(temp)
  )
}

plan(multisession, workers = workers)
all_results <- future_pmap(params, get_results)
plan(sequential)

# remove nulls
all_results_df <- bind_rows(all_results[!sapply(all_results, is.null)])

# save results
saveRDS(all_results_df, out_file)
print("Collection complete")