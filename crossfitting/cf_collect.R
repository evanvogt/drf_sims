##########
# title: collect up crossfitting comparison results
##########
# Reads the per-run files and applies the metric definitions in one streaming
# pass, so the flat metrics tibble is the only thing held in memory. The other
# studies write a large nested intermediate first; there is nothing to gain from
# that here because the per-run files carry no data or nuisance matrices.

# libraries
library(dplyr)
library(tidyr)
library(here)
library(future)
library(furrr)

# functions
source(here("crossfitting", "cf_metrics.R"))

# paths
res_path <- file.path(dirname(here()), "results", "crossfitting")
out_file <- file.path(res_path, "cf_metrics.RDS")

# parameters
workers <- 4

params <- expand.grid(scenario = c(1, 4, 6, 9),
                      stringsAsFactors = F)

get_metrics <- function(scenario) {
  folder <- file.path(res_path, paste0("scenario_", scenario))
  result_files <- list.files(folder, pattern = "^res_sim_\\d+\\.RDS$", full.names = TRUE)

  if (length(result_files) == 0) return(NULL)

  out <- bind_rows(lapply(result_files, function(res_file) {
    sim_res <- readRDS(res_file)
    run_metrics(sim_res, scenario)
  }))
  gc()
  out
}

plan(multisession, workers = workers)
all_metrics <- future_pmap(params, get_metrics)
plan(sequential)

metrics <- bind_rows(all_metrics[!sapply(all_metrics, is.null)])

# save
saveRDS(metrics, out_file)
print(paste0("Collection complete! ", nrow(metrics), " metric rows from ",
             n_distinct(metrics$run), " runs x ", n_distinct(metrics$scenario), " scenarios"))
