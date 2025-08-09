#####################
# title: collect all the failed job ids and rewrite array scripts for rerun
# date started: 24/04/2025
# date finished:
# author: Ellie Van Vogt
####################
# libraries
library(dplyr)
library(stringr)
library(purrr)
library(hms)
library(bench)

# paths
setwd("/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary")

# get all models
models <- c("CF", "DR_oracle", "DR_RF", "DR_semi_oracle", "DR_superlearner", "T_RF")
sample_sizes <- c("100", "250", "500", "1000")

for (model in models) {
  jobs_path <- file.path(model, "jobscripts")
  
  failed_files <- list.files(jobs_path, pattern = "^failed_\\d+\\.txt$", full.names = TRUE)
  
  if (length(failed_files) == 0) {
    message(paste("No failed files for", model, "- skipping"))
    next
  }
  
  for (failed_file in failed_files) {
    fname <- basename(failed_file)
    n_val <- str_extract(fname, "(?<=failed_)\\d+")
    
    print(paste("Processing", model, "sample size", n_val))
    
    failed_ids <- readLines(failed_file)
    num_failed <- length(failed_ids)
    
    if (num_failed == 0) {
      message(paste("No failed jobs for", model, "sample size", n_val, "- skipping"))
      next
    }
    
    # Read the original simulation array script template for this model
    original_script_path <- paste0(model, "/jobscripts/", model, "_", n_val, "_array.sh")
    
    if (!file.exists(original_script_path)) {
      message(paste("Original simulation script not found:", original_script_path, "- skipping"))
      next
    }
    
    old_script <- readLines(original_script_path)
    
    # Extract the old job parameters
    old_walltime <- as_hms(sub("^#PBS -l walltime=([0-9:]+).*", "\\1", grep("^#PBS -l walltime=", old_script, value = TRUE)))
    resource_line <- grep("^#PBS -l select=", old_script, value = TRUE)
    old_ncpus <- as.numeric(sub(".*ncpus=(\\d+).*", "\\1", resource_line))
    old_mem <- as.numeric(gsub("gb", "", sub(".*mem=([0-9a-zA-Z]+).*", "\\1", resource_line)))
    
    # Increase parameters by 50%
    new_walltime <- hms(as.numeric(old_walltime) * 1.5)
    new_ncpus <- ceiling(old_ncpus * 1.5)
    new_mem <- paste0(ceiling(old_mem * 1.5), "gb")
    
    # Create the rerun script
    script_path <- file.path(jobs_path, paste0(model, "_", n_val, "_rerun.sh"))
    
    new_script <- old_script %>%
      map_chr(~{
        if (str_detect(., "^#PBS -l walltime=")) {
          paste0("#PBS -l walltime=", new_walltime)
        } else if (str_detect(., "^#PBS -l select=")) {
          paste0("#PBS -l select=1:ncpus=", new_ncpus, ":mem=", new_mem)
        } else if (str_detect(., "^#PBS -J")) {
          paste0("#PBS -J 1-", num_failed)
        } else if (str_detect(., "^#PBS -N")) {
          paste0("#PBS -N ", model, "_", n_val, "_rerun")
        } else if (str_detect(., "^echo \"running")) {
          # Keep the original echo but indicate it's a rerun
          str_replace(., "running", "rerunning failed simulation")
        } else if (str_detect(., "^Rscript")) {
          # Keep the original Rscript call - it should work with the original job ID
          .
        } else {
          .
        }
      })
    
    # Find where to insert the new job ID lookup logic
    # Insert after the conda activation line
    conda_line_idx <- max(grep("conda activate", new_script))
    
    new_lines <- c(
      "",
      "# Get the original failed job ID from the consolidated list",
      paste0('cd "', normalizePath(jobs_path), '"'),
      paste0('original_job_id=$(sed -n "${PBS_ARRAY_INDEX}p" failed_', n_val, '.txt)'),
      "",
      "# Use the original job ID for the simulation",
      "# Override PBS_ARRAY_INDEX to use the original failed job ID",
      "export PBS_ARRAY_INDEX=$original_job_id",
      'echo "Rerunning original simulation job ID: $original_job_id"',
      ""
    )
    
    new_script <- append(new_script, new_lines, after = conda_line_idx)
    
    # Keep the original parameter grid and calculations since they're needed for the simulation
    # The simulation script will use the overridden PBS_ARRAY_INDEX
    
    writeLines(new_script, script_path)
    print(paste0("Rerun script for ", model, "_", n_val, " created with ", num_failed, " failed jobs"))
  }
}