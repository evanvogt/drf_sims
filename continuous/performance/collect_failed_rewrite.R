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
setwd("/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/continuous")

# get all models
models <- c("CF", "DR_oracle", "DR_RF", "T_RF")

for (model in models) {
  jobs_path <- file.path(model, "jobscripts")
  failed_files <- list.files(jobs_path, pattern = "^failed_\\d+_scenario_\\d+\\.txt$", full.names = TRUE)
  
  if (length(failed_files) == 0) {
    message(paste("No failed files for", model, "- skipping"))
    next
  }
  # group by sample size
  failed_info <- tibble(file = failed_files) %>%
    mutate(
      fname = basename(file),
      n = str_extract(fname, "(?<=failed_)\\d+"),
      scenario = str_extract(fname, "scenario_\\d+")
    )
  
  # combine files by sample size
  failed_info %>%
    group_by(n) %>%
    group_walk(~{
      n_val <- .y[["n"]]
      print(n_val)
      combined_path <- file.path(jobs_path, paste0("combined_failed_", n_val, ".txt"))
      all_lines <- .x$file %>% map(readLines) %>% unlist() %>% sort()
      num_failed <- length(all_lines)
  
      writeLines(as.character(all_lines), combined_path) # new file containing all the failed ids
      
      old_script <- readLines(paste0(model, "/jobscripts/", model, "_", n_val, "_array.sh"))
      
      # extract the old job parameters
      old_walltime <- as_hms(sub("^#PBS -l walltime=([0-9:]+).*", "\\1", grep("^#PBS -l walltime=", old_script, value = TRUE)))
      resource_line <- grep("^#PBS -l select=", old_script, value = TRUE)
      old_ncpus <- as.numeric(sub(".*ncpus=(\\d+).*", "\\1", resource_line))
      old_ompthreads <- as.numeric(sub(".*ompthreads=(\\d+).*", "\\1", resource_line))
      old_mem <- as.numeric(gsub("gb", "", sub(".*mem=([0-9a-zA-Z]+).*", "\\1", resource_line)))
      # increase parameters by 50%
      new_walltime <- hms(as.numeric(old_walltime)*1.5)
      
      new_ncpus <- ceiling(old_ncpus*1.5)
      new_ompthreads <- new_ncpus
      new_mem <- paste0(ceiling(old_mem*1.5), "gb")
      
      # creating the rerun script
      script_path <- file.path(jobs_path, paste0(model, "_", n_val, "_rerun.sh"))
      n_val <- .y[["n"]]
      new_script <- old_script %>%
        map_chr(~{
          if (str_detect(., "^#PBS -l walltime=")) paste0("#PBS -l walltime=", new_walltime)
          else if (str_detect(., "^#PBS -l select=")) paste0("#PBS -l select=1:ncpus=", new_ncpus,
                                                             ":ompthreads=", new_ompthreads,
                                                             ":mem=", new_mem)
          else if (str_detect(., "^#PBS -J")) paste0("#PBS -J ", "1-", num_failed)
          else if (str_detect(., "^#PBS -N")) paste0("#PBS -N ", model, "_", n_val, "_rerun")
          else if (str_detect(., "^sim_id=")) "sim_id=$(((jobid - 1) % 1000 + 1)) # 1-1000"
          else if (str_detect(., "^scen_id=")) "scen_id=$(((jobid - 1) / 1000))  # 0-10"
          else .
        })
      
      # add in new bits:
      combined_directory <- paste0("/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/continuous/", jobs_path)
      
      new_lines <- c("#get the failed job ids",
                     paste0('cd "', combined_directory, '"'),
                     paste0('jobid=$(sed -n "${PBS_ARRAY_INDEX}p" ', "combined_failed_", .y$n, ".txt)"))
      
      insert_place <- max(grep("scenarios=", new_script)) + 1
      new_script <- append(new_script, new_lines, after = insert_place)
      
      writeLines(new_script, script_path)
      print(paste0("rerun script for ", model, "_", n_val, " finished - check"))
    })
}
