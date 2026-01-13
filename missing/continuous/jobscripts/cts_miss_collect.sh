#!/bin/bash
#PBS -l walltime=01:00:00  
#PBS -l select=1:ncpus=2:ompthreads=2:mem=10gb
#PBS -N cts_miss_collect
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/missing/continuous/jobscripts/
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/missing/continuous/jobscripts/

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/missing/continuous"

# Run R script with parameters
Rscript cts_miss_collect.R 
