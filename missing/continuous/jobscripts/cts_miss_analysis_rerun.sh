#!/bin/bash
#PBS -l walltime=01:00:00  
#PBS -l select=1:ncpus=2:ompthreads=2:mem=5gb
#PBS -J 1-9
#PBS -N cts_miss_rerun
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/missing/continuous/logs_rerun/
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/missing/continuous/logs_rerun/

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# failed job ids
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/missing/continuous/jobscripts"
jobid=$(sed -n "${PBS_ARRAY_INDEX}p" failed_ids.txt)

echo "rerunning index: $jobid"

# Navigate to script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/missing/continuous"

# Run R script with parameters
Rscript cts_miss_analysis.R "$jobid"
