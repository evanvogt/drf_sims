#PBS -l walltime=01:00:00  
#PBS -l select=1:ncpus=10:ompthreads=10:mem=15gb
#PBS -N bin_results_collect
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/binary/all_methods/jobscripts
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/binary/all_methods/jobscripts

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/binary/all_methods"

# Run the R script for the assigned scenario and sample size
Rscript collect_results_bin.R
