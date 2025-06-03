#PBS -l walltime=3:00:00  
#PBS -l select=1:ncpus=30:ompthreads=30:mem=30gb
#PBS -J 1-44  
#PBS -N h_lasso_array
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/continuous/H_lasso/jobscripts/logs/
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/continuous/H_lasso/jobscripts/logs/

module purge
module add tools/prod
module add R/4.2.1-foss-2022a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate drf-env

# Define the scenarios and sample sizes
SCENARIOS=(1 2 3 4 5 6 7 8 9 10 11)
SAMPLESIZES=(250 500 1000 5000)

# Compute indices from the array job ID
IDX=$((PBS_ARRAY_INDEX - 1))
SCENARIO_IDX=$((IDX / 4))  # 4 sample sizes per scenario
SIZE_IDX=$((IDX % 4))

SCENARIO="scenario_${SCENARIOS[$SCENARIO_IDX]}"
N="${SAMPLESIZES[$SIZE_IDX]}"

# Navigate to the script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/continuous/H_lasso"

# Run the R script for the assigned scenario and sample size
Rscript hlasso.R "$SCENARIO" "$N"
