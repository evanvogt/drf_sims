#!/bin/bash
#PBS -l walltime=01:00:00  
#PBS -l select=1:ncpus=15:ompthreads=15:mem=10gb
#PBS -J 1-1200
#PBS -N cate_surv
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/competing_risk/sub_dist/logs
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/competing_risk/sub_dist/logs

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Define all parameter arrays
SAMPLE_SIZES=(250 500 1000)  
SCENARIOS=(1 2 3 4)

# Convert PBS_ARRAY_INDEX (1-based) to 0-based for calculations
IDX=$((PBS_ARRAY_INDEX - 1))

# Calculate parameters using integer division and modulo
SAMPLE_IDX=$((IDX / 400)) # 4 scenarios, 100 sims
REMAINING=$((IDX % 400))

SCENARIO_IDX=$((REMAINING / 100))
SIM_RUN=$((REMAINING % 100 + 1))

# Get actual values from arrays
SAMPLE_SIZE=${SAMPLE_SIZES[$SAMPLE_IDX]}
SCENARIO=${SCENARIOS[$SCENARIO_IDX]}

# Display the mapping (useful for debugging)
echo "PBS_ARRAY_INDEX: $PBS_ARRAY_INDEX"
echo "Sample Size: $SAMPLE_SIZE"
echo "Scenario: $SCENARIO" 
echo "Simulation Run: $SIM_RUN"

# Navigate to the script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/competing_risk/"

# Run the R script with all parameters
Rscript cate_surv_run.R "$SCENARIO" "$SAMPLE_SIZE" "$SIM_RUN"