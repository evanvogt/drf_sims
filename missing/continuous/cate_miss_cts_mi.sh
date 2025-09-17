#!/bin/bash
#PBS -l walltime=02:00:00  
#PBS -l select=1:ncpus=25:ompthreads=25:mem=10gb
#PBS -J 1-2400
#PBS -N missing_cts_mi
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/missing/continuous/logs_mi
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/missing/continuous/logs_mi

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Define all parameter arrays
SCENARIOS=(1 3 8 10)
MISS_TYPE=("prognostic" "predictive" "both")
MISS_PROP=(0.2 0.4)

# Calculate total combinations for verification
# 4 scenarios × 3 miss types × 2 miss props × 100 sims = 2400 total jobs

# Convert PBS_ARRAY_INDEX (1-based) to 0-based for calculations
IDX=$((PBS_ARRAY_INDEX - 1))

SCENARIO_IDX=$((IDX / 600))              # 600 = 3×2×100
REMAINING=$((IDX % 600))

MISS_TYPE_IDX=$((REMAINING / 200))       # 200 = 2×100
REMAINING=$((REMAINING % 200))

MISS_PROP_IDX=$((REMAINING / 100))       # 100 runs
SIM_RUN=$((REMAINING % 100 + 1))      

# Get actual values from arrays
SAMPLE_SIZE=1000 # fix at 1000 for now
SCENARIO=${SCENARIOS[$SCENARIO_IDX]}
MISS_TYPE_VAL=${MISS_TYPE[$MISS_TYPE_IDX]}
MISS_PROP_VAL=${MISS_PROP[$MISS_PROP_IDX]}

# Display the mapping (useful for debugging)
echo "PBS_ARRAY_INDEX: $PBS_ARRAY_INDEX"
echo "Sample Size: $SAMPLE_SIZE"
echo "Scenario: $SCENARIO" 
echo "Miss Type: $MISS_TYPE_VAL"
echo "Miss Proportion: $MISS_PROP_VAL"
echo "Simulation Run: $SIM_RUN"

# Navigate to the script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/missing/continuous"

# Run the R script with all parameters
Rscript cate_analysis_miss_cts.R "$SCENARIO" "$SAMPLE_SIZE" "$SIM_RUN" "$MISS_TYPE_VAL" "$MISS_PROP_VAL" "multiple_imputation"