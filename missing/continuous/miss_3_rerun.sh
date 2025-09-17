#!/bin/bash
#PBS -l walltime=01:00:00  
#PBS -l select=1:ncpus=20:ompthreads=20:mem=10gb
#PBS -J 1-1800
#PBS -N missing_cts_scen3
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/missing/continuous/logs
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/missing/continuous/logs

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Define parameter arrays
MISS_METHOD=("complete_cases" "mean_imputation" "multiple_imputation")  
MISS_TYPE=("prognostic" "predictive" "both")
MISS_PROP=(0.2 0.4)

# Scenario fixed to 3
SCENARIO=3

# Total jobs:
#get the failed job ids
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/missing/continuous"
jobid=$(sed -n "${PBS_ARRAY_INDEX}p" failed_ids.txt)

# Convert PBS_ARRAY_INDEX (1-based) to 0-based for calculations
IDX=$((jobid - 1))

MISS_METHOD_IDX=$((IDX / 600))          # 600 = 3×2×100
REMAINING=$((IDX % 600))

MISS_TYPE_IDX=$((REMAINING / 200))      # 200 = 2×100
REMAINING=$((REMAINING % 200))

MISS_PROP_IDX=$((REMAINING / 100))      # 100 sims per combo
SIM_RUN=$((REMAINING % 100 + 1))       

# Get actual values
SAMPLE_SIZE=1000
MISS_METHOD_VAL=${MISS_METHOD[$MISS_METHOD_IDX]}
MISS_TYPE_VAL=${MISS_TYPE[$MISS_TYPE_IDX]}
MISS_PROP_VAL=${MISS_PROP[$MISS_PROP_IDX]}

# Debug info
echo "PBS_ARRAY_INDEX: $PBS_ARRAY_INDEX"
echo "Sample Size: $SAMPLE_SIZE"
echo "Scenario: $SCENARIO"
echo "Miss Type: $MISS_TYPE_VAL"
echo "Miss Proportion: $MISS_PROP_VAL"
echo "Miss Method: $MISS_METHOD_VAL"
echo "Simulation Run: $SIM_RUN"

# Navigate to script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/missing/continuous"

# Run R script with parameters
Rscript cate_analysis_miss_cts.R "$SCENARIO" "$SAMPLE_SIZE" "$SIM_RUN" "$MISS_TYPE_VAL" "$MISS_PROP_VAL" "$MISS_METHOD_VAL"
