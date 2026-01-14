#PBS -l walltime=02:00:00  
#PBS -l select=1:ncpus=25:ompthreads=25:mem=10gb
#PBS -J 1-352
#PBS -N bin_rerun
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/binary/all_methods/logs
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/binary/all_methods/logs

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# scenario and sample size arrays
SAMPLE_SIZES=(100 250 500 1000)  
SCENARIOS=(1 3 8 10) 

#get the failed job ids
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/binary/all_methods"
jobid=$(sed -n "${PBS_ARRAY_INDEX}p" failed_ids.txt)

# Convert PBS_ARRAY_INDEX (1-based) to 0-based for calculations
IDX=$((jobid - 1))

# Calculate parameters using integer division and modulo
SAMPLE_IDX=$((IDX / 400))
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
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/binary/all_methods"

# Run the R script for the assigned scenario and sample size
Rscript cate_analysis_bin.R "$SCENARIO" "$SAMPLE_SIZE" "$SIM_RUN"
