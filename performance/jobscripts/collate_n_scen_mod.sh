#PBS -l walltime=1:00:00  
#PBS -l select=1:ncpus=5:mem=5gb
#PBS -J 1-160
#PBS -N col_n_sim_mod
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/performance/jobscripts/logs/
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/performance/jobscripts/logs/

module purge
module add tools/prod
module add R/4.2.1-foss-2022a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate drf-env

# parameter grid
SCENARIOS=(1 2 3 4 5 6 7 8 9 10)
SAMPLE_SIZES=(250 500 1000 5000)
MODELS=(CF DR_oracle DR_RF T_RF)

# numbers for dividing the index
TOTAL_SCENARIOS=${#SCENARIOS[@]}
TOTAL_SAMPLE_SIZES=${#SAMPLE_SIZES[@]}
TOTAL_MODELS=${#MODELS[@]}

IDX=$((PBS_ARRAY_INDEX - 1))  # Convert to zero-based index
SCEN_IDX=$((IDX / (TOTAL_SAMPLE_SIZES * TOTAL_MODELS)))
REMAINING=$((IDX % (TOTAL_SAMPLE_SIZES * TOTAL_MODELS)))
SAMPLE_IDX=$((REMAINING / TOTAL_MODELS))
MODEL_IDX=$((REMAINING % TOTAL_MODELS))

SCENARIO="scenario_${SCENARIOS[$SCEN_IDX]}"
SAMPLE_SIZE=${SAMPLE_SIZES[$SAMPLE_IDX]}
MODEL=${MODELS[$MODEL_IDX]}


# Navigate to the script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/performance"

echo "running ${SCENARIO}_${SAMPLE_SIZE} ${MODEL}"

# Run R script with parameters
Rscript collate_n_scen_mod.R "$SCENARIO" "$SAMPLE_SIZE" "$MODEL"