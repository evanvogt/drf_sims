#PBS -l walltime=00:30:00  
#PBS -l select=1:ncpus=5:mem=10gb
#PBS -J 1-24
#PBS -N res_comb_bin
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/performance/jobscripts/logs/
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/performance/jobscripts/logs/

module purge
module add tools/prod
module add R/4.2.1-foss-2022a
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate drf-env

# parameter grid (scenarios removed since R script handles all scenarios)
SAMPLE_SIZES=(100 250 500 1000)
MODELS=(CF DR_oracle DR_semi_oracle DR_superlearner DR_RF T_RF)

# numbers for dividing the index
TOTAL_SAMPLE_SIZES=${#SAMPLE_SIZES[@]}
TOTAL_MODELS=${#MODELS[@]}

IDX=$((PBS_ARRAY_INDEX - 1))  # Convert to zero-based index

# Calculate indices for sample size and model
SAMPLE_IDX=$((IDX / TOTAL_MODELS))
MODEL_IDX=$((IDX % TOTAL_MODELS))

SAMPLE_SIZE=${SAMPLE_SIZES[$SAMPLE_IDX]}
MODEL=${MODELS[$MODEL_IDX]}

# Navigate to the script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/performance"

echo "running all scenarios for sample size ${SAMPLE_SIZE} and model ${MODEL}"

# Run R script with parameters (no scenario parameter needed)
Rscript 01_collate_n_scen_mod.R "$SAMPLE_SIZE" "$MODEL"