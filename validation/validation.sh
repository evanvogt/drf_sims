#PBS -l walltime=00:10:00  
#PBS -l select=1:ncpus=5:ompthreads=5:mem=10gb
#PBS -J 1-300
#PBS -N validation_cts
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/validation/logs
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/validation/logs

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# arrays
INTERIM_PROPS=(0.25 0.5 0.75)

# Convert PBS_ARRAY_INDEX (1-based) to 0-based for calculations
IDX=$((PBS_ARRAY_INDEX - 1))

# Calculate parameters using integer division and modulo

PROP_IDX=$((IDX / 100)) 
SIM_RUN=$((IDX % 100 + 1))     

# Get actual values from arrays
INTERIM_PROP=${INTERIM_PROPS[$PROP_IDX]}

# Display the mapping (useful for debugging)
echo "PBS_ARRAY_INDEX: $PBS_ARRAY_INDEX"
echo "IDX: $IDX"
echo "PROP_IDX: $PROP_IDX"
echo "Interim proportion: $INTERIM_PROP"
echo "Simulation Run: $SIM_RUN"

# Navigate to the script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/validation"

# Run the R script for the assigned scenario and sample size
Rscript cate_validation_cts.R "$SIM_RUN" "$INTERIM_PROP" 
