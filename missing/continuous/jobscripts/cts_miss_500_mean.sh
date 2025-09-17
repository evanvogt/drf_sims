#!/bin/bash
#PBS -l walltime=03:00:00  
#PBS -l select=1:ncpus=15:ompthreads=15:mem=10gb
#PBS -J 1-10000
#PBS -N miss_cts_mean_500
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/missing/continuous/jobscripts/logs_500_mean
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/missing/continuous/jobscripts/logs_500_mean

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Parameter arrays
n=500
scenarios=(1 2 3 4 5)
miss_types=("prognostic" "predictive" "both")
miss_props=(0.2 0.4)


# 0-based index for calculation
id=$((PBS_ARRAY_INDEX - 1))

scen_id=$((id / 2000))
remainder=$((id % 2000))

miss_prop_id=$((remainder / 1000))
sim_id=$((remainder % 1000 + 1)) 


# Get actual values from arrays
scenario=${scenarios[$scen_id]}
miss_prop=${miss_props[$miss_prop_id]}

# Display the mapping (useful for debugging)
echo "PBS_ARRAY_INDEX: $PBS_ARRAY_INDEX"
echo "Sample Size: $n"
echo "Scenario: $scenario" 
echo "Miss Proportion: $miss_prop"
echo "Simulation Run: $sim_id"

for miss_type in "${miss_types[@]}"; do
    echo "Miss Type: $miss_type"
    # Navigate to the script directory
    cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/missing/continuous"
    
    # Run the R script with all parameters
    Rscript cts_miss_analysis.R "$scenario" "$n" "$sim_id" "$miss_type" "$miss_prop" "mean_imputation"
done