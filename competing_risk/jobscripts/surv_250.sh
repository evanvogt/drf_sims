#PBS -l walltime=01:00:00  
#PBS -l select=1:ncpus=10:ompthreads=10:mem=10gb
#PBS -J 1-7000
#PBS -N surv_250
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/competing_risk/jobscripts/logs_250
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/competing_risk/jobscripts/logs_250

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# scenario and sample sze arrays
n=250
scenarios=(1 2 3 4 5 6 7) 


# convert to 0-based for counting along arrays
id=$((PBS_ARRAY_INDEX - 1))

# calculate parameters
scen_id=$((id/1000))
sim_id=$((id % 1000 + 1))

# Get actual values from arrays
scenario=${scenarios[$scen_id]}

# Display the mapping
echo "PBS_ARRAY_INDEX: $PBS_ARRAY_INDEX"
echo "Sample Size: $n"
echo "Scenario: $scenario" 
echo "Simulation Run: $sim_id"

# Run script
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/competing_risk"

Rscript surv_analysis.R "$scenario" "$n" "$sim_id"
