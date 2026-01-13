#PBS -l walltime=02:00:00  
#PBS -l select=1:ncpus=10:ompthreads=10:mem=10gb
#PBS -J 1-12
#PBS -N cts_250_rerun
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/continuous/all_methods/jobscripts/logs_250
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/continuous/all_methods/jobscripts/logs_250

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# scenario and sample sze arrays
n=250
scenarios=(1 2 3 4 5 6 7 8 9 10) 

# read in failed ids
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/continuous/all_methods/jobscripts"
jobid=$(sed -n "${PBS_ARRAY_INDEX}p" failed_ids_250.txt)

# convert to 0-based for counting along arrays
id=$((jobid - 1))

# calculate parameters
scen_id=$((id/1000))
sim_id=$((id % 1000 + 1))

# Get actual values from arrays
scenario=${scenarios[$scen_id]}

# Display the mapping
echo "PBS_ARRAY_INDEX: $PBS_ARRAY_INDEX"
echo "Original ID: $jobid"
echo "Sample Size: $n"
echo "Scenario: $scenario" 
echo "Simulation Run: $sim_id"

# Run script
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/drf_sims/continuous/all_methods"

Rscript cts_analysis.R "$scenario" "$n" "$sim_id"
