#PBS -l walltime=00:30:00  
#PBS -l select=1:ncpus=10:ompthreads=10:mem=10gb
#PBS -J 1-1000 
#PBS -N CF_500_reruns
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/continuous/CF/jobscripts/logs_500/
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/continuous/CF/jobscripts/logs_500/

module purge
module add tools/prod
module add R/4.2.1-foss-2022a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate drf-env

#get the failed job ids
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/continuous/CF/jobscripts"
jobid=$(sed -n "${PBS_ARRAY_INDEX}p" failed_500_scenario_6.txt)


# scenarios
scenarios=(1 2 3 4 5 6 7 8 9 10)


# Compute indices from the array job ID
sim_id=$(((jobid - 1) % 1000 + 1)) # 1-1000
scen_id=$(((jobid - 1) / 1000))  # 0-10


scenario="scenario_${scenarios[$scen_id]}"
n="500"

echo "running: scenario_${scenario}_${n}, simulation: ${sim_id}"

# Navigate to the script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/continuous/CF"

# Run the R script for the assigned scenario and sample size
Rscript CF_sim.R "$scenario" "$n" "$sim_id"
