#PBS -l walltime=00:45:00  
#PBS -l select=1:ncpus=40:ompthreads=40:mem=30gb
#PBS -J 1-59
#PBS -N CF_5000_reruns
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/CF/jobscripts/logs_5000/
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/CF/jobscripts/logs_5000/

module purge
module add tools/prod
module add R/4.2.1-foss-2022a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate drf-env

# scenarios
scenarios=(1 2 3 4 5 6 7 8 9 10)

#get the failed job ids
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/CF/jobscripts"
jobid=$(sed -n "${PBS_ARRAY_INDEX}p" combined_failed_5000.txt)

# Compute indices from the array job ID
sim_id=$(((jobid-1) % 1000 + 1)) # 1-1000
scen_id=$(((jobid - 1) / 1000))  # 0-10


scenario="scenario_${scenarios[$scen_id]}"
n="5000"

echo "running: scenario_$scenario_$n, simulation: $sim_id"

# Navigate to the script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/CF"

# Run the R script for the assigned scenario and sample size
Rscript CF_sim.R "$scenario" "$n" "$sim_id"
