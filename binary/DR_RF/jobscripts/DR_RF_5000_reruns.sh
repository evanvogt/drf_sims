#PBS -l walltime=4:00:00  
#PBS -l select=1:ncpus=45:ompthreads=45:mem=30gb
#PBS -J 1-816
#PBS -N DR_RF_5000_rerun
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/DR_RF/jobscripts/logs_5000/
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/DR_RF/jobscripts/logs_5000/

module purge
module add tools/prod
module add R/4.2.1-foss-2022a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate drf-env

# scenarios
scenarios=(1 2 3 4 5 6 7 8 9 10)

#get the failed job ids
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/DR_RF/jobscripts"
jobid=$(sed -n "${PBS_ARRAY_INDEX}p" failed_5000_scenario_10.txt)

# Compute indices from the array job ID
sim_id=$(((jobid-1) % 1000 + 1)) # 1-1000
scen_id=$(((jobid - 1) / 1000))  # 0-10


scenario="scenario_${scenarios[$scen_id]}"
n="5000"

echo "running: scenario_${scenario}_${n}, simulation: $sim_id"

# Navigate to the script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/DR_RF"

# Run the R script for the assigned scenario and sample size
Rscript DR_RF_sim.R "$scenario" "$n" "$sim_id"
