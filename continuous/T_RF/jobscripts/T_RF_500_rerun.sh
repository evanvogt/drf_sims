#PBS -l walltime=00:22:30
#PBS -l select=1:ncpus=15:ompthreads=15:mem=8gb
#PBS -J 1-250
#PBS -N T_RF_500_rerun
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/continuous/T_RF/jobscripts/logs_500/
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/continuous/T_RF/jobscripts/logs_500/

module purge
module add tools/prod
module add R/4.2.1-foss-2022a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate drf-env

# scenarios
scenarios=(1 2 3 4 5 6 7 8 9 10)

#get the failed job ids
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/continuous/T_RF/jobscripts"
jobid=$(sed -n "${PBS_ARRAY_INDEX}p" combined_failed_500.txt)

# Compute indices from the array job ID
sim_id=$(((jobid - 1) % 1000 + 1)) # 1-1000
scen_id=$(((jobid - 1) / 1000))  # 0-10


scenario="scenario_${scenarios[$scen_id]}"
n="500"

echo "running: scenario_${scenario}_${n}, simulation: $sim_id"

# Navigate to the script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/continuous/T_RF"

# Run the R script for the assigned scenario and sample size
Rscript T_RF_sim.R "$scenario" "$n" "$sim_id"
