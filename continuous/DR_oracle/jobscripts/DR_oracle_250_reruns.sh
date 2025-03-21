#PBS -l walltime=00:45:00  
#PBS -l select=1:ncpus=10:ompthreads=10:mem=10gb
#PBS -J 1-368
#PBS -N DR_oracle_250_array
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/DR_oracle/jobscripts/logs_250/
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/DR_oracle/jobscripts/logs_250/

module purge
module add tools/prod
module add R/4.2.1-foss-2022a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate drf-env

# scenarios
scenarios=(1 2 3 4 5 6 7 8 9 10)

#failed job ids
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/DR_oracle/jobscripts"
jobid=$(sed -n "${PBS_ARRAY_INDEX}p" failed_250.txt)

# Compute indices from the array job ID
sim_id=$(((jobid-1) % 1000 + 1)) # 1-1000
scen_id=$(((jobid - 1) / 1000))  # 0-10


scenario="scenario_${scenarios[$scen_id]}"
n="250"

echo "rerunning: ${scenario}_${n}, simulation: $sim_id"

# Navigate to the script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/DR_oracle"

# Run the R script for the assigned scenario and sample size
Rscript DR_oracle_sim.R "$scenario" "$n" "$sim_id"
