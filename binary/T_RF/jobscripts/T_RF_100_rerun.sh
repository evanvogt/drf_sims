#PBS -l walltime=01:07:30
#PBS -l select=1:ncpus=15:mem=11gb
#PBS -J 1-200
#PBS -N T_RF_100_rerun
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/T_RF/jobscripts/logs_100/
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/T_RF/jobscripts/logs_100/

module purge
module add tools/prod
module add R/4.2.1-foss-2022a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate drf-env

# Get the original failed job ID from the consolidated list
cd "/rds/general/project/nihr_drf_simulations/live/scripts/binary/T_RF/jobscripts"
original_job_id=$(sed -n "${PBS_ARRAY_INDEX}p" failed_100.txt)

# Use the original job ID for the simulation
# Override PBS_ARRAY_INDEX to use the original failed job ID
export PBS_ARRAY_INDEX=$original_job_id
echo "Rerunning original simulation job ID: $original_job_id"


# scenarios
scenarios=(1 2 3 4 5 6 7 8 9 10)


# Compute indices from the array job ID
sim_id=$(((PBS_ARRAY_INDEX-1) % 1000 + 1)) # 1-1000
scen_id=$(((PBS_ARRAY_INDEX - 1) / 1000))  # 0-10


scenario="scenario_${scenarios[$scen_id]}"
n="100"

echo "rerunning failed simulation: scenario_${scenario}_${n}, simulation: $sim_id"

# Navigate to the script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/T_RF"

# Run the R script for the assigned scenario and sample size
Rscript T_RF_sim.R "$scenario" "$n" "$sim_id"
