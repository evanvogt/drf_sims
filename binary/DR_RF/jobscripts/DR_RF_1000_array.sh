#PBS -l walltime=3:00:00  
#PBS -l select=1:ncpus=20:ompthreads=20:mem=10gb
#PBS -J 1-10000
#PBS -N DR_RF_1000_array
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/DR_RF/jobscripts/logs_1000/
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/DR_RF/jobscripts/logs_1000/

module purge
module add tools/prod
module add R/4.2.1-foss-2022a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate drf-env

# scenarios
scenarios=(1 2 3 4 5 6 7 8 9 10)


# Compute indices from the array job ID
sim_id=$(((PBS_ARRAY_INDEX-1) % 1000 + 1)) # 1-1000
scen_id=$(((PBS_ARRAY_INDEX - 1) / 1000))  # 0-10


scenario="scenario_${scenarios[$scen_id]}"
n="1000"

echo "running: scenario_${scenario}_${n}, simulation: $sim_id"

# Navigate to the script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/DR_RF"

# Run the R script for the assigned scenario and sample size
Rscript DR_RF_sim.R "$scenario" "$n" "$sim_id"
