#PBS -l walltime=00:45:00  
#PBS -l select=1:ncpus=30:ompthreads=30:mem=20gb
#PBS -J 1-5
#PBS -N CF_1000_reruns
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/CF/jobscripts/logs_1000/
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/CF/jobscripts/logs_1000/

module purge
module add tools/prod
module add R/4.2.1-foss-2022a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate drf-env

failed=(4363 4364 4365 4387 4409)
id=${failed[$((PBS_ARRAY_INDEX - 1))]}
# scenarios
scenarios=(1 2 3 4 5 6 7 8 9 10)


# Compute indices from the array job ID
sim_id=$(((id-1) % 1000 + 1)) # 1-1000
scen_id=$(((id - 1) / 1000))  # 0-10


scenario="scenario_${scenarios[$scen_id]}"
n="1000"

echo "running: scenario_$scenario_$n, simulation: $sim_id"

# Navigate to the script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/CF"

# Run the R script for the assigned scenario and sample size
Rscript CF_sim.R "$scenario" "$n" "$sim_id"
