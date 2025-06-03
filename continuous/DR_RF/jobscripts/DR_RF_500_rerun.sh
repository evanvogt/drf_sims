#PBS -l walltime=00:15:00  
#PBS -l select=1:ncpus=10:ompthreads=10:mem=5gb
#PBS -J 1-1000
#PBS -N DR_RF_500_rerun
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/continuous/DR_RF/jobscripts/logs_500/
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/continuous/DR_RF/jobscripts/logs_500/

module purge
module add tools/prod
module add R/4.2.1-foss-2022a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate drf-env

# Compute indices from the array job ID
sim_id=$(((PBS_ARRAY_INDEX-1) % 1000 + 1)) # 1-1000

scenario="scenario_3"
n="500"

echo "running: ${scenario}_${n}, simulation: $sim_id"

# script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/continuous/DR_RF"

# Run the R script for the assigned scenario and sample size
Rscript DR_RF_sim.R "$scenario" "$n" "$sim_id"
