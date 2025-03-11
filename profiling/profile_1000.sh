#PBS -l walltime=3:00:00  
#PBS -l select=1:ncpus=50:ompthreads=50:mem=100gb
#PBS -N profile_1000

module purge
module add tools/prod
module add R/4.2.1-foss-2022a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate drf-env

# Navigate to the script directory
cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/profiling"

Rscript profile_1000.R