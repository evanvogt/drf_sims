#PBS -l walltime=2:00:00  
#PBS -l select=1:ncpus=10:mem=30gb
#PBS -N cate_corr_mse

module purge
module add tools/prod
module add R/4.2.1-foss-2022a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate drf-env

cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/continuous/performance"


Rscript cate_corr_mse.R
