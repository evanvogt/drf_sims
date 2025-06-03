#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=10:mem=20gb
#PBS -N collect_all_data_bin

module purge
module add tools/prod
module add R/4.2.1-foss-2022a

eval "$(~/miniforge3/bin/conda shell.bash hook)"

conda activate drf-env

cd /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/DGMs

Rscript all_data.R
