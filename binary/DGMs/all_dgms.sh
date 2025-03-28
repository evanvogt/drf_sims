#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=10:mem=20gb
#PBS -N all_dgms

module purge
module add tools/prod
module add R/4.2.1-foss-2022a

eval "$(~/miniforge3/bin/conda shell.bash hook)"

conda activate drf-env

for i in {1..10}; do

  cd /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/DGMs
  
  Rscript scenario_${i}.R

done