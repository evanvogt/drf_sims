#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=5:mem=1gb
#PBS -J 1-10
#PBS -N run_all_dgms
#PBS -o /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/DGMs/jobscripts/logs/
#PBS -e /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/DGMs/jobscripts/logs/


module purge
module add tools/prod
module add R/4.2.1-foss-2022a

eval "$(~/miniforge3/bin/conda shell.bash hook)"

conda activate drf-env

cd /rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/binary/DGMs

samplesizes=(100 250 500 1000)
  
Rscript scenario_${PBS_ARRAY_INDEX}.R "${samplesizes[@]}"