#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=60:mem=150gb
#PBS -N h_lasso_all

module purge

module add tools/prod

module add R/4.2.1-foss-2022a

eval "$(~/miniforge3/bin/conda shell.bash hook)"

conda activate drf-env

for i in {1..11}; do
  for j in 250 500 1000 5000; do

    cd "/rds/general/user/evanvogt/projects/nihr_drf_simulations/live/scripts/H_lasso"

    scenario="scenario_${i}"
    n="${j}"

    Rscript hlasso.R "$scenario" "$n"
  done
done

