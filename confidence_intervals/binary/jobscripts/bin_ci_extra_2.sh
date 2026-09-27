#!/bin/bash
# Runs 101-500 for scenarios 1-4 - grid rows 20001-52000 of
# bin_ci_config.R, split over bin_ci_extra_{1..4}.sh to stay within the
# 10,000-subjob array limit. This is part 2: rows 30001-40000.
# Resources are bin_ci_1.sh's. Needs logs_extra/ to exist.
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=5gb
#PBS -J 30001-40000%100
#PBS -N ci_bin_extra_2
#PBS -o logs_extra/
#PBS -e logs_extra/

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript bin_ci_analysis.R "$PBS_ARRAY_INDEX"
