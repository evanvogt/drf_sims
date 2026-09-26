#!/bin/bash
# Scenario 6 (main-study scenario 3) - grid rows 9901-12600 of
# cts_miss_config.R. Resources are cts_miss_1.sh's. Needs logs_extra/ to exist.
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=2gb
#PBS -J 9901-12600%190
#PBS -N cts_miss_extra
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
Rscript cts_miss_analysis.R "${PBS_ARRAY_INDEX}"
