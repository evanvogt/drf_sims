#!/bin/bash
# Runs 101-500 for scenarios 1-4 - grid rows 4001-10400 of
# bin_config.R. Resources are bin_1.sh's. Needs logs_extra/ to exist.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=5gb
#PBS -J 4001-10400%380
#PBS -N bin_extra
#PBS -o logs_extra/
#PBS -j oe


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript bin_analysis.R "$PBS_ARRAY_INDEX" 1 1
