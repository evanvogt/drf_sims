#!/bin/bash
#PBS -l walltime=00:05:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=3gb
#PBS -J 1-1100
#PBS -N cts_val_1
#PBS -o logs_1/
#PBS -e logs_1/

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script - PBS array index is a row number of study$grid (cts_val_config.R)
Rscript cts_val_analysis.R "$PBS_ARRAY_INDEX" 1 1
