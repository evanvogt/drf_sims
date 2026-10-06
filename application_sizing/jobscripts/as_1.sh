#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:ompthreads=8:mem=32gb
#PBS -J 1-201%50
#PBS -N as_1
#PBS -o logs/
#PBS -j oe

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters: grid row, 4 future workers x 2 grf threads
Rscript as_analysis.R "$PBS_ARRAY_INDEX" 4 2
