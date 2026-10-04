#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=4gb
#PBS -J 1-3000%200
#PBS -N se_1
#PBS -o logs_1/
#PBS -j oe


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript se_analysis.R "$PBS_ARRAY_INDEX"
