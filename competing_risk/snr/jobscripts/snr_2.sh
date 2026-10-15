#!/bin/bash
# n = 1000 and 2000, hte_mult 1. Locally (R 4.5.3) n = 2000 took about 4 min
# and peaked at 8gb, nearly all of it rfsrc (README.md, Resources).
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=16gb
#PBS -J 401-600
#PBS -N snr_2
#PBS -o logs_2/
#PBS -j oe


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript snr_analysis.R "$PBS_ARRAY_INDEX"
