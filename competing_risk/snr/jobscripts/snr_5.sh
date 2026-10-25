#!/bin/bash
# Scenarios 1-6, n = 1000 and 2000, hte_mult 1. Scenario 7 locally (R 4.5.3)
# took about 4 min and peaked at 8gb at n = 2000, nearly all of it rfsrc
# (README.md, Resources); 12gb leaves headroom for the other scenarios.
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=12gb
#PBS -J 2501-3700
#PBS -N snr_5
#PBS -o logs_5/
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
