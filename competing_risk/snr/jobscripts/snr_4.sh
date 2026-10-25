#!/bin/bash
# Scenarios 1-2 (n = 500, hte_mult 1) and 3-6 (n = 500, hte_mult 1 / 1.5 / 2 / 3).
# As snr_1.sh (scenario 7): about 30 s and 1.3gb locally.
#PBS -l walltime=00:15:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=4gb
#PBS -J 701-2500
#PBS -N snr_4
#PBS -o logs_4/
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
