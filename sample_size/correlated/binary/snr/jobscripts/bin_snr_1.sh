#!/bin/bash
# n = 500: scenario 1 at x1, scenarios 2-4 at hte_mult 1 / 1.5 / 2 / 3.
# Locally (R 4.5.3, 2 threads) about 10 s and well under 1gb.
#PBS -l walltime=00:15:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=2gb
#PBS -J 1-1300
#PBS -N bin_snr_1
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
Rscript bin_snr_analysis.R "$PBS_ARRAY_INDEX" 2
