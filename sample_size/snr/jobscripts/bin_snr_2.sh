#!/bin/bash
# n = 1000 / 2000 / 4000, scenarios 1-4 at hte_mult 1. Locally (R 4.5.3,
# 2 threads) n = 4000 took about 1.3 min and peaked under 1gb.
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=4gb
#PBS -J 1301-2500
#PBS -N bin_snr_2
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
Rscript bin_snr_analysis.R "$PBS_ARRAY_INDEX" 2
