#!/bin/bash
# n = 500, hte_mult 1 / 1.5 / 2 / 3. Locally (R 4.5.3) about 30 s and 1.3gb.
#PBS -l walltime=00:15:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=4gb
#PBS -J 1-400
#PBS -N snr_1
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
Rscript snr_analysis.R "$PBS_ARRAY_INDEX"
