#!/bin/bash
# The oracle add-on (snr_oracle_config.R) on every row of the grid. The array
# steps by 100, so subjob k gets PBS_ARRAY_INDEX = 100(k - 1) + 1 and fits rows
# PBS_ARRAY_INDEX .. +99: one cell, 43 subjobs. Each row is a readRDS and two
# grf forests, no RSF: locally (R 4.5.3) about 2 s a row at n = 500. The
# n = 4000 cells are the long ones; untimed, so walltime has headroom.
# Submit after the main study's check/rerun: a row with no main result file
# is skipped and left for snr_oracle_rerun.sh.
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=4gb
#PBS -J 1-4300:100
#PBS -N snr_oracle
#PBS -o logs_oracle/
#PBS -j oe


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters: first row, rows per subjob, grf threads
Rscript snr_oracle_analysis.R "$PBS_ARRAY_INDEX" 100 2
