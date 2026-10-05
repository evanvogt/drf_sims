#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=2:ompthreads=1:mem=3gb
#PBS -J 1-4800%100
#PBS -N cts_corr_rerun
#PBS -o logs_rerun/
#PBS -j oe


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

jobid=$(sed -n "${PBS_ARRAY_INDEX}p" "${PBS_O_WORKDIR}/failed_ids.txt")

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript cts_corr_analysis.R "$jobid" 1 1