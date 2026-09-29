#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=2gb
#PBS -J 1-947%100
#PBS -N cts_rerun
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
Rscript cts_analysis.R "$jobid" 1 1