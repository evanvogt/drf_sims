#!/bin/bash
# Sized as bin_snr_2.sh (the largest n), since a rerun list can mix every n.
# bin_snr_check.R rewrites -J to match failed_ids.txt.
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=3:ompthreads=2:mem=4gb
#PBS -J 1-11%100
#PBS -N bin_snr_rerun
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
Rscript bin_snr_analysis.R "$jobid" 2
