#!/bin/bash
# Sized for the largest cell (n = 4000, as snr_3.sh), since a rerun list can
# mix every n. study.R check rewrites -J to match failed_ids.txt.
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=3:ompthreads=2:mem=64gb
#PBS -J 1-208%100
#PBS -N snr_rerun
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
Rscript snr_analysis.R "$jobid"
