#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=3:ompthreads=2:mem=5gb
#PBS -J 1-2%100
#PBS -N se_rerun
#PBS -o logs_rerun/
#PBS -j oe

# -J and the resources are rewritten by se_check.R (check_failed()) to match
# failed_ids.txt - run that before submitting this.

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

jobid=$(sed -n "${PBS_ARRAY_INDEX}p" "${PBS_O_WORKDIR}/failed_ids.txt")

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript se_analysis.R "$jobid"
