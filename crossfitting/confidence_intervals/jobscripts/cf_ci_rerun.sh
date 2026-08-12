#!/bin/bash
# Reruns only the array indices listed in failed_ids.txt by cf_ci_check.R.
# -J and the resource request are rewritten by check_failed(); do not hand-edit.
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=10gb
#PBS -J 1-53
#PBS -N cf_ci_rerun
#PBS -o logs_ci_rerun/
#PBS -e logs_ci_rerun/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

jobid=$(sed -n "${PBS_ARRAY_INDEX}p" "${PBS_O_WORKDIR}/failed_ids.txt")

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

Rscript cf_ci_analysis.R "$jobid" 200 2 1
