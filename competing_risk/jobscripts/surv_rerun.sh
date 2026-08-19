#!/bin/bash
# Reruns only the array indices listed in failed_ids.txt by surv_check.R.
# -J and the resource request are rewritten by check_failed(); the values here
# are what it computes from surv_1.sh, so the first check leaves them alone.
# Needs logs_rerun/ to exist on the cluster - PBS aborts a job whose output
# directory is missing, and *logs*/ is gitignored so it is never checked out.
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=3:ompthreads=2:mem=20gb
#PBS -J 1-19%100
#PBS -N surv_rerun
#PBS -o logs_rerun/
#PBS -e logs_rerun/

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

jobid=$(sed -n "${PBS_ARRAY_INDEX}p" "${PBS_O_WORKDIR}/failed_ids.txt")

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript surv_analysis.R "$jobid"
