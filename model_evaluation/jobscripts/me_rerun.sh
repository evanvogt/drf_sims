#!/bin/bash
# Reruns only the array indices listed in failed_ids.txt by me_check.R.
# Set -J to 1-<number of lines in failed_ids.txt> before submitting.
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=2gb
#PBS -J 1-1
#PBS -N me_rerun
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

# Trailing args are workers/n_cores - keep these in sync with whatever
# me_profile_summary.R last wrote into me_1.sh.
Rscript me_analysis.R "$jobid" 1 1
