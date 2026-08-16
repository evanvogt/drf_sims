#!/bin/bash
# Reruns only the array indices listed in failed_ids.txt by me_check.R.
# -J and the resource request are rewritten by check_failed(); do not hand-edit.
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=3:ompthreads=2:mem=12gb
#PBS -J 1-3%2
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
