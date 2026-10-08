#!/bin/bash
# Reruns only the array indices listed in failed_ids_strat.txt.
#
# NOTE: me_check.R deliberately does NOT write that file automatically for
# this tree (it runs check_failed with write=FALSE) - study_strat's grid is
# all 600 rows, so any runs excluded from the main study would be reported
# missing forever and queued here forever. `Rscript me_check.R strategies`
# prints the indices that genuinely need resubmitting; put those in
# failed_ids_strat.txt by hand and set -J to match - keep the %N at 4 or
# below, as in me_strategies.sh.
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=10:ompthreads=1:mem=16gb
#PBS -J 1-2%1
#PBS -N me_strat_rerun
#PBS -o logs_strat_rerun/
#PBS -e logs_strat_rerun/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

jobid=$(sed -n "${PBS_ARRAY_INDEX}p" "${PBS_O_WORKDIR}/failed_ids_strat.txt")

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

Rscript me_strategies.R "$jobid" 10
