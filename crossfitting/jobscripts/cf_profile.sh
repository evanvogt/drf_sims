#!/bin/bash
# Deliberately over-provisioned so nothing fails for the wrong reason. ncpus=8
# covers the largest configuration in the sweep (workers=4 x grf_threads=2).
# ompthreads is left off the select line: cf_profile.R sets OMP_NUM_THREADS from
# the sweep grid before spawning its multisession workers.
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -J 1-36
#PBS -N cf_profile
#PBS -o logs_profile/
#PBS -e logs_profile/

mkdir -p "${PBS_O_WORKDIR}/logs_profile"

# background sampler for true peak memory - the multisession workers are separate
# processes, so R's own gc() cannot see them but PBS accounting can
USAGE_LOG="${PBS_O_WORKDIR}/logs_profile/usage_${PBS_ARRAY_INDEX}.txt"
(
  while true; do
    qstat -f "$PBS_JOBID" 2>/dev/null | grep -E 'resources_used\.(mem|vmem|cput|walltime)'
    sleep 30
  done
) > "$USAGE_LOG" 2>&1 &
SAMPLER=$!

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

Rscript cf_profile.R "$PBS_ARRAY_INDEX"
STATUS=$?

kill $SAMPLER 2>/dev/null

# final accounting
qstat -fx "$PBS_JOBID" 2>/dev/null | grep resources_used >> "$USAGE_LOG"

exit $STATUS
