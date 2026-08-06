#!/bin/bash
# Pilot: half-sample bootstrap CIs for the 5 crossfit-structured RF arms.
# 3 scenarios x 50 runs = 150 array jobs. PLACEHOLDER resources - the
# cf_profile.R/syrup sweep was deliberately skipped for this pilot (built for
# tuning a 2000-job production run, not worth it for 150 jobs). Before
# submitting the full array, run a small manual batch first (e.g.
# `qsub -J 1-5` below, or a couple of local `Rscript` timings) and check
# `qstat -fx <jobid> | grep resources_used`, then hand-set walltime/mem to
# ~2x what was observed - the bootstrap refits (~200 x V forests per arm x 5
# arms) dominate cost here, well beyond what cf_1.sh's point-estimate-only
# timings would suggest.
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=4gb
#PBS -J 1-150
#PBS -N cf_ci_1
#PBS -o logs_ci_1/
#PBS -e logs_ci_1/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters. Trailing args are CI_boot/workers/grf_threads
# - CI_boot=200 is the pilot's target; leave workers/grf_threads at 2 1 until
# resources above are confirmed from a small manual batch.
Rscript cf_ci_analysis.R "$PBS_ARRAY_INDEX" 200 2 1
