#!/bin/bash
# Pilot: half-sample bootstrap CIs for the 5 crossfit-structured RF arms.
# 3 scenarios x 50 runs = 150 array jobs. Resources below are set by
# cf_ci_profile.R / cf_ci_profile_summary.R (the same syrup sweep + writeback
# cf_1.sh uses, adapted for the CI pilot's CI_boot/workers/grf_threads trailing
# args - see crossfitting/README.md's "Sizing the CI pilot's array job"). If
# they still read as placeholders below, the sweep hasn't been run yet: do that
# before submitting the full array, rather than hand-guessing - the bootstrap
# refits (~200 x V forests per arm x 5 arms) dominate cost here, well beyond
# what cf_1.sh's point-estimate-only timings would suggest.
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

# Run R script with parameters. Trailing args are CI_boot/workers/grf_threads,
# overwritten by cf_ci_profile_summary.R alongside the #PBS -l lines above - kept
# in sync since both are written by the same sub() call there. Below is still the
# pre-profiling placeholder (workers=2 here vs. ncpus=1 above): do not hand-edit
# one without the other, run the sweep instead.
Rscript cf_ci_analysis.R "$PBS_ARRAY_INDEX" 200 2 1
