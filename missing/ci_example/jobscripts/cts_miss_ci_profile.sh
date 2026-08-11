#!/bin/bash
# Deliberately over-provisioned so nothing fails for the wrong reason. ncpus=20
# covers the largest inner configuration (workers1=3 x workers2=3 x
# grf_threads=2 = 18 cores) with headroom for the level-1 "requester"
# processes that stay alive while their own level-2 pool runs (see
# cts_miss_ci_profile.R's expected_procs comment). ompthreads=20 matches
# ncpus so parallelly::availableCores() sees all of them and
# plan(list(tweak(multisession, workers=3), ...)) doesn't fail; the profile
# script still sets OMP_NUM_THREADS from the sweep grid before spawning its
# workers, so per-worker thread counts are controlled there, not by
# ompthreads.
#
# #PBS -J 1-10, not the 500-row production grid: cts_miss_ci_profile.R uses a
# two-level outer(scenario x run)/inner(workers1 x workers2 x grf_threads x
# CI_boot) grid, not a flat one - one array index covers the OUTER 10-row
# grid and loops all 16 inner cells itself, reusing the one
# mice(m=50)-imputed dataset across them.
#
# Walltime is generous and largely a guess: unlike missing/continuous's
# cts_miss_profile.R (one run_all_cate_methods() call per inner cell), each
# cell here is up to 50 imputations x 4 CATE arms x 60 bootstrap draws, and
# the cheapest cell (workers1=1, workers2=1) has no parallelism to hide
# behind at all - if 36:00:00 proves too tight once real timings come in,
# widen it and resubmit rather than trust this number.
#
# Memory and CPU are measured inside R by syrup, so there is no background
# qstat sampler here and nothing depends on the scheduler.
#PBS -l walltime=36:00:00
#PBS -l select=1:ncpus=20:ompthreads=20:mem=40gb
#PBS -J 1-10
#PBS -N cts_miss_ci_profile
#PBS -o logs_ci_profile/
#PBS -e logs_ci_profile/

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

Rscript cts_miss_ci_profile.R "$PBS_ARRAY_INDEX"
