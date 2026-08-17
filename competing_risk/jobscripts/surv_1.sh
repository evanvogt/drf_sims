#!/bin/bash
# walltime trimmed 02:00:00 -> 01:00:00 after the crossfitting change removed the
# C(V,2) fold-pair loops. Hand-set from local timings, NOT a cf_profile.R-style
# sweep (this study has no profile script): three smoke-test replicates at
# n=500, V=10, 2 workers took 9.3, 9.9 and 10.7 min on R 4.5.3. One hour leaves
# ~6x headroom, which covers the cluster being slower and on R 4.3.2.
# mem=2gb. This comment used to say "left at 10gb", which the select= line below
# never matched. 2gb is now the measured answer rather than a guess: 1175 of the
# 1400 array jobs completed at this request, and the 225 that did not were an R
# error in pseudo_sl_t_split (see competing_risk/surv_failed_diagnose.R), not an
# OOM - they failed again unchanged at 20gb.
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=2gb
#PBS -J 1-1400%100
#PBS -N surv_1
#PBS -o logs_1/
#PBS -e logs_1/


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript surv_analysis.R "$PBS_ARRAY_INDEX"