#!/bin/bash
# Scenarios 1-6, n = 4000, hte_mult 1 (snr_3.sh is scenario 7). NOT run to completion
# locally: rfsrc died with an allocation failure at about 9.4gb on a 15gb laptop, after 8.6 min. Memory went
# 1.3 / 2.7 / 8.0gb at n = 500 / 1000 / 2000, so expect roughly 25-30gb here;
# 32gb and 2h are headroom, not a measurement. Submit two first
# (qsub -J 3701-3702 snr_6.sh), read their .o files, then qsub -J 3703-4300 snr_6.sh.
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=2:ompthreads=2:mem=32gb
#PBS -J 3701-4300
#PBS -N snr_6
#PBS -o logs_6/
#PBS -j oe


module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to script directory
cd "${PBS_O_WORKDIR}/.."

# Run R script with parameters
Rscript snr_analysis.R "$PBS_ARRAY_INDEX"
