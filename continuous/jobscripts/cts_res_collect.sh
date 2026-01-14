#PBS -l walltime=01:00:00  
#PBS -l select=1:ncpus=10:ompthreads=10:mem=15gb
#PBS -N cts_results_collect
#PBS -o logs_collect/
#PBS -e logs_collect/

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# Navigate to the script directory
cd "${PBS_O_WORKDIR}/.."

# Run the R script for the assigned scenario and sample size
Rscript collect_results_cts.R
