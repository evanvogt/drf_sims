#PBS -l walltime=02:00:00  
#PBS -l select=1:ncpus=20:ompthreads=20:mem=15gb
#PBS -J 1-10000
#PBS -N cts_500
#PBS -o logs_500/
#PBS -e logs_500/

module purge
module add tools/prod
module add R/4.3.2-gfbf-2023a

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate sim-env

# scenario and sample sze arrays
n=500
scenarios=(1 2 3 4 5 6 7 8 9 10) 

# convert to 0-based for counting along arrays
id=$((PBS_ARRAY_INDEX - 1))

# calculate parameters
scen_id=$((id/1000))
sim_id=$((id % 1000 + 1))

# Get actual values from arrays
scenario=${scenarios[$scen_id]}

# Display the mapping
echo "PBS_ARRAY_INDEX: $PBS_ARRAY_INDEX"
echo "Sample Size: $n"
echo "Scenario: $scenario" 
echo "Simulation Run: $sim_id"

# Run script
cd "${PBS_O_WORKDIR}/.."

Rscript cts_analysis.R "$scenario" "$n" "$sim_id"
