#!/bin/bash
#SBATCH --job-name=instability_run
#SBATCH --output=instability_run_%j.out
#SBATCH --error=instability_run_%j.err
#SBATCH --time=01:00:00           # Time limit hrs:min:sec
#SBATCH --ntasks=1                # Number of tasks (1)
#SBATCH --cpus-per-task=4         # Number of CPU cores per task (adjust as needed)
#SBATCH --mem=32G                  # Memory per node (adjust as needed)

# Activate conda environment
source /home/jboktor/.bashrc
source activate spatialomics

# Read command-line arguments
FOLDER_PATH=$1
K=$2
B=$4
SHAPE_FLAT=$5
NAME_DG=$6
OUTPUT_FILE=$7

# Run the Python script
srun python /central/groups/mthomson/jboktor/spatial_genomics/jess_2024-01-23/notebooks/python_scripts/run_instability.py $FOLDER_PATH $K $K $B $SHAPE_FLAT $NAME_DG $OUTPUT_FILE
