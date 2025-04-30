#!/bin/bash
#SBATCH --job-name=nmf_parallel
#SBATCH --output=nmf_parallel_%j.out
#SBATCH --error=nmf_parallel_%j.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --mem=64G

module load python/3.9  # Adjust the module load command as per your environment

# Activate your virtual environment if necessary
# source /path/to/your/virtualenv/bin/activate

# Run your Python script
srun python nmf_parallel.py
