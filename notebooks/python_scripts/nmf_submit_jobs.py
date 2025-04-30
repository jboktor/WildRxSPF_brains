import os
import subprocess
import numpy as np

K1 = 2
K2 = 4
numPatterns = np.arange(K1, K2+1)  # Example values for number of patterns
folder_path = '/central/groups/mthomson/jboktor/spatial_genomics/jess_2024-01-23/data/interim/osNMF/instability_runs'  # Define your folder path
B = 100  # Number of bootstrap replicates

slurm_template = """#!/bin/bash
#SBATCH --job-name=nmf_K{K}_b{b}
#SBATCH --output=nmf_K{K}_b{b}_%j.out
#SBATCH --error=nmf_K{K}_b{b}_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem=32G

# Activate your virtual environment if necessary
source /home/jboktor/.bashrc
source activate spatialomics

# Run your Python script
python /central/groups/mthomson/jboktor/spatial_genomics/jess_2024-01-23/notebooks/python_scripts/nmf_parallel.py {K} {b} {folder_path}
"""

for K in numPatterns:
    for b in range(B):
        # Check if the expected output file already exists
        expected_output_file = f"{folder_path}/K={K}/DG_PP_{b}.npz"
        if os.path.exists(expected_output_file):
            continue
        
        slurm_script = slurm_template.format(K=K, b=b, folder_path=folder_path)
        script_filename = f"nmf_slurm_K{K}_b{b}.sh"
        
        with open(script_filename, 'w') as f:
            f.write(slurm_script)
        
        # Submit the job
        subprocess.run(["sbatch", script_filename])
