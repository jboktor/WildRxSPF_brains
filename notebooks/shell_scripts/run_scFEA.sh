#!/bin/bash
#SBATCH --output="/resnick/groups/MazmanianLab/jboktor/WILDRxSPF_brains/.cluster_runs/scFEA/log_scFEA_%x_%j.log"
#SBATCH --job-name="scFEA"
#SBATCH --time=2-00:00:00     # walltime
#SBATCH --cpus-per-task=8  # number of cores
#SBATCH --partition=gpu
#SBATCH --gres=gpu:nvidia_h200:1 # number of GPUs
#SBATCH --ntasks=1
#SBATCH --mem-per-gpu=40G

log() {
  # Timestamped log lines (go to both Slurm log + our tee log)
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"
}

# load CUDA module
module load cuda/12.2.1-gcc-11.3.1-sdqrj2e
log "Loaded CUDA module"

# Activate conda environment
source "/home/${USER}/.bashrc"
source activate scFEA
log "Activated conda env (scFEA)"


scfea_gitpath="/resnick/groups/MazmanianLab/jboktor/git/scFEA"
scfea_input_path="/resnick/groups/MazmanianLab/jboktor/WILDRxSPF_brains/data/input/scFEA/snRNASeq"
scfea_output_path="/resnick/groups/MazmanianLab/jboktor/WILDRxSPF_brains/data/interim/scFEA/snRNASeq"

log "Running scFEA"

cd ${scfea_gitpath}
log "CWD after cd: $(pwd)"

test_file="snRNASeq_log2cpm_HYP.csv"
input_name="${test_file%.*}"
log "Test file: ${test_file}"
log "Input name: ${input_name}"

python src/scFEA.py \
--data_dir ${scfea_gitpath}/data \
--input_dir "${scfea_input_path}" \
--test_file "${test_file}" \
--moduleGene_file module_gene_complete_mouse_m168.csv \
--stoichiometry_matrix cmMat_complete_mouse_c70_m168.csv \
--sc_imputation True \
--output_flux_file "${scfea_output_path}/flux_${input_name}.csv" \
--output_balance_file "${scfea_output_path}/balance_${input_name}.csv"

log "scFEA completed!"