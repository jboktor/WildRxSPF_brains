#!/bin/bash
#SBATCH --output="/resnick/groups/MazmanianLab/jboktor/WILDRxSPF_brains/.cluster_runs/scFEA/log_scFEA_%x_%j.log"
#SBATCH --job-name="scFEA"
#SBATCH --time=7-00:00:00     # walltime
#SBATCH --cpus-per-task=24  # number of cores
#SBATCH --partition=gpu
#SBATCH --gres=gpu:nvidia_h200:1 # number of GPUs
#SBATCH --ntasks=1
#SBATCH --mem=600G

while getopts "i:o:t:" opt; do
  case $opt in
    i) scfea_input_path="$OPTARG";;
    o) scfea_output_path="$OPTARG";;
    t) input_file="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&2; exit 1;;
  esac
done


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
# scfea_input_path="/resnick/groups/MazmanianLab/jboktor/WILDRxSPF_brains/data/input/scFEA/snRNASeq"
# scfea_output_path="/resnick/groups/MazmanianLab/jboktor/WILDRxSPF_brains/data/interim/scFEA/snRNASeq"

log "Running scFEA"

cd ${scfea_gitpath}
log "CWD after cd: $(pwd)"

# HERE IS WHERE YOU CHANGE THE FILE YOU ARE RUNNING
# test_file="snRNASeq_log2cpm_AMY.csv"

log "Test file: ${input_file}"
input_name="${input_file%.*}"
log "Input name: ${input_name}"

export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

python src/scFEA.py \
--data_dir ${scfea_gitpath}/data \
--input_dir "${scfea_input_path}" \
--test_file "${input_file}" \
--moduleGene_file module_gene_complete_mouse_m168.csv \
--stoichiometry_matrix cmMat_complete_mouse_c70_m168.csv \
--sc_imputation True \
--output_flux_file "${scfea_output_path}/flux_${input_name}_run3.csv" \
--output_balance_file "${scfea_output_path}/balance_${input_name}_run3.csv"

log "scFEA completed!"