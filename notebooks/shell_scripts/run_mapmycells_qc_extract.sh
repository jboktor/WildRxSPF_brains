#!/bin/bash
# Flatten MapMyCells extended JSON -> per-cell mapping-quality parquet.
# The SN JSON is 1.25 GB (~25 GB RAM to parse), hence SLURM.
#
# Submit:   sbatch notebooks/shell_scripts/run_mapmycells_qc_extract.sh          # SN (default)
# Override: sbatch --export=ALL,MODALITY=seqFISH notebooks/shell_scripts/run_mapmycells_qc_extract.sh
# Output:   data/interim/MapMyCells_QC/mapmycells_scores_<MODALITY>_<date>.parquet
#SBATCH --job-name=mmc_qc_extract
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --partition=expansion
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/mapmycells_qc/log_%x_%j.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/mapmycells_qc/log_%x_%j.err"
set -euo pipefail

WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
PY="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python"
MODALITY="${MODALITY:-SN}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

log "start MODALITY=${MODALITY} host=$(hostname)"
"${PY}" "${WKDIR}/notebooks/python_scripts/extract_mapmycells_scores.py" --modality "${MODALITY}"
log "done"
