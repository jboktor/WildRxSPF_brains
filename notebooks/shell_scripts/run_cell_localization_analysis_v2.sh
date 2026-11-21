#!/bin/bash
# Submit with:
#   sbatch notebooks/shell_scripts/run_cell_localization_analysis_v2.sh
# Submit the 16 larger whole-slice jobs separately:
#   sbatch --array=146-161%4 --time=08:00:00 --cpus-per-task=4 --mem=96G \
#     notebooks/shell_scripts/run_cell_localization_analysis_v2.sh

#SBATCH --job-name="cell-neighborhoods-v2"
#SBATCH --partition=expansion
#SBATCH --array=1-145%24
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/cell_neighborhoods_mapmycells/log_cell_neighborhoods_%A_%a.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/cell_neighborhoods_mapmycells/log_cell_neighborhoods_%A_%a.err"
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END

set -euo pipefail

WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
RSCRIPT="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript"
SCRIPT="${WKDIR}/notebooks/R_scripts/run_cell_localization_analysis_v2.R"

log() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"
}

log "Starting array job ${SLURM_ARRAY_JOB_ID}, task ${SLURM_ARRAY_TASK_ID}"
log "Host: $(hostname)"
log "Rscript: ${RSCRIPT}"

cd "${WKDIR}"
"${RSCRIPT}" "${SCRIPT}" run-job "${SLURM_ARRAY_TASK_ID}"

log "Completed task ${SLURM_ARRAY_TASK_ID}"
