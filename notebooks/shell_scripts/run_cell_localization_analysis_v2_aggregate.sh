#!/bin/bash
# Submit only after the array job succeeds:
#   sbatch --dependency=afterok:<array_job_id> \
#     notebooks/shell_scripts/run_cell_localization_analysis_v2_aggregate.sh

#SBATCH --job-name="cell-neighborhoods-v2-aggregate"
#SBATCH --partition=expansion
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=96G
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/cell_neighborhoods_mapmycells/log_cell_neighborhoods_%x_%j.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/cell_neighborhoods_mapmycells/log_cell_neighborhoods_%x_%j.err"
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END

set -euo pipefail

WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
RSCRIPT="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript"
SCRIPT="${WKDIR}/notebooks/R_scripts/run_cell_localization_analysis_v2.R"

log() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"
}

log "Starting aggregation and primary visualizations"
log "Host: $(hostname)"
log "Rscript: ${RSCRIPT}"

cd "${WKDIR}"
"${RSCRIPT}" "${SCRIPT}" aggregate
"${RSCRIPT}" "${SCRIPT}" plot-confidence
"${RSCRIPT}" "${SCRIPT}" plot-whole-slice
"${RSCRIPT}" "${SCRIPT}" plot-manual-gates

log "Aggregation and primary visualizations completed"
