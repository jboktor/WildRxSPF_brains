#!/bin/bash
# Submit with:
#   sbatch notebooks/shell_scripts/run_cell_neighborhood_group_tests_v2.sh

#SBATCH --job-name="cell-neighborhoods-group-tests"
#SBATCH --partition=expansion
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/cell_neighborhoods_mapmycells/log_group_tests_%x_%j.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/cell_neighborhoods_mapmycells/log_group_tests_%x_%j.err"
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END

set -euo pipefail

WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
RSCRIPT="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript"
SCRIPT="${WKDIR}/notebooks/R_scripts/run_cell_neighborhood_group_tests_v2.R"

log() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"
}

log "Starting exploratory ROI-level SPF-WildR neighborhood tests"
log "Host: $(hostname)"
log "Rscript: ${RSCRIPT}"

cd "${WKDIR}"
"${RSCRIPT}" "${SCRIPT}" run

log "Exploratory group tests completed"
