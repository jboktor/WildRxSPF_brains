#!/bin/bash
#SBATCH --job-name="integrated_sn"
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/sn_integrated/log_integrated_sn_%j.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/sn_integrated/log_integrated_sn_%j.err"
#SBATCH --time=0-08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=expansion
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END
#
# Integrated MAST × pySCENIC × D-SPIN × popAlign analysis (stages gene → cell → animal → circuit → evidence).
# Reads saved outputs only; no upstream method is re-run.
# Submit:  sbatch notebooks/shell_scripts/run_integrated_analysis_sn.sh [stage] [run-tag]
#          stage defaults to "all" minus figures (figures are drawn separately, directly).

set -euo pipefail
WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
RSCRIPT="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript"
RUN_TAG="${2:-2026-10-06}"
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"; }
cd "${WKDIR}"
for st in ${1:-gene cell animal circuit evidence}; do
  log "Stage ${st}"
  "${RSCRIPT}" notebooks/R_scripts/run_integrated_analysis_sn.R --stage "${st}" --run-tag "${RUN_TAG}"
done
log "Done"
