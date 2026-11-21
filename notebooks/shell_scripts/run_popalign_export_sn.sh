#!/bin/bash
#SBATCH --job-name="popalign_export_sn"
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/popalign/log_popalign_export_sn_%j.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/popalign/log_popalign_export_sn_%j.err"
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=400G
#SBATCH --partition=expansion
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END
#
# Load the snRNA Seurat object once and write popAlign mtx inputs for every unit
# in data/interim/popalign/snRNASeq/manifest_<run-tag>.csv (class, subclass,
# supertype). Units whose inputs already exist are skipped (inputs are shared
# across run tags).
# Build the manifest first (light, no Seurat):
#   Rscript notebooks/R_scripts/export_popalign_inputs_sn.R --mode manifest --run-tag 2026-10-04
# Submit (arg 1 = run tag, default 2026-10-04):
#   sbatch notebooks/shell_scripts/run_popalign_export_sn.sh 2026-10-04

set -euo pipefail
WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
RSCRIPT="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript"
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"; }

cd "${WKDIR}"
log "Exporting popAlign inputs"
RUN_TAG="${1:-2026-10-04}"
"${RSCRIPT}" notebooks/R_scripts/export_popalign_inputs_sn.R --mode export --run-tag "${RUN_TAG}"
log "Done"
