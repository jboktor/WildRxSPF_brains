#!/bin/bash
#SBATCH --job-name="dspin_export"
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/dspin/log_dspin_export_%j.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/dspin/log_dspin_export_%j.err"
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --partition=expansion
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END
#
# Modes via first arg:
#   candidates   (default) — cell counts + meta from Seurat
#   rank         — join MAST/Augur → pilot_manifest.csv
#   export-pilot — write raw_counts.h5ad for interactive pilot slices
#   all          — candidates → rank → export-pilot

set -euo pipefail

WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
RSCRIPT="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript"
MODE="${1:-all}"

log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"; }

mkdir -p "${WKDIR}/.cluster_runs/dspin" "${WKDIR}/data/interim/dspin"
cd "${WKDIR}"

case "${MODE}" in
  candidates)
    log "=== candidates ==="
    "${RSCRIPT}" notebooks/R_scripts/export_dspin_h5ads.R --mode candidates
    ;;
  rank)
    log "=== rank pilots ==="
    "${RSCRIPT}" notebooks/R_scripts/rank_dspin_pilots.R
    ;;
  export-pilot)
    log "=== export-pilot (interactive) ==="
    "${RSCRIPT}" notebooks/R_scripts/export_dspin_h5ads.R --mode export-pilot --tier interactive
    ;;
  all)
    log "=== candidates ==="
    "${RSCRIPT}" notebooks/R_scripts/export_dspin_h5ads.R --mode candidates
    log "=== rank pilots ==="
    "${RSCRIPT}" notebooks/R_scripts/rank_dspin_pilots.R
    log "=== export-pilot (interactive) ==="
    "${RSCRIPT}" notebooks/R_scripts/export_dspin_h5ads.R --mode export-pilot --tier interactive
    ;;
  *)
    echo "Unknown mode: ${MODE}" >&2
    exit 1
    ;;
esac

log "=== Done (mode=${MODE}) ==="
