#!/bin/bash
#SBATCH --job-name="dspin_export_full"
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/dspin/log_dspin_export_full_%j.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/dspin/log_dspin_export_full_%j.err"
#SBATCH --time=0-18:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --partition=expansion
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END
#
# Rank all hard_pass (tissue × supertype_label) units + export raw_counts.h5ad.
# Paths: data/interim/dspin/{tissue}__{supertype_label}/  (snRNASeq only).
# Then HVG-prep every unit missing filtered.h5ad.

set -euo pipefail

WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
RSCRIPT="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript"
PY="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python"

log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"; }

mkdir -p "${WKDIR}/.cluster_runs/dspin" "${WKDIR}/data/interim/dspin"
cd "${WKDIR}"

log "=== rank → full_manifest.csv ==="
"${RSCRIPT}" notebooks/R_scripts/rank_dspin_pilots.R

log "=== export raw h5ads (tier=full; skips existing) ==="
"${RSCRIPT}" notebooks/R_scripts/export_dspin_h5ads.R \
  --mode export-pilot \
  --tier full \
  --manifest "${WKDIR}/data/interim/dspin/full_manifest.csv"

log "=== HVG prep all slices ==="
"${PY}" notebooks/python_scripts/dspin_prep_filtered.py \
  --manifest "${WKDIR}/data/interim/dspin/full_manifest.csv"

log "=== Done export+prep ==="
n_raw=$(ls -1 data/interim/dspin/*/raw_counts.h5ad 2>/dev/null | wc -l)
n_filt=$(ls -1 data/interim/dspin/*/filtered.h5ad 2>/dev/null | wc -l)
log "raw_counts.h5ad=${n_raw} filtered.h5ad=${n_filt}"
