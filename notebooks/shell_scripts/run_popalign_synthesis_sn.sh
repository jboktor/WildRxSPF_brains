#!/bin/bash
#SBATCH --job-name="popalign_synth_sn"
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/popalign/log_popalign_synthesis_sn_%j.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/popalign/log_popalign_synthesis_sn_%j.err"
#SBATCH --time=0-10:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=expansion
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END
#
# Aggregate per-unit popAlign results → animal-level stats, annotation, catalog, figures.
# Submit after the array finishes:
#   sbatch --dependency=afterany:<array_jobid> notebooks/shell_scripts/run_popalign_synthesis_sn.sh --run-tag 2026-10-04
# (afterany: synthesise whatever finished; failed units are listed in the status table)

set -euo pipefail
WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
RSCRIPT="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript"
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"; }
cd "${WKDIR}"
log "popAlign synthesis"
"${RSCRIPT}" notebooks/R_scripts/run_popalign_synthesis_sn.R "$@"
log "Done"
