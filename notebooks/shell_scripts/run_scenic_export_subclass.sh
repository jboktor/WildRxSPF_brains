#!/bin/bash
#SBATCH --job-name="scenic_export_split"
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/scenic/log_scenic_export_subclass_%x_%j.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/scenic/log_scenic_export_subclass_%x_%j.err"
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=expansion
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END

# Slice raw-UMI subclass looms out of the class looms (no Seurat load).
# Writes data/interim/scenic/looms_subclass/{*.loom,*_meta.csv,manifest_subclass.csv,
#   looms_subclass_{small,large}.txt}; existing looms are skipped, not rewritten.
#
# Submit:
#   sbatch notebooks/shell_scripts/run_scenic_export_subclass.sh                         # subclass
#   sbatch --export=ALL,LEVEL=supertype notebooks/shell_scripts/run_scenic_export_subclass.sh  # supertype
# Then submit the GRN arrays (run_pyscenic_subclass_array.sh header).

set -euo pipefail

WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
PY="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python"

log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"; }

cd "${WKDIR}"
LEVEL="${LEVEL:-subclass}"   # subclass | supertype (sbatch --export=ALL,LEVEL=supertype)
log "Export ${LEVEL} looms"
"${PY}" "${WKDIR}/notebooks/python_scripts/export_scenic_subclass_looms.py" --level "${LEVEL}" --max-cells 45000 --large-n 5000
log "DONE"
