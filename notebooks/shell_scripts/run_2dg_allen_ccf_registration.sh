#!/bin/bash
#
# Register one prepared 2DG Template slice to Allen CCF.
#
# Primary inputs are the scalar stacks under data/input/2DG_data/Template and
# data/input/2DG_data/WildRvSPF Intruder (not the older RGB overlay exports).
#
# Prepare once on the login node (quick, no registration):
#   /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python \
#     notebooks/python_scripts/register_2dg_to_allen_ccf.py prepare
#
# Submit all 66 slices:
#   mkdir -p /resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/2dg_allen_ccf
#   sbatch notebooks/shell_scripts/run_2dg_allen_ccf_registration.sh
#
# After every array task succeeds, aggregate explicitly:
#   /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python \
#     notebooks/python_scripts/register_2dg_to_allen_ccf.py aggregate
#   /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python \
#     notebooks/python_scripts/register_2dg_to_allen_ccf.py export-scalar-t
#
# Override one task for troubleshooting:
#   sbatch --array=30 notebooks/shell_scripts/run_2dg_allen_ccf_registration.sh

#SBATCH --job-name="2dg_ccf"
#SBATCH --partition=expansion
#SBATCH --array=0-65%12
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=12G
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/2dg_allen_ccf/log_2dg_ccf_%A_%a.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/2dg_allen_ccf/log_2dg_ccf_%A_%a.err"
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END

set -euo pipefail

WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
PYTHON="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python"
SCRIPT="${WKDIR}/notebooks/python_scripts/register_2dg_to_allen_ccf.py"
MANIFEST="${WKDIR}/data/interim/2DG_CCF_Integration/2026-11-23c/2dg_slice_manifest_2026-11-23c.csv"

log() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"
}

if [[ ! -f "${MANIFEST}" ]]; then
  log "ERROR: preparation manifest is missing: ${MANIFEST}"
  log "Run the documented prepare command before submitting the array."
  exit 1
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  log "ERROR: SLURM_ARRAY_TASK_ID is not set; submit this file with sbatch."
  exit 1
fi

cd "${WKDIR}"
log "Starting 2DG slice position ${SLURM_ARRAY_TASK_ID}"
log "Python: ${PYTHON}"
log "Script: ${SCRIPT}"

"${PYTHON}" "${SCRIPT}" register-slice \
  --slice-position "${SLURM_ARRAY_TASK_ID}"

log "Completed 2DG slice position ${SLURM_ARRAY_TASK_ID}"
