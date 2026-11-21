#!/bin/bash
#SBATCH --job-name="popalign_sn"
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/popalign/log_popalign_sn_%A_%a.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/popalign/log_popalign_sn_%A_%a.err"
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --partition=expansion
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --array=0-243%40
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END
#
# One popAlign run per manifest row. Run 2026-10-04 (244 units):
#   0-20 = class, 21-104 = subclass, 105-243 = supertype.
# Arg 1 = manifest (default: run 2026-10-04). Set --array to 0-(rows-1).
# Submit AFTER the export job:
#   sbatch --array=0-243%40 --dependency=afterok:<export_jobid> \
#     notebooks/shell_scripts/run_popalign_array_sn.sh data/interim/popalign/snRNASeq/manifest_2026-10-04.csv
# Re-run only failed / missing units (finished units are skipped via status.json):
#   sbatch --array=<idx,idx,...> notebooks/shell_scripts/run_popalign_array_sn.sh <manifest>

set -euo pipefail
WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
PY="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/popalign/bin/python"
export PYTHONPATH="/resnick/groups/mthomson/jboktor/spatial_genomics/popalign"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1   # popAlign parallelises with Pool
export MPLBACKEND=Agg
MANIFEST="${1:-${WKDIR}/data/interim/popalign/snRNASeq/manifest_2026-10-04.csv}"
[[ "${MANIFEST}" = /* ]] || MANIFEST="${WKDIR}/${MANIFEST}"
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"; }

log "Task ${SLURM_ARRAY_TASK_ID} on $(hostname) | manifest ${MANIFEST}"
"${PY}" "${WKDIR}/notebooks/python_scripts/popalign_run_unit_sn.py" \
  --manifest "${MANIFEST}" --index "${SLURM_ARRAY_TASK_ID}" --skip-done
log "Done"
