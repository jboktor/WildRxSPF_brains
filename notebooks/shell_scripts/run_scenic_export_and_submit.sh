#!/bin/bash
#SBATCH --job-name="scenic_export"
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/scenic/log_scenic_export_%j.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/scenic/log_scenic_export_%j.err"
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --partition=expansion
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END

# 1) Export raw-count looms for class splits (tissue-wide track discontinued —
#    whole-tissue GRN runs did not finish within 7-day walltime; see CLAS_01
#    subsampling approach instead for large splits)
# 2) sbatch one pySCENIC job per class loom

set -euo pipefail

WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
RSCRIPT="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript"
EXPORT="${WKDIR}/notebooks/R_scripts/export_scenic_looms.R"
DYNAMIC="${WKDIR}/notebooks/shell_scripts/run_pyscenic_dynamic.sh"
INPUT_DIR="${WKDIR}/data/input/scenic/snRNASeq"
SUBMIT_LIST="${WKDIR}/.cluster_runs/scenic/submitted_jobs.txt"

log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"; }

mkdir -p "${WKDIR}/.cluster_runs/scenic"
cd "${WKDIR}"

log "=== Export looms (mode=class; RNA counts = raw UMIs) ==="
"${RSCRIPT}" "${EXPORT}" --mode class

log "=== Submit class-level pySCENIC jobs ==="
: > "${SUBMIT_LIST}"
shopt -s nullglob
class_looms=("${INPUT_DIR}"/snRNASeq_*_CLAS_*.loom)
log "Found ${#class_looms[@]} class looms"
for loom in "${class_looms[@]}"; do
  base=$(basename "${loom}")
  stem=${base%.loom}
  meta="${INPUT_DIR}/${stem}_meta.csv"
  n_cells=0
  if [[ -f "${meta}" ]]; then
    # header + N rows
    n_cells=$(($(wc -l < "${meta}") - 1))
  fi
  # Scale resources by cell count (GRNBoost2 memory grows with n)
  mem="300G"; time_lim="3-00:00:00"; cpus=24
  if (( n_cells >= 50000 )); then
    mem="500G"; time_lim="5-00:00:00"; cpus=32
  elif (( n_cells >= 20000 )); then
    mem="400G"; time_lim="4-00:00:00"; cpus=28
  fi
  log "Submit ${base} n_cells=${n_cells} mem=${mem} cpus=${cpus} time=${time_lim}"
  jid=$(sbatch --parsable \
    --job-name="scenic_${stem}" \
    --mem="${mem}" --time="${time_lim}" --cpus-per-task="${cpus}" \
    "${DYNAMIC}" -i "${INPUT_DIR}" -t "${base}")
  echo "CLASS ${jid} ${base} n=${n_cells} mem=${mem}" | tee -a "${SUBMIT_LIST}"
done

log "=== Submitted job list written to ${SUBMIT_LIST} ==="
cat "${SUBMIT_LIST}"
log "DONE export + submit"
