#!/bin/bash
#SBATCH --job-name="pyscenic"
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/scenic/log_pyscenic_%x_%j.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/scenic/log_pyscenic_%x_%j.err"
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=300G
#SBATCH --partition=expansion
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL

# Unit job: one loom → pySCENIC grn/ctx/aucell (+ WildR vs SPF diff).
# Override resources on submit for tissue-scale looms, e.g.:
#   sbatch --mem=600G --time=7-00:00:00 --cpus-per-task=32 run_pyscenic_dynamic.sh -t snRNASeq_AMY.loom
#
# Usage:
#   sbatch run_pyscenic_dynamic.sh -t snRNASeq_AMY_CS20230722_CLAS_01.loom
#   sbatch run_pyscenic_dynamic.sh -i /path/to/looms -t snRNASeq_AMY.loom

while getopts "i:t:" opt; do
  case $opt in
    i) scenic_input_path="$OPTARG";;
    t) input_file="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&2; exit 1;;
  esac
done

WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
RSCRIPT="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript"
RUNNER="${WKDIR}/notebooks/R_scripts/run_scenic_one_split.R"

scenic_input_path="${scenic_input_path:-${WKDIR}/data/input/scenic/snRNASeq}"

log() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"
}

if [[ -z "${input_file:-}" ]]; then
  log "ERROR: pass -t <loom_basename.loom>"
  exit 1
fi

if [[ "${input_file}" = /* ]]; then
  loom_path="${input_file}"
else
  loom_path="${scenic_input_path}/${input_file}"
fi

if [[ ! -f "${loom_path}" ]]; then
  log "ERROR: loom not found: ${loom_path}"
  exit 1
fi

export SCENIC_N_WORKERS="${SLURM_CPUS_PER_TASK:-24}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

log "cpus=${SCENIC_N_WORKERS} mem_mb=${SLURM_MEM_PER_NODE:-?} job=${SLURM_JOB_ID:-local}"
log "Loom: ${loom_path}"

cd "${WKDIR}" || exit 1
"${RSCRIPT}" "${RUNNER}" --loom "${loom_path}"
status=$?
log "run_scenic_one_split.R exit status: ${status}"
exit ${status}
