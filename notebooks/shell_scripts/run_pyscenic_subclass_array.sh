#!/bin/bash
#SBATCH --job-name="pyscenic_subclass"
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/scenic/log_pyscenic_subclass_%A_%a.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/scenic/log_pyscenic_subclass_%A_%a.err"
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=expansion
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END

# Array worker: task i runs pySCENIC grn → ctx → aucell + the animal-level
# WildR vs SPF test (run_scenic_one_split.R) on line i+1 of a loom list.
# Outputs: data/interim/scenic/by_<level>/<stem>/, figures/SCENIC/snRNASeq_by_<level>/<stem>/
# (level inferred from the loom stem: _SUBC_ = subclass, _SUPT_ = supertype).
# Same commands for supertype with L=data/interim/scenic/looms_supertype and looms_supertype_*.txt.
#
# Two resource tiers (lists written by export_scenic_subclass_looms.py).
# Sizing from class runs: 49k nuclei used 85 GB / 2.4 h at 24 CPUs; 11k used 34 GB / 17 min.
#   L=data/interim/scenic/looms_subclass
#   N=$(grep -c . $L/looms_subclass_small.txt)
#   sbatch --array=0-$((N-1))%20 notebooks/shell_scripts/run_pyscenic_subclass_array.sh -l $L/looms_subclass_small.txt
#   N=$(grep -c . $L/looms_subclass_large.txt)
#   sbatch --array=0-$((N-1))%8 --cpus-per-task=24 --mem=160G --time=1-00:00:00 \
#     notebooks/shell_scripts/run_pyscenic_subclass_array.sh -l $L/looms_subclass_large.txt

set -euo pipefail

while getopts "l:" opt; do
  case $opt in
    l) list_file="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&2; exit 1;;
  esac
done

WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
RSCRIPT="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript"
RUNNER="${WKDIR}/notebooks/R_scripts/run_scenic_one_split.R"

log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"; }

: "${list_file:?pass -l <loom list file>}"
: "${SLURM_ARRAY_TASK_ID:?must run as an array task}"
loom=$(grep . "${list_file}" | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
[[ -n "${loom}" && -f "${loom}" ]] || { log "ERROR: no loom for task ${SLURM_ARRAY_TASK_ID} in ${list_file}"; exit 1; }

export SCENIC_N_WORKERS="${SLURM_CPUS_PER_TASK:-16}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

log "task=${SLURM_ARRAY_TASK_ID} job=${SLURM_ARRAY_JOB_ID:-?} cpus=${SCENIC_N_WORKERS} mem_mb=${SLURM_MEM_PER_NODE:-?}"
log "Loom: ${loom}"
cd "${WKDIR}"
# Rscript parses its file incrementally, so editing the runner while tasks are
# mid-run corrupts them (happened 2026-10-06: 27 subclass tasks died with a
# parse error). Run a private snapshot taken at task start instead.
snap="$(mktemp -d)/run_scenic_one_split.R"
cp "${RUNNER}" "${snap}"
"${RSCRIPT}" "${snap}" --loom "${loom}"
log "DONE $(basename "${loom}")"
