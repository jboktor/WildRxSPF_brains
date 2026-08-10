#!/bin/bash
#SBATCH --job-name="dspin_fit_arr"
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/dspin/log_dspin_fit_arr_%A_%a.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/dspin/log_dspin_fit_arr_%A_%a.err"
#SBATCH --time=0-06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --partition=expansion
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --array=0-139%25
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END
#
# Fit one full_manifest row per array task.
# Submit AFTER export+prep:
#   sbatch notebooks/shell_scripts/run_dspin_fit_array.sh
#
# Override array range if needed:
#   sbatch --array=0-139%20 notebooks/shell_scripts/run_dspin_fit_array.sh

set -euo pipefail

WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
PY="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python"
MANIFEST="${WKDIR}/data/interim/dspin/full_manifest.csv"
IDX="${SLURM_ARRAY_TASK_ID:-0}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${WKDIR}/.cluster_runs/dspin" "${WKDIR}/figures/DSPIN/snRNASeq"
cd "${WKDIR}"

echo "[$(date)] array task ${IDX} starting"
"${PY}" notebooks/python_scripts/dspin_fit_one.py \
  --manifest "${MANIFEST}" \
  --index "${IDX}" \
  --num-spin 20 \
  --num-repeat 10 \
  --skip-ok
echo "[$(date)] array task ${IDX} done"
