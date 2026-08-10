#!/bin/bash
#SBATCH --job-name="dspin_david"
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/dspin/log_dspin_david_%A_%a.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/dspin/log_dspin_david_%A_%a.err"
#SBATCH --time=0-04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=expansion
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --array=0-139%5
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END
#
# DAVID annotation for each full_manifest DSPIN unit (species=mouse).
# Requires DAVID-registered email in DAVID_TOKEN or .secrets/david_token
#
#   export DAVID_TOKEN='you@caltech.edu'   # registered at davidbioinformatics.nih.gov
#   sbatch --export=ALL,DAVID_TOKEN notebooks/shell_scripts/run_dspin_david_array.sh

set -euo pipefail

WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
PY="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python"
MANIFEST="${WKDIR}/data/interim/dspin/full_manifest.csv"
IDX="${SLURM_ARRAY_TASK_ID:-0}"
SECRET="${WKDIR}/.secrets/david_token"

# Fix mygene/httpx SSL on cluster nodes
export SSL_CERT_FILE="$("${PY}" -c 'import certifi; print(certifi.where())')"
export REQUESTS_CA_BUNDLE="${SSL_CERT_FILE}"
export CURL_CA_BUNDLE="${SSL_CERT_FILE}"

mkdir -p "${WKDIR}/.cluster_runs/dspin"
cd "${WKDIR}"

if [ -z "${DAVID_TOKEN:-}" ] && [ -f "${SECRET}" ]; then
  # shellcheck disable=SC2155
  export DAVID_TOKEN="$(head -n1 "${SECRET}" | tr -d '\r')"
fi
if [ -z "${DAVID_TOKEN:-}" ]; then
  echo "ERROR: DAVID_TOKEN unset and ${SECRET} missing" >&2
  exit 2
fi

echo "[$(date)] DAVID array task ${IDX}"
"${PY}" notebooks/python_scripts/dspin_david_one.py \
  --manifest "${MANIFEST}" \
  --index "${IDX}"
echo "[$(date)] array task ${IDX} done"
