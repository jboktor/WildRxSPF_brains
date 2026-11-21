#!/bin/bash
#SBATCH --job-name="seurat_gene_meta_sn"
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/popalign/log_seurat_gene_meta_sn_%j.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/popalign/log_seurat_gene_meta_sn_%j.err"
#SBATCH --time=0-02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=400G
#SBATCH --partition=expansion
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END
#
# Write the snRNA Seurat RNA feature metadata (gene_id, gene_name, gene_type, ...)
# to data/interim/seurat/snRNASeq/seurat_gene_metadata_2026-10-05.csv.
# Submit: sbatch notebooks/shell_scripts/run_export_seurat_gene_metadata_sn.sh

set -euo pipefail
WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
RSCRIPT="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript"
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"; }
cd "${WKDIR}"
log "Exporting Seurat gene metadata"
"${RSCRIPT}" notebooks/R_scripts/export_seurat_gene_metadata_sn.R
log "Done"
