#!/bin/bash
#SBATCH --job-name="dspin_migrate"
#SBATCH --output="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/dspin/log_dspin_migrate_%j.log"
#SBATCH --error="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/.cluster_runs/dspin/log_dspin_migrate_%j.err"
#SBATCH --time=0-06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=expansion
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL,END
#
# Rebuild label-keyed candidates from cell_meta, migrate dirs slug→label,
# rank full_manifest, finish HVG prep. Then submit fit array with afterok.

set -euo pipefail
WKDIR="/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
RSCRIPT="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript"
PY="/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python"
cd "${WKDIR}"
mkdir -p .cluster_runs/dspin

echo "[$(date)] rebuild candidates from cell_meta (label-keyed)"
"${RSCRIPT}" - <<'RS'
suppressPackageStartupMessages(library(tidyverse))
wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
out_dir <- file.path(wkdir, "data/interim/dspin")
meta <- read_csv(file.path(out_dir, "cell_meta_keep_cell.csv"), show_col_types = FALSE)
min_cells_per_microbiome <- 30L
min_cells_per_sample_id_floor <- 25L
min_samples_per_microbiome <- 2L
tissues <- c("AMY", "HYP")
cand <- meta %>%
  filter(keep_cell, tissue %in% tissues,
         !is.na(supertype_label), nzchar(as.character(supertype_label)),
         !is.na(supertype_name), nzchar(as.character(supertype_name))) %>%
  mutate(sample_id = as.character(sample),
         microbiome = as.character(microbiome),
         tissue = as.character(tissue),
         supertype_label = as.character(supertype_label),
         supertype_name = as.character(supertype_name)) %>%
  group_by(tissue, supertype_label, sample_id, microbiome) %>%
  summarise(n_cells = n(),
            supertype_name = dplyr::first(supertype_name),
            .groups = "drop") %>%
  group_by(tissue, supertype_label) %>%
  summarise(
    supertype_name = dplyr::first(supertype_name),
    n_SPF = sum(n_cells[microbiome == "SPF"]),
    n_WildR = sum(n_cells[microbiome == "WildR"]),
    min_cells_per_sample_id = min(n_cells),
    n_sample_SPF = n_distinct(sample_id[microbiome == "SPF"]),
    n_sample_WildR = n_distinct(sample_id[microbiome == "WildR"]),
    n_cells = sum(n_cells),
    .groups = "drop"
  ) %>%
  mutate(
    batch_control_ok = (n_SPF >= min_cells_per_sample_id_floor) &
      (n_WildR >= min_cells_per_sample_id_floor),
    hard_pass =
      (n_SPF >= min_cells_per_microbiome) &
      (n_WildR >= min_cells_per_microbiome) &
      (min_cells_per_sample_id >= min_cells_per_sample_id_floor) &
      (n_sample_SPF >= min_samples_per_microbiome) &
      (n_sample_WildR >= min_samples_per_microbiome) &
      batch_control_ok,
    unit_id = paste(tissue, supertype_label, sep = "__")
  )
write_csv(cand, file.path(out_dir, "candidate_cellcounts.csv"))
cat("units=", nrow(cand), " hard_pass=", sum(cand$hard_pass), "\n")
RS

echo "[$(date)] migrate dirs name-slug → label"
"${PY}" notebooks/python_scripts/dspin_migrate_to_labels.py

echo "[$(date)] rank → full_manifest"
"${RSCRIPT}" notebooks/R_scripts/rank_dspin_pilots.R

echo "[$(date)] HVG prep (skip exists)"
"${PY}" notebooks/python_scripts/dspin_prep_filtered.py \
  --manifest "${WKDIR}/data/interim/dspin/full_manifest.csv"

n_raw=$(ls -1 data/interim/dspin/*/raw_counts.h5ad 2>/dev/null | wc -l)
n_filt=$(ls -1 data/interim/dspin/*/filtered.h5ad 2>/dev/null | wc -l)
echo "[$(date)] raw=${n_raw} filtered=${n_filt} Done migrate+prep"
