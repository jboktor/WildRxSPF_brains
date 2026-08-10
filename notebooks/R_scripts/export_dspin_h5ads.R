#!/usr/bin/env Rscript
# Export DSPIN AnnData (raw UMI counts) + candidate cell-count tables.
#
# Modes:
#   --mode candidates     Load Seurat once; write candidate_cellcounts.csv +
#                         cell_meta_keep_cell.csv (AMY/HYP, keep_cell). No h5ad.
#   --mode export-pilot   Read manifest; export raw-count h5ad per
#                         (tissue, supertype_label) analysis unit.
#   --mode both           candidates then export-pilot.
#
# DSPIN metadata (written into adata.obs):
#   sample_id  <- sample   (one animal = one DSPIN condition / h)
#   if_control <- microbiome == "SPF"
#   batch      <- tissue
#
# Path key: tissue__{supertype_label}  (filesystem-safe; never use
# supertype_name — names contain "/" which create nested folders).
#
# Cells: keep_cell + tissue + type only. No extra per-cell QC subsetting
# (no min_genes_cell downfilter). Gene sparsity filter only.
#
# Usage:
#   Rscript notebooks/R_scripts/export_dspin_h5ads.R --mode candidates
#   Rscript notebooks/R_scripts/export_dspin_h5ads.R --mode export-pilot --tier full

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(glue)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1]]
}
mode <- get_arg("--mode", "candidates")
tier_filter <- get_arg("--tier", "interactive") # interactive | all | optional_sbatch | full | pool
manifest_arg <- get_arg("--manifest", NULL)
stopifnot(mode %in% c("candidates", "export-pilot", "both"))
stopifnot(tier_filter %in% c("interactive", "all", "optional_sbatch", "full", "pool"))

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
data_path <- file.path(wkdir, "data/interim")
out_dir <- file.path(data_path, "dspin")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

py_bin <- "/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python"
seurat_rds <- file.path(
  data_path, "seurat/snRNASeq",
  "seurat_obj_onlyparsefilters_celltyped_2025-12-27.rds"
)

min_cells_per_microbiome <- 30L
min_cells_per_sample_id_floor <- 25L
min_samples_per_microbiome <- 2L
min_cells_gene <- 10L
tissues <- c("AMY", "HYP")

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", paste0(...), "\n")
  flush.console()
}

# Directory stem = Allen label (safe). Never path-encode name.
unit_dir <- function(tissue, supertype_label) {
  file.path(out_dir, paste0(tissue, "__", as.character(supertype_label)))
}

assert_meta_cols <- function(md) {
  needed <- c(
    "sample", "microbiome", "tissue", "keep_cell",
    "supertype_label", "supertype_name", "class_label", "class_name"
  )
  missing <- setdiff(needed, colnames(md))
  if (length(missing)) {
    stop("Seurat meta missing columns: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

compute_candidate_table <- function(meta) {
  meta %>%
    filter(
      keep_cell,
      tissue %in% tissues,
      !is.na(supertype_label), nzchar(as.character(supertype_label)),
      !is.na(supertype_name), nzchar(as.character(supertype_name))
    ) %>%
    mutate(
      sample_id = as.character(sample),
      microbiome = as.character(microbiome),
      tissue = as.character(tissue),
      supertype_label = as.character(supertype_label),
      supertype_name = as.character(supertype_name)
    ) %>%
    group_by(tissue, supertype_label, sample_id, microbiome) %>%
    summarise(
      n_cells = n(),
      # name is display-only / MAST join; 1:1 with label
      supertype_name = dplyr::first(supertype_name),
      .groups = "drop"
    ) %>%
    {
      per_sample <- .
      per_sample %>%
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
    }
}

write_h5ad_from_counts <- function(counts, meta, h5ad_path, tmp_prefix) {
  counts <- as(counts, "dgCMatrix")
  i_bin <- paste0(tmp_prefix, ".i.bin")
  p_bin <- paste0(tmp_prefix, ".p.bin")
  x_bin <- paste0(tmp_prefix, ".x.bin")
  dims_txt <- paste0(tmp_prefix, ".dims.txt")
  genes_tsv <- paste0(tmp_prefix, ".genes.tsv")
  bars_tsv <- paste0(tmp_prefix, ".barcodes.tsv")
  meta_csv <- paste0(tmp_prefix, ".meta.csv")

  writeBin(as.integer(counts@i), i_bin, size = 4)
  writeBin(as.integer(counts@p), p_bin, size = 4)
  writeBin(as.double(counts@x), x_bin, size = 8)
  writeLines(as.character(c(nrow(counts), ncol(counts), length(counts@x))), dims_txt)
  write_tsv(tibble(gene = rownames(counts)), genes_tsv, col_names = FALSE)
  write_tsv(tibble(cell = colnames(counts)), bars_tsv, col_names = FALSE)
  write_csv(meta, meta_csv)

  py <- glue(
    '
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

prefix = Path("{tmp_prefix}")
dims = [int(x) for x in open(str(prefix) + ".dims.txt")]
nrow, ncol, nnz = dims
i = np.fromfile(str(prefix) + ".i.bin", dtype=np.int32)
p = np.fromfile(str(prefix) + ".p.bin", dtype=np.int32)
x = np.fromfile(str(prefix) + ".x.bin", dtype=np.float64)
counts = sparse.csc_matrix((x, i, p), shape=(nrow, ncol)).astype(np.float32).tocsr()
genes = pd.read_csv(str(prefix) + ".genes.tsv", header=None)[0].astype(str).tolist()
cells = pd.read_csv(str(prefix) + ".barcodes.tsv", header=None)[0].astype(str).tolist()
meta = pd.read_csv(str(prefix) + ".meta.csv").set_index("cell_id").loc[cells]
meta = meta.copy()
meta["sample_id"] = meta["sample"].astype(str)
meta["batch"] = meta["tissue"].astype(str)
meta["if_control"] = (meta["microbiome"].astype(str) == "SPF")
adata = ad.AnnData(X=counts.T, obs=meta, var=pd.DataFrame(index=genes))
adata.var_names_make_unique()
out = Path("{h5ad_path}")
out.parent.mkdir(parents=True, exist_ok=True)
if out.exists():
    out.unlink()
adata.write_h5ad(out, compression="gzip")
print("Wrote", out, "cells", adata.n_obs, "genes", adata.n_vars)
'
  )
  tmp_py <- tempfile(fileext = ".py")
  writeLines(py, tmp_py)
  status <- system2(py_bin, tmp_py)
  unlink(c(i_bin, p_bin, x_bin, dims_txt, genes_tsv, bars_tsv, meta_csv, tmp_py))
  if (status != 0) stop("anndata write failed for ", h5ad_path)
}

prepare_counts <- function(seurat_sub) {
  DefaultAssay(seurat_sub) <- "RNA"
  counts <- GetAssayData(seurat_sub, layer = "counts")

  # Gene sparsity only — do NOT drop cells beyond keep_cell + type + tissue.
  gene_keep <- Matrix::rowSums(counts > 0) >= min_cells_gene
  counts <- counts[gene_keep, , drop = FALSE]

  gene_symbols <- sub("-GRCm39$", "", rownames(counts))
  keep_sym <- !duplicated(gene_symbols) & nzchar(gene_symbols) & !grepl("^ENSMUSG", gene_symbols)
  counts <- counts[keep_sym, , drop = FALSE]
  rownames(counts) <- gene_symbols[keep_sym]

  # Align meta to remaining columns (all cells kept)
  seurat_sub <- seurat_sub[, colnames(counts)]

  meta <- seurat_sub@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    transmute(
      cell_id,
      sample = as.character(sample),
      microbiome = as.character(microbiome),
      tissue = as.character(tissue),
      class_label = as.character(class_label),
      class_name = as.character(class_name),
      subclass_label = as.character(subclass_label),
      subclass_name = as.character(subclass_name),
      supertype_label = as.character(supertype_label),
      supertype_name = as.character(supertype_name),
      sample_id = as.character(sample),
      batch = as.character(tissue),
      if_control = microbiome == "SPF"
    )

  list(counts = counts, meta = meta)
}

export_unit <- function(seurat_obj, tissue_i, label_i) {
  stem_dir <- unit_dir(tissue_i, label_i)
  h5ad_path <- file.path(stem_dir, "raw_counts.h5ad")
  meta_path <- file.path(stem_dir, "cells_meta.csv")

  cells <- colnames(seurat_obj)[
    seurat_obj$keep_cell &
      seurat_obj$tissue == tissue_i &
      seurat_obj$supertype_label == label_i
  ]
  if (length(cells) < 1) {
    log_msg("Skip empty unit ", tissue_i, " / ", label_i)
    return(NULL)
  }
  seurat_sub <- seurat_obj[, cells]
  n_spf <- sum(seurat_sub$microbiome == "SPF")
  n_wildr <- sum(seurat_sub$microbiome == "WildR")
  name_i <- as.character(unique(seurat_sub$supertype_name))[1]
  log_msg(
    "Exporting ", tissue_i, "__", label_i,
    " (", name_i, ") n=", ncol(seurat_sub),
    " SPF=", n_spf, " WildR=", n_wildr
  )

  if (file.exists(h5ad_path) && file.info(h5ad_path)$size > 1000) {
    log_msg("Exists, skip: ", h5ad_path)
    return(tibble(
      tissue = tissue_i,
      supertype_label = as.character(label_i),
      supertype_name = name_i,
      n_cells = ncol(seurat_sub),
      n_SPF = n_spf,
      n_WildR = n_wildr,
      h5ad = h5ad_path,
      status = "exists"
    ))
  }

  dir.create(stem_dir, showWarnings = FALSE, recursive = TRUE)
  prep <- prepare_counts(seurat_sub)
  write_csv(prep$meta, meta_path)
  write_h5ad_from_counts(
    prep$counts, prep$meta, h5ad_path,
    file.path(tempdir(), paste0(tissue_i, "_", gsub("[^A-Za-z0-9]+", "_", label_i)))
  )

  tibble(
    tissue = tissue_i,
    supertype_label = as.character(label_i),
    supertype_name = name_i,
    n_cells = ncol(prep$counts),
    n_genes = nrow(prep$counts),
    n_SPF = sum(prep$meta$microbiome == "SPF"),
    n_WildR = sum(prep$meta$microbiome == "WildR"),
    h5ad = h5ad_path,
    status = "wrote"
  )
}

load_seurat <- function() {
  log_msg("Loading Seurat (snRNASeq): ", seurat_rds)
  seurat_obj <- readRDS(seurat_rds)
  assert_meta_cols(seurat_obj@meta.data)
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj
}

run_candidates <- function(seurat_obj = NULL) {
  if (is.null(seurat_obj)) seurat_obj <- load_seurat()
  meta <- seurat_obj@meta.data %>% tibble::rownames_to_column("cell_id")
  meta_path <- file.path(out_dir, "cell_meta_keep_cell.csv")
  meta_keep <- meta %>%
    filter(keep_cell, tissue %in% tissues) %>%
    select(
      cell_id, sample, microbiome, tissue, keep_cell,
      any_of(c(
        "class_label", "class_name", "subclass_label", "subclass_name",
        "supertype_label", "supertype_name"
      ))
    )
  write_csv(meta_keep, meta_path)
  log_msg("Wrote ", meta_path, " rows=", nrow(meta_keep))

  cand <- compute_candidate_table(meta)
  cand_path <- file.path(out_dir, "candidate_cellcounts.csv")
  write_csv(cand, cand_path)
  log_msg(
    "Wrote ", cand_path, " units=", nrow(cand),
    " hard_pass=", sum(cand$hard_pass)
  )
  invisible(cand)
}

run_export_pilot <- function(seurat_obj = NULL) {
  man_path <- if (!is.null(manifest_arg)) {
    manifest_arg
  } else if (identical(tier_filter, "full") && file.exists(file.path(out_dir, "full_manifest.csv"))) {
    file.path(out_dir, "full_manifest.csv")
  } else {
    file.path(out_dir, "pilot_manifest.csv")
  }
  if (!file.exists(man_path)) {
    stop("Missing ", man_path, " — run rank_dspin_pilots.R first.")
  }
  man <- read_csv(man_path, show_col_types = FALSE)
  if (!"run_tier" %in% names(man)) {
    stop(man_path, " needs run_tier column")
  }
  if (!"supertype_label" %in% names(man)) {
    stop(man_path, " needs supertype_label (filesystem-safe path key)")
  }
  if (tier_filter != "all") {
    man <- man %>% filter(run_tier == tier_filter)
  }
  if (!nrow(man)) stop("No rows in ", man_path, " for tier=", tier_filter)
  log_msg("Exporting ", nrow(man), " units from ", man_path, " (tier=", tier_filter, ")")

  if (is.null(seurat_obj)) seurat_obj <- load_seurat()
  rows <- list()
  for (i in seq_len(nrow(man))) {
    row <- export_unit(seurat_obj, man$tissue[[i]], man$supertype_label[[i]])
    if (!is.null(row)) rows[[length(rows) + 1]] <- row
  }
  exp_man <- bind_rows(rows)
  out_name <- if (identical(tier_filter, "full") || grepl("full_manifest", basename(man_path))) {
    "export_full_manifest.csv"
  } else {
    "export_pilot_manifest.csv"
  }
  write_csv(exp_man, file.path(out_dir, out_name))
  log_msg("Export complete: ", nrow(exp_man), " rows → ", out_name)
  invisible(exp_man)
}

log_msg("MODE=", mode, " tier=", tier_filter)
seurat_obj <- NULL
if (mode %in% c("candidates", "both", "export-pilot")) {
  if (mode %in% c("candidates", "both")) {
    seurat_obj <- load_seurat()
    run_candidates(seurat_obj)
  }
  if (mode %in% c("export-pilot", "both")) {
    if (is.null(seurat_obj) && file.exists(file.path(out_dir, "pilot_manifest.csv"))) {
      seurat_obj <- load_seurat()
    }
    if (mode == "both" && !file.exists(file.path(out_dir, "pilot_manifest.csv"))) {
      log_msg("No pilot_manifest.csv yet — skipping export-pilot (rank first).")
    } else {
      run_export_pilot(seurat_obj)
    }
  }
}
log_msg("Done → ", out_dir)
