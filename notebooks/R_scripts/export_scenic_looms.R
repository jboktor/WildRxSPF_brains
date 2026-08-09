#!/usr/bin/env Rscript
# Export pySCENIC looms from Seurat RNA counts (raw UMIs).
#
# Modes:
#   --mode class   tissue × Allen class_label (Flux-style; default)
#   --mode tissue  whole tissue (AMY and HYP separately)
#   --mode both    class then tissue
#
# Input verification: GetAssayData(layer="counts") must be non-negative
# integer-like UMIs (pySCENIC / GRNBoost2 recommended input). We do NOT
# export the Seurat "data" layer.
#
# Outputs under data/input/scenic/snRNASeq/:
#   snRNASeq_{tissue}_{class_label}.loom (+ _meta.csv)   [class]
#   snRNASeq_{tissue}.loom (+ _meta.csv)                 [tissue]
#   manifest_class.csv / manifest_tissue.csv

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
mode <- get_arg("--mode", "both")
stopifnot(mode %in% c("class", "tissue", "both"))

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
data_path <- file.path(wkdir, "data/interim")
out_dir <- file.path(wkdir, "data/input/scenic/snRNASeq")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

py_bin <- "/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python"
seurat_rds <- file.path(
  data_path, "seurat/snRNASeq",
  "seurat_obj_onlyparsefilters_celltyped_2025-12-27.rds"
)

min_cells_total <- 100L
min_cells_per_microbiome <- 30L
min_cells_gene <- 10L
min_genes_cell <- 200L

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", paste0(...), "\n")
  flush.console()
}

assert_raw_counts <- function(counts, label) {
  x <- counts@x
  if (!length(x)) stop(label, ": empty counts matrix")
  if (any(x < 0, na.rm = TRUE)) stop(label, ": negative values in counts")
  frac_int <- mean(abs(x - round(x)) < 1e-6)
  if (frac_int < 0.99) {
    stop(
      label, ": counts layer is NOT integer-like (frac_int=",
      signif(frac_int, 3), "). Refusing to export — pySCENIC expects raw UMIs."
    )
  }
  lib <- Matrix::colSums(counts)
  log_msg(
    label, " raw-count check OK: genes=", nrow(counts),
    " cells=", ncol(counts),
    " frac_integer_nnz=", signif(frac_int, 4),
    " lib_median=", signif(median(lib), 4),
    " lib_max=", signif(max(lib), 4)
  )
  invisible(TRUE)
}

write_loom_from_counts <- function(counts, meta, loom_path, tmp_prefix) {
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
import loompy
from scipy import sparse

prefix = Path("{tmp_prefix}")
dims = [int(x) for x in open(str(prefix) + ".dims.txt")]
nrow, ncol, nnz = dims
i = np.fromfile(str(prefix) + ".i.bin", dtype=np.int32)
p = np.fromfile(str(prefix) + ".p.bin", dtype=np.int32)
x = np.fromfile(str(prefix) + ".x.bin", dtype=np.float64)
counts = sparse.csc_matrix((x, i, p), shape=(nrow, ncol)).astype(np.float32)
genes = pd.read_csv(str(prefix) + ".genes.tsv", header=None)[0].astype(str).tolist()
cells = pd.read_csv(str(prefix) + ".barcodes.tsv", header=None)[0].astype(str).tolist()
meta = pd.read_csv(str(prefix) + ".meta.csv").set_index("cell_id").loc[cells]
row_attrs = {{"Gene": np.array(genes, dtype=object)}}
col_attrs = {{"CellID": np.array(cells, dtype=object)}}
for col in ["sample", "microbiome", "tissue", "class_label", "class_name",
            "subclass_name", "supertype_name"]:
    if col in meta.columns:
        col_attrs[col] = meta[col].astype(str).to_numpy()
out = Path("{loom_path}")
if out.exists():
    out.unlink()
loompy.create(str(out), counts, row_attrs, col_attrs)
print("Wrote", out, "shape", counts.shape, "nnz", counts.nnz)
'
  )
  tmp_py <- tempfile(fileext = ".py")
  writeLines(py, tmp_py)
  status <- system2(py_bin, tmp_py)
  unlink(c(i_bin, p_bin, x_bin, dims_txt, genes_tsv, bars_tsv, meta_csv, tmp_py))
  if (status != 0) stop("loompy failed for ", loom_path)
}

prepare_counts <- function(seurat_sub, label) {
  DefaultAssay(seurat_sub) <- "RNA"
  counts <- GetAssayData(seurat_sub, layer = "counts")
  assert_raw_counts(counts, label)

  gene_keep <- Matrix::rowSums(counts > 0) >= min_cells_gene
  cell_keep <- Matrix::colSums(counts > 0) >= min_genes_cell
  counts <- counts[gene_keep, cell_keep, drop = FALSE]
  seurat_sub <- seurat_sub[, colnames(counts)]

  # Parse gene names are Symbol-GRCm39; cisTarget TF lists use bare symbols
  gene_symbols <- sub("-GRCm39$", "", rownames(counts))
  keep_sym <- !duplicated(gene_symbols) & nzchar(gene_symbols) & !grepl("^ENSMUSG", gene_symbols)
  counts <- counts[keep_sym, , drop = FALSE]
  rownames(counts) <- gene_symbols[keep_sym]

  meta <- seurat_sub@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    select(
      any_of(c(
        "cell_id", "sample", "microbiome", "tissue", "group",
        "class_label", "class_name", "subclass_name", "supertype_name",
        "broad_cell_group"
      ))
    )

  list(counts = counts, meta = meta, seurat = seurat_sub)
}

export_one <- function(seurat_sub, stem, split_level, extra_manifest = list()) {
  loom_path <- file.path(out_dir, glue("{stem}.loom"))
  meta_path <- file.path(out_dir, glue("{stem}_meta.csv"))

  n_spf <- sum(seurat_sub$microbiome == "SPF")
  n_wildr <- sum(seurat_sub$microbiome == "WildR")
  if (ncol(seurat_sub) < min_cells_total) {
    log_msg("Skip ", stem, " — n=", ncol(seurat_sub), " < ", min_cells_total)
    return(NULL)
  }
  if (n_spf < min_cells_per_microbiome || n_wildr < min_cells_per_microbiome) {
    log_msg("Skip ", stem, " — microbiome imbalance SPF=", n_spf, " WildR=", n_wildr)
    return(NULL)
  }

  if (file.exists(loom_path) && file.info(loom_path)$size > 1000) {
    log_msg("Exists, skip: ", loom_path)
    return(tibble(
      split_level = split_level,
      tissue = unique(as.character(seurat_sub$tissue))[1],
      class_label = extra_manifest$class_label %||% NA_character_,
      class_name = extra_manifest$class_name %||% NA_character_,
      n_cells = ncol(seurat_sub), n_SPF = n_spf, n_WildR = n_wildr,
      loom = basename(loom_path), meta = basename(meta_path),
      status = "exists"
    ))
  }

  log_msg("Exporting ", stem, " n=", ncol(seurat_sub),
          " SPF=", n_spf, " WildR=", n_wildr)
  prep <- prepare_counts(seurat_sub, stem)
  write_csv(prep$meta, meta_path)
  write_loom_from_counts(prep$counts, prep$meta, loom_path, file.path(tempdir(), stem))

  tibble(
    split_level = split_level,
    tissue = unique(as.character(prep$meta$tissue))[1],
    class_label = extra_manifest$class_label %||% NA_character_,
    class_name = if (!is.null(extra_manifest$class_name)) {
      extra_manifest$class_name
    } else if ("class_name" %in% names(prep$meta)) {
      unique(as.character(prep$meta$class_name))[1]
    } else {
      NA_character_
    },
    n_cells = ncol(prep$counts),
    n_genes = nrow(prep$counts),
    n_SPF = sum(prep$meta$microbiome == "SPF"),
    n_WildR = sum(prep$meta$microbiome == "WildR"),
    loom = basename(loom_path),
    meta = basename(meta_path),
    status = "wrote"
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

log_msg("MODE=", mode)
log_msg("Loading Seurat (raw counts from RNA counts layer)")
seurat_obj <- readRDS(seurat_rds)
seurat_obj <- subset(seurat_obj, subset = keep_cell)
DefaultAssay(seurat_obj) <- "RNA"

# Global raw-count sanity check on a small sample before long export
sample_cells <- sample(colnames(seurat_obj), min(200L, ncol(seurat_obj)))
assert_raw_counts(
  GetAssayData(seurat_obj, layer = "counts")[, sample_cells],
  "global_sample"
)

tissues <- c("AMY", "HYP")
manifest_class <- list()
manifest_tissue <- list()

if (mode %in% c("class", "both")) {
  for (ts in tissues) {
    seurat_ts <- subset(seurat_obj, subset = tissue == ts)
    celltype_counts <- table(seurat_ts$class_label)
    celltype_labels <- names(celltype_counts)[celltype_counts > min_cells_total]

    for (celltype in celltype_labels) {
      seurat_ct <- subset(seurat_ts, subset = class_label == celltype)
      stem <- glue("snRNASeq_{ts}_{celltype}")
      row <- export_one(
        seurat_ct, stem, "class",
        extra_manifest = list(
          class_label = celltype,
          class_name = unique(as.character(seurat_ct$class_name))[1]
        )
      )
      if (!is.null(row)) manifest_class[[length(manifest_class) + 1]] <- row
    }
  }
  df <- bind_rows(manifest_class)
  write_csv(df, file.path(out_dir, "manifest_class.csv"))
  # Keep legacy name for Flux-style notebook chunk
  write_csv(df, file.path(out_dir, "manifest.csv"))
  log_msg("Class manifest rows: ", nrow(df))
}

if (mode %in% c("tissue", "both")) {
  for (ts in tissues) {
    seurat_ts <- subset(seurat_obj, subset = tissue == ts)
    stem <- glue("snRNASeq_{ts}")
    row <- export_one(seurat_ts, stem, "tissue")
    if (!is.null(row)) manifest_tissue[[length(manifest_tissue) + 1]] <- row
  }
  df <- bind_rows(manifest_tissue)
  write_csv(df, file.path(out_dir, "manifest_tissue.csv"))
  log_msg("Tissue manifest rows: ", nrow(df))
}

log_msg("Export complete → ", out_dir)
