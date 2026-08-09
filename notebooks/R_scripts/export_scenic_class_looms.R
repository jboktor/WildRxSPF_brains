#!/usr/bin/env Rscript
# DEPRECATED wrapper — use export_scenic_looms.R --mode class|tissue|both
# (raw RNA counts; class + optional tissue splits).
#
# Export one pySCENIC loom per tissue × Allen class_label (FluxAnalysis pattern).
# Keeps both microbiome groups in each matrix so GRN is shared and WildR vs SPF
# AUCell contrasts are on the same regulons.
#
# Outputs:
#   data/input/scenic/snRNASeq/snRNASeq_{tissue}_{class_label}.loom
#   data/input/scenic/snRNASeq/snRNASeq_{tissue}_{class_label}_meta.csv
#   data/input/scenic/snRNASeq/manifest.csv

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(glue)
})

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
data_path <- file.path(wkdir, "data/interim")
out_dir <- file.path(wkdir, "data/input/scenic/snRNASeq")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

py_bin <- "/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python"
seurat_rds <- file.path(
  data_path, "seurat/snRNASeq",
  "seurat_obj_onlyparsefilters_celltyped_2025-12-27.rds"
)

min_cells_total <- 100L          # same floor as FluxAnalysis
min_cells_per_microbiome <- 30L  # need both groups for downstream tests
min_cells_gene <- 10L
min_genes_cell <- 200L

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", paste0(...), "\n")
  flush.console()
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

log_msg("Loading Seurat")
seurat_obj <- readRDS(seurat_rds)
seurat_obj <- subset(seurat_obj, subset = keep_cell)

tissues <- c("AMY", "HYP")
manifest <- list()

for (ts in tissues) {
  seurat_ts <- subset(seurat_obj, subset = tissue == ts)
  celltype_counts <- table(seurat_ts$class_label)
  celltype_labels <- names(celltype_counts)[celltype_counts > min_cells_total]

  for (celltype in celltype_labels) {
    stem <- glue("snRNASeq_{ts}_{celltype}")
    loom_path <- file.path(out_dir, glue("{stem}.loom"))
    meta_path <- file.path(out_dir, glue("{stem}_meta.csv"))

    seurat_ct <- subset(seurat_ts, subset = class_label == celltype)
    n_spf <- sum(seurat_ct$microbiome == "SPF")
    n_wildr <- sum(seurat_ct$microbiome == "WildR")
    if (ncol(seurat_ct) < min_cells_total) next
    if (n_spf < min_cells_per_microbiome || n_wildr < min_cells_per_microbiome) {
      log_msg("Skip ", stem, " — microbiome imbalance SPF=", n_spf, " WildR=", n_wildr)
      next
    }

    if (file.exists(loom_path) && file.info(loom_path)$size > 1000) {
      log_msg("Exists, skip: ", loom_path)
      manifest[[length(manifest) + 1]] <- tibble(
        tissue = ts, class_label = celltype,
        class_name = unique(as.character(seurat_ct$class_name))[1],
        n_cells = ncol(seurat_ct), n_SPF = n_spf, n_WildR = n_wildr,
        loom = basename(loom_path), meta = basename(meta_path),
        status = "exists"
      )
      next
    }

    log_msg("Exporting ", stem, " n=", ncol(seurat_ct),
            " SPF=", n_spf, " WildR=", n_wildr)

    counts <- GetAssayData(seurat_ct, layer = "counts")
    gene_keep <- Matrix::rowSums(counts > 0) >= min_cells_gene
    cell_keep <- Matrix::colSums(counts > 0) >= min_genes_cell
    counts <- counts[gene_keep, cell_keep, drop = FALSE]
    seurat_ct <- seurat_ct[, colnames(counts)]

    gene_symbols <- sub("-GRCm39$", "", rownames(counts))
    keep_sym <- !duplicated(gene_symbols) & nzchar(gene_symbols) & !grepl("^ENSMUSG", gene_symbols)
    counts <- counts[keep_sym, , drop = FALSE]
    rownames(counts) <- gene_symbols[keep_sym]

    meta <- seurat_ct@meta.data %>%
      tibble::rownames_to_column("cell_id") %>%
      select(
        cell_id, sample, microbiome, tissue, group,
        class_label, class_name, subclass_name, supertype_name, broad_cell_group
      )
    write_csv(meta, meta_path)

    tmp_prefix <- file.path(tempdir(), stem)
    write_loom_from_counts(counts, meta, loom_path, tmp_prefix)

    manifest[[length(manifest) + 1]] <- tibble(
      tissue = ts, class_label = celltype,
      class_name = unique(as.character(meta$class_name))[1],
      n_cells = ncol(counts), n_genes = nrow(counts),
      n_SPF = sum(meta$microbiome == "SPF"),
      n_WildR = sum(meta$microbiome == "WildR"),
      loom = basename(loom_path), meta = basename(meta_path),
      status = "wrote"
    )
  }
}

manifest_df <- bind_rows(manifest)
write_csv(manifest_df, file.path(out_dir, "manifest.csv"))
log_msg("Wrote manifest with ", nrow(manifest_df), " splits to ", out_dir)
print(manifest_df)
