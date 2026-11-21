#!/usr/bin/env Rscript
# Export popAlign inputs (10x-style mtx per analysis unit) from the snRNA-seq
# Seurat object, at three Allen taxonomy levels: tissue × class, tissue × subclass
# and tissue × supertype.
#
# Modes:
#   --mode manifest   Light. Reads data/interim/dspin/cell_meta_keep_cell.csv,
#                     applies the unit inclusion rule, writes
#                     data/interim/popalign/snRNASeq/manifest_<run-tag>.csv. No Seurat load.
#   --run-tag         Results folder / manifest suffix (default 2026-10-04). Each run
#                     tag gets its own results_<tag>/ so earlier runs are never overwritten.
#   Inputs are shared across runs: export skips units whose input folder is complete.
#   --mode export     Heavy (~400G RAM). Loads Seurat once and writes, per unit:
#                     matrix.mtx (genes × cells, raw UMI), genes.tsv, barcodes.tsv,
#                     metadata.csv (cell_barcode, sample_id, ... , cell_type).
#   --mode both       manifest then export.
#
# Unit rule (both levels): unit has >= min_cells_per_animal cells in EVERY one of
# the 12 animals of that tissue (6 SPF + 6 WildR). Animals = popAlign samples.
# Large units are subsampled to at most max_cells_per_animal cells per animal
# (uniform within animal, so per-animal composition is preserved).
#
# Path key: {level}/{tissue}__{label}  (Allen labels only; names contain "/").
#
# Usage:
#   Rscript notebooks/R_scripts/export_popalign_inputs_sn.R --mode manifest
#   sbatch notebooks/shell_scripts/run_popalign_export_sn.sh      # --mode export

suppressPackageStartupMessages({
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
mode <- get_arg("--mode", "manifest")
only_level <- get_arg("--level", NULL) # optional: class | subclass | supertype
run_tag <- get_arg("--run-tag", "2026-10-04")
stopifnot(mode %in% c("manifest", "export", "both"))

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
data_path <- file.path(wkdir, "data/interim")
pa_dir <- file.path(data_path, "popalign/snRNASeq")
input_root <- file.path(pa_dir, "inputs")
result_root <- file.path(pa_dir, glue("results_{run_tag}"))
manifest_csv <- file.path(pa_dir, glue("manifest_{run_tag}.csv"))
dir.create(input_root, showWarnings = FALSE, recursive = TRUE)

cell_meta_csv <- file.path(data_path, "dspin/cell_meta_keep_cell.csv")
seurat_rds <- file.path(
  data_path, "seurat/snRNASeq",
  "seurat_obj_onlyparsefilters_celltyped_2025-12-27.rds"
)

MIN_CELLS_ANIMAL <- 25L
MAX_CELLS_ANIMAL <- 2500L
min_cells_gene <- 10L
n_animals_expected <- 12L
seed <- 42L

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", paste0(...), "\n")
  flush.console()
}

level_cols <- c(class = "class_label", subclass = "subclass_label", supertype = "supertype_label")
name_cols <- c(class = "class_name", subclass = "subclass_name", supertype = "supertype_name")
# next-finer Allen level stored as metadata `cell_type` (annotation / purity only)
finer_cols <- c(class = "subclass_label", subclass = "supertype_label", supertype = "cluster_label")

# ── manifest ─────────────────────────────────────────────────────────────────
build_manifest <- function() {
  meta <- read_csv(cell_meta_csv, show_col_types = FALSE) %>%
    filter(keep_cell, tissue %in% c("AMY", "HYP"))
  n_animals <- meta %>% distinct(tissue, sample) %>% count(tissue, name = "n_animals_tissue")

  man <- map_dfr(names(level_cols), function(lv) {
    lab <- level_cols[[lv]]; nm <- name_cols[[lv]]
    meta %>%
      filter(!is.na(.data[[lab]])) %>%
      count(tissue, label = .data[[lab]], name_ = .data[[nm]], sample, microbiome, name = "n") %>%
      group_by(tissue, label) %>%
      summarise(
        name = first(name_),
        n_cells = sum(n),
        n_SPF = sum(n[microbiome == "SPF"]),
        n_WildR = sum(n[microbiome == "WildR"]),
        n_animals = n_distinct(sample),
        min_cells_per_animal = min(n),
        n_cells_exported = sum(pmin(n, MAX_CELLS_ANIMAL)),
        .groups = "drop"
      ) %>%
      mutate(level = lv)
  }) %>%
    left_join(n_animals, by = "tissue") %>%
    mutate(
      # a unit missing from an animal has n_animals < 12 → fails
      pass = n_animals == n_animals_tissue & n_animals_tissue == n_animals_expected &
        min_cells_per_animal >= MIN_CELLS_ANIMAL,
      unit_id = glue("{level}__{tissue}__{label}"),
      input_dir = file.path(input_root, level, glue("{tissue}__{label}")),
      output_dir = file.path(result_root, level, glue("{tissue}__{label}"))
    ) %>%
    mutate(level = factor(level, levels = names(level_cols))) %>%
    arrange(level, tissue, label) %>% mutate(level = as.character(level))

  write_csv(man, file.path(pa_dir, glue("candidate_units_{run_tag}.csv")))
  run_man <- man %>% filter(pass) %>%
    mutate(array_index = row_number() - 1L) %>%
    relocate(array_index, unit_id, level, tissue, label, name)
  write_csv(run_man, manifest_csv)
  log_msg("Candidates: ", nrow(man), " | passing: ", nrow(run_man))
  print(run_man %>% count(level, tissue))
  run_man
}

# ── export ───────────────────────────────────────────────────────────────────
write_unit <- function(counts_all, md_all, row) {
  lab_col <- level_cols[[row$level]]
  cells <- md_all %>%
    filter(tissue == row$tissue, .data[[lab_col]] == row$label)
  set.seed(seed)
  cells <- cells %>%
    group_by(sample) %>%
    group_modify(~ if (nrow(.x) > MAX_CELLS_ANIMAL) slice_sample(.x, n = MAX_CELLS_ANIMAL) else .x) %>%
    ungroup()

  m <- counts_all[, cells$cell_id, drop = FALSE]
  keep_gene <- Matrix::rowSums(m > 0) >= min_cells_gene
  m <- m[keep_gene, , drop = FALSE]

  dir.create(row$input_dir, showWarnings = FALSE, recursive = TRUE)
  Matrix::writeMM(m, file.path(row$input_dir, "matrix.mtx"))
  write_tsv(tibble(a = rownames(m), b = rownames(m)), file.path(row$input_dir, "genes.tsv"), col_names = FALSE)
  write_tsv(tibble(a = colnames(m)), file.path(row$input_dir, "barcodes.tsv"), col_names = FALSE)

  # cell_type = the next-finer Allen level within the unit (used only for
  # annotation / purity; not for the GMMs).
  ct_col <- finer_cols[[row$level]]
  if (!ct_col %in% names(cells)) ct_col <- "supertype_label"
  md_out <- cells %>%
    transmute(
      cell_barcode = cell_id, sample_id = sample, microbiome, tissue,
      across(any_of(c("class_label", "subclass_label", "supertype_label", "cluster_label"))),
      cell_type = .data[[ct_col]]
    )
  write_csv(md_out, file.path(row$input_dir, "metadata.csv"))
  tibble(unit_id = row$unit_id, n_cells = ncol(m), n_genes = nrow(m))
}

input_complete <- function(d) all(file.exists(file.path(d, c("matrix.mtx", "genes.tsv", "barcodes.tsv", "metadata.csv"))))

export_units <- function(man) {
  if (!is.null(only_level)) man <- man %>% filter(level == only_level)
  todo <- man %>% filter(!map_lgl(input_dir, input_complete))
  log_msg("Units: ", nrow(man), " | inputs already complete: ", nrow(man) - nrow(todo), " | to export: ", nrow(todo))
  if (!nrow(todo)) return(invisible(NULL))
  man <- todo
  suppressPackageStartupMessages(library(Seurat))
  log_msg("Loading Seurat: ", seurat_rds)
  so <- readRDS(seurat_rds)
  md <- so@meta.data %>% rownames_to_column("cell_id") %>%
    filter(keep_cell, tissue %in% c("AMY", "HYP")) %>%
    mutate(across(any_of(c("sample", "microbiome", "tissue", "class_label",
                           "subclass_label", "supertype_label", "cluster_label")), as.character))
  counts <- GetAssayData(so, assay = "RNA", layer = "counts")[, md$cell_id]
  rm(so); gc()

  sym <- sub("-GRCm39$", "", rownames(counts))
  keep <- !duplicated(sym) & nzchar(sym) & !grepl("^ENSMUSG", sym) &
    !grepl("^mt-", sym) & !grepl("^Rp[sl][0-9]", sym)   # popAlign's own filter is human-only
  counts <- counts[keep, ]
  rownames(counts) <- sym[keep]
  log_msg("Counts: ", nrow(counts), " genes × ", ncol(counts), " cells")

  res <- map_dfr(seq_len(nrow(man)), function(i) {
    row <- man[i, ]
    log_msg("[", i, "/", nrow(man), "] ", row$unit_id)
    write_unit(counts, md, row)
  })
  write_csv(left_join(man, res, by = "unit_id"), file.path(pa_dir, glue("export_manifest_{run_tag}.csv")))
  log_msg("Export complete → ", input_root)
}

man <- if (mode %in% c("manifest", "both")) build_manifest() else
  read_csv(manifest_csv, show_col_types = FALSE)
if (mode %in% c("export", "both")) export_units(man)
