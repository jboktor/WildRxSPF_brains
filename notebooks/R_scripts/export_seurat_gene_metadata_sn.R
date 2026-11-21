#!/usr/bin/env Rscript
# Export the snRNA-seq Seurat object's RNA feature metadata (meta.features) to CSV.
#
# Row names carry a genome tag plus Seurat de-duplication suffixes
# ("Arhgap26-GRCm39", "Arhgap26-GRCm39-1"): the suffixed rows are DIFFERENT
# Ensembl genes sharing a symbol. gene_id (Ensembl) is the unambiguous key;
# gene_name is the plain symbol. `export_symbol` reproduces the feature names
# used by the popAlign / D-SPIN exports (sub("-GRCm39$", "", rowname)) so those
# outputs can be mapped back to Ensembl IDs.
#
# Usage: sbatch notebooks/shell_scripts/run_export_seurat_gene_metadata_sn.sh

suppressPackageStartupMessages({ library(Seurat); library(tidyverse) })
wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
seurat_rds <- file.path(wkdir, "data/interim/seurat/snRNASeq/seurat_obj_onlyparsefilters_celltyped_2025-12-27.rds")
out_csv <- file.path(wkdir, "data/interim/seurat/snRNASeq/seurat_gene_metadata_2026-10-05.csv")
stopifnot(!file.exists(out_csv))   # never overwrite a data artifact

message(Sys.time(), " loading ", seurat_rds)
so <- readRDS(seurat_rds)
mf <- so[["RNA"]][[]]
message("meta.features columns: ", paste(colnames(mf), collapse = ", "), " | rows: ", nrow(mf))
gm <- mf %>%
  rownames_to_column("feature") %>%
  mutate(export_symbol = sub("-GRCm39$", "", feature),
         symbol_shared = duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE))
write_csv(gm, out_csv)
message(Sys.time(), " wrote ", out_csv, " (", nrow(gm), " features; ",
        sum(gm$symbol_shared), " share a gene_name with another Ensembl ID)")
