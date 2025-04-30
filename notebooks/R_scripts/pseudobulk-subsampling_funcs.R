library(tidyverse)
library(magrittr)
library(glue)
library(Seurat)
library(grid)
library(biomaRt)
library(batchtools)
library(strex)
# stats 
library(lme4)
library(broom)
# parallelization
library(BiocParallel)
library(future)
library(furrr)
library(foreach)
library(doParallel)

# plotting
library(RColorBrewer)
library(ggsci)
library(Matrix)
library(plotly)
library(DT)
library(gtsummary)
library(aplot)
library(patchwork)

# Load libraries
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(data.table)
library(DEGreport)

# setting paths
homedir <- "/central/groups/mthomson/jboktor"
wkdir <- glue("{homedir}/spatial_genomics/jess_2024-01-23")
source(glue("{wkdir}/notebooks/R_scripts/_misc_functions.R"))

subsamp_dir <- glue("{wkdir}/data/interim/DESeq2/subsampling/",
"gated-cell-types_12-samples_data"
)

count_df <- readRDS(
  glue("{subsamp_dir}/count_df_2025-01-05.rds")
)
sample_n_list <- readRDS(
  glue("{subsamp_dir}/sample_n_list_2025-01-05.rds")
)
core_deseq_meta_df <- readRDS(
  glue("{subsamp_dir}/core_deseq_meta_df_{Sys.Date()}.rds")
)


# Functions for analysis
subsample_celltype <- function(df, cell_id, subsamp_n, seed) {
  require(magrittr)
  require(dplyr)
  set.seed(seed)
  res <- list()

  # Counting cells per slice for the specified cell type
  slice_counts <- df %>%
    dplyr::filter(cluster_id == cell_id) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarize(n = n())

  # Filtering for slices with at least the specified number of cells
  valid_slices <- slice_counts %>%
    dplyr::filter(n >= subsamp_n) %>%
    dplyr::pull(sample_id)
  
  subsamp_df <- df %>%
    dplyr::filter(cluster_id == cell_id, sample_id %in% valid_slices) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::sample_n(size = subsamp_n, replace = FALSE) %>%
    dplyr::ungroup(sample_id)

  res[["subsamp_counts"]] <- subsamp_df %>%
    dplyr::select(-c(sample_id, cluster_id, group, run, roi, bregma))
  res[["subsamp_meta"]] <- subsamp_df %>%
    dplyr::select(c(sample_id))
  return(res)
}
