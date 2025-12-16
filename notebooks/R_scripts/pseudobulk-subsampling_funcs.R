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



run_deseq2_sensitivity_analysis <- function(wkdir, clust, nthreads) {
  
  require(glue)
  source(glue("{wkdir}/notebooks/R_scripts/pseudobulk-subsampling_funcs.R"))
  doParallel::registerDoParallel(cores = nthreads)

  #' Following variables are loaded from saved files
  # sample_n_list
  # count_df
  # core_deseq_meta_df

  deseq2_clust_dir <- glue(
    "{wkdir}/data/interim/DESeq2/subsampling/",
    "gated_12-samples_group-run/{clust}"
  )

  # Middle loop for subsampling numbers
  for (subsamp_n in sample_n_list[[clust]]) {
    # Check if all iterations for this cluster and subsample size already exist
    all_iterations_exist <- all(sapply(1:50, function(niter) {
      deseq2_dir <- glue("{deseq2_clust_dir}/n_{subsamp_n}/iter_{niter}")
      dir.exists(deseq2_dir)
    }))

    if (all_iterations_exist) {
      message("All iterations already exist for cluster ",
        clust, " and subsample size ", subsamp_n
      )
      next
    }

    # Parallel execution of iterations
    results <- foreach::foreach(niter = 1:50, .packages = c("glue", "DESeq2", "dplyr", "Matrix.utils")) %dopar% {
      deseq2_dir <- glue("{deseq2_clust_dir}/n_{subsamp_n}/iter_{niter}")
      if (dir.exists(deseq2_dir)) {
        return(NULL)
      }

      dir.create(deseq2_dir, showWarnings = FALSE, recursive = TRUE)

      # Subsampling cell type
      subsamp_data <- subsample_celltype(
        df = count_df,
        cell_id = clust,
        subsamp_n = subsamp_n,
        seed = 42 + niter  # Use different seed for each iteration
      )

      aggr_counts_subsamp <- Matrix.utils::aggregate.Matrix(
        subsamp_data$subsamp_counts,
        groupings = subsamp_data$subsamp_meta,
        fun = "sum"
      ) %>% t()

      deseq_meta <- dplyr::distinct(subsamp_data$subsamp_meta) %>% 
        dplyr::left_join(core_deseq_meta_df, by = c("sample_id")) %>% 
        dplyr::mutate_if(is.character, as.factor) %>%
        dplyr::distinct()

      # Create DESeq2 object        
      dds <- DESeq2::DESeqDataSetFromMatrix(aggr_counts_subsamp, 
                                    colData = deseq_meta, 
                                    design = ~ group + run)

      # RUNNING DESEQ2
      tryCatch({
        dds <- DESeq2::DESeq(dds)
        saveRDS(dds, glue("{deseq2_dir}/dds_{Sys.Date()}.rds"))

        # Generate results object
        res <- DESeq2::results(dds, 
                        name = "group_wildr_vs_spf",
                        alpha = 0.05,
                        contrast = c("group", "wildr", "spf"))
        res <- DESeq2::lfcShrink(dds, 
                          coef = "group_wildr_vs_spf",
                          res=res,
                          type = "apeglm")

        saveRDS(res, glue("{deseq2_dir}/res_{Sys.Date()}.rds"))

        return(TRUE)  # Indicate successful completion
      }, error = function(e) {
        message("Error occurred in iteration ", niter,
                " for cluster ", clust,
                " and subsample ", subsamp_n,
                ": ", e$message)
        return(FALSE)  # Indicate unsuccessful completion
      })
    }
  }
}

