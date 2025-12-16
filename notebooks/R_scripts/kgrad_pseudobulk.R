require(tidyverse)
require(magrittr)
require(glue)
require(Seurat)
require(grid)
require(strex)
require(broom)
# parallelization
require(BiocParallel)
require(future)
require(furrr)
require(foreach)
require(doParallel)
# plotting
require(Matrix)
# Load libraries
require(cowplot)
require(Matrix.utils)
require(edgeR)
require(Matrix)
require(reshape2)
require(S4Vectors)
require(SingleCellExperiment)
require(pheatmap)
require(apeglm)
require(png)
require(DESeq2)
require(data.table)
require(DEGreport)


# setting paths
homedir <- "/central/groups/MazmanianLab/joeB"
wkdir <- glue("{homedir}/WILDRxSPF_brains")
source(glue("{wkdir}/notebooks/R_scripts/_misc_functions.R"))


# Functions for analysis
subsample_celltype <- function(
    df, cell_id, subsamp_n, seed,
    meta_cols = c("sample_id", "cluster_id", "group")) {
    require(magrittr)
    require(dplyr)
    require(tidyselect)
    res <- list()

    # Counting cells per slice for the specified cell type
    slice_counts <- dplyr::filter(df, cluster_id == cell_id) %>%
        dplyr::group_by(sample_id) %>%
        dplyr::summarize(n = n())

    # Filtering for slices with at least the specified number of cells
    valid_slices <- dplyr::filter(slice_counts, n >= subsamp_n) %>%
        dplyr::pull(sample_id)

    set.seed(seed)
    subsamp_df <- df %>%
        dplyr::filter(cluster_id == cell_id, sample_id %in% valid_slices) %>%
        dplyr::group_by(sample_id) %>%
        dplyr::sample_n(size = subsamp_n, replace = FALSE) %>%
        dplyr::ungroup(sample_id)

    res[["subsamp_counts"]] <- subsamp_df %>% dplyr::select(-all_of(meta_cols))
    res[["subsamp_meta"]] <- subsamp_df %>%
        dplyr::select(c(sample_id))
    return(res)
}


kgrad_pseudobulk <- function(
    wkdir, clust, nthreads, sample_n_list_path, count_df_path, core_meta_path,
    deseq2_clust_dir) {
    require(glue)
    message(glue("Starting kgrad_pseudobulk analysis for cluster: {clust}"))
    source(glue("{wkdir}/notebooks/R_scripts/kgrad_pseudobulk.R"))
    dir.create(deseq2_clust_dir, showWarnings = FALSE, recursive = TRUE)

    message("Loading input data...")
    sample_n_list <- readRDS(sample_n_list_path)
    count_df <- readRDS(count_df_path)
    core_deseq_meta_df <- readRDS(core_meta_path)
    message(glue("Loaded data successfully. Number of samples in count_df: {ncol(count_df)}"))

    # Middle loop for subsampling numbers
    message(glue("Setting up parallel processing with {nthreads} threads"))
    doParallel::registerDoParallel(cores = nthreads)
    for (subsamp_n in sample_n_list[[ as.character(clust) ]]) {
        message(glue("\nProcessing subsample size: {subsamp_n}"))
        # Check if all iterations for this cluster and subsample size already exist
        all_iterations_exist <- all(sapply(1:50, function(niter) {
            deseq2_dir <- glue("{deseq2_clust_dir}/n_{subsamp_n}/iter_{niter}")
            dir.exists(deseq2_dir)
        }))

        if (all_iterations_exist) {
            message(
                "All iterations already exist for cluster ",
                clust, " and subsample size ", subsamp_n
            )
            next
        }

        message(glue("Starting parallel execution for {subsamp_n} cells"))
        # Parallel execution of iterations
        results <- foreach::foreach(
            niter = 1:50, .packages =
                c("glue", "DESeq2", "dplyr", "Matrix.utils")
        ) %dopar% {
            deseq2_dir <- glue("{deseq2_clust_dir}/n_{subsamp_n}/iter_{niter}")
            if (dir.exists(deseq2_dir)) {
                message(glue("Iteration {niter} already exists, skipping"))
                return(NULL)
            }

            message(glue("Starting iteration {niter}"))
            dir.create(deseq2_dir, showWarnings = FALSE, recursive = TRUE)

            # Subsampling cell type
            message(glue("Subsampling cells for iteration {niter}"))
            subsamp_data <- subsample_celltype(
                df = count_df,
                cell_id = clust,
                subsamp_n = subsamp_n,
                seed = 42 + niter # Use different seed for each iteration
            )

            message(glue("Aggregating counts for iteration {niter}"))
            aggr_counts_subsamp <- Matrix.utils::aggregate.Matrix(
                subsamp_data$subsamp_counts,
                groupings = subsamp_data$subsamp_meta,
                fun = "sum"
            ) %>% t()

            message(glue("Preparing metadata for iteration {niter}"))
            deseq_meta <- dplyr::distinct(subsamp_data$subsamp_meta) %>%
                dplyr::left_join(core_deseq_meta_df, by = c("sample_id")) %>%
                dplyr::mutate_if(is.character, as.factor) %>%
                dplyr::distinct()

            # Create DESeq2 object
            message(glue("Creating DESeq2 object for iteration {niter}"))
            dds <- DESeq2::DESeqDataSetFromMatrix(aggr_counts_subsamp,
                colData = deseq_meta,
                design = ~group
            )

            # RUNNING DESEQ2
            tryCatch(
                {
                    message(glue("Running DESeq2 analysis for iteration {niter}"))
                    dds <- DESeq2::DESeq(dds)
                    saveRDS(dds, glue("{deseq2_dir}/dds_{Sys.Date()}.rds"))

                    # Generate results object
                    message(glue("Generating results for iteration {niter}"))
                    res <- DESeq2::results(dds,
                        contrast = c("group", "WildR", "SPF") # WildR is Numerator
                    )
                    saveRDS(res, glue("{deseq2_dir}/res_{Sys.Date()}.rds"))

                    message(glue("Performing LFC shrinkage for iteration {niter}"))
                    res_lfcShrink <- DESeq2::lfcShrink(dds,
                        coef = "group_WildR_vs_SPF",
                        res = res,
                        type = "apeglm"
                    )
                    saveRDS(res_lfcShrink, glue("{deseq2_dir}/res_lfcShrink_{Sys.Date()}.rds"))

                    message(glue("Successfully completed iteration {niter}"))
                    return(TRUE) # Indicate successful completion
                },
                error = function(e) {
                    message(
                        "Error occurred in iteration ", niter,
                        " for cluster ", clust,
                        " and subsample ", subsamp_n,
                        ": ", e$message
                    )
                    return(FALSE) # Indicate unsuccessful completion
                }
            )
        }
        message(glue("Completed all iterations for subsample size {subsamp_n}"))
    }
    message(glue("Finished kgrad_pseudobulk analysis for cluster: {clust}"))
}

#' Function to aggregated interations of DESeq2 subsamples
process_celltype_results <- function(gated_celltype,
                                   structured_res_paths,
                                   complied_dir,
                                   analysis_date = Sys.Date()) {
  require(glue)
  require(magrittr)

  # Check if the file already exists
  output_file <- glue("{complied_dir}/{gated_celltype}_res_{analysis_date}.rds")
  if (file.exists(output_file)) {
    message(glue("File already exists for {gated_celltype}. Skipping..."))
    return(NULL)
  }
  message(gated_celltype)
  
  res_df <- structured_res_paths[[gated_celltype]] %>% 
    purrr::set_names(strex::str_after_nth(., "/", -4)) %>%
    purrr::map(
      ~ readRDS(.x) %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "gene_symbol")
    ) %>%
    dplyr::bind_rows(.id = "info") %>%
    tidyr::separate(info, 
      into = c("celltype", "subsamp_n", "iter", "date"), 
      sep = "/", 
      remove = FALSE
    ) %>%
    dplyr::mutate(subsamp_n = as.numeric(gsub("n_", "", subsamp_n))) %>%
    dplyr::mutate(subsamp_n = as.numeric(round(subsamp_n)))

  saveRDS(res_df, glue("{complied_dir}/{gated_celltype}_res_{analysis_date}.rds"))
}


#------------------------------------------------------------------------------
#                        Visualizing functions
#------------------------------------------------------------------------------
# Create line plots for p-value, adjusted p-value, and log2FC
plot_metrics <- function(data, y_var, y_label, title, yint) {
    data %>%
        ggplot(aes(
            x = subsamp_n, y = !!sym(y_var),
            color = HIT, alpha = HIT, group = gene_symbol
        )) +
        geom_line(data = data %>% filter(HIT == FALSE), linewidth = 0.5) +
        geom_line(data = data %>% filter(HIT == TRUE), linewidth = 0.5) +
        geom_hline(
            yintercept = yint,
            # yintercept = ifelse(y_var == "median_log2fc", 0, -log10(0.05)),
            linetype = "dashed", color = "black"
        ) +
        scale_color_manual(values = c("TRUE" = "darkred", "FALSE" = "grey")) +
        scale_alpha_manual(values = c("TRUE" = 0.8, "FALSE" = 0.5)) +
        labs(
            x = "Subsample Size",
            y = y_label,
            title = title,
            color = "HIT"
        ) +
        theme_bw() +
        theme(legend.position = "none")
}


