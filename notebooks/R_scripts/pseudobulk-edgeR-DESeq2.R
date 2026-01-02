
run_pseudobulk_de <- function(cell_type_column, group_column, tissue, celltype) {
    require(Seurat)
    require(glue)
    require(dplyr)
    require(magrittr)
    require(Matrix.utils)
    require(DESeq2)
    require(edgeR)
    
    get_time <- function() print(format(Sys.time(), "%Y-%m-%d_%H:%M:%S"))
    message(glue("{get_time()} Processing {celltype} {tissue}"))
    future::plan("multisession", workers = 4)
    # future::plan("sequential")
    
    # read in seurat object
    data_path <- "/resnick/groups/MazmanianLab/jboktor/WILDRxSPF_brains/data/interim"
    seur_path <- glue("{data_path}/seurat/snRNASeq/seurat_obj_onlyparsefilters_celltyped_2025-12-27.rds")
    seurat_obj <- readRDS(seur_path)
    seur_analysis <- seurat_obj %>% subset(subset = keep_cell)

    # create output directories
    edgeR_dir <- glue("{data_path}/edgeR/snRNASeq/{cell_type_column}/{celltype}_{tissue}")
    dir.create(edgeR_dir, showWarnings = FALSE, recursive = TRUE)
    deseq2_dir <- glue("{data_path}/DESeq2/snRNASeq/{cell_type_column}/{celltype}_{tissue}")
    dir.create(deseq2_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Subset CELLS for celltype and tissue
    message(glue("{get_time()} Filtering cells"))
    cells_to_keep <- seur_analysis@meta.data[[cell_type_column]] == celltype & 
                    seur_analysis@meta.data$tissue == tissue
    seur_analysis <- seur_analysis[, cells_to_keep]
    message(glue("{get_time()} Number of cells in seur: {ncol(seur_analysis)}"))
    
    # Subset GENES for protein coding only
    message(glue("{get_time()} Filtering genes for protein coding only"))
    gene_meta <- seur_analysis[["RNA"]][[]] %>% glimpse()
    mat <- SeuratObject::GetAssayData(seur_analysis, assay = "RNA", layer = "counts")
    gene_prevalence <- Matrix::rowSums(mat > 0) / ncol(mat)
    gene_thres <- gene_prevalence >= 0.1
    genes_to_test <- gene_meta[gene_thres, ] %>% 
        dplyr::filter(gene_type == "protein_coding") %>% 
        rownames()

    message(glue("{get_time()} Filtered from {nrow(seur_analysis)} to {length(genes_to_test)} genes"))
    seur_analysis <- seur_analysis[genes_to_test, ]
    
    # ------------------------------
    # Define pseudobulk grouping (by sample)
    # ------------------------------
    grouping <- seur_analysis@meta.data %>%
        as.data.frame() %>%
        dplyr::select(sample_id = sample) %>%
        glimpse()
    analysis_meta <- seur_analysis@meta.data %>% 
        dplyr::select(sample_id = sample, group = {{group_column}}) %>%
        dplyr::mutate_if(is.character, as.factor) %>%
        dplyr::distinct() %>% 
        tibble::remove_rownames() 
    
    # ------------------------------
    # Aggregate pseudobulk counts (sum across cells)
    # ------------------------------
    aggr_count_matrix <- Matrix.utils::aggregate.Matrix(
        t(seur_analysis@assays$RNA@counts),  # cells × genes
        groupings = grouping,
        fun = "sum"
        ) %>% t()  # genes × pseudobulk samples

    # ------------------------------
    # filter GENES using criteria (thresholded prevalence) / variance 
    # ------------------------------    
    min_counts <- 10
    min_samples <- 1
    keep_genes <- rowSums(aggr_count_matrix >= min_counts) >= min_samples
    aggr_count_matrix_filt <- aggr_count_matrix[keep_genes, ]

    # Filter by variance (drop bottom 30%)
    gene_vars <- apply(aggr_count_matrix_filt, 1, var)
    var_threshold <- quantile(gene_vars, 0.3)
    aggr_count_matrix_filt <- aggr_count_matrix_filt[gene_vars > var_threshold, ]

    message(glue("Filtered from {nrow(aggr_count_matrix)} to {nrow(aggr_count_matrix_filt)} genes"))

    # ------------------------------
    # Create edgeR DGEList object
    # ------------------------------
    message(glue("{get_time()} Creating edgeR DGEList object"))
    dge <- edgeR::DGEList(counts = aggr_count_matrix_filt)
    dge$samples$group <- analysis_meta$group
    dge$samples$sample_id <- analysis_meta$sample_id

    message(glue("{get_time()} Calculating TMM normalization factors"))
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    design <- model.matrix(~ group, data = analysis_meta)
    dge <- edgeR::estimateDisp(dge, design, robust = TRUE)

    # ------------------------------
    # Fit GLM with QL (quasi-likelihood) and test
    # ------------------------------
    message(glue("{get_time()} Fitting GLM with QL (quasi-likelihood)"))
    fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
    qlf <- edgeR::glmQLFTest(fit, coef = "groupWildR")

    res <- edgeR::topTags(qlf, n = Inf)$table
    res <- res %>% dplyr::arrange(FDR) %>% glimpse()

    # ------------------------------
    # Fit GLM with Vanilla model with LRT
    # ------------------------------
    message(glue("{get_time()} Fitting GLM with Vanilla model with LRT"))
    fit <- edgeR::glmFit(dge, design, robust = TRUE)
    lrt <- edgeR::glmLRT(fit, coef = "groupWildR")

    res <- edgeR::topTags(lrt, n = Inf)$table
    res <- res %>% dplyr::arrange(FDR) %>% glimpse()

    message(glue("{get_time()} Saving edgeR results"))
    saveRDS(res, glue("{edgeR_dir}/res_ql_{Sys.Date()}.rds"))
    saveRDS(res, glue("{edgeR_dir}/res_lrt_{Sys.Date()}.rds"))
    saveRDS(dge, glue("{edgeR_dir}/dge_{Sys.Date()}.rds"))

    # ------------------------------
    # Create DESeq2 object
    # ------------------------------
    message(glue("\n\nCreating DESeq2 object for iteration {celltype}_{tissue}"))
    dds <- DESeq2::DESeqDataSetFromMatrix(aggr_count_matrix_filt,
        colData = analysis_meta, design = ~ group)
    dds <- DESeq2::DESeq(dds)
    saveRDS(dds, glue("{deseq2_dir}/dds_{Sys.Date()}.rds"))

    # Extract and inspect size factors
    message(glue("{get_time()} Extracting and inspecting size factors"))
    size_factors <- DESeq2::sizeFactors(dds)
    analysis_meta$size_factor <- size_factors

    # Generate results object
    message(glue("{get_time()} Generating DESeq2 results"))
    res_deseq2 <- DESeq2::results(dds,
        contrast = c("group", "WildR", "SPF") # WildR is Numerator
    )
    message(glue("{get_time()} Performing LFC shrinkage"))
    res_deseq2_lfcShrink <- DESeq2::lfcShrink(dds,
        coef = "group_WildR_vs_SPF",
        res = res_deseq2,
        type = "apeglm"
    )
    message(glue("{get_time()} Saving DESeq2 results"))
    message(glue("\n\nSaving DESeq2 results for {celltype}_{tissue}"))
    saveRDS(res_deseq2, glue("{deseq2_dir}/res_{Sys.Date()}.rds"))
    saveRDS(res_deseq2_lfcShrink, glue("{deseq2_dir}/res_lfcShrink_{Sys.Date()}.rds")) 
}
