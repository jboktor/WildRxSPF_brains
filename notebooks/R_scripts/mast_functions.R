# Gene inclusion filter (expressed in at least 10% of cells)
# Cell type filter - must have at least 50 WildR and 50 SPF cells
# If there are at least 2k cells per group, randomly subsample to 2k cells per group

# "24 MY Glut"
# celltype <- "CS20230722_CLAS_24"
# tissue <- "HYP"
# cell_type_column <- "class_label"

run_mast <- function(cell_type_column, celltype, tissue)  {
    get_time <- function() print(format(Sys.time(), "%Y-%m-%d_%H:%M:%S"))

    require(MAST)
    require(dplyr)
    require(data.table)
    require(Seurat)
    require(glue)
    require(magrittr)

    wkdir <- "/resnick/groups/MazmanianLab/jboktor/WILDRxSPF_brains"
    data_path <- glue("{wkdir}/data/interim")
    seur_path <- glue("{data_path}/seurat/snRNASeq/seurat_obj_onlyparsefilters_celltyped_2025-12-27.rds")
    seurat_obj <- readRDS(seur_path)

    # filter out quality cells only
    seur_analysis <- seurat_obj %>% subset(subset = keep_cell)

    # Filter CELLS using metadata directly
    message(glue("{get_time()} Filtering cells"))
    cells_to_keep <- seur_analysis@meta.data[[cell_type_column]] == celltype & 
                    seur_analysis@meta.data$tissue == tissue
    seur_analysis <- seur_analysis[, cells_to_keep]
    message(glue("{get_time()} Number of cells in seur: {ncol(seur_analysis)}"))
    
    # check if either group has more than 2000 cells
    # subsample down to 2000 cells per group if necessary
    message(glue("{get_time()} Checking cell counts per group"))
    cells_per_group <- table(seur_analysis@meta.data$microbiome) %>% print()

    if (any(cells_per_group > 2000)) {
      groups_to_subsample <- which(cells_per_group > 2000) %>% names()
      message(glue("{get_time()} Groups to subsample: {groups_to_subsample}"))
      for (group in groups_to_subsample) {
        cells_in_group <- seur_analysis@meta.data$microbiome == group
        # selecting sampled cells from the currentgroup
        set.seed(42)
        cells_in_group_sampled <- seq_along(cells_in_group) %in% sample(which(cells_in_group), 2000)
        # selecting sampled cells + cells from the other group to subset
        selection_vector <- cells_in_group_sampled | !cells_in_group
        seur_analysis <- seur_analysis[, selection_vector]
      }
    } else {
      message(glue("{get_time()} No need to subsample down cells"))
    }
    message(glue("{get_time()} Number of cells in seur after subsampling: {ncol(seur_analysis)}"))
        
    # generate a log standardized nFeature_RNA score
    message(glue("{get_time()} Generating log standardized nFeature_RNA score"))
    seur_analysis$log_std_nFeature_RNA <- scale(log10(seur_analysis@meta.data$nFeature_RNA + 1))

    # cpm normalization
    message(glue("{get_time()} Normalizing data"))
    seur_analysis <- Seurat::NormalizeData(seur_analysis, 
                                normalization.method = "RC", 
                                scale.factor = 1e6,
                                verbose = TRUE)

    # filter GENES using criteria (protein coding and thresholded prevalence)
    message(glue("{get_time()} Filtering genes"))
    gene_meta <- seur_analysis[["RNA"]][[]] %>% glimpse()
    mat <- SeuratObject::GetAssayData(seur_analysis, assay = "RNA", layer = "counts")
    gene_prevalence <- Matrix::rowSums(mat > 0) / ncol(mat)
    gene_thres <- gene_prevalence >= 0.1
    genes_to_test <- gene_meta[gene_thres, ] %>% 
        dplyr::filter(gene_type == "protein_coding") %>% 
        rownames()

    message(glue("{get_time()} Filtered from {nrow(seur_analysis)} to {length(genes_to_test)} genes"))
    seur_analysis <- seur_analysis[genes_to_test, ]

    # pulling out the cpm matrix
    cpm_matrix <- SeuratObject::GetAssayData(seur_analysis, layer = "data")
    # Log2 scaling CPM gene counts
    message(glue("{get_time()} Log2 scaling CPM matrix"))
    log2_cpm_matrix <- log2(cpm_matrix + 1)

    # assembling the SingleCellAssay object
    sca <- MAST::FromMatrix(
        exprsArray = as.matrix(log2_cpm_matrix), 
        cData = as.data.frame(seur_analysis@meta.data)
    )
    # fit MAST models 
    message(glue("{get_time()} Fitting MAST models\n"))
    zlm <- MAST::zlm(~ microbiome + log_std_nFeature_RNA, sca)

    # LRT analysis 
    message(glue("{get_time()} Performing LRT analysis"))
    summary_res <- summary(zlm, doLRT = "microbiomeWildR")
    summary_res %>% glimpse()

    summary_resDt <- summary_res$datatable
        mast_res_df <- merge(
            summary_resDt[
                contrast=='microbiomeWildR' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values    
            summary_resDt[
                contrast=='microbiomeWildR' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], #logFC coefficients
                by='primerid')

    mast_res_df[,p_adjust:=p.adjust(`Pr(>Chisq)`, 'BH')]
    output_path = glue(
            "/resnick/groups/MazmanianLab/jboktor/WILDRxSPF_brains/data/interim/",
            "mast/snRNASeq/{cell_type_column}/",
            "{celltype}_{tissue}_{Sys.Date()}.rds")

    dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
    message(glue("{get_time()} Saving results to {output_path}"))
    saveRDS(mast_res_df, output_path)
}
