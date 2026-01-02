

#' Pseudobulk analysis with edgeRv4 --
#'  Spliting cells into NMF dims - 
#' for each osNMF dimension and then merging the results
#' 
#' 





get_dds_LRTresults <- function(clustx) {
  
  deseq2_dir <- glue("{wkdir}/figures/DESeq2/{Sys.Date()}/{clustx}")
  dir.create(deseq2_dir, showWarnings = FALSE, recursive = TRUE)
  
  print(clustx) # useful for debugging
  
  # Extract counts matrix and metadata for cluster x
  idx <- which(names(counts_ls) == clustx)
  cluster_counts <- counts_ls[[idx]]
  cluster_metadata <- metadata_ls[[idx]]
  
  # Print error message if sample names do not match
  if ( all(colnames(cluster_counts) != rownames(cluster_metadata)) ) {
    print("ERROR: sample names in counts matrix columns and metadata rows do not match!")
  }
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ group + run + bregma)
  dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)
  
  # Extract results  
  res_LRT <- results(dds_lrt, 
    name = "group_wildr_vs_spf",
    alpha = 0.05
    )
  # Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
  res_LRT <- lfcShrink(dds_lrt, 
    coef = "group_wildr_vs_spf",
    res=res_LRT,
    type = "ashr"
    )

  # Create a tibble for LRT results
  res_LRT_tb <- res_LRT %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>% 
    as_tibble()
  
  # Save all results
  if (!dir.exists(deseq2_dir)) { dir.create(deseq2_dir) }
  write.csv(res_LRT_tb,
            glue("{deseq2_dir}/{clustx}_LRT_all_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  # Subset to return genes with padj < 0.05
  sigLRT_genes <- res_LRT_tb %>% 
    filter(padj < 0.05)
  
  # Save significant results
  write.csv(sigLRT_genes,
    glue("{deseq2_dir}/{clustx}_LRT_signif_genes.csv"),
    quote = FALSE, 
    row.names = FALSE
    )
  
  # Transform counts for data visualization
  rld <- rlog(dds_lrt, blind = TRUE)
  
  # Extract the rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)

    # QC plots
  for (pca_group in c("group", "run", "cell_count")) {
    p_pca <- DESeq2::plotPCA(rld, ntop = 500, intgroup = pca_group)
    ggsave(
      glue("{deseq2_dir}/PCA_{pca_group}_{Sys.Date()}.png"),
      p_pca,
      width = 7, height = 5
    )
  }
  
  # Obtain rlog values for those significant genes
  cluster_rlog <- rld_mat[sigLRT_genes$gene, ]
  cluster_meta_sig <- cluster_metadata[which(rownames(cluster_metadata) %in% colnames(cluster_rlog)), ]
  
  # # Use the `degPatterns` function from DEGreport package to show gene clusters across sample groups
  # cluster_groups <- degPatterns(cluster_rlog, 
  #   metadata = cluster_meta_sig,
  #   time = "group", col = NULL
  #   )
  # ggsave(
  #   glue("{deseq2_dir}/{clustx}_LRT_DEgene_groups.png")
  # )
  
  # Save what is stored in the `df` component
  write.csv(cluster_groups$df,
    glue("{deseq2_dir}/{clustx}_LRT_DEgene_groups.csv"),
    quote = FALSE,
    row.names = FALSE
    )
  
  saveRDS(
    cluster_groups,
    glue("{deseq2_dir}/{clustx}_LRT_DEgene_groups.rds")
    )
  
  save(dds_lrt, cluster_groups, res_LRT, sigLRT_genes, 
       file = glue("{deseq2_dir}/{clustx}_all_LRTresults.Rdata")
       )
  
    # Subset the significant results
  # sig_res <- dplyr::filter(res_LRT_tb, padj < padj_cutoff) %>%
  #   dplyr::arrange(padj)

  # Volcano plot
  res_table_thres <- res_LRT_tb[!is.na(res_LRT_tb$padj), ] %>% 
    mutate(threshold = padj < padj_cutoff & abs(log2FoldChange) >= log2fc_cutoff)

  p_volcano <- ggplot(res_table_thres) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
    labs(x = "log2 fold change", y = "-log10 adjusted p-value") +
    scale_color_manual(values = c("grey60", "red3")) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.3), hjust = 0.5),
          axis.title = element_text(size = rel(1.15)))

  ggsave(
    glue("{deseq2_dir}/volcano_plot_{Sys.Date()}.png"),
    p_volcano,
    width = 4, height = 4
  )
}

# purrr::map(cluster_names, get_dds_LRTresults)






# Psuedo bulk analysis


# Load libraries
library(tidyverse)
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
library(RColorBrewer)
library(data.table)
library(glue)

merged_rois <- readRDS(
    glue("{wkdir}/data/interim/merged_roi_seurat_filtered_2024-07-22.rds")
)

# merged_rois$osNMF_2_postqc %>% hist()
# merged_rois$osNMF_2_postqc_Q10 %>% table()



# Extract raw counts and metadata to create SingleCellExperiment object
counts <- merged_rois@assays$SCT@counts %>% glimpse()
metadata <- merged_rois@meta.data %>% glimpse()

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(merged_rois@meta.data$singleR_labels)

# Create single cell experiment object
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = metadata
)

# ## Check the counts matrix
# dim(counts(sce))
# counts(sce)[1:6, 1:6]
# dim(colData(sce))
# head(colData(sce))

cluster_names <- unique(colData(sce)$singleR_labels)
cluster_names
length(cluster_names)

# Extract unique names of samples (= levels of sample_id factor variable)
sample_names <- unique(colData(sce)$slice_id)
sample_names

# Total number of samples
length(sample_names)


# Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
groups <- colData(sce)[, c("cluster_id", "slice_id")]
head(groups)

# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts <- aggregate.Matrix(
  t(counts(sce)), 
  groupings = groups, fun = "sum"
  ) %>% 
  t()

aggr_counts[1:6, 1:6]
dim(aggr_counts)


# Loop over all cell types to extract corresponding counts, and store information in a list
## Initiate empty list
counts_ls <- list()

for (i in 1:length(cluster_names)) {
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])

  ## Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
}

glimpse(counts_ls)


# Reminder: explore structure of metadata
head(colData(sce))

# Extract sample-level variables
metadata <- colData(sce) %>% 
  as.data.frame() %>% 
  dplyr::select(group, run, roi, bregma, slice_id) %>%
  distinct() %>%
  glimpse()

dim(metadata)
head(metadata)

rownames(metadata) <- metadata$slice_id
head(metadata)

# Number of cells per sample and cluster
t <- table(colData(sce)$slice_id,
           colData(sce)$cluster_id)
t[1:6, 1:6]

# Creating metadata list
## Initiate empty list
metadata_ls <- list()

for (i in 1:length(counts_ls)) {
  
    ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
    df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
    
    ## Use tstrsplit() to separate cluster (cell type) and sample IDs
    df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
    df$sample_id  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
    
    
    ## Retrieve cell count information for this cluster from global cell count table
    idx <- which(colnames(t) == unique(df$cluster_id))
    cell_counts <- t[, idx]
    
    ## Remove samples with zero cell contributing to the cluster
    cell_counts <- cell_counts[cell_counts > 0]
    
    ## Match order of cell_counts and sample_ids
    sample_order <- match(df$sample_id, names(cell_counts))
    cell_counts <- cell_counts[sample_order]
    
    ## Append cell_counts to data frame
    df$cell_count <- cell_counts
    
    
    ## Join data frame (capturing metadata specific to cluster) to generic metadata
    df %<>% dplyr::left_join(metadata, by = c("sample_id" = "slice_id"))
    
    ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
    rownames(df) <- df$cluster_sample_id
    
    ## Store complete metadata for cluster i in list
    metadata_ls[[i]] <- df
    names(metadata_ls)[i] <- unique(df$cluster_id)

  # Create DESeq2 object        
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ group + run + bregma)

  # Transform counts for data visualization
  rld <- rlog(dds, blind=TRUE)
}

# Explore the different components of the list
str(metadata_ls)


# ---------- DESEQ ANALYSIS
library(DEGreport)

# Double-check that both lists have same names
all(names(counts_ls) == names(metadata_ls))


