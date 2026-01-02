

# Cell typing --  clustermole marker genes
#______________________________________________________________________________
# devtools::install_github("igordot/clustermole")

library(clustermole)
merged_rois <- readRDS(
    glue("{wkdir}/data/interim/2024-04-12_merged_roi_seurat.rds")
)

cluster_markers_conserved <- readRDS(
    glue(
      "{wkdir}/data/interim/",
      "2024-04-12_cluster_markers_conserved_list.rds"
    )
)
markers <- clustermole_markers(species = "mm") %>%
  filter(species == "Mouse") %>% 
  filter(organ == "Brain") %>% 
  glimpse

markers_list <- split(x = markers$gene, f = markers$celltype_full)
names(markers_list)



# Note species identifier doesn't work properly for package,
# using hs instead of mm
cluster_alignments <- names(cluster_markers) %>%
  purrr::set_names() %>%
  purrr::map(
    ~ clustermole_overlaps(
      genes = cluster_markers_conserved[[.]] %>% 
      filter(spf_avg_log2FC > 0) %>% 
      pull(gene_symbol), species = "hs"
    ) %>%
      filter(species == "Mouse") %>%
      filter(organ == "Brain")
      # filter(db != "PanglaoDB")
        # slice_min(order_by = fdr, n = 1) %>% 
      # pull("celltype")
  )

cluster_alignments
#______________________________________________________________________________




library(tidyverse)


df <- data.frame(
  x=runif(100),
  y=runif(100)
)
df %>% ggplot(aes(x = x, y = y)) +
  geom_point()



merged_rois <- readRDS(
    glue("{wkdir}/data/interim/2024-01-23_merged_roi_umap_seurat.rds")
)

merged_rois

seqfish_obj <- BuildNicheAssay(
    object = merged_rois,
    fov = "RNA",
    group.by = "cluster_labels",
    niches.k = 5,
    neighbors.k = 30
)


seqfish_obj <- BuildNicheAssay(
    object = merged_rois,
    fov = "Spatial",
    group.by = "cluster_labels",
    niches.k = 5,
    neighbors.k = 30
)




coords <- data.frame(
  x = merged_rois$center_x,
  y = merged_rois$center_y
)

library(Seurat)


# Assuming merged_rois is your main Seurat object
# and it contains metadata with 'center_x' and 'center_y' for spatial coordinates

# First, create a data frame with the spatial coordinates
spatial_coords <- data.frame(
  x = merged_rois@meta.data$center_x,
  y = merged_rois@meta.data$center_y
)

# Assuming you have a matrix of gene expression data for the spatial assay
# The rows are genes, and the columns are cells/spots
# Replace this with your actual spatial expression data
spatial_expression_data <- <YourSpatialData>

# Create a new Assay object with the spatial expression data
spatial_assay <- CreateAssayObject(counts = spatial_expression_data)

# Add the spatial coordinates to the metadata of the new assay
spatial_assay <- AddMetaData(spatial_assay, metadata = spatial_coords)

# Add the spatial assay to your Seurat object
merged_rois[['spatial']] <- spatial_assay

# Now your Seurat object has a new assay named 'spatial' with the spatial coordinates in its metadata
# You can try to run BuildNicheAssay with the updated Seurat object



seqfish_obj <- BuildNicheAssay(
  object = roi1,
  fov = "RNA",  # Referencing the new spatial assay
#   group.by = "cluster_labels",
  niches.k = 5,
  neighbors.k = 30
)

# Note: This is a hypothetical code and may need adjustments based on the actual data and functions available in Seurat and BuildNicheAssay


names(roi1)
roi1 %>% str
roi1$x
roi1$y



roi$FindNeighbors.SCT.pca:Formal
roi1$SCT_snn_res.0.3



merged_cell_data <- tibble(
  label = roi1$label,
  x = roi1$x,
  y = roi1$y,
  area = roi1$area,
#   roi = roi1$roi,
#   group = roi1$group,
  cluster_labels = roi1$SCT_snn_res.0.3
)

merged_cell_data %>% 
  ggplot(aes(x=x, y=-y)) +
  geom_point(aes(color = cluster_labels), size = 0.1, alpha = 0.6) +
#   scale_color_manual(values = color_pal) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  coord_fixed() +
  labs(color = "Cell Type") +
  theme_minimal()




merged_rois

# we will use data from 2 technologies for the reference
brain_ref_raw <- readRDS(glue("{wkdir}/data/input/allen_mop_2020.rds"))
# pre-process dataset (without integration)
brain_ref <- NormalizeData(brain_ref)
brain_ref <- FindVariableFeatures(brain_ref)
brain_ref <- ScaleData(brain_ref)
brain_ref <- RunPCA(brain_ref)
brain_ref <- FindNeighbors(brain_ref, dims = 1:30)
brain_ref <- FindClusters(brain_ref)
brain_ref <- RunUMAP(brain_ref, dims = 1:30)

DimPlot(brain_ref, group.by = c("cluster"))


brain_ref <- IntegrateLayers(
  object = brain_ref,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)
brain_ref <- FindNeighbors(brain_ref, reduction = "integrated.cca", dims = 1:30)
brain_ref <- FindClusters(brain_ref)









# Extract metadata
metadata <- brain_ref_raw@meta.data

# Extract expression data for a specific gene, for example "GeneX"
gene_expression <- FetchData(brain_ref_raw, vars = "Lhfpl3")

# Combine metadata with gene expression
combined_data <- cbind(metadata, gene_expression)

# Use ggplot to create the plot
ggplot(combined_data, aes(x = class, y = Lhfpl3)) + 
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Expression of GeneX in Different Cell Types",
       x = "Cell Type",
       y = "Expression Level")








#______________________________________________________________________________
# Integration tutorial

library(Seurat)
library(SeuratData)
library(patchwork)

# install dataset
InstallData("ifnb")

# load dataset
ifnb <- LoadData("ifnb")
# split the RNA measurements into two layers one for control cells, one for stimulated cells

ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb

# run standard anlaysis workflow
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb)

ifnb <- RunUMAP(ifnb,
  dims = 1:30, reduction = "pca",
  reduction.name = "umap.unintegrated"
)
DimPlot(ifnb,
  reduction = "umap.unintegrated",
  group.by = c("stim", "seurat_annotations")
)


ifnb <- IntegrateLayers(
  object = ifnb, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

# re-join layers after integration
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")
# Visualization
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations")) #+
  theme(legend.position = "none")


Idents(ifnb) <- "seurat_annotations"
nk.markers <- FindConservedMarkers(
  ifnb,
  ident.1 = "NK", 
  grouping.var = "stim", 
  verbose = FALSE
)
nk.markers %>% glimpse


# NEEDS TO BE FIXED AND SET ORDER CORRECTLY
Idents(ifnb) <- factor(Idents(ifnb), levels = c("pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono",
    "B Activated", "B", "CD8 T", "NK", "T activated", "CD4 Naive T", "CD4 Memory T"))

markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
    "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
    "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
DotPlot(ifnb, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") +
    RotatedAxis()


library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

aggregate_ifnb <- AggregateExpression(ifnb,
  group.by = c("seurat_annotations", "stim"), return.seurat = TRUE
)
genes.to.label = c(
  "ISG15", "LY6E", "IFI6",
  "ISG20", "MX1", "IFIT2", 
  "IFIT1", "CXCL10", "CCL8"
)

p1 <- CellScatter(aggregate_ifnb, "CD14 Mono_CTRL", "CD14 Mono_STIM", highlight = genes.to.label)
p2 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)

p3 <- CellScatter(aggregate_ifnb, "CD4 Naive T_CTRL", "CD4 Naive T_STIM", highlight = genes.to.label)
p4 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE)

p2 + p4


ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")
Idents(ifnb) <- "celltype.stim"
b.interferon.response <- FindMarkers(ifnb, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
head(b.interferon.response, n = 15)

FeaturePlot(ifnb, features = c("CD3D", "GNLY", "IFI6"), split.by = "stim", max.cutoff = 3, cols = c("grey",
    "red"), reduction = "umap")

plots <- VlnPlot(ifnb, features = c("LYZ", "ISG15", "CXCL10"), split.by = "stim", group.by = "seurat_annotations",
    pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)


# ______________________________________________________________________________
# Integration attempt

# May need to update splitting methods (using pre-existing split which may not be right)
# ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
# ifnb

# run standard anlaysis workflow
# merged_rois <- merged_rois %>%
#   NormalizeData() %>%
#   FindVariableFeatures() %>%
#   ScaleData() %>%
#   RunPCA()

# merged_rois %<>% RunUMAP(
#   dims = 1:30, reduction = "pca",
#   reduction.name = "umap.unintegrated"
# )



# we will use data from 2 technologies for the reference
brain_ref_raw <- readRDS(glue("{wkdir}/data/input/allen_mop_2020.rds"))
# pre-process dataset (without integration)
brain_ref <- NormalizeData(brain_ref_raw)
brain_ref <- FindVariableFeatures(brain_ref)
brain_ref <- ScaleData(brain_ref)
brain_ref <- RunPCA(brain_ref)
brain_ref <- FindNeighbors(brain_ref, dims = 1:30)
brain_ref <- FindClusters(brain_ref)
brain_ref <- RunUMAP(brain_ref, dims = 1:30)





merged_rois <- readRDS(
    glue("{wkdir}/data/interim/2024-01-23_merged_roi_umap_seurat.rds")
)

# Join data and resplit by ref / data
merged_rois[["RNA"]] <- JoinLayers(merged_rois[["RNA"]])

# Merge data with reference
merged_ref <- Merge_Seurat_List(
  list(merged_rois, brain_ref),
  add.cell.ids = c("mouse_data", "mouse_ref"),
  merge.data = TRUE,
  project = "SeuratProject_Reference"
)

# jointly format objs
merged_ref %<>%
  UpdateSeuratObject() %>%
  SCTransform(assay = "Spatial") %>%
  RunPCA() %>%
  RunUMAP(dims = 1:4) %>%
  FindNeighbors(reduction = "pca", dims = 1:4) %>%
  FindClusters(resolution = 0.3)


merged_ref <- IntegrateLayers(
  object = merged_ref,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = TRUE
)

brain_ref <- FindNeighbors(brain_ref, reduction = "integrated.cca", dims = 1:30)
brain_ref <- FindClusters(brain_ref)



merged_rois %>% DimPlot(
  reduction = "umap",
  group.by = c("roi", "group", "seurat_clusters")
)


merged_rois <- IntegrateLayers(
  object = merged_rois, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

# re-join layers after integration
merged_rois[["RNA"]] <- JoinLayers(merged_rois[["RNA"]])

merged_rois %<>%
  FindNeighbors(reduction = "integrated.cca", dims = 1:30) %>%
  FindClusters(resolution = 1) %>%
  RunUMAP(dims = 1:30, reduction = "integrated.cca")

# Visualization
DimPlot(merged_rois, reduction = "umap", group.by = c("stim", "seurat_annotations")) #+
  theme(legend.position = "none")


Idents(merged_rois) <- "seurat_annotations"
nk.markers <- FindConservedMarkers(
  merged_rois,
  ident.1 = "NK", 
  grouping.var = "stim", 
  verbose = FALSE
)
nk.markers %>% glimpse


# # NEEDS TO BE FIXED AND SET ORDER CORRECTLY
# Idents(merged_rois) <- factor(Idents(merged_rois), levels = c("pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono",
#     "B Activated", "B", "CD8 T", "NK", "T activated", "CD4 Naive T", "CD4 Memory T"))

# markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
#     "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
#     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
# DotPlot(merged_rois, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") +
#     RotatedAxis()


# library(ggplot2)
# library(cowplot)
# theme_set(theme_cowplot())

# aggregate_merged_rois <- AggregateExpression(merged_rois,
#   group.by = c("seurat_annotations", "stim"), return.seurat = TRUE
# )
# genes.to.label = c(
#   "ISG15", "LY6E", "IFI6",
#   "ISG20", "MX1", "IFIT2",
#   "IFIT1", "CXCL10", "CCL8"
# )

# p1 <- CellScatter(aggregate_merged_rois, "CD14 Mono_CTRL", "CD14 Mono_STIM", highlight = genes.to.label)
# p2 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)

# p3 <- CellScatter(aggregate_merged_rois, "CD4 Naive T_CTRL", "CD4 Naive T_STIM", highlight = genes.to.label)
# p4 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE)

# p2 + p4


# merged_rois$celltype.stim <- paste(merged_rois$seurat_annotations, merged_rois$stim, sep = "_")
# Idents(merged_rois) <- "celltype.stim"
# b.interferon.response <- FindMarkers(merged_rois, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
# head(b.interferon.response, n = 15)

# FeaturePlot(merged_rois, features = c("CD3D", "GNLY", "IFI6"), split.by = "stim", max.cutoff = 3, cols = c("grey",
#     "red"), reduction = "umap")

# plots <- VlnPlot(merged_rois, features = c("LYZ", "ISG15", "CXCL10"), split.by = "stim", group.by = "seurat_annotations",
#     pt.size = 0, combine = FALSE)
# wrap_plots(plots = plots, ncol = 1)


# #______________________________________________________________________________
# # Integration attempt







#______________________________________________________________________________

merged_rois <- readRDS(
    glue("{wkdir}/data/interim/2024-04-12_merged_roi_seurat.rds")
)

cluster_markers_conserved <- unique(merged_rois$seurat_clusters) %>%
  purrr::set_names() %>%
  purrr::map(~ Seurat::FindConservedMarkers(merged_rois,
    ident.1 = ., grouping.var = "group",
    assay = "SCT"
  ) %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene_symbol"))

saveRDS(
    cluster_markers_conserved,
    glue(
      "{wkdir}/data/interim/",
      "{Sys.Date()}_cluster_markers_conserved_list.rds"
    )
)


# p <- FeaturePlot(merged_rois, features = c("OLIG1"), min.cutoff = "q10")
# ggsave(p, file = "OLIG1.png")
# Idents(merged_rois) 



# %>% 
#   as.data.frame() %>% 
#   rownames_to_column(var = "gene_symbol") %>% 

# cluster_markers_df <- cluster_markers %>% bind_rows(.id = "cluster number")

# cluster_markers[[1]] %>%
marker_list_data <- names(cluster_markers) %>%
  purrr::set_names() %>%
  purrr::map(
    ~ cluster_markers[[.]] %>%
      slice_min(order_by = p_val, n = 100) %>%
      slice_max(order_by = abs(avg_log2FC), n = 50)
  )

# Chat GPT generated list
# cell_type_list <- list(
#   `0` = "Microglia",
#   `1` = "Microglia",
#   `2` = "Oligodendrocyte precursor cells",
#   `3` = "Microglia",
#   `4` = "Endothelial",
#   `5` = "Oligodendrocyte precursor cells",
#   `6` = "Oligodendrocytes",
#   `7` = "Oligodendrocytes",
#   `8` = "Endothelial",
#   `9` = "Oligodendrocytes",
#   `10` = "Interneurons",
#   `11` = "Macrophages",
#   `12` = "Endothelial",
#   `13` = "Microglia",
#   `14` = "Smooth muscle cell or Pericyte",
#   `15` = "Oligodendrocytes",
#   `16` = "Ependymal cell"
# )

# Manual list viewed from cluster markers
cell_type_list <- list(
  `0` = "Interneurons",
  `1` = "Interneurons",
  `2` = "Oligodendrocyte precursor cells",
  `3` = "Astrocytes",
  `4` = "Interneurons",
  `5` = "Astrocytes",
  `6` = "Interneurons",
  `7` = "Oligodendrocytes",
  `8` = "Endothelial cell",
  `9` = "Interneurons",
  `10` = "Interneurons",
  `11` = "Astrocytes",
  `12` = "Endothelial cell",
  `13` = "Microglia",
  `14` = "Bergmann glia",
  `15` = "Interneurons",
  `16` = "Neural stem cell"
)


merged_rois <- readRDS(
    glue("{wkdir}/data/interim/2024-04-12_merged_roi_seurat.rds")
)

# relabel cluster identity using chatgpt generated cluster cell type assignment 
merged_rois$cluster_labels <-
  cell_type_list[merged_rois$seurat_clusters] %>% 
  unlist() %>% 
  unname() %>% 
  as.factor()



p_umap_group <- Seurat::DimPlot(
  object = merged_rois,
  reduction = "umap",
  group.by = "cluster_labels",
  label = TRUE,
  label.size = 3,
  repel = TRUE,
  label.box = TRUE,
  label.color = "black",
  split.by = 'group',
  raster=FALSE,
  alpha = 0.7
) + 
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) +
  NoLegend() +
  scale_color_manual(values = color_pal) +
  scale_fill_manual(values = color_pal)

ggsave(
  glue("{wkdir}/figures/merged_group_UMAP_{Sys.Date()}.png"),
  p_umap_group,
  width = 9, height = 5
)




















# ______________________________________________________




# adding new metadata column for group specific cell type labels
merged_rois$singleR_labels_group <- paste0(
  merged_rois@meta.data$singleR_labels,
  "__",
  merged_rois@meta.data$group
)

# temporarily setting cluster labels to singleR labels
Idents(merged_rois) <- merged_rois@meta.data$singleR_labels_group


# Explore differentially expressed genes between 
# SPF and WildR groups within clusters

# options(future.globals.maxSize = 10 * 1024^3)  # Increase to 10 GiB
# future::plan("multisession", workers = 32)
singler_cluster_expression_comp <- 
  unique(merged_rois@meta.data$singleR_labels) %>%
  purrr::set_names() %>%
  # furrr::future_map(
  purrr::map(
    ~ Seurat::FindMarkers(
      merged_rois,
      ident.1 = glue("{.}__spf"),
      ident.2 = glue("{.}__wildr"),
      min.pct = 0.25
    ) %>%
      as.data.frame() %>%
      rownames_to_column(var = "gene_symbol") %>%
      mutate(SingleR_cluster = .x) %>%
      filter(p_val_adj <= 0.1)
  ) %>%
  bind_rows() %>%
  glimpse()
  

saveRDS(
  singler_cluster_expression_comp,
  glue(
    "{wkdir}/data/interim/",
    "{Sys.Date()}_singler_cluster_expression_comp.rds"
  )
)





singler_cluster_expression_comp <- readRDS(
  glue(
    "{wkdir}/data/interim/",
    "2024-04-15_singler_cluster_expression_comp.rds"
  )
)

singler_cluster_expression_comp %>% glimpse

singler_cluster_expression_comp %>%
  group_by(SingleR_cluster) %>%
  summarize(n = n()) %>% 
  arrange(-n) %>% 
  View


singler_cluster_expression_comp %>% glimpse

# Plotting L2FC overview
# get order of clusters with positive changes
cluster_order <- singler_cluster_expression_comp %>%
  mutate(avg_log2FC = if_else(avg_log2FC > 0, avg_log2FC, 0)) %>%
  group_by(SingleR_cluster) %>%
  summarize(sum_log2FC = sum(avg_log2FC)) %>% 
  arrange(sum_log2FC) %>% 
  pull(SingleR_cluster)


p_cluster_exp_survey <- singler_cluster_expression_comp %>%
  mutate(SingleR_cluster = factor(SingleR_cluster, levels = cluster_order)) %>%
    ggplot(aes(
      x = avg_log2FC,
      y = SingleR_cluster,
      color = gene_symbol
    )) +
    geom_bar(stat = "identity", position = "stack", fill = "white") +
     labs(y = NULL, x = expression(log[2]*"Fold Change WildR/SPF")) +
      theme_bw() +
      scale_color_viridis_d(option = "turbo") +
      theme(legend.position = "none")
    

ggsave(
  glue("{wkdir}/figures/celltype_log2FC_stacked-bars_{Sys.Date()}.png"),
  p_cluster_exp_survey,
  width = 7, height = 7
)










library(clusterProfiler)
library(org.Mm.eg.db)

gene_symbols <- singler_cluster_expression_comp$gene_symbol %>% unique()
gene_symbol_map <- data.frame(
  gene_symbol = gene_symbols,
  gene_symbol_formatted = tools::toTitleCase(tolower(gene_symbols))
)

deg_df <- singler_cluster_expression_comp %<>%
  left_join(gene_symbol_map, by = "gene_symbol") %>%
  mutate(direction = if_else(
    avg_log2FC > 0, "WildR Upregulated", "WildR Downregulated"
  ))

singler_cluster_expression_comp %>% glimpse



# format params for GO enrichment
params_1 <- expand.grid(
  unique(deg_df$SingleR_cluster),
  c( "WildR Upregulated", "WildR Downregulated")
) %>% 
dplyr::rename(
  SingleR_cluster = Var1,
  direction = Var2
)
params_2 <- expand.grid(
  c("BP", "CC", "MF"),
  c( "WildR Upregulated", "WildR Downregulated")
) %>% 
dplyr::rename(
  Ontology = Var1,
  direction = Var2
)

go_params <- params_1 %>% 
  full_join(params_2, by = "direction", relationship = "many-to-many") %>% 
  glimpse()



future::plan("multisession", workers = 16)
enrichGO_res <-
  furrr::future_pmap(
    list(
      go_params$Ontology,
      go_params$direction,
      go_params$SingleR_cluster
    ),
    ~ enrichGO(
      gene = deg_df %>%
        filter(direction == ..2 & SingleR_cluster == ..3) %>%
        pull(gene_symbol_formatted),
      OrgDb = org.Mm.eg.db,
      keyType = "SYMBOL",
      ont = ..1,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    ) %>%
      as.data.frame() %>%
      mutate(Ontology = ..1, direction = ..2, SingleR_cluster = ..3)
  ) %>%
  bind_rows() %>% 
  glimpse()

saveRDS(
  enrichGO_res,
  glue(
    "{wkdir}/data/interim/",
    "{Sys.Date()}_enrichGO_res.rds"
  )
)


enrichGO_res <- readRDS(
  glue(
    "{wkdir}/data/interim/",
    "2024-04-15_enrichGO_res.rds"
  )
)


enrichGO_res %>% head()

enrichGO_res %>% View
 
for (direc in c("WildR Upregulated", "WildR Downregulated")) {
  for (cluster in unique(enrichGO_res$SingleR_cluster)) {
    message(glue("Plotting: {cluster} - {direc}"))
    p <- enrichGO_res %>%
      filter(direction == direc) %>%
      filter(SingleR_cluster == cluster) %>%
      filter(qvalue <= 0.05) %>%
      slice_min(qvalue, n = 30) %>%
      slice_max(Count, n = 20, with_ties = FALSE) %>%
      mutate(Description = glue("{Description}:{Ontology}")) %>% 
      ggplot(aes(x = Count, y = fct_reorder(Description, Count))) +
      geom_col(width = 0.4, fill = "red") +
      theme_bw() +
      labs(
        x = "Gene Count", y = NULL,
        title = glue("{cluster}: {direc}")
      )

    ggsave(
      glue(
        "{wkdir}/figures/GO/",
        "GO_enrichment_{cluster}_{direc}_{Sys.Date()}.png"
        ),
      p,
      width = 7, height = 5.5
    )
  }
}


enrichGO_res %>% glimpse
enrichGO_res %>%
  filter(qvalue <= 0.05) %>%
  group_by(SingleR_cluster, direction) %>%
  summarize(n = n()) %>%
  View

enrichGO_res %>%
  filter(qvalue <= 0.05) %>%
  filter(SingleR_cluster == "34 Immune") %>%
  View()

  # group_by(SingleR_cluster, direction) %>%
  # summarize(n = n()) %>%
  table()


Idents(merged_rois) <- merged_rois@meta.data$singleR_labels

RidgePlot(merged_rois, features = c("CX3CR1")) +
  scale_fill_manual(values = abca_color_pal[["cluster_colors"]]) +
  theme_minimal() +
  theme(legend.position = "none")

FeaturePlot(merged_rois, features = c("CX3CR1", "SALL1", "MOG"))


merged_rois

merged_rois %>%
  rownames() %>%
  keep(grepl("CX", .))

"CX3CR1"

"OXTR"

"IGF1"
"IGF1R"
"IL1B"

# "iba1"
"SALL1"

"cx3cr1"

"tmem119"



p <- enrichGO_res %>%
  filter(direction == "WildR Upregulated") %>%
  filter(SingleR_cluster == "28 CB GABA") %>%
  filter(qvalue <= 0.05) %>%
  slice_min(qvalue, n = 30) %>%
  slice_max(Count, n = 20, with_ties = FALSE) %>%
  ggplot(aes(x = Count, y = fct_reorder(Description, Count))) +
  geom_col(width = 0.4, fill = "red") +
  theme_bw() +
  labs(
    x = "Gene Count", y = NULL,
    title = "28 CB GABA: WildR Upregulated"
  )

  ggsave(
    glue("{wkdir}/figures/GO/GO_enrichment_28_CB_GABA_{Sys.Date()}.png"),
    p,
    width = 7, height = 5.5
  )





# key columns are center_x, center_y, slice_id, 
# subclass_label_transfer, spatial_modules_level_1

merged_rois <- readRDS(
  glue(
    "{wkdir}/data/interim/",
    "2024-04-14_merged_roi_seurat_singleR-annot.rds")
)

col_pal <- readRDS(glue("{wkdir}/data/interim/abca_color_pal.rds"))
col_pal %>% names
# merged_rois

merged_rois@meta.data %>% glimpse

cell_metadata_df <- as.data.frame(merged_rois@meta.data) %>%
  mutate(subclass_label_transfer = singleR_labels, 
  spatial_modules_level_1 = NA) %>%
  mutate(sample_id = slice_id, nCount_RNA = nCount_RNA) %>%
  mutate(class_color = purrr::map_chr(
    singleR_labels,
    ~ col_pal$class_colors[.x]
  )) %>%
  glimpse()
  

# cell_metadata_df %>%
#   filter(sample_id == "roi2_run1") %>%
#   ggplot(aes(center_x, -center_y, color = nCount_RNA)) +
#   geom_point() +
#   scale_color_viridis() +
#   theme_bw()


write.csv(
  cell_metadata_df,
  glue("{wkdir}/data/interim/registration/cell-metadata_{Sys.Date()}.csv"),
  quote = FALSE
)




#--------------------------------------------
# geneNMF analysis

library(GeneNMF)

merged_rois <- readRDS(
  glue(
    "{wkdir}/data/interim/",
    "2024-06-09_merged_roi_seurat_singleR-annot.rds"
  )
)
ndim <- 15
merged_rois <- runNMF(merged_rois, k = ndim, assay="SCT")

merged_rois <- RunUMAP(merged_rois,
  reduction = "NMF",
  dims = 1:ndim,
  reduction.name = "NMF_UMAP",
  reduction.key = "nmfUMAP_"
)

saveRDS(
  merged_rois, 
  glue(
    "{wkdir}/data/interim/",
    "{Sys.Date()}_merged_roi_seurat_singleR-annot_geneNMF.rds"
  )
)

merged_rois@meta.data %>% glimpse

DimPlot(merged_rois,
  reduction = "NMF_UMAP",
  group.by = "singleR_labels_refined",
  label = TRUE) +
  theme(
    aspect.ratio = 1, axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
    ) + 
  ggtitle("NMF UMAP") +
  NoLegend()




merged_rois@reductions %>% glimpse










# Plotting NMF DEGs

# WildR SPF DEGs - via NMF dimensions









# Read in STalign results =------------------------------------

sta_res <- read.csv(
  glue("{wkdir}/data/interim/STalign_registration/roi2_run2.csv")
) 

sta_df <- sta_res %>%
  mutate(
    X = glue("roi2_run2_{X}"),
    # convert back into pixel coordinates
    center_x = x / 0.106,
    center_y = y/ 0.106
    ) %>%
  # left_join(cell_metadata_df %>% rownames_to_column(var = "X"), by = "X") %>%
  left_join(cell_metadata_df, by = c("center_x", "center_y")) %>%
  glimpse()

sta_df %>%
  arrange(desc(nCount_RNA)) %>%
  ggplot(aes(center_x, -center_y, color = nCount_RNA)) +
  geom_point(size = 0.2) +
  scale_color_viridis() +
  theme_bw()

sta_df %>%
  ggplot(aes(x, -y, color = nCount_RNA)) +
  geom_point() +
  scale_color_viridis() +
  theme_bw()


cell_metadata_df %>%
  filter(sample_id == "roi2_run2") %>%
  dim()

sta_df %>% dim()

sta_df %>% glimpse
  drop_na(sample_id) %>%
  dim()


sta_df %>%
  ggplot(aes(x, -y, color = acronym)) +
  geom_point(size = 0.8) +
  # scale_color_viridis() +
  theme_bw() +
  coord_fixed() +
  theme(legend.position = "none")

# INTERACTIVE GPU SESS
# srun --job-name "InteractiveGPUJob" --partition=gpu --gres=gpu:1 --cpus-per-task 4 --mem 50G --time 1:00:00 --pty bash
# module load cuda/11.8.0-gcc-13.2.0-hqcjtct # Make sure this matches the version of CUDA you need
# mamba activate spatialomics








# Visualizing instability analysis results
instability_df <- read.csv(
  glue("{wkdir}/data/interim/osNMF/instability_runs/instabilityDG_PP.csv"),
  header = FALSE
) %>%
  dplyr::rename(
    "K" = "V1",
    "instability" = "V2",
    "instability_std" = "V3",
  ) %>%
  distinct() %>% 
  glimpse()

instability_df %>%
  ggplot(aes(x = K, y = instability)) +
  # geom_point() +
  geom_pointrange(aes(
    ymin=instability-instability_std,
    ymax=instability+instability_std),
  color="blue", fill="white", shape=22) +
  geom_line() +
  theme_bw()





# ______________________________________________________



 

# lintr::use_lintr(type = "tidyverse")


#' points_1 and points_2 are two sets of points in 2D space that define a boundary
#' In this case point_1 should always be the upper left-hand corner point 
#' and point_2 should be the bottom right-hand corner point. 
#' 
#' Examples: 
#' point_1 = list(boundary1 = c(x1, y1), boundary2 = c(x1, y1))
#' point_2 = list(boundary1 = c(x2, y2), boundary2 = c(x2, y2))
#' 
#' This function annotates a datafrmae with columns center_x and center_y as within a boundary
points_1[[bound_n]][1]
points_1[[bound_n]][2]

points_2[[bound_n]][1]
points_2[[bound_n]][2]

points_1 = list(
  "boundary1" = c(2000, 500), 
  "boundary2" = c(2500, 3800)
  )
points_2 = list(
  "boundary1" = c(3000, 0), 
  "boundary2" = c(3000, 3500)
  )





































# popAlign Analysis and viz

# meta_df %>% glimpse
# merged_rois <- readRDS(
#     glue("{wkdir}/data/interim/",
#     "merged_roi_seurat_filtered_staNMF-updated_2024-07-24.rds")
# )


merged_rois@meta.data$slice_id %>% unique() 

# Calculate delta bregma values between slices
merged_rois@meta.data %>% glimpse()
# ref_bregma <- -1.75

bregma_df <- merged_rois@meta.data %>% 
  dplyr::select(slice_id, bregma) %>%
  dplyr::distinct() %>%
  mutate(bregma_delta = bregma - ref_bregma) %>%
  glimpse()


# Adding osNMF Quantile metadata to seurat object

popAlign_output_root <- glue(
  "{wkdir}/data/interim/popAlign_results/2024-08-07_SCT-counts"
)

delta_stats_df <- list.files(
  popAlign_output_root,
  recursive = TRUE,
  pattern = "delta_stats.csv"
) %>%
  purrr::set_names() %>%
  purrr::map(
    ~ read.csv(glue("{popAlign_output_root}/{.}"))
  ) %>%
    bind_rows(.id = "path") %>%
      mutate(nmf = strex::str_before_first(path, "/")) %>%
      dplyr::select(-c(origidx, path)) %>%
      filter(nmf != "testruns", 
      # orderedsamples != "SPF"
      ) %>% 
      glimpse()

mu_mat <- delta_stats_df %>%
  filter(orderedsamples != "SPF run2 roi2") %>%
  dplyr::group_by(nmf, cell_type) %>%
  dplyr::summarize(mean_delta_mu = mean(mean_delta_mu)) %>%
  pivot_wider(
    names_from = "cell_type",
    values_from = "mean_delta_mu", values_fill = 0
  ) %>%
  tibble::column_to_rownames(var = "nmf") %>%
    as.matrix() %>%
    glimpse()

nmf_order <- seriate_matrix_rows(mu_mat)
cell_order <- seriate_matrix_rows(t(mu_mat))
# print(names(mu_mat)) == trimws(names(mu_mat))

delta_stats_df_long <- delta_stats_df %>%
  pivot_longer(
    cols = c(mean_delta_w, mean_delta_mu, mean_delta_cov),
    names_to = "mean_type",
    values_to = "mean_value"
  ) %>%
  pivot_longer(
    cols = c(pvals_w, pvals_mu, pvals_cov),
    names_to = "pval_type",
    values_to = "pval_value"
  ) %>%
    mutate(mean_type = case_when(
      mean_type == "mean_delta_w" ~ "Delta * omega",
      mean_type == "mean_delta_mu" ~ "Delta * mu",
      mean_type == "mean_delta_cov" ~ "Delta * Sigma"
    )) %>%
      mutate(
        nmf = factor(nmf, levels = nmf_order),
        cell_type = factor(cell_type, levels = cell_order)
      ) %>%
    mutate(group = strex::str_before_first(orderedsamples, " ")) %>%
    left_join(bregma_df, by = c("orderedsamples" = "slice_id")) %>%
  glimpse()

metrics <- c("mean_delta_mu", "mean_delta_w", "mean_delta_cov")
cell_types <- unique(delta_stats_df$cell_type)
nmf_dims <- unique(delta_stats_df$nmf)

# > cell = "20 MB GABA"
# > nmf = "osNMF_9"

stats_df <- tibble()
for (nmf_dim in nmf_dims) {
  for (metric in metrics) {
    for (cell in cell_types) {
      
      pop_stat_fig_dir <- glue("{wkdir}/figures/popAlign/delta_stats/{Sys.Date()}")
      dir.create(pop_stat_fig_dir, showWarnings = FALSE, recursive = TRUE)
      
      df <- delta_stats_df %>%
        filter(cell_type == cell & nmf == nmf_dim) %>%
        mutate(group = strex::str_before_first(orderedsamples, " ")) %>%
        left_join(bregma_df, by = c("orderedsamples" = "slice_id")) %>%
        filter(orderedsamples != "SPF run2 roi2") %>%
        mutate(run = strex::str_after_first(orderedsamples, " ") %>% 
          strex::str_before_first(., " "))
      
      if (nrow(df) == 0) next
      message(glue("Plotting: {metric} - {cell} - {nmf_dim} "))

      # p <- df %>%
      #   ggplot(aes(x = group, y = !!sym(metric))) +
      #   geom_boxplot() +
      #   geom_point(aes(color = bregma_delta), size = 2, position = position_jitter(width=0.2)) +
      #   scale_color_viridis_c(option = "magma") +
      #   labs(x = NULL, y = NULL, ) +
      #   theme_bw()
      
      # dftst <- df %>% rename(test_metric = metric)
      fm <- formula(glue("{metric} ~ group + run + bregma_delta"))
      lm_stats <- lm(fm, data = df) %>% broom::tidy() %>%
        mutate(cell_type = cell, nmf_dim = nmf_dim, metric_name = metric)
      stats_df %<>% bind_rows(lm_stats)

      # ggsave(
      #   glue("{pop_stat_fig_dir}/",
      #   "Ref_SPF-run4-roi4_{nmf_dim}_{cell}.png"),
      #   p,
      #   width = 3, height = 3
      # )
    }
  }
}


stats_df %>% glimpse()
stats_df %>% View

stats_df %>%
  filter(term == "groupWILDR") %>% View


# stats_df %>%
#   filter(metric_name == "mean_delta_mu") %>%
#   filter(term == "groupWILDR") %>% View


p_pop_lm_stats <- stats_df %>%
  # filter(term == "groupWILDR") %>%
  filter(term != "(Intercept)") %>%
  filter(p_value < 0.05) %>%
  ggplot(aes(x = estimate, y = nmf_dim )) +
  geom_col(aes(fill = cell_type), position = "identity", width = 0.5) +
  facet_grid(term ~ metric_name, scales = "free", space = "free_y") +
  scale_fill_manual(values = abca_color_pal[["cluster_colors"]]) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw()

ggsave(
  glue("{wkdir}/figures/popAlign/delts_lm_stats_{Sys.Date()}.png"),
  p_pop_lm_stats,
  width = 11, height = 10
)

# delta_stats_df_long %>%
#   dplyr::select(-c(pval_type, pval_value)) %>%
#   group_by(nmf, cell_type, mean_type) %>%
#   # nest() %>%
#   nest(data = c(orderedsamples, group, mean_type, mean_value)) %>%
#   filter(cell_type == "01 IT-ET Glut" & nmf == "osNMF_10" & mean_type == "Delta * mu") %>%
#   unnest(cols = c(data))%>% print(n = 100)
#   # nest(data = c(orderedsamples, nmf)) %>% print(n = 100) %>%
#   # nest(data = c(cell_type, nmf)) %>% print(n = 100) %>%
#   # filter(orderedsamples == "SPF run1 roi1" & mean_type == "Delta * mu") %>%
#   # unnest(data)
#   # nest(data = c(cell_type, nmf, mean_type, mean_value)) %>% print(n = 100)
# filter()


# delta_stats_df_long %>%
#   dplyr::select(-c(pval_type, pval_value)) %>%
#   group_by(nmf, cell_type, mean_type) %>%
#   nest(data = c(orderedsamples, group, mean_type, mean_value)) %>%
#   filter(cell_type == "01 IT-ET Glut" & nmf == "osNMF_10") %>%
#   unnest



delta_stats_df_long_summary <- delta_stats_df_long %>%
  filter(orderedsamples != "SPF run2 roi2") %>%
  group_by(group, cell_type, nmf, mean_type) %>%
  summarize(
    mean = mean(mean_value, na.rm = TRUE),
    sd = sd(mean_value, na.rm = TRUE),
    n = n()
  ) %>%
  glimpse()


p_pop_csum <- delta_stats_df_long_summary %>%
  # mutate(orderedsamples = gsub(" ", "_", orderedsamples)) %>%
  ggplot(aes(x = mean, y = nmf)) +
  geom_bar(stat = "identity", width = 0.7,
    aes(fill = cell_type)) +
  facet_grid(orderedsamples ~ mean_type, scales = "free", space = "free",
  labeller = label_parsed) +
  scale_fill_manual(values = abca_color_pal[["cluster_colors"]]) +
  labs(x = NULL, y = NULL) +
  theme_bw()

ggsave(
  glue("{wkdir}/figures/popAlign/GMM_NMF_csum_delta_stats_{Sys.Date()}.png"),
  p_pop_csum,
  width = 11, height = 10
)




all_delta_stats_df_long %>% 
  filter(ref_slice == "SPF run1 roi1" & ref_slice !=  test_slice) %>% View
  glimpse




# exploring data












# # for multiple NMF dims as y axis
# allvall_stats %>%
#   janitor::clean_names() %>%
#     mutate(metric = case_when(
#     metric == "mean_delta_w" ~ "Delta * omega",
#     metric == "mean_delta_mu" ~ "Delta * mu",
#     metric == "mean_delta_cov" ~ "Delta * Sigma"
#   )) %>%
#   filter(
#     # term == "test_groupwildr",
#      p_value < 0.05
#      ) %>%
#   ggplot(aes(x = value, y = nmf_dim)) +
#   geom_col(aes(fill = cell_type), position = "identity", width = 0.5) +
#   facet_grid(term ~ metric, labeller = label_parsed, 
#     scales = "free", space = "free_y") +
#     scale_fill_manual(values = abca_color_pal[["class_colors"]]) +
#   theme_bw()


# filtered_df %>%
#   mutate(sample_pair = 
#   purrr::map2_chr(ref_slice, test_slice, 
#     ~paste(sort(c(.x, .y)), collapse = "_"))
#   ) %>%
#   # dplyr::select(ref_slice, test_slice, sample_pair) %>%
#   View














# All against all popAlign delta stats
library(strex)
library(nlme)


all_delta_stats_paths <- list.files(
  glue("{wkdir}/data/interim/popAlign_results/2024-08-08_SCT-counts-nmf"),
  recursive = TRUE,
  pattern = "delta_stats.csv",
  full.names = TRUE
) %>%
  keep(~ grepl("multicomp", .))

all_delta_stats_df_raw <- all_delta_stats_paths %>%
  purrr::set_names() %>%
  purrr::map(
    ~ read.csv(.x)
  ) %>%
  bind_rows(.id = "filepath") %>%
  glimpse()

meta_df <- merged_rois@meta.data %>% glimpse()
ref_meta <- meta_df %>%
  tibble() %>%
  select(
    ref_slice = slice_id, ref_group = group, 
    ref_run = run, ref_roi = roi, ref_bregma = bregma
    ) %>%
  distinct()

test_meta <- meta_df %>%
  tibble() %>%
  select(
    test_slice = slice_id, test_group = group, 
    test_run = run, test_roi = roi, test_bregma = bregma
    ) %>%
  distinct()

all_delta_stats_df <- all_delta_stats_df_raw %>%
  mutate(
    ref_slice = str_before_last(filepath, "/") %>% str_after_last("/"), 
    nmf_dim = str_before_last(filepath, "/multicomp") %>% str_after_last("/")
    ) %>%
  tibble() %>%
  dplyr::rename(test_slice = orderedsamples) %>%
  left_join(ref_meta, by = c("ref_slice" )) %>%
  left_join(test_meta, by = c("test_slice" )) %>%
  dplyr::mutate(bregma_delta = ref_bregma - test_bregma) %>%
  mutate(pair_type = case_when(
    ref_group == test_group & ref_group == "spf" ~ "Within-SPF",
    ref_group == test_group & ref_group == "wildr" ~ "Within-WildR",
    ref_group == "wildr" ~ "Between-WildR-to-SPF",
    ref_group == "spf" ~ "Between-SPF-to-WildRSPF",
    TRUE ~ "Error"
  )) %>%
  mutate(sample_pair = purrr::map2_chr(ref_slice, test_slice, 
    ~paste(sort(c(.x, .y)), collapse = "__"))
  ) %>%
  mutate(sample_pair_asym = glue("{ref_slice}__{test_slice}")) %>%
  glimpse()



all_delta_stats_df_long <- all_delta_stats_df %>% 
  dplyr::select(-contains("pvals"), -origidx) %>%
  pivot_longer(
    cols = c(mean_delta_w, mean_delta_mu, mean_delta_cov),
    names_to = "metric",
    values_to = "metric_value"
  ) %>%
    mutate(metric_formatted = case_when(
      metric == "mean_delta_w" ~ "Delta * omega",
      metric == "mean_delta_mu" ~ "Delta * mu",
      metric == "mean_delta_cov" ~ "Delta * Sigma"
  )) %>%
  glimpse()


sp_order <- all_delta_stats_df_long %>% 
  filter(test_slice != ref_slice) %>%
  filter(metric == "mean_delta_cov") %>%
  group_by(sample_pair_asym) %>%
  summarize(mean_metric = mean(metric_value)) %>%
  arrange(mean_metric) %>%
  pull(sample_pair_asym)

p_dist <- all_delta_stats_df_long %>%
  filter(test_slice != ref_slice) %>%
  mutate(sample_pair_asym = factor(sample_pair_asym, levels = sp_order, ordered  = TRUE)) %>%
  ggplot(aes(y=sample_pair_asym, x=metric_value)) +
  geom_point(aes(color = cell_type),
    position = position_jitter(width = 0.05), alpha = 0.4, size = 0.8) +
  geom_boxplot(alpha = 0.3, outlier.alpha = 0) +
  facet_grid(pair_type ~ metric_formatted, labeller = label_parsed, 
    scales = "free", space = "free_y") +
  scale_color_manual(values = abca_color_pal[["class_colors"]]) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  theme_bw()

ggsave(
  glue("{wkdir}/figures/popAlign/delta_stats_pairwise_grouped_NMF-split_{Sys.Date()}.png"),
  p_dist,
  width = 17, height = 25
)



allvall_stats <- tibble()
model_list <- list()


# Loop through parameters, dimensions, and cell types
for (param in c("mean_delta_w", "mean_delta_mu", "mean_delta_cov")) {
  for (nd in unique(all_delta_stats_df$nmf_dim)) {
    nmf_dim_df <- all_delta_stats_df %>%
      filter(nmf_dim == nd) %>%
      filter(ref_group != test_group)
    for (ct in unique(nmf_dim_df$cell_type)) {
      message("Fitting model ", param, " ", nd, " ", ct)
      filtered_df <- nmf_dim_df %>% 
        filter(cell_type == ct)
      
      tryCatch({
        # Fit the nlme model
        model_list[[nd]][[ct]][[param]] <- nlme::lme(fixed = as.formula(glue("{param} ~ test_group + bregma_delta")),
                                                     random = ~ 1 | ref_run/ref_slice,
                                                     #  random = ~ 1 | ref_run + sample_pair,
                                                     data = filtered_df,
                                                     method = "REML")
        
        # Extract the model summary and convert to a tidy format
        res <- summary(model_list[[nd]][[ct]][[param]])$tTable %>%
          as.data.frame() %>%
          rownames_to_column("term") %>%
          as_tibble() %>%
          mutate(nmf_dim = nd, cell_type = ct, metric = param)
        
        # Combine results
        allvall_stats %<>% bind_rows(res)
      }, error = function(e) {
        message("Error fitting model for param: ", param, ", dim: ", nd, ", cell type: ", ct)
        message("Error: ", e$message)
      })
    }
  }
}


# # Loop through parameters, dimensions, and cell types
# for (param in c("mean_delta_w", "mean_delta_mu", "mean_delta_cov")) {
#   for (nd in unique(all_delta_stats_df$nmf_dim)) {
#     nmf_dim_df <- all_delta_stats_df %>%
#       filter(nmf_dim == nd) %>%
#       filter(ref_group != test_group)
#     for (ct in unique(nmf_dim_df$cell_type)) {
#       message("Fitting model ", param, " ", nd, " ", ct)
#       filtered_df <- nmf_dim_df %>% 
#         filter(cell_type == ct)
      
#       # Fit the nlme model
#       model_list[[nd]][[ct]][[param]] <- nlme::lme(fixed = as.formula(glue("{param} ~ test_group + bregma_delta")),
#                                              random = ~ 1 | ref_run/ref_slice,
#                                             #  random = ~ 1 | ref_run + sample_pair,
#                                              data = filtered_df,
#                                              method = "REML")
      
#       # Extract the model summary and convert to a tidy format
#       res <- summary(model_list[[nd]][[ct]][[param]])$tTable %>%
#         as.data.frame() %>%
#         rownames_to_column("term") %>%
#         as_tibble() %>%
#         mutate(nmf_dim = nd, cell_type = ct, metric = param)
      
#       # Combine results
#       allvall_stats %<>% bind_rows(res)
#     }
#   }
# }

allvall_stats %<>%
  janitor::clean_names() %>%
    mutate(metric = case_when(
    metric == "mean_delta_w" ~ "beta (Delta * omega)",
    metric == "mean_delta_mu" ~ "beta (Delta * mu)",
    metric == "mean_delta_cov" ~ "beta (Delta * Sigma)"
  ))

# saveRDS(
#   allvall_stats,
#   glue(
#     "{wkdir}/data/interim/popAlign_results/2024-08-01_SCT-counts/",
#     "whole_slice/lmm_stats/",
#     "whole-slice_lmm_stats_{Sys.Date()}.rds"
#     )
# )


# allvall_stats <- readRDS(
#   glue(
#     "{wkdir}/data/interim/popAlign_results/2024-08-01_SCT-counts/",
#     "whole_slice/lmm_stats/",
#     "whole-slice_lmm_stats_2024-08-07.rds"
#     )
# )
# allvall_stats %>% glimpse()

# print significant results in a table
allvall_stats %>% 
  filter(grepl("test_group", term)) %>% 
  filter(p_value < 0.05) %>%
  DT::datatable(
    options = list(
      pageLength = 10,
      scrollX = TRUE
    )
  )

allvall_stats$term %>% unique()

# plotting significant results
p_lmm_sig_bars <- allvall_stats %>%
  janitor::clean_names() %>%
  filter(
    term == "test_groupwildr",
    # term != "(Intercept)",
     p_value <= 0.05
     ) %>%
  ggplot(aes(x = value, y = cell_type)) +
  geom_col(aes(fill = cell_type), position = "identity", width = 0.5) +
  facet_grid(nmf_dim ~ metric, labeller = label_parsed, 
    scales = "free", space = "free_y") +
    scale_fill_manual(values = abca_color_pal[["class_colors"]]) +
    labs(x = NULL, y = NULL) +
  theme_bw()

ggsave(
  glue("{wkdir}/figures/popAlign/delta_stats_lm_sig_bars_NMF-split_{Sys.Date()}.png"),
  p_lmm_sig_bars,
  width = 8, height = 10
)









#______________________________________________________________________________

# Making PCA plots of cell types 
#______________________________________________________________________________


seur_objs <- SplitObject(merged_rois, split.by = "singleR_labels")
celltypes <- unique(merged_rois@meta.data$singleR_labels)
seur_objs <- seur_objs %>% purrr::map(~ RunPCA(.))

saveRDS(
  seur_objs,
    glue(
      "{wkdir}/data/interim/",
      "merged_roi_seurat_filtered_staNMF-updated",
      "_celltype-list_{Sys.Date()}.rds")
)


seur_objs <- readRDS(
    glue(
      "{wkdir}/data/interim/",
      "merged_roi_seurat_filtered_staNMF-updated",
      "_celltype-list_2024-08-12.rds"
      )
)


# Functions for plotting
dimplot_manual <- function(embed_data, group_var, 
                          dim_x, dim_y, contour = FALSE){
  p_contour <- embed_data %>%
    ggplot(aes(x = !!sym(dim_x), y = !!sym(dim_y), color = !!sym(group_var))) +
    theme_bw()
  if (contour) {
    p_contour  <- p_contour + geom_density_2d(bins = 15)
  } else{
    p_contour <- p_contour + geom_point(size = 0.3)
  }
  return(p_contour)
}

seur_pca_elbowplot <- function(seur_obj){
  # Creating a PCA Elbow plot
  pca_sd <- seur_obj[["pca"]]@stdev
  pct <- (pca_sd / sum(pca_sd)) * 100
  # cumulative percents for each PC
  cumu <- cumsum(pct)

  # # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 80 & pct < 5)[1]
  # last point where change of % of variation is more than 0.1%.
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

  # Create a data frame for plotting
  plot_df <- data.frame(pct = pct, 
            cumu = cumu, 
            rank = 1:length(pct))

  # Elbow plot to visualize 
  ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > co2)) + 
    geom_text() + 
    geom_vline(xintercept = 90, color = "grey") + 
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
    scale_color_manual(values = c("red", "black")) +
    labs(x = "Cumulative Percent of Variation",
        y = "Percent of Variation") +
    guides(color = FALSE) +
    theme_bw()
}


for (celltype in celltypes) {
  message("\nPlotting: ", celltype)
  # subset_obj <- merged_rois %>% subset(singleR_labels == celltype)
  # subset_obj <- seur_objs[[celltype]]
  # seur_objs[[celltype]] <- RunPCA(seur_objs[[celltype]])
  # seur_objs[[celltype]] <- RunUMAP(seur_objs[[celltype]], dims = 1:10)  # Adjust dims as needed

  # # merging embedding data and metadatae
  # print(seur_objs[[celltype]][["pca"]]@stdev)
  embed_data <- seur_objs[[celltype]]@meta.data %>%
    bind_cols(Embeddings(seur_objs[[celltype]], reduction = "pca"))

  # p_contour_group_pca <- 
  #   dimplot_manual(embed_data, "group", "PC_1", "PC_2", contour = TRUE)

  p_elbow <- seur_pca_elbowplot(seur_objs[[celltype]]) %>% print()
  p_pt_group_pca_1 <- 
    dimplot_manual(embed_data, "group", "PC_1", "PC_2", contour = FALSE)
  p_pt_group_pca_2 <- 
    dimplot_manual(embed_data, "group", "PC_3", "PC_4", contour = FALSE)
  p_pt_group_pca_3 <- 
    dimplot_manual(embed_data, "group", "PC_5", "PC_6", contour = FALSE)

  cmbd_group <- 
    (p_elbow + p_pt_group_pca_1) / (p_pt_group_pca_2 + p_pt_group_pca_3) + 
    plot_layout(guides = "collect") +
    plot_annotation(title = glue("{celltype}"))

  ggsave(
    glue("{wkdir}/figures/SingleR/DimRed/PCA_{celltype}_{Sys.Date()}.png"),
    cmbd_group,
    width = 12, height = 12
  )
}


#______________________________________________________________________________






#______________________________________________________________________________








#______________________________________________________________________________

# SingleR cell type gene exp. Markers
#______________________________________________________________________________


# temporarily setting cluster labels to singleR labels
Idents(merged_rois) <- merged_rois@meta.data$singleR_labels

# Explore differentially expressed genes between SingleR Cell labels

# options(future.globals.maxSize = 10 * 1024^3)  # Increase to 10 GiB
# future::plan("multisession", workers = 32)
singler_celltype_markers <- 
  unique(merged_rois@meta.data$singleR_labels) %>%
  purrr::set_names() %>%
  # furrr::future_map(
  purrr::map(
    ~ Seurat::FindMarkers(
      merged_rois,
      ident.1 = .,
      ident.2 = NULL
    ) %>%
      as.data.frame() %>%
      rownames_to_column(var = "gene_symbol") %>%
      mutate(SingleR_cluster = .x) #%>%
      # filter(p_val_adj <= 0.1)
  ) %>%
  bind_rows() %>%
  glimpse()


saveRDS(
  singler_celltype_markers, 
  glue("{wkdir}/data/interim/")
)

singler_celltype_markers %>% glimpse
singler_celltype_markers %>%
  select(-p_val) %>%
  filter(SingleR_cluster == "12 HY GABA") %>%
  View

singler_celltype_markers$SingleR_cluster %>% unique

singler_celltype_markers %>%
  select(-p_val) %>%
  filter(SingleR_cluster == "18 TH Glut") %>%
  View





#______________________________________________________________________________

# For genes of interest, produce gene exp. summary and distrib plots
#______________________________________________________________________________


cois <- c(
  "OTP",
  "TMEM132D",
  "GPX3",
  "C1QA",
  "ALDH1L1",
  "STAT3",
  "P2RY12",
  "CALCR",
  "GABRQ",
  "GPR101",
  "SELPLG",
  "STX1A",
  "SLC1A3",
  "RASAL1"
)

library(data.table)

count_df <- merged_rois@assays$SCT@counts

# Print summary for each gene in cois
for (gene in cois) {
  if (gene %in% rownames(count_df)) {
    cat("\nSummary for", gene, ":\n")
    print(summary(count_df[gene,]))
  } else {
    cat("\nGene", gene, "not found in the dataset.\n")
  }
}


#______________________________________________________________________________






meta_df_gated <- readRDS(
  glue(
    "{wkdir}/data/interim/",
    "gated-seurat-metadata_2024-09-17.rds"
    )
  )

# create a new column with cell types and gated regions
meta_df_gated %<>%
  mutate(final_gate = case_when(
    final_gate != "" ~ glue(" {final_gate}"),
    TRUE ~ final_gate
  )) %>% 
  mutate(gated_cell_labels = glue("{singleR_labels}{final_gate}")) %>% 
  glimpse()

meta_df_gated$gated_cell_labels %>% table()

# meta_df_gated %>% glimpse()
# meta_df %>% glimpse()
merged_rois$final_gate <- meta_df_gated$final_gate
merged_rois$gated_cell_labels <- meta_df_gated$gated_cell_labels
merged_rois@meta.data %>% glimpse()


amyg_regions <- c("AMYG_CeMe", "AMYG_CT", "AMYG_BLA", "AMYG")
relab_df <- data.frame()

for (region in amyg_regions) {
  temp_df <- merged_rois@meta.data %>% 
    filter(grepl(region, gated_cell_labels)) %>%
    pull(singleR_labels) %>%
    table() %>%
    as.data.frame() %>% 
    mutate(region = region, 
    perc = Freq / sum(Freq)
    )
  
  relab_df %<>% bind_rows(temp_df)

}



relab_df %>%
  ggplot(aes(x = region, y = perc, fill = .)) +
  scale_fill_manual(values = abca_color_pal[["cluster_colors"]]) +
  geom_col()




seurat_obj <- merged_rois
gene <- "C1QA"

# Function to create a boxplot for a single gene expression across cell types
plot_gene_expression <- function(seurat_obj, gene, assay = "SCT") {
  # Ensure the gene exists in the dataset
  if (!gene %in% rownames(seurat_obj)) {
    stop(paste("Gene", gene, "not found in the dataset"))
  }

  # Extract expression data and metadata
  expression_data <- FetchData(seurat_obj, 
    vars = c(gene, "singleR_labels"), 
    layer = "counts", 
    assay = assay
  )
  
  # Create the plot
  ggplot(expression_data, aes(x = singleR_labels, y = .data[[gene]], fill = singleR_labels)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
    scale_fill_manual(values = abca_color_pal[["cluster_colors"]]) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = paste("Expression of", gene, "across cell types"),
         x = "Cell Type",
         y = paste(gene, "Expression (SCT)"))
}

# Example usage (replace 'YourGeneOfInterest' with an actual gene name)
p <- plot_gene_expression(merged_rois, "C1QA")
ggsave(
  glue("{wkdir}/figures/SingleR/GeneExp/C1QA_expression_boxplot_{Sys.Date()}.png"),
  p,
  width = 12, height = 8
)








#-----------------------------------------------------------------------------
# Testing for Piriform specific markers

meta_df_gated <- readRDS(
  glue(
    "{wkdir}/data/interim/",
    "gated-seurat-metadata_2024-09-17.rds"
    )
)

# create a new column with cell types and gated regions
meta_df_gated %<>%
  mutate(final_gate = case_when(
    final_gate != "" ~ glue(" {final_gate}"),
    TRUE ~ final_gate
  )) %>% 
  mutate(amyg_gate = case_when(
    grepl("AMYG", final_gate) ~ "Amygdala",
    grepl("Piriform", final_gate) ~ "Piriform",
    TRUE ~ ""
  )) %>% 
  mutate(gated_cell_labels = glue("{singleR_labels}{final_gate}")) %>% 
  mutate(gated_amyg_cell_labels = glue("{singleR_labels} {amyg_gate}")) %>% 
  glimpse()

# Add gated_amyg_cell_labels to the Seurat object
merged_rois$gated_amyg_cell_labels <- meta_df_gated$gated_amyg_cell_labels


# Display tables and glimpse for verification
meta_df_gated$amyg_gate %>% table()
meta_df_gated$gated_cell_labels %>% table()


# Identify neuronal cell types with > 500 cells total in all piriform samples
neuronal_cell_types <- meta_df_gated$gated_amyg_cell_labels %>% 
  table() %>% 
  as.data.frame() %>% 
  dplyr::rename("cell_type" = ".", "count" = "Freq") %>% 
  filter(grepl("GABA|Glut", cell_type)) %>% 
  filter(grepl("Piriform", cell_type)) %>% 
  filter(count > 500) %>% 
  pull(cell_type)  %>% 
  as.character() %>% 
  gsub(" Piriform", "", .)

neuronal_cell_types


# Define bad slice_ids to exclude
merged_rois_filt <- subset(merged_rois,
  slice_id %nin% c(
    "WILDR run1 roi3", "SPF run2 roi3",
    "SPF run3 roi1", "WILDR run3 roi2"
  )
)

merged_rois_filt@meta.data %>% glimpse()

# Set the Idents of the Seurat object to gated_amyg_cell_labels
Idents(merged_rois_filt) <- "gated_amyg_cell_labels"

# Verify the new identity classes
table(Idents(merged_rois_filt))

# Perform differential expression analysis for each valid neuronal cell type
amyg_vs_piriform_markers <- 
  neuronal_cell_types %>%
  purrr::set_names() %>%
  purrr::map(
    ~ Seurat::FindMarkers(
      merged_rois_filt,
      ident.1 = paste0(.x, " Amygdala"),
      ident.2 = paste0(.x, " Piriform"),
      verbose = FALSE
    ) %>%
      as.data.frame() %>%
      rownames_to_column(var = "gene_symbol") %>%
      mutate(cell_type = .x,
             comparison = "Amygdala vs Piriform")
  ) %>%
  bind_rows() %>%
  glimpse()

# Display the results
print(amyg_vs_piriform_markers)

amyg_vs_piriform_markers %>% glimpse()

amyg_vs_piriform_markers %>% 
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(data = amyg_vs_piriform_markers %>% filter(p_val_adj > 0.05 | abs(avg_log2FC) < ),
    color = "grey", alpha = 0.8) +
  geom_point(data = amyg_vs_piriform_markers %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.2),
    color = "darkred", alpha = 0.8) +
  facet_wrap(~ cell_type) +
  theme_bw()




#-----------------------------------------------------------------------------
# Adding new gates
#-----------------------------------------------------------------------------

annot_paths <- list.files(
  glue("{wkdir}/data/input/JAG_annotations_v2"),
  full.names = TRUE
)

annot_paths_r2 <- list.files(
  glue("{wkdir}/data/input/JAG_annotations_2024-10-15"),
  full.names = TRUE
)

annot_paths <- c(annot_paths, annot_paths_r2)

# Quality Check for file names
annot_types <- annot_paths %>% basename(.) %>% 
      str_before_last("_") %>% 
      str_after_first("_")

annot_types %>% table()

clean_gate_tries <- function(gate_list){
  unique_numbers <- gate_list %>%
    na.omit() %>% 
    unique() %>% 
    strsplit(",") %>%
    unlist() %>%
    as.numeric() %>%
    unique() %>%
    sort()
  }

gated_res_list <- annot_paths %>% 
  set_names(basename(.)) %>% 
  map(~ {
    annotation_type <- basename(.) %>% 
      str_before_last("_") %>% 
      str_after_first("_")
    
    df <- readRDS(.)
    clean_gates <- df$gated %>% clean_gate_tries()
    max_gate <- max(clean_gates, na.rm = TRUE)
    
    df %>%
      mutate(final_gate = grepl(max_gate, gated)) %>% 
      filter(final_gate) %>%
      mutate(!!sym(annotation_type) := TRUE)
  })

# Prepare the gated_res_list dataframes
gated_res_list_trim <- gated_res_list %>%
  map(~ .x %>% 
    select(starts_with("gated-")) %>% 
    mutate(cell_id = rownames(.))
  )

# compile a list of all cell ids for each gating type
agg_annots <- list()
for (ii in unique(annot_types)) {
  annot_type_paths <- names(gated_res_list_trim) %>% 
    keep(grepl(ii, .))
  agg_annots[[ii]] <- 
    gated_res_list_trim[annot_type_paths] %>% bind_rows()
}

# Join the gated results to meta_df
meta_df_updated <- meta_df %>%
  mutate(cell_id = rownames(.))

# now joining each aggregated list 
for (df in agg_annots) {
  meta_df_updated %<>% full_join(df)
}
meta_df_updated %>% glimpse()
meta_df_updated$slice_id %>% unique

meta_df_updated <- meta_df_updated %>% 
  mutate(gated_count = rowSums(select(., starts_with("gated-")), na.rm = TRUE)) %>%
  glimpse()

# counting the frequency of cell type labeling 
meta_df_updated$gated_count %>% table()

p_gate_count <- meta_df_updated %>% 
  arrange(slice_id) %>% 
  mutate(gated_count = as.character(gated_count)) %>% 
  ggplot(aes(center_x, center_y)) +
  geom_point(aes(color = gated_count ), size = 0.3, alpha = 0.5) +
  theme_bw() +
  coord_fixed() +
  scale_y_reverse() +
  scale_color_manual(values = c("0" = "lightgray", "1" = "red", "2" = "green")) + # set custom colors
  facet_wrap(~slice_id, ncol = 4)

ggsave(
  filename = glue("{wkdir}/figures/gating/gated_count_scatter_{Sys.Date()}.png"),
  plot = p_gate_count,
  width = 20,
  height = 20
)

meta_df_updated %>% glimpse()

# Strategy for prioritizing cells grouped in multiple gates
meta_df_gated <- meta_df_updated %>%
  mutate(across(starts_with("gated-"), 
    ~ ifelse(. == TRUE, cur_column(), NA))) %>%
  unite("temp_gate", all_of(unique(annot_types)), remove = FALSE, sep = ",", na.rm = TRUE) %>% 
  mutate(final_gate = case_when(
    is.na(temp_gate) | temp_gate == "" ~ "",
    grepl("gated-RT", temp_gate) ~ "RT",
    grepl("gated-Piriform", temp_gate) ~ "Piriform",
    grepl("gated-AMYG-BLA", temp_gate) ~ "AMYG_BLA",
    grepl("gated-AMYG-Cor_Tran", temp_gate) ~ "AMYG_CT",
    grepl("gated-AMYG-CeMe", temp_gate) ~ "AMYG_CeMe",
    grepl("gated-caudate", temp_gate) ~ "CAUD",
    grepl("gated-thalamus", temp_gate) ~ "THAL",
    grepl("gated-hippocampus", temp_gate) ~ "HIPP",
    grepl("gated-sscortex", temp_gate) ~ "SSC",
    grepl("gated-hyp", temp_gate) ~ "HYP",
    TRUE ~ "ERROR"
    )
  )

meta_df_gated$temp_gate %>% table()
meta_df_gated$final_gate %>% table()

# Visualizing manually identified regions
p_gate_ids <- meta_df_gated %>% 
  arrange(slice_id) %>% 
  ggplot(aes(center_x, center_y)) +
  geom_point(data = filter(meta_df_gated, final_gate == ""),
    size = 0.3, alpha = 0.5, color = "grey") +
  geom_point(data = filter(meta_df_gated, final_gate != ""),
    aes(color = final_gate ), size = 0.3, alpha = 0.5) +
  theme_bw() +
  coord_fixed() +
  scale_y_reverse() +
  scale_color_d3() +
  facet_wrap(~slice_id, ncol = 4) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(
  filename = glue("{wkdir}/figures/gating/gated_labels_scatter_{Sys.Date()}.png"),
  plot = p_gate_ids,
  width = 20,
  height = 20
)

meta_df_gated %>% glimpse

# Saving results
saveRDS(
  meta_df_gated,
  glue(
    "{wkdir}/data/interim/",
    "gated-seurat-metadata_{Sys.Date()}.rds"
    )
  )





