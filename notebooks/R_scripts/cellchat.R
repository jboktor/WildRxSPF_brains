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

library(CellChat)

# setting paths
homedir <- "/central/groups/mthomson/jboktor"
wkdir <- glue("{homedir}/spatial_genomics/jess_2024-01-23")
source(glue("{wkdir}/notebooks/R_scripts/_misc_functions.R"))

abca_color_pal <- readRDS(glue("{wkdir}/data/interim/abca_color_pal.rds"))
merged_rois <- readRDS(
    glue(
      "{wkdir}/data/interim/",
      "merged_roi_seurat_filtered_staNMF-updated_2024-07-24.rds"
    )
)
meta_df <- merged_rois@meta.data %>% glimpse()



# cellchat <- createCellChat(object = merged_rois, 
#   group.by = "singleR_labels_refined", 
#   assay = "RNA"
# )







# Prepare input data for CelChat analysis
data.input1 = Seurat::GetAssayData(seu1, slot = "data", assay = "SCT") # normalized data matrix
data.input2 = Seurat::GetAssayData(seu2, slot = "data", assay = "SCT") 

genes.common <- intersect(rownames(data.input1), rownames(data.input2))
colnames(data.input1) <- paste0("A1_", colnames(data.input1))
colnames(data.input2) <- paste0("A2_", colnames(data.input2))
data.input <- cbind(data.input1[genes.common, ], data.input2[genes.common, ])

# define the meta data
# a column named `samples` should be provided for spatial transcriptomics analysis, which is useful for analyzing cell-cell communication by aggregating multiple samples/replicates. Of note, for comparison analysis across different conditions, users still need to create a CellChat object seperately for each condition.  
meta1 = data.frame(labels = Idents(seu1), samples = "A1") # manually create a dataframe consisting of the cell labels
meta2 = data.frame(labels = Idents(seu2), samples = "A2") 

meta <- rbind(meta1, meta2)
rownames(meta) <- colnames(data.input)
# a factor level should be defined for the `meta$labels` and `meta$samples`
meta$labels <- factor(meta$labels, levels = levels(Idents(seu1)))
meta$samples <- factor(meta$samples, levels = c("A1", "A2"))
unique(meta$labels) # check the cell labels
#>  [1] Myofibroblasts        Cycling Cells         Undifferentiated     
#>  [4] Glial                 DC2                   B-Cells              
#>  [7] Goblets               T-Cells               Plasma Cells         
#> [10] Colonocytes           Crypt Top Colonocytes ILCs                 
#> [13] Mast Cells            BEST4+/OTOP2+ Cell    NK                   
#> [16] Endothelial 2         Enteroendocrines      Stromal 3            
#> [19] Stromal 1             Stromal 2             Pericytes            
#> [22] Endothelial 1         DC1                   Macrophages&Monocytes
#> [25] Stromal 4            
#> 25 Levels: B-Cells BEST4+/OTOP2+ Cell Colonocytes ... Undifferentiated
unique(meta$samples) # check the sample labels
#> [1] A1 A2
#> Levels: A1 A2

# load spatial transcriptomics information
# Spatial locations of spots from full (NOT high/low) resolution images are required. For 10X Visium, this information is in `tissue_positions.csv`. 
spatial.locs1 = Seurat::GetTissueCoordinates(seu1, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs2 = Seurat::GetTissueCoordinates(seu2, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs <- rbind(spatial.locs1, spatial.locs2)
rownames(spatial.locs) <- colnames(data.input)

# Scale factors of spatial coordinates
# For 10X Visium, the conversion factor of converting spatial coordinates from Pixels to Micrometers can be computed as the ratio of the theoretical spot size (i.e., 65um) over the number of pixels that span the diameter of a theoretical spot size in the full-resolution image (i.e., 'spot_diameter_fullres' in pixels in the 'scalefactors_json.json' file). 
# Of note, the 'spot_diameter_fullres' factor is different from the `spot` in Seurat object and thus users still need to get the value from the original json file. 
scalefactors1 = jsonlite::fromJSON(txt = file.path("/Users/suoqinjin/Library/CloudStorage/OneDrive-Personal/works/CellChat/tutorial/spatial_imaging_data-intestinalA1", 'scalefactors_json.json'))
spot.size = 65 # the theoretical spot size (um) in 10X Visium
conversion.factor1 = spot.size/scalefactors1$spot_diameter_fullres
spatial.factors1 = data.frame(ratio = conversion.factor1, tol = spot.size/2)

scalefactors2 = jsonlite::fromJSON(txt = file.path("/Users/suoqinjin/Library/CloudStorage/OneDrive-Personal/works/CellChat/tutorial/spatial_imaging_data-intestinalA2", 'scalefactors_json.json'))
conversion.factor2 = spot.size/scalefactors2$spot_diameter_fullres
spatial.factors2 = data.frame(ratio = conversion.factor2, tol = spot.size/2)

spatial.factors <- rbind(spatial.factors1, spatial.factors2)
rownames(spatial.factors) <- c("A1", "A2")



### Create a CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
#> [1] "Create a CellChat object from a data matrix"
#> Create a CellChat object from spatial transcriptomics data... 
#> Set cell identities for the new CellChat object 
#> The cell groups used for CellChat analysis are  B-Cells, BEST4+/OTOP2+ Cell, Colonocytes, Crypt Top Colonocytes, Cycling Cells, DC1, DC2, Endothelial 1, Endothelial 2, Enteroendocrines, Glial, Goblets, ILCs, Macrophages&Monocytes, Mast Cells, Myofibroblasts, NK, Pericytes, Plasma Cells, Stromal 1, Stromal 2, Stromal 3, Stromal 4, T-Cells, Undifferentiated
cellchat
#> An object of class CellChat created from a single dataset 
#>  15609 genes.
#>  4965 cells. 
#> CellChat analysis of spatial data! The input spatial locations are 
#>                       x_cent y_cent
#> A1_AAACAAGTATCTCCCA-1   4372   5303
#> A1_AAACAGAGCGACTCCT-1   1753   4960
#> A1_AAACATTTCCCGGATT-1   5173   5097
#> A1_AAACCACTACACAGAT-1    949   5919
#> A1_AAACCCGAACGAAATC-1   4006   5846
#> A1_AAACCGGAAATGTTAA-1   4660   6225


### Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# set the used database in the object
cellchat@DB <- CellChatDB.use




# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 422
 
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
#> [1] 20.40818



# Part II: Inference of cell-cell communication network
# Compute the communication probability and infer cellular communication network
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 100)
#> truncatedMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on spatial transcriptomics data without distance values as constraints of the computed communication probability <<< [2024-02-22 12:06:15.985272]"
#> Molecules of the input L-R pairs are diffusible. Run CellChat in a diffusion manner based on the `interaction.range`.
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2024-02-22 12:09:10.370361]"



cellchat <- filterCommunication(cellchat, min.cells = 10)
#> The cell-cell communication related with the following cell groups are excluded due to the few number of cells:  Stromal 4 !     0.0% interactions are removed!


# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)


# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
#> [1] 188.3582

# We can also visualize the aggregated cell-cell communication network. For example, showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot or heatmap plot.
ptm = Sys.time()

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")



