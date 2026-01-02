


striatum_bams

#   seurat_obj <- 
striatum_bams_reproc <- striatum_bams %>%
    SCTransform(assay = "RNA") %>%
    RunPCA(npcs = 4, features = rownames(striatum_bams)) %>%
    RunUMAP(dims = 1:4)

striatum_bams_reproc


p_striatum_umap2 <- FeaturePlot(striatum_bams_reproc,
  features = c("H2-Aa", "Cd74", "Mrc1", "Flt1", "Pf4", "Cd163", "Gpnmb"),
  split.by = "condition",
  order = TRUE,
  pt.size	= 3,
  reduction = "umap"
)

ggsave(
  glue("{wkdir}/figures/temp/umap_gene_highlights_STRIATUM-BAMs.png"),
  p_striatum_umap2,
  width = 14,
  height = 21
)



enrichGO_res <- readRDS(
  glue(
    "{wkdir}/data/input/temp/",
    "BAM_EnrichGO_res.rds"
  )
)

enrichGO_res %>% glimpse


enrichGO_res %>%
  filter(grepl("regulation of leukocyte migration", Description)) %>%
  View


xlsx::write.xlsx(enrichGO_res, 
  glue("{wkdir}/data/processed/BAM_EnrichGO_res.xlsx"),
  sheetName = "GO-results",
  col.names = TRUE, row.names = TRUE, append = FALSE
)

