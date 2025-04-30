#______________________________________________________________________________
### ATTEMPT AT SENSITIVITY / SUBSAMPLING ANALYSIS OF HYP DATA
#______________________________________________________________________________

library(foreach)
library(doParallel)

abca_color_pal <- readRDS(glue("{wkdir}/data/interim/abca_color_pal.rds"))
merged_rois <- readRDS(
    glue(
      "{wkdir}/data/interim/",
      "merged_roi_seurat_filtered_staNMF-updated_2024-07-24.rds"
    )
)
meta_df <- merged_rois@meta.data %>% glimpse()

# ADDING GATING METADATA TO SEURAT OBJ

gated_suffixes <- c(
  "RT",
  "Piriform",
  "AMYG_BLA", 
  "AMYG_CT",
  "AMYG_CeMe",
  "CAUD",
  "THAL", 
  "HIPP",
  "SSC",
  "HYP"
)

meta_df_gated <- readRDS(
  glue(
    "{wkdir}/data/interim/",
    "gated-seurat-metadata_2024-12-31.rds"
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

merged_rois$final_gate <- meta_df_gated$final_gate
merged_rois$gated_cell_labels <- meta_df_gated$gated_cell_labels
merged_rois@meta.data %>% glimpse

#______________________________________

meta_df <- merged_rois@meta.data %>% as.data.table() %>% glimpse

# Analysis to determine if there is enrichment in cell types within gated regions
nested_cell_abund <- meta_df %>% 
  filter(final_gate != "") %>%
  group_by(slice_id, group, run, bregma, final_gate, singleR_labels) %>% 
  summarize(n = n(), .groups = "drop") %>%
  group_by(slice_id, final_gate) %>%
  mutate(rel_abundance = n / sum(n)) %>%
  ungroup() %>%
  group_by(final_gate, singleR_labels) %>% 
  nest()

gated_cell_type_prev_models_df <- nested_cell_abund %>% 
  mutate(count_model = purrr::map(data, possibly(function(df) lm(n ~ group + run + bregma, data = df), NULL))) %>% 
  mutate(rel_abund_model = purrr::map(data, possibly(function(df) lm(rel_abundance ~ group + run + bregma, data = df), NULL))) %>% 
  print(n = Inf)


# Extract model results for count models
count_model_results <- gated_cell_type_prev_models_df %>%
  filter(count_model != "NULL") %>%
  mutate(count_tidy = purrr::map(count_model, possibly(broom::tidy, NULL))) %>%
  select(final_gate, singleR_labels, count_tidy) %>%
  unnest(count_tidy) %>%
  mutate(model_type = "Cell Counts") %>% 
  glimpse()

# Extract model results for relative abundance models  
rel_abund_model_results <- gated_cell_type_prev_models_df %>%
  mutate(rel_abund_tidy = purrr::map(rel_abund_model, possibly(broom::tidy, NULL))) %>%
  select(final_gate, singleR_labels, rel_abund_tidy) %>%
  unnest(rel_abund_tidy) %>%
  mutate(model_type = "Relative Abundance") %>% 
  glimpse()

rel_abund_model_results %>% 
  filter(term == "groupwildr") %>% 
  View()

count_model_results %>% 
  filter(term == "groupwildr") %>% 
  View()

bothmods <- bind_rows(count_model_results, rel_abund_model_results)

# Combine results and create visualization
bothmods %>%
  mutate(colname = glue("{final_gate}_{singleR_labels}")) %>%
  mutate(signif = p.value <= 0.05) %>%
  # filter(term == "groupwildr") %>%
  ggplot(aes(
    x = estimate,
    y = fct_reorder(colname, estimate)
  )) +
  geom_bar(stat = "identity", width = 0.5, aes(fill = signif)) +
  geom_errorbar(aes(xmin = estimate - std.error, xmax = estimate + std.error),
                position = position_dodge(width = 0.8), width = 0.2) +
  theme_bw() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(
    x = "Beta Coefficient WildR / SPF",
    y = NULL,
    fill = "Significant"
  ) +
  facet_grid(model_type ~ term, space = "free", scales = "free") +
  theme(legend.position = "top")


p_count_mod <- bothmods %>%
  mutate(colname = glue("{final_gate}_{singleR_labels}")) %>%
  mutate(signif = p.value <= 0.05) %>%
  filter(term == "groupwildr") %>%
  filter(model_type == "Cell Counts") %>%
  ggplot(aes(
    x = estimate,
    y = fct_reorder(colname, estimate)
  )) +
  geom_bar(stat = "identity", width = 0.5, aes(fill = signif)) +
  theme_bw() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(
    x = "Beta Coefficient WildR / SPF",
    y = NULL,
    fill = "Significant"
  ) +
  facet_grid(model_type ~ term, space = "free", scales = "free") +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 6)) +
  scale_y_discrete(breaks = function(x) {
    # Only show labels for significant results
    sig_labels <- bothmods %>%
      filter(term == "groupwildr",
             model_type == "Cell Counts",
             p.value <= 0.05) %>%
      mutate(colname = glue("{final_gate}_{singleR_labels}")) %>%
      pull(colname)
    x[x %in% sig_labels]
  })

ggsave(
  glue("{wkdir}/figures/Celltype_prop_LM/",
    "celltype_counts_lm_results-counts-GATED_{Sys.Date()}.png"),
  p_count_mod,
  width = 4, height = 6
)

p_relab_mod <- bothmods %>%
  mutate(colname = glue("{final_gate}_{singleR_labels}")) %>%
  mutate(signif = p.value <= 0.05) %>%
  filter(term == "groupwildr") %>%
  filter(model_type == "Relative Abundance") %>%
  ggplot(aes(
    x = estimate,
    y = fct_reorder(colname, estimate)
  )) +
  geom_bar(stat = "identity", width = 0.5, aes(fill = signif)) +
  # geom_errorbar(aes(xmin = estimate - std.error, xmax = estimate + std.error),
  #               position = position_dodge(width = 0.8), width = 0.2) +
  theme_bw() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(
    x = "Beta Coefficient WildR / SPF",
    y = NULL,
    fill = "Significant"
  ) +
  facet_grid(model_type ~ term, space = "free", scales = "free") +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 6)) +
  scale_y_discrete(breaks = function(x) {
    # Only show labels for significant results
    sig_labels <- bothmods %>%
      filter(term == "groupwildr",
             model_type == "Relative Abundance",
             p.value <= 0.05) %>%
      mutate(colname = glue("{final_gate}_{singleR_labels}")) %>%
      pull(colname)
    x[x %in% sig_labels]
  })

ggsave(
  glue("{wkdir}/figures/Celltype_prop_LM/",
    "celltype_counts_lm_results-relab_GATED_{Sys.Date()}.png"),
  p_relab_mod,
  width = 4, height = 4
)




# mtcars_nested <- mtcars_nested %>% 
#   mutate(model = purrr::map(data, function(df) lm(mpg ~ wt, data = df)))
# mtcars_nested

#______________________________________


merged_rois_filt <- subset(merged_rois, slice_id %nin% c(
    "WILDR run1 roi3", "SPF run2 roi3",
    "SPF run3 roi1", "WILDR run3 roi2"
  )
)

# Checking cell labels
merged_rois_filt@meta.data$gated_cell_labels %>% table()
merged_rois_filt@meta.data$singleR_labels %>% table()

# Get cell counts per label
label_counts <- table(merged_rois_filt$gated_cell_labels)
# Get labels that have at least 500 cells
keep_labels <- names(label_counts[label_counts >= 200])
# Filter the Seurat object
merged_rois_filt <- subset(merged_rois_filt, gated_cell_labels %in% keep_labels)
merged_rois_filt

# Save filtered seurat object
saveRDS(merged_rois_filt, glue(
  "{wkdir}/data/interim/",
  "merged_roi_seurat_filtered_staNMF-updated",
  "_12-sample-filtered-gated-cell-types_2025-01-05.rds"), 
  compress = FALSE
)


# 763,451 to 759,645

counts <- merged_rois_filt@assays$SCT@counts %>% glimpse()
metadata <- merged_rois_filt@meta.data %>% glimpse()

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(merged_rois_filt@meta.data$gated_cell_labels)

# Create single cell experiment object
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = metadata
)

essential_meta_cols <- c(
  # "cluster_id" = "singleR_labels",
  "sample_id" = "slice_id",
  "cluster_id",
  "group", "run", "roi", "bregma"
  )

meta_df_grps <- dplyr::select(metadata, essential_meta_cols) %>% glimpse()
core_deseq_meta_df <- meta_df_grps %>% 
  dplyr::select(-cluster_id) %>%
  distinct() %>%
  glimpse()

count_df <- t(counts(sce)) %>% 
  as.data.table() %>%
  bind_cols(meta_df_grps) %>%
  glimpse()

count_df$cluster_id %>% table()

cois <- unique(metadata$cluster_id)

# c(
#   "12 HY GABA", "14 HY Glut", "13 CNU-HYa Glut", "09 CNU-LGE GABA"
# )

# Printing min/ max cell counts across slices
celltype_summary <- celltype_min  <- celltype_max <- list()

for (celltype in unique(metadata$cluster_id)) {
  message("\n", celltype)
  cellstats <- metadata %>%
    # NOTE SWITCHING TO GATED CELL TYPES
    filter(gated_cell_labels == celltype) %>%
    group_by(slice_id) %>%
    summarize(n = n()) %>%
    pull(n) %>%
    summary()

  cellstats %>% print()
  celltype_min[[celltype]] <- cellstats[[1]]
  celltype_max[[celltype]] <- cellstats[[6]]
  celltype_summary[[celltype]] <- cellstats
}

# 12 HY GABA
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     547    1455    1781    1875    2251    3621 

# 12 HY GABA
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1409    1681    2020    2136    2340    3621 


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



#' generate a list of 10 evenly dispersed values between the 
#' min and max cell type counts per slice 

sample_n_list <- list()
for (ct in cois){
  if ( is.na(celltype_summary[[ct]][[1]]) ) next
  if ( celltype_summary[[ct]][[3]] < 50 ) next
  min_val <- 10
  max_val <- celltype_summary[[ct]][[3]] # first quartile
  sample_n_list[[ct]] <- seq(min_val, max_val, length.out = 15) %>% print
}

# stats of cell types to be analyzed
cois <- names(sample_n_list)
for (celltype in cois) {
  message("\n", celltype)
  cellstats <- metadata %>%
    # NOTE SWITCHING TO GATED CELL TYPES
    filter(gated_cell_labels == celltype) %>%
    group_by(slice_id) %>%
    summarize(n = n()) %>%
    pull(n) %>%
    summary()

  cellstats %>% print()
  celltype_min[[celltype]] <- cellstats[[1]]
  celltype_max[[celltype]] <- cellstats[[6]]
  celltype_summary[[celltype]] <- cellstats
}

# remap cell types of interest by those passing abundance threshold
# a median of at least 50 across all slices

# Saving count data for future use
saveRDS(
  count_df,
  glue(
    "{wkdir}/data/interim/DESeq2/subsampling/",
    "gated-cell-types_12-samples_data/count_df_{Sys.Date()}.rds"
  )
)

# Saving sample n list for future use
saveRDS(
  sample_n_list,
  glue(
    "{wkdir}/data/interim/DESeq2/subsampling/",
    "gated-cell-types_12-samples_data/sample_n_list_{Sys.Date()}.rds"
  )
)

saveRDS(
  core_deseq_meta_df,
  glue(
    "{wkdir}/data/interim/DESeq2/subsampling/",
    "gated-cell-types_12-samples_data/core_deseq_meta_df_{Sys.Date()}.rds"
  )
)


#______________________________________________________________________________
# Manual single node loop below ---
#______________________________________________________________________________

# Set up parallel backend
registerDoParallel(cores = 8)
# cois_short <- c("12 HY GABA")

# Outer loop for cell types
for (clust in cois) {

  deseq2_clust_dir <- glue(
    "{wkdir}/data/interim/DESeq2/subsampling/",
    "gated-cell-types_12-samples_latest/{clust}"
  )

  # Middle loop for subsampling numbers
  for (subsamp_n in sample_n_list[[clust]]) {
    # Check if all iterations for this cluster and subsample size already exist
    all_iterations_exist <- all(sapply(1:50, function(niter) {
      deseq2_dir <- glue(
        "{deseq2_clust_dir}/n_{subsamp_n}/iter_{niter}"
      )
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

# Stop the parallel backend
stopImplicitCluster()


#______________________________________________________________________________
# Slurm batchtools loop  ---
#______________________________________________________________________________

run_deseq2_sensitivity_analysis <- function(wkdir, clust, nthreads) {
  
  require(glue)
  source(glue("{wkdir}/notebooks/R_scripts/pseudobulk-subsampling_funcs.R"))
  registerDoParallel(cores = nthreads)

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


sample_n_list_gated <- names(sample_n_list) %>% 
  str_detect(paste(gated_suffixes, collapse = "|")) %>% 
  sample_n_list[.]
names(sample_n_list_gated)

batchtools_params <- data.frame("clust" = names(sample_n_list_gated)) %>%
  mutate(nthreads = 4, wkdir = wkdir) %>% 
  glimpse()

# configure registry ----
cluster_run <- glue("{get_time()}_DESeq2_sensitivity_analysis")
message("\n\nRUNNING:  ", cluster_run, "\n")
breg <- makeRegistry(
  file.dir = glue("{wkdir}/.cluster_runs/", cluster_run),
  seed = 42
)
breg$cluster.functions <- batchtools::makeClusterFunctionsSlurm(
  template = glue("{wkdir}/batchtools_templates/batchtools.slurm_mamba-spatialomics.tmpl"),
  scheduler.latency = 0.1,
  fs.latency = 1
)

# Submit Jobs ----
jobs <- batchMap(
  fun = run_deseq2_sensitivity_analysis,
  args = batchtools_params,
  reg = breg
)
# jobs[, chunk := chunk(job.id, chunk.size = 1)]
# print(jobs[, .N, by = chunk])

submitJobs(jobs,
    resources = list(
        walltime = 120, # min (2hrs)
        memory = "50G", # MB (100 GB)
        ncpus = 4,
        max.concurrent.jobs = 9999
    )
)


#______________________________________________________________________________
# Compile results ---
#______________________________________________________________________________

subsamp_dir <- glue(
    "{wkdir}/data/interim/DESeq2/subsampling/",
    "gated_12-samples_group-run"
  )

res_paths <- list.files(
  subsamp_dir,
  pattern = "res_*",
  full.names = TRUE,
  recursive = TRUE
)

structured_res_paths <- res_paths %>%
  set_names(
    map_chr(., function(path) {
      path %>%
        strex::str_after_nth("/", -4) %>%
        strex::str_before_first("/")
    })
  ) %>%
  split(names(.)) %>%
  map(unname)


complied_dir <- glue(
    "{wkdir}/data/interim/DESeq2/subsampling/",
    "gated_12-samples_group-run_compiled"
  )
dir.create(complied_dir, showWarnings = FALSE, recursive = TRUE)


#' Function to aggregated interations of DESeq2 subsamples
process_celltype_results <- function(gated_celltype,
                                   structured_res_paths,
                                   complied_dir) {
  require(glue)
  require(magrittr)

  # Check if the file already exists
  output_file <- glue("{complied_dir}/{gated_celltype}_res_{Sys.Date()}.rds")
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

  saveRDS(res_df, glue("{complied_dir}/{gated_celltype}_res_{Sys.Date()}.rds"))
}

future::plan("multisession", workers = 28)
names(structured_res_paths) %>% 
  furrr::future_walk(
    ~ process_celltype_results(.x, structured_res_paths, complied_dir),
    .progress = TRUE
  )


# for (gated_celltype in names(structured_res_paths)) {

#   # Check if the file already exists
#   output_file <- glue("{complied_dir}/{gated_celltype}_res_{Sys.Date()}.rds")
#   if (file.exists(output_file)) {
#     message(glue("File already exists for {gated_celltype}. Skipping..."))
#     next
#   }
#   message(gated_celltype)
  
#   res_df <- structured_res_paths[[gated_celltype]] %>% 
#     set_names(strex::str_after_nth(., "/", -4)) %>%
#     purrr::map(
#       ~ readRDS(.x) %>%
#         as.data.frame() %>%
#         rownames_to_column(var = "gene_symbol")
#     ) %>%
#     bind_rows(.id = "info") %>%
#     separate(info, 
#       into = c("celltype", "subsamp_n", "iter", "date"), 
#       sep = "/", 
#       remove = FALSE
#     ) %>%
#     mutate(subsamp_n = as.numeric(gsub("n_", "", subsamp_n))) %>%
#     mutate(subsamp_n = as.numeric(round(subsamp_n)))

#   saveRDS(res_df, glue("{complied_dir}/{gated_celltype}_res_{Sys.Date()}.rds"))
# }

  # Create line plots for p-value, adjusted p-value, and log2FC
  plot_metrics <- function(data, y_var, y_label, title) {
    data %>%
      arrange(HIT) %>%
      ggplot(aes(x = subsamp_n, y = !!sym(y_var), 
        color = HIT, alpha = HIT, group = gene_symbol
      )) +
      geom_line() +
      geom_hline(yintercept = ifelse(y_var == "median_log2fc", 0, -log10(0.05)),
        linetype = "dashed", color = "red"
      ) +
      scale_color_manual(values = c("TRUE" = "darkred", "FALSE" = "grey")) +
      scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4)) +
      labs(
        x = "Subsample Size",
        y = y_label,
        title = title,
        color = "HIT"
      ) +
      theme_bw() +
      theme(legend.position = "none")
  }




for (celltype_gate in names(structured_res_paths)) {
  message(celltype_gate)

  res <- readRDS(
    glue("{complied_dir}/{celltype_gate}_res_{Sys.Date()}.rds")
  )
  res %>% glimpse

  summary_stats <- res %>%
    group_by(gene_symbol, subsamp_n) %>%
    summarize(
      median_pvalue = median(pvalue, na.rm = TRUE),
      mean_pvalue = mean(pvalue, na.rm = TRUE),
      sd_pvalue = sd(pvalue, na.rm = TRUE),
      median_adj_pvalue = median(padj, na.rm = TRUE),
      mean_adj_pvalue = mean(padj, na.rm = TRUE),
      sd_adj_pvalue = sd(padj, na.rm = TRUE),
      median_log2fc = median(log2FoldChange, na.rm = TRUE),
      mean_log2fc = mean(log2FoldChange, na.rm = TRUE),
      sd_log2fc = sd(log2FoldChange, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    group_by(gene_symbol) %>%
    mutate(HIT = any(median_pvalue < 0.05)) %>%
    ungroup() %>% 
    mutate(
      neg_log_pval = -log10(median_pvalue),
      neg_log_adj_pval = -log10(median_adj_pvalue)
    )

  plots <- list(
    p_median_pvalue = plot_metrics(summary_stats, 
      "neg_log_pval",
      "-log10(Median P-Value)",
      paste(celltype_gate, "\nMedian P-Value vs Subsample Size")
    ),
    p_median_adj_pvalue = plot_metrics(summary_stats,
      "neg_log_adj_pval",
      "-log10(Median Adjusted P-Value)",
      paste(celltype_gate, "\nMedian Adjusted P-Value vs Subsample Size")
    ),
    p_median_log2fc = plot_metrics(summary_stats,
      "median_log2fc",
      "Median log2FC",
      paste(celltype_gate, "\nMedian log2FC vs Subsample Size")
    )
  )

  # Stack the plots vertically using aplot
  stacked_plot <- aplot::plot_list(
    plots$p_median_pvalue,
    plots$p_median_adj_pvalue,
    plots$p_median_log2fc,
    ncol = 1,
    heights = c(1, 1, 1)
  )

  # Save the stacked plot
  ggsave(
    filename = glue(
      "{wkdir}/figures/pseudobulk-subsampling/12_sample-run-group_{Sys.Date()}/",
      "{celltype_gate}_subsamp_size_summary_{Sys.Date()}.png"
    ),
    plot = stacked_plot,
    width = 10,
    height = 12,
    dpi = 300
  )
}



# filter for genes with significant median p-values
pseudobulk_stat_summary_df <- tibble()
for (celltype_gate in names(structured_res_paths)) {
  message(celltype_gate)

  res <- readRDS(
    glue("{complied_dir}/{celltype_gate}_res_{Sys.Date()}.rds")
  )

  summary_stats <- res %>%
    group_by(gene_symbol, subsamp_n) %>%
    summarize(
      median_pvalue = median(pvalue, na.rm = TRUE),
      mean_pvalue = mean(pvalue, na.rm = TRUE),
      sd_pvalue = sd(pvalue, na.rm = TRUE),
      median_adj_pvalue = median(padj, na.rm = TRUE),
      mean_adj_pvalue = mean(padj, na.rm = TRUE),
      sd_adj_pvalue = sd(padj, na.rm = TRUE),
      median_log2fc = median(log2FoldChange, na.rm = TRUE),
      mean_log2fc = mean(log2FoldChange, na.rm = TRUE),
      sd_log2fc = sd(log2FoldChange, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    group_by(gene_symbol) %>%
    mutate(HIT = any(median_pvalue < 0.05)) %>%
    ungroup() %>% 
    mutate(region = celltype_gate)

  pseudobulk_stat_summary_df %<>% bind_rows(summary_stats)
}

gated_12-samples_group-run_misc_data

saveRDS(
  pseudobulk_stat_summary_df,
  glue("{wkdir}/data/interim/DESeq2/subsampling/",
    "gated_12-samples_group-run_misc-data/",
    "pseudobulk_stat_summary_df_{Sys.Date()}.rds"
  )
)



pseudobulk_stat_summary_df <- readRDS(
  glue("{wkdir}/data/interim/DESeq2/subsampling/",
    "gated_12-samples_group-run_misc-data/",
    "pseudobulk_stat_summary_df_2025-01-05.rds"
  )
)


pseudobulk_stat_summary_df %>% glimpse
#  create a summary df with the mean of the median p-values and log2fc
pseudobulk_stat_summary_df_avg <- pseudobulk_stat_summary_df %>%
  group_by(region, gene_symbol) %>%
  summarize(
    mean_median_pvalue = mean(median_pvalue, na.rm = TRUE),
    mean_median_log2fc = mean(median_log2fc, na.rm = TRUE),
    .groups = 'drop'
  )


# Create volcano plots for each region
volcano_plots <- pseudobulk_stat_summary_df_avg %>%
  group_by(region) %>%
  group_map(~ {
    plot <- ggplot(.x, 
      aes(x = mean_median_log2fc, y = -log10(mean_median_pvalue))
    ) +
      geom_point(
        aes(color = (abs(mean_median_log2fc) > 0.2 & 
          -log10(mean_median_pvalue) > -log10(0.05))
        ),
        alpha = 0.6
      ) +
      geom_hline(
        yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      geom_vline(
        xintercept = c(-0.2, 0.2), linetype = "dashed", color = "red") +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
      labs(
        title = paste("Volcano Plot -", .y),
        x = "Mean Median log2FC",
        y = "-log10(Mean Median p-value)"
      ) +
      theme_bw() +
      theme(legend.position = "none") +
      ggrepel::geom_text_repel(
        data = ~ .x %>% 
          filter(abs(mean_median_log2fc) > 0.2 & 
          -log10(mean_median_pvalue) > -log10(0.05)) %>%
          mutate(text_label = glue("{gene_symbol}\n{region}")),
        aes(label = text_label),
        size = 3,
        box.padding = 0.5,
        point.padding = 0.1,
        segment.color = "grey50"
      )

    setNames(list(plot), .y)
  }, .keep = TRUE) %>%
  unlist(recursive = FALSE)


# Save individual volcano plots
for (i in seq_along(volcano_plots)) {
  region_name <- names(volcano_plots)[i]
  ggsave(
    filename = glue(
      "{wkdir}/figures/pseudobulk-subsampling/",
      "12_sample-run-group/volcano_plots/",
      "volcano_plot_{region_name}_{Sys.Date()}.png"),
    plot = volcano_plots[[i]],
    width = 10,
    height = 8,
    dpi = 300
  )
}


volcano_plot_main <- pseudobulk_stat_summary_df_avg %>%
  mutate(text_label = glue("{gene_symbol}\n{region}")) %>%
  ggplot(aes(x = mean_median_log2fc, y = -log10(mean_median_pvalue))) +
  geom_point(alpha = 0.6, color = "grey") +
  geom_point(
    data = . %>% filter(abs(mean_median_log2fc) > 0.2 & -log10(mean_median_pvalue) > -log10(0.05)),
    aes(color = gene_symbol),
    alpha = 1
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "red"
  ) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "red") +
  scale_color_manual(values = ggsci::pal_igv()(n_distinct(pseudobulk_stat_summary_df_avg$gene_symbol[abs(pseudobulk_stat_summary_df_avg$mean_median_log2fc) > 0.2 & -log10(pseudobulk_stat_summary_df_avg$mean_median_pvalue) > -log10(0.05)]))) +
  labs(
    title = "Volcano Plot - All Regions",
    x = "Mean Median log2FC",
    y = "-log10(Mean Median p-value)"
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  ggrepel::geom_text_repel(
    data = . %>% filter(abs(mean_median_log2fc) > 0.2 & -log10(mean_median_pvalue) > -log10(0.05)),
    aes(label = text_label),
    size = 1.5,
    max.overlaps = 20,
    box.padding = 0.2,
    point.padding = 0.1,
    segment.color = "grey50"
  )


# Save the combined plot
ggsave(
  filename = glue("{wkdir}/figures/pseudobulk-subsampling/12_sample-run-group/",
    "combined_volcano_plot_{Sys.Date()}.png"
  ),
  plot = volcano_plot_main,
  width = 10,
  height = 7,
  dpi = 300
)



volcano_plot_tissue <- pseudobulk_stat_summary_df_avg %>%
  mutate(text_label = glue("{gene_symbol}\n{region}")) %>%
  ggplot(aes(x = mean_median_log2fc, y = -log10(mean_median_pvalue))) +
  geom_point(alpha = 0.6, color = "grey") +
  geom_point(
    data = . %>% filter(abs(mean_median_log2fc) > 0.2 & -log10(mean_median_pvalue) > -log10(0.05)),
    aes(color = region),
    alpha = 1
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "red"
  ) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "red") +
  scale_color_manual(values = 
    ggsci::pal_igv()(n_distinct(pseudobulk_stat_summary_df_avg$region[abs(pseudobulk_stat_summary_df_avg$mean_median_log2fc) > 0.2 & -log10(pseudobulk_stat_summary_df_avg$mean_median_pvalue) > -log10(0.05)]))) +
  labs(
    title = "Volcano Plot - All Regions",
    x = "Mean Median log2FC",
    y = "-log10(Mean Median p-value)"
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  ggrepel::geom_text_repel(
    data = . %>% filter(abs(mean_median_log2fc) > 0.2 & -log10(mean_median_pvalue) > -log10(0.05)),
    aes(label = text_label),
    size = 1.5,
    max.overlaps = 20,
    box.padding = 0.2,
    point.padding = 0.1,
    segment.color = "grey50"
  )


# Save the combined plot
ggsave(
  filename = glue("{wkdir}/figures/pseudobulk-subsampling/12_sample-run-group/",
    "combined_volcano_plot-tissuecol_{Sys.Date()}.png"
  ),
  plot = volcano_plot_tissue,
  width = 10,
  height = 7,
  dpi = 300
)

p_region_bars <- pseudobulk_stat_summary_df_avg %>% 
  filter(abs(mean_median_log2fc) > 0.2 & mean_median_pvalue < 0.05) %>%
  ggplot(aes(y = region, x = mean_median_log2fc, fill = gene_symbol)) +
  geom_bar(stat = "identity", position = "stack") +
  ggsci::scale_fill_igv() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Region", x = "Mean Median Log2FC", fill = "Gene")

ggsave(
  filename = glue(
    "{wkdir}/figures/pseudobulk-subsampling/12_sample-run-group/",
    "region_bars_{Sys.Date()}.png"
  ),
  plot = p_region_bars,
  width = 8,
  height = 6,
  dpi = 300
)

p_region_bars_location <- pseudobulk_stat_summary_df_avg %>% 
  mutate(location = strex::str_after_last(region, " ")) %>%
  filter(abs(mean_median_log2fc) > 0.2 & mean_median_pvalue < 0.05) %>%
  ggplot(aes(y = region, x = mean_median_log2fc, fill = gene_symbol)) +
  geom_bar(stat = "identity", position = "stack") +
  ggsci::scale_fill_igv() + 
  facet_grid(location~., scales = "free", space = "free") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 0)
  ) +
  labs(y = "Region", x = "Mean Median Log2FC", fill = "Gene")

ggsave(
  filename = glue(
    "{wkdir}/figures/pseudobulk-subsampling/12_sample-run-group/",
    "region_bars_location_{Sys.Date()}.png"
  ),
  plot = p_region_bars_location,
  width = 8,
  height = 6,
  dpi = 300
)




# Filter hits_df_avg for genes with mean_median_pvalue < 0.05
significant_genes <- hits_df_avg %>%
  filter(mean_median_pvalue < 0.05) %>% 
  arrange(mean_median_pvalue)

gene <- significant_genes$gene_symbol[1]
region <- significant_genes$region[1]

# Function to create plots for each gene
create_gene_plots <- function(gene, region) {
  # Read the full data from the correct RDS file
  gene_data <- readRDS(glue("{complied_dir}/{region}_res_2024-09-19.rds")) %>% 
    filter(gene_symbol == gene)
  
  # Common theme for all plots
  common_theme <- theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Plot for p-value
  p_pvalue <- gene_data %>%
    ggplot(aes(x = as.numeric(subsamp_n), y = -log10(pvalue))) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(x = "Subsample Size", y = "-log10(p-value)") +
    common_theme +
    scale_x_continuous(breaks = unique(gene_data$subsamp_n))

  # Plot for adjusted p-value
  p_padj <- gene_data %>%
    ggplot(aes(x = as.numeric(subsamp_n), y = -log10(padj))) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(x = "Subsample Size", y = "-log10(adj. p-value)") +
    common_theme +
    scale_x_continuous(breaks = unique(gene_data$subsamp_n))

  # Plot for log2 fold change
  p_log2fc <- gene_data %>%
    ggplot(aes(x = as.numeric(subsamp_n), y = log2FoldChange)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    labs(x = "Subsample Size", y = "log2 Fold Change") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_continuous(breaks = unique(gene_data$subsamp_n))

  # Combine plots using patchwork
  combined_plot <- (p_pvalue / p_padj / p_log2fc) +
    plot_layout(heights = c(1, 1, 1.2)) +
    plot_annotation(
      title = glue("{gene} - {region}"),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    )

  # Save the combined plot
  ggsave(
    filename = glue(
      "{wkdir}/figures/pseudobulk-subsampling/12_sample-run-group/top_hits/",
      "{gene}_{region}_pval_padj_log2fc_{Sys.Date()}.png"
    ),
    plot = combined_plot,
    width = 6,
    height = 8,
    dpi = 300
  )
}

# Apply the function to each significant gene
purrr::walk2(
  significant_genes$gene_symbol,
  significant_genes$region,
  ~create_gene_plots(.x, .y)
)





#_________
# Seurat Object needs to have gated cell types metadata column

# Function to plot gene epxression on spatial data
plot_xy_gene_exp <- function(seur, gene) {
  gene = toupper(gene)
  meta_df <- merged_rois@meta.data
  meta_df <- bind_cols(meta_df, 
    genename = merged_rois@assays$SCT@data[gene, ]) %>% 
    dplyr::mutate(slice_id = 
      glue("{toupper(group)} {run} {roi}")) %>%
    arrange(group, run, genename)
    
  p <- meta_df %>%
    ggplot(aes(x = center_x, y = center_y)) +
    geom_point(aes(color = genename), size = 0.2, alpha = 0.4) +
    coord_fixed() +
    facet_wrap(~slice_id, ncol = 4) +
    labs(color = gene) +
    scale_color_viridis_c(option = "magma") +
    theme_bw()
  return(p)
}

# Function to plot gene epxression on spatial data
plot_celltype_xy_gene_exp <- function(seur, gene, cell_filter, cell_type_column = "singleR_labels") {
  # Check if the cell_type_column exists in the metadata
  if (!cell_type_column %in% colnames(seur@meta.data)) {
    stop(glue("Column '{cell_type_column}' not found in the metadata."))
  }
  
  gene = toupper(gene)
  meta_df <- seur@meta.data
  meta_df <- bind_cols(meta_df, 
    genename = seur@assays$SCT@data[gene, ]) %>%
    mutate(celltype_of_interest = 
      ifelse(!!sym(cell_type_column) == cell_filter, TRUE, FALSE)) %>%
    arrange(group, run, genename)

  p <- meta_df %>%
    ggplot(aes(x = center_x, y = center_y)) +
    geom_point(data = dplyr::filter(meta_df, !celltype_of_interest),
        color = "lightgrey", size = 0.3, alpha = 0.5) +
    geom_point(data = dplyr::filter(meta_df, celltype_of_interest),
        aes(color = genename), size = 0.3, alpha = 0.8) +
    facet_wrap(~slice_id, ncol = 4) +
    coord_fixed() +
    scale_y_reverse() +
    labs(color = gene, title = glue("Plotting {gene} in {cell_filter}")) +
    scale_color_viridis_c(option = "viridis") +
    theme_bw()
  return(p)
}






pseudobulk_stat_summary_df_avg %>% glimpse()
plot_deg_df <- pseudobulk_stat_summary_df_avg %>% 
  filter(abs(mean_median_log2fc) > 0.2 & 
    -log10(mean_median_pvalue) > -log10(0.05))  %>% 
  select(cell_type = region, gene = gene_symbol) %>%
  distinct() %>% 
  glimpse()

gate_list <- c("AMYG_", "Piriform", "RT")

# plot_deg_df <- deg_all_df %>%
#   filter(threshold) 

for (ii in 1:nrow(plot_deg_df)){
  gid <- plot_deg_df$gene[ii]
  coi <- plot_deg_df$cell_type[ii]
  # message(glue("Plotting: {gid} - {coi}"))

  figout_all <- glue("{wkdir}/figures/gene_exp_scatter/2024-08-23/",
    "{gid}_expression_SCT_2024-08-23.png"
  )
  figout_celltype <- glue("{wkdir}/figures/gene_exp_scatter/2024-09-21/",
    "{gid}_{coi}_expression_SCT_2024-08-23.png"
  )

  # if (!file.exists(figout_all)) {
  #   p <- plot_xy_gene_exp(merged_rois, gid)
  #   ggsave(filename = figout_all, p,
  #     width = 16, height = 16
  #   )
  # }

  if (!file.exists(figout_celltype)) {

    # check if the cell type is in the metadata
    if (any(str_detect(coi, gate_list))) {
      print(glue("Plotting: GATED REGION {gid} - {coi}"))
      p2 <- plot_celltype_xy_gene_exp(merged_rois, gid, coi, 
        cell_type_column = "gated_cell_labels"
      )
      ggsave(filename = figout_celltype, p2,
        width = 16, height = 16
      )
    } else {
      print(glue("Plotting: UNGATED REGION {gid} - {coi}"))
      p2 <- plot_celltype_xy_gene_exp(merged_rois, gid, coi, 
        cell_type_column = "singleR_labels"
      )
      ggsave(filename = figout_celltype, p2,
        width = 16, height = 16
      )
    }

  }
}






any(str_detect(coi, gate_list))











gois <- c(
  # "CAMK2N1",
  "CALM2"
  # "SLC1A2",
  # "SLC17A6",
  # "APP",
  # "GAP43",
  # "APOE",
  # "CALM1",
  # "CALM3",
  # "SNAP25",
  # "ATP1A3"
)

# for (gid in gois){
#   message(glue("Plotting: {gid}"))
#   figout <- glue("{wkdir}/figures/gene_exp_scatter/{gid}_expression_SCT_{Sys.Date()}.png")
#   # if (file.exists(figout)) next
  
#   p <- plot_xy_gene_exp(merged_rois, gid)
#   ggsave(
#     filename = figout,
#     p,
#     width = 22, height = 20
#   )
# }










# full_res_df <- res_paths %>%
#   keep(grepl(" RT", .)) %>%
#   # keep(grepl("AMYG_|Piriform| RT", .)) %>%
#   set_names(strex::str_after_nth(., "/", -4)) %>%
#   purrr::map(
#     ~ readRDS(.x) %>%
#       as.data.frame() %>%
#       rownames_to_column(var = "gene_symbol") %>%
#       glimpse()
#   ) %>%
#   bind_rows(.id = "info") %>%
#   separate(info, 
#     into = c("celltype", "subsamp_n", "iter", "date"), 
#     sep = "/", 
#     remove = FALSE
#   ) %>%
#   mutate(subsamp_n = as.numeric(gsub("n_", "", subsamp_n))) %>%
#   mutate(subsamp_n = as.numeric(round(subsamp_n))) %>%
#   glimpse()


# # Group the full res dataframe and count significant genes
# grouped_counts <- full_res_df %>%
#   group_by(info, subsamp_n) %>%
#   summarize(
#     sig_genes = sum(padj <= 0.05 & abs(log2FoldChange) >= 0.2, na.rm = TRUE),
#     .groups = 'drop'
#   ) %>%
#   mutate(subsamp_n = as.numeric(str_extract(subsamp_n, "\\d+")))

# # Create the ggplot2 scatterplot
# p_sig_genes <- ggplot(grouped_counts, aes(x = subsamp_n, y = sig_genes)) +
#   geom_point(position = position_jitter(width = 0.2, height = 0)) +
#   # geom_boxplot() +
#   # scale_x_log10() +
#   labs(
#     x = "Subsample Size",
#     y = "Number of Significant Genes",
#     title = "Significant Genes vs Subsample Size",
#     color = "Sample Info"
#   ) +
#   theme_bw() +
#   theme(legend.position = "none")

# ggsave(
#   filename = glue("{wkdir}/figures/deseq2_sensitivity/00_sig_genes_subsamp_size_{Sys.Date()}.png"),
#   plot = p_sig_genes,
#   width = 8,
#   height = 6,
#   dpi = 300
# )




# full_res_df %>%
#   filter(celltype == "12 HY GABA") %>%
#   filter(padj <= 0.05) %>%
#   DT::datatable()


# full_res_df %>%
#   filter(celltype == "12 HY GABA") %>%
#   filter(gene_symbol == "OTP") %>%
#   filter(subsamp_n == 1047) %>%
#   # filter(padj <= 0.05) %>%
#   DT::datatable()

# gois <- c(
#   "OTP",
#   "TMEM132D",
#   "GPX3",
#   "C1QA",
#   "ALDH1L1",
#   "STAT3",
#   "P2RY12",
#   "CALCR",
#   "GABRQ",
#   "GPR101",
#   "SELPLG",
#   "STX1A",
#   "SLC1A3",
#   "RASAL1",
#   "PCP4L1",
#   "SCN4B"
# )


# for (gene in gois) {
#   # Plot for p-value
#   p_pvalue <- full_res_df %>% 
#     filter(gene_symbol == gene) %>%
#     ggplot(aes(x = factor(subsamp_n), y = -log10(pvalue))) +
#     geom_point(size = 0.5, alpha = 0.5, position = position_jitter(width = 0.1, height = 0)) +
#     geom_boxplot(alpha = 0.2, outlier.shape = NA) +
#     geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
#     labs(
#       x = "Subsample Size",
#       y = "-log10(p-value)"
#     ) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))

#   # Plot for log2 fold change
#   p_log2fc <- full_res_df %>% 
#     filter(gene_symbol == gene) %>%
#     ggplot(aes(x = factor(subsamp_n), y = log2FoldChange)) +
#     geom_point(size = 0.5, alpha = 0.5, position = position_jitter(width = 0.1, height = 0)) +
#     geom_boxplot(alpha = 0.2, outlier.shape = NA) +
#     scale_y_log10() +
#     # geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
#     labs(
#       x = "Subsample Size",
#       y = "log2 Fold Change"
#     ) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))

#   # Combine plots using patchwork
#   combined_plot <- p_pvalue / p_log2fc +
#     plot_layout(ncol = 1, heights = c(1, 1)) +
#     plot_annotation(title = gene)

#   # Save the combined plot
#   ggsave(
#     filename = glue("{wkdir}/figures/deseq2_sensitivity/run-group_{Sys.Date()}/{gene}_pval_log2fc.png"),
#     plot = combined_plot,
#     width = 10,
#     height = 10,
#     dpi = 300
#   )
# }

