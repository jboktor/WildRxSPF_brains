
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
    glue("{wkdir}/data/interim/",
        "gated-seurat-metadata_2024-12-31.rds")
    )

meta_df_gated %>% glimpse()

# create a new column with cell types and gated regions & convert pixels to micrometers* 0.106
mini_meta <- meta_df_gated %>%
  mutate(x_microns = center_x * 0.106, y_microns = center_y * 0.106) %>%
  dplyr::select(cell_id, slice_id, group, run, bregma, final_gate, singleR_labels, x_microns, y_microns) %>%
  glimpse()

# mini_meta %>% glimpse
# cell_type <- "34 Immune"
# gate <- "Piriform"
# sid <- "WILDR run4 roi3"



calculate_cell_neighborhoods <- function(sid, cell_type, gate, wkdir) {
    require(glue)
    require(magrittr)
    require(dplyr)
    require(data.table)
    require(tibble)
    require(tidyr)
    require(tidyverse)

    start_time <- Sys.time()
    message(glue("Starting at {start_time}"))
    message(glue("Processing {sid} - {cell_type} in {gate}"))

    output_dir <- glue(
        "{wkdir}/data/interim/cell_neighborhoods/",
        "{sid}/{gate}"
        )
    output_path <- glue("{output_dir}/cell_neighborhood_stats_{sid}_{gate}_{cell_type}.rds")
    
    if (file.exists(output_path)) {
        message("Output file already exists, skipping...")
        return(NULL)
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    message(glue("{Sys.time()}: Loading metadata..."))
    meta_df_gated <- readRDS(
        glue("{wkdir}/data/interim/",
            "gated-seurat-metadata_2024-12-31.rds")
        )
    # create a new column with cell types and gated regions
    # convert pixels to micrometers* 0.106
    df <- meta_df_gated %>%
        mutate(x_microns = center_x * 0.106,
               y_microns = center_y * 0.106) %>%
        dplyr::select(cell_id, slice_id, group, run, bregma,
                     final_gate, singleR_labels, x_microns, y_microns)
    
    # subset dataset
    message(glue("{Sys.time()}: Subsetting data..."))
    cell_subset_df <- df %>%
        dplyr::filter(final_gate == gate) %>%
        dplyr::filter(slice_id == sid)

    message(glue("{Sys.time()}: Calculating distances for {nrow(cell_subset_df)} cells..."))
    dist_matrix <- sqrt(
        (outer(cell_subset_df$x_microns, cell_subset_df$x_microns, "-"))^2 + 
        (outer(cell_subset_df$y_microns, cell_subset_df$y_microns, "-"))^2
    )
    dist_dt <- data.table::as.data.table(dist_matrix)
    colnames(dist_dt) <- rownames(dist_dt) <- cell_subset_df$cell_id

    message(glue("{Sys.time()}: Reformatting distance matrix and adding metadata..."))
    dist_dt_long <- dist_dt %>% 
        tibble::rownames_to_column("cell1_id") %>%
        tidyr::pivot_longer(
            cols = -cell1_id,
            names_to = "cell2_id", 
            values_to = "distance"
        ) %>%
        # Join metadata for both cells
        dplyr::left_join(
            cell_subset_df %>% dplyr::select(cell_id, singleR_labels, run, group),
            by = c("cell1_id" = "cell_id")
        ) %>%
        dplyr::rename(cell1_type = singleR_labels) %>%
        dplyr::left_join(
            cell_subset_df %>% dplyr::select(cell_id, singleR_labels),
            by = c("cell2_id" = "cell_id")
        ) %>%
        dplyr::rename(cell2_type = singleR_labels)

    message(glue("{Sys.time()}: Calculating neighborhood statistics..."))
    distb_df <- dist_dt_long %>% 
        dplyr::filter(cell1_type == cell_type) %>% 
        dplyr::filter(cell1_id !=  cell2_id) %>%
        dplyr::filter(distance <= 60) %>% 
        dplyr::group_by(cell1_id, cell1_type, cell2_type) %>%
        dplyr::summarize(n = n()) %>% 
        dplyr::group_by(cell1_id) %>% 
        dplyr::mutate(rel_abund = n / sum(n))

    message(glue("{Sys.time()}: Saving results..."))
    saveRDS(distb_df, output_path)
    
    end_time <- Sys.time()
    elapsed <- round(difftime(end_time, start_time, units="mins"), 2)
    message(glue("Done! Completed in {elapsed} minutes"))
}


batchtools_params <- expand_grid(
    sid = unique(mini_meta$slice_id),
    cell_type = unique(mini_meta$singleR_labels),
    gate = gated_suffixes,
    wkdir = wkdir
) %>% 
    dplyr::mutate(
        output_path = glue(
            "{wkdir}/data/interim/cell_neighborhoods/{sid}/{gate}/",
            "cell_neighborhood_stats_{sid}_{gate}_{cell_type}.rds")
    ) %>%
    dplyr::filter(!file.exists(output_path)) %>% 
    dplyr::select(-c(output_path)) %>% 
    glimpse()


# future::plan("multisession", workers = 30)
# purrr::pmap(params, calculate_cell_neighborhoods, .progress = TRUE)


# configure registry ----
cluster_run <- glue("{get_time()}_cell_neighborhoods")
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
  fun = calculate_cell_neighborhoods,
  args = batchtools_params, #%>% filter(cell_type == "01 IT-ET Glut", gate == "SSC"),
  reg = breg
)

jobs[, chunk := chunk(job.id, chunk.size = 50)]
print(jobs[, .N, by = chunk])

submitJobs(jobs,
    resources = list(
        walltime = 120, # min (1hrs)
        memory = "400GB", # MB (100 GB)
        ncpus = 2,
        max.concurrent.jobs = 9999
    )
)


#------------------------------------------------------------------------------
# load_cell_neighborhood_data ----
#------------------------------------------------------------------------------

load_cell_neighborhood_data <- function(res_paths, wkdir) {
    res_paths %>%
        purrr::set_names() %>%
        furrr::future_map_dfr(function(path, wkdir) {
            readRDS(glue("{wkdir}/data/interim/cell_neighborhoods/{path}")) %>%
                mutate(file_path = path)
        }, wkdir = wkdir, .progress = TRUE) %>%
        data.table::as.data.table() %>%
        dplyr::mutate(
            slice_id = strex::str_before_first(file_path, "/"),
            cell_type = strex::str_after_last(file_path, "_") %>%
                gsub(".rds", "", .),
            gating = strex::str_before_nth(file_path, "/", 2) %>%
                strex::str_after_last(., "/"),
            group = strex::str_before_first(slice_id, " "),
            run = strex::str_after_first(slice_id, " ") %>%
                strex::str_before_first(., " ")
        )
}

poor_quality_slices = c(
  "SPF run3 roi4", 
  "WILDR run1 roi3", 
  "SPF run2 roi3", 
  "WILDR run3 roi2"
  )

# library(furrr)

all_paths <- list.files(
    glue("{wkdir}/data/interim/cell_neighborhoods/"),
    recursive = TRUE
    )
length(all_paths)

cell_neighbors_df <- load_cell_neighborhood_data(all_paths, wkdir) %>% 
    glimpse()

saveRDS(
    cell_neighbors_df,
    glue("{wkdir}/data/interim/cell_neighborhoods_summary/",
        "cell_neighborhood_stats_all_slices.rds"),
    compress = FALSE
)


#------------------------------------------------------------------------------
# Load data ----
#------------------------------------------------------------------------------

library(broom.mixed)

cell_neighbors_df <- readRDS(
    glue("{wkdir}/data/interim/cell_neighborhoods_summary/",
        "cell_neighborhood_stats_all_slices.rds")
)

cell_neighbors_df %<>% filter(slice_id %nin% poor_quality_slices)
cell_neighbors_df %>% glimpse()
cell_neighbors_df$gating %>% table

# fit models for each cell type and gating 
models_df <- cell_neighbors_df %>%
    nest(data = -c(gating, cell1_type, cell2_type)) %>%
    mutate(model = purrr::map(data,
        ~ purrr::possibly(
            ~ nlme::lme(log1p(rel_abund) ~ group, 
                random = ~1|run/slice_id, data = na.omit(.x)),
            otherwise = NULL
        )(.))
    ) %>%
    glimpse()

# models_df$model[[1]] %>% summary()
# models_df$model[[1]] %>% broom.mixed::tidy()

models_summary_df <- models_df %>% 
    mutate(model_summary = purrr::map(model, ~ broom.mixed::tidy(.x))) %>%
    unnest(model_summary) %>%
    janitor::clean_names() %>%
    glimpse()

saveRDS(
    models_summary_df,
    glue("{wkdir}/data/interim/cell_neighborhoods_summary/",
        "cell_neighborhood_stats_all_slices_models.rds")
)

# Visualizing model results ----

models_summary_df <- readRDS(
    glue("{wkdir}/data/interim/cell_neighborhoods_summary/",
        "cell_neighborhood_stats_all_slices_models.rds")
)
models_summary_df %>% glimpse()

models_summary_df %>% 
    dplyr::select(-c(data, model)) %>%
    filter(term == "groupWILDR") %>% 
    mutate(sig = p_value <= 0.05) %>% 
    filter(sig) %>% 
    View()

p_mod_summary <- models_summary_df %>%
    dplyr::select(-c(data, model)) %>%
    filter(!is.na(p_value)) %>%
    mutate(sig = p_value <= 0.05) %>%
    filter(term == "groupWILDR") %>%
    # Arrange data so significant points are plotted last
    arrange(sig) %>%
    ggplot(aes(x=estimate, y= fct_reorder(cell2_type, p_value))) +
    geom_pointrange(aes(xmin = estimate - std_error, xmax = estimate + std_error, color = sig), 
        size = 0.4) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    facet_grid(cell1_type ~ gating, scales = "free", space = "free") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(angle = 0)
    ) +
    labs(y = "Cell Type", x = "Hierarchical Linear Mixed Model Estimate")

ggsave(
    glue("{wkdir}/figures/cell_neighborhoods/misc/",
        "All_rel_abund_lmm_models.png"), 
    p_mod_summary, width = 16, height = 50, limitsize = FALSE
)



p_dist <- cell_neighbors_df %>% 
    # filter(gating == "Piriform") %>% 
    ggplot(aes(x = fct_reorder(cell2_type, rel_abund), y = rel_abund)) +
    geom_boxplot(aes(fill = group)) +
    facet_wrap(~cell1_type) +
    theme_bw() +
    scale_fill_d3() +
    facet_wrap(~gating, ncol = 1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Relative Abundance of Immune Cells in Piriform Cortex",
         x = "Cell Type",
         y = "Relative Abundance")

ggsave(glue("{wkdir}/figures/cell_neighborhoods/Piriform/",
    "piriform_immune_rel_abund.png"), 
    p_dist, width = 10, height = 24
)









#------------------------------------------------------------------------------
# Visualizing data distributions ----
#------------------------------------------------------------------------------

ssc_itet_paths <- list.files(
    glue("{wkdir}/data/interim/cell_neighborhoods/"),
    recursive = TRUE,
    pattern = "SSC_01 IT-ET Glut"
    )
length(ssc_itet_paths)

future::plan("multisession", workers = 8)
ssc_itet_df <- load_cell_neighborhood_data(ssc_itet_paths, wkdir) %>% 
    filter(slice_id %nin% poor_quality_slices) %>%
    glimpse()

# Visualizing data distributions 
p_bxdist <- ssc_itet_df %>% 
    ggplot(aes(y = cell2_type, x = rel_abund)) +
    geom_point(aes(color = group), alpha = 0.2, position = position_jitterdodge(), size = 0.5) +
    geom_boxplot(aes(fill = group), alpha = 0.2, outlier.shape = NA, color = "black") +
    scale_color_d3(guide = "none") +
    scale_fill_d3() +
    # scale_x_continuous(transform = scales::log_trans()) +
    scale_x_log10() +
    theme_bw() +
    labs(x = "Relative Abundance", y = "Cell Type", fill = "Group")

ggsave(
    glue("{wkdir}/figures/cell_neighborhoods/misc/",
        "SSC_01 IT-ET Glut_rel_abund_dist.png"), 
    p_bxdist, width = 8, height = 6
)

p_deepdive <- ssc_itet_df %>% 
    filter(cell2_type == "34 Immune") %>% 
    ggplot(aes(x = rel_abund, y= slice_id)) +
    geom_point(aes(color = group), alpha = 0.2, size = 0.5, position = position_jitter(height = 0.2, width = 0)) +
    geom_boxplot(alpha = 0.2, outlier.shape = NA, color = "black") +
    scale_x_log10() +
    scale_color_d3() +
    facet_grid(row = vars(run), scales = "free_y", space = "free_y") +
    theme_bw() +
    labs(x = "Relative Abundance", y = "Slice ID") +
    theme(legend.position = "none")

ggsave(
    glue("{wkdir}/figures/cell_neighborhoods/misc/",
        "SSC_01 IT-ET Glut x 34 Immune_rel_abund_dist.png"), 
    p_deepdive, width = 6, height = 4
)

# filter out poor quality slices
models_df <- ssc_itet_df %>%
    nest(data = -c(gating, cell1_type, cell2_type)) %>% 
    mutate(model = purrr::map(data, 
        ~ purrr::possibly(
            ~ nlme::lme(log1p(rel_abund) ~ group, random = ~1|run/slice_id, data = na.omit(.x)),
            otherwise = NULL
        )(.))
    ) %>%
    glimpse()

models_summary_df <- models_df %>% 
  mutate(model_summary = purrr::map(model, ~ broom.mixed::tidy(.x))) %>%
  unnest(model_summary) %>%
  janitor::clean_names() %>%
  dplyr::select(-c(data, model)) %>%
  glimpse()

models_summary_df %>% 
    filter(term == "groupWILDR") %>% 
    View





