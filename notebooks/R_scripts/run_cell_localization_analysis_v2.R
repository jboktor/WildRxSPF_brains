#!/usr/bin/env Rscript

# Cell-neighborhood analysis V2 for seqFISH.
#
# This is the MapMyCells-based successor to cell-localization-analysis.R.
# It keeps the original 60-micron neighborhood definition while:
#   1. using subclass and supertype MapMyCells labels from the latest Seurat object;
#   2. analyzing both whole slices and manually drawn final_gate regions;
#   3. using bootstrap >= 0.5 as the primary label set and all labels as sensitivity;
#   4. finding radius neighbors once per slice/region rather than rebuilding a
#      dense all-by-all distance matrix for every query cell type; and
#   5. producing descriptive ROI-level summaries only. WildR-vs-SPF models require
#      an ROI-to-animal map and are intentionally not fitted here.
#
# Run with the spatialomics R binary:
#   Rscript run_cell_localization_analysis_v2.R prepare
#   Rscript run_cell_localization_analysis_v2.R run-job <job_id>
#   Rscript run_cell_localization_analysis_v2.R aggregate
#   Rscript run_cell_localization_analysis_v2.R status
#
# For a SLURM array, call `run-job ${SLURM_ARRAY_TASK_ID}` after `prepare`.

suppressPackageStartupMessages({
  library(BiocNeighbors)
  library(data.table)
  library(glue)
  library(SeuratObject)
})

WKDIR <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
ANALYSIS_DATE <- "2026-11-06"
RADIUS_UM <- 60
MICRONS_PER_PIXEL <- 0.106
PRIMARY_BOOTSTRAP_MIN <- 0.5

LATEST_SEURAT <- file.path(
  WKDIR,
  "data/interim/seurat/seqFISH/merged_roi_seurat_prenorm_2026-02-16.rds"
)
GATED_METADATA <- file.path(
  WKDIR,
  "data/interim/gated-seurat-metadata_2024-12-31.rds"
)

ANALYSIS_DIR <- file.path(
  WKDIR,
  "data/interim/cell_neighborhoods_mapmycells",
  ANALYSIS_DATE
)
LOG_DIR <- file.path(WKDIR, ".cluster_runs/cell_neighborhoods_mapmycells")
METADATA_DIR <- file.path(ANALYSIS_DIR, "slice_metadata")
JOB_DIR <- file.path(ANALYSIS_DIR, "jobs")
MANIFEST_PATH <- file.path(
  ANALYSIS_DIR,
  glue("cell_neighborhood_jobs_{ANALYSIS_DATE}.csv")
)
PREP_QC_PATH <- file.path(
  ANALYSIS_DIR,
  glue("cell_neighborhood_prepare_qc_{ANALYSIS_DATE}.csv")
)
SUMMARY_PATH <- file.path(
  ANALYSIS_DIR,
  glue("cell_neighborhood_roi_summary_{ANALYSIS_DATE}.rds")
)
JOB_QC_PATH <- file.path(
  ANALYSIS_DIR,
  glue("cell_neighborhood_job_qc_{ANALYSIS_DATE}.csv")
)
SUMMARY_TABLE_PATH <- file.path(
  WKDIR,
  "tables/CellNeighborhoods",
  glue("CellNeighborhoods_MapMyCells_roi_summary_{ANALYSIS_DATE}.csv.gz")
)
FIGURE_DIR <- file.path(WKDIR, "figures/CellNeighborhoods")
CONFIDENCE_FIGURE <- file.path(
  FIGURE_DIR,
  "CellNeighborhoods_MapMyCells_confidence_retention.png"
)
WHOLE_SLICE_FIGURE <- file.path(
  FIGURE_DIR,
  "CellNeighborhoods_MapMyCells_whole_slice_subclass_heatmap.png"
)
MANUAL_GATE_FIGURE <- file.path(
  FIGURE_DIR,
  "CellNeighborhoods_MapMyCells_manual_gate_supertype_self_neighborhood.png"
)

LEVELS <- c("subclass", "supertype")
LABEL_SETS <- c("primary_confident", "all_labels_sensitivity")

log_message <- function(...) {
  message(glue("[{format(Sys.time(), '%Y-%m-%d %H:%M:%S')}] "), ...)
}

stop_if_exists <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing)) {
    stop(
      "Refusing to overwrite existing artifact(s):\n",
      paste(existing, collapse = "\n"),
      call. = FALSE
    )
  }
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "-", x)
}

require_columns <- function(x, columns, object_name) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    stop(
      object_name,
      " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

prepare_inputs <- function() {
  dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(MANIFEST_PATH) && file.exists(PREP_QC_PATH)) {
    log_message("Preparation artifacts already exist; leaving them unchanged.")
    log_message("Manifest: ", MANIFEST_PATH)
    return(invisible(fread(MANIFEST_PATH)))
  }
  if (file.exists(MANIFEST_PATH) || file.exists(PREP_QC_PATH)) {
    stop(
      "Preparation is partial. Refusing to overwrite it. Inspect ",
      ANALYSIS_DIR,
      call. = FALSE
    )
  }

  dir.create(METADATA_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(JOB_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(SUMMARY_TABLE_PATH), recursive = TRUE, showWarnings = FALSE)

  log_message("Loading latest seqFISH Seurat object: ", LATEST_SEURAT)
  seurat_obj <- readRDS(LATEST_SEURAT)
  spatial_meta <- as.data.table(seurat_obj@meta.data)
  rm(seurat_obj)
  invisible(gc())

  mapmycells_columns <- c(
    "cell_id", "center_x", "center_y", "roi", "run", "block", "group", "bregma",
    "class_label", "class_name", "class_bootstrapping_probability",
    "subclass_label", "subclass_name", "subclass_bootstrapping_probability",
    "supertype_label", "supertype_name", "supertype_bootstrapping_probability"
  )
  require_columns(spatial_meta, mapmycells_columns, "Latest seqFISH metadata")
  spatial_meta <- spatial_meta[, ..mapmycells_columns]

  if (anyDuplicated(spatial_meta$cell_id)) {
    stop("Latest seqFISH metadata has duplicated cell_id values.", call. = FALSE)
  }

  spatial_meta[, slice_id := paste(group, run, roi)]
  spatial_meta[, `:=`(
    x_microns = center_x * MICRONS_PER_PIXEL,
    y_microns = center_y * MICRONS_PER_PIXEL
  )]

  log_message("Joining manually drawn final_gate calls: ", GATED_METADATA)
  gate_meta <- as.data.table(readRDS(GATED_METADATA))
  require_columns(
    gate_meta,
    c("cell_id", "center_x", "center_y", "run", "group", "bregma", "final_gate"),
    "Gated metadata"
  )
  if (anyDuplicated(gate_meta$cell_id)) {
    stop("Gated metadata has duplicated cell_id values.", call. = FALSE)
  }

  match_index <- match(gate_meta$cell_id, spatial_meta$cell_id)
  if (anyNA(match_index)) {
    stop("Some gated-metadata cells are absent from the latest Seurat object.", call. = FALSE)
  }
  for (field in c("center_x", "center_y", "run", "group", "bregma")) {
    old_value <- gate_meta[[field]]
    new_value <- spatial_meta[[field]][match_index]
    equal <- if (is.numeric(old_value)) {
      isTRUE(all.equal(old_value, new_value, check.attributes = FALSE))
    } else {
      identical(as.character(old_value), as.character(new_value))
    }
    if (!equal) {
      stop("Join validation failed for field: ", field, call. = FALSE)
    }
  }

  gate_join <- gate_meta[, .(
    cell_id,
    final_gate = trimws(as.character(final_gate))
  )]
  gate_join[final_gate == "", final_gate := NA_character_]
  spatial_meta[gate_join, final_gate := i.final_gate, on = "cell_id"]
  spatial_meta[, manual_gate_source_match := cell_id %chin% gate_join$cell_id]

  slice_ids <- sort(unique(spatial_meta$slice_id))
  slice_paths <- file.path(
    METADATA_DIR,
    glue("seqfish_mapmycells_{safe_name(slice_ids)}_{ANALYSIS_DATE}.rds")
  )
  slice_files_exist <- file.exists(slice_paths)
  if (all(slice_files_exist)) {
    log_message(
      "All ", length(slice_ids),
      " per-slice metadata files already exist; reusing them unchanged."
    )
  } else if (any(slice_files_exist)) {
    stop(
      "Only some per-slice metadata files exist. Refusing to overwrite or ",
      "continue from an ambiguous partial preparation in ",
      METADATA_DIR,
      call. = FALSE
    )
  } else {
    log_message("Writing ", length(slice_ids), " per-slice metadata files.")
    for (i in seq_along(slice_ids)) {
      saveRDS(
        spatial_meta[slice_id == slice_ids[[i]]],
        slice_paths[[i]],
        compress = FALSE
      )
    }
  }

  slice_lookup <- data.table(slice_id = slice_ids, metadata_path = slice_paths)
  whole_slice_jobs <- spatial_meta[
    ,
    .(n_cells = .N),
    by = .(slice_id, group, run, roi, bregma)
  ][
    slice_lookup,
    metadata_path := i.metadata_path,
    on = "slice_id"
  ][
    ,
    `:=`(scope = "whole_slice", region = "WholeSlice")
  ]

  manual_gate_jobs <- spatial_meta[
    !is.na(final_gate),
    .(n_cells = .N),
    by = .(slice_id, group, run, roi, bregma, region = final_gate)
  ][
    slice_lookup,
    metadata_path := i.metadata_path,
    on = "slice_id"
  ][
    ,
    scope := "manual_gate"
  ]

  manifest <- rbindlist(
    list(whole_slice_jobs, manual_gate_jobs),
    use.names = TRUE,
    fill = TRUE
  )
  setorder(manifest, scope, slice_id, region)
  manifest[, job_id := .I]
  manifest[, output_path := file.path(
    JOB_DIR,
    paste0(
      "cell_neighborhood_", sprintf("%03d", job_id), "_",
      safe_name(slice_id), "_", safe_name(region), "_",
      ANALYSIS_DATE, ".rds"
    )
  )]
  setcolorder(
    manifest,
    c(
      "job_id", "scope", "region", "slice_id", "group", "run", "roi",
      "bregma", "n_cells", "metadata_path", "output_path"
    )
  )

  prep_qc <- rbindlist(lapply(LEVELS, function(level) {
    probability <- spatial_meta[[paste0(level, "_bootstrapping_probability")]]
    manual_gate_qc <- spatial_meta[
      !is.na(final_gate),
      .(
        n_cells = .N,
        n_primary_confident = sum(
          get(paste0(level, "_bootstrapping_probability")) >=
            PRIMARY_BOOTSTRAP_MIN
        ),
        fraction_primary_confident = mean(
          get(paste0(level, "_bootstrapping_probability")) >=
            PRIMARY_BOOTSTRAP_MIN
        )
      ),
      by = .(region = final_gate)
    ]
    manual_gate_qc[, `:=`(scope = "manual_gate", level = level)]
    rbindlist(list(
      data.table(
        scope = "whole_slice",
        region = "WholeSlice",
        level = level,
        n_cells = nrow(spatial_meta),
        n_primary_confident = sum(probability >= PRIMARY_BOOTSTRAP_MIN),
        fraction_primary_confident = mean(probability >= PRIMARY_BOOTSTRAP_MIN)
      ),
      manual_gate_qc
    ), use.names = TRUE)
  }))
  setcolorder(
    prep_qc,
    c(
      "scope", "region", "level", "n_cells", "n_primary_confident",
      "fraction_primary_confident"
    )
  )

  fwrite(manifest, MANIFEST_PATH)
  fwrite(prep_qc, PREP_QC_PATH)
  log_message("Prepared ", nrow(manifest), " jobs.")
  log_message("Manifest: ", MANIFEST_PATH)
  invisible(manifest)
}

summarize_level <- function(cells, edges, level, label_set) {
  label_column <- paste0(level, "_label")
  name_column <- paste0(level, "_name")
  probability_column <- paste0(level, "_bootstrapping_probability")

  labels <- as.character(cells[[label_column]])
  names_by_cell <- as.character(cells[[name_column]])
  probabilities <- cells[[probability_column]]

  eligible <- !is.na(labels) & nzchar(labels)
  if (label_set == "primary_confident") {
    eligible <- eligible &
      !is.na(probabilities) &
      probabilities >= PRIMARY_BOOTSTRAP_MIN
  }

  level_edges <- edges[eligible[from] & eligible[to]]
  eligible_cells <- which(eligible)
  query_counts <- data.table(
    query_label = labels[eligible_cells],
    query_name = names_by_cell[eligible_cells]
  )[
    ,
    .(n_query_total = .N),
    by = .(query_label, query_name)
  ]
  neighbor_counts <- data.table(
    neighbor_label = labels[eligible_cells],
    neighbor_name = names_by_cell[eligible_cells]
  )[
    ,
    .(n_neighbor_region = .N),
    by = .(neighbor_label, neighbor_name)
  ]

  if (!nrow(level_edges)) {
    return(list(
      cell_pairs = data.table(),
      roi_summary = data.table(),
      qc = data.table(
        level = level,
        label_set = label_set,
        n_eligible_cells = length(eligible_cells),
        n_directed_edges = 0L,
        n_cell_type_pairs = 0L
      )
    ))
  }

  level_edges[, `:=`(
    query_label = labels[from],
    query_name = names_by_cell[from],
    neighbor_label = labels[to],
    neighbor_name = names_by_cell[to]
  )]

  cell_pairs <- level_edges[
    ,
    .(n = .N),
    by = .(
      cell_index = from,
      query_label,
      query_name,
      neighbor_label,
      neighbor_name
    )
  ]
  neighbor_totals <- level_edges[
    ,
    .(n_total_neighbors = .N),
    by = .(cell_index = from)
  ]
  cell_pairs[
    neighbor_totals,
    n_total_neighbors := i.n_total_neighbors,
    on = "cell_index"
  ]
  cell_pairs[, rel_abund := n / n_total_neighbors]
  cell_pairs[, cell_id := cells$cell_id[cell_index]]

  query_with_neighbors <- unique(
    cell_pairs[, .(cell_index, query_label, query_name)]
  )[
    ,
    .(n_query_with_neighbors = .N),
    by = .(query_label, query_name)
  ]

  roi_summary <- cell_pairs[
    ,
    .(
      n_query_with_neighbor_type = uniqueN(cell_index),
      total_neighbor_occurrences = sum(n),
      sum_rel_abund = sum(rel_abund)
    ),
    by = .(query_label, query_name, neighbor_label, neighbor_name)
  ]
  roi_summary[
    query_counts,
    n_query_total := i.n_query_total,
    on = .(query_label, query_name)
  ]
  roi_summary[
    query_with_neighbors,
    n_query_with_neighbors := i.n_query_with_neighbors,
    on = .(query_label, query_name)
  ]
  roi_summary[
    neighbor_counts,
    n_neighbor_region := i.n_neighbor_region,
    on = .(neighbor_label, neighbor_name)
  ]
  roi_summary[, `:=`(
    prevalence = n_query_with_neighbor_type / n_query_with_neighbors,
    mean_neighbor_count = total_neighbor_occurrences / n_query_with_neighbors,
    mean_rel_abund = sum_rel_abund / n_query_with_neighbors,
    mean_region_normalized_count =
      total_neighbor_occurrences /
      n_query_with_neighbors /
      n_neighbor_region
  )]

  cell_pairs[, `:=`(level = level, label_set = label_set)]
  setcolorder(
    cell_pairs,
    c(
      "level", "label_set", "cell_id", "query_label", "query_name",
      "neighbor_label", "neighbor_name", "n", "n_total_neighbors", "rel_abund"
    )
  )
  cell_pairs[, cell_index := NULL]
  roi_summary[, `:=`(level = level, label_set = label_set)]

  list(
    cell_pairs = cell_pairs,
    roi_summary = roi_summary,
    qc = data.table(
      level = level,
      label_set = label_set,
      n_eligible_cells = length(eligible_cells),
      n_directed_edges = nrow(level_edges),
      n_cell_type_pairs = nrow(roi_summary)
    )
  )
}

calculate_cell_neighborhoods_v2 <- function(job) {
  start_time <- Sys.time()
  log_message(
    "Starting job ", job$job_id, ": ", job$slice_id, " / ",
    job$scope, " / ", job$region
  )

  cells <- as.data.table(readRDS(job$metadata_path))
  if (job$scope == "manual_gate") {
    cells <- cells[final_gate == job$region]
  } else if (job$scope != "whole_slice") {
    stop("Unknown scope: ", job$scope, call. = FALSE)
  }
  if (nrow(cells) != job$n_cells) {
    stop(
      "Manifest expected ", job$n_cells, " cells but found ", nrow(cells),
      call. = FALSE
    )
  }

  coordinates <- as.matrix(cells[, .(x_microns, y_microns)])
  neighbor_result <- BiocNeighbors::findNeighbors(
    coordinates,
    threshold = RADIUS_UM,
    get.distance = FALSE
  )
  edges <- data.table(
    from = rep.int(seq_along(neighbor_result$index), lengths(neighbor_result$index)),
    to = unlist(neighbor_result$index, use.names = FALSE)
  )
  edges <- edges[from != to]

  results <- lapply(LEVELS, function(level) {
    lapply(LABEL_SETS, function(label_set) {
      summarize_level(cells, edges, level, label_set)
    })
  })
  results <- unlist(results, recursive = FALSE)

  cell_pairs <- rbindlist(
    lapply(results, `[[`, "cell_pairs"),
    use.names = TRUE,
    fill = TRUE
  )
  roi_summary <- rbindlist(
    lapply(results, `[[`, "roi_summary"),
    use.names = TRUE,
    fill = TRUE
  )
  qc <- rbindlist(lapply(results, `[[`, "qc"), use.names = TRUE, fill = TRUE)

  job_columns <- c(
    "job_id", "scope", "region", "slice_id", "group", "run", "roi", "bregma"
  )
  for (column in job_columns) {
    roi_summary[, (column) := job[[column]]]
    qc[, (column) := job[[column]]]
  }
  setcolorder(roi_summary, c(job_columns, setdiff(names(roi_summary), job_columns)))
  setcolorder(qc, c(job_columns, setdiff(names(qc), job_columns)))

  elapsed_minutes <- as.numeric(
    difftime(Sys.time(), start_time, units = "mins")
  )
  qc[, `:=`(
    n_region_cells = nrow(cells),
    n_spatial_edges = nrow(edges),
    radius_um = RADIUS_UM,
    elapsed_minutes = elapsed_minutes
  )]

  list(
    parameters = list(
      radius_um = RADIUS_UM,
      microns_per_pixel = MICRONS_PER_PIXEL,
      primary_bootstrap_min = PRIMARY_BOOTSTRAP_MIN,
      levels = LEVELS,
      label_sets = LABEL_SETS
    ),
    job = as.list(job),
    cell_pairs = cell_pairs,
    roi_summary = roi_summary,
    qc = qc
  )
}

run_job <- function(job_id) {
  if (!file.exists(MANIFEST_PATH)) {
    stop("Run `prepare` first. Missing manifest: ", MANIFEST_PATH, call. = FALSE)
  }
  manifest <- fread(MANIFEST_PATH)
  requested_job_id <- as.integer(job_id)
  if (is.na(requested_job_id)) {
    stop("job_id must be numeric: ", job_id, call. = FALSE)
  }
  job <- manifest[job_id == requested_job_id]
  if (nrow(job) != 1L) {
    stop("Unknown or duplicated job_id: ", job_id, call. = FALSE)
  }
  if (file.exists(job$output_path)) {
    log_message("Output already exists; leaving it unchanged: ", job$output_path)
    return(invisible(job$output_path))
  }

  result <- calculate_cell_neighborhoods_v2(job)
  saveRDS(result, job$output_path)
  log_message("Saved: ", job$output_path)
  invisible(job$output_path)
}

aggregate_jobs <- function() {
  if (!file.exists(MANIFEST_PATH)) {
    stop("Run `prepare` first. Missing manifest: ", MANIFEST_PATH, call. = FALSE)
  }
  manifest <- fread(MANIFEST_PATH)
  missing <- manifest[!file.exists(output_path)]
  if (nrow(missing)) {
    stop(
      nrow(missing),
      " job output(s) are missing. First missing IDs: ",
      paste(head(missing$job_id, 20), collapse = ", "),
      call. = FALSE
    )
  }
  stop_if_exists(c(SUMMARY_PATH, JOB_QC_PATH, SUMMARY_TABLE_PATH))

  log_message("Aggregating ", nrow(manifest), " completed jobs.")
  roi_summaries <- vector("list", nrow(manifest))
  job_qcs <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    result <- readRDS(manifest$output_path[[i]])
    roi_summaries[[i]] <- result$roi_summary
    job_qcs[[i]] <- result$qc
    rm(result)
  }
  roi_summary <- rbindlist(roi_summaries, use.names = TRUE, fill = TRUE)
  job_qc <- rbindlist(job_qcs, use.names = TRUE, fill = TRUE)
  rm(roi_summaries, job_qcs)
  invisible(gc())

  saveRDS(roi_summary, SUMMARY_PATH, compress = FALSE)
  fwrite(job_qc, JOB_QC_PATH)
  fwrite(roi_summary, SUMMARY_TABLE_PATH, compress = "gzip")
  log_message("ROI summary: ", SUMMARY_PATH)
  log_message("Table: ", SUMMARY_TABLE_PATH)
  log_message(
    "No WildR-vs-SPF model was fitted; an ROI-to-animal map is required."
  )
  invisible(roi_summary)
}

plot_confidence_retention <- function() {
  if (!file.exists(PREP_QC_PATH)) {
    stop("Missing preparation QC table: ", PREP_QC_PATH, call. = FALSE)
  }
  stop_if_exists(CONFIDENCE_FIGURE)
  suppressPackageStartupMessages({
    library(ggplot2)
    library(scales)
  })
  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

  qc <- fread(PREP_QC_PATH)
  qc[, region_display := fifelse(scope == "whole_slice", "Whole slice", region)]
  region_order <- qc[
    ,
    .(retained = mean(fraction_primary_confident)),
    by = region_display
  ][order(retained), region_display]
  qc[, region_display := factor(region_display, levels = region_order)]

  plot <- ggplot(
    qc,
    aes(
      x = region_display,
      y = fraction_primary_confident,
      color = level,
      shape = level
    )
  ) +
    geom_point(
      position = position_dodge(width = 0.45),
      size = 3
    ) +
    coord_flip() +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    labs(
      x = NULL,
      y = "Cells retained at bootstrap probability >= 0.5",
      color = "MapMyCells level",
      shape = "MapMyCells level",
      title = "MapMyCells confidence retention by spatial scope"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")

  ggsave(CONFIDENCE_FIGURE, plot, width = 9, height = 6.5, dpi = 200)
  log_message("Saved: ", CONFIDENCE_FIGURE)
  invisible(CONFIDENCE_FIGURE)
}

read_roi_summary <- function() {
  if (!file.exists(SUMMARY_PATH)) {
    stop("Missing aggregated ROI summary: ", SUMMARY_PATH, call. = FALSE)
  }
  as.data.table(readRDS(SUMMARY_PATH))
}

plot_whole_slice_subclass_heatmap <- function() {
  stop_if_exists(WHOLE_SLICE_FIGURE)
  suppressPackageStartupMessages(library(ggplot2))
  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

  summary <- read_roi_summary()[
    scope == "whole_slice" &
      level == "subclass" &
      label_set == "primary_confident"
  ]
  query_sizes <- unique(
    summary[, .(slice_id, query_label, query_name, n_query_total)]
  )[
    ,
    .(n_cells = sum(n_query_total)),
    by = .(query_label, query_name)
  ][order(-n_cells)]
  top_labels <- head(query_sizes$query_label, 15)
  label_names <- query_sizes[
    query_label %chin% top_labels,
    setNames(query_name, query_label)
  ]

  heatmap <- summary[
    query_label %chin% top_labels & neighbor_label %chin% top_labels,
    .(
      mean_rel_abund = weighted.mean(
        mean_rel_abund,
        w = n_query_with_neighbors,
        na.rm = TRUE
      )
    ),
    by = .(query_label, neighbor_label)
  ]
  grid <- CJ(query_label = top_labels, neighbor_label = top_labels)
  heatmap <- heatmap[grid, on = .(query_label, neighbor_label)]
  heatmap[is.na(mean_rel_abund), mean_rel_abund := 0]
  heatmap[, `:=`(
    query_name = label_names[query_label],
    neighbor_name = label_names[neighbor_label]
  )]
  ordered_names <- query_sizes[
    query_label %chin% top_labels,
    query_name
  ]
  heatmap[, `:=`(
    query_name = factor(query_name, levels = rev(ordered_names)),
    neighbor_name = factor(neighbor_name, levels = ordered_names)
  )]

  plot <- ggplot(
    heatmap,
    aes(x = neighbor_name, y = query_name, fill = mean_rel_abund)
  ) +
    geom_tile() +
    scale_fill_viridis_c(trans = "sqrt") +
    labs(
      x = "Neighbor subclass",
      y = "Query subclass",
      fill = "Mean local\nrelative abundance",
      title = "Whole-slice neighborhoods of the 15 most abundant subclasses",
      subtitle = "Primary MapMyCells labels (bootstrap probability >= 0.5)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  ggsave(WHOLE_SLICE_FIGURE, plot, width = 12, height = 10, dpi = 200)
  log_message("Saved: ", WHOLE_SLICE_FIGURE)
  invisible(WHOLE_SLICE_FIGURE)
}

plot_manual_gate_supertype_self_neighborhood <- function() {
  stop_if_exists(MANUAL_GATE_FIGURE)
  suppressPackageStartupMessages(library(ggplot2))
  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

  summary <- read_roi_summary()[
    scope == "manual_gate" &
      level == "supertype" &
      label_set == "primary_confident" &
      query_label == neighbor_label
  ]
  query_sizes <- unique(
    summary[, .(slice_id, region, query_label, query_name, n_query_total)]
  )[
    ,
    .(n_cells = sum(n_query_total)),
    by = .(query_label, query_name)
  ][order(-n_cells)]
  top_labels <- head(query_sizes$query_label, 20)
  ordered_names <- query_sizes[
    query_label %chin% top_labels,
    query_name
  ]

  heatmap <- summary[
    query_label %chin% top_labels,
    .(
      mean_self_rel_abund = weighted.mean(
        mean_rel_abund,
        w = n_query_with_neighbors,
        na.rm = TRUE
      )
    ),
    by = .(region, query_label, query_name)
  ]
  heatmap[, query_name := factor(query_name, levels = rev(ordered_names))]

  plot <- ggplot(
    heatmap,
    aes(x = region, y = query_name, fill = mean_self_rel_abund)
  ) +
    geom_tile() +
    scale_fill_viridis_c(trans = "sqrt", na.value = "white") +
    labs(
      x = "Manually drawn anatomical region",
      y = "Supertype",
      fill = "Mean same-type\nrelative abundance",
      title = "Local self-neighborhoods across manually gated regions",
      subtitle = "Twenty most abundant supertypes; primary confident labels"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  ggsave(MANUAL_GATE_FIGURE, plot, width = 12, height = 10, dpi = 200)
  log_message("Saved: ", MANUAL_GATE_FIGURE)
  invisible(MANUAL_GATE_FIGURE)
}

report_status <- function() {
  cat("Analysis directory:", ANALYSIS_DIR, "\n")
  cat("Manifest:", if (file.exists(MANIFEST_PATH)) "present" else "missing", "\n")
  if (!file.exists(MANIFEST_PATH)) {
    return(invisible(NULL))
  }
  manifest <- fread(MANIFEST_PATH)
  complete <- file.exists(manifest$output_path)
  cat("Jobs:", sum(complete), "/", nrow(manifest), "complete\n")
  if (any(!complete)) {
    cat("Pending job IDs:", paste(manifest$job_id[!complete], collapse = ","), "\n")
  }
  cat("Aggregate:", if (file.exists(SUMMARY_PATH)) "present" else "missing", "\n")
  invisible(manifest)
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  run_cell_localization_analysis_v2.R prepare",
      "  run_cell_localization_analysis_v2.R run-job <job_id>",
      "  run_cell_localization_analysis_v2.R aggregate",
      "  run_cell_localization_analysis_v2.R plot-confidence",
      "  run_cell_localization_analysis_v2.R plot-whole-slice",
      "  run_cell_localization_analysis_v2.R plot-manual-gates",
      "  run_cell_localization_analysis_v2.R status",
      sep = "\n"
    ),
    "\n"
  )
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || args[[1]] %in% c("-h", "--help", "help")) {
    usage()
    return(invisible(NULL))
  }

  action <- args[[1]]
  if (action == "prepare") {
    prepare_inputs()
  } else if (action == "run-job") {
    if (length(args) < 2L) {
      stop("run-job requires a numeric job_id.", call. = FALSE)
    }
    run_job(args[[2]])
  } else if (action == "aggregate") {
    aggregate_jobs()
  } else if (action == "plot-confidence") {
    plot_confidence_retention()
  } else if (action == "plot-whole-slice") {
    plot_whole_slice_subclass_heatmap()
  } else if (action == "plot-manual-gates") {
    plot_manual_gate_supertype_self_neighborhood()
  } else if (action == "status") {
    report_status()
  } else {
    usage()
    stop("Unknown action: ", action, call. = FALSE)
  }
}

if (sys.nframe() == 0L) {
  main()
}
