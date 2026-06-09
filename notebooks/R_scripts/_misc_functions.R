# Miscellanous functions
#' This script contains function that are used
#' across multiple scripts notebooks 

# ______________________________________________________________________________
#                     Utility Functions
# ______________________________________________________________________________

# inverse of logical statement
`%nin%` <- Negate(`%in%`)

#' Function for executing shell command in R with STDOUT and STEDRR control
#' specifying TRUE for stdout_path redirects the output of the command 
#' into an R object
shell_do <- function(command_string, stdout_path = "", stderr_path = "") {
  inputs <- unlist(stringr::str_split(command_string, " "))
  system2(
    command = inputs[1],
    args = inputs[-1],
    stdout = stdout_path,
    stderr = stderr_path
  )
}


# Function to run any shell command using slurm and future.batchtools
slurm_shell_do <- function(cmd,
                           jobname = glue("slurm-shell-{get_time()}"),
                           working_dir = wkdir,
                           template_path = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
                           memory = "1G",
                           ncpus = 1,
                           walltime = 3600) {
  require(magrittr)
  require(future)
  require(future.batchtools)
  # Initiate future.batchtools backend for parallel processing
  future::plan(
    future.batchtools::batchtools_slurm,
    template = template_path,
    resources = list(
      name = jobname,
      memory = memory,
      ncpus = ncpus,
      walltime = walltime
    )
  )
  job %<-% shell_do(cmd)
}


get_time <- function(){
  print(format(Sys.time(), "%Y-%m-%d_%H:%M:%S"))
}

# function to count ongoing slurm jobs
count_slurm_jobs <- function(params = c("-u", "jboktor")) {
  queue <- system2(
    command = "squeue",
    args = params,
    stdout = TRUE
  )
  length(queue) - 1
}

#' sleep until a condition is true
wait_until <- function(conditional, interval = 2) {
  # Keep looping until the condition is met
  while (!conditional()) {
    Sys.sleep(interval)
  }
}

check_slurm_overload <- function(njobs = 9999) {
  wait_until(function() {
    count_slurm_jobs() < njobs
  })
}

#' Function to rapidly extract the number of files matching a search string
#' while integrating a list of sample names.
#' Input:
#' search_pattern   - a string with glue syntax {.}, not yet glued
#'                      to indicate placement of sampleIDs
#' name_list        - list of names
#' nworkders        - number of threads to use for
#'                    parallel processing
#' Output:
#' A named list of samples with the number of files
#' matching the search critera for each sample
search_file_n <- function(search_pattern, name_list, nworkers = 2, ...) {
  plan(multisession, workers = nworkers)
  file_n <- name_list %>%
    purrr::set_names() %>%
    purrr::map(~glue(search_pattern)) %>%
    furrr::future_map(
      ~ system2(
        command = "ls",
        args = c(., "|", "wc", "-l"),
        stdout = TRUE
      ) %>% as.numeric(),
      .progress = TRUE
    )
  return(file_n)
  }

# chunking function from https://stackoverflow.com/a/16275428
chunk_func <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))

has_error_message <- function(stderr_file) {
  lines <- readLines(stderr_file)
  for (line in lines) {
    if (grepl("error|ERR|quota exceeded", line)) {
      return(TRUE)
    }
  }
  return(FALSE)
}

seriate_matrix_rows <- function(mat,
                                seriate_method = "HC_average",
                                dist_method = "euclidean",
                                nthreads = 8) {
  order <- mat %>%
    parallelDist::parDist(method = dist_method, threads = nthreads) %>%
    # stats::dist(method = dist_method) %>%
    seriation::seriate(method = seriate_method) %>%
    seriation::get_order()
  ranked_order <- rownames(mat)[order]
  return(ranked_order)
}


get_unique_gene_counts <- function(seurat_obj) {
    counts <- SeuratObject::GetAssayData(seurat_obj, layer = "counts")
    clean_genes <- strex::str_before_first(rownames(counts), "-GRCm39")
    clean_counts <- rowsum(counts, group = clean_genes)
    # Sanity check
    # sum(counts[grep("^Gnai3", rownames(counts)), ])
    # sum(clean_counts["Gnai3", ])
    return(clean_counts)
}


# ______________________________________________________________________________
#                     Gene Expression Extraction
# ______________________________________________________________________________

#' Extract gene expression for filtered cells with metadata (memory- and speed-optimized)
#'
#' Returns a data frame of metadata for cells matching celltype (and optionally tissue),
#' with gene expression (raw counts or log2cpm) appended. Avoids Seurat subsetting;
#' uses nCount_RNA when available for log2cpm to skip colSums.
#'
#' Uses column-before-row subsetting for dgCMatrix (column-major) and LayerData
#' with features/cells when Assay5 (avoids loading full matrix). Pass pre-fetched
#' counts via \code{counts} when calling repeatedly for many genes to avoid
#' redundant fetches.
#'
#' @param seurat_obj Seurat object
#' @param celltype_column Name of metadata column for cell type (e.g. "supertype_name")
#' @param celltype_value Value to filter (e.g. "0568 VMH Nr5a1 Glut_3")
#' @param gene Gene symbol (e.g. "Cntn6"); handles -GRCm39 suffix
#' @param tissue Optional tissue value to filter (e.g. "HYP"); NULL skips tissue filter
#' @param return_format "log2cpm" or "raw"
#' @param tissue_column Name of tissue column in metadata (default "tissue")
#' @param assay Assay name; default uses DefaultAssay(seurat_obj)
#' @param counts Optional pre-fetched count matrix (e.g. from GetAssayData); avoids
#'   redundant fetch when calling repeatedly for many genes
#' @return Data frame: metadata for filtered cells + column {gene}_{return_format}
get_gene_expression_df <- function(seurat_obj,
                                  celltype_column,
                                  celltype_value,
                                  gene,
                                  tissue = NULL,
                                  return_format = c("log2cpm", "raw"),
                                  tissue_column = "tissue",
                                  assay = NULL,
                                  counts = NULL) {
  return_format <- match.arg(return_format)
  if (is.null(assay)) assay <- Seurat::DefaultAssay(seurat_obj)

  meta <- seurat_obj@meta.data
  cells_to_keep <- meta[[celltype_column]] == celltype_value
  if (!is.null(tissue)) {
    cells_to_keep <- cells_to_keep & (meta[[tissue_column]] == tissue)
  }

  n_cells <- sum(cells_to_keep)
  if (n_cells == 0) {
    warning("No cells match the filter criteria (celltype: ", celltype_value,
            if (!is.null(tissue)) paste0(", tissue: ", tissue) else "", ")")
    meta_out <- meta[integer(0), , drop = FALSE]
    meta_out[[paste0(gene, "_", return_format)]] <- numeric(0)
    return(meta_out)
  }

  cells_idx <- which(cells_to_keep)
  cell_ids <- colnames(seurat_obj)[cells_idx]
  use_layer_subset <- FALSE
  use_ncount <- "nCount_RNA" %in% names(meta)

  # Resolve gene to feature name (try exact match first, then grep for -GRCm39 etc.)
  if (!is.null(counts)) {
    feature_names <- rownames(counts)
  } else if (inherits(seurat_obj[[assay]], "Assay5")) {
    feature_names <- SeuratObject::Features(seurat_obj, assay = assay)
  } else {
    feature_names <- rownames(Seurat::GetAssayData(seurat_obj, assay = assay, layer = "counts"))
  }
  gene_feat_idx <- match(gene, feature_names)
  if (is.na(gene_feat_idx)) {
    gene_feat_idx <- match(paste0(gene, "-GRCm39"), feature_names)
  }
  if (is.na(gene_feat_idx)) {
    gene_feat_name <- grep(paste0("^", gene, "(-|$)"), feature_names, value = TRUE)[1]
  } else {
    gene_feat_name <- feature_names[gene_feat_idx]
  }
  if (is.na(gene_feat_name)) {
    stop("Gene '", gene, "' not found in assay '", assay, "'.")
  }

  # Fetch data: LayerData subset (Assay5) when possible, else single GetAssayData + column-before-row
  if (is.null(counts)) {
    assay_obj <- seurat_obj[[assay]]
    is_assay5 <- inherits(assay_obj, "Assay5")
    can_use_subset <- is_assay5 && (return_format == "raw" || use_ncount)
    if (can_use_subset) {
      tryCatch({
        counts_sub <- SeuratObject::LayerData(seurat_obj, layer = "counts", assay = assay,
                                              features = gene_feat_name, cells = cell_ids)
        gene_counts <- as.vector(counts_sub)
        if (length(gene_counts) == n_cells) {
          use_layer_subset <- TRUE
        }
      }, error = function(e) NULL)
    }
  }

  if (!use_layer_subset) {
    if (is.null(counts)) {
      counts <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = "counts")
    }
    counts_sub <- counts[, cells_idx, drop = FALSE]
    gene_counts <- as.vector(counts_sub[gene_feat_name, ])
  }

  if (return_format == "log2cpm") {
    if (use_ncount) {
      lib_size <- meta$nCount_RNA[cells_to_keep]
    } else {
      lib_size <- Matrix::colSums(counts_sub)
    }
    cpm <- (gene_counts / lib_size) * 1e6
    gene_values <- log2(cpm + 1)
  } else {
    gene_values <- gene_counts
  }

  meta_out <- meta[cells_to_keep, , drop = FALSE]
  meta_out[[paste0(gene, "_", return_format)]] <- gene_values
  meta_out
}


#' Extract gene expression for multiple genes (memory- and speed-optimized)
#'
#' Like \code{get_gene_expression_df} but accepts a vector of gene names and returns
#' a single data frame with metadata plus one column per gene (e.g. \code{Cntn6_log2cpm},
#' \code{Gabrg3_log2cpm}). Fetches the count matrix once for all genes to avoid
#' redundant I/O.
#'
#' @param seurat_obj Seurat object
#' @param celltype_column Name of metadata column for cell type (e.g. "supertype_name")
#' @param celltype_value Value to filter (e.g. "0568 VMH Nr5a1 Glut_3")
#' @param genes Character vector of gene symbols (e.g. c("Cntn6", "Gabrg3")); handles -GRCm39 suffix
#' @param tissue Optional tissue value to filter (e.g. "HYP"); NULL skips tissue filter
#' @param return_format "log2cpm" or "raw"
#' @param tissue_column Name of tissue column in metadata (default "tissue")
#' @param assay Assay name; default uses DefaultAssay(seurat_obj)
#' @return Data frame: metadata for filtered cells + columns {gene}_{return_format} for each gene
get_gene_expression_df_multi <- function(seurat_obj,
                                         celltype_column,
                                         celltype_value,
                                         genes,
                                         tissue = NULL,
                                         return_format = c("log2cpm", "raw"),
                                         tissue_column = "tissue",
                                         assay = NULL) {
  return_format <- match.arg(return_format)
  if (is.null(assay)) assay <- Seurat::DefaultAssay(seurat_obj)

  meta <- seurat_obj@meta.data
  cells_to_keep <- meta[[celltype_column]] == celltype_value
  if (!is.null(tissue)) {
    cells_to_keep <- cells_to_keep & (meta[[tissue_column]] == tissue)
  }

  n_cells <- sum(cells_to_keep)
  if (n_cells == 0) {
    warning("No cells match the filter criteria (celltype: ", celltype_value,
            if (!is.null(tissue)) paste0(", tissue: ", tissue) else "", ")")
    meta_out <- meta[integer(0), , drop = FALSE]
    for (g in genes) meta_out[[paste0(g, "_", return_format)]] <- numeric(0)
    return(meta_out)
  }

  cells_idx <- which(cells_to_keep)
  cell_ids <- colnames(seurat_obj)[cells_idx]
  use_ncount <- "nCount_RNA" %in% names(meta)

  # Resolve all genes to feature names
  assay_obj <- seurat_obj[[assay]]
  is_assay5 <- inherits(assay_obj, "Assay5")
  if (is_assay5) {
    feature_names <- SeuratObject::Features(seurat_obj, assay = assay)
  } else {
    feature_names <- rownames(Seurat::GetAssayData(seurat_obj, assay = assay, layer = "counts"))
  }

  resolve_gene <- function(gene) {
    gene_feat_idx <- match(gene, feature_names)
    if (is.na(gene_feat_idx)) gene_feat_idx <- match(paste0(gene, "-GRCm39"), feature_names)
    if (is.na(gene_feat_idx)) {
      gene_feat_name <- grep(paste0("^", gene, "(-|$)"), feature_names, value = TRUE)[1]
    } else {
      gene_feat_name <- feature_names[gene_feat_idx]
    }
    gene_feat_name
  }

  gene_feat_names <- vapply(genes, resolve_gene, character(1), USE.NAMES = FALSE)
  missing <- is.na(gene_feat_names)
  if (any(missing)) {
    stop("Gene(s) not found in assay '", assay, "': ",
         paste(genes[missing], collapse = ", "))
  }

  # Fetch counts for all genes at once
  can_use_subset <- is_assay5 && (return_format == "raw" || use_ncount)
  counts_sub <- NULL
  if (can_use_subset) {
    tryCatch({
      counts_sub <- SeuratObject::LayerData(seurat_obj, layer = "counts", assay = assay,
                                            features = gene_feat_names, cells = cell_ids)
      if (nrow(counts_sub) != length(genes) || ncol(counts_sub) != n_cells) {
        counts_sub <- NULL
      }
    }, error = function(e) counts_sub <- NULL)
  }

  counts_full <- NULL
  if (is.null(counts_sub)) {
    counts_full <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = "counts")
    counts_sub <- counts_full[gene_feat_names, cells_idx, drop = FALSE]
  }

  # Compute expression values
  if (return_format == "log2cpm") {
    if (use_ncount) {
      lib_size <- meta$nCount_RNA[cells_to_keep]
    } else {
      lib_size <- Matrix::colSums(
        if (!is.null(counts_full)) counts_full[, cells_idx, drop = FALSE]
        else Seurat::GetAssayData(seurat_obj, assay = assay, layer = "counts")[, cells_idx, drop = FALSE]
      )
    }
    cpm <- (as.matrix(counts_sub) / rep(lib_size, each = nrow(counts_sub))) * 1e6
    gene_values <- log2(cpm + 1)
  } else {
    gene_values <- as.matrix(counts_sub)
  }

  meta_out <- meta[cells_to_keep, , drop = FALSE]
  for (i in seq_along(genes)) {
    meta_out[[paste0(genes[i], "_", return_format)]] <- gene_values[i, ]
  }
  meta_out
}


#' Compute KS test statistics for all sample pairs
#'
#' For each row in sample_pairs (with sample1, sample2 columns), computes the
#' two-sample Kolmogorov-Smirnov statistic comparing the expression distribution
#' between the two samples.
#'
#' @param gene_df Data frame with expression values and sample IDs (from get_gene_expression_df)
#' @param sample_pairs Data frame with sample1 and sample2 columns (unique pairs)
#' @param expr_col Name of expression column, e.g. "Cntn6_raw" or "Cntn6_log2cpm". If NULL, uses gene + return_format.
#' @param gene Gene symbol when expr_col is NULL (e.g. "Cntn6")
#' @param return_format "log2cpm" or "raw" when expr_col is NULL (default "log2cpm")
#' @param sample_col Name of sample column in gene_df (default "sample")
#' @return sample_pairs with ks_stat column added
compute_ks_stats <- function(gene_df, sample_pairs, expr_col = NULL, gene = NULL, return_format = c("log2cpm", "raw"), sample_col = "sample") {
  if (is.null(expr_col)) {
    if (is.null(gene)) stop("Provide either expr_col or gene.")
    return_format <- match.arg(return_format)
    expr_col <- paste0(gene, "_", return_format)
  }
  if (!expr_col %in% names(gene_df)) {
    stop("Column '", expr_col, "' not found in gene_df. Available: ",
         paste(intersect(names(gene_df), grep("_(raw|log2cpm)$", names(gene_df), value = TRUE)), collapse = ", "))
  }
  sample_pairs %>%
    mutate(ks_stat = purrr::map2_dbl(
      sample1, sample2,
      ~ stats::ks.test(
        gene_df[[expr_col]][gene_df[[sample_col]] == .x],
        gene_df[[expr_col]][gene_df[[sample_col]] == .y]
      )$statistic
    ))
}


#' Get top UP and DOWN regulated genes from DEG results (DESeq2, edgeR, or MAST)
#'
#' Auto-detects format from column names and returns the top n most significant
#' genes in each direction for volcano plot highlighting.
#'
#' @param deg_df DESeqResults, DGEExact/GLM, or data frame from DESeq2, edgeR, or MAST
#' @param n Number of genes per direction (default 10)
#' @param strip_suffix Remove -GRCm39 from gene names (default TRUE)
#' @param p_threshold Optional; filter to p < threshold before ranking
#' @return List with \code{up} and \code{down} character vectors of gene names
get_top_deg_genes <- function(deg_df, n = 10, strip_suffix = TRUE, p_threshold = NULL) {
  df <- as.data.frame(deg_df)

  if (!"gene" %in% names(df) && "primerid" %in% names(df)) {
    df$gene <- df$primerid
  } else if (!"gene" %in% names(df)) {
    df$gene <- rownames(df)
  }

  lfc_col <- NULL
  for (c in c("log2FoldChange", "logFC", "coef")) {
    if (c %in% names(df)) { lfc_col <- c; break }
  }
  pval_col <- NULL
  for (c in c("padj", "FDR", "p_adjust")) {
    if (c %in% names(df)) { pval_col <- c; break }
  }
  if (is.null(lfc_col) || is.null(pval_col)) {
    stop("Could not detect LFC column (log2FoldChange/logFC/coef) or p-value column (padj/FDR/p_adjust) in deg_df.")
  }

  keep <- !is.na(df[[lfc_col]]) & !is.na(df[[pval_col]])
  if (!is.null(p_threshold)) {
    keep <- keep & (df[[pval_col]] < p_threshold)
  }
  df <- df[keep, , drop = FALSE]
  df <- df[order(df[[pval_col]]), , drop = FALSE]

  up_idx <- df[[lfc_col]] > 0.25
  down_idx <- df[[lfc_col]] < -0.25
  genes_up <- head(df$gene[up_idx], n)
  genes_down <- head(df$gene[down_idx], n)

  if (strip_suffix) {
    genes_up <- gsub("-GRCm39$", "", genes_up)
    genes_down <- gsub("-GRCm39$", "", genes_down)
  }

  list(up = as.character(genes_up), down = as.character(genes_down))
}


#' Volcano plot ggpacket template for DEG results
#'
#' Reusable template (points, thresholds, ggrepel labels, colors, theme) for
#' DESeq2, edgeR, or MAST volcano plots. Requires \code{ggpackets} and \code{ggrepel}.
#'
#' @param label_data Data frame of genes to label (must have gene_label, sig, and LFC column)
#' @param lfc_col Name of LFC column for nudge direction ("log2FoldChange", "logFC", or "coef")
#' @param lfc_threshold LFC threshold for vertical lines and legend (default 0.25)
#' @param label_sig_only If TRUE, only label significant genes (default FALSE)
#' @param ylim_max Optional max for y-axis; NULL skips (default NULL)
#' @return ggpacket to add to ggplot
#' @seealso https://dgkf.github.io/ggpackets/
ggpk_volcano_deg <- function(label_data,
                             lfc_col = "log2FoldChange",
                             lfc_threshold = 0.25,
                             label_sig_only = FALSE,
                             ylim_max = NULL) {
  if (label_sig_only && "sig" %in% names(label_data)) {
    label_data <- label_data[label_data$sig %in% TRUE, , drop = FALSE]
  }
  nudge_x_vals <- if (nrow(label_data) > 0) {
    ifelse(label_data[[lfc_col]] > 0, 0.1, -0.1)
  } else {
    numeric(0)
  }
  sig_packet <- ggpackets::ggpacket() +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "gray50") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = gene_label, color = sig),
      show.legend = FALSE,
      size = 3,
      max.overlaps = Inf,
      force = 50,
      force_pull = 0.1,
      nudge_x = nudge_x_vals,
      nudge_y = 0.1,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey40",
      segment.size = 0.3,
      segment.alpha = 0.8,
      direction = "both",
      segment.curvature = 0.2,
      segment.ncp = 3,
      segment.angle = 15,
      min.segment.length = 0
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "gray70", "TRUE" = "#E41A1C"),
      labels = c("NS", glue::glue("FDR < 0.05, |LFC| > {lfc_threshold}")),
      na.translate = FALSE
    ) +
    ggplot2::labs(
      x = expression(log[2] ~ "fold change"),
      y = expression(-log[10] ~ "(adjusted p-value)"),
      color = NULL
    ) +
    ggplot2::theme(legend.position = "top")
  if (!is.null(ylim_max)) {
    sig_packet <- sig_packet +
      ggplot2::scale_y_continuous(limits = c(0, ylim_max))
  }
  sig_packet
}


# ______________________________________________________________________________
#                     Gene Expression Plotting
# ______________________________________________________________________________

#' Plot gene expression by sample and microbiome (two-panel violin + boxplot)
#'
#' Produces a combined plot: sample-level (left) and microbiome-level (right)
#' with Wilcoxon p-value. Expects output from get_gene_expression_df.
#'
#' @param gene_expr_df Data frame from get_gene_expression_df (must have sample, microbiome, and {gene}_{return_format} columns)
#' @param gene Gene symbol (e.g. "Cntn6")
#' @param return_format "log2cpm" or "raw" (must match column in df)
#' @param tissue Optional; used in plot title
#' @param celltype_value Optional; used in plot title
#' @param microbiome_col Name of microbiome column (default "microbiome")
#' @param sample_col Name of sample column (default "sample")
#' @param fill_colors Named vector for fill (default SPF/WildR)
#' @return Patchwork object (two-panel plot)
plot_gene_expression <- function(gene_expr_df,
                                 gene,
                                 return_format = c("log2cpm", "raw"),
                                 tissue = NULL,
                                 celltype_value = NULL,
                                 microbiome_col = "microbiome",
                                 sample_col = "sample",
                                 fill_colors = c(SPF = "#1D78B4", WildR = "#FF7F0F")) {
  return_format <- match.arg(return_format)
  expr_col_log2 <- paste0(gene, "_log2cpm")
  expr_col_raw <- paste0(gene, "_raw")
  if (expr_col_log2 %in% names(gene_expr_df)) {
    expr_col <- expr_col_log2
    return_format <- "log2cpm"
  } else if (expr_col_raw %in% names(gene_expr_df)) {
    expr_col <- expr_col_raw
    return_format <- "raw"
  } else {
    stop("Neither '", expr_col_log2, "' nor '", expr_col_raw, "' found in gene_expr_df.")
  }

  y_lab <- if (return_format == "log2cpm") {
    bquote(.(gene) ~ log[2] ~ "(CPM + 1)")
  } else {
    paste0(gene, " counts")
  }

  title_parts <- c(tissue, celltype_value)
  title <- if (length(title_parts) > 0 && any(!sapply(title_parts, is.null))) {
    glue::glue("{gene} expression in {paste(Filter(Negate(is.null), title_parts), collapse = ' ')}")
  } else {
    glue::glue("{gene} expression")
  }

  p_sample <- gene_expr_df %>%
    ggplot2::ggplot(ggplot2::aes(
      x = forcats::fct_reorder(.data[[sample_col]], as.numeric(.data[[microbiome_col]])),
      y = .data[[expr_col]]
    )) +
    ggplot2::geom_violin(scale = "width", alpha = 0.3, trim = FALSE, ggplot2::aes(fill = .data[[microbiome_col]])) +
    ggplot2::geom_jitter(width = 0.2, height = 0, alpha = 0.4, size = 0.5, shape = 21) +
    ggplot2::geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.8,
      position = ggplot2::position_dodge(width = 0.9), ggplot2::aes(fill = .data[[microbiome_col]])) +
    ggplot2::labs(x = NULL, y = y_lab) +
    ggplot2::scale_fill_manual(values = fill_colors) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  p_microbiome <- gene_expr_df %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[[microbiome_col]], y = .data[[expr_col]], fill = .data[[microbiome_col]])) +
    ggplot2::geom_violin(scale = "width", alpha = 0.4, trim = FALSE) +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8,
      position = ggplot2::position_dodge(width = 0.9)) +
    ggplot2::labs(x = NULL, y = y_lab) +
    ggplot2::scale_fill_manual(values = fill_colors) +
    ggplot2::theme_bw() +
    ggpubr::stat_compare_means(
      # method = "t.test",
      method = "wilcox.test",
      label = "p.format",
      label.x.npc = 0.2
      # label.y.npc = "top",
      # hjust = 0,
      # vjust = 1
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  p_sample + p_microbiome + patchwork::plot_layout(nrow = 1, widths = c(4, 1)) +
    patchwork::plot_annotation(title = title)
}
