#!/usr/bin/env Rscript

# Exploratory SPF-versus-WildR tests for seqFISH cell neighborhoods.
#
# The unit analyzed here is the ROI section, not the cell. These results are
# explicitly exploratory because the repository does not contain an ROI-to-animal
# mapping. If multiple sections come from one animal, uncertainty is understated.
#
# Primary outcome: mean local relative abundance within 60 microns.
# Model: asin(sqrt(outcome)) ~ imaging run + microbiome group, fitted with
# limma empirical-Bayes moderation. Raw-scale run-adjusted differences and CIs
# are reported for interpretation. A run-stratified exact permutation test
# evaluates each family of neighborhood pairs jointly.
#
# Usage:
#   <spatialomics>/bin/Rscript run_cell_neighborhood_group_tests_v2.R self-test
#   <spatialomics>/bin/Rscript run_cell_neighborhood_group_tests_v2.R run

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(glue)
  library(limma)
  library(scales)
})

WKDIR <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
ANALYSIS_DATE <- "2026-11-06"
RUN_VERSION <- "r3"
ANALYSIS_DIR <- file.path(
  WKDIR,
  "data/interim/cell_neighborhoods_mapmycells",
  ANALYSIS_DATE
)
INPUT_CSV <- file.path(
  WKDIR,
  "tables/CellNeighborhoods",
  glue("CellNeighborhoods_MapMyCells_roi_summary_{ANALYSIS_DATE}.csv.gz")
)
MANIFEST_CSV <- file.path(
  ANALYSIS_DIR,
  glue("cell_neighborhood_jobs_{ANALYSIS_DATE}.csv")
)

RESULT_RDS <- file.path(
  ANALYSIS_DIR,
  glue(
    "cell_neighborhood_spf_wildr_roi_exploratory_",
    "{RUN_VERSION}_{ANALYSIS_DATE}.rds"
  )
)
PAIRWISE_CSV <- file.path(
  WKDIR,
  "tables/CellNeighborhoods",
  glue(
    "CellNeighborhoods_SPF-WildR_pairwise_roi_exploratory_",
    "{RUN_VERSION}_{ANALYSIS_DATE}.csv.gz"
  )
)
GLOBAL_CSV <- file.path(
  WKDIR,
  "tables/CellNeighborhoods",
  glue(
    "CellNeighborhoods_SPF-WildR_global_roi_exploratory_",
    "{RUN_VERSION}_{ANALYSIS_DATE}.csv"
  )
)
SENSITIVITY_CSV <- file.path(
  WKDIR,
  "tables/CellNeighborhoods",
  glue(
    "CellNeighborhoods_SPF-WildR_confidence_sensitivity_",
    "{RUN_VERSION}_{ANALYSIS_DATE}.csv"
  )
)
ELIGIBILITY_CSV <- file.path(
  WKDIR,
  "tables/CellNeighborhoods",
  glue(
    "CellNeighborhoods_SPF-WildR_eligibility_",
    "{RUN_VERSION}_{ANALYSIS_DATE}.csv"
  )
)
TOP_HITS_CSV <- file.path(
  WKDIR,
  "tables/CellNeighborhoods",
  glue(
    "CellNeighborhoods_SPF-WildR_top_pairs_",
    "{RUN_VERSION}_{ANALYSIS_DATE}.csv"
  )
)

FIGURE_DIR <- file.path(WKDIR, "figures/CellNeighborhoods")
GLOBAL_FIGURE <- file.path(
  FIGURE_DIR,
  glue("CellNeighborhoods_SPF-WildR_global_permutation_{RUN_VERSION}.png")
)
EFFECT_FIGURE <- file.path(
  FIGURE_DIR,
  glue("CellNeighborhoods_SPF-WildR_top_effects_{RUN_VERSION}.png")
)
SENSITIVITY_FIGURE <- file.path(
  FIGURE_DIR,
  glue("CellNeighborhoods_SPF-WildR_confidence_sensitivity_{RUN_VERSION}.png")
)

MIN_QUERY_CELLS <- 20L
MIN_NEIGHBOR_CELLS <- 20L
MIN_GROUP_ROIS <- 3L
MIN_RUNS <- 3L
GLOBAL_T2_CAP <- 100
POOR_QUALITY_SLICES <- c(
  "spf run3 roi4",
  "wildr run1 roi3",
  "spf run2 roi3",
  "wildr run3 roi2"
)

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

required_columns <- function(x, columns, object_name) {
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

family_label <- function(scope, region, level, label_set) {
  paste(scope, region, level, label_set, sep = "__")
}

build_sample_metadata <- function(manifest) {
  samples <- unique(
    manifest[
      ,
      .(scope, region, slice_id, group, run, roi, bregma)
    ]
  )
  samples[, `:=`(
    group = factor(tolower(group), levels = c("spf", "wildr")),
    run = factor(run, levels = paste0("run", 1:4)),
    slice_id = tolower(slice_id)
  )]
  samples <- samples[!slice_id %chin% POOR_QUALITY_SLICES]
  setorder(samples, scope, region, run, group, roi)
  samples[
    ,
    sample_index := seq_len(.N),
    by = .(scope, region)
  ]
  samples
}

build_family_matrix <- function(family_rows, sample_meta) {
  pair_columns <- c(
    "query_label", "query_name", "neighbor_label", "neighbor_name"
  )

  pairs <- unique(family_rows[, ..pair_columns])
  pairs[, feature_id := .I]

  query_counts <- unique(
    family_rows[
      ,
      .(
        slice_id,
        query_label,
        query_name,
        n_query_total,
        n_query_with_neighbors
      )
    ]
  )
  neighbor_counts <- unique(
    family_rows[
      ,
      .(
        slice_id,
        neighbor_label,
        neighbor_name,
        n_neighbor_region
      )
    ]
  )

  observed <- pairs[
    family_rows,
    on = pair_columns,
    nomatch = 0L
  ][
    ,
    .(feature_id, slice_id, mean_rel_abund)
  ]

  sample_indices <- sample_meta$sample_index
  grid <- pairs[
    ,
    .(sample_index = sample_indices),
    by = pair_columns
  ]
  grid[
    pairs,
    feature_id := i.feature_id,
    on = pair_columns
  ]
  grid[
    sample_meta,
    `:=`(
      slice_id = i.slice_id,
      group = i.group,
      run = i.run
    ),
    on = "sample_index"
  ]
  grid[
    query_counts,
    `:=`(
      n_query_total = i.n_query_total,
      n_query_with_neighbors = i.n_query_with_neighbors
    ),
    on = .(slice_id, query_label, query_name)
  ]
  grid[
    neighbor_counts,
    n_neighbor_region := i.n_neighbor_region,
    on = .(slice_id, neighbor_label, neighbor_name)
  ]
  grid[
    observed,
    mean_rel_abund := i.mean_rel_abund,
    on = .(feature_id, slice_id)
  ]

  grid[
    ,
    eligible := !is.na(n_query_total) &
      !is.na(n_query_with_neighbors) &
      !is.na(n_neighbor_region) &
      n_query_total >= MIN_QUERY_CELLS &
      n_query_with_neighbors >= MIN_QUERY_CELLS &
      n_neighbor_region >= MIN_NEIGHBOR_CELLS
  ]
  grid[eligible == TRUE & is.na(mean_rel_abund), mean_rel_abund := 0]
  grid[eligible == FALSE, mean_rel_abund := NA_real_]

  available_by_group <- sample_meta[, .N, by = group]
  min_spf <- max(
    MIN_GROUP_ROIS,
    ceiling(available_by_group[group == "spf", N] / 2)
  )
  min_wildr <- max(
    MIN_GROUP_ROIS,
    ceiling(available_by_group[group == "wildr", N] / 2)
  )

  eligibility <- grid[
    eligible == TRUE,
    .(
      n_spf = sum(group == "spf"),
      n_wildr = sum(group == "wildr"),
      n_runs = uniqueN(run)
    ),
    by = feature_id
  ]
  eligibility[
    ,
    tested := n_spf >= min_spf &
      n_wildr >= min_wildr &
      n_runs >= MIN_RUNS
  ]
  tested_ids <- eligibility[tested == TRUE, feature_id]

  if (!length(tested_ids)) {
    return(NULL)
  }

  model_grid <- grid[feature_id %in% tested_ids]
  wide <- dcast(
    model_grid,
    feature_id ~ sample_index,
    value.var = "mean_rel_abund"
  )
  setorder(wide, feature_id)
  feature_ids <- wide$feature_id
  raw_matrix <- as.matrix(wide[, -1L])
  storage.mode(raw_matrix) <- "double"
  colnames(raw_matrix) <- as.character(sample_meta$sample_index)
  transformed_matrix <- asin(sqrt(pmin(pmax(raw_matrix, 0), 1)))

  design <- make_design(sample_meta)
  estimable <- vapply(
    seq_len(nrow(raw_matrix)),
    function(i) {
      observed <- !is.na(raw_matrix[i, ])
      sum(observed) > ncol(design) &&
        qr(design[observed, , drop = FALSE])$rank == ncol(design)
    },
    logical(1)
  )
  raw_matrix <- raw_matrix[estimable, , drop = FALSE]
  transformed_matrix <- transformed_matrix[estimable, , drop = FALSE]
  feature_ids <- feature_ids[estimable]
  if (!length(feature_ids)) {
    return(NULL)
  }

  list(
    raw_matrix = raw_matrix,
    transformed_matrix = transformed_matrix,
    feature_meta = pairs[
      data.table(feature_id = feature_ids),
      on = "feature_id"
    ],
    eligibility = eligibility[
      data.table(feature_id = feature_ids),
      on = "feature_id"
    ],
    sample_meta = copy(sample_meta),
    min_spf = min_spf,
    min_wildr = min_wildr
  )
}

make_design <- function(sample_meta, group_vector = NULL) {
  group_numeric <- if (is.null(group_vector)) {
    as.integer(sample_meta$group == "wildr")
  } else {
    as.integer(group_vector)
  }
  run_design <- model.matrix(~ run, data = sample_meta)
  design <- cbind(run_design, group_wildr = group_numeric)
  design
}

fit_pairwise_models <- function(family_data) {
  design <- make_design(family_data$sample_meta)
  if (qr(design)$rank < ncol(design)) {
    stop("Observed family design is rank-deficient.", call. = FALSE)
  }
  coef_index <- which(colnames(design) == "group_wildr")

  robust <- requireNamespace("statmod", quietly = TRUE)
  transformed_fit <- eBayes(
    lmFit(family_data$transformed_matrix, design),
    robust = robust,
    trend = TRUE
  )
  raw_fit <- eBayes(
    lmFit(family_data$raw_matrix, design),
    robust = robust,
    trend = TRUE
  )

  transformed_table <- topTable(
    transformed_fit,
    coef = coef_index,
    number = Inf,
    sort.by = "none",
    adjust.method = "BH"
  )
  transformed_table <- as.data.table(
    transformed_table,
    keep.rownames = "matrix_row"
  )
  transformed_table[, matrix_row := as.integer(matrix_row)]

  raw_se <- raw_fit$stdev.unscaled[, coef_index] * sqrt(raw_fit$s2.post)
  raw_df <- raw_fit$df.total
  raw_critical <- qt(0.975, df = raw_df)
  raw_delta <- raw_fit$coefficients[, coef_index]

  result <- cbind(
    family_data$feature_meta,
    family_data$eligibility[, .(n_spf, n_wildr, n_runs)],
    transformed_table[
      ,
      .(
        estimate_asin_sqrt = logFC,
        average_asin_sqrt = AveExpr,
        moderated_t = t,
        p_value = P.Value,
        fdr_family = adj.P.Val,
        log_odds = B
      )
    ],
    data.table(
      delta_raw_adjusted = raw_delta,
      delta_raw_ci_low = raw_delta - raw_critical * raw_se,
      delta_raw_ci_high = raw_delta + raw_critical * raw_se
    )
  )

  groups <- family_data$sample_meta$group
  result[, `:=`(
    mean_spf_raw = rowMeans(
      family_data$raw_matrix[, groups == "spf", drop = FALSE],
      na.rm = TRUE
    ),
    mean_wildr_raw = rowMeans(
      family_data$raw_matrix[, groups == "wildr", drop = FALSE],
      na.rm = TRUE
    )
  )]
  result[]
}

enumerate_run_stratified_groups <- function(sample_meta) {
  run_indices <- split(seq_len(nrow(sample_meta)), sample_meta$run, drop = TRUE)
  run_options <- lapply(run_indices, function(indices) {
    n_wildr <- sum(sample_meta$group[indices] == "wildr")
    combinations <- if (n_wildr == 0L) {
      matrix(integer(), nrow = 0L, ncol = 1L)
    } else if (n_wildr == length(indices)) {
      matrix(indices, ncol = 1L)
    } else {
      combn(indices, n_wildr)
    }
    lapply(seq_len(ncol(combinations)), function(i) combinations[, i])
  })

  option_grid <- do.call(
    expand.grid,
    c(
      lapply(run_options, seq_along),
      stringsAsFactors = FALSE
    )
  )
  permutations <- matrix(
    0L,
    nrow = nrow(option_grid),
    ncol = nrow(sample_meta)
  )
  for (i in seq_len(nrow(option_grid))) {
    for (j in seq_along(run_options)) {
      selected <- run_options[[j]][[option_grid[i, j]]]
      permutations[i, selected] <- 1L
    }
  }
  permutations
}

global_t2_statistic <- function(y, design) {
  if (qr(design)$rank < ncol(design)) {
    return(NA_real_)
  }
  coef_index <- which(colnames(design) == "group_wildr")
  fit <- suppressWarnings(lmFit(y, design))
  t_value <- fit$coefficients[, coef_index] /
    (fit$stdev.unscaled[, coef_index] * fit$sigma)
  mean(pmin(t_value^2, GLOBAL_T2_CAP), na.rm = TRUE)
}

run_global_permutation <- function(family_data) {
  permutations <- enumerate_run_stratified_groups(family_data$sample_meta)
  observed_group <- as.integer(family_data$sample_meta$group == "wildr")
  observed_design <- make_design(family_data$sample_meta, observed_group)
  observed_statistic <- global_t2_statistic(
    family_data$transformed_matrix,
    observed_design
  )

  null_statistics <- vapply(
    seq_len(nrow(permutations)),
    function(i) {
      design <- make_design(family_data$sample_meta, permutations[i, ])
      global_t2_statistic(family_data$transformed_matrix, design)
    },
    numeric(1)
  )
  valid_null <- is.finite(null_statistics)
  p_value <- mean(
    null_statistics[valid_null] >= observed_statistic - sqrt(.Machine$double.eps)
  )

  data.table(
    n_tested_pairs = nrow(family_data$transformed_matrix),
    n_permutations = sum(valid_null),
    global_mean_t2 = observed_statistic,
    permutation_p = p_value,
    null_t2_median = median(null_statistics[valid_null]),
    null_t2_95 = quantile(null_statistics[valid_null], 0.95, names = FALSE)
  )
}

analyze_family <- function(
    all_rows,
    all_samples,
    scope_value,
    region_value,
    level_value,
    label_set_value
) {
  log_message(
    "Family: ", scope_value, " / ", region_value, " / ",
    level_value, " / ", label_set_value
  )

  family_rows <- all_rows[
    scope == scope_value &
      region == region_value &
      level == level_value &
      label_set == label_set_value
  ]
  sample_meta <- all_samples[
    scope == scope_value & region == region_value
  ]
  if (!nrow(family_rows) || !nrow(sample_meta)) {
    return(NULL)
  }
  setorder(sample_meta, run, group, roi)
  sample_meta[, sample_index := seq_len(.N)]

  available_groups <- sample_meta[, .N, by = group]
  if (
    nrow(available_groups) < 2L ||
      any(available_groups$N < MIN_GROUP_ROIS) ||
      uniqueN(sample_meta$run) < MIN_RUNS
  ) {
    log_message("Skipping family: insufficient ROI sections or runs.")
    return(NULL)
  }

  family_data <- build_family_matrix(family_rows, sample_meta)
  if (is.null(family_data)) {
    log_message("Skipping family: no eligible neighborhood pairs.")
    return(NULL)
  }

  pairwise <- fit_pairwise_models(family_data)
  global <- run_global_permutation(family_data)
  family_columns <- list(
    scope = scope_value,
    region = region_value,
    level = level_value,
    label_set = label_set_value
  )
  for (column in names(family_columns)) {
    pairwise[, (column) := family_columns[[column]]]
    global[, (column) := family_columns[[column]]]
  }
  setcolorder(
    pairwise,
    c(
      "scope", "region", "level", "label_set",
      setdiff(names(pairwise), c("scope", "region", "level", "label_set"))
    )
  )
  setcolorder(
    global,
    c(
      "scope", "region", "level", "label_set",
      setdiff(names(global), c("scope", "region", "level", "label_set"))
    )
  )

  eligibility <- data.table(
    scope = scope_value,
    region = region_value,
    level = level_value,
    label_set = label_set_value,
    n_available_spf = sum(sample_meta$group == "spf"),
    n_available_wildr = sum(sample_meta$group == "wildr"),
    n_available_runs = uniqueN(sample_meta$run),
    min_required_spf = family_data$min_spf,
    min_required_wildr = family_data$min_wildr,
    n_tested_pairs = nrow(pairwise)
  )

  list(pairwise = pairwise, global = global, eligibility = eligibility)
}

make_sensitivity_summary <- function(pairwise) {
  join_columns <- c(
    "scope", "region", "level",
    "query_label", "query_name", "neighbor_label", "neighbor_name"
  )
  primary <- pairwise[
    label_set == "primary_confident",
    c(
      as.list(.SD),
      list(
        delta_primary = delta_raw_adjusted,
        p_primary = p_value,
        fdr_primary = fdr_family
      )
    ),
    .SDcols = join_columns
  ]
  all_labels <- pairwise[
    label_set == "all_labels_sensitivity",
    c(
      as.list(.SD),
      list(
        delta_all_labels = delta_raw_adjusted,
        p_all_labels = p_value,
        fdr_all_labels = fdr_family
      )
    ),
    .SDcols = join_columns
  ]
  joined <- merge(primary, all_labels, by = join_columns)
  summary <- joined[
    ,
    .(
      n_matched_pairs = .N,
      spearman_effect = cor(
        delta_primary,
        delta_all_labels,
        method = "spearman",
        use = "complete.obs"
      ),
      pearson_effect = cor(
        delta_primary,
        delta_all_labels,
        method = "pearson",
        use = "complete.obs"
      ),
      sign_agreement = mean(
        sign(delta_primary) == sign(delta_all_labels),
        na.rm = TRUE
      )
    ),
    by = .(scope, region, level)
  ]
  list(joined = joined, summary = summary)
}

plot_global_results <- function(global) {
  primary <- global[label_set == "primary_confident"]
  primary[, region_display := ifelse(scope == "whole_slice", "Whole slice", region)]
  primary[, evidence := -log10(pmax(permutation_p, .Machine$double.xmin))]

  plot <- ggplot(
    primary,
    aes(
      x = evidence,
      y = reorder(region_display, evidence),
      shape = level,
      alpha = global_fdr < 0.10
    )
  ) +
    geom_vline(xintercept = -log10(0.05), linetype = 2) +
    geom_point(size = 3) +
    facet_grid(scope ~ level, scales = "free_y", space = "free_y") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.45), guide = "none") +
    labs(
      x = expression(-log[10]("run-stratified permutation p")),
      y = NULL,
      shape = "MapMyCells level",
      title = "Global SPF-WildR neighborhood-composition tests",
      subtitle = "ROI-section analysis; filled points have global BH-FDR < 0.10"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")

  ggsave(GLOBAL_FIGURE, plot, width = 11, height = 8, dpi = 200)
}

plot_top_effects <- function(pairwise) {
  top <- pairwise[
    label_set == "primary_confident" &
      is.finite(delta_raw_adjusted) &
      is.finite(p_value)
  ][
    order(p_value)
  ]
  top <- head(top, 30L)
  top[, pair_display := paste0(
    query_name, " -> ", neighbor_name, " [", region, "]"
  )]
  top[, pair_display := factor(pair_display, levels = rev(pair_display))]

  plot <- ggplot(
    top,
    aes(
      x = delta_raw_adjusted,
      y = pair_display,
      xmin = delta_raw_ci_low,
      xmax = delta_raw_ci_high,
      shape = level,
      alpha = fdr_family < 0.10
    )
  ) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbarh(height = 0) +
    geom_point(size = 2.5) +
    facet_grid(level ~ scope, scales = "free_y", space = "free_y") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.45), guide = "none") +
    labs(
      x = "Run-adjusted WildR - SPF difference in mean local relative abundance",
      y = NULL,
      shape = "MapMyCells level",
      title = "Largest exploratory SPF-WildR neighborhood signals",
      subtitle = "Thirty smallest moderated p-values; filled points have family FDR < 0.10"
    ) +
    theme_bw(base_size = 9) +
    theme(legend.position = "top")

  ggsave(EFFECT_FIGURE, plot, width = 14, height = 12, dpi = 200)
}

plot_sensitivity_results <- function(sensitivity_summary) {
  plot_data <- copy(sensitivity_summary)
  plot_data[, region_display := ifelse(scope == "whole_slice", "Whole slice", region)]
  plot_data[, region_display := factor(
    region_display,
    levels = unique(region_display[order(spearman_effect)])
  )]

  plot <- ggplot(
    plot_data,
    aes(
      x = spearman_effect,
      y = region_display,
      shape = level,
      size = n_matched_pairs
    )
  ) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point(alpha = 0.8) +
    facet_grid(scope ~ level, scales = "free_y", space = "free_y") +
    scale_x_continuous(limits = c(-1, 1)) +
    scale_size_area(max_size = 6) +
    labs(
      x = "Spearman correlation: primary vs all-label effect estimates",
      y = NULL,
      shape = "MapMyCells level",
      size = "Matched pairs",
      title = "Confidence-filter sensitivity of SPF-WildR effects"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")

  ggsave(SENSITIVITY_FIGURE, plot, width = 11, height = 8, dpi = 200)
}

write_top_hits_table <- function(pairwise, overwrite = FALSE) {
  if (!overwrite) {
    stop_if_exists(TOP_HITS_CSV)
  }
  top_hits <- pairwise[
    label_set == "primary_confident" &
      fdr_family < 0.10
  ][
    order(fdr_family, p_value)
  ][
    ,
    .(
      scope,
      region,
      level,
      query_name,
      neighbor_name,
      n_spf,
      n_wildr,
      mean_spf_raw,
      mean_wildr_raw,
      delta_raw_adjusted,
      delta_raw_ci_low,
      delta_raw_ci_high,
      p_value,
      fdr_family,
      fdr_studywide
    )
  ]
  fwrite(top_hits, TOP_HITS_CSV)
  invisible(top_hits)
}

run_analysis <- function() {
  outputs <- c(
    RESULT_RDS,
    PAIRWISE_CSV,
    GLOBAL_CSV,
    SENSITIVITY_CSV,
    ELIGIBILITY_CSV,
    TOP_HITS_CSV,
    GLOBAL_FIGURE,
    EFFECT_FIGURE,
    SENSITIVITY_FIGURE
  )
  stop_if_exists(outputs)
  dir.create(dirname(PAIRWISE_CSV), recursive = TRUE, showWarnings = FALSE)
  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

  log_message("Loading ROI summaries: ", INPUT_CSV)
  columns <- c(
    "scope", "region", "slice_id", "group", "run", "roi", "bregma",
    "level", "label_set", "query_label", "query_name", "neighbor_label",
    "neighbor_name", "n_query_total", "n_query_with_neighbors",
    "n_neighbor_region", "mean_rel_abund"
  )
  neighborhood_rows <- fread(INPUT_CSV, select = columns)
  neighborhood_rows[, `:=`(
    slice_id = tolower(slice_id),
    group = tolower(group)
  )]
  neighborhood_rows <- neighborhood_rows[
    !slice_id %chin% POOR_QUALITY_SLICES
  ]

  manifest <- fread(MANIFEST_CSV)
  manifest[, `:=`(
    slice_id = tolower(slice_id),
    group = tolower(group)
  )]
  sample_metadata <- build_sample_metadata(manifest)

  families <- unique(
    neighborhood_rows[, .(scope, region, level, label_set)]
  )
  setorder(families, label_set, scope, region, level)
  results <- vector("list", nrow(families))
  for (i in seq_len(nrow(families))) {
    results[[i]] <- analyze_family(
      neighborhood_rows,
      sample_metadata,
      families$scope[[i]],
      families$region[[i]],
      families$level[[i]],
      families$label_set[[i]]
    )
    invisible(gc())
  }
  results <- Filter(Negate(is.null), results)

  pairwise <- rbindlist(lapply(results, `[[`, "pairwise"), use.names = TRUE)
  global <- rbindlist(lapply(results, `[[`, "global"), use.names = TRUE)
  eligibility <- rbindlist(
    lapply(results, `[[`, "eligibility"),
    use.names = TRUE
  )
  pairwise[
    ,
    fdr_studywide := p.adjust(p_value, method = "BH"),
    by = label_set
  ]
  global[, global_fdr := p.adjust(permutation_p, method = "BH"), by = label_set]

  sensitivity <- make_sensitivity_summary(pairwise)

  saveRDS(
    list(
      parameters = list(
        analysis_unit = "ROI section (exploratory; animal mapping unavailable)",
        outcome = "mean_rel_abund",
        transform = "asin(sqrt(mean_rel_abund))",
        design = "~ run + group",
        min_query_cells = MIN_QUERY_CELLS,
        min_neighbor_cells = MIN_NEIGHBOR_CELLS,
        min_group_rois = MIN_GROUP_ROIS,
        min_runs = MIN_RUNS,
        global_t2_cap = GLOBAL_T2_CAP,
        poor_quality_slices = POOR_QUALITY_SLICES
      ),
      pairwise = pairwise,
      global = global,
      sensitivity_pairs = sensitivity$joined,
      sensitivity_summary = sensitivity$summary,
      eligibility = eligibility
    ),
    RESULT_RDS,
    compress = FALSE
  )
  fwrite(pairwise, PAIRWISE_CSV, compress = "gzip")
  fwrite(global, GLOBAL_CSV)
  fwrite(sensitivity$summary, SENSITIVITY_CSV)
  fwrite(eligibility, ELIGIBILITY_CSV)
  write_top_hits_table(pairwise)

  plot_global_results(global)
  plot_top_effects(pairwise)
  plot_sensitivity_results(sensitivity$summary)

  log_message("Pairwise tests: ", format(nrow(pairwise), big.mark = ","))
  log_message("Global families: ", nrow(global))
  log_message("Pairwise FDR < 0.10: ", sum(pairwise$fdr_family < 0.10, na.rm = TRUE))
  log_message(
    "Study-wide pairwise FDR < 0.10: ",
    sum(pairwise$fdr_studywide < 0.10, na.rm = TRUE)
  )
  log_message("Global FDR < 0.10: ", sum(global$global_fdr < 0.10, na.rm = TRUE))
  log_message("Saved: ", RESULT_RDS)
}

write_display_tables <- function() {
  if (!file.exists(RESULT_RDS)) {
    stop("Missing result artifact: ", RESULT_RDS, call. = FALSE)
  }
  result <- readRDS(RESULT_RDS)
  hits <- write_top_hits_table(result$pairwise)
  log_message("Display hits: ", nrow(hits))
  log_message("Saved: ", TOP_HITS_CSV)
}

self_test <- function() {
  sample_meta <- data.table(
    slice_id = paste0("s", 1:8),
    group = factor(
      rep(c("spf", "wildr"), 4),
      levels = c("spf", "wildr")
    ),
    run = factor(rep(paste0("run", 1:4), each = 2)),
    roi = paste0("roi", rep(1:2, 4)),
    sample_index = 1:8
  )
  permutations <- enumerate_run_stratified_groups(sample_meta)
  stopifnot(nrow(permutations) == 16L)
  stopifnot(all(rowSums(permutations) == 4L))

  set.seed(42)
  y <- matrix(runif(5 * 8, 0, 0.2), nrow = 5)
  y[1, sample_meta$group == "wildr"] <-
    y[1, sample_meta$group == "wildr"] + 0.5
  y <- asin(sqrt(y))
  observed <- global_t2_statistic(y, make_design(sample_meta))
  stopifnot(is.finite(observed), observed > 0)
  cat("SELF_TEST_OK\n")
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  action <- if (length(args)) args[[1]] else "help"
  if (action == "self-test") {
    self_test()
  } else if (action == "run") {
    run_analysis()
  } else if (action == "display-tables") {
    write_display_tables()
  } else {
    cat(
      "Usage:\n",
      "  run_cell_neighborhood_group_tests_v2.R self-test\n",
      "  run_cell_neighborhood_group_tests_v2.R run\n",
      "  run_cell_neighborhood_group_tests_v2.R display-tables\n",
      sep = ""
    )
  }
}

if (sys.nframe() == 0L) {
  main()
}
