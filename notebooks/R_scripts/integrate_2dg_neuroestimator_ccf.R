#!/usr/bin/env Rscript

# Join the latest seqFISH Allen CCF coordinates, spatial NEUROeSTIMator scores,
# and the accepted 2DG Allen-region summaries.
#
# Run with:
#   /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript \
#     notebooks/R_scripts/integrate_2dg_neuroestimator_ccf.R integrate
#
# This script is descriptive by design. The 2DG image is a group-level contrast,
# and the repository does not contain a verified seqFISH section-to-animal key.
# It therefore does not fit cell-level or section-level inferential models.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

WKDIR <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
ANALYSIS_DATE <- "2026-11-11"
EFFECT_ANALYSIS_DATE <- "2026-11-11"

CCF_PATH <- file.path(
  WKDIR,
  "data/interim/registration/Allen_CCF_rebuild/ccf3d.csv"
)
NEURO_PATH <- file.path(
  WKDIR,
  "data/interim/neuroestimator/seqFISH_meta_2025-05-01.rds"
)
DG_REGION_PATH <- file.path(
  WKDIR,
  "tables/2DG_CCF_Integration",
  paste0("2DG_CCF_region_summary_strict_", ANALYSIS_DATE, ".csv")
)
REGISTRATION_SUMMARY_PATH <- file.path(
  WKDIR,
  "data/interim/2DG_CCF_Integration",
  ANALYSIS_DATE,
  paste0("2dg_registration_summary_strict_", ANALYSIS_DATE, ".csv")
)
ANNOTATION_NAME_PATH <- file.path(
  WKDIR,
  "tables/2DG_CCF_Integration",
  paste0("Allen_CCF_annotation_names_hierarchy_", ANALYSIS_DATE, ".csv")
)

OUTPUT_DIR <- file.path(
  WKDIR,
  "data/interim/2DG_CCF_Integration",
  ANALYSIS_DATE
)
TABLE_DIR <- file.path(WKDIR, "tables/2DG_CCF_Integration")
FIGURE_DIR <- file.path(WKDIR, "figures/2DG_CCF_Integration")

JOINED_PATH <- file.path(
  OUTPUT_DIR,
  paste0("seqFISH_neuroestimator_ccf_strict_", ANALYSIS_DATE, ".rds")
)
SEQFISH_REGION_PATH <- file.path(
  TABLE_DIR,
  paste0("SeqFISH_NEUROeSTIMator_CCF_region_summary_strict_", ANALYSIS_DATE, ".csv.gz")
)
CONCORDANCE_PATH <- file.path(
  TABLE_DIR,
  paste0("2DG_SeqFISH_CCF_region_concordance_strict_", ANALYSIS_DATE, ".csv")
)
COVERAGE_PATH <- file.path(
  TABLE_DIR,
  paste0("2DG_SeqFISH_CCF_join_coverage_strict_", ANALYSIS_DATE, ".csv")
)
CONCORDANCE_FIGURE <- file.path(
  FIGURE_DIR,
  paste0("2DG_SeqFISH_CCF_region_concordance_strict_", ANALYSIS_DATE, ".png")
)
COVERAGE_FIGURE <- file.path(
  FIGURE_DIR,
  paste0("2DG_SeqFISH_CCF_join_coverage_strict_", ANALYSIS_DATE, ".png")
)
DG_OVERLAY_PATH <- file.path(
  OUTPUT_DIR,
  paste0(
    "2dg_overlay_points_ccf_scalar-t_",
    EFFECT_ANALYSIS_DATE,
    ".csv.gz"
  )
)
GATE_BRIDGE_PATH <- file.path(
  TABLE_DIR,
  paste0(
    "2DG_SeqFISH_annotation_gate_bridge_scalar-t_",
    EFFECT_ANALYSIS_DATE,
    ".csv.gz"
  )
)
GATED_EFFECT_PATH <- file.path(
  TABLE_DIR,
  paste0(
    "2DG_NEUROeSTIMator_gated_effects_scalar-t_",
    EFFECT_ANALYSIS_DATE,
    ".csv"
  )
)
GATED_CORRELATION_PATH <- file.path(
  TABLE_DIR,
  paste0(
    "2DG_NEUROeSTIMator_gated_effect_correlations_scalar-t_",
    EFFECT_ANALYSIS_DATE,
    ".csv"
  )
)
GATED_RANKING_FIGURE <- file.path(
  FIGURE_DIR,
  paste0(
    "2DG_NEUROeSTIMator_gated_effect_correlation_ranking_scalar-t_",
    EFFECT_ANALYSIS_DATE,
    ".png"
  )
)
GATED_BEST_FIGURE <- file.path(
  FIGURE_DIR,
  paste0(
    "2DG_NEUROeSTIMator_gated_effect_best_scatter_scalar-t_",
    EFFECT_ANALYSIS_DATE,
    ".png"
  )
)

log_message <- function(...) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
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

require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(
      "Required input(s) are missing:\n",
      paste(missing, collapse = "\n"),
      call. = FALSE
    )
  }
}

require_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    stop(
      label,
      " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

weighted_mean_safe <- function(x, weight) {
  keep <- is.finite(x) & is.finite(weight) & weight > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * weight[keep]) / sum(weight[keep])
}

hedges_g <- function(wildr, spf) {
  wildr <- wildr[is.finite(wildr)]
  spf <- spf[is.finite(spf)]
  n_wildr <- length(wildr)
  n_spf <- length(spf)
  if (n_wildr < 2L || n_spf < 2L) return(NA_real_)
  degrees_freedom <- n_wildr + n_spf - 2
  pooled_variance <- (
    (n_wildr - 1) * stats::var(wildr) +
      (n_spf - 1) * stats::var(spf)
  ) / degrees_freedom
  if (!is.finite(pooled_variance) || pooled_variance <= 0) return(NA_real_)
  cohen_d <- (mean(wildr) - mean(spf)) / sqrt(pooled_variance)
  correction <- 1 - 3 / (4 * degrees_freedom - 1)
  correction * cohen_d
}

integrate <- function() {
  require_files(c(
    CCF_PATH,
    NEURO_PATH,
    DG_REGION_PATH,
    REGISTRATION_SUMMARY_PATH,
    ANNOTATION_NAME_PATH
  ))
  outputs <- c(
    JOINED_PATH,
    SEQFISH_REGION_PATH,
    CONCORDANCE_PATH,
    COVERAGE_PATH,
    CONCORDANCE_FIGURE,
    COVERAGE_FIGURE
  )
  stop_if_exists(outputs)
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

  log_message("Loading NEUROeSTIMator metadata: ", NEURO_PATH)
  neuro <- as.data.table(readRDS(NEURO_PATH))
  require_columns(
    neuro,
    c("cell_id", "predicted_activity", "group", "roi", "run", "bregma"),
    "NEUROeSTIMator metadata"
  )
  neuro <- neuro[
    ,
    .(
      cell_id = as.character(cell_id),
      predicted_activity = as.numeric(predicted_activity),
      group = as.character(group),
      roi = as.character(roi),
      run = as.character(run),
      bregma = as.numeric(bregma),
      final_gate = if ("final_gate" %in% names(neuro)) as.character(final_gate) else NA_character_
    )
  ]
  if (anyDuplicated(neuro$cell_id)) {
    stop("NEUROeSTIMator cell_id is not unique.", call. = FALSE)
  }
  neuro[, sample_id := paste(roi, run, sep = "_")]

  log_message("Reading selected CCF columns: ", CCF_PATH)
  ccf <- fread(
    CCF_PATH,
    select = c("cell_id", "ccfx", "ccfy", "ccfz", "annotation", "sample_id")
  )
  require_columns(
    ccf,
    c("cell_id", "ccfx", "ccfy", "ccfz", "annotation"),
    "CCF table"
  )
  ccf[, cell_id := as.character(cell_id)]
  if (anyDuplicated(ccf$cell_id)) {
    stop("CCF cell_id is not unique.", call. = FALSE)
  }
  setnames(ccf, "sample_id", "ccf_sample_id")

  joined <- merge(neuro, ccf, by = "cell_id", all.x = TRUE, sort = FALSE)
  joined[, registered_ccf := !is.na(annotation)]
  joined[, sample_id_matches_ccf := sample_id == ccf_sample_id]

  matched <- joined[registered_ccf == TRUE]
  seqfish_by_sample <- matched[
    ,
    .(
      n_cells = .N,
      mean_predicted_activity = safe_mean(predicted_activity),
      median_predicted_activity = if (all(is.na(predicted_activity))) {
        NA_real_
      } else {
        median(predicted_activity, na.rm = TRUE)
      }
    ),
    by = .(annotation = as.integer(annotation), group, sample_id)
  ]
  seqfish_region <- seqfish_by_sample[
    ,
    .(
      n_samples = uniqueN(sample_id),
      n_cells = sum(n_cells),
      mean_predicted_activity = safe_mean(mean_predicted_activity),
      median_sample_activity = if (all(is.na(mean_predicted_activity))) {
        NA_real_
      } else {
        median(mean_predicted_activity, na.rm = TRUE)
      }
    ),
    by = .(annotation, group)
  ]
  annotation_names <- fread(ANNOTATION_NAME_PATH)
  require_columns(
    annotation_names,
    c("annotation", "annotation_name"),
    "Allen annotation-name map"
  )
  annotation_names[, annotation := as.integer(annotation)]
  seqfish_region <- merge(
    seqfish_region,
    annotation_names,
    by = "annotation",
    all.x = TRUE
  )

  # A descriptive WildR-minus-SPF contrast of section means. This is not an
  # animal-level test and intentionally carries no p-value.
  activity_wide <- dcast(
    seqfish_region,
    annotation ~ group,
    value.var = c("mean_predicted_activity", "n_samples", "n_cells")
  )
  wildr_column <- grep(
    "^mean_predicted_activity_.*wild",
    names(activity_wide),
    value = TRUE,
    ignore.case = TRUE
  )
  spf_column <- grep(
    "^mean_predicted_activity_.*spf",
    names(activity_wide),
    value = TRUE,
    ignore.case = TRUE
  )
  if (length(wildr_column) != 1L || length(spf_column) != 1L) {
    stop(
      "Could not identify exactly one WildR and one SPF activity column. Found: ",
      paste(names(activity_wide), collapse = ", "),
      call. = FALSE
    )
  }
  activity_wide[
    ,
    seqfish_wildr_minus_spf_descriptive :=
      get(wildr_column) - get(spf_column)
  ]

  dg <- fread(DG_REGION_PATH)
  require_columns(
    dg,
    c(
      "annotation",
      "annotation_name",
      "n_slices",
      "n_overlay_pixels",
      "mean_overlay_signed_relative"
    ),
    "2DG region summary"
  )
  dg[, annotation := as.integer(annotation)]
  dg[, annotation_name_export := annotation_name]
  dg[, annotation_name := NULL]
  dg <- merge(dg, annotation_names, by = "annotation", all.x = TRUE)
  concordance <- merge(dg, activity_wide, by = "annotation", all = FALSE)
  concordance[
    ,
    interpretation_status :=
      "descriptive_only_2DG_group_map_and_no_verified_seqFISH_animal_key"
  ]
  concordance[
    ,
    t_value_status := "scalar_tiff_y_and_email_formula"
  ]

  registration <- fread(REGISTRATION_SUMMARY_PATH)
  coverage <- data.table(
    metric = c(
      "NEUROeSTIMator cells",
      "Cells joined to current CCF",
      "Cells with matching derived/current sample ID",
      "Accepted 2DG slices",
      "2DG Allen regions",
      "Shared Allen regions"
    ),
    value = c(
      nrow(joined),
      sum(joined$registered_ccf),
      sum(joined$sample_id_matches_ccf, na.rm = TRUE),
      sum(registration$accepted_final),
      uniqueN(dg$annotation),
      uniqueN(concordance$annotation)
    )
  )

  saveRDS(joined, JOINED_PATH, compress = "xz")
  fwrite(seqfish_region, SEQFISH_REGION_PATH)
  fwrite(concordance, CONCORDANCE_PATH)
  fwrite(coverage, COVERAGE_PATH)

  plot_coverage(coverage)
  plot_concordance(concordance)
  log_message("Joined object: ", JOINED_PATH)
  log_message("Shared-region table: ", CONCORDANCE_PATH)
}

plot_coverage <- function(coverage) {
  plot_data <- copy(coverage)
  plot_data[, metric := factor(metric, levels = rev(metric))]
  p <- ggplot(plot_data, aes(x = value, y = metric)) +
    geom_col(fill = "#FF6C0C", width = 0.72) +
    geom_text(
      aes(label = format(value, big.mark = ",", scientific = FALSE)),
      hjust = -0.08,
      size = 3.4
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      title = "Cross-modality Allen CCF coverage",
      x = "Count",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank())
  ggsave(COVERAGE_FIGURE, p, width = 9, height = 4.8, dpi = 220)
}

plot_concordance <- function(concordance) {
  x_column <- if (
    "mean_t_value" %in% names(concordance) &&
      any(is.finite(concordance$mean_t_value))
  ) {
    "mean_t_value"
  } else {
    "mean_overlay_signed_relative"
  }
  plot_data <- concordance[
    is.finite(get(x_column)) &
      is.finite(seqfish_wildr_minus_spf_descriptive)
  ]
  if (!nrow(plot_data)) {
    stop("No finite shared-region values are available for plotting.", call. = FALSE)
  }
  plot_data[, two_dg_x := get(x_column)]
  plot_data[
    ,
    label_rank := frank(
      -abs(two_dg_x * seqfish_wildr_minus_spf_descriptive),
      ties.method = "first"
    )
  ]
  p <- ggplot(
    plot_data,
    aes(
      x = two_dg_x,
      y = seqfish_wildr_minus_spf_descriptive,
      size = pmax(n_overlay_pixels, 1)
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_point(alpha = 0.58, color = "#FF6C0C") +
    geom_text(
      data = plot_data[label_rank <= 12],
      aes(label = annotation_name),
      size = 2.5,
      check_overlap = TRUE,
      vjust = -0.8,
      show.legend = FALSE
    ) +
    scale_size_area(max_size = 8, guide = "none") +
    labs(
      title = "Descriptive regional alignment of 2DG and NEUROeSTIMator",
      subtitle = "No inferential p-value: 2DG is a group map and seqFISH animal identity is unverified",
      x = if (identical(x_column, "mean_t_value")) {
        "WildRvSPF Intruder mean t (scalar TIFF)"
      } else {
        "Intruder contrast overlay (signed relative display magnitude)"
      },
      y = "NEUROeSTIMator WildR − SPF (mean section activity)"
    ) +
    theme_minimal(base_size = 12)
  ggsave(CONCORDANCE_FIGURE, p, width = 9, height = 6.5, dpi = 220)
}

calculate_gated_effect_correlations <- function() {
  require_files(c(JOINED_PATH, DG_OVERLAY_PATH, DG_REGION_PATH))
  outputs <- c(
    GATE_BRIDGE_PATH,
    GATED_EFFECT_PATH,
    GATED_CORRELATION_PATH,
    GATED_RANKING_FIGURE,
    GATED_BEST_FIGURE
  )
  stop_if_exists(outputs)
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

  log_message("Loading joined seqFISH/CCF object: ", JOINED_PATH)
  joined <- as.data.table(readRDS(JOINED_PATH))
  require_columns(
    joined,
    c(
      "predicted_activity",
      "group",
      "sample_id",
      "final_gate",
      "annotation",
      "registered_ccf"
    ),
    "Joined seqFISH/CCF object"
  )
  joined[, final_gate := trimws(as.character(final_gate))]
  joined[, group := tolower(as.character(group))]

  # Section-level means are the least pseudoreplicated summary currently
  # possible. They are not animal-level replicates because the project lacks a
  # verified section-to-animal key.
  section_gate <- joined[
    !is.na(final_gate) & final_gate != "" & is.finite(predicted_activity),
    .(
      n_cells = .N,
      mean_predicted_activity = mean(predicted_activity)
    ),
    by = .(final_gate, group, sample_id)
  ]
  neuro_effect <- section_gate[
    ,
    {
      wildr <- mean_predicted_activity[group == "wildr"]
      spf <- mean_predicted_activity[group == "spf"]
      list(
        n_sections_wildr = length(wildr),
        n_sections_spf = length(spf),
        n_cells = sum(n_cells),
        neuro_wildr_mean = safe_mean(wildr),
        neuro_spf_mean = safe_mean(spf),
        neuro_raw_difference = safe_mean(wildr) - safe_mean(spf),
        neuro_hedges_g = hedges_g(wildr, spf)
      )
    },
    by = final_gate
  ]

  # Estimate P(gate | Allen annotation) using all registered seqFISH cells.
  # Ungated cells remain in the denominator so an annotation represented only
  # weakly by a manual gate cannot donate all of its 2DG signal to that gate.
  bridge_source <- joined[
    registered_ccf == TRUE & !is.na(annotation) & annotation != 0
  ]
  bridge_source[
    ,
    bridge_gate := fifelse(
      !is.na(final_gate) & final_gate != "",
      final_gate,
      "__UNGATED__"
    )
  ]
  bridge <- bridge_source[
    ,
    .(n_seqfish_cells = .N),
    by = .(annotation = as.integer(annotation), final_gate = bridge_gate)
  ]
  bridge[
    ,
    annotation_total_cells := sum(n_seqfish_cells),
    by = annotation
  ]
  bridge[
    ,
    p_gate_given_annotation := n_seqfish_cells / annotation_total_cells
  ]
  bridge_gate <- bridge[final_gate != "__UNGATED__"]
  fwrite(bridge_gate, GATE_BRIDGE_PATH)

  log_message("Loading CCF-mapped 2DG overlay: ", DG_OVERLAY_PATH)
  overlay <- fread(
    DG_OVERLAY_PATH,
    select = c("annotation", "t_value")
  )
  overlay[, annotation := as.integer(annotation)]
  overlay_allocated <- merge(
    overlay,
    bridge_gate[
      ,
      .(annotation, final_gate, p_gate_given_annotation)
    ],
    by = "annotation",
    allow.cartesian = TRUE
  )
  overlay_allocated[
    ,
    positive := t_value > 0
  ]
  overlay_allocated[
    ,
    negative := t_value < 0
  ]
  dg_overlay_effect <- overlay_allocated[
    ,
    .(
      allocated_overlay_pixels = sum(p_gate_given_annotation),
      two_dg_mean_t = weighted_mean_safe(
        t_value,
        p_gate_given_annotation
      ),
      two_dg_mean_absolute_t = weighted_mean_safe(
        abs(t_value),
        p_gate_given_annotation
      ),
      two_dg_t_sign_balance = (
        sum(p_gate_given_annotation[positive]) -
          sum(p_gate_given_annotation[negative])
      ) / sum(p_gate_given_annotation),
      two_dg_positive_mean_t = weighted_mean_safe(
        t_value[positive],
        p_gate_given_annotation[positive]
      ),
      two_dg_negative_mean_t = weighted_mean_safe(
        t_value[negative],
        p_gate_given_annotation[negative]
      ),
      weighted_signed_signal = sum(
        p_gate_given_annotation * t_value
      ),
      weighted_absolute_signal = sum(
        p_gate_given_annotation * abs(t_value)
      )
    ),
    by = final_gate
  ]

  dg_regions <- fread(
    DG_REGION_PATH,
    select = c("annotation", "n_template_pixels", "n_overlay_pixels")
  )
  dg_regions[, annotation := as.integer(annotation)]
  dg_template_allocated <- merge(
    dg_regions,
    bridge_gate[
      ,
      .(annotation, final_gate, p_gate_given_annotation)
    ],
    by = "annotation",
    allow.cartesian = TRUE
  )
  dg_denominator <- dg_template_allocated[
    ,
    .(
      allocated_template_pixels = sum(
        p_gate_given_annotation * n_template_pixels
      ),
      allocated_overlay_pixels_from_regions = sum(
        p_gate_given_annotation * n_overlay_pixels
      )
    ),
    by = final_gate
  ]
  dg_effect <- merge(
    dg_overlay_effect,
    dg_denominator,
    by = "final_gate",
    all = TRUE
  )
  dg_effect[
    ,
    two_dg_signed_t_density :=
      weighted_signed_signal / allocated_template_pixels
  ]
  dg_effect[
    ,
    two_dg_absolute_t_density :=
      weighted_absolute_signal / allocated_template_pixels
  ]
  dg_effect[
    ,
    two_dg_overlay_density :=
      allocated_overlay_pixels_from_regions / allocated_template_pixels
  ]

  gated_effect <- merge(neuro_effect, dg_effect, by = "final_gate", all = FALSE)
  gated_effect[
    ,
    interpretation_status :=
      "exploratory_section_level_NEUROeSTIMator_and_scalar_2DG_t"
  ]
  gated_effect[
    ,
    t_value_status := "scalar_tiff_y_and_email_formula"
  ]
  fwrite(gated_effect, GATED_EFFECT_PATH)

  metric_columns <- grep("^two_dg_", names(gated_effect), value = TRUE)
  neuro_columns <- c("neuro_hedges_g", "neuro_raw_difference")
  correlation_rows <- list()
  set.seed(42)
  row_index <- 1L
  for (neuro_column in neuro_columns) {
    for (metric_column in metric_columns) {
      pair <- gated_effect[
        is.finite(get(neuro_column)) & is.finite(get(metric_column)),
        .(
          final_gate,
          neuro_effect = get(neuro_column),
          two_dg_effect = get(metric_column)
        )
      ]
      if (nrow(pair) < 5L || uniqueN(pair$two_dg_effect) < 3L) next
      pearson <- cor(
        pair$two_dg_effect,
        pair$neuro_effect,
        method = "pearson"
      )
      spearman <- cor(
        pair$two_dg_effect,
        pair$neuro_effect,
        method = "spearman"
      )
      permuted <- replicate(
        10000,
        cor(
          sample(pair$two_dg_effect),
          pair$neuro_effect,
          method = "spearman"
        )
      )
      permutation_p <- (
        1 + sum(abs(permuted) >= abs(spearman), na.rm = TRUE)
      ) / (1 + sum(is.finite(permuted)))
      loo <- vapply(
        seq_len(nrow(pair)),
        function(index) {
          cor(
            pair$two_dg_effect[-index],
            pair$neuro_effect[-index],
            method = "spearman"
          )
        },
        numeric(1)
      )
      correlation_rows[[row_index]] <- data.table(
        neuro_effect = neuro_column,
        two_dg_effect = metric_column,
        n_gates = nrow(pair),
        pearson_r = pearson,
        spearman_rho = spearman,
        permutation_p = permutation_p,
        loo_spearman_min = min(loo, na.rm = TRUE),
        loo_spearman_median = median(loo, na.rm = TRUE),
        loo_spearman_max = max(loo, na.rm = TRUE)
      )
      row_index <- row_index + 1L
    }
  }
  correlations <- rbindlist(correlation_rows)
  correlations[
    ,
    permutation_fdr := p.adjust(permutation_p, method = "BH"),
    by = neuro_effect
  ]
  correlations[, abs_spearman_rho := abs(spearman_rho)]
  setorder(
    correlations,
    neuro_effect,
    -abs_spearman_rho
  )
  correlations[
    ,
    rank_within_neuro_effect := seq_len(.N),
    by = neuro_effect
  ]
  fwrite(correlations, GATED_CORRELATION_PATH)

  plot_gated_correlation_ranking(correlations)
  plot_best_gated_correlation(correlations, gated_effect)
  log_message("Gated effects: ", GATED_EFFECT_PATH)
  log_message("Gated correlation ranking: ", GATED_CORRELATION_PATH)
}

plot_gated_correlation_ranking <- function(correlations) {
  ranking <- correlations[neuro_effect == "neuro_hedges_g"]
  ranking[, effect_label := gsub("^two_dg_", "", two_dg_effect)]
  ranking[, effect_label := gsub("_", " ", effect_label)]
  ranking[
    ,
    effect_label := factor(
      effect_label,
      levels = effect_label[order(abs(spearman_rho))]
    )
  ]
  plot_data <- melt(
    ranking,
    id.vars = c("effect_label", "n_gates"),
    measure.vars = c("spearman_rho", "pearson_r"),
    variable.name = "correlation",
    value.name = "coefficient"
  )
  plot_data[
    ,
    correlation := factor(
      correlation,
      levels = c("spearman_rho", "pearson_r"),
      labels = c("Spearman", "Pearson")
    )
  ]
  p <- ggplot(
    plot_data,
    aes(x = coefficient, y = effect_label, fill = correlation)
  ) +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "grey70") +
    geom_col(position = position_dodge(width = 0.75), width = 0.68) +
    scale_x_continuous(limits = c(-1, 1)) +
    labs(
      title = "Which scalar 2DG t effect best tracks gated NEUROeSTIMator effects?",
      subtitle = "Correlations across manual gates; WildRvSPF Intruder scalar y converted with the supplied t formula",
      x = "Correlation coefficient",
      y = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank())
  ggsave(GATED_RANKING_FIGURE, p, width = 10, height = 6.5, dpi = 220)
}

plot_best_gated_correlation <- function(correlations, gated_effect) {
  best <- correlations[
    neuro_effect == "neuro_hedges_g"
  ][order(-abs(spearman_rho))][1]
  metric <- best$two_dg_effect
  plot_data <- gated_effect[
    is.finite(neuro_hedges_g) & is.finite(get(metric)),
    .(
      final_gate,
      neuro_hedges_g,
      two_dg_effect = get(metric),
      n_sections_wildr,
      n_sections_spf
    )
  ]
  metric_label <- gsub("_", " ", sub("^two_dg_", "", metric))
  p <- ggplot(
    plot_data,
    aes(x = two_dg_effect, y = neuro_hedges_g)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, color = "grey45") +
    geom_point(size = 3, color = "#FF6C0C") +
    geom_text(
      aes(label = final_gate),
      vjust = -0.8,
      size = 3,
      check_overlap = TRUE
    ) +
    labs(
      title = paste("Best gated-region scalar 2DG metric:", metric_label),
      subtitle = sprintf(
        "Across %d gates: Spearman rho = %.2f; permutation p = %.3f",
        best$n_gates,
        best$spearman_rho,
        best$permutation_p
      ),
      x = paste("Decoded 2DG", metric_label),
      y = "NEUROeSTIMator WildR − SPF Hedges' g"
    ) +
    theme_minimal(base_size = 12)
  ggsave(GATED_BEST_FIGURE, p, width = 8.5, height = 6.5, dpi = 220)
}

report_status <- function() {
  paths <- c(
    registration_summary = REGISTRATION_SUMMARY_PATH,
    dg_region_summary = DG_REGION_PATH,
    joined_object = JOINED_PATH,
    concordance = CONCORDANCE_PATH,
    coverage = COVERAGE_PATH,
    gated_effects = GATED_EFFECT_PATH,
    gated_correlations = GATED_CORRELATION_PATH
  )
  print(data.table(
    artifact = names(paths),
    path = unname(paths),
    exists = file.exists(paths)
  ))
}

usage <- function() {
  cat(
    "Usage:\n",
    "  integrate_2dg_neuroestimator_ccf.R integrate\n",
    "  integrate_2dg_neuroestimator_ccf.R gated-effect-correlation\n",
    "  integrate_2dg_neuroestimator_ccf.R status\n",
    sep = ""
  )
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) {
    usage()
    quit(status = 1)
  }
  if (identical(args[[1]], "integrate")) {
    integrate()
  } else if (identical(args[[1]], "gated-effect-correlation")) {
    calculate_gated_effect_correlations()
  } else if (identical(args[[1]], "status")) {
    report_status()
  } else {
    usage()
    quit(status = 1)
  }
}

if (sys.nframe() == 0L) {
  main()
}
