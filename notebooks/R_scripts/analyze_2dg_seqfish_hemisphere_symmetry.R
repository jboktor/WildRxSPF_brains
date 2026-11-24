#!/usr/bin/env Rscript

# Hemisphere-aware bilateral summaries + symmetry signatures for
# 2DG (Baseline / Intruder / I−B) and seqFISH NEUROeSTIMator.
#
# Follows analysis_md/2DG_seqFISH_hemisphere_symmetry_plan.md §5–7 (steps C–F).
# Does not overwrite Allen rebuild tables; writes dated add-ons only.
#
# /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript \
#   notebooks/R_scripts/analyze_2dg_seqfish_hemisphere_symmetry.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

WKDIR <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
SOURCE_DATE <- "2026-11-23c"
CONTROL_DATE <- "2026-11-23c"
ANALYSIS_DATE <- "2026-11-23c"

SOURCE_DIR <- file.path(WKDIR, "data/interim/2DG_CCF_Integration", SOURCE_DATE)
CONTROL_DIR <- file.path(WKDIR, "data/interim/2DG_CCF_Integration", CONTROL_DATE)
OUT_DIR <- file.path(WKDIR, "data/interim/2DG_CCF_Integration", ANALYSIS_DATE)
FIGURE_DIR <- file.path(WKDIR, "figures/2DG_CCF_Integration")
TABLE_DIR <- file.path(WKDIR, "tables/2DG_CCF_Integration")

JOINED <- file.path(SOURCE_DIR, paste0("seqFISH_neuroestimator_ccf_strict_", SOURCE_DATE, ".rds"))
HEMISPHERE_PROVENANCE <- file.path(
  TABLE_DIR, "seqFISH_section_hemisphere_provenance_2026-11-14.csv"
)
LIMBIC_SET <- file.path(
  TABLE_DIR, paste0("2DG_prespecified_limbic_hypothalamic_annotations_", CONTROL_DATE, ".csv")
)

INTRUDER_OVERLAY <- file.path(
  SOURCE_DIR, paste0("2dg_overlay_points_ccf_scalar-t_", SOURCE_DATE, ".csv.gz")
)
BASELINE_OVERLAY <- file.path(
  CONTROL_DIR, paste0("2dg_overlay_points_ccf_baseline_wildr_vs_spf_", CONTROL_DATE, ".csv.gz")
)
WILDR_IB_OVERLAY <- file.path(
  CONTROL_DIR, paste0("2dg_overlay_points_ccf_wildr_intruder_vs_baseline_", CONTROL_DATE, ".csv.gz")
)
SPF_IB_OVERLAY <- file.path(
  CONTROL_DIR, paste0("2dg_overlay_points_ccf_spf_intruder_vs_baseline_", CONTROL_DATE, ".csv.gz")
)

SEQ_HEMI_RDS <- file.path(OUT_DIR, paste0("seqFISH_neuroestimator_ccf_with_hemisphere_", ANALYSIS_DATE, ".rds"))
SEQ_HEMI_CSV <- file.path(TABLE_DIR, paste0("seqFISH_CCF_hemisphere_cell_summary_", ANALYSIS_DATE, ".csv.gz"))
TWODG_BILAT <- file.path(TABLE_DIR, paste0("2DG_region_bilateral_summaries_", ANALYSIS_DATE, ".csv"))
NE_BILAT <- file.path(TABLE_DIR, paste0("SeqFISH_NE_region_bilateral_summaries_", ANALYSIS_DATE, ".csv"))
SYMM_TABLE <- file.path(TABLE_DIR, paste0("2DG_NE_hemisphere_symmetry_classes_", ANALYSIS_DATE, ".csv"))
SYMM_CONFIG <- file.path(OUT_DIR, paste0("hemisphere_symmetry_thresholds_", ANALYSIS_DATE, ".json"))
FIG_CLASSES <- file.path(FIGURE_DIR, paste0("2DG_NE_hemisphere_symmetry_classes_", ANALYSIS_DATE, ".png"))
FIG_LAT <- file.path(FIGURE_DIR, paste0("2DG_NE_lateralization_concordance_", ANALYSIS_DATE, ".png"))
FIG_COUNTS <- file.path(FIGURE_DIR, paste0("2DG_NE_symmetry_class_counts_", ANALYSIS_DATE, ".png"))

ML_MIDLINE_UM <- 5700
MIDLINE_BAND_UM <- 250

# Locked before inspecting class counts (plan §5.1 / §9).
TAU_2DG <- 0.75
TAU0_2DG <- 0.25
DELTA_REL <- 0.35
TAU_NE <- 0.004
TAU0_NE <- 0.0015
MIN_2DG_VOXELS <- 25L
MIN_NE_CELLS_SIDE <- 50L
EPS <- 1e-6

log_message <- function(...) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
}

stop_if_exists <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing)) {
    stop("Refusing to overwrite:\n", paste(existing, collapse = "\n"), call. = FALSE)
  }
}

hemisphere_from_ml <- function(ml) {
  fifelse(
    abs(ml - ML_MIDLINE_UM) <= MIDLINE_BAND_UM,
    "midline_band",
    fifelse(ml < ML_MIDLINE_UM, "left", "right")
  )
}

apply_anatomical_ml <- function(dt, sample_col = "sample_id") {
  prov <- fread(HEMISPHERE_PROVENANCE)
  prov[, sample_id := as.character(sample_id)]
  prov[, reflect := as.logical(reflect)]
  out <- copy(dt)
  out[, (sample_col) := as.character(get(sample_col))]
  out <- merge(
    out,
    unique(prov[, .(sample_id, reflect)]),
    by.x = sample_col,
    by.y = "sample_id",
    all.x = TRUE
  )
  out[is.na(reflect), reflect := FALSE]
  out[, ccfz_registered := ccfz]
  out[reflect == TRUE, ccfz_anatomical := (2 * ML_MIDLINE_UM) - ccfz_registered]
  out[reflect == FALSE, ccfz_anatomical := ccfz_registered]
  out[, anatomical_hemisphere := hemisphere_from_ml(ccfz_anatomical)]
  out[, registration_hemisphere := hemisphere_from_ml(ccfz_registered)]
  out
}

classify_symmetry <- function(a, b, tau, tau0, delta = DELTA_REL, ok = TRUE) {
  if (!isTRUE(ok) || !is.finite(a) || !is.finite(b)) return("indeterminate")
  aa <- abs(a)
  bb <- abs(b)
  if (aa < tau0 && bb < tau0) return("indeterminate")
  if (sign(a) == sign(b) && aa >= tau && bb >= tau) {
    rel <- abs(a - b) / (aa + bb + EPS)
    if (rel <= delta) return("symmetric")
    return("indeterminate") # same sign but unequal magnitude
  }
  if ((aa >= tau && bb <= tau0) || (bb >= tau && aa <= tau0)) return("unilateral")
  if (sign(a) == -sign(b) && aa >= tau && bb >= tau) return("opposing")
  "indeterminate"
}

lateralization_index <- function(left, right) {
  (right - left) / (abs(right) + abs(left) + EPS)
}

summarize_2dg_map <- function(path, map_name) {
  dt <- fread(
    path,
    select = c("ccfz", "t_value", "annotation", "annotation_name")
  )
  dt <- dt[is.finite(t_value) & is.finite(ccfz) & !is.na(annotation) & annotation != 0]
  dt[, hemisphere := hemisphere_from_ml(ccfz)]
  dt <- dt[hemisphere %in% c("left", "right")]
  by_hemi <- dt[
    ,
    .(
      n_voxels = .N,
      mean_t = mean(t_value),
      median_t = median(t_value),
      mean_abs_t = mean(abs(t_value))
    ),
    by = .(annotation = as.integer(annotation), annotation_name, hemisphere)
  ]
  wide <- dcast(
    by_hemi,
    annotation + annotation_name ~ hemisphere,
    value.var = c("n_voxels", "mean_t", "median_t", "mean_abs_t")
  )
  for (col in c(
    "n_voxels_left", "n_voxels_right",
    "mean_t_left", "mean_t_right",
    "median_t_left", "median_t_right",
    "mean_abs_t_left", "mean_abs_t_right"
  )) {
    if (!col %in% names(wide)) wide[, (col) := NA_real_]
  }
  wide[, map_family := map_name]
  wide[, n_voxels_bilateral := n_voxels_left + n_voxels_right]
  wide[, mean_t_bilateral := {
    wL <- n_voxels_left
    wR <- n_voxels_right
    ((wL * mean_t_left) + (wR * mean_t_right)) / pmax(wL + wR, 1)
  }]
  wide[, mean_t_L_minus_R := mean_t_left - mean_t_right]
  wide[, lateralization_index := lateralization_index(mean_t_left, mean_t_right)]
  wide[, qc_ok := n_voxels_left >= MIN_2DG_VOXELS & n_voxels_right >= MIN_2DG_VOXELS]
  wide[
    ,
    symmetry_class := mapply(
      classify_symmetry,
      mean_t_left,
      mean_t_right,
      MoreArgs = list(tau = TAU_2DG, tau0 = TAU0_2DG),
      ok = qc_ok
    )
  ]
  wide
}

summarize_ne <- function(cells) {
  dt <- cells[
    anatomical_hemisphere %in% c("left", "right") &
      is.finite(predicted_activity) &
      !is.na(annotation) &
      annotation != 0
  ]
  # Per hemisphere × group means, then WildR − SPF
  side_group <- dt[
    ,
    .(
      n_cells = .N,
      mean_activity = mean(predicted_activity),
      n_samples = uniqueN(sample_id)
    ),
    by = .(
      annotation = as.integer(annotation),
      annotation_name,
      anatomical_hemisphere,
      group = as.character(group)
    )
  ]
  # Keep only spf/wildr
  side_group <- side_group[group %in% c("spf", "wildr", "SPF", "WildR", "wildR")]
  side_group[, group := tolower(group)]
  side_group[group == "wildr", group := "wildr"]
  wide_g <- dcast(
    side_group,
    annotation + annotation_name + anatomical_hemisphere ~ group,
    value.var = c("n_cells", "mean_activity", "n_samples")
  )
  for (needed in c(
    "n_cells_spf", "n_cells_wildr",
    "mean_activity_spf", "mean_activity_wildr",
    "n_samples_spf", "n_samples_wildr"
  )) {
    if (!needed %in% names(wide_g)) wide_g[, (needed) := NA_real_]
  }
  wide_g[, ne_delta := mean_activity_wildr - mean_activity_spf]
  wide_g[, n_cells_total := n_cells_spf + n_cells_wildr]
  hemi <- dcast(
    wide_g,
    annotation + annotation_name ~ anatomical_hemisphere,
    value.var = c(
      "ne_delta", "n_cells_total", "n_cells_spf", "n_cells_wildr",
      "mean_activity_spf", "mean_activity_wildr"
    )
  )
  for (col in names(hemi)) {
    if (is.numeric(hemi[[col]])) hemi[is.na(get(col)), (col) := NA_real_]
  }
  # Ensure expected columns
  expect <- c(
    "ne_delta_left", "ne_delta_right",
    "n_cells_total_left", "n_cells_total_right"
  )
  for (col in expect) if (!col %in% names(hemi)) hemi[, (col) := NA_real_]
  hemi[, lateralization_index := lateralization_index(ne_delta_left, ne_delta_right)]
  hemi[, qc_ok := n_cells_total_left >= MIN_NE_CELLS_SIDE &
    n_cells_total_right >= MIN_NE_CELLS_SIDE &
    is.finite(ne_delta_left) & is.finite(ne_delta_right)]
  hemi[
    ,
    symmetry_class := mapply(
      classify_symmetry,
      ne_delta_left,
      ne_delta_right,
      MoreArgs = list(tau = TAU_NE, tau0 = TAU0_NE),
      ok = qc_ok
    )
  ]
  # Also available-side only rows for hemibrain honesty
  avail <- wide_g[
    ,
    .(
      annotation,
      annotation_name,
      anatomical_hemisphere,
      ne_delta,
      n_cells_total,
      n_cells_spf,
      n_cells_wildr
    )
  ]
  list(bilateral = hemi, available_side = avail)
}

build_concordance <- function(twodg, ne_bilat, limbic) {
  intr <- twodg[map_family == "intruder_wildr_vs_spf"]
  base <- twodg[map_family == "baseline_wildr_vs_spf"]
  rename_existing <- function(dt, old, new) {
    keep <- old %in% names(dt)
    if (any(keep)) setnames(dt, old[keep], new[keep])
    dt
  }
  intr <- rename_existing(
    intr,
    c("mean_t_left", "mean_t_right", "lateralization_index", "symmetry_class", "qc_ok",
      "n_voxels_left", "n_voxels_right", "mean_t_bilateral", "mean_t_L_minus_R"),
    c("intruder_mean_t_left", "intruder_mean_t_right", "intruder_lat_index",
      "twodg_intruder_symmetry_class", "intruder_qc_ok",
      "intruder_n_voxels_left", "intruder_n_voxels_right",
      "intruder_mean_t_bilateral", "intruder_mean_t_L_minus_R")
  )
  base <- rename_existing(
    base,
    c("mean_t_left", "mean_t_right", "lateralization_index", "symmetry_class", "qc_ok",
      "n_voxels_left", "n_voxels_right", "mean_t_bilateral", "mean_t_L_minus_R"),
    c("baseline_mean_t_left", "baseline_mean_t_right", "baseline_lat_index",
      "twodg_baseline_symmetry_class", "baseline_qc_ok",
      "baseline_n_voxels_left", "baseline_n_voxels_right",
      "baseline_mean_t_bilateral", "baseline_mean_t_L_minus_R")
  )
  ne <- copy(ne_bilat)
  ne <- rename_existing(
    ne,
    c("lateralization_index", "symmetry_class", "qc_ok",
      "n_cells_total_left", "n_cells_total_right"),
    c("ne_lat_index", "ne_symmetry_class", "ne_qc_ok",
      "ne_n_cells_left", "ne_n_cells_right")
  )
  out <- merge(
    intr[, .(
      annotation, annotation_name,
      intruder_mean_t_left, intruder_mean_t_right, intruder_mean_t_bilateral,
      intruder_mean_t_L_minus_R, intruder_lat_index, twodg_intruder_symmetry_class,
      intruder_qc_ok, intruder_n_voxels_left, intruder_n_voxels_right
    )],
    base[, .(
      annotation,
      baseline_mean_t_left, baseline_mean_t_right, baseline_mean_t_bilateral,
      baseline_mean_t_L_minus_R, baseline_lat_index, twodg_baseline_symmetry_class,
      baseline_qc_ok, baseline_n_voxels_left, baseline_n_voxels_right
    )],
    by = "annotation",
    all = TRUE
  )
  out <- merge(
    out,
    ne[, .(
      annotation, annotation_name_ne = annotation_name,
      ne_delta_left, ne_delta_right, ne_lat_index, ne_symmetry_class,
      ne_qc_ok, ne_n_cells_left, ne_n_cells_right
    )],
    by = "annotation",
    all = TRUE
  )
  out[is.na(annotation_name) & !is.na(annotation_name_ne), annotation_name := annotation_name_ne]
  out[, annotation_name_ne := NULL]
  out <- merge(
    out,
    limbic[, .(annotation = as.integer(annotation), in_prespecified_set)],
    by = "annotation",
    all.x = TRUE
  )
  out[is.na(in_prespecified_set), in_prespecified_set := FALSE]
  out[, baseline_vs_intruder_class_agree :=
    twodg_baseline_symmetry_class == twodg_intruder_symmetry_class &
      !is.na(twodg_baseline_symmetry_class) &
      !is.na(twodg_intruder_symmetry_class) &
      twodg_baseline_symmetry_class != "indeterminate"]
  out[, intruder_vs_ne_class_agree :=
    twodg_intruder_symmetry_class == ne_symmetry_class &
      !is.na(twodg_intruder_symmetry_class) &
      !is.na(ne_symmetry_class) &
      twodg_intruder_symmetry_class != "indeterminate"]
  out
}

plot_class_counts <- function(symm, path) {
  long <- rbindlist(list(
    symm[, .(class = twodg_baseline_symmetry_class, map = "Baseline 2DG")],
    symm[, .(class = twodg_intruder_symmetry_class, map = "Intruder 2DG")],
    symm[, .(class = ne_symmetry_class, map = "NEUROeSTIMator")]
  ))
  long <- long[!is.na(class)]
  counts <- long[, .N, by = .(map, class)]
  counts[, class := factor(
    class,
    levels = c("symmetric", "unilateral", "opposing", "indeterminate")
  )]
  p <- ggplot(counts, aes(x = class, y = N, fill = class)) +
    geom_col(width = 0.72) +
    facet_wrap(~map, nrow = 1) +
    scale_fill_manual(
      values = c(
        symmetric = "#4C78A8",
        unilateral = "#F58518",
        opposing = "#E45756",
        indeterminate = "#B0B0B0"
      )
    ) +
    labs(
      title = "Hemisphere symmetry classes by map",
      subtitle = paste0(
        "2DG τ=", TAU_2DG, " / τ0=", TAU0_2DG,
        "; NE τ=", TAU_NE, " / τ0=", TAU0_NE,
        "; δ=", DELTA_REL, "; midline±", MIDLINE_BAND_UM, " µm excluded from L/R"
      ),
      x = NULL,
      y = "Allen regions",
      fill = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "#f2f2f2", colour = NA),
      axis.text.x = element_text(angle = 25, hjust = 1)
    )
  ggsave(path, p, width = 11, height = 4.2, dpi = 200, bg = "white")
}

plot_lat_concordance <- function(symm, path) {
  d <- symm[is.finite(intruder_lat_index) & is.finite(ne_lat_index) &
    is.finite(baseline_lat_index)]
  d[, set_label := fifelse(in_prespecified_set, "Limbic/hypothalamic set", "Other")]
  p1 <- ggplot(d, aes(intruder_lat_index, ne_lat_index, colour = set_label)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_point(alpha = 0.75, size = 2.2) +
    geom_smooth(method = "lm", se = TRUE, colour = "#222222", fill = "#cccccc", linewidth = 0.5) +
    scale_colour_manual(values = c("Limbic/hypothalamic set" = "#FF6C0C", Other = "#4C78A8")) +
    labs(
      title = "Lateralization concordance",
      subtitle = "Index = (R − L) / (|R| + |L|); positive = right-biased WildR−SPF effect",
      x = "Intruder 2DG lat. index",
      y = "NEUROeSTIMator lat. index",
      colour = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  p2 <- ggplot(d, aes(baseline_lat_index, ne_lat_index, colour = set_label)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_point(alpha = 0.75, size = 2.2) +
    geom_smooth(method = "lm", se = TRUE, colour = "#222222", fill = "#cccccc", linewidth = 0.5) +
    scale_colour_manual(values = c("Limbic/hypothalamic set" = "#FF6C0C", Other = "#4C78A8")) +
    labs(
      x = "Baseline 2DG lat. index",
      y = "NEUROeSTIMator lat. index",
      colour = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  # simple side-by-side via cowplot if available, else save two pages as one with patchwork-free grid
  png(path, width = 1400, height = 620, res = 140)
  print(p1, vp = grid::viewport(x = 0.25, width = 0.48))
  print(p2, vp = grid::viewport(x = 0.75, width = 0.48))
  dev.off()
}

plot_class_strip <- function(symm, path) {
  keep <- symm[
    twodg_intruder_symmetry_class %in% c("symmetric", "unilateral", "opposing") |
      ne_symmetry_class %in% c("symmetric", "unilateral", "opposing")
  ]
  if (!nrow(keep)) {
    keep <- symm[1:min(30, nrow(symm))]
  }
  long <- rbindlist(list(
    keep[, .(annotation_name, map = "Baseline", class = twodg_baseline_symmetry_class, in_prespecified_set)],
    keep[, .(annotation_name, map = "Intruder", class = twodg_intruder_symmetry_class, in_prespecified_set)],
    keep[, .(annotation_name, map = "NE", class = ne_symmetry_class, in_prespecified_set)]
  ))
  long <- long[!is.na(class)]
  # Limit to readable set: limbic + non-indeterminate
  focus <- unique(long[
    in_prespecified_set == TRUE | class != "indeterminate",
    annotation_name
  ])
  long <- long[annotation_name %in% focus]
  if (length(focus) > 40) {
    focus <- unique(long[in_prespecified_set == TRUE | class %in% c("opposing", "unilateral", "symmetric"), annotation_name])
    long <- long[annotation_name %in% head(focus, 40)]
  }
  long[, class := factor(class, levels = c("symmetric", "unilateral", "opposing", "indeterminate"))]
  long[, map := factor(map, levels = c("Baseline", "Intruder", "NE"))]
  p <- ggplot(long, aes(map, annotation_name, fill = class)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_manual(
      values = c(
        symmetric = "#4C78A8",
        unilateral = "#F58518",
        opposing = "#E45756",
        indeterminate = "#E8E8E8"
      )
    ) +
    labs(
      title = "Per-region hemisphere signature strip",
      subtitle = "Selected limbic/non-indeterminate regions; full table on disk",
      x = NULL,
      y = NULL,
      fill = "Class"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 7),
      legend.position = "bottom"
    )
  ggsave(path, p, width = 8.5, height = max(5, 0.18 * uniqueN(long$annotation_name) + 2), dpi = 200, bg = "white")
}

main <- function() {
  outputs <- c(
    SEQ_HEMI_RDS, SEQ_HEMI_CSV, TWODG_BILAT, NE_BILAT, SYMM_TABLE, SYMM_CONFIG,
    FIG_CLASSES, FIG_LAT, FIG_COUNTS
  )
  stop_if_exists(outputs)
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

  cfg <- sprintf(
    paste0(
      "{\n",
      "  \"analysis_date\": \"%s\",\n",
      "  \"ml_midline_um\": %d,\n",
      "  \"midline_band_um\": %d,\n",
      "  \"tau_2dg\": %s,\n",
      "  \"tau0_2dg\": %s,\n",
      "  \"tau_ne\": %s,\n",
      "  \"tau0_ne\": %s,\n",
      "  \"delta_rel\": %s,\n",
      "  \"min_2dg_voxels_per_side\": %d,\n",
      "  \"min_ne_cells_per_side\": %d,\n",
      "  \"notes\": \"Thresholds locked before class inspection per hemisphere plan. NE deltas are descriptive section/gate-level, not animal-level.\"\n",
      "}\n"
    ),
    ANALYSIS_DATE, ML_MIDLINE_UM, MIDLINE_BAND_UM,
    TAU_2DG, TAU0_2DG, TAU_NE, TAU0_NE, DELTA_REL,
    MIN_2DG_VOXELS, MIN_NE_CELLS_SIDE
  )
  writeLines(cfg, SYMM_CONFIG)

  log_message("Loading seqFISH join + applying anatomical ML")
  joined <- as.data.table(readRDS(JOINED))
  cells <- joined[
    registered_ccf == TRUE & is.finite(ccfx) & is.finite(ccfz)
  ][
    ,
    .(
      cell_id,
      predicted_activity,
      group = as.character(group),
      final_gate = as.character(final_gate),
      sample_id = as.character(sample_id),
      ccfx,
      ccfy,
      ccfz,
      annotation = as.integer(annotation)
    )
  ]
  cells <- apply_anatomical_ml(cells)
  names_map <- unique(fread(INTRUDER_OVERLAY, select = c("annotation", "annotation_name")))
  names_map[, annotation := as.integer(annotation)]
  cells <- merge(cells, names_map, by = "annotation", all.x = TRUE)
  saveRDS(cells, SEQ_HEMI_RDS, compress = "xz")
  fwrite(cells, SEQ_HEMI_CSV)

  log_message(
    "seqFISH anatomical hemispheres: ",
    paste(names(table(cells$anatomical_hemisphere)), table(cells$anatomical_hemisphere), sep = "=", collapse = ", ")
  )

  log_message("2DG bilateral region summaries")
  twodg <- rbindlist(list(
    summarize_2dg_map(INTRUDER_OVERLAY, "intruder_wildr_vs_spf"),
    summarize_2dg_map(BASELINE_OVERLAY, "baseline_wildr_vs_spf"),
    summarize_2dg_map(WILDR_IB_OVERLAY, "wildr_intruder_vs_baseline"),
    summarize_2dg_map(SPF_IB_OVERLAY, "spf_intruder_vs_baseline")
  ), fill = TRUE)
  fwrite(twodg, TWODG_BILAT)

  log_message("NE bilateral / available-side summaries")
  ne <- summarize_ne(cells)
  ne_avail_path <- file.path(
    TABLE_DIR, paste0("SeqFISH_NE_region_available_side_", ANALYSIS_DATE, ".csv")
  )
  if (file.exists(ne_avail_path)) {
    stop("Refusing to overwrite:\n", ne_avail_path, call. = FALSE)
  }
  fwrite(ne$bilateral, NE_BILAT)
  fwrite(ne$available_side, ne_avail_path)

  limbic <- if (file.exists(LIMBIC_SET)) {
    fread(LIMBIC_SET)
  } else {
    data.table(annotation = integer(), in_prespecified_set = logical())
  }

  log_message("Cross-modality symmetry concordance")
  symm <- build_concordance(twodg, ne$bilateral, limbic)
  fwrite(symm, SYMM_TABLE)

  log_message("Figures")
  plot_class_counts(symm, FIG_COUNTS)
  plot_class_strip(symm, FIG_CLASSES)
  if (requireNamespace("grid", quietly = TRUE)) {
    tryCatch(plot_lat_concordance(symm, FIG_LAT), error = function(e) {
      log_message("lat concordance plot failed: ", conditionMessage(e))
      # fallback single panel
      d <- symm[is.finite(intruder_lat_index) & is.finite(ne_lat_index)]
      p <- ggplot(d, aes(intruder_lat_index, ne_lat_index)) +
        geom_point(alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE) +
        theme_bw() +
        labs(title = "Intruder vs NE lateralization index", x = "Intruder lat.", y = "NE lat.")
      ggsave(FIG_LAT, p, width = 6.5, height = 5, dpi = 200, bg = "white")
    })
  }

  # Console summary
  cat("\nClass counts (Intruder 2DG):\n")
  print(symm[, .N, by = twodg_intruder_symmetry_class])
  cat("\nClass counts (Baseline 2DG):\n")
  print(symm[, .N, by = twodg_baseline_symmetry_class])
  cat("\nClass counts (NE):\n")
  print(symm[, .N, by = ne_symmetry_class])
  ok <- symm[is.finite(intruder_lat_index) & is.finite(ne_lat_index)]
  if (nrow(ok) > 5) {
    cat(
      "\nSpearman lat(Intruder) vs lat(NE): ",
      cor(ok$intruder_lat_index, ok$ne_lat_index, method = "spearman"),
      "\nSpearman lat(Baseline) vs lat(NE): ",
      cor(ok$baseline_lat_index, ok$ne_lat_index, method = "spearman", use = "complete.obs"),
      "\n",
      sep = ""
    )
  }
  log_message("Done: ", ANALYSIS_DATE)
}

main()
