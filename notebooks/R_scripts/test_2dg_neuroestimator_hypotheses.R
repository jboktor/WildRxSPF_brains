#!/usr/bin/env Rscript

# Pre-specified H1–H5 battery for 2DG × seqFISH CCF co-embedding.
# Group-map / section-level only — not animal-level inference.
#
# /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript \
#   notebooks/R_scripts/test_2dg_neuroestimator_hypotheses.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

WKDIR <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
SOURCE_DATE <- "2026-11-11"
ANALYSIS_DATE <- "2026-11-12"
N_PERM <- 10000L
ALPHA_H3 <- 0.10
SET_SEED <- 42L

TABLE_DIR <- file.path(WKDIR, "tables/2DG_CCF_Integration")
FIGURE_DIR <- file.path(WKDIR, "figures/2DG_CCF_Integration")
SOURCE_DIR <- file.path(WKDIR, "data/interim/2DG_CCF_Integration", SOURCE_DATE)
OUT_DIR <- file.path(WKDIR, "data/interim/2DG_CCF_Integration", ANALYSIS_DATE)

REGISTRATION_PATH <- file.path(
  SOURCE_DIR, paste0("2dg_registration_summary_strict_", SOURCE_DATE, ".csv")
)
ANNOTATION_PATH <- file.path(
  TABLE_DIR, paste0("Allen_CCF_annotation_names_hierarchy_", SOURCE_DATE, ".csv")
)
INTRUDER_REGION_PATH <- file.path(
  TABLE_DIR, paste0("2DG_CCF_region_summary_strict_", SOURCE_DATE, ".csv")
)
INTRUDER_OVERLAY_PATH <- file.path(
  SOURCE_DIR, paste0("2dg_overlay_points_ccf_scalar-t_", SOURCE_DATE, ".csv.gz")
)
BASELINE_REGION_PATH <- file.path(
  TABLE_DIR, paste0("2DG_CCF_region_summary_baseline_wildr_vs_spf_", ANALYSIS_DATE, ".csv")
)
WILDR_IB_REGION_PATH <- file.path(
  TABLE_DIR, paste0("2DG_CCF_region_summary_wildr_intruder_vs_baseline_", ANALYSIS_DATE, ".csv")
)
SPF_IB_REGION_PATH <- file.path(
  TABLE_DIR, paste0("2DG_CCF_region_summary_spf_intruder_vs_baseline_", ANALYSIS_DATE, ".csv")
)
GATE_BRIDGE_PATH <- file.path(
  TABLE_DIR, paste0("2DG_SeqFISH_annotation_gate_bridge_scalar-t_", SOURCE_DATE, ".csv.gz")
)
GATED_EFFECT_PATH <- file.path(
  TABLE_DIR, paste0("2DG_NEUROeSTIMator_gated_effects_scalar-t_", SOURCE_DATE, ".csv")
)
JOINED_PATH <- file.path(
  SOURCE_DIR, paste0("seqFISH_neuroestimator_ccf_strict_", SOURCE_DATE, ".rds")
)
ROI_PATH <- file.path(
  TABLE_DIR, paste0("2DG_ROI_seed_ccf_overlap_", ANALYSIS_DATE, ".csv")
)

BATTERY_PATH <- file.path(
  TABLE_DIR, paste0("2DG_hypothesis_battery_controlled_", ANALYSIS_DATE, ".csv")
)
LIMBIC_SET_PATH <- file.path(
  TABLE_DIR, paste0("2DG_prespecified_limbic_hypothalamic_annotations_", ANALYSIS_DATE, ".csv")
)
GATE_CONTROL_PATH <- file.path(
  TABLE_DIR, paste0("2DG_gate_control_triad_controlled_", ANALYSIS_DATE, ".csv")
)
H2_ENRICH_PATH <- file.path(
  TABLE_DIR, paste0("2DG_H2_region_set_enrichment_controlled_", ANALYSIS_DATE, ".csv")
)
BASELINE_OVERLAY_PATH <- file.path(
  OUT_DIR, paste0("2dg_overlay_points_ccf_baseline_wildr_vs_spf_", ANALYSIS_DATE, ".csv.gz")
)
WILDR_IB_OVERLAY_PATH <- file.path(
  OUT_DIR, paste0("2dg_overlay_points_ccf_wildr_intruder_vs_baseline_", ANALYSIS_DATE, ".csv.gz")
)
SPF_IB_OVERLAY_PATH <- file.path(
  OUT_DIR, paste0("2dg_overlay_points_ccf_spf_intruder_vs_baseline_", ANALYSIS_DATE, ".csv.gz")
)

log_message <- function(...) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
}

stop_if_exists <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing)) {
    stop("Refusing to overwrite:\n", paste(existing, collapse = "\n"), call. = FALSE)
  }
}

require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing:\n", paste(missing, collapse = "\n"), call. = FALSE)
  }
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

weighted_mean_safe <- function(x, weight) {
  keep <- is.finite(x) & is.finite(weight) & weight > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * weight[keep]) / sum(weight[keep])
}

spearman_perm_p <- function(x, y, n_perm = N_PERM) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 5L || uniqueN(x) < 3L || uniqueN(y) < 3L) {
    return(list(rho = NA_real_, p = NA_real_, n = length(x)))
  }
  rho <- suppressWarnings(cor(x, y, method = "spearman"))
  set.seed(SET_SEED)
  null <- replicate(n_perm, suppressWarnings(cor(sample(x), y, method = "spearman")))
  p <- (1 + sum(abs(null) >= abs(rho), na.rm = TRUE)) / (1 + sum(is.finite(null)))
  list(rho = rho, p = p, n = length(x), null_mean = mean(null, na.rm = TRUE))
}

allocate_positive_mean_t_from_overlay <- function(overlay_path, bridge) {
  overlay <- fread(overlay_path, select = c("annotation", "t_value"))
  overlay[, annotation := as.integer(annotation)]
  merged <- merge(
    overlay,
    bridge[, .(annotation, final_gate, p_gate_given_annotation)],
    by = "annotation",
    allow.cartesian = TRUE
  )
  merged[, positive := t_value > 0]
  merged[
    ,
    .(
      two_dg_positive_mean_t = weighted_mean_safe(
        t_value[positive],
        p_gate_given_annotation[positive]
      ),
      two_dg_mean_t = weighted_mean_safe(t_value, p_gate_given_annotation),
      allocated_overlay_pixels = sum(p_gate_given_annotation)
    ),
    by = final_gate
  ]
}

build_limbic_set <- function(annotation_names) {
  pattern <- paste0(
    "(?i)amygdala|amygdalar|hypothal|bed nuclei of the stria|lateral septal|",
    "medial septal|nucleus accumbens|paraventricular hypothalamic|",
    "medial preoptic|infralimbic|prelimbic|anterior cingulate"
  )
  set <- annotation_names[grepl(pattern, annotation_name)]
  set[, in_prespecified_set := TRUE]
  set
}

test_h1 <- function(registration) {
  accepted <- registration[accepted_final == TRUE]
  continuity_fail <- registration[
    accepted_slice_fit == TRUE & accepted_final == FALSE
  ]
  blank <- registration[rejection_reason == "blank_or_empty_template"]
  estimate <- sprintf(
    "accepted_final=%d/%d; mean_dice=%.3f; continuity_rejects=%d; blank=%d",
    nrow(accepted),
    nrow(registration),
    mean(accepted$dice),
    nrow(continuity_fail),
    nrow(blank)
  )
  decision <- if (nrow(accepted) >= 40L && mean(accepted$dice) >= 0.9) {
    "supported"
  } else if (nrow(accepted) >= 20L) {
    "inconclusive"
  } else {
    "not_supported"
  }
  data.table(
    hypothesis_id = "H1",
    metric = "accepted_slice_dice_and_continuity_qc",
    estimate = estimate,
    null_summary = "QC claim: continuity gate retained; no chance Dice model",
    p_perm = NA_real_,
    fdr_or_predeclared_alpha = NA_real_,
    control_result = sprintf(
      "%d high-Dice fits rejected only by AP continuity",
      nrow(continuity_fail)
    ),
    decision = decision
  )
}

test_h2 <- function(intruder_regions, baseline_regions, limbic_set) {
  limbic_ids <- limbic_set$annotation
  score_map <- function(regions, label) {
    dt <- copy(regions)
    dt[, annotation := as.integer(annotation)]
    dt[, in_set := annotation %in% limbic_ids]
    # Region-level mean absolute t among regions with any overlay mass.
    dt <- dt[is.finite(mean_absolute_t) | is.finite(mean_t_value)]
    if (!"mean_absolute_t" %in% names(dt) || all(!is.finite(dt$mean_absolute_t))) {
      dt[, mean_absolute_t := abs(mean_t_value)]
    }
    observed <- safe_mean(dt[in_set == TRUE, mean_absolute_t]) -
      safe_mean(dt[in_set == FALSE, mean_absolute_t])
    set.seed(SET_SEED)
    null <- replicate(N_PERM, {
      shuffled <- sample(dt$in_set)
      safe_mean(dt$mean_absolute_t[shuffled]) - safe_mean(dt$mean_absolute_t[!shuffled])
    })
    p <- (1 + sum(null >= observed, na.rm = TRUE)) / (1 + sum(is.finite(null)))
    list(
      label = label,
      observed = observed,
      p = p,
      null_mean = mean(null, na.rm = TRUE),
      n_in = sum(dt$in_set),
      n_out = sum(!dt$in_set)
    )
  }
  intr <- score_map(intruder_regions, "intruder")
  base <- score_map(baseline_regions, "baseline")
  # Supported if Intruder enrichment significant AND stronger than Baseline.
  decision <- if (is.finite(intr$p) && intr$p < 0.05 && intr$observed > base$observed) {
    "supported"
  } else if (is.finite(intr$p) && intr$p < 0.05) {
    "inconclusive"
  } else {
    "not_supported"
  }
  enrich <- data.table(
    map_family = c(intr$label, base$label),
    delta_mean_abs_t = c(intr$observed, base$observed),
    p_perm = c(intr$p, base$p),
    null_mean = c(intr$null_mean, base$null_mean),
    n_in_set = c(intr$n_in, base$n_out * 0 + intr$n_in),
    n_out_set = c(intr$n_out, base$n_out)
  )
  enrich[map_family == "baseline", `:=`(n_in_set = base$n_in, n_out_set = base$n_out)]
  fwrite(enrich, H2_ENRICH_PATH)
  data.table(
    hypothesis_id = "H2",
    metric = "limbic_hypothalamic_mean_abs_t_minus_rest",
    estimate = sprintf(
      "intruder_delta=%.4f (p=%.4g); baseline_delta=%.4f (p=%.4g)",
      intr$observed, intr$p, base$observed, base$p
    ),
    null_summary = "annotation-label permutation within map (10k)",
    p_perm = intr$p,
    fdr_or_predeclared_alpha = NA_real_,
    control_result = sprintf(
      "baseline_delta=%.4f; intruder_stronger_than_baseline=%s",
      base$observed,
      intr$observed > base$observed
    ),
    decision = decision
  )
}

test_h3 <- function(gate_triad) {
  intr <- gate_triad[!is.na(neuro_hedges_g) & !is.na(intruder_positive_mean_t)]
  base <- gate_triad[!is.na(neuro_hedges_g) & !is.na(baseline_positive_mean_t)]
  intr_stat <- spearman_perm_p(intr$intruder_positive_mean_t, intr$neuro_hedges_g)
  base_stat <- spearman_perm_p(base$baseline_positive_mean_t, base$neuro_hedges_g)
  decision <- if (is.finite(intr_stat$p) && intr_stat$p < ALPHA_H3) {
    if (!is.finite(base_stat$rho) || abs(intr_stat$rho) > abs(base_stat$rho)) {
      "supported"
    } else {
      "inconclusive"
    }
  } else {
    "not_supported"
  }
  data.table(
    hypothesis_id = "H3",
    metric = "spearman_intruder_positive_mean_t_vs_neuro_hedges_g",
    estimate = sprintf(
      "rho=%.3f; n_gates=%d; baseline_rho=%.3f",
      intr_stat$rho, intr_stat$n, base_stat$rho
    ),
    null_summary = "gate-label permutation of 2DG scores (10k)",
    p_perm = intr_stat$p,
    fdr_or_predeclared_alpha = ALPHA_H3,
    control_result = sprintf(
      "baseline_rho=%.3f; baseline_p=%.4g; collapses_or_weaker=%s",
      base_stat$rho,
      base_stat$p,
      is.finite(base_stat$rho) && abs(base_stat$rho) >= abs(intr_stat$rho)
    ),
    decision = decision
  )
}

test_h4 <- function(gate_triad) {
  dt <- gate_triad[
    is.finite(neuro_hedges_g) &
      is.finite(intruder_positive_mean_t) &
      is.finite(wildr_ib_positive_mean_t) &
      is.finite(spf_ib_positive_mean_t)
  ]
  intr <- spearman_perm_p(dt$intruder_positive_mean_t, dt$neuro_hedges_g)
  wildr <- spearman_perm_p(dt$wildr_ib_positive_mean_t, dt$neuro_hedges_g)
  spf <- spearman_perm_p(dt$spf_ib_positive_mean_t, dt$neuro_hedges_g)
  control_max <- max(abs(c(wildr$rho, spf$rho)), na.rm = TRUE)
  # Bootstrap gates for confidence on |rho_intr| - control_max
  set.seed(SET_SEED)
  boot <- replicate(2000, {
    idx <- sample.int(nrow(dt), replace = TRUE)
    r_i <- suppressWarnings(cor(
      dt$intruder_positive_mean_t[idx],
      dt$neuro_hedges_g[idx],
      method = "spearman"
    ))
    r_w <- suppressWarnings(cor(
      dt$wildr_ib_positive_mean_t[idx],
      dt$neuro_hedges_g[idx],
      method = "spearman"
    ))
    r_s <- suppressWarnings(cor(
      dt$spf_ib_positive_mean_t[idx],
      dt$neuro_hedges_g[idx],
      method = "spearman"
    ))
    abs(r_i) - max(abs(c(r_w, r_s)), na.rm = TRUE)
  })
  p_boot <- mean(boot <= 0, na.rm = TRUE)
  decision <- if (is.finite(intr$rho) && abs(intr$rho) > control_max && p_boot < 0.5) {
    "supported"
  } else if (control_max >= abs(intr$rho)) {
    "not_supported"
  } else {
    "inconclusive"
  }
  data.table(
    hypothesis_id = "H4",
    metric = "abs_rho_intruder_minus_max_abs_rho_IB",
    estimate = sprintf(
      "|rho_intr|=%.3f; |rho_wildr_IB|=%.3f; |rho_spf_IB|=%.3f; delta=%.3f",
      abs(intr$rho), abs(wildr$rho), abs(spf$rho), abs(intr$rho) - control_max
    ),
    null_summary = "bootstrap over gates (2000): P(delta<=0)",
    p_perm = p_boot,
    fdr_or_predeclared_alpha = NA_real_,
    control_result = sprintf(
      "max_|rho|_IB=%.3f; IB_explains_as_well_or_better=%s",
      control_max,
      control_max >= abs(intr$rho)
    ),
    decision = decision
  )
}

test_h5 <- function(roi) {
  usable <- roi[slice_accepted == TRUE & is.finite(iou_target)]
  if (!nrow(usable)) {
    return(data.table(
      hypothesis_id = "H5",
      metric = "roi_seed_iou_vs_target_parcels",
      estimate = "no accepted seed slices available",
      null_summary = "slice-level anatomical correspondence QC",
      p_perm = NA_real_,
      fdr_or_predeclared_alpha = NA_real_,
      control_result = "several seed slices failed continuity QC",
      decision = "inconclusive"
    ))
  }
  mean_iou <- mean(usable$iou_target)
  n_pos <- sum(usable$n_target_pixels > 0)
  decision <- if (mean_iou >= 0.05 && n_pos >= 3L) {
    "supported"
  } else if (n_pos >= 1L) {
    "inconclusive"
  } else {
    "not_supported"
  }
  data.table(
    hypothesis_id = "H5",
    metric = "roi_seed_iou_vs_target_parcels",
    estimate = sprintf(
      "mean_iou=%.4f across %d usable seeds; n_with_target_overlap=%d/%d",
      mean_iou, nrow(usable), n_pos, nrow(usable)
    ),
    null_summary = "QC threshold: mean IoU>=0.05 and >=3 ROIs with overlap",
    p_perm = NA_real_,
    fdr_or_predeclared_alpha = NA_real_,
    control_result = paste(
      usable[, sprintf("%s:iou=%.3f", roi, iou_target)],
      collapse = "; "
    ),
    decision = decision
  )
}

plot_decision_strip <- function(battery) {
  path <- file.path(
    FIGURE_DIR, paste0("2DG_hypothesis_decision_strip_controlled_", ANALYSIS_DATE, ".png")
  )
  if (file.exists(path)) return(invisible(path))
  plot_data <- copy(battery)
  plot_data[, hypothesis_id := factor(hypothesis_id, levels = rev(unique(hypothesis_id)))]
  plot_data[, decision := factor(
    decision,
    levels = c("supported", "inconclusive", "not_supported")
  )]
  pal <- c(
    supported = "#1B9E77",
    inconclusive = "#D95F02",
    not_supported = "#7570B3"
  )
  p <- ggplot(plot_data, aes(x = 1, y = hypothesis_id, fill = decision)) +
    geom_tile(width = 0.7, height = 0.7) +
    geom_text(aes(label = decision), color = "white", size = 3.6, fontface = "bold") +
    scale_fill_manual(values = pal, drop = FALSE) +
    labs(
      title = "Pre-specified hypothesis decisions",
      subtitle = "Only H1–H5 are promoted as findings",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  ggsave(path, p, width = 7.5, height = 4.2, dpi = 220)
  invisible(path)
}

plot_h2_enrichment <- function() {
  path <- file.path(
    FIGURE_DIR, paste0("2DG_H2_region_set_enrichment_controlled_", ANALYSIS_DATE, ".png")
  )
  if (file.exists(path) || !file.exists(H2_ENRICH_PATH)) return(invisible(NULL))
  enrich <- fread(H2_ENRICH_PATH)
  p <- ggplot(enrich, aes(x = map_family, y = delta_mean_abs_t, fill = map_family)) +
    geom_col(width = 0.55) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50") +
    geom_text(
      aes(label = sprintf("p=%.3g", p_perm)),
      vjust = ifelse(enrich$delta_mean_abs_t >= 0, -0.4, 1.3),
      size = 3.3
    ) +
    scale_fill_manual(values = c(intruder = "#FF6C0C", baseline = "#4C78A8")) +
    labs(
      title = "H2: limbic/hypothalamic |t| enrichment vs rest of brain",
      subtitle = "Positive = higher mean |t| inside the pre-specified annotation set",
      x = NULL,
      y = "Δ mean |t| (set − complement)"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  ggsave(path, p, width = 7.2, height = 4.5, dpi = 220)
}

plot_control_triad <- function(gate_triad) {
  path <- file.path(
    FIGURE_DIR, paste0("2DG_gate_control_triad_controlled_", ANALYSIS_DATE, ".png")
  )
  if (file.exists(path)) return(invisible(path))
  long <- melt(
    gate_triad,
    id.vars = c("final_gate", "neuro_hedges_g"),
    measure.vars = c(
      "intruder_positive_mean_t",
      "baseline_positive_mean_t",
      "wildr_ib_positive_mean_t",
      "spf_ib_positive_mean_t"
    ),
    variable.name = "map_metric",
    value.name = "two_dg_score"
  )
  long[, map_metric := factor(
    map_metric,
    levels = c(
      "intruder_positive_mean_t",
      "baseline_positive_mean_t",
      "wildr_ib_positive_mean_t",
      "spf_ib_positive_mean_t"
    ),
    labels = c(
      "WildR vs SPF Intruder",
      "WildR vs SPF Baseline",
      "WildR Intruder−Baseline",
      "SPF Intruder−Baseline"
    )
  )]
  p <- ggplot(long, aes(x = two_dg_score, y = neuro_hedges_g)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_point(color = "#FF6C0C", size = 2.4, alpha = 0.85) +
    geom_text(aes(label = final_gate), size = 2.3, vjust = -0.7, check_overlap = TRUE) +
    facet_wrap(~map_metric, scales = "free_x") +
    labs(
      title = "Control triad: gate-allocated positive mean t vs NEUROeSTIMator g",
      x = "Gate-allocated positive mean t",
      y = "NEUROeSTIMator Hedges' g (WildR − SPF)"
    ) +
    theme_minimal(base_size = 11)
  ggsave(path, p, width = 10.5, height = 7.2, dpi = 220)
}

plot_registration_qc <- function(registration) {
  path <- file.path(
    FIGURE_DIR, paste0("2DG_CCF_registration_qc_status_", ANALYSIS_DATE, ".png")
  )
  if (file.exists(path)) return(invisible(path))
  dt <- copy(registration)
  dt[, status := fifelse(
    accepted_final == TRUE,
    "accepted_final",
    fifelse(accepted_slice_fit == TRUE, "continuity_rejected", "fit_failed")
  )]
  pal <- c(
    accepted_final = "#1B9E77",
    continuity_rejected = "#D95F02",
    fit_failed = "#999999"
  )
  p <- ggplot(dt, aes(x = source_number, y = dice, color = status)) +
    geom_hline(yintercept = 0.45, linetype = "dashed", color = "grey50") +
    geom_point(size = 2.2) +
    scale_color_manual(values = pal) +
    labs(
      title = "H1 registration QC: Template → Allen CCF",
      subtitle = "Continuity rejects keep high Dice without inventing an AP trajectory",
      x = "Source slice number",
      y = "Tissue Dice",
      color = NULL
    ) +
    theme_minimal(base_size = 12)
  ggsave(path, p, width = 9.5, height = 4.4, dpi = 220)
}

main <- function() {
  require_files(c(
    REGISTRATION_PATH,
    ANNOTATION_PATH,
    INTRUDER_REGION_PATH,
    BASELINE_REGION_PATH,
    INTRUDER_OVERLAY_PATH,
    BASELINE_OVERLAY_PATH,
    WILDR_IB_OVERLAY_PATH,
    SPF_IB_OVERLAY_PATH,
    GATE_BRIDGE_PATH,
    GATED_EFFECT_PATH,
    ROI_PATH
  ))
  outputs <- c(
    BATTERY_PATH,
    GATE_CONTROL_PATH,
    H2_ENRICH_PATH,
    file.path(FIGURE_DIR, paste0("2DG_hypothesis_decision_strip_controlled_", ANALYSIS_DATE, ".png")),
    file.path(FIGURE_DIR, paste0("2DG_H2_region_set_enrichment_controlled_", ANALYSIS_DATE, ".png")),
    file.path(FIGURE_DIR, paste0("2DG_gate_control_triad_controlled_", ANALYSIS_DATE, ".png")),
    file.path(FIGURE_DIR, paste0("2DG_CCF_registration_qc_status_", ANALYSIS_DATE, ".png"))
  )
  # Limbic set may already exist from the first battery pass.
  if (!file.exists(LIMBIC_SET_PATH)) {
    outputs <- c(outputs, LIMBIC_SET_PATH)
  }
  stop_if_exists(outputs)
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

  registration <- fread(REGISTRATION_PATH)
  annotations <- fread(ANNOTATION_PATH)
  annotations[, annotation := as.integer(annotation)]
  if (file.exists(LIMBIC_SET_PATH)) {
    limbic_set <- fread(LIMBIC_SET_PATH)
  } else {
    limbic_set <- build_limbic_set(annotations)
    fwrite(limbic_set, LIMBIC_SET_PATH)
  }

  intruder_regions <- fread(INTRUDER_REGION_PATH)
  intruder_regions[, annotation := as.integer(annotation)]
  # H2 uses |t| on significant overlay voxels, not tissue-wide near-zero means.
  intruder_overlay_abs <- fread(
    INTRUDER_OVERLAY_PATH,
    select = c("annotation", "t_value")
  )[
    ,
    .(mean_absolute_t = mean(abs(t_value))),
    by = .(annotation = as.integer(annotation))
  ]
  intruder_regions[, mean_absolute_t := NULL]
  intruder_regions <- merge(
    intruder_regions,
    intruder_overlay_abs,
    by = "annotation",
    all.x = TRUE
  )
  baseline_regions <- fread(BASELINE_REGION_PATH)
  bridge <- fread(GATE_BRIDGE_PATH)
  bridge[, annotation := as.integer(annotation)]
  neuro <- fread(GATED_EFFECT_PATH)[
    ,
    .(final_gate, neuro_hedges_g, neuro_raw_difference)
  ]

  # Pixel-level allocation matches the validated Intruder gated-effect metric.
  intr_gate <- allocate_positive_mean_t_from_overlay(INTRUDER_OVERLAY_PATH, bridge)
  setnames(
    intr_gate,
    c("two_dg_positive_mean_t", "two_dg_mean_t", "allocated_overlay_pixels"),
    c("intruder_positive_mean_t", "intruder_mean_t", "intruder_allocated_overlay")
  )
  base_gate <- allocate_positive_mean_t_from_overlay(BASELINE_OVERLAY_PATH, bridge)
  setnames(
    base_gate,
    c("two_dg_positive_mean_t", "two_dg_mean_t", "allocated_overlay_pixels"),
    c("baseline_positive_mean_t", "baseline_mean_t", "baseline_allocated_overlay")
  )
  wildr_gate <- allocate_positive_mean_t_from_overlay(WILDR_IB_OVERLAY_PATH, bridge)
  setnames(
    wildr_gate,
    c("two_dg_positive_mean_t", "two_dg_mean_t", "allocated_overlay_pixels"),
    c("wildr_ib_positive_mean_t", "wildr_ib_mean_t", "wildr_ib_allocated_overlay")
  )
  spf_gate <- allocate_positive_mean_t_from_overlay(SPF_IB_OVERLAY_PATH, bridge)
  setnames(
    spf_gate,
    c("two_dg_positive_mean_t", "two_dg_mean_t", "allocated_overlay_pixels"),
    c("spf_ib_positive_mean_t", "spf_ib_mean_t", "spf_ib_allocated_overlay")
  )

  gate_triad <- Reduce(
    function(x, y) merge(x, y, by = "final_gate", all = TRUE),
    list(neuro, intr_gate, base_gate, wildr_gate, spf_gate)
  )
  fwrite(gate_triad, GATE_CONTROL_PATH)

  roi <- fread(ROI_PATH)
  battery <- rbindlist(list(
    test_h1(registration),
    test_h2(intruder_regions, baseline_regions, limbic_set),
    test_h3(gate_triad),
    test_h4(gate_triad),
    test_h5(roi)
  ), fill = TRUE)
  fwrite(battery, BATTERY_PATH)

  plot_decision_strip(battery)
  plot_h2_enrichment()
  plot_control_triad(gate_triad)
  plot_registration_qc(registration)

  log_message("Hypothesis battery: ", BATTERY_PATH)
  print(battery[, .(hypothesis_id, decision, estimate)])
}

main()
