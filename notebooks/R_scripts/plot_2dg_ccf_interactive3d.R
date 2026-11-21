#!/usr/bin/env Rscript

# Upgraded interactive 2DG × seqFISH CCF viewers:
#   1) Tall multi-map 3D explorer (Baseline / Intruder / NE / CCF / AP overlap)
#   2) Animated contrast transition (Baseline → Intruder → NE)
#   3) Per-region NE vs Baseline/Intruder co-variation (interactive + static)
#
# Polarity (all WildR−SPF signed maps):
#   positive / blue  = WildR higher
#   negative / red   = SPF higher
#
# /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript \
#   notebooks/R_scripts/plot_2dg_ccf_interactive3d.R

suppressPackageStartupMessages({
  library(data.table)
  library(plotly)
  library(htmlwidgets)
  library(htmltools)
  library(ggplot2)
})

WKDIR <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
SOURCE_DATE <- "2026-11-11"
CONTROL_DATE <- "2026-11-12"
ANALYSIS_DATE <- "2026-11-15"

SOURCE_DIR <- file.path(WKDIR, "data/interim/2DG_CCF_Integration", SOURCE_DATE)
CONTROL_DIR <- file.path(WKDIR, "data/interim/2DG_CCF_Integration", CONTROL_DATE)
OUT_DIR <- file.path(WKDIR, "data/interim/2DG_CCF_Integration", ANALYSIS_DATE)
FIGURE_DIR <- file.path(WKDIR, "figures/2DG_CCF_Integration")
TABLE_DIR <- file.path(WKDIR, "tables/2DG_CCF_Integration")

INTRUDER_OVERLAY <- file.path(
  SOURCE_DIR, paste0("2dg_overlay_points_ccf_scalar-t_", SOURCE_DATE, ".csv.gz")
)
BASELINE_OVERLAY <- file.path(
  CONTROL_DIR, paste0("2dg_overlay_points_ccf_baseline_wildr_vs_spf_", CONTROL_DATE, ".csv.gz")
)
BORDER_POINTS <- file.path(
  CONTROL_DIR, paste0("allen_ccf_border_points_subsampled_", CONTROL_DATE, ".csv.gz")
)
JOINED <- file.path(
  SOURCE_DIR, paste0("seqFISH_neuroestimator_ccf_strict_", SOURCE_DATE, ".rds")
)
CONCORDANCE <- file.path(
  TABLE_DIR, paste0("2DG_SeqFISH_CCF_region_concordance_strict_", SOURCE_DATE, ".csv")
)
BASELINE_REGIONS <- file.path(
  TABLE_DIR, paste0("2DG_CCF_region_summary_baseline_wildr_vs_spf_", CONTROL_DATE, ".csv")
)
INTRUDER_REGIONS <- file.path(
  TABLE_DIR, paste0("2DG_CCF_region_summary_strict_", SOURCE_DATE, ".csv")
)
LIMBIC_SET <- file.path(
  TABLE_DIR, paste0("2DG_prespecified_limbic_hypothalamic_annotations_", CONTROL_DATE, ".csv")
)
SHELL_MESH_DIR <- file.path(WKDIR, "data/interim/2DG_CCF_Integration", "2026-11-14d")
SHELL_VERTS <- file.path(SHELL_MESH_DIR, "allen_ccf_shell_vertices.csv.gz")
SHELL_FACES <- file.path(SHELL_MESH_DIR, "allen_ccf_shell_faces.csv.gz")
SHELL_EDGES <- file.path(SHELL_MESH_DIR, "allen_ccf_shell_edges.csv.gz")
HEMISPHERE_PROVENANCE <- file.path(
  TABLE_DIR, "seqFISH_section_hemisphere_provenance_2026-11-14.csv"
)
ML_MIDLINE_UM <- 5700

WIDGET_3D_HTML <- file.path(FIGURE_DIR, paste0("2DG_CCF_interactive3d_", ANALYSIS_DATE, ".html"))
WIDGET_3D_RDS <- file.path(OUT_DIR, paste0("2DG_CCF_interactive3d_widget_", ANALYSIS_DATE, ".rds"))
TRANSITION_HTML <- file.path(FIGURE_DIR, paste0("2DG_contrast_transition_", ANALYSIS_DATE, ".html"))
TRANSITION_RDS <- file.path(OUT_DIR, paste0("2DG_contrast_transition_widget_", ANALYSIS_DATE, ".rds"))
TRANSITION_PAGE <- file.path(FIGURE_DIR, paste0("2DG_contrast_transition_fullscreen_", ANALYSIS_DATE, ".html"))
COVAR_HTML <- file.path(FIGURE_DIR, paste0("2DG_region_covariation_", ANALYSIS_DATE, ".html"))
COVAR_RDS <- file.path(OUT_DIR, paste0("2DG_region_covariation_widget_", ANALYSIS_DATE, ".rds"))
COVAR_TABLE <- file.path(TABLE_DIR, paste0("2DG_region_covariation_", ANALYSIS_DATE, ".csv"))
COVAR_PNG <- file.path(FIGURE_DIR, paste0("2DG_region_covariation_", ANALYSIS_DATE, ".png"))
COLOR_NOTE <- file.path(OUT_DIR, paste0("2DG_viewer_color_sampling_note_", ANALYSIS_DATE, ".json"))

N_MAP_LOW <- 8000L
N_MAP_MED <- 20000L
N_MAP_HIGH <- 35000L
N_SEQ_LOW <- 5000L
N_SEQ_MED <- 12000L
N_SEQ_HIGH <- 25000L
SET_SEED <- 42L
# Near-fullscreen vertical canvas; scene uses the full plot domain (no bottom pin).
PLOT_HEIGHT_3D <- 1600L
PLOT_HEIGHT_TRANSITION <- 1600L
PLOT_WIDTH <- 1500L

log_message <- function(...) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
}

stop_if_exists <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing)) {
    stop("Refusing to overwrite:\n", paste(existing, collapse = "\n"), call. = FALSE)
  }
}

sample_n <- function(dt, n) {
  if (nrow(dt) <= n) return(dt)
  dt[sample.int(.N, n)]
}

load_map_points <- function(path, n_high) {
  dt <- fread(path, select = c("ccfx", "ccfy", "ccfz", "t_value", "annotation", "annotation_name"))
  dt <- dt[is.finite(t_value) & is.finite(ccfx)]
  dt[, t_clip := pmax(pmin(t_value, 4), -4)]
  sample_n(dt, n_high)
}

t_colorbar <- function(title_prefix, x = 1.02) {
  list(
    title = list(
      text = paste0(
        title_prefix, "<br>",
        "<span style='font-size:11px'>blue = WildR↑&nbsp;&nbsp;|&nbsp;&nbsp;red = SPF↑</span>"
      ),
      side = "right"
    ),
    tickvals = c(-4, -2, 0, 2, 4),
    ticktext = c("SPF↑ −4", "−2", "0", "+2", "WildR↑ +4"),
    len = 0.72,
    y = 0.45,
    x = x,
    thickness = 16
  )
}

ne_colorbar <- function(x = 1.12) {
  list(
    title = list(
      text = "NEUROeSTIMator<br>predicted activity<br><span style='font-size:11px'>Viridis: low → high</span>",
      side = "right"
    ),
    len = 0.55,
    y = 0.35,
    x = x,
    thickness = 14
  )
}

scene_layout <- function() {
  # Full-plot domain: prior y=c(0.02,0.78) pinned the brain to the bottom of tall canvases.
  list(
    domain = list(x = c(0.0, 0.92), y = c(0.0, 1.0)),
    xaxis = list(
      title = "CCF AP (µm)",
      backgroundcolor = "#ffffff",
      gridcolor = "#e6e6e6",
      showbackground = TRUE,
      zeroline = FALSE
    ),
    yaxis = list(
      title = "CCF ML (µm)",
      backgroundcolor = "#ffffff",
      gridcolor = "#e6e6e6",
      showbackground = TRUE,
      zeroline = FALSE
    ),
    zaxis = list(
      title = "CCF DV (µm)",
      autorange = "reversed",
      backgroundcolor = "#ffffff",
      gridcolor = "#e6e6e6",
      showbackground = TRUE,
      zeroline = FALSE
    ),
    aspectmode = "data",
    camera = list(
      eye = list(x = 1.55, y = 1.45, z = 0.85),
      center = list(x = 0, y = 0, z = 0),
      up = list(x = 0, y = 0, z = 1)
    ),
    bgcolor = "#ffffff"
  )
}

apply_anatomical_ml <- function(dt, sample_col = "sample_id") {
  # Un-mirror sections that were reflected into the atlas medial half during
  # seqFISH registration (workbook column `reflect`). Keeps registered ML in
  # ccfz_registered; plotting uses ccfz as anatomical ML.
  if (!file.exists(HEMISPHERE_PROVENANCE)) {
    stop("Missing hemisphere provenance: ", HEMISPHERE_PROVENANCE, call. = FALSE)
  }
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
  out[reflect == TRUE, ccfz := (2 * ML_MIDLINE_UM) - ccfz_registered]
  out[, anatomical_hemisphere := fifelse(
    abs(ccfz - ML_MIDLINE_UM) <= 250,
    "midline_band",
    fifelse(ccfz < ML_MIDLINE_UM, "left", "right")
  )]
  out
}

load_ccf_shell <- function() {
  if (!file.exists(SHELL_VERTS) || !file.exists(SHELL_FACES)) {
    stop("CCF shell mesh missing. Run export_allen_ccf_shell_mesh.py first.", call. = FALSE)
  }
  verts <- fread(SHELL_VERTS)
  faces <- fread(SHELL_FACES)
  edges <- if (file.exists(SHELL_EDGES)) fread(SHELL_EDGES) else NULL
  # Subsample edges for a sketched outline without exploding widget size.
  if (!is.null(edges) && nrow(edges) > 12000L) {
    set.seed(SET_SEED)
    edges <- edges[sample.int(.N, 12000L)]
  }
  list(verts = verts, faces = faces, edges = edges)
}

add_ccf_shell_traces <- function(p, shell, visible = FALSE, opacity = 0.16) {
  verts <- shell$verts
  faces <- shell$faces
  p <- add_trace(
    p,
    type = "mesh3d",
    x = verts$ap,
    y = verts$ml,
    z = verts$dv,
    i = faces$i,
    j = faces$j,
    k = faces$k,
    name = "CCF outline shell",
    opacity = opacity,
    color = "#d9dde3",
    flatshading = FALSE,
    lighting = list(
      ambient = 0.78,
      diffuse = 0.55,
      specular = 0.05,
      roughness = 0.9,
      fresnel = 0.05
    ),
    lightposition = list(x = 1e4, y = 1e4, z = 0),
    hoverinfo = "skip",
    visible = visible,
    showlegend = TRUE
  )
  if (!is.null(shell$edges) && nrow(shell$edges) > 0) {
    e <- shell$edges
    # Line segments with NA breaks
    x <- as.numeric(rbind(verts$ap[e$a + 1L], verts$ap[e$b + 1L], NA_real_))
    y <- as.numeric(rbind(verts$ml[e$a + 1L], verts$ml[e$b + 1L], NA_real_))
    z <- as.numeric(rbind(verts$dv[e$a + 1L], verts$dv[e$b + 1L], NA_real_))
    p <- add_trace(
      p,
      type = "scatter3d",
      mode = "lines",
      x = x,
      y = y,
      z = z,
      name = "CCF outline wires",
      line = list(color = "rgba(40,44,52,0.35)", width = 1.1),
      hoverinfo = "skip",
      visible = visible,
      showlegend = TRUE
    )
  }
  p
}

build_explorer <- function(intruder, baseline, borders, seqfish, shell, ap_min, ap_max) {
  set.seed(SET_SEED)
  intr_low <- sample_n(intruder, N_MAP_LOW)
  intr_med <- sample_n(intruder, N_MAP_MED)
  intr_high <- sample_n(intruder, N_MAP_HIGH)
  base_med <- sample_n(baseline, N_MAP_MED)
  base_high <- sample_n(baseline, N_MAP_HIGH)
  seq_low <- sample_n(seqfish, N_SEQ_LOW)
  seq_med <- sample_n(seqfish, N_SEQ_MED)
  seq_high <- sample_n(seqfish, N_SEQ_HIGH)

  intr_overlap <- intr_med[ccfx >= ap_min & ccfx <= ap_max]
  base_overlap <- base_med[ccfx >= ap_min & ccfx <= ap_max]
  seq_overlap <- seq_med[ccfx >= ap_min & ccfx <= ap_max]

  add_map_trace <- function(p, data, name, visible, showscale = FALSE) {
    add_trace(
      p,
      data = data,
      x = ~ccfx,
      y = ~ccfz,
      z = ~ccfy,
      type = "scatter3d",
      mode = "markers",
      name = name,
      marker = list(
        size = 2.2,
        opacity = 0.55,
        color = ~t_clip,
        colorscale = "RdBu",
        reversescale = TRUE,
        cmin = -4,
        cmax = 4,
        showscale = showscale,
        colorbar = if (showscale) t_colorbar("2DG t (WildR − SPF)") else NULL
      ),
      text = ~sprintf(
        "%s<br>t=%.2f (WildR−SPF)<br>%s<br>%s",
        name,
        t_value,
        ifelse(t_value >= 0, "WildR higher", "SPF higher"),
        annotation_name
      ),
      hoverinfo = "text",
      visible = visible
    )
  }

  add_seq_trace <- function(p, data, name, visible, showscale = FALSE) {
    add_trace(
      p,
      data = data,
      x = ~ccfx,
      y = ~ccfz,
      z = ~ccfy,
      type = "scatter3d",
      mode = "markers",
      name = name,
      marker = list(
        size = 2.0,
        opacity = 0.40,
        color = ~predicted_activity,
        colorscale = "Viridis",
        showscale = showscale,
        colorbar = if (showscale) ne_colorbar() else NULL
      ),
      text = ~sprintf(
        "seqFISH %s<br>NE activity=%.3f<br>gate=%s<br>hemi=%s<br>reflect=%s<br>ML anat=%.0f (reg=%.0f)",
        group,
        predicted_activity,
        final_gate,
        anatomical_hemisphere,
        reflect,
        ccfz,
        ccfz_registered
      ),
      hoverinfo = "text",
      visible = visible
    )
  }

  p <- plot_ly(width = PLOT_WIDTH, height = PLOT_HEIGHT_3D)
  p <- add_map_trace(p, intr_med, "Intruder WildR−SPF (med)", TRUE, TRUE)
  p <- add_map_trace(p, base_med, "Baseline WildR−SPF (med)", FALSE, FALSE)
  p <- add_map_trace(p, intr_high, "Intruder WildR−SPF (high dens)", FALSE, FALSE)
  p <- add_map_trace(p, base_high, "Baseline WildR−SPF (high dens)", FALSE, FALSE)
  p <- add_map_trace(p, intr_low, "Intruder WildR−SPF (low dens)", FALSE, FALSE)
  p <- add_trace(
    p,
    data = borders,
    x = ~ccfx,
    y = ~ccfz,
    z = ~ccfy,
    type = "scatter3d",
    mode = "markers",
    name = "CCF region borders",
    marker = list(size = 1.3, opacity = 0.16, color = "#333333"),
    hoverinfo = "skip",
    visible = "legendonly"
  )
  p <- add_seq_trace(p, seq_med, "seqFISH NE (med)", FALSE, TRUE)
  p <- add_seq_trace(p, seq_high, "seqFISH NE (high dens)", FALSE, FALSE)
  p <- add_seq_trace(p, seq_low, "seqFISH NE (low dens)", FALSE, FALSE)
  p <- add_map_trace(p, intr_overlap, "Intruder ∩ seqFISH AP", FALSE, FALSE)
  p <- add_map_trace(p, base_overlap, "Baseline ∩ seqFISH AP", FALSE, FALSE)
  p <- add_seq_trace(p, seq_overlap, "seqFISH in AP window", FALSE, FALSE)
  # Traces 13–14: translucent Allen CCF shell + sketched wires
  p <- add_ccf_shell_traces(p, shell, visible = FALSE, opacity = 0.16)

  # 14 traces total
  n_tr <- 14L
  vis <- function(..., shell_on = FALSE, points_on = FALSE) {
    vals <- rep(FALSE, n_tr)
    idx <- c(...)
    if (length(idx)) vals[idx] <- TRUE
    if (points_on) vals[6] <- TRUE
    if (shell_on) {
      vals[13] <- TRUE
      vals[14] <- TRUE
    }
    vals
  }
  # 0-based indices for restyle on mesh/wires
  mesh_idx <- list(12, 13)

  p <- layout(
    p,
    title = list(
      text = paste0(
        "CCF explorer — 2DG heatmaps + seqFISH<br>",
        "<sup>Polarity: blue = WildR higher, red = SPF higher. ",
        "seqFISH ML is anatomical (reflect=TRUE sections un-mirrored about 5700 µm). ",
        "Anatomy: CCF points or translucent Allen outline shell with α controls. ",
        "AP window ", ap_min, "–", ap_max, " µm</sup>"
      ),
      x = 0.02,
      xanchor = "left",
      y = 0.98
    ),
    scene = scene_layout(),
    margin = list(l = 0, r = 80, b = 30, t = 110),
    legend = list(orientation = "h", y = -0.04, x = 0, font = list(size = 11)),
    paper_bgcolor = "#ffffff",
    plot_bgcolor = "#ffffff",
    updatemenus = list(
      list(
        type = "buttons",
        direction = "right",
        x = 0.0,
        xanchor = "left",
        y = 1.12,
        yanchor = "top",
        pad = list(r = 4, t = 2),
        bgcolor = "#f3f3f3",
        buttons = list(
          list(
            label = "Intruder",
            method = "update",
            args = list(list(visible = vis(1)))
          ),
          list(
            label = "Baseline",
            method = "update",
            args = list(list(visible = vis(2)))
          ),
          list(
            label = "Baseline + Intruder",
            method = "update",
            args = list(list(visible = vis(1, 2)))
          ),
          list(
            label = "Intruder + CCF points",
            method = "update",
            args = list(list(visible = vis(1, points_on = TRUE)))
          ),
          list(
            label = "Intruder + CCF outline",
            method = "update",
            args = list(list(visible = vis(1, shell_on = TRUE)))
          ),
          list(
            label = "seqFISH AP overlap",
            method = "update",
            args = list(list(visible = vis(10, 11, 12, points_on = TRUE)))
          ),
          list(
            label = "All maps + outline",
            method = "update",
            args = list(list(visible = vis(1, 2, 7, shell_on = TRUE)))
          )
        )
      ),
      list(
        type = "buttons",
        direction = "right",
        x = 0.0,
        xanchor = "left",
        y = 1.07,
        yanchor = "top",
        pad = list(r = 4, t = 2),
        bgcolor = "#f3f3f3",
        buttons = list(
          list(label = "Density low", method = "update", args = list(list(visible = vis(5, points_on = TRUE)))),
          list(label = "Density med", method = "update", args = list(list(visible = vis(1, points_on = TRUE)))),
          list(label = "Density high", method = "update", args = list(list(visible = vis(3, points_on = TRUE)))),
          list(label = "Larger points", method = "restyle", args = list(list(marker.size = 3.6))),
          list(label = "Default point size", method = "restyle", args = list(list(marker.size = 2.2)))
        )
      ),
      list(
        type = "buttons",
        direction = "right",
        x = 0.0,
        xanchor = "left",
        y = 1.02,
        yanchor = "top",
        pad = list(r = 4, t = 2),
        bgcolor = "#f3f3f3",
        buttons = list(
          list(
            label = "Shell α 0.08",
            method = "restyle",
            args = list(list(opacity = 0.08), list(12))
          ),
          list(
            label = "Shell α 0.16",
            method = "restyle",
            args = list(list(opacity = 0.16), list(12))
          ),
          list(
            label = "Shell α 0.28",
            method = "restyle",
            args = list(list(opacity = 0.28), list(12))
          ),
          list(
            label = "Shell α 0.45",
            method = "restyle",
            args = list(list(opacity = 0.45), list(12))
          ),
          list(
            label = "Wires on",
            method = "restyle",
            args = list(list(visible = TRUE), list(13))
          ),
          list(
            label = "Wires off",
            method = "restyle",
            args = list(list(visible = FALSE), list(13))
          )
        )
      )
    )
  )
  config(p, displaylogo = FALSE, responsive = TRUE)
}

build_transition <- function(intruder, baseline, region_centroids, shell) {
  set.seed(SET_SEED)
  base <- sample_n(baseline, N_MAP_MED)
  intr <- sample_n(intruder, N_MAP_MED)
  ne <- copy(region_centroids)
  ne[, t_clip := pmax(pmin(seqfish_wildr_minus_spf_descriptive * 8, 4), -4)]
  ne[, t_value := seqfish_wildr_minus_spf_descriptive]

  add_frame <- function(p, dt, label, size = 2.2, opacity = 0.55) {
    add_trace(
      p,
      data = dt,
      x = ~ccfx,
      y = ~ccfz,
      z = ~ccfy,
      type = "scatter3d",
      mode = "markers",
      frame = label,
      marker = list(
        size = size,
        opacity = opacity,
        color = ~t_clip,
        colorscale = "RdBu",
        reversescale = TRUE,
        cmin = -4,
        cmax = 4,
        colorbar = t_colorbar("Signed effect (WildR − SPF)")
      ),
      text = ~sprintf(
        "%s<br>%.3f<br>%s<br>%s",
        label,
        t_value,
        ifelse(t_value >= 0, "WildR higher", "SPF higher"),
        annotation_name
      ),
      hoverinfo = "text",
      showlegend = FALSE
    )
  }
  p <- plot_ly(width = PLOT_WIDTH, height = PLOT_HEIGHT_TRANSITION)
  p <- add_frame(p, base, "Baseline WildR−SPF")
  p <- add_frame(p, intr, "Intruder WildR−SPF")
  p <- add_frame(p, ne, "NEUROeSTIMator WildR−SPF", size = 5.5, opacity = 0.9)
  # Static translucent shell behind all frames
  p <- add_ccf_shell_traces(p, shell, visible = TRUE, opacity = 0.12)
  # 0-based indices of mesh / wire traces (last two added)
  n_data <- length(p$x$data)
  mesh_i <- n_data - 2L
  wire_i <- n_data - 1L
  p <- layout(
    p,
    title = list(
      text = paste0(
        "Contrast transition — Baseline → Intruder → NEUROeSTIMator<br>",
        "<sup>All frames: blue = WildR higher, red = SPF higher (NE Δ ×8). ",
        "Allen CCF outline shell shown with α controls below.</sup>"
      ),
      x = 0.02,
      xanchor = "left"
    ),
    scene = scene_layout(),
    margin = list(l = 0, r = 80, b = 60, t = 90),
    paper_bgcolor = "#ffffff",
    plot_bgcolor = "#ffffff",
    updatemenus = list(
      list(
        type = "buttons",
        direction = "right",
        x = 0.0,
        xanchor = "left",
        y = 1.06,
        bgcolor = "#f3f3f3",
        font = list(color = "#222222"),
        buttons = list(
          list(label = "Shell α 0.08", method = "restyle", args = list(list(opacity = 0.08), list(mesh_i))),
          list(label = "Shell α 0.12", method = "restyle", args = list(list(opacity = 0.12), list(mesh_i))),
          list(label = "Shell α 0.25", method = "restyle", args = list(list(opacity = 0.25), list(mesh_i))),
          list(label = "Shell α 0.40", method = "restyle", args = list(list(opacity = 0.40), list(mesh_i))),
          list(
            label = "Hide shell",
            method = "restyle",
            args = list(list(visible = FALSE), list(mesh_i, wire_i))
          ),
          list(
            label = "Show shell",
            method = "restyle",
            args = list(list(visible = TRUE), list(mesh_i, wire_i))
          )
        )
      )
    )
  )
  p <- animation_opts(p, frame = 1200, transition = 800, redraw = TRUE)
  p <- animation_slider(p, currentvalue = list(prefix = "", font = list(size = 13)))
  p <- animation_button(p, label = "Play / hold cycles", x = 0.02, y = 0)
  config(p, displaylogo = FALSE, responsive = TRUE)
}

spearman_label <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5L) return("ρ = NA")
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  sprintf("Spearman ρ = %.3f (p = %.3g, n = %d)", unname(ct$estimate), ct$p.value, sum(ok))
}

build_covariation <- function(covar) {
  covar <- copy(covar)
  covar[, set_label := ifelse(in_prespecified_set %in% TRUE, "Limbic/hypothalamic set", "Other regions")]
  covar[, size_cells := pmax(7, 3 + 4 * log10(pmax(n_cells_total, 1)))]

  mk <- function(xcol, xlab, name, color) {
    plot_ly(
      covar,
      x = as.formula(paste0("~", xcol)),
      y = ~seqfish_wildr_minus_spf_descriptive,
      type = "scatter",
      mode = "markers",
      name = name,
      color = ~set_label,
      colors = c("Limbic/hypothalamic set" = "#FF6C0C", "Other regions" = color),
      marker = list(size = ~size_cells, opacity = 0.75, line = list(width = 0.4, color = "#333")),
      text = ~sprintf(
        "%s<br>%s<br>2DG mean t=%.2f<br>NE WildR−SPF=%.4f<br>cells=%d; overlay px=%d",
        annotation_name,
        set_label,
        get(xcol),
        seqfish_wildr_minus_spf_descriptive,
        n_cells_total,
        n_overlay_pixels
      ),
      hoverinfo = "text",
      showlegend = FALSE
    ) %>%
      layout(
        xaxis = list(
          title = paste0(xlab, "<br><sup>→ WildR higher  |  SPF higher ←</sup>"),
          zeroline = TRUE,
          zerolinecolor = "#888"
        ),
        yaxis = list(
          title = "NEUROeSTIMator WildR − SPF<br><sup>→ WildR higher activity</sup>",
          zeroline = TRUE,
          zerolinecolor = "#888"
        ),
        annotations = list(
          list(
            text = spearman_label(
              covar[[xcol]],
              covar$seqfish_wildr_minus_spf_descriptive
            ),
            xref = "paper",
            yref = "paper",
            x = 0.02,
            y = 1.08,
            showarrow = FALSE,
            align = "left",
            font = list(size = 12)
          )
        )
      )
  }

  p1 <- mk("intruder_mean_t", "Intruder 2DG mean t (WildR − SPF)", "vs Intruder", "#9ecae1")
  p2 <- mk("baseline_mean_t", "Baseline 2DG mean t (WildR − SPF)", "vs Baseline", "#a1d99b")
  p3 <- plot_ly(
    covar,
    x = ~baseline_mean_t,
    y = ~intruder_mean_t,
    type = "scatter",
    mode = "markers",
    name = "Baseline vs Intruder",
    marker = list(
      size = ~size_cells,
      opacity = 0.8,
      color = ~seqfish_wildr_minus_spf_descriptive,
      colorscale = "RdBu",
      reversescale = TRUE,
      colorbar = list(
        title = list(
          text = "NE Δ<br><span style='font-size:11px'>blue=WildR↑</span>"
        ),
        tickvals = NULL
      ),
      cmin = -max(abs(covar$seqfish_wildr_minus_spf_descriptive), na.rm = TRUE),
      cmax = max(abs(covar$seqfish_wildr_minus_spf_descriptive), na.rm = TRUE),
      line = list(width = 0.4, color = "#333")
    ),
    text = ~sprintf(
      "%s<br>Baseline t=%.2f<br>Intruder t=%.2f<br>NE Δ=%.4f",
      annotation_name,
      baseline_mean_t,
      intruder_mean_t,
      seqfish_wildr_minus_spf_descriptive
    ),
    hoverinfo = "text",
    showlegend = FALSE
  ) %>%
    layout(
      xaxis = list(title = "Baseline 2DG mean t (WildR − SPF)", zeroline = TRUE),
      yaxis = list(title = "Intruder 2DG mean t (WildR − SPF)", zeroline = TRUE),
      annotations = list(list(
        text = spearman_label(covar$baseline_mean_t, covar$intruder_mean_t),
        xref = "paper", yref = "paper", x = 0.02, y = 1.08,
        showarrow = FALSE, font = list(size = 12)
      ))
    )

  subplot(p1, p2, p3, nrows = 1, shareY = FALSE, titleX = TRUE, titleY = TRUE, margin = 0.06) %>%
    layout(
      title = list(
        text = paste0(
          "Per-region co-variation — NEUROeSTIMator vs 2DG effects<br>",
          "<sup>Filtered Allen parcels with ≥100 seqFISH cells and ≥10 Intruder overlay pixels. ",
          "Orange = pre-specified limbic/hypothalamic set. Marker size ∝ log10(cells).</sup>"
        ),
        x = 0.02,
        xanchor = "left",
        y = 0.98
      ),
      height = 720,
      width = PLOT_WIDTH,
      showlegend = TRUE,
      legend = list(orientation = "h", y = -0.12),
      margin = list(t = 90, b = 80, l = 70, r = 40)
    ) %>%
    config(displaylogo = FALSE, responsive = TRUE)
}

build_covariation_static <- function(covar, path) {
  covar <- copy(covar)
  covar[, set_label := ifelse(in_prespecified_set %in% TRUE, "Limbic/hypothalamic set", "Other regions")]
  long <- rbindlist(list(
    covar[, .(
      annotation_name, set_label, n_cells_total,
      x = intruder_mean_t,
      y = seqfish_wildr_minus_spf_descriptive,
      panel = "NE Δ vs Intruder 2DG t"
    )],
    covar[, .(
      annotation_name, set_label, n_cells_total,
      x = baseline_mean_t,
      y = seqfish_wildr_minus_spf_descriptive,
      panel = "NE Δ vs Baseline 2DG t"
    )],
    covar[, .(
      annotation_name, set_label, n_cells_total,
      x = baseline_mean_t,
      y = intruder_mean_t,
      panel = "Intruder t vs Baseline t"
    )]
  ))
  long[, panel := factor(
    panel,
    levels = c(
      "NE Δ vs Intruder 2DG t",
      "NE Δ vs Baseline 2DG t",
      "Intruder t vs Baseline t"
    )
  )]
  stats <- long[
    ,
    {
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 5L) {
        .(label = "ρ = NA")
      } else {
        ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
        .(label = sprintf("ρ = %.3f\np = %.2g\nn = %d", unname(ct$estimate), ct$p.value, sum(ok)))
      }
    },
    by = panel
  ]
  p <- ggplot(long, aes(x = x, y = y)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_point(
      aes(colour = set_label, size = n_cells_total),
      alpha = 0.75
    ) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.6, colour = "#222222", fill = "#cccccc") +
    geom_text(
      data = stats,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = 1.2,
      size = 3.2
    ) +
    facet_wrap(~panel, scales = "free", nrow = 1) +
    scale_colour_manual(
      values = c(
        "Limbic/hypothalamic set" = "#FF6C0C",
        "Other regions" = "#4C78A8"
      )
    ) +
    scale_size_area(max_size = 8, breaks = c(100, 500, 2000, 8000)) +
    labs(
      title = "Per-region co-variation: NEUROeSTIMator WildR−SPF vs 2DG effects",
      subtitle = paste0(
        "Positive on either axis = WildR higher. Filtered: ≥100 cells and ≥10 Intruder overlay pixels. ",
        "Orange = pre-specified limbic/hypothalamic set."
      ),
      x = "2DG regional mean t (WildR − SPF)  |  right panel x = Baseline",
      y = "NEUROeSTIMator WildR − SPF  |  right panel y = Intruder t",
      colour = NULL,
      size = "seqFISH cells"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "#f2f2f2", colour = NA),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9)
    )
  ggsave(path, p, width = 14.5, height = 5.2, dpi = 200, bg = "white")
}

write_fullscreen_page <- function(widget, path) {
  saveWidget(widget, path, selfcontained = TRUE, title = "2DG contrast transition")
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  inject <- paste0(
    "<style>",
    "html,body{margin:0;padding:0;background:#ffffff;color:#222;}",
    ".twodg-bar{position:sticky;top:0;z-index:50;padding:10px 16px;",
    "background:#ffffff;border-bottom:1px solid #ddd;display:flex;",
    "gap:12px;align-items:center;font-family:Helvetica,Arial,sans-serif;}",
    ".twodg-bar button{background:#FF6C0C;color:#111;border:0;padding:10px 16px;",
    "border-radius:8px;font-weight:700;cursor:pointer;}",
    ".js-plotly-plot,.plotly{min-height:calc(100vh - 60px) !important;}",
    "</style>",
    "<div class='twodg-bar'><strong>Contrast transition</strong>",
    "<span>Baseline → Intruder → NEUROeSTIMator. Blue = WildR higher; red = SPF higher.</span>",
    "<button type='button' onclick=\"(function(){",
    "var el=document.querySelector('.js-plotly-plot');",
    "if(!el||!window.Plotly)return;",
    "Plotly.animate(el);",
    "})()\">Play transition</button></div>"
  )
  if (grepl("<body[^>]*>", html, ignore.case = TRUE)) {
    html <- sub("(<body[^>]*>)", paste0("\\1\n", inject), html, ignore.case = TRUE)
  } else {
    html <- paste0(inject, "\n", html)
  }
  writeLines(html, path)
}

region_centroids_from_joined <- function(joined, concordance) {
  cells <- joined[
    registered_ccf == TRUE &
      is.finite(ccfx) &
      !is.na(annotation) &
      annotation != 0
  ][
    ,
    .(
      ccfx,
      ccfy,
      ccfz,
      annotation = as.integer(annotation),
      sample_id = as.character(sample_id)
    )
  ]
  cells <- apply_anatomical_ml(cells)
  cents <- cells[
    ,
    .(
      ccfx = mean(ccfx),
      ccfy = mean(ccfy),
      ccfz = mean(ccfz),
      n_cells = .N,
      n_left = sum(anatomical_hemisphere == "left"),
      n_right = sum(anatomical_hemisphere == "right")
    ),
    by = .(annotation)
  ]
  merge(
    cents,
    concordance[
      ,
      .(
        annotation = as.integer(annotation),
        annotation_name,
        seqfish_wildr_minus_spf_descriptive
      )
    ],
    by = "annotation",
    all = FALSE
  )[is.finite(seqfish_wildr_minus_spf_descriptive)]
}

build_covariation_table <- function() {
  conc <- fread(CONCORDANCE)
  base <- fread(BASELINE_REGIONS)
  intr <- fread(INTRUDER_REGIONS)
  limbic <- if (file.exists(LIMBIC_SET)) {
    fread(LIMBIC_SET)[, .(annotation = as.integer(annotation), in_prespecified_set)]
  } else {
    data.table(annotation = integer(), in_prespecified_set = logical())
  }
  conc[, annotation := as.integer(annotation)]
  base[, annotation := as.integer(annotation)]
  intr[, annotation := as.integer(annotation)]
  out <- merge(
    conc[
      ,
      .(
        annotation,
        annotation_name,
        seqfish_wildr_minus_spf_descriptive,
        n_cells_spf,
        n_cells_wildr,
        n_overlay_pixels,
        intruder_mean_t = mean_t_value
      )
    ],
    base[, .(annotation, baseline_mean_t = mean_t_value)],
    by = "annotation",
    all = FALSE
  )
  if (any(!is.finite(out$intruder_mean_t))) {
    out <- merge(
      out,
      intr[, .(annotation, intruder_mean_t_fill = mean_t_value)],
      by = "annotation",
      all.x = TRUE
    )
    out[
      !is.finite(intruder_mean_t) & is.finite(intruder_mean_t_fill),
      intruder_mean_t := intruder_mean_t_fill
    ]
    out[, intruder_mean_t_fill := NULL]
  }
  out <- merge(out, limbic, by = "annotation", all.x = TRUE)
  out[is.na(in_prespecified_set), in_prespecified_set := FALSE]
  out[, n_cells_total := n_cells_spf + n_cells_wildr]
  out <- out[
    is.finite(intruder_mean_t) &
      is.finite(baseline_mean_t) &
      is.finite(seqfish_wildr_minus_spf_descriptive)
  ]
  # Analysis-quality filter for displayed co-variation
  out[n_cells_total >= 100L & n_overlay_pixels >= 10L]
}

main <- function() {
  outputs <- c(
    WIDGET_3D_HTML, WIDGET_3D_RDS,
    TRANSITION_HTML, TRANSITION_RDS, TRANSITION_PAGE,
    COVAR_HTML, COVAR_RDS, COVAR_TABLE, COVAR_PNG, COLOR_NOTE
  )
  stop_if_exists(outputs)
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
  set.seed(SET_SEED)

  log_message("Loading map point clouds")
  intruder <- load_map_points(INTRUDER_OVERLAY, N_MAP_HIGH)
  baseline <- load_map_points(BASELINE_OVERLAY, N_MAP_HIGH)
  borders <- fread(BORDER_POINTS)
  log_message("Loading Allen CCF outline shell mesh")
  shell <- load_ccf_shell()

  log_message("Loading seqFISH (anatomical ML: un-mirror reflect=TRUE sections)")
  joined <- as.data.table(readRDS(JOINED))
  seqfish <- joined[
    registered_ccf == TRUE &
      is.finite(ccfx) &
      is.finite(predicted_activity)
  ][
    ,
    .(
      ccfx,
      ccfy,
      ccfz,
      predicted_activity,
      group = as.character(group),
      final_gate = as.character(final_gate),
      annotation = as.integer(annotation),
      sample_id = as.character(sample_id)
    )
  ]
  seqfish <- apply_anatomical_ml(seqfish)
  log_message(
    "seqFISH anatomical hemispheres: ",
    paste(names(table(seqfish$anatomical_hemisphere)), table(seqfish$anatomical_hemisphere), sep = "=", collapse = ", ")
  )
  ap_min <- as.numeric(min(seqfish$ccfx, na.rm = TRUE))
  ap_max <- as.numeric(max(seqfish$ccfx, na.rm = TRUE))

  concordance <- fread(CONCORDANCE)
  centroids <- region_centroids_from_joined(joined, concordance)
  covar <- build_covariation_table()
  fwrite(covar, COVAR_TABLE)

  note <- sprintf(
    paste0(
      "{\n",
      "  \"analysis_date\": \"%s\",\n",
      "  \"polarity\": {\n",
      "    \"signed_maps\": \"WildR minus SPF\",\n",
      "    \"positive_or_blue\": \"WildR higher\",\n",
      "    \"negative_or_red\": \"SPF higher\"\n",
      "  },\n",
      "  \"color_model\": {\n",
      "    \"twodg\": \"signed t on randomly subsampled significant overlay voxels; RdBu diverging\",\n",
      "    \"seqfish\": \"NEUROeSTIMator predicted_activity on randomly subsampled registered cells; Viridis continuous; NOT stratified by cell type\",\n",
      "    \"ne_difference_frame\": \"Allen-region centroids colored by descriptive WildR-SPF mean predicted_activity\"\n",
      "  },\n",
      "  \"sampling\": {\"map_low\": %d, \"map_med\": %d, \"map_high\": %d, \"seq_low\": %d, \"seq_med\": %d, \"seq_high\": %d, \"seed\": %d},\n",
      "  \"seqfish_ap_window_um\": [%.1f, %.1f],\n",
      "  \"covariation_filter\": {\"min_cells\": 100, \"min_intruder_overlay_pixels\": 10, \"n_regions\": %d}\n",
      "}\n"
    ),
    ANALYSIS_DATE,
    N_MAP_LOW, N_MAP_MED, N_MAP_HIGH, N_SEQ_LOW, N_SEQ_MED, N_SEQ_HIGH, SET_SEED,
    ap_min, ap_max, nrow(covar)
  )
  writeLines(note, COLOR_NOTE)

  log_message("Building 3D explorer")
  explorer <- build_explorer(intruder, baseline, borders, seqfish, shell, ap_min, ap_max)
  saveWidget(explorer, WIDGET_3D_HTML, selfcontained = TRUE, title = "2DG CCF explorer")
  saveRDS(explorer, WIDGET_3D_RDS, compress = "xz")

  log_message("Building contrast transition")
  transition <- build_transition(intruder, baseline, centroids, shell)
  saveWidget(transition, TRANSITION_HTML, selfcontained = TRUE, title = "2DG contrast transition")
  saveRDS(transition, TRANSITION_RDS, compress = "xz")
  write_fullscreen_page(transition, TRANSITION_PAGE)

  log_message("Building region co-variation")
  covar_plot <- build_covariation(covar)
  saveWidget(covar_plot, COVAR_HTML, selfcontained = TRUE, title = "2DG region covariation")
  saveRDS(covar_plot, COVAR_RDS, compress = "xz")
  build_covariation_static(covar, COVAR_PNG)

  log_message("Done: ", ANALYSIS_DATE, " regions_in_covariation=", nrow(covar))
}

main()
