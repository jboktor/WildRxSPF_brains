#!/usr/bin/env Rscript

# Build a fixed-viewport, modular 2DG × seqFISH Allen CCF explorer.
#
# Every control is orthogonal:
#   - Intruder/Baseline share one unambiguous contrast selector;
#   - seqFISH remains an independent layer;
#   - density and AP scope apply symmetrically to both 2DG maps;
#   - whole-brain context, domain, subdomain, surface, and wireframe are
#     independent of the selected data layer;
#   - camera presets do not alter layer visibility.
#
# Run with:
# /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript \
#   notebooks/R_scripts/plot_2dg_ccf_modular_viewer.R

suppressPackageStartupMessages({
  library(data.table)
  library(plotly)
  library(htmlwidgets)
  library(jsonlite)
})

WKDIR <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
SOURCE_DATE <- "2026-11-23c"
CONTROL_DATE <- "2026-11-23c"
ANALYSIS_DATE <- "2026-11-23c"
MESH_DATE <- "2026-11-21e"

SOURCE_DIR <- file.path(WKDIR, "data/interim/2DG_CCF_Integration", SOURCE_DATE)
CONTROL_DIR <- file.path(WKDIR, "data/interim/2DG_CCF_Integration", CONTROL_DATE)
OUT_DIR <- file.path(WKDIR, "data/interim/2DG_CCF_Integration", ANALYSIS_DATE)
MESH_DIR <- file.path(WKDIR, "data/interim/2DG_CCF_Integration", MESH_DATE)
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
HEMISPHERE_PROVENANCE <- file.path(
  TABLE_DIR, "seqFISH_section_hemisphere_provenance_2026-11-14.csv"
)
MESH_MANIFEST <- file.path(
  MESH_DIR, paste0("allen_ccf_region_mesh_manifest_", MESH_DATE, ".json")
)
MESH_VERTICES <- file.path(
  MESH_DIR, paste0("allen_ccf_region_mesh_vertices_", MESH_DATE, ".csv.gz")
)
MESH_FACES <- file.path(
  MESH_DIR, paste0("allen_ccf_region_mesh_faces_", MESH_DATE, ".csv.gz")
)
MESH_EDGES <- file.path(
  MESH_DIR, paste0("allen_ccf_region_mesh_edges_", MESH_DATE, ".csv.gz")
)

VIEWER_HTML <- file.path(
  FIGURE_DIR, paste0("2DG_CCF_interactive3d_", ANALYSIS_DATE, ".html")
)
VIEWER_NOTE <- file.path(
  OUT_DIR, paste0("2DG_viewer_color_sampling_note_", ANALYSIS_DATE, ".json")
)

N_MAP <- c(low = 8000L, medium = 20000L, high = 35000L)
N_SEQ <- c(low = 5000L, medium = 12000L, high = 25000L)
SET_SEED <- 42L
ML_MIDLINE_UM <- 5700

log_message <- function(...) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
}

stop_if_exists <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing)) {
    stop("Refusing to overwrite:\n", paste(existing, collapse = "\n"), call. = FALSE)
  }
}

nested_samples <- function(dt, counts) {
  if (!nrow(dt)) stop("Cannot sample an empty table.", call. = FALSE)
  order <- sample.int(nrow(dt), min(nrow(dt), max(counts)))
  sampled <- dt[order]
  lapply(counts, function(n) sampled[seq_len(min(n, nrow(sampled)))])
}

load_map_points <- function(path) {
  dt <- fread(
    path,
    select = c("ccfx", "ccfy", "ccfz", "t_value", "annotation", "annotation_name")
  )
  dt <- dt[
    is.finite(t_value) &
      is.finite(ccfx) &
      is.finite(ccfy) &
      is.finite(ccfz)
  ]
  if (!nrow(dt)) stop("No finite 2DG points in ", path, call. = FALSE)
  dt[, t_clip := pmax(pmin(t_value, 4), -4)]
  dt
}

apply_anatomical_ml <- function(dt, sample_col = "sample_id") {
  if (!file.exists(HEMISPHERE_PROVENANCE)) {
    stop("Missing hemisphere provenance: ", HEMISPHERE_PROVENANCE, call. = FALSE)
  }
  prov <- fread(HEMISPHERE_PROVENANCE)
  prov[, sample_id := as.character(sample_id)]
  prov[, reflect := as.logical(reflect)]
  if (anyDuplicated(prov$sample_id)) {
    stop("Duplicate sample_id values in hemisphere provenance.", call. = FALSE)
  }
  if (anyNA(prov$sample_id) || anyNA(prov$reflect)) {
    stop("Missing sample_id/reflect values in hemisphere provenance.", call. = FALSE)
  }
  out <- copy(dt)
  out[, (sample_col) := as.character(get(sample_col))]
  missing_samples <- setdiff(unique(out[[sample_col]]), prov$sample_id)
  if (length(missing_samples)) {
    stop(
      "Hemisphere provenance is incomplete for: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }
  out <- merge(
    out,
    prov[, .(sample_id, reflect)],
    by.x = sample_col,
    by.y = "sample_id",
    all.x = TRUE
  )
  out[, ccfz_registered := ccfz]
  out[reflect == TRUE, ccfz := (2 * ML_MIDLINE_UM) - ccfz_registered]
  out
}

load_mesh_catalog <- function() {
  required <- c(MESH_MANIFEST, MESH_VERTICES, MESH_FACES, MESH_EDGES)
  missing <- required[!file.exists(required)]
  if (length(missing)) {
    stop(
      "Missing modular mesh assets. Run export_2dg_ccf_region_meshes.py:\n",
      paste(missing, collapse = "\n"),
      call. = FALSE
    )
  }
  manifest <- fromJSON(MESH_MANIFEST)
  if (!identical(as.character(manifest$analysis_date), MESH_DATE)) {
    stop("Mesh manifest date does not match MESH_DATE.", call. = FALSE)
  }
  mesh <- list(
    regions = as.data.table(manifest$regions),
    vertices = fread(MESH_VERTICES),
    faces = fread(MESH_FACES),
    edges = fread(MESH_EDGES)
  )
  required_region_columns <- c("key", "label", "group", "level", "color", "n_vertices", "n_faces")
  if (!all(required_region_columns %in% names(mesh$regions))) {
    stop("Mesh manifest region schema is incomplete.", call. = FALSE)
  }
  if (anyDuplicated(mesh$regions$key)) {
    stop("Mesh manifest has duplicate region keys.", call. = FALSE)
  }
  region_sets <- lapply(mesh[c("vertices", "faces", "edges")], function(x) unique(x$region))
  expected <- sort(mesh$regions$key)
  if (!all(vapply(region_sets, function(x) identical(sort(x), expected), logical(1)))) {
    stop("Mesh tables and manifest have different region sets.", call. = FALSE)
  }
  for (key in expected) {
    region_key <- key
    n_vertices <- mesh$regions[key == region_key, n_vertices][[1]]
    face_index <- unlist(mesh$faces[region == region_key, .(i, j, k)], use.names = FALSE)
    if (!length(face_index) || min(face_index) < 0 || max(face_index) >= n_vertices) {
      stop("Out-of-range mesh face index for ", key, call. = FALSE)
    }
  }
  mesh
}

scene_layout <- function() {
  list(
    domain = list(x = c(0, 1), y = c(0, 1)),
    xaxis = list(
      title = "CCF AP (µm)",
      backgroundcolor = "#fbfbfc",
      gridcolor = "#eceef1",
      showbackground = TRUE,
      zeroline = FALSE
    ),
    yaxis = list(
      title = "CCF ML (µm)",
      backgroundcolor = "#fbfbfc",
      gridcolor = "#eceef1",
      showbackground = TRUE,
      zeroline = FALSE
    ),
    zaxis = list(
      title = "CCF DV (µm)",
      autorange = "reversed",
      backgroundcolor = "#fbfbfc",
      gridcolor = "#eceef1",
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

add_2dg_trace <- function(p, data, dataset, density, scope, initially_visible = FALSE) {
  label <- if (dataset == "intruder") "Intruder WildR−SPF" else "Baseline WildR−SPF"
  add_trace(
    p,
    data = data,
    x = ~ccfx,
    y = ~ccfz,
    z = ~ccfy,
    type = "scatter3d",
    mode = "markers",
    name = paste(label, density, scope),
    marker = list(
      size = 2.2,
      opacity = 0.58,
      color = ~t_clip,
      colorscale = "RdBu",
      reversescale = FALSE,
      cmin = -4,
      cmax = 4,
      showscale = FALSE
    ),
    text = ~sprintf(
      "%s<br>t = %.2f (WildR−SPF)<br>%s<br>%s",
      label,
      t_value,
      ifelse(t_value >= 0, "WildR higher", "SPF higher"),
      annotation_name
    ),
    hoverinfo = "text",
    visible = initially_visible,
    showlegend = FALSE,
    meta = list(
      role = "data",
      dataset = dataset,
      density = density,
      scope = scope,
      base_size = 2.2
    )
  )
}

add_seqfish_trace <- function(
  p,
  data,
  density,
  ne_min,
  ne_max,
  initially_visible = FALSE
) {
  add_trace(
    p,
    data = data,
    x = ~ccfx,
    y = ~ccfz,
    z = ~ccfy,
    type = "scatter3d",
    mode = "markers",
    name = paste("seqFISH NE", density),
    marker = list(
      size = 2,
      opacity = 0.42,
      color = ~predicted_activity,
      colorscale = "Viridis",
      cmin = ne_min,
      cmax = ne_max,
      showscale = FALSE
    ),
    text = ~sprintf(
      "seqFISH %s<br>NE activity = %.3f<br>gate = %s<br>anatomical ML = %.0f µm",
      group,
      predicted_activity,
      final_gate,
      ccfz
    ),
    hoverinfo = "text",
    visible = initially_visible,
    showlegend = FALSE,
    meta = list(
      role = "data",
      dataset = "seqfish",
      density = density,
      scope = "both",
      base_size = 2
    )
  )
}

add_border_trace <- function(p, borders) {
  add_trace(
    p,
    data = borders,
    x = ~ccfx,
    y = ~ccfz,
    z = ~ccfy,
    type = "scatter3d",
    mode = "markers",
    name = "Allen parcel borders",
    marker = list(size = 1.1, opacity = 0.14, color = "#3e434a"),
    hoverinfo = "skip",
    visible = FALSE,
    showlegend = FALSE,
    meta = list(role = "border")
  )
}

add_region_traces <- function(p, mesh) {
  for (row_index in seq_len(nrow(mesh$regions))) {
    info <- mesh$regions[row_index]
    key <- info$key[[1]]
    verts <- mesh$vertices[region == key]
    faces <- mesh$faces[region == key]
    edges <- mesh$edges[region == key]
    if (!nrow(verts) || !nrow(faces)) next
    edge_cap <- switch(
      info$level[[1]],
      root = 3000L,
      domain = 1500L,
      subdomain = 500L,
      500L
    )
    if (nrow(edges) > edge_cap) edges <- edges[seq_len(edge_cap)]
    initially_root <- key == "root"
    p <- add_trace(
      p,
      type = "mesh3d",
      x = verts$ap,
      y = verts$ml,
      z = verts$dv,
      i = faces$i,
      j = faces$j,
      k = faces$k,
      name = paste(info$label[[1]], "surface"),
      vertexcolor = I(rep(info$color[[1]], nrow(verts))),
      opacity = if (initially_root) 0.09 else 0.34,
      flatshading = FALSE,
      lighting = list(
        ambient = 0.72,
        diffuse = 0.72,
        specular = 0.12,
        roughness = 0.76,
        fresnel = 0.08
      ),
      lightposition = list(x = 10000, y = 14000, z = -3000),
      hovertemplate = paste0(info$label[[1]], "<extra></extra>"),
      visible = initially_root,
      showlegend = FALSE,
      meta = list(
        role = "anatomy",
        kind = "surface",
        key = key,
        level = info$level[[1]],
        group = info$group[[1]],
        base_opacity = if (initially_root) 0.09 else 0.34
      )
    )
    if (nrow(edges)) {
      x <- as.numeric(rbind(edges$ap1, edges$ap2, NA_real_))
      y <- as.numeric(rbind(edges$ml1, edges$ml2, NA_real_))
      z <- as.numeric(rbind(edges$dv1, edges$dv2, NA_real_))
      p <- add_trace(
        p,
        type = "scatter3d",
        mode = "lines",
        x = x,
        y = y,
        z = z,
        name = paste(info$label[[1]], "wireframe"),
        line = list(color = info$color[[1]], width = 1.15),
        hoverinfo = "skip",
        visible = FALSE,
        showlegend = FALSE,
        meta = list(
          role = "anatomy",
          kind = "wire",
          key = key,
          level = info$level[[1]],
          group = info$group[[1]],
          base_opacity = 0.55
        )
      )
    }
  }
  p
}

build_viewer <- function(
  intruder,
  baseline,
  borders,
  seqfish,
  mesh,
  ap_min,
  ap_max,
  ne_min,
  ne_max
) {
  set.seed(SET_SEED)
  p <- plot_ly(
    width = NULL,
    height = NULL,
    colors = rep("#7f8791", 256)
  )

  intruder_samples <- nested_samples(intruder, N_MAP)
  baseline_samples <- nested_samples(baseline, N_MAP)
  seqfish_samples <- nested_samples(seqfish, N_SEQ)

  for (density in names(N_MAP)) {
    intr_sample <- intruder_samples[[density]]
    base_sample <- baseline_samples[[density]]
    p <- add_2dg_trace(
      p,
      intr_sample,
      "intruder",
      density,
      "full",
      initially_visible = density == "medium"
    )
    p <- add_2dg_trace(p, base_sample, "baseline", density, "full")
    p <- add_2dg_trace(
      p,
      intr_sample[ccfx >= ap_min & ccfx <= ap_max],
      "intruder",
      density,
      "ap"
    )
    p <- add_2dg_trace(
      p,
      base_sample[ccfx >= ap_min & ccfx <= ap_max],
      "baseline",
      density,
      "ap"
    )
  }
  for (density in names(N_SEQ)) {
    p <- add_seqfish_trace(
      p,
      seqfish_samples[[density]],
      density,
      ne_min,
      ne_max
    )
  }
  p <- add_border_trace(p, borders)
  p <- add_region_traces(p, mesh)

  p <- layout(
    p,
    scene = scene_layout(),
    margin = list(l = 0, r = 8, b = 0, t = 4),
    showlegend = FALSE,
    paper_bgcolor = "#ffffff",
    plot_bgcolor = "#ffffff",
    hoverlabel = list(
      bgcolor = "#ffffff",
      bordercolor = "#d8dce2",
      font = list(color = "#17191c", size = 12)
    )
  )
  config(
    p,
    displaylogo = FALSE,
    responsive = TRUE,
    scrollZoom = TRUE,
    modeBarButtonsToRemove = c(
      "toImage", "sendDataToCloud", "select3d", "lasso3d", "resetCameraLastSave3d"
    )
  )
}

region_catalog <- function(regions) {
  domains <- regions[level == "domain", .(
    key,
    label,
    color
  )]
  subdomains <- regions[level == "subdomain", .(
    key,
    label,
    group,
    color
  )]
  list(domains = domains, subdomains = subdomains)
}

viewer_controls_html <- function(catalog, ap_min, ap_max, ne_min, ne_max) {
  catalog_json <- toJSON(catalog, dataframe = "rows", auto_unbox = TRUE)
  paste0(
    "<header class='viewer-head'>",
    "<div><h1>2DG × seqFISH · Allen CCF explorer</h1>",
    "<span>Independent data and anatomy layers</span></div>",
    "<div class='head-actions'><div class='viewer-help'>Drag to rotate · wheel to zoom · hover for values</div>",
    "<button id='panel-toggle' type='button' aria-expanded='true'>Hide controls</button></div>",
    "</header>",
    "<aside class='viewer-controls' aria-label='Viewer controls'>",
    "<section class='data-controls'><h2>Data layers</h2>",
    "<label>2DG contrast<select id='twodg-layer'>",
    "<option value='off'>None</option>",
    "<option value='intruder' selected>Intruder WildR−SPF</option>",
    "<option value='baseline'>Baseline WildR−SPF</option>",
    "</select></label>",
    "<label class='check'><input id='layer-seqfish' type='checkbox'>seqFISH NE activity</label>",
    "<div class='two-col'>",
    "<label>Density<select id='density'><option value='low'>Low</option>",
    "<option value='medium' selected>Medium</option><option value='high'>High</option></select></label>",
    "<label>Spatial scope<select id='scope'><option value='full'>Full 2DG volume</option>",
    "<option value='ap'>seqFISH AP band</option></select></label></div>",
    "<label>Point size<input id='point-size' type='range' min='0.7' max='1.8' step='0.1' value='1'></label>",
    "</section>",
    "<section class='anatomy-controls'><h2>Anatomical context</h2>",
    "<label class='check'><input id='root-context' type='checkbox' checked>Whole-brain context</label>",
    "<label>Large domain<select id='domain'></select></label>",
    "<label>Subdomain / region<select id='subdomain'><option value='off'>None</option></select></label>",
    "<div class='check-grid'>",
    "<label class='check'><input id='domain-on' type='checkbox' checked>Selected domain</label>",
    "<label class='check'><input id='surface-on' type='checkbox' checked>Surface</label>",
    "<label class='check'><input id='wire-on' type='checkbox'>Wire mesh</label>",
    "<label class='check'><input id='border-on' type='checkbox'>Parcel points</label>",
    "</div>",
    "<label>Region opacity<input id='region-opacity' type='range' min='0.08' max='0.85' step='0.01' value='0.34'></label>",
    "<div class='region-key'><span id='region-swatch'></span><span id='region-name'>Amygdala</span></div>",
    "</section>",
    "<section class='view-controls'><h2>View</h2>",
    "<div class='two-col'><label>Camera<select id='camera'>",
    "<option value='three-quarter'>Three-quarter</option>",
    "<option value='sagittal'>Sagittal</option>",
    "<option value='frontal'>Frontal</option>",
    "<option value='dorsal'>Dorsal</option></select></label>",
    "<button id='reset-view' type='button'>Reset view</button></div>",
    "</section>",
    "<section class='legend-block color-controls'><h2>Color</h2>",
    "<div id='legend-2dg'><div class='legend-title'>2DG t · WildR − SPF</div>",
    "<div class='gradient effect'></div><div class='legend-labels'><span>SPF↑ −4</span><span>0</span><span>WildR↑ +4</span></div></div>",
    "<div id='legend-ne' class='hidden'><div class='legend-title'>NEUROeSTIMator activity</div>",
    "<div class='gradient viridis'></div><div class='legend-labels'><span>",
    sprintf("%.2f", ne_min), "</span><span>", sprintf("%.2f", ne_max), "</span></div></div>",
    "</section>",
    "<p class='scope-note'>AP-matched scope: ", round(ap_min), "–", round(ap_max),
    " µm. Blue 2DG points indicate WildR-higher t values; red indicates SPF-higher.</p>",
    "</aside>",
    "<aside class='floating-legend' aria-label='Color legend'>",
    "<div id='floating-legend-2dg'><div class='legend-title'>2DG t · WildR − SPF</div>",
    "<div class='gradient effect'></div><div class='legend-labels'><span>SPF↑ −4</span><span>0</span><span>WildR↑ +4</span></div></div>",
    "<div id='floating-legend-ne' class='hidden'><div class='legend-title'>NEUROeSTIMator activity</div>",
    "<div class='gradient viridis'></div><div class='legend-labels'><span>",
    sprintf("%.2f", ne_min), "</span><span>", sprintf("%.2f", ne_max), "</span></div></div>",
    "</aside>",
    "<script>window.REGION_CATALOG=", catalog_json, ";</script>"
  )
}

viewer_css <- function() {
  paste0(
    "<style>",
    ":root{--ink:#17191c;--muted:#66707b;--line:#dfe3e8;--soft:#f6f7f8;--accent:#ff6c0c;}",
    "*{box-sizing:border-box}html,body{margin:0!important;padding:0!important;width:100%;height:100%;",
    "overflow:hidden!important;background:#fff;color:var(--ink);font-family:Inter,Helvetica,Arial,sans-serif}",
    ".viewer-head{position:fixed;z-index:50;inset:0 0 auto 0;height:60px;padding:10px 18px;",
    "display:flex;justify-content:space-between;align-items:center;background:rgba(255,255,255,.97);",
    "border-bottom:1px solid var(--line)}",
    ".viewer-head h1{margin:0;font-size:18px;line-height:1.15;letter-spacing:-.015em}.viewer-head span,.viewer-help{font-size:12px;color:var(--muted)}",
    ".head-actions{display:flex;align-items:center;gap:12px}#panel-toggle{border:1px solid #cfd5dc;border-radius:7px;",
    "background:#fff;color:var(--ink);padding:7px 10px;font-size:12px;font-weight:650;cursor:pointer}",
    ".viewer-controls{position:fixed;z-index:45;left:0;top:60px;bottom:0;width:318px;padding:12px 14px 14px;",
    "overflow:hidden;background:#fff;border-right:1px solid var(--line);transition:transform .18s ease;",
    "display:flex;flex-direction:column}",
    ".viewer-controls section{flex:0 0 auto;padding:0 0 10px;margin:0 0 10px;border-bottom:1px solid var(--line)}",
    ".viewer-controls .anatomy-controls{flex:1 1 auto;min-height:120px;overflow-x:hidden;overflow-y:auto}",
    ".viewer-controls h2{font-size:12px;text-transform:uppercase;letter-spacing:.08em;margin:0 0 8px;color:#4d5661}",
    ".viewer-controls label{display:block;font-size:12.5px;color:#4d5661;margin:7px 0}",
    ".viewer-controls select,.viewer-controls button{display:block;width:100%;margin-top:4px;border:1px solid #cfd5dc;",
    "border-radius:7px;background:#fff;color:var(--ink);padding:7px 8px;font-size:13px}",
    ".viewer-controls button{align-self:end;cursor:pointer;font-weight:650}.viewer-controls button:hover{border-color:var(--accent)}",
    ".check{display:flex!important;gap:7px;align-items:center;color:var(--ink)!important}.check input{accent-color:var(--accent)}",
    ".check-grid{display:grid;grid-template-columns:1fr 1fr}.two-col{display:grid;grid-template-columns:1fr 1fr;gap:8px}",
    "input[type=range]{width:100%;accent-color:var(--accent)}",
    ".region-key{display:flex;align-items:center;gap:8px;margin-top:8px;font-size:12.5px;font-weight:650}",
    "#region-swatch{width:12px;height:12px;border-radius:3px;background:#d1495b;border:1px solid rgba(0,0,0,.12)}",
    ".legend-block{display:grid;gap:10px}.legend-title{font-size:12px;margin-bottom:4px}.gradient{height:9px;border-radius:99px;border:1px solid #d8dce2}",
    ".gradient.effect{background:linear-gradient(90deg,#b2182b,#f7f7f7 50%,#2166ac)}",
    ".gradient.viridis{background:linear-gradient(90deg,#440154,#31688e,#35b779,#fde725)}",
    ".legend-labels{display:flex;justify-content:space-between;font-size:11px;color:var(--muted);margin-top:2px}",
    ".scope-note{font-size:11px;line-height:1.45;color:var(--muted);margin:0}.hidden{display:none}",
    ".floating-legend{display:none;position:fixed;z-index:48;left:14px;bottom:14px;width:290px;",
    "gap:10px;padding:11px 12px;background:rgba(255,255,255,.96);border:1px solid var(--line);",
    "border-radius:9px;box-shadow:0 4px 16px rgba(30,40,55,.10)}",
    "#htmlwidget_container{position:fixed!important;left:318px!important;right:0!important;top:60px!important;bottom:0!important;",
    "width:auto!important;height:auto!important;overflow:hidden!important}",
    "#htmlwidget_container>.html-widget,.plotly.html-widget,.js-plotly-plot{width:100%!important;height:100%!important}",
    "body.controls-collapsed .viewer-controls{transform:translateX(-100%)}",
    "body.controls-collapsed #htmlwidget_container{left:0!important}",
    "body.controls-collapsed .floating-legend{display:grid}",
    "@media(max-width:820px){.viewer-controls{width:min(318px,86vw)}.two-col{grid-template-columns:1fr}",
    ".viewer-help{display:none}.viewer-head{padding:8px 10px}.viewer-head h1{font-size:16px}}",
    "</style>"
  )
}

viewer_javascript <- function() {
  paste0(
    "<script>",
    "(function(){",
    "const $=id=>document.getElementById(id);",
    "let initialized=false;",
    "function wait(){if(initialized)return;const gd=document.querySelector('.js-plotly-plot');",
    "if(!gd||!window.Plotly){setTimeout(wait,60);return;}initialized=true;init(gd);}",
    "function init(gd){",
    "const C=window.REGION_CATALOG;",
    "const domain=$('domain'),sub=$('subdomain');",
    "domain.innerHTML='<option value=\"off\">None</option>'+C.domains.map(d=>`<option value=\"${d.key}\" data-color=\"${d.color}\">${d.label}</option>`).join('');",
    "const amy=C.domains.findIndex(d=>d.label==='Amygdala');domain.selectedIndex=amy>=0?amy+1:0;",
    "function selectedDomain(){return C.domains.find(d=>d.key===domain.value)}",
    "function rebuildSub(){const d=selectedDomain();const rows=d?C.subdomains.filter(s=>s.group===d.label):[];",
    "sub.innerHTML='<option value=\"off\">None</option>'+rows.map(s=>`<option value=\"${s.key}\" data-color=\"${s.color}\">${s.label}</option>`).join('');",
    "$('region-name').textContent=d?d.label:'No selected region';$('region-swatch').style.background=d?d.color:'#c8c8cf';}",
    "function camera(which){const eyes={",
    "'three-quarter':{x:1.55,y:1.45,z:.85},sagittal:{x:.05,y:2.25,z:.15},frontal:{x:2.25,y:.05,z:.15},dorsal:{x:.05,y:.05,z:2.35}};",
    "Plotly.relayout(gd,{'scene.camera':{eye:eyes[which],center:{x:0,y:0,z:0},up:{x:0,y:0,z:1}}});}",
    "function apply(){",
    "const twodg=$('twodg-layer').value;",
    "const layers={intruder:twodg==='intruder',baseline:twodg==='baseline',seqfish:$('layer-seqfish').checked};",
    "const density=$('density').value,scope=$('scope').value,domainKey=domain.value,subKey=sub.value;",
    "const rootOn=$('root-context').checked,domainOn=$('domain-on').checked&&domainKey!=='off';",
    "const surfaceOn=$('surface-on').checked,wireOn=$('wire-on').checked,borderOn=$('border-on').checked;",
    "const vis=[],anatomyIdx=[],anatomyOpacity=[];const pointScale=Number($('point-size').value),regionAlpha=Number($('region-opacity').value);",
    "gd.data.forEach((tr,i)=>{const m=tr.meta||{};let on=false;",
    "if(m.role==='data'){on=!!layers[m.dataset]&&m.density===density&&(m.dataset==='seqfish'||m.scope===scope);}",
    "else if(m.role==='border'){on=borderOn;}",
    "else if(m.role==='anatomy'){const kindOn=m.kind==='surface'?surfaceOn:wireOn;",
    "on=kindOn&&((m.key==='root'&&rootOn)||(domainOn&&m.key===domainKey)||(subKey!=='off'&&m.key===subKey));",
    "anatomyIdx.push(i);let a=m.key==='root' ? 0.09 : (m.key===subKey ? Math.min(0.95,regionAlpha+0.28) : regionAlpha);",
    "if(m.kind==='wire')a=Math.min(.9,a+.25);anatomyOpacity.push(a);}",
    "vis.push(on);});",
    "Plotly.restyle(gd,{visible:vis});",
    "if(anatomyIdx.length)Plotly.restyle(gd,{opacity:anatomyOpacity},anatomyIdx);",
    "const dataIdx=[];const sizes=[];gd.data.forEach((tr,i)=>{if((tr.meta||{}).role==='data'){dataIdx.push(i);sizes.push((tr.meta.base_size||2)*pointScale);}});",
    "if(dataIdx.length)Plotly.restyle(gd,{'marker.size':sizes},dataIdx);",
    "[$('legend-2dg'),$('floating-legend-2dg')].forEach(el=>el.classList.toggle('hidden',twodg==='off'));",
    "[$('legend-ne'),$('floating-legend-ne')].forEach(el=>el.classList.toggle('hidden',!layers.seqfish));",
    "const active=subKey!=='off'?C.subdomains.find(s=>s.key===subKey):(domainOn?selectedDomain():null);",
    "$('region-name').textContent=active?active.label:'No selected region';$('region-swatch').style.background=active?active.color:'#c8c8cf';",
    "setTimeout(()=>Plotly.Plots.resize(gd),0);}",
    "rebuildSub();",
    "['twodg-layer','layer-seqfish','density','scope','point-size','root-context','domain-on','surface-on','wire-on','border-on','region-opacity']",
    ".forEach(id=>$(id).addEventListener('input',apply));",
    "domain.addEventListener('change',()=>{rebuildSub();apply();});sub.addEventListener('change',apply);",
    "$('camera').addEventListener('change',e=>camera(e.target.value));",
    "$('reset-view').addEventListener('click',()=>{camera('three-quarter');$('camera').value='three-quarter';});",
    "$('panel-toggle').addEventListener('click',()=>{const collapsed=document.body.classList.toggle('controls-collapsed');",
    "$('panel-toggle').textContent=collapsed?'Show controls':'Hide controls';$('panel-toggle').setAttribute('aria-expanded',String(!collapsed));",
    "setTimeout(()=>Plotly.Plots.resize(gd),190);});",
    "window.addEventListener('resize',()=>Plotly.Plots.resize(gd));apply();camera('three-quarter');",
    "}",
    "window.addEventListener('load',wait);setTimeout(wait,80);",
    "})();",
    "</script>"
  )
}

validate_viewer <- function(widget) {
  built <- plotly_build(widget)
  meta <- lapply(built$x$data, function(trace) trace$meta)
  role <- vapply(meta, function(x) as.character(x$role), character(1))
  role_counts <- table(role)
  if (!identical(
    as.integer(role_counts[c("anatomy", "border", "data")]),
    c(180L, 1L, 15L)
  )) {
    stop("Unexpected viewer trace-role counts.", call. = FALSE)
  }
  data_meta <- meta[role == "data"]
  key <- vapply(meta, function(x) if (is.null(x$key)) "" else as.character(x$key), character(1))
  required_keys <- c("root", "domain__amygdala", "domain__hypothalamus", "region__PVH")
  if (!all(required_keys %in% key)) {
    stop("Required anatomy traces are missing.", call. = FALSE)
  }
  combinations <- vapply(
    data_meta,
    function(x) paste(x$dataset, x$density, x$scope, sep = ":"),
    character(1)
  )
  expected <- c(
    as.vector(outer(
      c("intruder", "baseline"),
      as.vector(outer(names(N_MAP), c("full", "ap"), paste, sep = ":")),
      paste,
      sep = ":"
    )),
    paste("seqfish", names(N_SEQ), "both", sep = ":")
  )
  if (!setequal(combinations, expected)) {
    stop("Viewer data trace matrix is incomplete.", call. = FALSE)
  }
  invisible(TRUE)
}

atomic_move_new <- function(temp_path, final_path) {
  if (file.exists(final_path)) {
    stop("Refusing to overwrite: ", final_path, call. = FALSE)
  }
  if (!file.rename(temp_path, final_path)) {
    stop("Atomic move failed: ", temp_path, " -> ", final_path, call. = FALSE)
  }
}

write_viewer_page <- function(
  widget,
  catalog,
  ap_min,
  ap_max,
  ne_min,
  ne_max,
  path
) {
  temp_path <- file.path(
    dirname(path),
    paste0(".", tools::file_path_sans_ext(basename(path)), ".tmp-", Sys.getpid(), ".html")
  )
  if (file.exists(temp_path)) {
    stop("Temporary viewer path already exists: ", temp_path, call. = FALSE)
  }
  saveWidget(
    widget,
    temp_path,
    selfcontained = TRUE,
    title = "2DG CCF modular explorer"
  )
  html <- paste(readLines(temp_path, warn = FALSE), collapse = "\n")
  head_inject <- viewer_css()
  body_inject <- viewer_controls_html(catalog, ap_min, ap_max, ne_min, ne_max)
  js_inject <- viewer_javascript()
  html <- sub("</head>", paste0(head_inject, "\n</head>"), html, ignore.case = TRUE)
  html <- sub("(<body[^>]*>)", paste0("\\1\n", body_inject), html, ignore.case = TRUE)
  html <- sub("</body>", paste0(js_inject, "\n</body>"), html, ignore.case = TRUE)
  writeLines(html, temp_path)
  atomic_move_new(temp_path, path)
}

main <- function() {
  outputs <- c(VIEWER_HTML, VIEWER_NOTE)
  stop_if_exists(outputs)
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
  set.seed(SET_SEED)

  log_message("Loading 2DG map points")
  intruder <- load_map_points(INTRUDER_OVERLAY)
  baseline <- load_map_points(BASELINE_OVERLAY)
  borders <- fread(BORDER_POINTS)

  log_message("Loading seqFISH points with anatomical ML coordinates")
  joined <- as.data.table(readRDS(JOINED))
  seqfish <- joined[
    registered_ccf == TRUE &
      is.finite(ccfx) &
      is.finite(ccfy) &
      is.finite(ccfz) &
      is.finite(predicted_activity),
    .(
      ccfx,
      ccfy,
      ccfz,
      predicted_activity,
      group = as.character(group),
      final_gate = as.character(final_gate),
      sample_id = as.character(sample_id)
    )
  ]
  if (!nrow(seqfish)) {
    stop("No finite registered seqFISH points are available.", call. = FALSE)
  }
  seqfish <- apply_anatomical_ml(seqfish)
  ap_min <- as.numeric(min(seqfish$ccfx, na.rm = TRUE))
  ap_max <- as.numeric(max(seqfish$ccfx, na.rm = TRUE))
  ne_limits <- as.numeric(
    quantile(seqfish$predicted_activity, c(0.01, 0.99), na.rm = TRUE)
  )
  if (length(ne_limits) != 2L || any(!is.finite(ne_limits)) || ne_limits[1] >= ne_limits[2]) {
    stop("Invalid global seqFISH color limits.", call. = FALSE)
  }
  ne_min <- ne_limits[[1]]
  ne_max <- ne_limits[[2]]

  log_message("Loading curated Allen domain and subdomain meshes")
  mesh <- load_mesh_catalog()
  catalog <- region_catalog(mesh$regions)

  log_message("Building modular viewer")
  viewer <- build_viewer(
    intruder,
    baseline,
    borders,
    seqfish,
    mesh,
    ap_min,
    ap_max,
    ne_min,
    ne_max
  )
  validate_viewer(viewer)
  write_viewer_page(
    viewer,
    catalog,
    ap_min,
    ap_max,
    ne_min,
    ne_max,
    VIEWER_HTML
  )

  note <- list(
    analysis_date = ANALYSIS_DATE,
    polarity = list(
      signed_maps = "WildR minus SPF",
      positive_or_blue = "WildR higher",
      negative_or_red = "SPF higher"
    ),
    controls = list(
      two_dg_selector = c("None", "Intruder", "Baseline"),
      seqfish_independent = TRUE,
      shared_density = names(N_MAP),
      shared_2dg_scope = c("full", "seqFISH AP band"),
      anatomy = c(
        "whole-brain context",
        "selected domain visibility",
        "subdomain",
        "surface",
        "wire mesh",
        "parcel points"
      )
    ),
    sampling = list(map = as.list(N_MAP), seqfish = as.list(N_SEQ), seed = SET_SEED),
    seqfish_ap_window_um = c(ap_min, ap_max),
    seqfish_color_limits_1_99pct = c(ne_min, ne_max),
    mesh_manifest = MESH_MANIFEST,
    n_mesh_regions = nrow(mesh$regions)
  )
  note_temp <- file.path(
    dirname(VIEWER_NOTE),
    paste0(".", basename(VIEWER_NOTE), ".tmp-", Sys.getpid())
  )
  write_json(note, note_temp, pretty = TRUE, auto_unbox = TRUE)
  atomic_move_new(note_temp, VIEWER_NOTE)
  log_message("Wrote ", VIEWER_HTML)
}

main()
