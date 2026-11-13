#!/usr/bin/env Rscript

# Build an interactive Plotly 3D viewer: Intruder t, optional CCF borders,
# optional seqFISH NEUROeSTIMator cells.
#
# /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript \
#   notebooks/R_scripts/plot_2dg_ccf_interactive3d.R

suppressPackageStartupMessages({
  library(data.table)
  library(plotly)
  library(htmlwidgets)
})

WKDIR <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
SOURCE_DATE <- "2026-11-11"
ANALYSIS_DATE <- "2026-11-12"

SOURCE_DIR <- file.path(WKDIR, "data/interim/2DG_CCF_Integration", SOURCE_DATE)
OUT_DIR <- file.path(WKDIR, "data/interim/2DG_CCF_Integration", ANALYSIS_DATE)
FIGURE_DIR <- file.path(WKDIR, "figures/2DG_CCF_Integration")

INTRUDER_OVERLAY <- file.path(
  SOURCE_DIR, paste0("2dg_overlay_points_ccf_scalar-t_", SOURCE_DATE, ".csv.gz")
)
BORDER_POINTS <- file.path(
  OUT_DIR, paste0("allen_ccf_border_points_subsampled_", ANALYSIS_DATE, ".csv.gz")
)
JOINED <- file.path(
  SOURCE_DIR, paste0("seqFISH_neuroestimator_ccf_strict_", SOURCE_DATE, ".rds")
)

WIDGET_HTML <- file.path(
  FIGURE_DIR, paste0("2DG_CCF_interactive3d_", ANALYSIS_DATE, ".html")
)
WIDGET_RDS <- file.path(
  OUT_DIR, paste0("2DG_CCF_interactive3d_widget_", ANALYSIS_DATE, ".rds")
)

N_INTRUDER <- 22000L
N_SEQFISH <- 12000L
SET_SEED <- 42L

log_message <- function(...) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
}

stop_if_exists <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing)) {
    stop("Refusing to overwrite:\n", paste(existing, collapse = "\n"), call. = FALSE)
  }
}

main <- function() {
  stop_if_exists(c(WIDGET_HTML, WIDGET_RDS))
  if (!file.exists(BORDER_POINTS)) {
    stop("Border point cloud missing. Run export_ccf_border_points.py first.", call. = FALSE)
  }
  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  set.seed(SET_SEED)

  log_message("Loading Intruder overlay")
  intruder <- fread(
    INTRUDER_OVERLAY,
    select = c("ccfx", "ccfy", "ccfz", "t_value")
  )
  intruder <- intruder[is.finite(t_value) & is.finite(ccfx)]
  if (nrow(intruder) > N_INTRUDER) {
    intruder <- intruder[sample.int(.N, N_INTRUDER)]
  }
  intruder[, t_clip := pmax(pmin(t_value, 4), -4)]

  log_message("Loading CCF border points")
  borders <- fread(BORDER_POINTS)

  log_message("Loading seqFISH subsample")
  joined <- as.data.table(readRDS(JOINED))
  seqfish <- joined[
    registered_ccf == TRUE &
      is.finite(ccfx) &
      is.finite(predicted_activity)
  ][
    ,
    .(ccfx, ccfy, ccfz, predicted_activity, group)
  ]
  if (nrow(seqfish) > N_SEQFISH) {
    seqfish <- seqfish[sample.int(.N, N_SEQFISH)]
  }

  scene <- list(
    xaxis = list(title = "CCF AP (µm)"),
    yaxis = list(title = "CCF ML (µm)"),
    zaxis = list(title = "CCF DV (µm)", autorange = "reversed"),
    aspectmode = "data"
  )

  p <- plot_ly()
  p <- add_trace(
    p,
    data = intruder,
    x = ~ccfx,
    y = ~ccfz,
    z = ~ccfy,
    type = "scatter3d",
    mode = "markers",
    name = "Intruder t",
    marker = list(
      size = 1.6,
      opacity = 0.55,
      color = ~t_clip,
      colorscale = "RdBu",
      reversescale = TRUE,
      cmin = -4,
      cmax = 4,
      colorbar = list(title = "t")
    ),
    text = ~sprintf("t=%.2f", t_value),
    hoverinfo = "text+name",
    visible = TRUE
  )
  p <- add_trace(
    p,
    data = borders,
    x = ~ccfx,
    y = ~ccfz,
    z = ~ccfy,
    type = "scatter3d",
    mode = "markers",
    name = "CCF borders",
    marker = list(size = 1.1, opacity = 0.18, color = "#444444"),
    hoverinfo = "skip",
    visible = "legendonly"
  )
  p <- add_trace(
    p,
    data = seqfish,
    x = ~ccfx,
    y = ~ccfz,
    z = ~ccfy,
    type = "scatter3d",
    mode = "markers",
    name = "seqFISH NEUROeSTIMator",
    marker = list(
      size = 1.4,
      opacity = 0.35,
      color = ~predicted_activity,
      colorscale = "Viridis",
      colorbar = list(title = "NE", x = 1.12)
    ),
    text = ~sprintf("%s<br>NE=%.3f", group, predicted_activity),
    hoverinfo = "text+name",
    visible = "legendonly"
  )

  p <- layout(
    p,
    title = list(
      text = "Interactive CCF co-embedding: Intruder t ± CCF ± seqFISH",
      y = 0.98
    ),
    scene = scene,
    legend = list(orientation = "h", y = -0.05),
    margin = list(l = 0, r = 0, b = 40, t = 50),
    updatemenus = list(
      list(
        type = "buttons",
        direction = "right",
        x = 0.02,
        y = 1.12,
        buttons = list(
          list(
            label = "Intruder",
            method = "update",
            args = list(list(visible = c(TRUE, FALSE, FALSE)))
          ),
          list(
            label = "+ CCF",
            method = "update",
            args = list(list(visible = c(TRUE, TRUE, FALSE)))
          ),
          list(
            label = "+ seqFISH",
            method = "update",
            args = list(list(visible = c(TRUE, FALSE, TRUE)))
          ),
          list(
            label = "All",
            method = "update",
            args = list(list(visible = c(TRUE, TRUE, TRUE)))
          )
        )
      )
    )
  )
  p <- config(p, displaylogo = FALSE, responsive = TRUE)

  log_message("Saving widget HTML: ", WIDGET_HTML)
  saveWidget(p, WIDGET_HTML, selfcontained = TRUE, title = "2DG CCF interactive 3D")
  saveRDS(p, WIDGET_RDS, compress = "xz")
  log_message("Done")
}

main()
