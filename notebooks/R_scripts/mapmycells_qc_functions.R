# MapMyCells mapping-quality QC: WildR vs SPF per cell type
# ---------------------------------------------------------------------------
# Used by notebooks/R_scripts/run_mapmycells_qc.R and
# notebooks/0X_MapMyCells_Mapping_QC.qmd.
#
# Inputs : data/interim/MapMyCells_QC/mapmycells_scores_{SN,seqFISH}_<date>.parquet
#          (python_scripts/extract_mapmycells_scores.py)
# Design : replicate = animal (snRNA, analysed within region) or ROI section
#          (seqFISH; ROI -> animal key not recorded, see notebook caveat).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(glue)
  library(ggplot2)
  library(patchwork)
})

# Group colours: WildR orange / SPF blue (repo convention; validated pair,
# CVD dE >= 20 and >= 3:1 contrast on white).
mmc_group_cols <- c(WildR = "#E8590C", SPF = "#2E6FA8")
# Diverging fill for effect sizes: SPF-higher (blue) <-> 0 (gray) <-> WildR-higher (orange)
mmc_div_cols <- c(low = "#2E6FA8", mid = "#F0EFEC", high = "#E8590C")

mmc_levels <- c("class", "subclass", "supertype", "cluster")

# Endpoints: replicate-level summaries of per-cell scores at the level (L)
# that defines the cell type.
mmc_endpoints <- tibble::tribble(
  ~endpoint,       ~cell_metric,    ~label,                                   ~role,
  "mean_boot",     "boot_prob",     "Bootstrap probability (mean)",           "primary",
  "frac_low_conf", "low_conf",      "Fraction low-confidence (boot < 0.5)",   "primary",
  "mean_corr",     "avg_corr",      "Correlation to reference (mean)",        "primary"
)
# Not used as an endpoint: correlation margin over the runner-up. MapMyCells
# only reports runner-ups when >= 1 bootstrap iteration disagreed, so the
# margin is missing for every cell with boot_prob == 1 (77% SN, 24% seqFISH):
# missing-not-at-random by construction.

theme_mmc <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#ECECEC", linewidth = 0.3),
      axis.text = element_text(colour = "#555555"),
      axis.title = element_text(colour = "#555555"),
      plot.title = element_text(face = "bold", colour = "#0A0A0A"),
      strip.text = element_text(face = "bold", colour = "#0A0A0A", hjust = 0),
      legend.position = "top",
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

# Figure sizes (inches). Every figure is saved 12 in wide and displayed at the
# full ~1,220 px content width (.fig-wide cards), so a 13 pt base font renders
# at ~16-18 px and all figures share one text scale. Height is the only knob.
mmc_fig_w <- 12
mmc_fig_dpi <- 200
mmc_fig_h <- c(
  ecdf = 8, depth = 4.6, landscape = 12, heatmap_class = 10, heatmap = 13,
  forest = 10, pvalue = 9.5, superplot = 6.5
)

#' Save a QC figure at the shared width / dpi.
mmc_save <- function(p, fig_dir, file, kind) {
  ggsave(file.path(fig_dir, file), p, width = mmc_fig_w, height = mmc_fig_h[[kind]],
         dpi = mmc_fig_dpi, bg = "white")
  invisible(file.path(fig_dir, file))
}

#' Load one modality's parquet and harmonise to a common long-ish frame.
#' Adds `stratum` (SN AMY / SN HYP / seqFISH), `replicate`, `log10_counts`.
load_mmc_scores <- function(path, modality, min_counts = 0) {
  df <- arrow::read_parquet(path)
  df <- df %>%
    mutate(
      modality = modality,
      group = factor(group, levels = c("SPF", "WildR")),
      log10_counts = log10(pmax(n_counts, 1))
    )
  if (modality == "SN") {
    df <- df %>% mutate(stratum = paste("SN", region), replicate = animal, run = NA_character_)
  } else {
    df <- df %>% mutate(stratum = "seqFISH", replicate = sample_id, region = NA_character_)
  }
  n0 <- nrow(df)
  df <- df %>% filter(n_counts >= min_counts)
  message(glue("{modality}: {n0 - nrow(df)} of {n0} cells dropped by min_counts = {min_counts}"))
  df
}

#' Per-cell metrics for the level that defines the cell type.
cell_metrics_at <- function(df, level) {
  df %>%
    transmute(
      modality, stratum, region, replicate, sample_id, run, group, log10_counts,
      level = level,
      cell_type = .data[[glue("{level}_name")]],
      class_name,
      boot_prob = .data[[glue("{level}_boot_prob")]],
      avg_corr = .data[[glue("{level}_avg_corr")]],
      low_conf = as.numeric(.data[[glue("{level}_boot_prob")]] < 0.5)
    )
}

#' Replicate-level (pseudobulk) summaries per cell type.
summarise_replicates <- function(cm) {
  cm %>%
    group_by(modality, stratum, level, class_name, cell_type, replicate, run, group) %>%
    summarise(
      n_cells = n(),
      median_log10_counts = median(log10_counts),
      mean_boot = mean(boot_prob),
      frac_low_conf = mean(low_conf),
      mean_corr = mean(avg_corr),
      .groups = "drop"
    )
}

#' Cell types testable in a stratum: >= min_cells in >= min_rep_frac of the
#' replicates of *each* group (SN: 3 of 6 animals; seqFISH: 4 of 8 ROIs).
eligible_types <- function(rep_df, min_cells = 20, min_rep_frac = 0.5) {
  n_tot <- rep_df %>% distinct(stratum, group, replicate) %>%
    count(stratum, group, name = "n_total")
  rep_df %>%
    filter(n_cells >= min_cells) %>%
    count(modality, stratum, level, cell_type, group, name = "n_rep") %>%
    left_join(n_tot, by = c("stratum", "group")) %>%
    mutate(ok = n_rep >= ceiling(min_rep_frac * n_total)) %>%
    group_by(modality, stratum, level, cell_type) %>%
    filter(n() == 2, all(ok)) %>%
    ungroup() %>%
    select(-n_total, -ok) %>%
    pivot_wider(names_from = group, values_from = n_rep, values_fill = 0,
                names_prefix = "n_rep_")
}

cliffs_delta <- function(x, y) {
  # P(x > y) - P(x < y); x = WildR, y = SPF
  d <- outer(x, y, "-")
  (sum(d > 0) - sum(d < 0)) / length(d)
}

hedges_g <- function(x, y) {
  nx <- length(x); ny <- length(y)
  sp <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  j <- 1 - 3 / (4 * (nx + ny) - 9)
  j * (mean(x) - mean(y)) / sp
}

#' Primary test: linear model on replicate means (pseudobulk). seqFISH ROIs are
#' blocked by imaging run (each run = 2 SPF + 2 WildR; run 1 is a low-depth
#' batch), so the model is value ~ group + run (classical pooled-variance t for snRNA, balanced 6 v 6).
#' Robustness: exact Wilcoxon on replicate means (seqFISH: after removing run
#' means) + Cliff's delta. BH-FDR within modality x stratum x level x endpoint.
test_replicate_level <- function(rep_df, elig, min_cells = 20) {
  rep_ok <- rep_df %>%
    filter(n_cells >= min_cells) %>%
    semi_join(elig, by = c("modality", "stratum", "level", "cell_type"))

  rep_ok %>%
    pivot_longer(all_of(mmc_endpoints$endpoint), names_to = "endpoint", values_to = "value") %>%
    group_by(modality, stratum, level, class_name, cell_type, endpoint) %>%
    group_modify(function(d, key) {
      blocked <- key$stratum == "seqFISH" && n_distinct(d$run) > 1
      f <- if (blocked) value ~ group + run else value ~ group
      fit <- tryCatch(lm(f, data = d), error = function(e) NULL)
      co <- if (!is.null(fit)) summary(fit)$coefficients else NULL
      ok <- !is.null(co) && "groupWildR" %in% rownames(co) && fit$df.residual > 0
      ci <- if (ok) suppressWarnings(confint(fit, "groupWildR")) else c(NA, NA)
      # run-centred values for the rank-based checks
      v <- if (blocked) d$value - ave(d$value, d$run) else d$value
      x <- v[d$group == "WildR"]; y <- v[d$group == "SPF"]
      wx <- suppressWarnings(tryCatch(wilcox.test(x, y, exact = TRUE),
                                      error = function(e) NULL))
      tibble(
        model = if (blocked) "value ~ group + run" else "value ~ group",
        n_rep_WildR = length(x), n_rep_SPF = length(y),
        n_cells_WildR = sum(d$n_cells[d$group == "WildR"]),
        n_cells_SPF = sum(d$n_cells[d$group == "SPF"]),
        mean_WildR = mean(d$value[d$group == "WildR"]),
        mean_SPF = mean(d$value[d$group == "SPF"]),
        delta = if (ok) co["groupWildR", "Estimate"] else NA_real_,
        ci_low = ci[1], ci_high = ci[2],
        df_resid = if (ok) fit$df.residual else NA_real_,
        p_lm = if (ok) co["groupWildR", "Pr(>|t|)"] else NA_real_,
        hedges_g = hedges_g(x, y),
        cliffs_delta = cliffs_delta(x, y),
        p_wilcox = if (is.null(wx)) NA_real_ else wx$p.value
      )
    }) %>%
    ungroup() %>%
    group_by(modality, stratum, level, endpoint) %>%
    mutate(fdr_lm = p.adjust(p_lm, "BH"),
           fdr_wilcox = p.adjust(p_wilcox, "BH")) %>%
    ungroup()
}

#' Depth-adjusted sensitivity: cell-level random-intercept LMM (nlme; lme4 is
#' ABI-broken against Matrix in the spatialomics env)
#'   metric ~ group (+ run) + log10_counts + (1 | replicate)
#' nlme assigns between-replicate (containment) df to group, i.e.
#' n_replicates - (# between-replicate fixed effects), so cells are not
#' treated as independent replicates.
test_lmm_depth <- function(cm, elig, min_cells = 20,
                           metrics = c(mean_boot = "boot_prob", frac_low_conf = "low_conf",
                                       mean_corr = "avg_corr")) {
  keys <- elig %>% select(modality, stratum, level, cell_type)
  cm_ok <- cm %>%
    semi_join(keys, by = c("modality", "stratum", "level", "cell_type")) %>%
    group_by(modality, stratum, level, cell_type, replicate) %>%
    filter(n() >= min_cells) %>%
    ungroup()

  na_row <- function(ep) tibble(endpoint = ep, delta_adj = NA_real_, se_adj = NA_real_,
                                df_adj = NA_real_, p_adj_depth = NA_real_,
                                depth_slope = NA_real_)
  cm_ok %>%
    group_by(modality, stratum, level, cell_type) %>%
    group_modify(function(d, key) {
      n_run <- if (key$stratum == "seqFISH") n_distinct(d$run) else 1
      rhs <- if (n_run > 1) "group + run" else "group"
      d$log10_counts_c <- d$log10_counts - mean(d$log10_counts)
      purrr::map_dfr(names(metrics), function(ep) {
        if (var(d[[metrics[[ep]]]]) == 0) return(na_row(ep))
        f <- as.formula(glue("{metrics[[ep]]} ~ {rhs} + log10_counts_c"))
        fit <- tryCatch(suppressWarnings(nlme::lme(
          f, random = ~ 1 | replicate, data = d, method = "REML",
          control = nlme::lmeControl(opt = "optim"))), error = function(e) NULL)
        if (is.null(fit)) return(na_row(ep))
        tt <- summary(fit)$tTable
        tibble(endpoint = ep,
               delta_adj = tt["groupWildR", "Value"],
               se_adj = tt["groupWildR", "Std.Error"],
               df_adj = tt["groupWildR", "DF"],
               p_adj_depth = tt["groupWildR", "p-value"],
               depth_slope = tt["log10_counts_c", "Value"])
      })
    }) %>%
    ungroup() %>%
    group_by(modality, stratum, level, endpoint) %>%
    mutate(fdr_adj_depth = p.adjust(p_adj_depth, "BH")) %>%
    ungroup()
}

#' Global (all-cell) replicate summary for the overview figure.
summarise_global <- function(df) {
  df %>%
    group_by(modality, stratum, replicate, run, group) %>%
    summarise(
      n_cells = n(),
      median_log10_counts = median(log10_counts),
      across(all_of(glue("{mmc_levels}_boot_prob")), mean, .names = "mean_{.col}"),
      across(all_of(glue("{mmc_levels}_avg_corr")), mean, .names = "mean_{.col}"),
      .groups = "drop"
    )
}

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

#' Fig 1A: per-replicate ECDF of bootstrap probability at each level.
plot_ecdf_levels <- function(df_all) {
  long <- df_all %>%
    select(stratum, replicate, group, all_of(glue("{mmc_levels}_boot_prob"))) %>%
    pivot_longer(-c(stratum, replicate, group), names_to = "level", values_to = "boot_prob") %>%
    mutate(level = factor(sub("_boot_prob", "", level), levels = mmc_levels))
  ggplot(long, aes(boot_prob, colour = group, group = replicate)) +
    stat_ecdf(linewidth = 0.35, alpha = 0.8, pad = FALSE) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "#A5A5A5", linewidth = 0.3) +
    facet_grid(stratum ~ level) +
    scale_colour_manual(values = mmc_group_cols, name = NULL) +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    labs(x = "Bootstrap probability", y = "Cumulative fraction of cells",
         title = "Mapping confidence by taxonomy level (one line per replicate)") +
    theme_mmc()
}

#' Fig 1B: replicate mean subclass correlation vs sequencing / probe depth.
plot_depth_vs_quality <- function(glob_df) {
  ggplot(glob_df, aes(median_log10_counts, mean_subclass_avg_corr, colour = group)) +
    geom_point(size = 2.8, stroke = 0.6, shape = 21, fill = "white") +
    geom_point(size = 1.6) +
    facet_wrap(~stratum, scales = "free", nrow = 1) +
    scale_colour_manual(values = mmc_group_cols, name = NULL) +
    labs(x = "Median log10 counts per cell (replicate)",
         y = "Mean subclass correlation to reference",
         title = "Depth drives absolute fit; groups overlap") +
    theme_mmc()
}

#' Fig 2: landscape heatmap -- pooled mean bootstrap probability per cell type
#' per stratum (tested types only; blank = not testable in that stratum).
#' Rows = union of the top_n most abundant types in each stratum, grouped by class.
plot_quality_landscape <- function(rep_df, elig, top_n = 15) {
  d <- rep_df %>%
    semi_join(elig, by = c("modality", "stratum", "level", "cell_type")) %>%
    group_by(stratum, class_name, cell_type) %>%
    summarise(mean_boot = weighted.mean(mean_boot, n_cells),
              n_cells = sum(n_cells), .groups = "drop")
  keep <- d %>% group_by(stratum) %>% slice_max(n_cells, n = top_n) %>% pull(cell_type)
  d <- d %>% filter(cell_type %in% keep) %>%
    mutate(cell_type = factor(cell_type, levels = rev(sort(unique(cell_type)))),
           lab = sprintf("%.2f", mean_boot))
  ggplot(d, aes(stratum, cell_type, fill = mean_boot)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = lab, colour = mean_boot > 0.8), size = 3.6, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#0A0A0A")) +
    scale_fill_gradient(low = "#FBE7D9", high = "#9C3600", limits = c(NA, 1),
                        name = "Mean bootstrap\nprobability") +
    labs(x = NULL, y = NULL,
         title = "Mapping confidence by cell type and modality",
         subtitle = glue("Top {top_n} most abundant tested types per stratum; blank = not testable there")) +
    theme_mmc() + theme(legend.position = "right", panel.grid = element_blank())
}

#' Fig 3: summary heatmap of Cliff's delta (WildR vs SPF) per cell type x endpoint,
#' with FDR markers. Rows = cell types tested in >=1 stratum.
plot_effect_heatmap <- function(res, endpoints = c("mean_boot", "frac_low_conf", "mean_corr"),
                                fdr_col = "fdr_lm", max_types = 60) {
  d <- res %>% filter(endpoint %in% endpoints)
  top <- d %>% group_by(cell_type) %>%
    summarise(min_fdr = min(.data[[fdr_col]], na.rm = TRUE),
              n = sum(n_cells_WildR + n_cells_SPF)) %>%
    arrange(min_fdr, desc(n)) %>% slice_head(n = max_types)
  d <- d %>% filter(cell_type %in% top$cell_type) %>%
    left_join(mmc_endpoints %>% select(endpoint, label), by = "endpoint") %>%
    mutate(
      label = factor(label, levels = mmc_endpoints$label),
      sig = case_when(.data[[fdr_col]] < 0.01 ~ "**", .data[[fdr_col]] < 0.05 ~ "*", TRUE ~ ""),
      cell_type = factor(cell_type, levels = rev(sort(unique(cell_type))))
    )
  ggplot(d, aes(stratum, cell_type, fill = cliffs_delta)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sig), size = 4.2, colour = "#0A0A0A", vjust = 0.75) +
    facet_grid(~label, labeller = label_wrap_gen(22)) +
    scale_fill_gradient2(low = mmc_div_cols[["low"]], mid = mmc_div_cols[["mid"]],
                         high = mmc_div_cols[["high"]], midpoint = 0, limits = c(-1, 1),
                         name = "Cliff's δ\n(+ = WildR higher)") +
    labs(x = NULL, y = NULL,
         title = "WildR vs SPF mapping quality per cell type",
         subtitle = "Replicate-level test; * BH-FDR < 0.05, ** < 0.01 (none occur)") +
    theme_mmc(12) +
    theme(legend.position = "right", panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
}

#' Fig 3B: forest plot of the replicate-level mean difference (WildR - SPF) for one endpoint.
plot_effect_forest <- function(res, endpoint_sel = "mean_corr", fdr_col = "fdr_lm",
                               max_types = 40) {
  d <- res %>% filter(endpoint == endpoint_sel)
  top <- d %>% group_by(cell_type) %>%
    summarise(n = sum(n_cells_WildR + n_cells_SPF)) %>% slice_max(n, n = max_types)
  d <- d %>% filter(cell_type %in% top$cell_type) %>%
    mutate(sig = .data[[fdr_col]] < 0.05,
           cell_type = factor(cell_type, levels = rev(sort(unique(cell_type)))))
  lab <- mmc_endpoints$label[mmc_endpoints$endpoint == endpoint_sel]
  ggplot(d, aes(delta, cell_type)) +
    geom_vline(xintercept = 0, colour = "#A5A5A5", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, colour = "#767676",
                   linewidth = 0.4) +
    geom_point(aes(shape = sig), size = 2.4, colour = "#0A0A0A", fill = "#0A0A0A") +
    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 19),
                       labels = c(`FALSE` = "FDR ≥ 0.05", `TRUE` = "FDR < 0.05"), name = NULL) +
    facet_grid(~stratum, scales = "free_x") +
    labs(x = glue("Δ {lab} (WildR − SPF), 95% CI"), y = NULL,
         title = glue("Effect sizes: {lab}")) +
    theme_mmc()
}

#' Fig 4: single-distribution "superplot" for one cell type in one stratum.
#' (A) all cells at full scale (violin + subsampled jitter) with replicate means
#' overlaid; (B) replicate means only, zoomed, with group mean +/- 95% CI -- the
#' level at which the test is done. Stats (replicate LM, Wilcoxon, Cliff's delta,
#' depth-adjusted LMM) in the caption line.
plot_superplot <- function(cm, res, lmm, cell_type_sel, stratum_sel, level_sel,
                           metric = "avg_corr", endpoint_sel = "mean_corr",
                           max_cells_per_rep = 400, seed = 1, note = NULL,
                           min_cells = 20) {
  # only replicates that entered the test (>= min_cells of this type)
  d <- cm %>% filter(cell_type == cell_type_sel, stratum == stratum_sel, level == level_sel) %>%
    group_by(replicate) %>% filter(n() >= min_cells) %>% ungroup()
  rep_means <- d %>% group_by(group, replicate, run) %>%
    summarise(value = mean(.data[[metric]]), n = n(), .groups = "drop")
  set.seed(seed)
  d_pts <- d %>% group_by(replicate) %>% slice_sample(n = max_cells_per_rep) %>% ungroup()
  r <- res %>% filter(cell_type == cell_type_sel, stratum == stratum_sel,
                      level == level_sel, endpoint == endpoint_sel)
  a <- lmm %>% filter(cell_type == cell_type_sel, stratum == stratum_sel,
                      level == level_sel, endpoint == endpoint_sel)
  rep_unit <- if (stratum_sel == "seqFISH") "ROIs" else "animals"
  stat_txt <- glue(
    "Replicate-level LM ({r$model}; {r$n_rep_WildR} WildR vs {r$n_rep_SPF} SPF {rep_unit}): ",
    "\u0394 = {signif(r$delta, 2)} [95% CI {signif(r$ci_low, 2)}, {signif(r$ci_high, 2)}], ",
    "p = {signif(r$p_lm, 2)}, BH-FDR = {signif(r$fdr_lm, 2)}\n",
    "Wilcoxon p = {signif(r$p_wilcox, 2)}, Cliff's \u03b4 = {round(r$cliffs_delta, 2)} | ",
    "Depth-adjusted LMM: \u0394 = {signif(a$delta_adj, 2)}, p = {signif(a$p_adj_depth, 2)} (df = {a$df_adj})"
  )
  ylab <- c(avg_corr = "Correlation to reference", boot_prob = "Bootstrap probability")[[metric]]

  pa <- ggplot(d, aes(group, .data[[metric]], colour = group)) +
    geom_violin(aes(fill = group), alpha = 0.12, colour = NA, scale = "width", width = 0.8) +
    geom_jitter(data = d_pts, width = 0.18, height = 0, size = 0.25, alpha = 0.25) +
    geom_point(data = rep_means, aes(y = value), size = 3, shape = 21, fill = "white",
               stroke = 1, position = position_jitter(width = 0.06, seed = seed)) +
    scale_colour_manual(values = mmc_group_cols, guide = "none") +
    scale_fill_manual(values = mmc_group_cols, guide = "none") +
    labs(x = NULL, y = ylab, title = "All cells (open circles = replicate means)") +
    theme_mmc()

  pb <- ggplot(rep_means, aes(group, value, colour = group)) +
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.12,
                 colour = "#0A0A0A", linewidth = 0.5) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4, colour = "#0A0A0A",
                 linewidth = 0.35) +
    geom_point(size = 3, alpha = 0.9, position = position_jitter(width = 0.06, seed = seed)) +
    scale_colour_manual(values = mmc_group_cols, guide = "none") +
    labs(x = NULL, y = glue("Replicate mean {tolower(ylab)}"),
         title = "Replicate means \u00b1 95% CI (unit of test)") +
    theme_mmc()

  (pa | pb) + plot_layout(widths = c(1.3, 1)) +
    plot_annotation(title = glue("{cell_type_sel} \u2014 {stratum_sel} ({level_sel})"),
                    subtitle = note,
                    caption = stat_txt,
                    theme = theme(plot.title = element_text(face = "bold", size = 16),
                                  plot.subtitle = element_text(size = 12, colour = "#555555"),
                                  plot.caption = element_text(hjust = 0, size = 11,
                                                              colour = "#555555")))
}

#' Null check: p-value histograms of the primary replicate-level test, one
#' panel per taxonomy level x stratum (endpoints pooled; they are correlated,
#' so pooled counts are not independent). Under no group effect these are
#' ~uniform; a spike near 0 indicates signal.
plot_pvalue_hist <- function(res_all) {
  d <- res_all %>% filter(!is.na(p_lm)) %>%
    mutate(level = factor(level, levels = mmc_levels))
  exp_line <- d %>% count(level, stratum) %>% mutate(expected = n / 20)
  ggplot(d, aes(p_lm)) +
    geom_histogram(breaks = seq(0, 1, 0.05), fill = "#767676", colour = "white",
                   linewidth = 0.3) +
    geom_hline(data = exp_line, aes(yintercept = expected), linetype = "dashed",
               colour = "#E8590C", linewidth = 0.5) +
    facet_grid(level ~ stratum, scales = "free_y") +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    labs(x = "Replicate-level p-value (WildR vs SPF), all endpoints",
         y = "Tests", title = "p-value distributions vs the uniform null (dashed)") +
    theme_mmc()
}

#' Precision of the null: typical 95% CI half-width and the largest effect
#' still compatible with the data (90th percentile of max |CI bound|).
summarise_precision <- function(res_all) {
  res_all %>%
    filter(!is.na(ci_low)) %>%
    group_by(level, stratum, endpoint) %>%
    summarise(
      n_types = n(),
      n_nominal_p05 = sum(p_lm < 0.05, na.rm = TRUE),
      expected_nominal = round(0.05 * n(), 1),
      n_fdr05 = sum(fdr_lm < 0.05, na.rm = TRUE),
      n_robust = sum(robust_hit, na.rm = TRUE),
      median_abs_delta = median(abs(delta)),
      median_ci_halfwidth = median((ci_high - ci_low) / 2),
      p90_max_abs_ci = quantile(pmax(abs(ci_low), abs(ci_high)), 0.9),
      median_group_mean = median((mean_WildR + mean_SPF) / 2),
      r_delta_vs_depth_adj = cor(delta, delta_adj, use = "complete.obs"),
      .groups = "drop"
    )
}
