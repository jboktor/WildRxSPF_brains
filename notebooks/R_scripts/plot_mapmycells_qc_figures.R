#!/usr/bin/env Rscript
# Redraw every MapMyCells QC figure from saved results -- no statistics are
# recomputed. Use after changing plot code or figure sizes.
#
# Run (spatialomics env, ~2 min, ~6 GB RAM for the per-cell parquet):
#   /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript \
#     notebooks/R_scripts/plot_mapmycells_qc_figures.R
#
# Inputs : data/interim/MapMyCells_QC/mapmycells_qc_results_<RESULTS_TAG>.rds
#          data/interim/MapMyCells_QC/mapmycells_scores_{SN,seqFISH}_<SCORES_DATE>.parquet
# Outputs: figures/MapMyCells_QC/*.png (same file names the driver writes)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(glue); library(ggplot2); library(patchwork)
})

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
source(file.path(wkdir, "notebooks/R_scripts/mapmycells_qc_functions.R"))
log_msg <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), glue(...), "\n", sep = "")

results_tag <- Sys.getenv("RESULTS_TAG", "all-levels_2026-10-02")
scores_date <- Sys.getenv("SCORES_DATE", "2026-10-02")
int_dir <- file.path(wkdir, "data/interim/MapMyCells_QC")
fig_dir <- file.path(wkdir, "figures/MapMyCells_QC")

res <- readRDS(file.path(int_dir, glue("mapmycells_qc_results_{results_tag}.rds")))
p <- res$params
levels_done <- intersect(mmc_levels, names(res))

sn <- load_mmc_scores(file.path(int_dir, glue("mapmycells_scores_SN_{scores_date}.parquet")),
                      "SN", p$min_counts_SN)
fish <- load_mmc_scores(file.path(int_dir, glue("mapmycells_scores_seqFISH_{scores_date}.parquet")),
                        "seqFISH", p$min_counts_FISH)
common <- intersect(names(sn), names(fish))
all_cells <- bind_rows(sn[common], fish[common]) %>%
  mutate(stratum = factor(stratum, levels = c("SN AMY", "SN HYP", "seqFISH")))
rm(sn, fish); invisible(gc())

mmc_save(plot_ecdf_levels(all_cells), fig_dir, "MapMyCells_QC_ecdf_levels.png", "ecdf")
mmc_save(plot_depth_vs_quality(res$global), fig_dir, "MapMyCells_QC_depth_vs_quality.png", "depth")

for (lvl in levels_done) {
  r <- res[[lvl]]
  mmc_save(plot_quality_landscape(r$replicate, r$eligible, top_n = 15), fig_dir,
           glue("MapMyCells_QC_landscape_{lvl}.png"), "landscape")
  mmc_save(plot_effect_heatmap(r$tests, max_types = 50), fig_dir,
           glue("MapMyCells_QC_effect-heatmap_{lvl}.png"),
           if (lvl == "class") "heatmap_class" else "heatmap")
  for (ep in c("mean_corr", "mean_boot"))
    mmc_save(plot_effect_forest(r$tests, ep, max_types = 35), fig_dir,
             glue("MapMyCells_QC_forest_{lvl}_{ep}.png"), "forest")

  pick <- r$superplot_pick
  cm <- cell_metrics_at(all_cells, lvl)
  lmm <- r$tests %>% select(modality, stratum, level, cell_type, endpoint,
                            delta_adj, p_adj_depth, df_adj)
  sp <- plot_superplot(
    cm, r$tests, lmm, pick$cell_type, as.character(pick$stratum), lvl,
    metric = c(mean_corr = "avg_corr", mean_boot = "boot_prob")[[pick$endpoint]],
    endpoint_sel = pick$endpoint, min_cells = p$min_cells,
    note = if (isTRUE(pick$is_robust_hit)) "Strongest robust hit (FDR < 0.05, depth-adjusted p < 0.05)" else
      "No robust hit at this level: strongest nominal test shown as an example")
  mmc_save(sp, fig_dir, glue("MapMyCells_QC_superplot_{lvl}_top.png"), "superplot")
  rm(cm); invisible(gc())
  log_msg("{lvl}: figures redrawn")
}

res_all <- bind_rows(lapply(levels_done, function(l) res[[l]]$tests))
mmc_save(plot_pvalue_hist(res_all), fig_dir, "MapMyCells_QC_pvalue_hist.png", "pvalue")
log_msg("all figures redrawn in {fig_dir}")
