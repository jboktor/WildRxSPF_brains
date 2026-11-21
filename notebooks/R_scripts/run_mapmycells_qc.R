#!/usr/bin/env Rscript
# MapMyCells mapping-quality QC — WildR vs SPF, per cell type, SN + seqFISH.
#
# Run (spatialomics env; full run is ~20-40 min single-core -> SLURM):
#   sbatch notebooks/shell_scripts/run_mapmycells_qc.sh
# Quick check (SN AMY subclass only, writes to $QC_OUT_ROOT if set):
#   QC_QUICK=1 QC_OUT_ROOT=/tmp/... Rscript notebooks/R_scripts/run_mapmycells_qc.R
#
# Inputs : data/interim/MapMyCells_QC/mapmycells_scores_{SN,seqFISH}_<SCORES_DATE>.parquet
# Outputs: tables/MapMyCells_QC/*.csv, figures/MapMyCells_QC/*.png,
#          data/interim/MapMyCells_QC/mapmycells_qc_results_<date>.rds

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(glue); library(ggplot2)
  library(patchwork)
})

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
source(file.path(wkdir, "notebooks/R_scripts/mapmycells_qc_functions.R"))

log_msg <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), glue(...), "\n", sep = "")

# ---- parameters ------------------------------------------------------------
params <- list(
  scores_date      = Sys.getenv("SCORES_DATE", "2026-10-02"),
  levels_tested    = mmc_levels,   # class, subclass, supertype, cluster
  min_counts_SN    = 0,    # Parse output already requires >= 224 genes / nucleus
  min_counts_FISH  = 10,   # drops empty/near-empty segmentation masks only (~0.05%)
  min_cells        = 20,   # cells per replicate for that replicate to count
  min_rep_frac     = 0.5,  # >= half the replicates of each group meet min_cells
                           # (SN 3/6 animals, seqFISH 4/8 ROIs)
  run_tag          = Sys.getenv("RUN_TAG", "all-levels"),
  run_date         = Sys.Date()
)
quick <- Sys.getenv("QC_QUICK") == "1"
out_root <- Sys.getenv("QC_OUT_ROOT", wkdir)
fig_dir <- file.path(out_root, "figures/MapMyCells_QC")
tbl_dir <- file.path(out_root, "tables/MapMyCells_QC")
int_dir <- file.path(out_root, "data/interim/MapMyCells_QC")
for (d in c(fig_dir, tbl_dir, int_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
res_rds <- file.path(int_dir, glue("mapmycells_qc_results_{params$run_tag}_{params$run_date}.rds"))
if (file.exists(res_rds)) stop(glue("refusing to overwrite {res_rds}"))
log_msg("quick = {quick}; outputs under {out_root}")

# ---- load ------------------------------------------------------------------
score_dir <- file.path(wkdir, "data/interim/MapMyCells_QC")
sn <- load_mmc_scores(glue("{score_dir}/mapmycells_scores_SN_{params$scores_date}.parquet"),
                      "SN", params$min_counts_SN)
fish <- load_mmc_scores(glue("{score_dir}/mapmycells_scores_seqFISH_{params$scores_date}.parquet"),
                        "seqFISH", params$min_counts_FISH)
common <- intersect(names(sn), names(fish))
all_cells <- bind_rows(sn[common], fish[common]) %>%
  mutate(stratum = factor(stratum, levels = c("SN AMY", "SN HYP", "seqFISH")))
rm(sn, fish); gc()
if (quick) {
  all_cells <- all_cells %>% filter(stratum == "SN AMY")
  params$levels_tested <- "subclass"
}
log_msg("{nrow(all_cells)} cells; replicates: ",
        paste(capture.output(print(count(distinct(all_cells, stratum, replicate, group), stratum, group))),
              collapse = "\n"))

# ---- global summaries + overview figures ------------------------------------
glob <- summarise_global(all_cells)
write_csv(glob, file.path(tbl_dir, "MapMyCells_QC_global_replicate_summary.csv"))

mmc_save(plot_ecdf_levels(all_cells), fig_dir, "MapMyCells_QC_ecdf_levels.png", "ecdf")
mmc_save(plot_depth_vs_quality(glob), fig_dir, "MapMyCells_QC_depth_vs_quality.png", "depth")

# ---- per-level tests --------------------------------------------------------
results <- list(params = params, global = glob)
for (lvl in params$levels_tested) {
  log_msg("level = {lvl}: computing per-cell metrics")
  cm <- cell_metrics_at(all_cells, lvl)
  rep_df <- summarise_replicates(cm)
  elig <- eligible_types(rep_df, params$min_cells, params$min_rep_frac)
  log_msg("{lvl}: eligible cell types per stratum: ",
          paste(capture.output(print(count(elig, stratum))), collapse = "\n"))

  t0 <- Sys.time()
  rep_res <- test_replicate_level(rep_df, elig, params$min_cells)
  log_msg("{lvl}: replicate-level tests done ({nrow(rep_res)} rows)")
  lmm_res <- test_lmm_depth(cm, elig, params$min_cells)
  log_msg("{lvl}: LMM done in {round(difftime(Sys.time(), t0, units = 'mins'), 1)} min")

  res <- rep_res %>%
    left_join(lmm_res, by = c("modality", "stratum", "level", "cell_type", "endpoint")) %>%
    mutate(robust_hit = fdr_lm < 0.05 & p_adj_depth < 0.05 & sign(delta) == sign(delta_adj)) %>%
    arrange(stratum, endpoint, p_lm)

  write_csv(rep_df, file.path(tbl_dir, glue("MapMyCells_QC_replicate_summary_{lvl}.csv")))
  write_csv(elig, file.path(tbl_dir, glue("MapMyCells_QC_eligible_types_{lvl}.csv")))
  write_csv(res, file.path(tbl_dir, glue("MapMyCells_QC_tests_{lvl}.csv")))

  mmc_save(plot_quality_landscape(rep_df, elig, top_n = 15), fig_dir,
           glue("MapMyCells_QC_landscape_{lvl}.png"), "landscape")
  mmc_save(plot_effect_heatmap(res, max_types = 50), fig_dir,
           glue("MapMyCells_QC_effect-heatmap_{lvl}.png"),
           if (lvl == "class") "heatmap_class" else "heatmap")
  for (ep in c("mean_corr", "mean_boot"))
    mmc_save(plot_effect_forest(res, ep, max_types = 35), fig_dir,
             glue("MapMyCells_QC_forest_{lvl}_{ep}.png"), "forest")

  # Single-distribution example: strongest robust hit, else the smallest-p test
  # among fully replicated cell types (all replicates of both groups present).
  hits <- res %>% filter(robust_hit, endpoint %in% c("mean_corr", "mean_boot"))
  n_full <- glob %>% count(stratum, group) %>%
    pivot_wider(names_from = group, values_from = n, names_prefix = "full_")
  pick <- if (nrow(hits)) slice_min(hits, fdr_lm, n = 1, with_ties = FALSE) else
    res %>% filter(endpoint %in% c("mean_corr", "mean_boot")) %>%
      left_join(n_full, by = "stratum") %>%
      filter(n_rep_WildR == full_WildR, n_rep_SPF == full_SPF) %>%
      slice_min(p_lm, n = 1, with_ties = FALSE)
  metric <- c(mean_corr = "avg_corr", mean_boot = "boot_prob")[[pick$endpoint]]
  p_sp <- plot_superplot(cm, res, lmm_res, pick$cell_type, as.character(pick$stratum), lvl,
                         metric = metric, endpoint_sel = pick$endpoint,
                         min_cells = params$min_cells,
                         note = if (nrow(hits)) "Strongest robust hit (FDR < 0.05, depth-adjusted p < 0.05)" else
                           "No robust hit at this level: strongest nominal test shown as an example")
  mmc_save(p_sp, fig_dir, glue("MapMyCells_QC_superplot_{lvl}_top.png"), "superplot")

  results[[lvl]] <- list(tests = res, replicate = rep_df, eligible = elig,
                         superplot_pick = pick %>% mutate(is_robust_hit = nrow(hits) > 0))
  rm(cm); gc()
}

# ---- cross-level null checks -------------------------------------------------
res_all <- bind_rows(lapply(params$levels_tested, function(l) results[[l]]$tests))
prec <- summarise_precision(res_all)
write_csv(prec, file.path(tbl_dir, "MapMyCells_QC_precision_summary.csv"))
mmc_save(plot_pvalue_hist(res_all), fig_dir, "MapMyCells_QC_pvalue_hist.png", "pvalue")
results$precision <- prec

saveRDS(results, res_rds)
log_msg("wrote {res_rds}")
