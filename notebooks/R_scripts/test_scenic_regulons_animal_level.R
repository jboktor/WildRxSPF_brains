#!/usr/bin/env Rscript
# Animal-level (pseudobulk) WildR vs SPF test of pySCENIC regulon activity,
# aggregated across all tissue × class, subclass or supertype splits.
#
# Why: the original per-split tests were nucleus-level Wilcoxon tests, which
# treat nuclei as independent replicates (p < 1e-40 for ΔAUC ~ 1e-3). The unit
# of replication is the animal (6 WildR vs 6 SPF per tissue). Since 2026-10-06
# run_scenic_one_split.R runs the same animal-level test per split; this script
# recomputes it from saved AUCell matrices so every split uses one code path
# (scenic_animal_level_helpers.R) and adds cross-split corrections.
#
# Method (no GRN re-fit):
#   1. per split, join auc_mtx.csv to the split metadata (sample -> animal)
#   2. per animal × regulon: mean AUCell (animals with >= SCENIC_MIN_CELLS nuclei)
#   3. splits need >= SCENIC_MIN_ANIMALS animals per group, else underpowered
#   4. pooled-variance t on animal means (== lm(mean ~ group)), Hedges g,
#      Cliff's delta, exact Wilcoxon, complete-separation count
#   5. BH-FDR within split (primary) and across all tests (global)
#   6. class level only: join the original nucleus-level result for comparison
#   7. cross-split TF consensus (per-animal mean of within-split z), per tissue
#
# Outputs (date-stamped; <sfx> = "" for class, "_subclass" / "_supertype" otherwise):
#   data/interim/scenic/animal_level/regulon_auc_per_animal<sfx>_<date>.csv
#   data/interim/scenic/animal_level/regulon_tf_consensus_per_animal<sfx>_<date>.csv
#   tables/SCENIC/SCENIC_animal_level_tests<sfx>_<date>.csv
#   tables/SCENIC/SCENIC_animal_level_split_summary<sfx>_<date>.csv
#   tables/SCENIC/SCENIC_tf_consensus_animal_level<sfx>_<date>.csv
#
# Run:
#   /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript \
#     notebooks/R_scripts/test_scenic_regulons_animal_level.R [--level class|subclass|supertype]

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(glue)
})

args <- commandArgs(trailingOnly = TRUE)
level <- if ("--level" %in% args) args[match("--level", args) + 1] else "class"
stopifnot(level %in% c("class", "subclass", "supertype"))

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
by_class <- file.path(wkdir, "data/interim/scenic", paste0("by_", level))
meta_dir <- file.path(wkdir, if (level == "class") "data/input/scenic/snRNASeq" else glue("data/interim/scenic/looms_{level}"))
out_interim <- file.path(wkdir, "data/interim/scenic/animal_level")
out_tbl <- file.path(wkdir, "tables/SCENIC")
dir.create(out_interim, showWarnings = FALSE, recursive = TRUE)
dir.create(out_tbl, showWarnings = FALSE, recursive = TRUE)

source(file.path(wkdir, "notebooks/R_scripts/scenic_animal_level_helpers.R"))
MIN_CELLS <- SCENIC_MIN_CELLS
MIN_ANIMALS <- SCENIC_MIN_ANIMALS
run_date <- Sys.Date()
sfx <- if (level == "class") "" else glue("_{level}")
`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

log_msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", paste0(...), "\n")

summaries <- list.files(by_class, pattern = "^run_summary\\.json$", recursive = TRUE, full.names = TRUE)
log_msg("Splits with run_summary.json: ", length(summaries))

per_animal_all <- list()
tests_all <- list()
split_info <- list()

for (sj in summaries) {
  s <- jsonlite::fromJSON(sj)
  d <- dirname(sj)
  auc_f <- file.path(d, "auc_mtx.csv")
  meta_f <- file.path(meta_dir, glue("{s$stem}_meta.csv"))
  if (!file.exists(auc_f) || !file.exists(meta_f)) {
    log_msg("SKIP ", s$stem, " (missing auc or meta)")
    next
  }
  auc <- fread(auc_f)
  setnames(auc, 1, "cell_id")
  meta <- fread(meta_f, select = c("cell_id", "sample", "microbiome"))
  auc <- merge(meta, auc, by = "cell_id")
  n_unmatched <- s$n_cells - nrow(auc)
  auc[, animal := animal_from_sample(sample)]

  reg_cols <- setdiff(names(auc), c("cell_id", "sample", "microbiome", "animal"))
  # Per-animal means (all animals kept here; filter for testing below)
  pa <- auc[, c(list(n_cells = .N), lapply(.SD, mean)), by = .(animal, microbiome), .SDcols = reg_cols]
  pa_long <- melt(pa, id.vars = c("animal", "microbiome", "n_cells"),
                  variable.name = "regulon", value.name = "mean_auc")
  pa_long[, `:=`(stem = s$stem, tissue = s$tissue, class_label = s$class_label,
                 class_name = s$class_name, subclass_name = s$subclass_name %||% NA_character_,
                 supertype_name = s$supertype_name %||% NA_character_)]
  per_animal_all[[s$stem]] <- pa_long

  keep <- pa[n_cells >= MIN_CELLS]
  n_w <- sum(keep$microbiome == "WildR"); n_s <- sum(keep$microbiome == "SPF")
  tested <- n_w >= MIN_ANIMALS && n_s >= MIN_ANIMALS
  split_info[[s$stem]] <- tibble(
    stem = s$stem, tissue = s$tissue, class_label = s$class_label, class_name = s$class_name,
    subclass_name = s$subclass_name %||% NA_character_,
    supertype_name = s$supertype_name %||% NA_character_,
    n_cells = nrow(auc), n_cells_unmatched = n_unmatched, n_regulons = length(reg_cols),
    n_animals_WildR_total = sum(pa$microbiome == "WildR"),
    n_animals_SPF_total = sum(pa$microbiome == "SPF"),
    n_animals_WildR_ok = n_w, n_animals_SPF_ok = n_s,
    min_cells_per_animal = min(pa$n_cells), median_cells_per_animal = median(pa$n_cells),
    tested = tested
  )
  if (!tested) {
    log_msg("UNDERPOWERED ", s$stem, " (animals >= ", MIN_CELLS, " cells: WildR ", n_w, ", SPF ", n_s, ")")
    next
  }
  m <- t(as.matrix(keep[, ..reg_cols]))
  colnames(m) <- keep$animal
  res <- row_tests(m, keep$microbiome) %>%
    mutate(stem = s$stem, tissue = s$tissue, class_label = s$class_label,
           class_name = s$class_name, subclass_name = s$subclass_name %||% NA_character_,
           supertype_name = s$supertype_name %||% NA_character_,
           n_cells = nrow(auc),
           fdr_split = p.adjust(p_animal, "BH"))
  tests_all[[s$stem]] <- res
  log_msg(s$stem, ": ", length(reg_cols), " regulons, animals ", n_w, "v", n_s,
          ", nominal p<0.05: ", sum(res$p_animal < 0.05, na.rm = TRUE),
          ", FDR<0.10: ", sum(res$fdr_split < 0.10, na.rm = TRUE))
}

per_animal <- rbindlist(per_animal_all)
fwrite(per_animal, file.path(out_interim, glue("regulon_auc_per_animal{sfx}_{run_date}.csv")))

# Nucleus-level results (original per-split Wilcoxon; class level only) for comparison
cell_level <- if (level == "class") {
  read_csv(file.path(wkdir, "data/processed/scenic_by_class/diff_regulons_all_splits.csv"),
           show_col_types = FALSE) %>%
    transmute(stem = paste0("snRNASeq_", tissue, "_", class_label),
              regulon = str_replace(regulon, "\\.pos$", "(+)") %>% str_replace("\\.neg$", "(-)"),
              delta_cell = delta, p_cell = p, padj_cell = padj)
} else {
  tibble(stem = character(), regulon = character(), delta_cell = numeric(),
         p_cell = numeric(), padj_cell = numeric())
}

tests <- bind_rows(tests_all) %>%
  mutate(fdr_global = p.adjust(p_animal, "BH"),
         tf = tf_from_regulon(regulon),
         direction = if_else(delta_animal > 0, "WildR higher", "SPF higher")) %>%
  left_join(cell_level, by = c("stem", "regulon")) %>%
  mutate(cell_level_hit = !is.na(padj_cell) & padj_cell < 0.05) %>%
  relocate(tissue, class_name, regulon, tf, direction, delta_animal, ci_low, ci_high,
           hedges_g, cliffs_delta, animals_separated, p_animal, fdr_split, fdr_global,
           p_wilcox_animal, cell_level_hit, delta_cell, padj_cell) %>%
  arrange(p_animal)

n_join <- sum(!is.na(tests$p_cell))
log_msg("Cell-level join matched ", n_join, " / ", nrow(tests), " animal-level tests")

write_csv(tests, file.path(out_tbl, glue("SCENIC_animal_level_tests{sfx}_{run_date}.csv")))

split_summary <- bind_rows(split_info) %>%
  left_join(
    tests %>% group_by(stem) %>% summarise(
      n_tested = n(),
      n_nominal_p05 = sum(p_animal < 0.05, na.rm = TRUE),
      n_expected_p05 = round(0.05 * n()),
      n_fdr10 = sum(fdr_split < 0.10, na.rm = TRUE),
      n_cell_level_hits = sum(cell_level_hit),
      n_cell_hits_nominal_animal = sum(cell_level_hit & p_animal < 0.05, na.rm = TRUE),
      n_cell_hits_concordant_sign = sum(cell_level_hit & sign(delta_cell) == sign(delta_animal), na.rm = TRUE),
      .groups = "drop"),
    by = "stem") %>%
  arrange(tissue, desc(n_cells))
write_csv(split_summary, file.path(out_tbl, glue("SCENIC_animal_level_split_summary{sfx}_{run_date}.csv")))

# ---- Cross-class consensus per TF (animal-level) ----
# Regulons are inferred separately per split, so raw AUC is not comparable
# across splits. Instead: within each split, z-score the 12 animal means of a
# TF's regulon; average each animal's z across the splits where that regulon
# was called (>= MIN_CELLS nuclei); 6 v 6 test on the per-animal consensus.
# Splits share animals, so this is one test per TF × tissue (not per split).
tested_stems <- unique(tests$stem)
consensus_input <- per_animal[stem %in% tested_stems & n_cells >= MIN_CELLS]
consensus_input[, tf := tf_from_regulon(as.character(regulon))]
consensus_input <- consensus_input[grepl("\\(\\+\\)$", regulon)]
consensus_input[, z := {
  s <- sd(mean_auc); if (is.na(s) || s == 0) NA_real_ else (mean_auc - mean(mean_auc)) / s
}, by = .(stem, regulon)]
consensus_animal <- consensus_input[!is.na(z), .(z_consensus = mean(z), n_splits = uniqueN(stem)),
                                    by = .(tissue, tf, animal, microbiome)]
fwrite(consensus_animal, file.path(out_interim, glue("regulon_tf_consensus_per_animal{sfx}_{run_date}.csv")))

consensus <- consensus_animal %>%
  group_by(tissue, tf) %>%
  filter(max(n_splits) >= 3, sum(microbiome == "WildR") >= MIN_ANIMALS,
         sum(microbiome == "SPF") >= MIN_ANIMALS) %>%
  group_modify(function(d, key) {
    m <- matrix(d$z_consensus, nrow = 1, dimnames = list(key$tf, d$animal))
    row_tests(m, d$microbiome) %>%
      select(n_animals_WildR, n_animals_SPF, delta_z = delta_animal, ci_low, ci_high,
             hedges_g, cliffs_delta, animals_separated, p_animal, p_wilcox_animal) %>%
      mutate(n_splits = max(d$n_splits))
  }) %>%
  ungroup() %>%
  group_by(tissue) %>%
  mutate(fdr_tissue = p.adjust(p_animal, "BH")) %>%
  ungroup() %>%
  left_join(
    tests %>% filter(cell_level_hit) %>% group_by(tissue, tf) %>%
      summarise(n_cell_hit_splits = n(),
                n_cell_hit_WildR_up = sum(delta_cell > 0),
                n_cell_hit_SPF_up = sum(delta_cell < 0), .groups = "drop"),
    by = c("tissue", "tf")) %>%
  mutate(across(starts_with("n_cell_hit"), ~ coalesce(.x, 0L))) %>%
  arrange(p_animal)
write_csv(consensus, file.path(out_tbl, glue("SCENIC_tf_consensus_animal_level{sfx}_{run_date}.csv")))
log_msg("Consensus TF tests: ", nrow(consensus), "; nominal p<0.05: ",
        sum(consensus$p_animal < 0.05), "; FDR_tissue<0.10: ", sum(consensus$fdr_tissue < 0.10))

log_msg("Tests: ", nrow(tests)," across ", n_distinct(tests$stem), " splits; nominal p<0.05: ",
        sum(tests$p_animal < 0.05, na.rm = TRUE), " (expected under null ≈ ",
        round(0.05 * nrow(tests)), "); FDR_split<0.10: ", sum(tests$fdr_split < 0.10, na.rm = TRUE),
        "; FDR_global<0.10: ", sum(tests$fdr_global < 0.10, na.rm = TRUE))
log_msg("Wrote ", out_tbl)
