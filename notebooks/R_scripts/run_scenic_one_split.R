#!/usr/bin/env Rscript
# Run pySCENIC grn → ctx → aucell for one loom (tissue, class or subclass split),
# then test WildR vs SPF regulon activity at the ANIMAL level (unit of
# replication; see scenic_animal_level_helpers.R). Steps whose outputs exist
# are skipped, so re-running on a finished split only redoes the test.
#
# Split level is inferred from the loom stem:
#   snRNASeq_{AMY|HYP}                  -> tissue   (data/interim/scenic/by_tissue)
#   snRNASeq_{tissue}_<class_label>     -> class    (data/interim/scenic/by_class)
#   snRNASeq_{tissue}_SUBC_<id>         -> subclass (data/interim/scenic/by_subclass)
#   snRNASeq_{tissue}_SUPT_<id>         -> supertype (data/interim/scenic/by_supertype)
# Metadata is read from <stem>_meta.csv next to the loom.
#
# Outputs per split (out_dir): adjacencies.tsv, regulons.csv, auc_mtx.csv,
#   regulon_auc_per_animal.csv, diff_regulons_WildR_vs_SPF_animal.csv,
#   run_summary.json. (Older runs also hold the superseded nucleus-level
#   diff_regulons_WildR_vs_SPF.csv; it is left untouched, and an older
#   run_summary.json is copied to run_summary_celllevel.json before update.)
#
# Usage:
#   Rscript run_scenic_one_split.R --loom snRNASeq_AMY_CS20230722_CLAS_01.loom
#   Rscript run_scenic_one_split.R --loom /full/path/to/snRNASeq_AMY_SUBC_064.loom

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(tictoc)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1]]
}

loom_arg <- get_arg("--loom")
if (is.null(loom_arg)) stop("Provide --loom <file.loom>")

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
input_dir <- file.path(wkdir, "data/input/scenic/snRNASeq")
db_dir <- file.path(wkdir, "data/input/scenic/cisTarget")

py_bin <- "/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python"
pyscenic_bin <- "/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/pyscenic"
n_workers <- as.integer(Sys.getenv("SCENIC_N_WORKERS", "16"))

loom_path <- if (grepl("^/", loom_arg)) loom_arg else file.path(input_dir, loom_arg)
stopifnot(file.exists(loom_path))
stem <- sub("\\.loom$", "", basename(loom_path))

source(file.path(wkdir, "notebooks/R_scripts/scenic_animal_level_helpers.R"))

split_level <- dplyr::case_when(
  grepl("^snRNASeq_(AMY|HYP)$", stem) ~ "tissue",
  grepl("_SUBC_", stem) ~ "subclass",
  grepl("_SUPT_", stem) ~ "supertype",
  TRUE ~ "class"
)
is_tissue <- split_level == "tissue"
out_root <- file.path(wkdir, "data/interim/scenic", paste0("by_", split_level))
fig_root <- file.path(wkdir, "figures/SCENIC", paste0("snRNASeq_by_", split_level))

out_dir <- file.path(out_root, stem)
fig_dir <- file.path(fig_root, stem)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

adj_tsv <- file.path(out_dir, "adjacencies.tsv")
regulons_csv <- file.path(out_dir, "regulons.csv")
auc_csv <- file.path(out_dir, "auc_mtx.csv")
diff_csv <- file.path(out_dir, "diff_regulons_WildR_vs_SPF_animal.csv")
per_animal_csv <- file.path(out_dir, "regulon_auc_per_animal.csv")
summary_json <- file.path(out_dir, "run_summary.json")

tf_list <- file.path(db_dir, "allTFs_mm.txt")
motif_tbl <- file.path(db_dir, "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl")
rank_dbs <- list.files(db_dir, pattern = "genes_vs_motifs\\.rankings\\.feather$", full.names = TRUE)
stopifnot(file.exists(tf_list), file.exists(motif_tbl), length(rank_dbs) >= 1)

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", paste0(...), "\n")
  flush.console()
}

run_cli <- function(cmd_args, label) {
  log_msg("START ", label)
  tic(label)
  status <- system2(pyscenic_bin, args = cmd_args)
  toc()
  if (status != 0) stop(label, " failed with status ", status)
  log_msg("DONE ", label)
}

log_msg("Split: ", stem, " level=", split_level, " workers=", n_workers)
log_msg("Loom: ", loom_path)

if (!file.exists(adj_tsv) || file.info(adj_tsv)$size < 1000) {
  run_cli(
    c("grn", "--num_workers", as.character(n_workers), "-o", adj_tsv, loom_path, tf_list),
    "pyscenic grn"
  )
} else {
  log_msg("Skip grn")
}

if (!file.exists(regulons_csv) || file.info(regulons_csv)$size < 1000) {
  run_cli(
    c(
      "ctx", adj_tsv, rank_dbs,
      "--annotations_fname", motif_tbl,
      "--expression_mtx_fname", loom_path,
      "--mode", "custom_multiprocessing",
      "--output", regulons_csv,
      "--num_workers", as.character(n_workers),
      "--mask_dropouts"
    ),
    "pyscenic ctx"
  )
} else {
  log_msg("Skip ctx")
}

if (!file.exists(auc_csv) || file.info(auc_csv)$size < 1000) {
  run_cli(
    c("aucell", loom_path, regulons_csv, "-o", auc_csv, "--num_workers", as.character(n_workers)),
    "pyscenic aucell"
  )
} else {
  log_msg("Skip aucell")
}

# ---- Downstream: WildR vs SPF at the animal level ----
meta_path <- file.path(dirname(loom_path), paste0(stem, "_meta.csv"))
stopifnot(file.exists(meta_path))
meta <- read_csv(meta_path, show_col_types = FALSE)
first1 <- function(col) if (col %in% names(meta)) unique(as.character(meta[[col]]))[1] else NA_character_
split_info <- list(
  tissue = first1("tissue"),
  class_label = if (!is_tissue) first1("class_label") else NA_character_,
  class_name = if (!is_tissue) first1("class_name") else NA_character_,
  subclass_name = if (split_level %in% c("subclass", "supertype")) first1("subclass_name") else NA_character_,
  supertype_name = if (split_level == "supertype") first1("supertype_name") else NA_character_
)

auc <- fread(auc_csv)
am <- animal_means(auc, meta)
log_msg("AUCell: ", nrow(auc), " nuclei (", am$n_cells_matched, " matched to meta), ",
        length(am$regulons), " regulons, ", nrow(am$wide), " animals")
fwrite(am$per_animal[, `:=`(stem = stem, split_level = split_level)], per_animal_csv)

tst <- animal_level_test(am)
if (is.null(tst$res)) {
  log_msg("UNDERPOWERED for animal-level test (", tst$reason, "); writing empty diff table")
  diff_reg <- tibble(regulon = character(), p_animal = numeric(), fdr_split = numeric())
} else {
  diff_reg <- tst$res %>%
    mutate(stem = stem, split_level = split_level, !!!split_info,
           n_cells = am$n_cells_matched) %>%
    arrange(p_animal)
}
write_csv(diff_reg, diff_csv)

if (nrow(diff_reg)) {
  group_cols <- c(WildR = "#E8590C", SPF = "#2E6FA8")
  vd <- diff_reg %>%
    mutate(flag = case_when(fdr_split < 0.10 ~ "FDR < 0.10", p_animal < 0.05 ~ "p < 0.05", TRUE ~ "n.s."),
           lab = if_else(rank(p_animal, ties.method = "first") <= 8, regulon, NA_character_))
  p_volc <- ggplot(vd, aes(hedges_g, -log10(p_animal))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#A5A5A5") +
    geom_point(aes(colour = flag), alpha = 0.8, size = 1.8) +
    ggrepel::geom_text_repel(aes(label = lab), size = 3.4, fontface = "italic", seed = 1,
                             min.segment.length = 0, max.overlaps = Inf) +
    scale_colour_manual(values = c("FDR < 0.10" = "#E8590C", "p < 0.05" = "#27455C", "n.s." = "#D4D4D4"),
                        name = NULL) +
    labs(title = glue("{stem} ({split_level})"),
         subtitle = glue("Animal-level WildR − SPF ({tst$n_WildR} vs {tst$n_SPF} animals, {comma_n} nuclei)",
                         comma_n = format(am$n_cells_matched, big.mark = ",")),
         x = "Hedges g on animal means (WildR − SPF)", y = "−log10 animal-level p") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(), legend.position = "top",
          plot.title = element_text(face = "bold"), plot.title.position = "plot",
          plot.background = element_rect(fill = "white", colour = NA))
  ggsave(file.path(fig_dir, "diff_regulon_volcano_animal.png"), p_volc,
         width = 9, height = 6.5, dpi = 200, bg = "white")
}

# Keep the earlier (nucleus-level) summary before replacing it
if (file.exists(summary_json)) {
  old <- jsonlite::fromJSON(summary_json)
  if (is.null(old$test)) file.copy(summary_json, file.path(out_dir, "run_summary_celllevel.json"), overwrite = FALSE)
}

summary <- c(
  list(stem = stem, split_level = split_level),
  split_info,
  list(
    n_cells = am$n_cells_matched,
    n_regulons = length(am$regulons),
    test = "animal_level_pooled_t",
    min_cells_per_animal = SCENIC_MIN_CELLS,
    n_animals_WildR = tst$n_WildR, n_animals_SPF = tst$n_SPF,
    tested = !is.null(tst$res),
    n_nominal_p05 = sum(diff_reg$p_animal < 0.05, na.rm = TRUE),
    n_fdr10 = sum(diff_reg$fdr_split < 0.10, na.rm = TRUE),
    auc_csv = auc_csv, regulons_csv = regulons_csv, diff_csv = diff_csv,
    per_animal_csv = per_animal_csv,
    input_matrix = "raw_counts_UMI",
    date = as.character(Sys.Date())
  )
)
jsonlite::write_json(summary, summary_json, auto_unbox = TRUE, pretty = TRUE, na = "null")
log_msg("DONE split ", stem, " regulons=", length(am$regulons),
        " animals=", tst$n_WildR, "v", tst$n_SPF,
        " p<0.05=", summary$n_nominal_p05, " FDR<0.10=", summary$n_fdr10)
