#!/usr/bin/env Rscript
# Run pySCENIC grn → ctx → aucell for one loom (class or tissue split).
#
# Usage:
#   Rscript run_scenic_one_split.R --loom snRNASeq_AMY_CS20230722_CLAS_01.loom
#   Rscript run_scenic_one_split.R --loom snRNASeq_AMY.loom
#   Rscript run_scenic_one_split.R --loom /full/path/to/file.loom

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(tictoc)
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

# snRNASeq_AMY.loom / snRNASeq_HYP.loom → tissue; else class
is_tissue <- grepl("^snRNASeq_(AMY|HYP)$", stem)
split_level <- if (is_tissue) "tissue" else "class"
out_root <- file.path(wkdir, "data/interim/scenic", if (is_tissue) "by_tissue" else "by_class")
fig_root <- file.path(
  wkdir, "figures/SCENIC",
  if (is_tissue) "snRNASeq_by_tissue" else "snRNASeq_by_class"
)

out_dir <- file.path(out_root, stem)
fig_dir <- file.path(fig_root, stem)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

adj_tsv <- file.path(out_dir, "adjacencies.tsv")
regulons_csv <- file.path(out_dir, "regulons.csv")
auc_csv <- file.path(out_dir, "auc_mtx.csv")
diff_csv <- file.path(out_dir, "diff_regulons_WildR_vs_SPF.csv")
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

# ---- Downstream: WildR vs SPF ----
suppressPackageStartupMessages({
  library(ggrepel)
})

meta_path <- file.path(input_dir, paste0(stem, "_meta.csv"))
meta <- read_csv(meta_path, show_col_types = FALSE)
auc <- read.csv(auc_csv, row.names = 1, check.names = FALSE)
auc <- as.matrix(auc)
if (mean(rownames(auc) %in% meta$cell_id) > 0.5) {
  cell_ids <- rownames(auc)
  auc_mat <- t(auc)
} else {
  cell_ids <- colnames(auc)
  auc_mat <- auc
}
rownames(auc_mat) <- make.unique(gsub("\\(\\+\\)", ".pos", gsub("\\(-\\)", ".neg", rownames(auc_mat))))
rownames(auc_mat) <- make.unique(gsub("[^A-Za-z0-9.-]", ".", rownames(auc_mat)))

meta <- meta %>% filter(cell_id %in% cell_ids)
common <- intersect(cell_ids, meta$cell_id)
auc_mat <- auc_mat[, common, drop = FALSE]
meta <- meta %>% filter(cell_id %in% common) %>% arrange(match(cell_id, common))

# Primary test: overall WildR vs SPF within this split
diff_reg <- as.data.frame(t(auc_mat)) %>%
  mutate(
    cell_id = common,
    microbiome = meta$microbiome[match(cell_id, meta$cell_id)]
  ) %>%
  pivot_longer(cols = all_of(rownames(auc_mat)), names_to = "regulon", values_to = "auc") %>%
  group_by(regulon) %>%
  summarise(
    n_SPF = sum(microbiome == "SPF"),
    n_WildR = sum(microbiome == "WildR"),
    mean_SPF = mean(auc[microbiome == "SPF"]),
    mean_WildR = mean(auc[microbiome == "WildR"]),
    delta = mean_WildR - mean_SPF,
    p = tryCatch(wilcox.test(auc ~ microbiome)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    padj = p.adjust(p, method = "BH"),
    tissue = unique(meta$tissue)[1],
    split_level = split_level,
    class_label = if ("class_label" %in% names(meta) && !is_tissue) {
      unique(meta$class_label)[1]
    } else {
      NA_character_
    },
    class_name = if ("class_name" %in% names(meta) && !is_tissue) {
      unique(as.character(meta$class_name))[1]
    } else {
      NA_character_
    }
  ) %>%
  arrange(padj, p)

write_csv(diff_reg, diff_csv)

# Tissue runs: also test within class_label (exploratory; same tissue-wide regulons)
if (is_tissue && "class_label" %in% names(meta)) {
  diff_by_class <- as.data.frame(t(auc_mat)) %>%
    mutate(
      cell_id = common,
      microbiome = meta$microbiome[match(cell_id, meta$cell_id)],
      class_label = meta$class_label[match(cell_id, meta$cell_id)],
      class_name = meta$class_name[match(cell_id, meta$cell_id)]
    ) %>%
    pivot_longer(cols = all_of(rownames(auc_mat)), names_to = "regulon", values_to = "auc") %>%
    group_by(class_label, class_name, regulon) %>%
    summarise(
      n_SPF = sum(microbiome == "SPF"),
      n_WildR = sum(microbiome == "WildR"),
      mean_SPF = mean(auc[microbiome == "SPF"]),
      mean_WildR = mean(auc[microbiome == "WildR"]),
      delta = mean_WildR - mean_SPF,
      p = tryCatch(
        if (sum(microbiome == "SPF") >= 10 && sum(microbiome == "WildR") >= 10) {
          wilcox.test(auc ~ microbiome)$p.value
        } else {
          NA_real_
        },
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      padj = p.adjust(p, method = "BH"),
      tissue = unique(meta$tissue)[1],
      split_level = "tissue_within_class"
    ) %>%
    arrange(padj, p)
  write_csv(diff_by_class, file.path(out_dir, "diff_regulons_WildR_vs_SPF_by_class.csv"))
}

p_volc <- diff_reg %>%
  mutate(
    hit = !is.na(padj) & padj < 0.05,
    label = if_else(hit & abs(delta) >= quantile(abs(delta), 0.8, na.rm = TRUE), regulon, NA_character_)
  ) %>%
  ggplot(aes(x = delta, y = -log10(p))) +
  geom_point(aes(color = hit), alpha = 0.7, size = 1.2) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#D62728")) +
  labs(
    title = glue("WildR − SPF regulon AUC — {stem} ({split_level})"),
    x = "Δ mean AUC (WildR − SPF)", y = "−log10(p)"
  ) +
  theme_bw()
ggsave(file.path(fig_dir, "diff_regulon_volcano.png"), p_volc, width = 8, height = 6, dpi = 150)

summary <- list(
  stem = stem,
  split_level = split_level,
  tissue = unique(meta$tissue)[1],
  class_label = if (!is_tissue && "class_label" %in% names(meta)) unique(meta$class_label)[1] else NA,
  class_name = if (!is_tissue && "class_name" %in% names(meta)) unique(as.character(meta$class_name))[1] else NA,
  n_cells = length(common),
  n_regulons = nrow(auc_mat),
  n_sig_padj05 = sum(diff_reg$padj < 0.05, na.rm = TRUE),
  auc_csv = auc_csv,
  regulons_csv = regulons_csv,
  diff_csv = diff_csv,
  input_matrix = "raw_counts_UMI"
)
jsonlite::write_json(summary, summary_json, auto_unbox = TRUE, pretty = TRUE)
log_msg("DONE split ", stem, " regulons=", nrow(auc_mat),
        " sig_padj05=", summary$n_sig_padj05)
