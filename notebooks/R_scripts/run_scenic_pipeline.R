#!/usr/bin/env Rscript
# End-to-end pySCENIC for WildR vs SPF snRNA-seq.
#
# Usage:
#   Rscript run_scenic_pipeline.R              # pilot (~5k cells)
#   Rscript run_scenic_pipeline.R --full       # all keep_cell nuclei (no subsample)
#
# Env overrides: SCENIC_N_WORKERS, SCENIC_MIN_CELLS, SCENIC_MIN_GENES
#
# Pilot  -> data/interim/scenic, data/processed/scenic, figures/SCENIC/snRNASeq
# Full   -> data/interim/scenic_full, data/processed/scenic_full, figures/SCENIC/snRNASeq_full

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(glue)
  library(patchwork)
  library(scico)
  library(ComplexHeatmap)
  library(circlize)
  library(ggrepel)
  library(tictoc)
  library(grid)
  library(jsonlite)
})

theme_set(theme_bw())
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
full_run <- "--full" %in% args || identical(Sys.getenv("SCENIC_FULL"), "1")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
interim <- file.path(wkdir, "data/interim")
input_db <- file.path(wkdir, "data/input/scenic/cisTarget")
tag <- if (full_run) "scenic_full" else "scenic"
interim_scenic <- file.path(wkdir, "data/interim", tag)
processed <- file.path(wkdir, "data/processed", tag)
fig_dir <- file.path(
  wkdir, "figures/SCENIC",
  if (full_run) "snRNASeq_full" else "snRNASeq"
)
cluster_log_dir <- file.path(wkdir, ".cluster_runs", tag)

for (d in c(input_db, interim_scenic, processed, fig_dir, cluster_log_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# Call spatialomics binaries directly (avoid reticulate conda metadata issues)
py_bin <- "/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python"
pyscenic_bin <- "/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/pyscenic"

# Knobs
if (full_run) {
  max_cells_per_group <- NULL
  target_n_cells <- NULL
  min_genes <- as.integer(Sys.getenv("SCENIC_MIN_GENES", "200"))
  min_cells <- as.integer(Sys.getenv("SCENIC_MIN_CELLS", "50")) # stricter for full
  n_workers <- as.integer(Sys.getenv("SCENIC_N_WORKERS", "24"))
} else {
  max_cells_per_group <- 15L
  target_n_cells <- 5000L
  min_genes <- as.integer(Sys.getenv("SCENIC_MIN_GENES", "200"))
  min_cells <- as.integer(Sys.getenv("SCENIC_MIN_CELLS", "10"))
  n_workers <- as.integer(Sys.getenv("SCENIC_N_WORKERS", "6"))
}

seurat_rds <- file.path(
  interim, "seurat/snRNASeq",
  "seurat_obj_onlyparsefilters_celltyped_2025-12-27.rds"
)

expr_loom <- file.path(interim_scenic, "expr_mat.loom")
counts_mtx <- file.path(interim_scenic, "counts.mtx")
genes_tsv <- file.path(interim_scenic, "genes.tsv")
barcodes_tsv <- file.path(interim_scenic, "barcodes.tsv")
meta_csv <- file.path(interim_scenic, "cell_meta.csv")
meta_rds <- file.path(interim_scenic, "cell_meta.rds")
adj_tsv <- file.path(interim_scenic, "adjacencies.tsv")
sparse_prefix <- file.path(interim_scenic, "counts_csc")

regulons_csv <- file.path(processed, "regulons.csv")
auc_csv <- file.path(processed, "auc_mtx.csv")
auc_assay_rds <- file.path(processed, "seurat_with_scenic_auc.rds")
auc_mat_rds <- file.path(processed, "auc_mat.rds")
diff_csv <- file.path(processed, glue("diff_regulons_WildR_vs_SPF_{Sys.Date()}.csv"))
rss_csv <- file.path(processed, glue("regulon_rss_{Sys.Date()}.csv"))
summary_json <- file.path(processed, "run_summary.json")
conclusions_txt <- file.path(processed, "conclusions.txt")
pipeline_log <- file.path(interim_scenic, "pipeline.log")

log_msg <- function(...) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", paste0(...))
  cat(line, "\n")
  cat(line, "\n", file = pipeline_log, append = TRUE)
  flush.console()
}

log_msg("MODE: ", if (full_run) "FULL (all keep_cell, no subsample)" else "PILOT")
log_msg("interim=", interim_scenic)
log_msg("processed=", processed)
log_msg("figures=", fig_dir)
log_msg("n_workers=", n_workers, " min_cells=", min_cells, " min_genes=", min_genes)

# ---------------------------------------------------------------------------
# 1. Download cisTarget resources
# ---------------------------------------------------------------------------
db_urls <- c(
  "https://resources.aertslab.org/cistarget/tf_lists/allTFs_mm.txt",
  "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl",
  "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
  "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
)
options(timeout = max(12000, getOption("timeout")))
for (u in db_urls) {
  dest <- file.path(input_db, basename(u))
  if (file.exists(dest) && file.info(dest)$size > 1e6) {
    log_msg("DB exists: ", dest)
    next
  }
  # small files (TF list) may be << 1e6
  if (file.exists(dest) && grepl("allTFs|\\.tbl$", basename(dest)) && file.info(dest)$size > 1000) {
    log_msg("DB exists: ", dest)
    next
  }
  log_msg("Downloading ", u)
  download.file(u, destfile = dest, mode = "wb", quiet = FALSE)
}

tf_list <- file.path(input_db, "allTFs_mm.txt")
motif_tbl <- file.path(input_db, "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl")
rank_dbs <- list.files(input_db, pattern = "genes_vs_motifs\\.rankings\\.feather$", full.names = TRUE)
stopifnot(file.exists(tf_list), file.exists(motif_tbl), length(rank_dbs) >= 1)

# ---------------------------------------------------------------------------
# 2. Load Seurat, subsample, filter
# ---------------------------------------------------------------------------
# Resume-friendly: if AUC already exists with barcodes, reuse that cell set
resume_auc <- file.exists(auc_csv) && file.info(auc_csv)$size > 1000 &&
  file.exists(barcodes_tsv)

log_msg("Loading Seurat object")
seurat_obj <- readRDS(seurat_rds)
seurat_obj <- subset(seurat_obj, subset = keep_cell)

if (resume_auc) {
  log_msg("Resuming from existing AUC / barcodes — skipping matrix + loom rebuild")
  keep_ids <- readLines(barcodes_tsv)
  missing <- setdiff(keep_ids, colnames(seurat_obj))
  if (length(missing)) {
    stop("Barcodes not found in Seurat object: ", length(missing), " missing")
  }
  seurat_obj <- seurat_obj[, keep_ids]
  if (file.exists(meta_rds)) {
    cell_meta <- readRDS(meta_rds)
  } else {
    cell_meta <- seurat_obj@meta.data %>%
      tibble::rownames_to_column("cell_id") %>%
      select(
        cell_id, sample, microbiome, tissue, group,
        class_name, subclass_name, supertype_name, broad_cell_group
      )
    saveRDS(cell_meta, meta_rds)
    write_csv(cell_meta, meta_csv)
  }
  log_msg("Resumed cells: ", ncol(seurat_obj))
  counts <- NULL
} else {
  if (!full_run) {
    keep_cells <- seurat_obj@meta.data %>%
      tibble::rownames_to_column("cell_id") %>%
      group_by(microbiome, tissue, subclass_name) %>%
      group_modify(~ slice_sample(.x, n = min(nrow(.x), max_cells_per_group))) %>%
      ungroup()
    if (!is.null(target_n_cells) && nrow(keep_cells) > target_n_cells) {
      keep_cells <- keep_cells %>% slice_sample(n = target_n_cells)
    }
    seurat_obj <- seurat_obj[, keep_cells$cell_id]
    log_msg("Subsampled cells: ", ncol(seurat_obj))
  } else {
    log_msg("FULL run — keeping all post-QC cells: ", ncol(seurat_obj))
  }

  DefaultAssay(seurat_obj) <- "RNA"
  counts <- GetAssayData(seurat_obj, layer = "counts")
  gene_keep <- Matrix::rowSums(counts > 0) >= min_cells
  cell_keep <- Matrix::colSums(counts > 0) >= min_genes
  counts <- counts[gene_keep, cell_keep, drop = FALSE]
  seurat_obj <- seurat_obj[, colnames(counts)]

  # Parse gene names are Symbol-GRCm39; cisTarget TF lists use bare symbols
  gene_symbols <- sub("-GRCm39$", "", rownames(counts))
  keep_sym <- !duplicated(gene_symbols) & nzchar(gene_symbols) & !grepl("^ENSMUSG", gene_symbols)
  counts <- counts[keep_sym, , drop = FALSE]
  rownames(counts) <- gene_symbols[keep_sym]
  log_msg(
    "Filtered matrix: ", nrow(counts), " genes x ", ncol(counts), " cells ",
    "(symbols stripped of -GRCm39)"
  )

  tf_names <- readLines(tf_list)
  tf_overlap <- sum(rownames(counts) %in% tf_names)
  log_msg("TF overlap with allTFs_mm.txt: ", tf_overlap, " / ", length(tf_names))
  if (tf_overlap < 50) {
    stop("Too few TFs overlap gene symbols; check naming (got ", tf_overlap, ")")
  }

  cell_meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    select(
      cell_id, sample, microbiome, tissue, group,
      class_name, subclass_name, supertype_name, broad_cell_group
    )
  saveRDS(cell_meta, meta_rds)
  write_csv(cell_meta, meta_csv)
  write_tsv(tibble(gene = rownames(counts)), genes_tsv, col_names = FALSE)
  write_tsv(tibble(cell = colnames(counts)), barcodes_tsv, col_names = FALSE)

  # ---------------------------------------------------------------------------
  # 3. Build loom (Python) — sparse CSC bins for full; MTX OK for pilot
  # ---------------------------------------------------------------------------
  log_msg("Building loom")
  if (full_run) {
    counts_csc <- as(counts, "dgCMatrix")
    writeBin(as.integer(counts_csc@i), paste0(sparse_prefix, ".i.bin"), size = 4)
    writeBin(as.integer(counts_csc@p), paste0(sparse_prefix, ".p.bin"), size = 4)
    writeBin(as.double(counts_csc@x), paste0(sparse_prefix, ".x.bin"), size = 8)
    writeLines(
      as.character(c(nrow(counts_csc), ncol(counts_csc), length(counts_csc@x))),
      paste0(sparse_prefix, ".dims.txt")
    )
    rm(counts_csc)
    gc()
    py_code <- glue(
      '
from pathlib import Path
import numpy as np
import pandas as pd
import loompy
from scipy import sparse

scenic_dir = Path("{interim_scenic}")
prefix = scenic_dir / "counts_csc"
dims = [int(x) for x in open(str(prefix) + ".dims.txt")]
nrow, ncol, nnz = dims
i = np.fromfile(str(prefix) + ".i.bin", dtype=np.int32)
p = np.fromfile(str(prefix) + ".p.bin", dtype=np.int32)
x = np.fromfile(str(prefix) + ".x.bin", dtype=np.float64)
counts = sparse.csc_matrix((x, i, p), shape=(nrow, ncol)).astype(np.float32)
genes = pd.read_csv(scenic_dir / "genes.tsv", header=None)[0].astype(str).tolist()
cells = pd.read_csv(scenic_dir / "barcodes.tsv", header=None)[0].astype(str).tolist()
meta = pd.read_csv(scenic_dir / "cell_meta.csv").set_index("cell_id").loc[cells]
assert counts.shape == (len(genes), len(cells)), (counts.shape, len(genes), len(cells))
row_attrs = {{"Gene": np.array(genes, dtype=object)}}
col_attrs = {{"CellID": np.array(cells, dtype=object)}}
for col in ["sample", "microbiome", "tissue", "subclass_name", "class_name", "supertype_name"]:
    if col in meta.columns:
        col_attrs[col] = meta[col].astype(str).to_numpy()
out = scenic_dir / "expr_mat.loom"
if out.exists():
    out.unlink()
loompy.create(str(out), counts, row_attrs, col_attrs)
print("Wrote", out, "shape", counts.shape, "nnz", counts.nnz)
'
    )
  } else {
    unlink(c(expr_loom, counts_mtx))
    Matrix::writeMM(counts, file = counts_mtx)
    py_code <- glue(
      '
from pathlib import Path
import numpy as np
import pandas as pd
import loompy
from scipy.io import mmread

scenic_dir = Path("{interim_scenic}")
counts = mmread(scenic_dir / "counts.mtx").tocsc().astype(np.float32)
genes = pd.read_csv(scenic_dir / "genes.tsv", header=None)[0].astype(str).tolist()
cells = pd.read_csv(scenic_dir / "barcodes.tsv", header=None)[0].astype(str).tolist()
meta = pd.read_csv(scenic_dir / "cell_meta.csv").set_index("cell_id").loc[cells]
assert counts.shape == (len(genes), len(cells)), (counts.shape, len(genes), len(cells))
row_attrs = {{"Gene": np.array(genes, dtype=object)}}
col_attrs = {{"CellID": np.array(cells, dtype=object)}}
for col in ["sample", "microbiome", "tissue", "subclass_name", "class_name", "supertype_name"]:
    if col in meta.columns:
        col_attrs[col] = meta[col].astype(str).to_numpy()
out = scenic_dir / "expr_mat.loom"
if out.exists():
    out.unlink()
loompy.create(str(out), counts, row_attrs, col_attrs)
print("Wrote", out, "shape", counts.shape)
'
    )
  }
  tmp_py <- tempfile(fileext = ".py")
  writeLines(py_code, tmp_py)
  stopifnot(system2(py_bin, tmp_py) == 0)
}

# ---------------------------------------------------------------------------
# 4. pySCENIC CLI
# ---------------------------------------------------------------------------
run_cli <- function(args, label) {
  log_msg("START ", label, ": ", paste(args, collapse = " "))
  tic(label)
  status <- system2(pyscenic_bin, args = args)
  toc()
  if (status != 0) stop(label, " failed with status ", status)
  log_msg("DONE ", label)
}

if (!file.exists(adj_tsv) || file.info(adj_tsv)$size < 1000) {
  run_cli(
    c("grn", "--num_workers", as.character(n_workers), "-o", adj_tsv, expr_loom, tf_list),
    "pyscenic grn"
  )
} else {
  log_msg("Skip grn; exists ", adj_tsv)
}

if (!file.exists(regulons_csv) || file.info(regulons_csv)$size < 1000) {
  run_cli(
    c(
      "ctx", adj_tsv, rank_dbs,
      "--annotations_fname", motif_tbl,
      "--expression_mtx_fname", expr_loom,
      "--mode", "custom_multiprocessing",
      "--output", regulons_csv,
      "--num_workers", as.character(n_workers),
      "--mask_dropouts"
    ),
    "pyscenic ctx"
  )
} else {
  log_msg("Skip ctx; exists ", regulons_csv)
}

if (!file.exists(auc_csv) || file.info(auc_csv)$size < 1000) {
  run_cli(
    c(
      "aucell", expr_loom, regulons_csv,
      "-o", auc_csv,
      "--num_workers", as.character(n_workers)
    ),
    "pyscenic aucell"
  )
} else {
  log_msg("Skip aucell; exists ", auc_csv)
}

# ---------------------------------------------------------------------------
# 5. Import AUC into Seurat + plots
# ---------------------------------------------------------------------------
log_msg("Importing AUC")
tmp <- read.csv(auc_csv, row.names = 1, check.names = FALSE)
tmp <- as.matrix(tmp)
if (mean(rownames(tmp) %in% colnames(seurat_obj)) > 0.5) {
  auc_mat <- t(tmp) # regulons x cells
} else if (mean(colnames(tmp) %in% colnames(seurat_obj)) > 0.5) {
  auc_mat <- tmp
} else {
  stop("Could not match AUC matrix cells to Seurat colnames")
}
# Seurat-safe regulon names: Alx1(+) -> Alx1.pos
rownames(auc_mat) <- make.unique(
  gsub("\\(\\+\\)", ".pos", gsub("\\(-\\)", ".neg", rownames(auc_mat)))
)
rownames(auc_mat) <- make.unique(gsub("[^A-Za-z0-9.-]", ".", rownames(auc_mat)))
storage.mode(auc_mat) <- "numeric"

common <- intersect(colnames(seurat_obj), colnames(auc_mat))
log_msg("AUC matrix: ", nrow(auc_mat), " regulons x ", length(common), " matched cells")
stopifnot(length(common) > 0, nrow(auc_mat) > 0)
seurat_obj <- seurat_obj[, common]
auc_mat <- auc_mat[, common, drop = FALSE]

# For full cohort, build a SCENIC-only object (drop RNA counts — too large to serialize)
if (full_run) {
  meta_keep <- seurat_obj@meta.data[common, , drop = FALSE]
  seurat_obj <- CreateSeuratObject(
    counts = auc_mat,
    meta.data = meta_keep,
    assay = "SCENIC",
    project = "SCENIC_full"
  )
} else {
  seurat_obj[["SCENIC"]] <- CreateAssayObject(counts = auc_mat)
}

# AUCell scores are already enrichment metrics — put them in data and scale manually
seurat_obj <- SetAssayData(seurat_obj, assay = "SCENIC", layer = "data", new.data = auc_mat)
auc_scaled <- t(scale(t(auc_mat)))
auc_scaled[is.na(auc_scaled)] <- 0
seurat_obj <- SetAssayData(seurat_obj, assay = "SCENIC", layer = "scale.data", new.data = auc_scaled)
DefaultAssay(seurat_obj) <- "SCENIC"

VariableFeatures(seurat_obj, assay = "SCENIC") <- rownames(auc_mat)
seurat_obj <- RunPCA(
  seurat_obj, assay = "SCENIC", features = rownames(auc_mat),
  reduction.name = "scenic_pca", verbose = FALSE, npcs = min(50, nrow(auc_mat) - 1L)
)
npc <- min(30, ncol(Embeddings(seurat_obj, "scenic_pca")))
seurat_obj <- RunUMAP(
  seurat_obj, reduction = "scenic_pca", dims = 1:npc,
  reduction.name = "scenic_umap", verbose = FALSE
)
saveRDS(auc_mat, auc_mat_rds)
saveRDS(seurat_obj, auc_assay_rds)
log_msg("Saved AUC matrix + Seurat: ", auc_assay_rds)

rv <- apply(auc_scaled, 1, var)
top_reg <- names(sort(rv, decreasing = TRUE))[seq_len(min(50, length(rv)))]
top_subclasses <- seurat_obj@meta.data %>%
  count(subclass_name, sort = TRUE) %>%
  slice_head(n = 12) %>%
  pull(subclass_name)

# ---- Heatmap ----
log_msg("Plot heatmap")
pb <- as.data.frame(t(auc_scaled[top_reg, , drop = FALSE])) %>%
  mutate(
    subclass_name = seurat_obj$subclass_name,
    microbiome = seurat_obj$microbiome,
    tissue = seurat_obj$tissue
  ) %>%
  group_by(tissue, subclass_name, microbiome) %>%
  summarise(across(all_of(top_reg), mean), .groups = "drop")

mat <- pb %>%
  mutate(label = paste(tissue, subclass_name, microbiome, sep = " | ")) %>%
  column_to_rownames("label") %>%
  select(all_of(top_reg)) %>%
  as.matrix() %>%
  t()

ann <- pb %>%
  transmute(
    label = paste(tissue, subclass_name, microbiome, sep = " | "),
    tissue, microbiome, subclass_name
  ) %>%
  column_to_rownames("label")

ha <- HeatmapAnnotation(
  tissue = ann[colnames(mat), "tissue"],
  microbiome = ann[colnames(mat), "microbiome"],
  col = list(
    microbiome = c(SPF = "#4C78A8", WildR = "#F58518"),
    tissue = c(AMY = "#54A24B", HYP = "#B279A2")
  )
)
ht <- Heatmap(
  mat, name = "scaled AUC", top_annotation = ha,
  show_column_names = FALSE, cluster_columns = TRUE,
  column_split = ann[colnames(mat), "microbiome"],
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  row_names_gp = gpar(fontsize = 7)
)
heatmap_png <- file.path(fig_dir, "regulon_heatmap_top50.png")
png(heatmap_png, width = 14, height = 10, units = "in", res = 150)
draw(ht)
dev.off()

# ---- UMAP ----
log_msg("Plot UMAP")
p_umap_ct <- DimPlot(
  seurat_obj, reduction = "scenic_umap", group.by = "subclass_name",
  label = TRUE, repel = TRUE, label.size = 2
) + ggtitle("SCENIC AUC UMAP — subclass") + NoLegend()
p_umap_mb <- DimPlot(
  seurat_obj, reduction = "scenic_umap", group.by = "microbiome"
) + ggtitle("SCENIC AUC UMAP — microbiome")
example_regs <- top_reg[seq_len(min(4, length(top_reg)))]
p_feat <- FeaturePlot(
  seurat_obj, features = example_regs, reduction = "scenic_umap",
  ncol = 2, order = TRUE
) & scale_color_scico(palette = "lajolla")
p_umap <- (p_umap_ct | p_umap_mb) / p_feat
umap_png <- file.path(fig_dir, "scenic_umap_feature.png")
ggsave(umap_png, p_umap, width = 14, height = 12, dpi = 150)

# ---- Dotplot ----
log_msg("Plot dotplot")
cells_use <- colnames(seurat_obj)[seurat_obj$subclass_name %in% top_subclasses]
seu_dot <- seurat_obj[, cells_use]
seu_dot$subclass_mb <- paste(seu_dot$subclass_name, seu_dot$microbiome, sep = " | ")
p_dot <- DotPlot(
  seu_dot,
  features = top_reg[seq_len(min(25, length(top_reg)))],
  group.by = "subclass_mb", assay = "SCENIC"
) + RotatedAxis() +
  ggtitle("Top variable regulons — activity by subclass × microbiome") +
  theme(axis.text.x = element_text(size = 7))
dot_png <- file.path(fig_dir, "regulon_dotplot.png")
ggsave(dot_png, p_dot, width = 14, height = 8, dpi = 150)

# ---- RSS ----
log_msg("Plot RSS")
calc_rss <- function(auc_mat, cell_labels) {
  labels <- as.character(cell_labels)
  uct <- unique(labels)
  rss <- matrix(
    NA_real_, nrow = nrow(auc_mat), ncol = length(uct),
    dimnames = list(rownames(auc_mat), uct)
  )
  for (ct in uct) {
    in_ct <- labels == ct
    q <- as.numeric(in_ct)
    q <- q / sum(q)
    for (i in seq_len(nrow(auc_mat))) {
      p <- as.numeric(auc_mat[i, ])
      p[p < 0] <- 0
      if (sum(p) == 0) next
      p <- p / sum(p)
      m <- 0.5 * (p + q)
      js <- 0.5 * sum(ifelse(p > 0, p * log2(p / m), 0)) +
        0.5 * sum(ifelse(q > 0, q * log2(q / m), 0))
      rss[i, ct] <- 1 - sqrt(js)
    }
  }
  rss
}
auc_raw <- as.matrix(GetAssayData(seurat_obj, assay = "SCENIC", layer = "data"))
rss <- calc_rss(auc_raw, seurat_obj$subclass_name)
rss_long <- as.data.frame(rss) %>%
  rownames_to_column("regulon") %>%
  pivot_longer(-regulon, names_to = "subclass_name", values_to = "rss") %>%
  filter(subclass_name %in% top_subclasses) %>%
  group_by(subclass_name) %>%
  slice_max(rss, n = 5, with_ties = FALSE) %>%
  ungroup()
write_csv(rss_long, rss_csv)
p_rss <- rss_long %>%
  mutate(regulon = fct_reorder(regulon, rss)) %>%
  ggplot(aes(x = rss, y = regulon, fill = subclass_name)) +
  geom_col() +
  facet_wrap(~subclass_name, scales = "free_y") +
  guides(fill = "none") +
  labs(
    title = "Top cell-type–specific regulons (RSS)",
    x = "Regulon Specificity Score", y = NULL
  )
rss_png <- file.path(fig_dir, "regulon_rss.png")
ggsave(rss_png, p_rss, width = 14, height = 10, dpi = 150)

# ---- Differential regulons ----
log_msg("Differential regulons WildR vs SPF")
auc_t <- as.data.frame(t(auc_raw)) %>%
  mutate(
    microbiome = seurat_obj$microbiome,
    tissue = seurat_obj$tissue,
    subclass_name = seurat_obj$subclass_name,
    sample = seurat_obj$sample
  )
diff_reg <- auc_t %>%
  pivot_longer(cols = all_of(rownames(auc_raw)), names_to = "regulon", values_to = "auc") %>%
  group_by(tissue, subclass_name, regulon) %>%
  summarise(
    n_SPF = sum(microbiome == "SPF"),
    n_WildR = sum(microbiome == "WildR"),
    mean_SPF = mean(auc[microbiome == "SPF"]),
    mean_WildR = mean(auc[microbiome == "WildR"]),
    delta = mean_WildR - mean_SPF,
    p = tryCatch(wilcox.test(auc ~ microbiome)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  filter(n_SPF >= 10, n_WildR >= 10) %>%
  mutate(padj = p.adjust(p, method = "BH"))
write_csv(diff_reg, diff_csv)

p_volc <- diff_reg %>%
  mutate(
    hit = !is.na(padj) & padj < 0.05 & abs(delta) > quantile(abs(delta), 0.9, na.rm = TRUE),
    label = if_else(hit, regulon, NA_character_)
  ) %>%
  ggplot(aes(x = delta, y = -log10(p))) +
  geom_point(aes(color = hit), alpha = 0.5, size = 0.8) +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 30) +
  facet_wrap(~tissue) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#D62728")) +
  labs(
    title = "Differential regulon activity (WildR − SPF)",
    x = "Δ mean AUC (WildR − SPF)", y = "−log10(p) Wilcoxon"
  )
volc_png <- file.path(fig_dir, "diff_regulon_volcano.png")
ggsave(volc_png, p_volc, width = 12, height = 6, dpi = 150)

# Composition bar for subsample
p_comp <- cell_meta %>%
  count(tissue, microbiome, class_name) %>%
  ggplot(aes(x = class_name, y = n, fill = microbiome)) +
  geom_col(position = "dodge") +
  facet_wrap(~tissue, scales = "free_x") +
  scale_fill_manual(values = c(SPF = "#4C78A8", WildR = "#F58518")) +
  coord_flip() +
  labs(title = "Cells in SCENIC pilot subsample", x = NULL, y = "n cells")
comp_png <- file.path(fig_dir, "subsample_composition.png")
ggsave(comp_png, p_comp, width = 12, height = 8, dpi = 150)

# ---------------------------------------------------------------------------
# 6. Summary + conclusions
# ---------------------------------------------------------------------------
top_hits <- diff_reg %>% filter(!is.na(padj), padj < 0.05) %>% arrange(padj)
n_hits <- nrow(top_hits)
by_tissue <- top_hits %>% count(tissue, name = "n_sig")
top5 <- paste(head(top_hits$regulon, 5), collapse = ", ")
if (!nzchar(top5)) top5 <- "(none at padj<0.05 in this pilot)"

conclusion_txt <- glue(
  "SCENIC pilot run summary ({Sys.Date()})\n",
  "=====================================\n",
  "Cells: {ncol(seurat_obj)}; genes: {if (is.null(counts)) NA else nrow(counts)}; regulons: {nrow(auc_raw)}\n",
  "Subsample: max {max_cells_per_group}/microbiome×tissue×subclass, capped at {target_n_cells}\n",
  "Significant WildR vs SPF regulon tests (padj<0.05): {n_hits}\n",
  "By tissue: {paste(by_tissue$tissue, by_tissue$n_sig, sep='=', collapse=', ')}\n",
  "Example top differential regulons: {top5}\n\n",
  "Figure reads:\n",
  "1. subsample_composition — confirms balanced WildR/SPF and AMY/HYP coverage in the pilot.\n",
  "2. scenic_umap_feature — AUC embedding recovers subclass structure; microbiome mixes within types,\n",
  "   consistent with state (not identity) shifts.\n",
  "3. regulon_heatmap_top50 / regulon_dotplot — high-variance regulons mark cell-type programs;\n",
  "   within-type WildR vs SPF blocks highlight candidate state TFs.\n",
  "4. regulon_rss — top RSS regulons should enrich for known lineage TFs; use as identity sanity check.\n",
  "5. diff_regulon_volcano — microbiome-associated AUC changes; prioritize overlap with DEG notebooks.\n\n",
  "Caveat: this is a {ncol(seurat_obj)}-cell pilot. Re-run with larger n (or max_cells_per_group=NULL)\n",
  "on a compute node before final biological claims.\n"
)
writeLines(conclusion_txt, conclusions_txt)

summary <- list(
  date = as.character(Sys.Date()),
  n_cells = ncol(seurat_obj),
  n_genes = if (is.null(counts)) {
    if (file.exists(genes_tsv)) length(readLines(genes_tsv)) else NA_integer_
  } else {
    nrow(counts)
  },
  n_regulons = nrow(auc_raw),
  n_sig_diff = n_hits,
  figures = list(
    composition = comp_png,
    heatmap = heatmap_png,
    umap = umap_png,
    dotplot = dot_png,
    rss = rss_png,
    volcano = volc_png
  ),
  tables = list(
    diff = diff_csv,
    rss = rss_csv,
    seurat = auc_assay_rds,
    conclusions = conclusions_txt
  )
)
jsonlite::write_json(summary, summary_json, auto_unbox = TRUE, pretty = TRUE)
log_msg(conclusion_txt)
log_msg("PIPELINE COMPLETE")
