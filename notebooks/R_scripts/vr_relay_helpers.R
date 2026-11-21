# Shared conventions for VR-1…VR-4 olfactory–amygdalar chemosensory relay notebooks.
# Source from each notebook appendix. Never re-run MAST / D-SPIN / pySCENIC GRN.

vr_project_root <- function() {
  "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
}

vr_paths <- function(wk = vr_project_root()) {
  list(
    wk = wk,
    data = file.path(wk, "data/interim"),
    mast_signif = file.path(wk, "data/interim/mast/snRNASeq/supertype_mast_signif_hits_2026-01-04.rds"),
    mast_full_dir = file.path(wk, "data/interim/mast/snRNASeq/supertype_label"),
    dspin = file.path(wk, "data/interim/dspin"),
    stars = file.path(wk, "tables/DSPIN/star_programs.csv"),
    scenic_by_class = file.path(wk, "data/interim/scenic/by_class"),
    scenic_cell_meta = file.path(wk, "data/interim/scenic/cell_meta.csv"),
    seurat_sn = file.path(
      wk, "data/interim/seurat/snRNASeq",
      "seurat_obj_onlyparsefilters_celltyped_2025-12-27.rds"
    ),
    droplist = file.path(wk, "data/interim/seurat/snRNASeq/celltype_droplist_2025-12-30.rds"),
    spatial_seurat = file.path(
      wk, "data/interim/seurat/seqFISH",
      "merged_roi_seurat_prenorm_2026-02-16.rds"
    ),
    spatial_meta = file.path(wk, "data/interim/merged_roi_seurat_metadata_2024-07-15.rds"),
    mapmycells = file.path(wk, "data/interim/MapMyCells/output.csv"),
    ccf_rebuild = file.path(wk, "data/interim/registration/Allen_CCF_rebuild"),
    inventory = file.path(wk, "analysis_md/00_input_inventory.md")
  )
}

# Permanently excluded from signatures / module scores / rankings (§0.4)
VR_EXCLUDED_GENES <- c("ENSMUSG00000095041", "Cmss1", "Lars2", "Mid1")

strip_grcm39 <- function(x) {
  sub("-GRCm39$", "", as.character(x))
}

animal_id_from_sample <- function(sample_id) {
  # SPF_3D4_AMY / WildR_11_HYP → animal key without tissue suffix
  gsub("_(AMY|HYP)$", "", as.character(sample_id))
}

group_from_sample <- function(sample_id) {
  dplyr::case_when(
    grepl("^WildR", sample_id) ~ "WildR",
    grepl("^SPF", sample_id) ~ "SPF",
    TRUE ~ NA_character_
  )
}

# Ordinal chain positions — 3G is parallel, never in 1–4 regressions
vr_chain_map <- function() {
  tibble::tribble(
    ~supertype_pattern, ~position, ~node, ~limb,
    "^0039 |^0040 |^0041 |^0042 |L2/3 IT PIR-ENTl", 1L, "Cortical input", "glut",
    "^0058 |^0059 |COAp", 2L, "Cortical amygdala", "glut",
    "^0055 |^0056 |^0057 |^0507 |MEA Slc17a7|MEA-COA-BMA", 3L, "MeA glutamatergic", "glut",
    "^0338 |^0355 |^0359 |MEA-BST", "3G", "MeA GABAergic output", "gaba",
    "^0060 |^0061 |LA-BLA-BMA-PA", 4L, "Basolateral complex", "glut"
  )
}

assign_chain_position <- function(supertype_name) {
  nm <- as.character(supertype_name)
  pos <- rep(NA_character_, length(nm))
  # Order matters: more specific first where needed
  rules <- list(
    list(pat = "0076|CA3 Glut", pos = "EXCL_HIPPO"),
    list(pat = "MEA-BST|0338 |0355 |0359 ", pos = "3G"),
    list(pat = "0055 |0056 |0057 |0507 |MEA Slc17a7|MEA-COA-BMA", pos = "3"),
    list(pat = "0058 |0059 |COAp", pos = "2"),
    list(pat = "0039 |0040 |0041 |0042 |PIR-ENTl", pos = "1"),
    list(pat = "0060 |0061 |LA-BLA-BMA-PA", pos = "4")
  )
  for (r in rules) {
    hit <- is.na(pos) & grepl(r$pat, nm)
    pos[hit] <- r$pos
  }
  pos
}

load_mast_hits <- function(padj_cut = 0.1, paths = vr_paths()) {
  m <- readRDS(paths$mast_signif)
  tibble::as_tibble(m) %>%
    dplyr::mutate(
      gene = strip_grcm39(primerid),
      direction = dplyr::if_else(coef > 0, "WildR-UP", "SPF-UP"),
      chain_position = assign_chain_position(supertype_name),
      excluded_gene = gene %in% VR_EXCLUDED_GENES
    ) %>%
    dplyr::filter(p_adjust < padj_cut)
}

# Cliff's delta + bootstrap CI (WildR vs SPF); n_boot default 10000
cliffs_delta <- function(x_wildr, x_spf, n_boot = 10000L, seed = 1L) {
  x_wildr <- as.numeric(x_wildr)
  x_spf <- as.numeric(x_spf)
  x_wildr <- x_wildr[is.finite(x_wildr)]
  x_spf <- x_spf[is.finite(x_spf)]
  delta_one <- function(a, b) {
    if (!length(a) || !length(b)) return(NA_real_)
    comp <- outer(a, b, "-")
    (mean(comp > 0) - mean(comp < 0))
  }
  d <- delta_one(x_wildr, x_spf)
  set.seed(seed)
  boots <- replicate(n_boot, {
    delta_one(sample(x_wildr, replace = TRUE), sample(x_spf, replace = TRUE))
  })
  list(
    delta = d,
    ci_low = as.numeric(stats::quantile(boots, 0.025, na.rm = TRUE)),
    ci_high = as.numeric(stats::quantile(boots, 0.975, na.rm = TRUE)),
    n_wildr = length(x_wildr),
    n_spf = length(x_spf)
  )
}

complete_separation <- function(values, groups, higher = "WildR") {
  # How many of the higher group exceed all of the other group
  g <- as.character(groups)
  v <- as.numeric(values)
  hi <- v[g == higher]
  lo <- v[g != higher]
  if (!length(hi) || !length(lo)) {
    return(list(n_exceed_all = NA_integer_, n_hi = length(hi), n_lo = length(lo),
                complete = FALSE))
  }
  n_ex <- sum(hi > max(lo, na.rm = TRUE))
  list(
    n_exceed_all = as.integer(n_ex),
    n_hi = length(hi),
    n_lo = length(lo),
    complete = identical(as.integer(n_ex), length(hi))
  )
}

spearman_ci <- function(x, y, n_boot = 5000L, seed = 1L) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3) {
    return(list(rho = NA_real_, ci_low = NA_real_, ci_high = NA_real_, n = length(x)))
  }
  rho <- suppressWarnings(stats::cor(x, y, method = "spearman"))
  set.seed(seed)
  boots <- replicate(n_boot, {
    i <- sample.int(length(x), replace = TRUE)
    suppressWarnings(stats::cor(x[i], y[i], method = "spearman"))
  })
  list(
    rho = as.numeric(rho),
    ci_low = as.numeric(stats::quantile(boots, 0.025, na.rm = TRUE)),
    ci_high = as.numeric(stats::quantile(boots, 0.975, na.rm = TRUE)),
    n = length(x)
  )
}

wilson_ci <- function(k, n, conf = 0.95) {
  if (n <= 0) return(c(est = NA_real_, low = NA_real_, high = NA_real_))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- k / n
  den <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / den
  half <- (z / den) * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(est = p, low = max(0, centre - half), high = min(1, centre + half))
}

read_dspin_unit_h <- function(tissue, supertype_label, program_id,
                              which = c("absolute", "relative"),
                              paths = vr_paths()) {
  which <- match.arg(which)
  unit <- file.path(paths$dspin, paste0(tissue, "__", supertype_label))
  f <- if (which == "absolute") {
    file.path(unit, "responses.csv")
  } else {
    file.path(unit, "relative_responses.csv")
  }
  if (!file.exists(f)) return(NULL)
  df <- readr::read_csv(f, show_col_types = FALSE) %>%
    dplyr::filter(program_id == .env$program_id)
  hcol <- if (which == "absolute") "h" else "relative_h"
  dplyr::transmute(
    df,
    sample_id,
    animal = animal_id_from_sample(sample_id),
    group = group_from_sample(sample_id),
    tissue = tissue,
    supertype_label = supertype_label,
    program_id = program_id,
    h = .data[[hcol]]
  )
}

gaba_output_programs <- function() {
  tibble::tribble(
    ~supertype_label, ~supertype_name, ~program_id,
    "CS20230722_SUPT_0338", "0338 MEA-BST Sox6 Gaba_2", "P10",
    "CS20230722_SUPT_0355", "0355 MEA-BST Lhx6 Nr2e1 Gaba_2", "P5",
    "CS20230722_SUPT_0359", "0359 MEA-BST Lhx6 Nfib Gaba_3", "P6"
  )
}

# Gene-set score: mean expression of positives minus negatives, control-matched
# via detection-rate bins (scanpy.tl.score_genes style approximation in R).
score_genes_matched <- function(expr_mat, genes_pos, genes_neg = character(),
                                n_bins = 25L, n_ctrl = 50L, seed = 1L) {
  # expr_mat: genes × cells (dgCMatrix or matrix), rownames = gene symbols
  stopifnot(!is.null(rownames(expr_mat)))
  genes_pos <- intersect(genes_pos, rownames(expr_mat))
  genes_neg <- intersect(genes_neg, rownames(expr_mat))
  if (!length(genes_pos) && !length(genes_neg)) {
    return(rep(NA_real_, ncol(expr_mat)))
  }
  # detection rate + mean expression for matching
  if (inherits(expr_mat, "dgCMatrix") || inherits(expr_mat, "Matrix")) {
    det <- Matrix::rowMeans(expr_mat > 0)
    mu <- Matrix::rowMeans(expr_mat)
  } else {
    det <- rowMeans(expr_mat > 0)
    mu <- rowMeans(expr_mat)
  }
  set.seed(seed)
  score_one_set <- function(genes) {
    if (!length(genes)) return(rep(0, ncol(expr_mat)))
    bins <- cut(rank(det, ties.method = "average") / length(det),
                breaks = n_bins, include.lowest = TRUE, labels = FALSE)
    names(bins) <- names(det)
    ctrl <- unique(unlist(lapply(genes, function(g) {
      b <- bins[[g]]
      pool <- setdiff(names(bins)[bins == b & !(names(bins) %in% genes)], genes)
      if (!length(pool)) return(character())
      sample(pool, min(n_ctrl, length(pool)))
    })))
    if (inherits(expr_mat, "Matrix")) {
      s_g <- Matrix::colMeans(expr_mat[genes, , drop = FALSE])
      s_c <- if (length(ctrl)) Matrix::colMeans(expr_mat[ctrl, , drop = FALSE]) else 0
    } else {
      s_g <- colMeans(expr_mat[genes, , drop = FALSE])
      s_c <- if (length(ctrl)) colMeans(expr_mat[ctrl, , drop = FALSE]) else 0
    }
    as.numeric(s_g - s_c)
  }
  score_one_set(genes_pos) - score_one_set(genes_neg)
}

jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (!length(a) && !length(b)) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

write_md <- function(path, lines) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
}

vr_brand_ggplot <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#FFF3EB", colour = NA),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      legend.position = "top"
    )
}

VR_COLORS <- c(WildR = "#FF6C0C", SPF = "#27455C")

# Gene scoring via Python subprocess (reticulate conda probe is broken on this cluster)
vr_python <- function() {
  "/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python"
}

vr_score_py <- function(wk = vr_project_root()) {
  file.path(wk, "notebooks/python_scripts/score_genes_h5ad.py")
}

score_unit_via_python <- function(supertype_label, genes_pos, genes_neg = character(),
                                  paths = vr_paths()) {
  h5 <- file.path(paths$dspin, paste0("AMY__", supertype_label), "raw_counts.h5ad")
  meta_f <- file.path(paths$dspin, paste0("AMY__", supertype_label), "cells_meta.csv")
  if (!file.exists(h5) || !file.exists(meta_f)) return(NULL)
  tmp <- tempfile("vr_score_")
  pos_f <- paste0(tmp, "_pos.txt")
  neg_f <- paste0(tmp, "_neg.txt")
  out_f <- paste0(tmp, "_out.csv")
  on.exit(unlink(c(pos_f, neg_f, out_f)), add = TRUE)
  writeLines(as.character(genes_pos), pos_f)
  if (length(genes_neg)) {
    writeLines(as.character(genes_neg), neg_f)
    neg_arg <- neg_f
  } else {
    neg_arg <- "-"
  }
  rc <- system2(
    vr_python(),
    c(vr_score_py(paths$wk), h5, meta_f, pos_f, neg_arg, out_f),
    stdout = TRUE, stderr = TRUE
  )
  if (!file.exists(out_f)) {
    warning("score failed for ", supertype_label, ": ", paste(rc, collapse = " | "))
    return(NULL)
  }
  scores <- readr::read_csv(out_f, show_col_types = FALSE)
  meta <- readr::read_csv(meta_f, show_col_types = FALSE)
  n <- min(nrow(scores), nrow(meta))
  tibble::tibble(
    cell_id = as.character(scores$cell_id[seq_len(n)]),
    sample_id = as.character(scores$sample_id[seq_len(n)]),
    animal = animal_id_from_sample(scores$sample_id[seq_len(n)]),
    group = group_from_sample(scores$sample_id[seq_len(n)]),
    microbiome = meta$microbiome[seq_len(n)],
    score_pos_only = as.numeric(scores$score_pos_only[seq_len(n)]),
    score_signed = as.numeric(scores$score_signed[seq_len(n)]),
    score = as.numeric(scores$score_pos_only[seq_len(n)]),
    genes_used = length(genes_pos),
    panel = NA_character_,
    supertype_label = supertype_label,
    median_umi_proxy = NA_real_
  )
}
