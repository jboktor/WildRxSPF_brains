# Helpers for the integrated MAST × pySCENIC × D-SPIN × popAlign analysis
# (notebooks/0X_Integrated-Analysis_SN.qmd; plan: analysis_md/integrated_analysis_plan.md).
#
# Scaffold: popAlign units (run 2026-10-04, synthesis 2026-10-05), tissue × Allen level × label.
# Unit of replication: the animal (12 per tissue, 6 WildR / 6 SPF). Gene identity: export
# symbols (Seurat row name minus "-GRCm39"), shared by popAlign, MAST, D-SPIN and pySCENIC.

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(Matrix)
  library(data.table)
})

int_paths <- function(run_tag = "2026-10-06", wk = "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains") {
  pa_tag <- "2026-10-04"; synth <- "run_2026-10-04_synth-2026-10-05"
  list(
    wk = wk, run_tag = run_tag,
    pa_manifest = file.path(wk, glue("data/interim/popalign/snRNASeq/manifest_{pa_tag}.csv")),
    pa_tbl = file.path(wk, "tables/popAlign/snRNASeq", synth),
    mast_dir = file.path(wk, "data/interim/mast/snRNASeq"),
    dspin_dir = file.path(wk, "data/interim/dspin"),
    scenic_dir = file.path(wk, "data/interim/scenic/by_class"),
    scenic_diff = file.path(wk, "data/processed/scenic_by_class/diff_regulons_all_splits.csv"),
    vr1 = file.path(wk, "tables/VR/VR1/VR1_per_animal_scores.csv"),
    gene_meta = file.path(wk, "data/interim/seurat/snRNASeq/seurat_gene_metadata_2026-10-05.csv"),
    out = file.path(wk, "data/interim/sn_integrated", glue("run_{run_tag}")),
    tbl = file.path(wk, "tables/SN_Integrated", glue("run_{run_tag}")),
    fig = file.path(wk, "figures/SN_Integrated", glue("run_{run_tag}"))
  )
}

# MeA → BST → medial hypothalamus network (same definition as the popAlign report)
circuit_pattern <- "MEA|BST|PVH|DMH|AHN|TMd|MPN|MPO|VMH|SBPV"
in_circuit <- function(name) str_detect(name, circuit_pattern)

log_msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), paste0(...))

# ── exact permutation, vectorised over features (animals × features) ─────────
perm_combos <- function(n = 12, k = 6) combn(n, k)
# two-sided exact p for the difference in group means, for every column of X
perm_p_matrix <- function(X, is_wildr, combos = NULL) {
  X <- as.matrix(X); n <- nrow(X); k <- sum(is_wildr)
  if (k < 3 || n - k < 3) return(rep(NA_real_, ncol(X)))
  cb <- if (!is.null(combos) && n == 12 && k == 6) combos else combn(n, k)
  ind <- matrix(0, ncol(cb), n); ind[cbind(rep(seq_len(ncol(cb)), each = k), as.vector(cb))] <- 1
  tot <- colSums(X)
  s_in <- ind %*% X                                   # combos × features
  null <- s_in / k - (matrix(tot, nrow(s_in), ncol(X), byrow = TRUE) - s_in) / (n - k)
  obs <- colMeans(X[is_wildr, , drop = FALSE]) - colMeans(X[!is_wildr, , drop = FALSE])
  colMeans(abs(null) >= matrix(abs(obs), nrow(null), ncol(X), byrow = TRUE) - 1e-12)
}
welch_p <- function(x, g) tryCatch(t.test(x[g], x[!g])$p.value, error = function(e) NA_real_)

# ── loaders ──────────────────────────────────────────────────────────────────
load_mast_unit <- function(P, level, tissue, label) {
  f <- file.path(P$mast_dir, glue("{level}_label"), glue("{label}_{tissue}_2026-01-04.rds"))
  if (!file.exists(f)) return(NULL)
  readRDS(f) %>% as_tibble() %>%
    transmute(gene = sub("-GRCm39$", "", primerid), coef, p = `Pr(>Chisq)`, padj = p_adjust) %>%
    filter(!is.na(coef))
}

load_dspin_unit <- function(P, tissue, label) {
  d <- file.path(P$dspin_dir, glue("{tissue}__{label}"))
  if (!file.exists(file.path(d, "relative_responses.csv"))) return(NULL)
  meta <- read_csv(file.path(d, "sample_meta.csv"), show_col_types = FALSE)
  rr <- read_csv(file.path(d, "relative_responses.csv"), show_col_types = FALSE) %>%
    left_join(meta %>% dplyr::select(sample_id, microbiome), by = "sample_id")
  genes <- read_csv(file.path(d, "programs_genes.csv"), show_col_types = FALSE) %>%
    group_by(program_id) %>% summarise(genes = list(gene), .groups = "drop")
  tests <- rr %>% group_by(program_id) %>%
    summarise(effect = mean(relative_h[microbiome == "WildR"]) - mean(relative_h[microbiome == "SPF"]),
              p_welch = welch_p(relative_h, microbiome == "WildR"),
              star = all(relative_h[microbiome == "WildR"] > 0) | all(relative_h[microbiome == "WildR"] < 0),
              .groups = "drop")
  list(responses = rr, genes = genes, tests = left_join(tests, genes, by = "program_id"))
}

# pySCENIC regulons of one class split: activating regulons, union of targets over motifs
load_scenic_regulons <- function(P, tissue, class_label) {
  f <- file.path(P$scenic_dir, glue("snRNASeq_{tissue}_{class_label}"), "regulons.csv")
  if (!file.exists(f)) return(NULL)
  r <- data.table::fread(f, skip = 3, header = FALSE, sep = ",", quote = "\"", fill = TRUE,
                         col.names = c("TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue",
                                       "OrthologousIdentity", "Annotation", "Context", "TargetGenes", "RankAtMax"))
  r <- r[grepl("activating", Context)]
  r[, targets := str_match_all(TargetGenes, "\\('([^']+)'")]
  r %>% as_tibble() %>% group_by(TF) %>%
    summarise(targets = list(unique(unlist(map(targets, ~ .x[, 2])))), .groups = "drop") %>%
    mutate(regulon = glue("{TF}(+)"), regulon_dot = glue("{TF}.pos"))
}

# AUCell matrix of one class split (nuclei × regulons), as data.table keyed by cell
load_scenic_auc <- function(P, tissue, class_label) {
  f <- file.path(P$scenic_dir, glue("snRNASeq_{tissue}_{class_label}"), "auc_mtx.csv")
  if (!file.exists(f)) return(NULL)
  a <- data.table::fread(f); setnames(a, 1, "cell"); setkey(a, cell); a
}

# popAlign unit inputs: raw counts (genes × nuclei), cells with component assignment
load_pa_counts <- function(input_dir) {
  m <- Matrix::readMM(file.path(input_dir, "matrix.mtx")); m <- as(m, "CsparseMatrix")
  rownames(m) <- read_tsv(file.path(input_dir, "genes.tsv"), col_names = FALSE, show_col_types = FALSE)$X1
  colnames(m) <- read_tsv(file.path(input_dir, "barcodes.tsv"), col_names = FALSE, show_col_types = FALSE)$X1
  m
}

# standardised WildR − SPF difference of every oNMF feature, per component
# (same statistic as popAlign's state-shift test; animals × features centroids)
component_feature_shifts <- function(cents) {
  fcols <- grep("^F[0-9]+$", names(cents), value = TRUE)
  cents %>% group_by(comp) %>% group_modify(function(d, key) {
    X <- as.matrix(d[, fcols]); g <- d$microbiome == "WildR"
    if (sum(g) < 3 || sum(!g) < 3) return(tibble())
    sd_w <- sqrt((apply(X[g, , drop = FALSE], 2, var) * (sum(g) - 1) +
                    apply(X[!g, , drop = FALSE], 2, var) * (sum(!g) - 1)) / (nrow(X) - 2))
    tibble(feature = fcols, d = (colMeans(X[g, , drop = FALSE]) - colMeans(X[!g, , drop = FALSE])) / sd_w)
  }) %>% ungroup() %>% filter(is.finite(d))
}

hyper_p <- function(set, hits, universe) {
  set <- intersect(set, universe); hits <- intersect(hits, universe)
  k <- length(intersect(set, hits))
  if (!length(set) || !length(hits)) return(c(k = k, p = NA))
  c(k = k, p = phyper(k - 1, length(hits), length(universe) - length(hits), length(set), lower.tail = FALSE))
}

# signed z of a gene set's mean MAST coefficient vs random sets of the same size
set_direction_z <- function(genes, coef_tbl) {
  x <- coef_tbl$coef[match(genes, coef_tbl$gene)]; x <- x[!is.na(x)]
  N <- nrow(coef_tbl); n <- length(x)
  if (n < 5) return(NA_real_)
  se <- sd(coef_tbl$coef) / sqrt(n) * sqrt((N - n) / (N - 1))
  (mean(x) - mean(coef_tbl$coef)) / se
}
