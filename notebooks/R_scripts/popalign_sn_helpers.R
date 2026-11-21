# Helpers for the snRNA-seq popAlign synthesis (0X_popAlign_SN.qmd).
# Levels: tissue × class, tissue × subclass, tissue × supertype.
# Statistics are animal-level throughout: 12 animals per tissue (6 SPF, 6 WildR),
# exact permutation over all choose(12, 6) = 924 group relabelings. The smallest
# attainable two-sided p is 2/924 ≈ 0.0022, so effect sizes and complete
# separation are reported alongside BH q-values, plus a continuous quasi-binomial
# GLM p as a secondary test.

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
})

pa_wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
pa_levels <- c("class", "subclass", "supertype")
# next-finer Allen label used for component purity, per level
pa_finer <- c(class = "subclass_label", subclass = "supertype_label", supertype = "cluster_label")

pa_paths <- function(run_tag = "2026-10-04", wk = pa_wkdir) {
  base <- file.path(wk, "data/interim/popalign/snRNASeq")
  list(
    run_tag = run_tag,
    base = base,
    manifest = file.path(base, glue("manifest_{run_tag}.csv")),
    results = file.path(base, glue("results_{run_tag}")),
    synth = file.path(base, glue("synthesis_{run_tag}")),
    fig = file.path(wk, "figures/popAlign/snRNASeq", glue("run_{run_tag}")),
    tbl = file.path(wk, "tables/popAlign/snRNASeq", glue("run_{run_tag}")),
    dspin = file.path(wk, "data/interim/dspin"),
    # Seurat RNA meta.features: feature, gene_id (Ensembl), gene_name, gene_type, export_symbol
    gene_meta = file.path(wk, "data/interim/seurat/snRNASeq/seurat_gene_metadata_2026-10-05.csv"),
    mast = file.path(wk, "data/interim/mast/snRNASeq",
                     glue("{pa_levels}_mast_signif_hits_2026-01-04.rds")) %>% set_names(pa_levels)
  )
}

pa_colors <- c(WildR = "#FF6C0C", SPF = "#27455C")

# ── loading ──────────────────────────────────────────────────────────────────
read_unit_csv <- function(man, file, ...) {
  rd <- function(f) read_csv(f, show_col_types = FALSE, ...)  # dots captured here, not in a ~ lambda
  man %>%
    mutate(path = file.path(output_dir, file)) %>%
    filter(file.exists(path)) %>%
    dplyr::select(unit_id, level, tissue, label, name, path) %>%
    mutate(data = map(path, rd)) %>%
    dplyr::select(-path) %>%
    # header-only files (e.g. no markers in K = 1 units) parse as all-character
    # columns and would clash with typed columns from other units: drop them
    filter(map_int(data, nrow) > 0) %>%
    unnest(data)
}

read_status <- function(man) {
  num <- function(st, key) map_dbl(st, ~ as.numeric(.x$steps[[key]] %||% NA))
  lgl <- function(st, key) map_lgl(st, ~ as.logical(.x$steps[[key]] %||% NA))
  man %>%
    mutate(status_f = file.path(output_dir, "status.json")) %>%
    mutate(st = map(status_f, ~ if (file.exists(.x)) jsonlite::fromJSON(.x) else NULL),
           status = map_chr(st, ~ .x$status %||% "missing"),
           beta = num(st, "beta"),
           atlas_K = num(st, "atlas_K"),
           atlas_pool = num(st, "atlas_pool"),
           boot_K_agreement = num(st, "boot_K_agreement"),
           n_features = num(st, "n_features"),
           n_features_at_bound = lgl(st, "n_features_at_bound"),
           atlas_K_at_upper_bound = lgl(st, "atlas_K_at_upper_bound"),
           plot_deltas = map_chr(st, ~ substr(.x$steps$plot_deltas %||% NA_character_, 1, 40)),
           runtime_min = map_dbl(st, ~ .x$runtime_min %||% NA_real_)) %>%
    dplyr::select(-st, -status_f)
}

# ── exact permutation machinery (animal level) ───────────────────────────────
perm_index <- function(n = 12, k = 6) combn(n, k)

perm_diff_test <- function(x, grp, combos = NULL) {
  # x: one value per animal; grp: "WildR"/"SPF". Two-sided exact p for diff of means.
  ok <- !is.na(x)
  x <- x[ok]; g <- grp[ok] == "WildR"
  n <- length(x); k <- sum(g)
  if (k < 3 || (n - k) < 3) return(tibble(diff = NA_real_, p_perm = NA_real_))
  obs <- mean(x[g]) - mean(x[!g])
  cb <- if (!is.null(combos) && n == 12 && k == 6) combos else combn(n, k)
  tot <- sum(x)
  sums <- colSums(matrix(x[cb], nrow = k))
  null <- sums / k - (tot - sums) / (n - k)
  tibble(diff = obs, p_perm = mean(abs(null) >= abs(obs) - 1e-12))
}

# multivariate version: distance between group mean vectors (rows = animals)
perm_dist_test <- function(X, grp, combos = NULL) {
  g <- grp == "WildR"
  if (sum(g) < 3 || sum(!g) < 3) return(tibble(stat = NA_real_, p_perm = NA_real_))
  stat <- function(gi) sqrt(sum((colMeans(X[gi, , drop = FALSE]) - colMeans(X[!gi, , drop = FALSE]))^2))
  obs <- stat(g)
  cb <- if (!is.null(combos) && nrow(X) == 12 && sum(g) == 6) combos else combn(nrow(X), sum(g))
  null <- apply(cb, 2, function(ix) stat(seq_len(nrow(X)) %in% ix))
  tibble(stat = obs, p_perm = mean(null >= obs - 1e-12))
}

hedges_g <- function(x, grp) {
  a <- x[grp == "WildR"]; b <- x[grp == "SPF"]
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  if (length(a) < 2 || length(b) < 2) return(NA_real_)
  sp <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2))
  J <- 1 - 3 / (4 * (length(a) + length(b)) - 9)
  J * (mean(a) - mean(b)) / sp
}

# number of animals of the "higher" group that exceed every animal of the other
# group (6 = complete separation)
separation <- function(x, grp) {
  a <- x[grp == "WildR"]; b <- x[grp == "SPF"]
  max(sum(a > max(b)), sum(b > max(a)), sum(a < min(b)), sum(b < min(a)))
}

# ── 1. abundance (Δw) ────────────────────────────────────────────────────────
# (a) unit omnibus: CLR-transformed per-animal fraction vectors, distance between
#     group means, exact permutation; BH across units within level.
# (b) component tests on logit fractions. In K = 2 units the two fractions are
#     complementary, so only the minority component is tested (the other is
#     marked redundant and excluded from BH).
# adepth: per-animal unit median log10 UMI (animal_depth()); adds a depth-adjusted
# quasi-binomial p so depth-associated components can still be judged on
# whether the microbiome effect survives the ~5% WildR/SPF depth difference.
test_abundance <- function(fracs, adepth = NULL, eps = 1e-3) {
  combos <- perm_index()
  if (!is.null(adepth)) fracs <- fracs %>% left_join(adepth, by = c("unit_id", "sample_id"))
  if (!"animal_log10_umi" %in% names(fracs)) fracs$animal_log10_umi <- NA_real_
  multi <- fracs %>% group_by(unit_id) %>% filter(n_distinct(comp) >= 2) %>% ungroup()

  omni <- multi %>%
    mutate(f = pmax(soft_frac, eps)) %>%
    group_by(unit_id, level, tissue, label, name, sample_id, microbiome) %>%
    mutate(clr = log(f) - mean(log(f))) %>% ungroup() %>%
    dplyr::select(unit_id, level, tissue, label, name, sample_id, microbiome, comp, clr) %>%
    pivot_wider(names_from = comp, values_from = clr, names_prefix = "c") %>%
    group_by(unit_id, level, tissue, label, name) %>%
    group_modify(function(d, key) {
      X <- as.matrix(d %>% dplyr::select(starts_with("c")))
      X <- X[, colSums(is.na(X)) == 0, drop = FALSE]
      perm_dist_test(X, d$microbiome, combos) %>% rename(omnibus_stat = stat, p_omnibus = p_perm)
    }) %>% ungroup() %>%
    group_by(level) %>% mutate(q_omnibus = p.adjust(p_omnibus, "BH")) %>% ungroup()

  comp <- multi %>%
    group_by(unit_id) %>%
    mutate(K = n_distinct(comp),
           minority = comp == comp[which.min(ave(soft_frac, comp, FUN = mean))]) %>%
    ungroup() %>%
    mutate(logit_frac = qlogis(pmin(pmax(soft_frac, eps), 1 - eps)),
           soft_n = round(soft_frac * n_cells)) %>%
    group_by(unit_id, level, tissue, label, name, comp, K, minority) %>%
    group_modify(function(d, key) {
      pt <- perm_diff_test(d$logit_frac, d$microbiome, combos)
      glm_p <- function(f) tryCatch({
        fit <- glm(f, family = quasibinomial, data = d)
        coef(summary(fit))["microbiomeWildR", "Pr(>|t|)"]
      }, error = function(e) NA_real_)
      p_glm <- glm_p(cbind(soft_n, n_cells - soft_n) ~ microbiome)
      p_glm_depth <- if (all(is.na(d$animal_log10_umi))) NA_real_ else
        glm_p(cbind(soft_n, n_cells - soft_n) ~ microbiome + animal_log10_umi)
      tibble(
        mean_frac_SPF = mean(d$soft_frac[d$microbiome == "SPF"]),
        mean_frac_WildR = mean(d$soft_frac[d$microbiome == "WildR"]),
        log2fc_frac = log2((mean_frac_WildR + eps) / (mean_frac_SPF + eps)),
        delta_logit = pt$diff, g = hedges_g(d$logit_frac, d$microbiome),
        separation = separation(d$logit_frac, d$microbiome),
        p_perm = pt$p_perm, p_glm = p_glm, p_glm_depth = p_glm_depth
      )
    }) %>% ungroup() %>%
    mutate(redundant = K == 2 & !minority) %>%
    group_by(level) %>%
    mutate(q_perm = replace(rep(NA_real_, n()), !redundant, p.adjust(p_perm[!redundant], "BH")),
           q_glm = replace(rep(NA_real_, n()), !redundant, p.adjust(p_glm[!redundant], "BH"))) %>%
    ungroup() %>%
    left_join(omni %>% dplyr::select(unit_id, p_omnibus, q_omnibus), by = "unit_id")

  list(unit = omni, comp = comp)
}

# ── 2. state shift (Δμ) within each component ────────────────────────────────
# Euclidean distance between WildR and SPF mean centroids after standardising each
# oNMF feature by its pooled within-group SD; exact permutation. For K = 1 units
# this is a whole-population state shift.
test_state_shift <- function(cents, min_eff_n = 5) {
  combos <- perm_index()
  cents %>%
    filter(eff_n >= min_eff_n) %>%
    group_by(unit_id, level, tissue, label, name, comp) %>%
    group_modify(function(d, key) {
      fcols <- grep("^F[0-9]+$", names(d), value = TRUE)
      X <- as.matrix(d[, fcols]); X <- X[, colSums(!is.na(X)) == nrow(X), drop = FALSE]
      g <- d$microbiome == "WildR"
      if (sum(g) < 3 || sum(!g) < 3 || ncol(X) == 0) return(tibble(shift = NA_real_, p_perm = NA_real_, n_animals = nrow(d)))
      sd_w <- sqrt((apply(X[g, , drop = FALSE], 2, var) * (sum(g) - 1) +
                      apply(X[!g, , drop = FALSE], 2, var) * (sum(!g) - 1)) / (nrow(X) - 2))
      sd_w[!is.finite(sd_w) | sd_w == 0] <- NA
      Z <- sweep(X, 2, sd_w, "/"); Z <- Z[, !is.na(sd_w), drop = FALSE]
      pt <- perm_dist_test(Z, d$microbiome, combos)
      dvec <- colMeans(Z[g, , drop = FALSE]) - colMeans(Z[!g, , drop = FALSE])
      tibble(shift = pt$stat, p_perm = pt$p_perm, n_animals = nrow(Z),
             top_shift_feature = names(dvec)[which.max(abs(dvec))],
             top_shift_value = dvec[which.max(abs(dvec))])
    }) %>%
    ungroup() %>%
    group_by(level) %>% mutate(q_perm = p.adjust(p_perm, "BH")) %>% ungroup()
}

# ── 3. native popAlign alignment summaries ───────────────────────────────────
summarise_alignments <- function(al) {
  grp <- al %>% mutate(microbiome = if_else(str_detect(sample_id, "^SPF"), "SPF", "WildR"))
  n_rep <- max(al$rep) + 1
  per_animal <- grp %>%
    group_by(unit_id, sample_id, microbiome, atlas_comp) %>%
    # popAlign's JeffreyDiv returns log10(JD); average on that scale
    summarise(present_frac = n_distinct(rep) / n_rep, mean_log10_jd = mean(jd), .groups = "drop")
  combos <- perm_index()
  per_animal %>%
    group_by(unit_id, atlas_comp) %>%
    group_modify(function(d, key) {
      pt <- perm_diff_test(d$mean_log10_jd, d$microbiome, combos)
      tibble(n_animals_aligned = n_distinct(d$sample_id),
             present_SPF = mean(d$present_frac[d$microbiome == "SPF"]),
             present_WildR = mean(d$present_frac[d$microbiome == "WildR"]),
             log10_jd_SPF = median(d$mean_log10_jd[d$microbiome == "SPF"]),
             log10_jd_WildR = median(d$mean_log10_jd[d$microbiome == "WildR"]),
             delta_log10_jd = pt$diff, p_perm_jd = pt$p_perm)
    }) %>% ungroup() %>%
    rename(comp = atlas_comp)
}

# ── 4. Allen composition of components (interpretability anchor) ─────────────
allen_composition <- function(man) {
  map_dfr(seq_len(nrow(man)), function(i) {
    r <- man[i, ]; f <- file.path(r$output_dir, "cells.csv.gz")
    if (!file.exists(f)) return(NULL)
    lab_col <- pa_finer[[r$level]]
    cells <- read_csv(f, col_select = any_of(c("comp", lab_col)), show_col_types = FALSE)
    if (!lab_col %in% names(cells)) return(NULL)
    cells %>%
      count(comp, allen = .data[[lab_col]]) %>%
      group_by(comp) %>%
      mutate(p = n / sum(n)) %>%
      summarise(dominant_allen = allen[which.max(p)],
                purity = max(p),
                entropy_bits = -sum(p * log2(p)),
                n_allen_10pct = sum(p >= 0.10),
                .groups = "drop") %>%
      mutate(unit_id = r$unit_id, allen_level = sub("_label$", "", lab_col))
  })
}

# ── 5. functional annotation (GO BP on component markers) ────────────────────
# markers are specificity-ranked by the runner (log2FC, Δpct, AUC >= 0.6)
component_marker_sets <- function(markers, n_top = 50, min_auc = 0.6) {
  markers %>%
    filter(auc >= min_auc) %>%
    group_by(unit_id, comp) %>%
    slice_max(auc, n = n_top, with_ties = FALSE) %>%
    summarise(genes = list(gene), n_markers = n(), .groups = "drop")
}

go_empty <- tibble(unit_id = character(), comp = integer(), ID = character(), Description = character(),
                   GeneRatio = character(), BgRatio = character(), p.adjust = double(), geneID = character())

# GO is run on Ensembl IDs (gene_id from the Seurat feature metadata): symbols are
# ambiguous here because 288 features share a gene_name with another Ensembl ID.
# `universe_by_unit` must hold Ensembl IDs; readable = TRUE reports symbols.
go_enrich_components <- function(sets, universe_by_unit, ont = "BP", min_markers = 10) {
  suppressPackageStartupMessages({ library(clusterProfiler); library(org.Mm.eg.db) })
  res <- sets %>%
    filter(n_markers >= min_markers) %>%
    mutate(go = map2(ids, unit_id, function(g, u) {
      res <- tryCatch(
        enrichGO(g, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = ont,
                 universe = universe_by_unit[[u]], pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                 minGSSize = 10, maxGSSize = 500, readable = TRUE),
        error = function(e) NULL)
      if (is.null(res) || nrow(as.data.frame(res)) == 0) return(tibble())
      as_tibble(as.data.frame(res)) %>% dplyr::select(ID, Description, GeneRatio, BgRatio, p.adjust, geneID) %>%
        slice_head(n = 15)
    })) %>%
    dplyr::select(unit_id, comp, go) %>%
    unnest(go)
  bind_rows(go_empty, res)   # keep the schema even when nothing is enriched
}

# ── 6. overlap with D-SPIN programs and MAST hits ────────────────────────────
dspin_overlap <- function(sets, man, paths = pa_paths()) {
  st <- man %>% filter(level == "supertype")   # D-SPIN units are supertypes
  empty <- tibble(unit_id = character(), comp = integer(), best_dspin_program = character(),
                  dspin_overlap = integer(), dspin_jaccard = double())
  res <- map_dfr(seq_len(nrow(st)), function(i) {
    r <- st[i, ]
    f <- file.path(paths$dspin, glue("{r$tissue}__{r$label}"), "programs_genes.csv")
    if (!file.exists(f)) return(NULL)
    prog <- read_csv(f, show_col_types = FALSE) %>% group_by(program_id) %>% summarise(pg = list(gene))
    s <- sets %>% filter(unit_id == r$unit_id)
    if (!nrow(s)) return(NULL)
    tidyr::expand_grid(s %>% dplyr::select(unit_id, comp, genes), prog) %>%
      mutate(n_overlap = map2_int(genes, pg, ~ length(intersect(.x, .y))),
             jaccard = map2_dbl(genes, pg, ~ length(intersect(.x, .y)) / length(union(.x, .y)))) %>%
      dplyr::select(unit_id, comp, program_id, n_overlap, jaccard)
  })
  if (!nrow(res)) return(empty)
  res %>%
    group_by(unit_id, comp) %>% slice_max(jaccard, n = 1, with_ties = FALSE) %>% ungroup() %>%
    rename(best_dspin_program = program_id, dspin_jaccard = jaccard, dspin_overlap = n_overlap)
}

mast_overlap <- function(sets, man, universe_by_unit, paths = pa_paths(), padj = 0.1) {
  load_hits <- function(f) {
    if (!file.exists(f)) return(NULL)
    readRDS(f) %>% as_tibble() %>% filter(p_adjust < padj) %>%
      transmute(tissue, label = as.character(cell_type), gene = sub("-GRCm39$", "", primerid), coef)
  }
  hits <- map_dfr(paths$mast, load_hits)
  empty <- tibble(unit_id = character(), comp = integer(), n_mast_hits_unit = integer(),
                  n_markers_mast = integer(), p_mast_enrich = double(), mast_genes = character())
  if (!nrow(sets)) return(empty)
  res <- sets %>%
    left_join(man %>% dplyr::select(unit_id, tissue, label), by = "unit_id") %>%
    mutate(res = pmap(list(genes, tissue, label, unit_id), function(g, t, l, u) {
      h <- hits %>% filter(tissue == t, label == l)
      U <- universe_by_unit[[u]]
      hg <- intersect(h$gene, U); g <- intersect(g, U)
      k <- length(intersect(g, hg))
      p <- if (length(hg) && length(g)) phyper(k - 1, length(hg), length(U) - length(hg), length(g), lower.tail = FALSE) else NA
      tibble(n_mast_hits_unit = length(hg), n_markers_mast = k, p_mast_enrich = p,
             mast_genes = paste(intersect(g, hg), collapse = ";"))
    })) %>%
    dplyr::select(unit_id, comp, res) %>% unnest(res)
  bind_rows(empty, res)
}

# ── 7. cross-unit meta-programs ──────────────────────────────────────────────
meta_programs <- function(sets, min_jaccard = 0.15, resolution = 1) {
  suppressPackageStartupMessages(library(igraph))
  s <- sets %>% mutate(node = glue("{unit_id}::{comp}"))
  empty_genes <- tibble(meta_program = integer(), n_members = integer(), gene = character(),
                        n_with_gene = integer(), frac_members = double())
  G <- s$genes; n <- length(G)
  ed <- list()
  if (n >= 2) for (i in seq_len(n - 1)) for (j in (i + 1):n) {
    jac <- length(intersect(G[[i]], G[[j]])) / length(union(G[[i]], G[[j]]))
    if (jac >= min_jaccard) ed[[length(ed) + 1]] <- c(i, j, jac)
  }
  if (!length(ed)) return(list(membership = s %>% mutate(meta_program = NA_integer_), genes = empty_genes))
  ed <- do.call(rbind, ed)
  g <- graph_from_data_frame(data.frame(from = s$node[ed[, 1]], to = s$node[ed[, 2]], weight = ed[, 3]),
                             directed = FALSE, vertices = data.frame(name = s$node))
  set.seed(42)
  cl <- cluster_leiden(g, objective_function = "modularity", weights = E(g)$weight, resolution_parameter = resolution)
  mem <- s %>% mutate(meta_program = as.integer(membership(cl)[node])) %>%
    add_count(meta_program, name = "mp_size") %>%
    mutate(meta_program = if_else(mp_size >= 3, meta_program, NA_integer_))
  mp_genes <- mem %>% filter(!is.na(meta_program)) %>%
    dplyr::select(meta_program, unit_id, gene = genes) %>% unnest(gene) %>%
    group_by(meta_program) %>% mutate(n_members = n_distinct(unit_id)) %>%
    count(meta_program, n_members, gene, name = "n_with_gene") %>%
    mutate(frac_members = n_with_gene / n_members) %>%
    filter(frac_members >= 0.3) %>% arrange(meta_program, desc(frac_members)) %>% ungroup()
  list(membership = mem, genes = bind_rows(empty_genes, mp_genes), graph = g)
}

# ── 8. coarse ↔ fine consistency (shared nuclei) ─────────────────────────────
# For each coarse-level component: how its nuclei distribute over the components
# of the finer-level units (class → subclass, subclass → supertype).
cross_level_map <- function(man, coarse = "class", fine = "subclass") {
  fine_col <- c(subclass = "subclass_label", supertype = "supertype_label")[[fine]]
  cu <- man %>% filter(level == coarse)
  map_dfr(seq_len(nrow(cu)), function(i) {
    r <- cu[i, ]; f <- file.path(r$output_dir, "cells.csv.gz")
    if (!file.exists(f)) return(NULL)
    cc <- read_csv(f, col_select = any_of(c("cell_barcode", "comp", fine_col)), show_col_types = FALSE)
    if (!fine_col %in% names(cc)) return(NULL)
    fu <- man %>% filter(level == fine, tissue == r$tissue, label %in% unique(cc[[fine_col]]))
    sc <- map_dfr(seq_len(nrow(fu)), function(j) {
      fj <- file.path(fu$output_dir[j], "cells.csv.gz")
      if (!file.exists(fj)) return(NULL)
      read_csv(fj, col_select = c(cell_barcode, comp), show_col_types = FALSE) %>%
        transmute(cell_barcode, fine_unit = fu$unit_id[j], fine_comp = comp)
    })
    if (!nrow(sc)) return(NULL)
    inner_join(cc, sc, by = "cell_barcode") %>%
      count(coarse_comp = comp, fine_unit, fine_comp) %>%
      group_by(coarse_comp) %>% mutate(p = n / sum(n)) %>% ungroup() %>%
      mutate(coarse_unit = r$unit_id, coarse_level = coarse, fine_level = fine)
  })
}

# per-animal median depth within each unit (covariate for test_abundance)
animal_depth <- function(man) {
  map_dfr(seq_len(nrow(man)), function(i) {
    r <- man[i, ]; f <- file.path(r$output_dir, "cells.csv.gz")
    if (!file.exists(f)) return(NULL)
    cc <- read_csv(f, col_select = any_of(c("sample_id", "log10_umi")), show_col_types = FALSE)
    if (!"log10_umi" %in% names(cc)) return(NULL)
    cc %>% group_by(sample_id) %>% summarise(animal_log10_umi = median(log10_umi), .groups = "drop") %>%
      mutate(unit_id = r$unit_id)
  })
}

# ── 9. sequencing-depth QC of components ─────────────────────────────────────
# snRNA components can split nuclei by capture depth (long-gene signature) rather
# than biology. Per unit: eta^2 = share of per-nucleus log10 UMI variance explained
# by component; per component: median log10 UMI shift vs the rest of the unit.
# Uses the runner's per-nucleus log10_umi (raw counts).
depth_qc <- function(man) {
  map_dfr(seq_len(nrow(man)), function(i) {
    r <- man[i, ]; f_cells <- file.path(r$output_dir, "cells.csv.gz")
    if (!file.exists(f_cells)) return(NULL)
    cc <- read_csv(f_cells, col_select = any_of(c("comp", "log10_umi", "sample_id", "microbiome")),
                   show_col_types = FALSE)
    if (!"log10_umi" %in% names(cc)) return(NULL)
    mu <- mean(cc$log10_umi)
    eta2 <- if (n_distinct(cc$comp) > 1)
      sum(tapply(cc$log10_umi, cc$comp, function(v) length(v) * (mean(v) - mu)^2)) / sum((cc$log10_umi - mu)^2) else 0
    cc %>% group_by(comp) %>%
      summarise(median_log10_umi = median(log10_umi), .groups = "drop") %>%
      mutate(delta_log10_umi = map_dbl(comp, function(j) {
        rest <- cc$log10_umi[cc$comp != j]
        if (length(rest)) median(cc$log10_umi[cc$comp == j]) - median(rest) else 0
      }),
      unit_id = r$unit_id, unit_depth_eta2 = eta2)
  })
}


# ── 10. gene identity: export symbols → Ensembl IDs / clean symbols ──────────
load_gene_map <- function(paths = pa_paths()) {
  read_csv(paths$gene_meta, show_col_types = FALSE) %>%
    dplyr::select(export_symbol, gene_id, gene_name, any_of("gene_type")) %>%
    distinct(export_symbol, .keep_all = TRUE)
}

# add Ensembl IDs (`ids`) and clean display symbols (`names`) to marker sets
annotate_sets <- function(sets, gmap) {
  lk_id <- setNames(gmap$gene_id, gmap$export_symbol)
  lk_nm <- setNames(gmap$gene_name, gmap$export_symbol)
  sets %>% mutate(ids = map(genes, ~ unname(na.omit(lk_id[.x]))),
                  names = map(genes, ~ unname(coalesce(lk_nm[.x], .x))))
}

# ── 11. cross-lineage contamination flag ─────────────────────────────────────
# A component whose top markers carry ANOTHER lineage's canonical genes is more
# likely ambient RNA / doublets than a state of the unit's own cell type (e.g.
# myelin genes in hypothalamic neuronal units; neuronal genes in OPC units).
contam_markers <- list(
  oligodendrocyte = c("Plp1", "Mbp", "Mobp", "Mog", "Mag", "Cldn11", "Opalin", "Ermn", "St18"),
  # Lhfpl3 / Megf11 dropped: also expressed by some neurons (false OPC flags in GABA units)
  OPC = c("Pdgfra", "Cspg4", "Gpr17", "Olig1", "Olig2", "Sox10"),
  astrocyte = c("Slc1a2", "Slc1a3", "Atp1a2", "Slc6a11", "Gja1", "Aqp4", "Gfap", "Aldh1l1", "Sparcl1"),
  microglia = c("Hexb", "Inpp5d", "Cx3cr1", "P2ry12", "Csf1r", "Tmem119", "C1qa", "C1qb"),
  endothelial = c("Flt1", "Cldn5", "Ly6c1", "Pecam1", "Vwf"),
  neuronal = c("Syt1", "Snap25", "Rbfox3", "Snhg11", "Celf4", "Rian", "Stmn2", "Syn2")
)

# lineage(s) that a unit legitimately expresses, from its Allen name
unit_lineages <- function(name) {
  n <- tolower(name)
  out <- character()
  if (grepl("oligo|mfol|mol|cop", n)) out <- c(out, "oligodendrocyte", "OPC")
  if (grepl("opc", n)) out <- c(out, "OPC", "oligodendrocyte")
  if (grepl("astro|epen|tanycyte", n)) out <- c(out, "astrocyte")
  if (grepl("immune|microglia|bam|macrophage", n)) out <- c(out, "microglia")
  if (grepl("vascular|endo|peri|smc|vlmc", n)) out <- c(out, "endothelial")
  # neuronal if it names a transmitter type, or if no glial/vascular lineage matched
  if (grepl("glut|gaba|dopa|chol|sero|glyc|hist", n) || length(out) == 0) out <- c(out, "neuronal")
  unique(out)
}

contamination_flags <- function(sets, man, min_genes = 2) {
  lin <- man %>% distinct(unit_id, name) %>% mutate(own = map(name, unit_lineages))
  sets %>%
    left_join(lin, by = "unit_id") %>%
    mutate(res = map2(names, own, function(g, own) {
      hits <- map_int(contam_markers, ~ length(intersect(g, .x)))
      foreign <- hits[!names(hits) %in% own & hits >= min_genes]
      tibble(contamination_flag = length(foreign) > 0,
             contamination_type = if (length(foreign)) paste(names(sort(foreign, decreasing = TRUE)), collapse = "+") else NA_character_,
             contamination_genes = if (length(foreign)) paste(intersect(g, unlist(contam_markers[names(foreign)])), collapse = ";") else NA_character_)
    })) %>%
    dplyr::select(unit_id, comp, res) %>% unnest(res)
}
