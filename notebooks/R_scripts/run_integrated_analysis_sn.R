#!/usr/bin/env Rscript
# Integrated MAST × pySCENIC × D-SPIN × popAlign analysis (snRNA-seq, WildR vs SPF).
# Plan: analysis_md/integrated_analysis_plan.md · report: notebooks/0X_Integrated-Analysis_SN.qmd
#
# Stages (each writes to tables/SN_Integrated/run_<tag>/ and data/interim/sn_integrated/run_<tag>/):
#   gene     gene-level bridge: popAlign shifted program vs MAST / D-SPIN / pySCENIC (all levels)
#   cell     cell-level bridge (supertype units): per animal × subpopulation regulon AUC,
#            pseudobulk limma-voom (unit, per component, heterogeneity), D-SPIN program activity
#   animal   animal-level integration: leave-one-animal-out modality scores, concordance, VR link
#   circuit  MeA → BST → medial hypothalamus enrichment, every modality
#   evidence integrated evidence table per supertype unit and per popAlign lead
#   figures  report figures
#   all      everything in order
#
# Usage: Rscript notebooks/R_scripts/run_integrated_analysis_sn.R --stage all [--run-tag 2026-10-06]
#        sbatch notebooks/shell_scripts/run_integrated_analysis_sn.sh

suppressPackageStartupMessages({ library(tidyverse); library(glue); library(edgeR); library(limma) })
source("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/notebooks/R_scripts/integrated_sn_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[[i + 1]] }
stage <- get_arg("--stage", "all")
only_units <- get_arg("--units", NULL)          # testing: comma-separated unit_ids
P <- int_paths(get_arg("--run-tag", "2026-10-06"))
if (!is.null(get_arg("--out"))) { o <- get_arg("--out"); P$out <- file.path(o, "interim"); P$tbl <- file.path(o, "tables"); P$fig <- file.path(o, "figures") }
walk(c(P$out, P$tbl, P$fig), dir.create, recursive = TRUE, showWarnings = FALSE)

man <- read_csv(P$pa_manifest, show_col_types = FALSE)
catalog <- read_csv(file.path(P$pa_tbl, "popAlign_subpopulation_catalog.csv"), show_col_types = FALSE)
scenic_diff <- read_csv(P$scenic_diff, show_col_types = FALSE)
if (!is.null(only_units)) man <- man %>% filter(unit_id %in% str_split(only_units, ",")[[1]])
unit_class <- read_csv(file.path(P$pa_tbl, "popAlign_subpopulation_catalog.csv"), show_col_types = FALSE) %>%
  distinct(unit_id, level, tissue, label, name)
# parent class of each unit (from its nuclei)
parent_class <- function(u) {
  read_csv(file.path(u$output_dir, "cells.csv.gz"), col_select = "class_label", n_max = 1, show_col_types = FALSE)$class_label
}
combos <- perm_combos()

# ════════════════════════════════════════════════════════════════════════════
# Stage 1 · gene-level bridge
# ════════════════════════════════════════════════════════════════════════════
stage_gene <- function() {
  ok_comps <- catalog %>% filter(robust, !contamination_flag) %>% dplyr::select(unit_id, comp, p_state, name)
  units <- man %>% semi_join(ok_comps, by = "unit_id")
  log_msg("Gene bridge: ", nrow(units), " units")
  scen_cache <- list()
  res <- map_dfr(seq_len(nrow(units)), function(i) {
    u <- units[i, ]
    cents <- read_csv(file.path(u$output_dir, "animal_component_centroids.csv"), show_col_types = FALSE)
    sh <- component_feature_shifts(cents)
    W <- read_csv(file.path(u$output_dir, "onmf_W.csv.gz"), show_col_types = FALSE)
    mast <- load_mast_unit(P, u$level, u$tissue, u$label)
    dsp <- if (u$level == "supertype") load_dspin_unit(P, u$tissue, u$label) else NULL
    cls <- parent_class(u); key <- glue("{u$tissue}_{cls}")
    if (is.null(scen_cache[[key]])) scen_cache[[key]] <<- load_scenic_regulons(P, u$tissue, cls) %||% tibble()
    reg <- scen_cache[[key]]
    U <- W$gene; if (!is.null(mast)) U <- intersect(U, mast$gene)
    mastU <- if (!is.null(mast)) mast %>% filter(gene %in% U) else NULL
    degs <- if (!is.null(mastU)) mastU$gene[mastU$padj < 0.1] else character()
    top_genes <- function(f) W %>% arrange(desc(.data[[f]])) %>% filter(gene %in% U) %>% slice_head(n = 50) %>% pull(gene)
    comps <- ok_comps %>% filter(unit_id == u$unit_id)
    map_dfr(comps$comp, function(cc) {
      s <- sh %>% filter(comp == cc); if (!nrow(s)) return(NULL)
      fs <- s$feature[which.max(abs(s$d))]; ds <- s$d[which.max(abs(s$d))]
      fc <- s$feature[which.min(abs(s$d))]; dc <- s$d[which.min(abs(s$d))]
      out <- tibble(unit_id = u$unit_id, level = u$level, tissue = u$tissue, name = u$name, comp = cc,
                    p_state = comps$p_state[comps$comp == cc], circuit = in_circuit(u$name),
                    shifted_feature = fs, d_shift = ds, control_feature = fc, d_control = dc,
                    n_universe = length(U), n_deg = length(degs))
      for (arm in c("shift", "control")) {
        g <- top_genes(if (arm == "shift") fs else fc); dd <- if (arm == "shift") ds else dc
        # MAST: DEG overlap and signed direction agreement
        if (!is.null(mastU)) {
          hp <- hyper_p(g, degs, U)
          out[[glue("mast_overlap_{arm}")]] <- hp[["k"]]; out[[glue("mast_p_{arm}")]] <- hp[["p"]]
          out[[glue("mast_dirz_{arm}")]] <- sign(dd) * set_direction_z(g, mastU)
        }
        # D-SPIN: best-overlapping program, its response and direction
        if (!is.null(dsp)) {
          ov <- dsp$tests %>% mutate(hp = map(genes, ~ hyper_p(g, .x, U)), k = map_dbl(hp, "k"), p = map_dbl(hp, "p")) %>%
            arrange(p, desc(k))
          b <- ov[1, ]
          out[[glue("dspin_best_{arm}")]] <- b$program_id; out[[glue("dspin_overlap_p_{arm}")]] <- b$p
          out[[glue("dspin_best_resp_p_{arm}")]] <- b$p_welch; out[[glue("dspin_best_star_{arm}")]] <- b$star
          out[[glue("dspin_dir_concord_{arm}")]] <- sign(b$effect) == sign(dd)
          out[[glue("dspin_any_responsive_overlap_{arm}")]] <- any(ov$p < 0.05 & ov$p_welch < 0.05, na.rm = TRUE)
        }
        # pySCENIC: regulon target enrichment of the program genes (BH within component × arm)
        if (nrow(reg)) {
          en <- reg %>% mutate(hp = map(targets, ~ hyper_p(g, .x, U)), k = map_dbl(hp, "k"), p = map_dbl(hp, "p")) %>%
            mutate(q = p.adjust(p, "BH")) %>% arrange(p)
          top <- en[1, ]
          dif <- scenic_diff %>% filter(tissue == u$tissue, class_label == cls, regulon == as.character(top$regulon_dot))
          out[[glue("scenic_n_enriched_{arm}")]] <- sum(en$q < 0.05, na.rm = TRUE)
          out[[glue("scenic_top_{arm}")]] <- as.character(top$TF); out[[glue("scenic_top_q_{arm}")]] <- top$q
          out[[glue("scenic_top_diff_padj_{arm}")]] <- if (nrow(dif)) dif$padj[1] else NA_real_
          out[[glue("scenic_top_diff_delta_{arm}")]] <- if (nrow(dif)) dif$delta[1] else NA_real_
        }
      }
      out
    })
  })
  write_csv(res, file.path(P$tbl, "integrated_gene_bridge.csv"))
  log_msg("Gene bridge rows: ", nrow(res))
  invisible(res)
}

# ════════════════════════════════════════════════════════════════════════════
# Stage 2 · cell-level bridge (supertype units)
# ════════════════════════════════════════════════════════════════════════════
voom_fit <- function(counts, design) {
  y <- DGEList(counts); keep <- filterByExpr(y, design = design); y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y); v <- voom(y, design); list(v = v, fit = eBayes(lmFit(v, design)))
}
pb_sum <- function(cm, groups) {   # genes × nuclei → genes × groups (sum)
  f <- factor(groups); ind <- sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), nlevels(f)))
  out <- as.matrix(cm %*% ind); colnames(out) <- levels(f); out
}

stage_cell <- function() {
  units <- man %>% filter(level == "supertype")
  units$class_label <- map_chr(seq_len(nrow(units)), ~ parent_class(units[.x, ]))
  log_msg("Cell bridge: ", nrow(units), " supertype units in ", n_distinct(paste(units$tissue, units$class_label)), " pySCENIC splits")
  out <- list(reg_unit = list(), reg_comp = list(), pb_unit = list(), pb_genes = list(), pb_comp = list(),
              het = list(), dsp_comp = list(), dsp_val = list(), animal_auc = list(), logcpm = list())
  for (sp in split(units, paste(units$tissue, units$class_label))) {
    auc <- load_scenic_auc(P, sp$tissue[1], sp$class_label[1])
    for (i in seq_len(nrow(sp))) {
      u <- sp[i, ]; uid <- u$unit_id
      cells <- read_csv(file.path(u$output_dir, "cells.csv.gz"), col_select = c(cell_barcode, sample_id, microbiome, comp), show_col_types = FALSE)
      animals <- cells %>% distinct(sample_id, microbiome) %>% arrange(microbiome, sample_id)
      gW <- animals$microbiome == "WildR"
      robust_comps <- catalog %>% filter(unit_id == uid, robust, !contamination_flag) %>% pull(comp)
      # ── regulon AUC per animal (unit) and per animal × component ──
      if (!is.null(auc)) {
        a <- auc[J(cells$cell_barcode), nomatch = 0L]; reg_cols <- setdiff(names(a), "cell")
        a <- as_tibble(a) %>% left_join(cells, by = c("cell" = "cell_barcode"))
        cov <- nrow(a) / nrow(cells)
        am <- a %>% group_by(sample_id) %>% summarise(across(all_of(reg_cols), mean), n = n(), .groups = "drop")
        am <- animals %>% left_join(am, by = "sample_id")
        if (all(!is.na(am$n))) {
          X <- as.matrix(am[, reg_cols])
          out$animal_auc[[uid]] <- am %>% mutate(unit_id = uid)
          out$reg_unit[[uid]] <- tibble(unit_id = uid, regulon = reg_cols,
            effect = colMeans(X[gW, , drop = FALSE]) - colMeans(X[!gW, , drop = FALSE]),
            p_welch = apply(X, 2, welch_p, g = gW), p_perm = perm_p_matrix(X, gW, combos), coverage = cov)
        }
        for (cc in sort(unique(cells$comp))) {
          ac <- a %>% filter(comp == cc) %>% group_by(sample_id) %>% filter(n() >= 5) %>%
            summarise(across(all_of(reg_cols), mean), .groups = "drop") %>% left_join(animals, by = "sample_id")
          g <- ac$microbiome == "WildR"
          if (sum(g) >= 4 && sum(!g) >= 4) {
            X <- as.matrix(ac[, reg_cols])
            out$reg_comp[[glue("{uid}::{cc}")]] <- tibble(unit_id = uid, comp = cc, robust = cc %in% robust_comps, regulon = reg_cols,
              effect = colMeans(X[g, , drop = FALSE]) - colMeans(X[!g, , drop = FALSE]),
              p_perm = perm_p_matrix(X, g, combos), n_animals = nrow(ac))
          }
        }
      }
      # ── pseudobulk: unit level ──
      cm <- load_pa_counts(u$input_dir)[, cells$cell_barcode]
      pbu <- pb_sum(cm, cells$sample_id)[, animals$sample_id]
      des <- model.matrix(~ factor(animals$microbiome, levels = c("SPF", "WildR")))
      vf <- voom_fit(pbu, des); tt <- topTable(vf$fit, coef = 2, number = Inf, sort.by = "none") %>% rownames_to_column("gene")
      out$logcpm[[uid]] <- vf$v$E
      mast <- load_mast_unit(P, "supertype", u$tissue, u$label)
      mj <- if (!is.null(mast)) inner_join(tt, mast, by = "gene") else NULL
      mdeg <- if (!is.null(mj)) mj %>% filter(padj < 0.1) else NULL
      out$pb_genes[[uid]] <- tt %>% transmute(unit_id = uid, gene, logFC, t, p = P.Value, padj = adj.P.Val)
      out$pb_unit[[uid]] <- tibble(unit_id = uid, n_genes = nrow(tt), n_fdr10 = sum(tt$adj.P.Val < 0.1), n_p01 = sum(tt$P.Value < 0.01),
        excess_p01 = sum(tt$P.Value < 0.01) / (0.01 * nrow(tt)), min_p = min(tt$P.Value),
        n_mast_deg = if (!is.null(mdeg)) nrow(mdeg) else NA_integer_,
        mast_deg_pb_confirmed = if (!is.null(mdeg) && nrow(mdeg)) mean(mdeg$P.Value < 0.05 & sign(mdeg$logFC) == sign(mdeg$coef)) else NA_real_,
        rho_mast_pb = if (!is.null(mj) && nrow(mj) > 10) cor(mj$coef, mj$logFC, method = "spearman") else NA_real_)
      # ── pseudobulk: per component, and heterogeneity of the microbiome effect across components ──
      cnt <- cells %>% count(sample_id, comp)
      elig <- cnt %>% filter(n >= 5) %>% left_join(animals, by = "sample_id") %>% group_by(comp) %>%
        filter(sum(microbiome == "WildR") >= 4, sum(microbiome == "SPF") >= 4) %>% ungroup()
      if (nrow(elig)) {
        key <- paste(cells$sample_id, cells$comp, sep = "::")
        pbc <- pb_sum(cm, key); smp <- elig %>% mutate(col = paste(sample_id, comp, sep = "::")) %>% filter(col %in% colnames(pbc))
        pbc <- pbc[, smp$col, drop = FALSE]
        dsp <- load_dspin_unit(P, u$tissue, u$label)
        for (cc in unique(smp$comp)) {
          s <- smp %>% filter(comp == cc); dd <- model.matrix(~ factor(s$microbiome, levels = c("SPF", "WildR")))
          vc <- tryCatch(voom_fit(pbc[, s$col, drop = FALSE], dd), error = function(e) NULL); if (is.null(vc)) next
          tc <- topTable(vc$fit, coef = 2, number = Inf, sort.by = "none")
          out$pb_comp[[glue("{uid}::{cc}")]] <- tibble(unit_id = uid, comp = cc, robust = cc %in% robust_comps, n_animals = nrow(s),
            n_genes = nrow(tc), n_fdr10 = sum(tc$adj.P.Val < 0.1), n_p01 = sum(tc$P.Value < 0.01), excess_p01 = sum(tc$P.Value < 0.01) / (0.01 * nrow(tc)))
          # D-SPIN program activity within the component (mean z-scored logCPM of program genes)
          if (!is.null(dsp)) {
            E <- vc$v$E; Z <- t(scale(t(E))); gw <- s$microbiome == "WildR"
            ps <- dsp$genes %>% mutate(score = map(genes, ~ { gg <- intersect(.x, rownames(Z)); if (length(gg) >= 3) colMeans(Z[gg, , drop = FALSE]) else NULL })) %>%
              filter(!map_lgl(score, is.null))
            if (nrow(ps)) {
              S <- do.call(cbind, ps$score); colnames(S) <- ps$program_id
              out$dsp_comp[[glue("{uid}::{cc}")]] <- tibble(unit_id = uid, comp = cc, robust = cc %in% robust_comps, program_id = ps$program_id,
                effect = colMeans(S[gw, , drop = FALSE]) - colMeans(S[!gw, , drop = FALSE]), p_perm = perm_p_matrix(S, gw, combos))
            }
          }
        }
        if (n_distinct(smp$comp) >= 2) {
          f_comp <- factor(smp$comp); wild <- as.integer(smp$microbiome == "WildR")
          dh <- model.matrix(~ 0 + f_comp + f_comp:wild); colnames(dh) <- make.names(colnames(dh))
          y <- DGEList(pbc); keep <- filterByExpr(y, design = dh); y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
          het <- tryCatch({
            v <- voom(y, dh); dc <- duplicateCorrelation(v, dh, block = smp$sample_id)
            fit <- lmFit(v, dh, block = smp$sample_id, correlation = dc$consensus)
            ic <- grep("wild", colnames(dh), value = TRUE)
            cm_ <- makeContrasts(contrasts = paste(ic[-1], "-", ic[1]), levels = dh)
            ft <- topTable(eBayes(contrasts.fit(fit, cm_)), number = Inf, sort.by = "none")
            tibble(unit_id = uid, n_comps = nlevels(f_comp), n_genes = nrow(ft), n_het_fdr10 = sum(ft$adj.P.Val < 0.1),
                   n_het_p01 = sum(ft$P.Value < 0.01), excess_het_p01 = sum(ft$P.Value < 0.01) / (0.01 * nrow(ft)),
                   frac_mast_deg_heterogeneous = if (!is.null(mdeg) && nrow(mdeg)) mean(ft[intersect(mdeg$gene, rownames(ft)), "P.Value"] < 0.05) else NA_real_,
                   dup_corr = dc$consensus)
          }, error = function(e) NULL)
          out$het[[uid]] <- het
        }
      }
      # ── D-SPIN validation: per-animal pseudobulk program score vs D-SPIN relative h ──
      dsp <- load_dspin_unit(P, u$tissue, u$label)
      if (!is.null(dsp)) {
        Z <- t(scale(t(vf$v$E)))
        out$dsp_val[[uid]] <- dsp$genes %>% mutate(rho = map2_dbl(genes, program_id, function(gg, pid) {
          gg <- intersect(gg, rownames(Z)); if (length(gg) < 3) return(NA_real_)
          sc <- colMeans(Z[gg, , drop = FALSE]); h <- dsp$responses %>% filter(program_id == pid)
          cor(sc[h$sample_id], h$relative_h, method = "spearman")
        })) %>% transmute(unit_id = uid, program_id, rho)
      }
      log_msg("  ", uid, " done (", nrow(cells), " nuclei)")
    }
    rm(auc); gc()
  }
  save_tbl <- function(x, f) { d <- bind_rows(x); write_csv(d, file.path(P$tbl, f)); d }
  save_tbl(out$reg_unit, "integrated_regulon_animal_unit.csv"); save_tbl(out$reg_comp, "integrated_regulon_animal_component.csv")
  save_tbl(out$pb_unit, "integrated_pseudobulk_unit_summary.csv"); save_tbl(out$pb_comp, "integrated_pseudobulk_component_summary.csv")
  save_tbl(out$het, "integrated_heterogeneity_summary.csv"); save_tbl(out$dsp_comp, "integrated_dspin_program_component.csv")
  save_tbl(out$dsp_val, "integrated_dspin_validation.csv")
  saveRDS(bind_rows(out$pb_genes), file.path(P$out, "pseudobulk_unit_genes.rds"))
  saveRDS(out$animal_auc, file.path(P$out, "animal_regulon_auc.rds")); saveRDS(out$logcpm, file.path(P$out, "pseudobulk_unit_logcpm.rds"))
  log_msg("Cell bridge complete")
}

# ════════════════════════════════════════════════════════════════════════════
# Stage 3 · animal-level integration (leave-one-animal-out modality scores)
# ════════════════════════════════════════════════════════════════════════════
# For a matrix animals × features: z-score each feature across animals, then for every
# animal take the WildR − SPF direction learned from the OTHER 11 animals and score the
# held-out animal along it. The animal never informs its own direction.
loao_score <- function(X, is_wildr) {
  X <- as.matrix(X); X <- X[, apply(X, 2, function(v) sd(v, na.rm = TRUE) > 0), drop = FALSE]
  if (!ncol(X)) return(rep(NA_real_, nrow(X)))
  Z <- scale(X)
  sapply(seq_len(nrow(Z)), function(a) {
    w <- colMeans(Z[-a, ][is_wildr[-a], , drop = FALSE]) - colMeans(Z[-a, ][!is_wildr[-a], , drop = FALSE])
    sum(Z[a, ] * sign(w) * abs(w)) / sum(abs(w))
  })
}

stage_animal <- function() {
  auc <- readRDS(file.path(P$out, "animal_regulon_auc.rds")); lc <- readRDS(file.path(P$out, "pseudobulk_unit_logcpm.rds"))
  sup <- man %>% filter(level == "supertype") %>% mutate(circuit = in_circuit(name))
  scope_scores <- function(scope_units, tis) {
    us <- scope_units %>% filter(tissue == tis)
    per_unit <- map_dfr(seq_len(nrow(us)), function(i) {
      u <- us[i, ]; res <- list()
      # popAlign: robust clean components, state centroids
      rc <- catalog %>% filter(unit_id == u$unit_id, robust, !contamination_flag) %>% pull(comp)
      if (length(rc)) {
        ce <- read_csv(file.path(u$output_dir, "animal_component_centroids.csv"), show_col_types = FALSE) %>% filter(comp %in% rc)
        sc <- ce %>% group_by(comp) %>% filter(n() == 12) %>% group_modify(function(d, k) {
          d <- arrange(d, microbiome, sample_id); f <- grep("^F[0-9]+$", names(d), value = TRUE)
          tibble(sample_id = d$sample_id, microbiome = d$microbiome, s = loao_score(d[, f], d$microbiome == "WildR"))
        }) %>% ungroup() %>% group_by(sample_id, microbiome) %>% summarise(s = mean(s), .groups = "drop")
        if (nrow(sc)) res$popAlign <- sc
      }
      # D-SPIN programs
      dsp <- load_dspin_unit(P, u$tissue, u$label)
      if (!is.null(dsp)) {
        w <- dsp$responses %>% dplyr::select(sample_id, microbiome, program_id, relative_h) %>% pivot_wider(names_from = program_id, values_from = relative_h) %>% arrange(microbiome, sample_id)
        res$`D-SPIN` <- tibble(sample_id = w$sample_id, microbiome = w$microbiome, s = loao_score(w[, -(1:2)], w$microbiome == "WildR"))
      }
      # pySCENIC regulon AUC (animal means)
      a <- auc[[u$unit_id]]
      if (!is.null(a)) {
        a <- arrange(a, microbiome, sample_id); rcol <- setdiff(names(a), c("sample_id", "microbiome", "n", "unit_id"))
        res$pySCENIC <- tibble(sample_id = a$sample_id, microbiome = a$microbiome, s = loao_score(a[, rcol], a$microbiome == "WildR"))
      }
      # pseudobulk expression: LOAO top-50 genes by Welch t, signed (MAST-equivalent signature)
      E <- lc[[u$unit_id]]
      if (!is.null(E)) {
        mic <- if_else(str_detect(colnames(E), "^WildR"), "WildR", "SPF"); g <- mic == "WildR"
        Z <- t(scale(t(E))); Z <- Z[is.finite(rowSums(Z)), , drop = FALSE]
        s_pb <- sapply(seq_len(ncol(Z)), function(a) {
          Zo <- Z[, -a]; go <- g[-a]
          t <- (rowMeans(Zo[, go]) - rowMeans(Zo[, !go])) / sqrt(apply(Zo[, go], 1, var) / sum(go) + apply(Zo[, !go], 1, var) / sum(!go))
          top <- order(-abs(t))[1:min(50, length(t))]
          mean(Z[top, a] * sign(t[top]))
        })
        res$Expression <- tibble(sample_id = colnames(E), microbiome = mic, s = s_pb)
      }
      imap_dfr(res, ~ .x %>% mutate(modality = .y, unit_id = u$unit_id))
    })
    if (!nrow(per_unit)) return(tibble())
    per_unit %>% group_by(modality, unit_id) %>% mutate(s = as.numeric(scale(s))) %>% ungroup() %>%
      group_by(modality, sample_id, microbiome) %>% summarise(score = mean(s, na.rm = TRUE), n_units = n_distinct(unit_id), .groups = "drop")
  }
  scores <- map_dfr(c("AMY", "HYP"), function(tis) bind_rows(
    scope_scores(sup %>% filter(circuit), tis) %>% mutate(scope = "network"),
    scope_scores(sup %>% filter(!circuit), tis) %>% mutate(scope = "other")) %>% mutate(tissue = tis))
  # VR relay scores (AMY): animal id without tissue suffix
  vr <- read_csv(P$vr1, show_col_types = FALSE) %>%
    transmute(sample_id = paste0(animal, "_AMY"), microbiome = group, `VR input` = input_score, `VR output` = output_score_mean, `VR D` = D) %>%
    pivot_longer(starts_with("VR"), names_to = "modality", values_to = "score") %>% mutate(scope = "network", tissue = "AMY", n_units = NA)
  scores <- bind_rows(scores, vr)
  write_csv(scores, file.path(P$tbl, "integrated_animal_scores.csv"))

  # LOAO prediction: does the direction learned from 11 animals place the held-out one correctly?
  # Leave-one-out scores are biased AGAINST the held-out animal's label under the null (removing an
  # animal pulls its group mean away from it), so only correct placement counts: one-sided test of
  # AUC > 0.5 against all 924 relabelings. The true null sits below 0.5, so this is conservative.
  auc_of <- function(score, w) { r <- rank(score); (sum(r[w]) - sum(w) * (sum(w) + 1) / 2) / (sum(w) * sum(!w)) }
  pred <- scores %>% filter(!str_detect(modality, "^VR")) %>% group_by(tissue, scope, modality) %>%
    summarise(n = n(), accuracy = mean(sign(score) == if_else(microbiome == "WildR", 1, -1)),
              auc = auc_of(score, microbiome == "WildR"),
              p_one_sided = { w <- microbiome == "WildR"; obs <- auc_of(score, w)
                null <- apply(combos, 2, function(ix) { ww <- seq_along(score) %in% ix; auc_of(score, ww) })
                mean(null >= obs - 1e-12) },
              .groups = "drop")
  write_csv(pred, file.path(P$tbl, "integrated_animal_loao_prediction.csv"))

  # cross-modality coupling: all animals, and within group (residual after removing group mean)
  coup <- scores %>% group_by(tissue, scope) %>% group_modify(function(d, k) {
    w <- d %>% dplyr::select(sample_id, microbiome, modality, score) %>% pivot_wider(names_from = modality, values_from = score) %>% drop_na()
    mods <- setdiff(names(w), c("sample_id", "microbiome")); if (length(mods) < 2) return(tibble())
    R <- w %>% group_by(microbiome) %>% mutate(across(all_of(mods), ~ .x - mean(.x))) %>% ungroup()
    pairs <- combn(mods, 2, simplify = FALSE)
    map_dfr(pairs, function(pr) {
      rho_all <- cor(w[[pr[1]]], w[[pr[2]]], method = "spearman")
      x <- R[[pr[1]]]; y <- R[[pr[2]]]; rho_w <- cor(x, y, method = "spearman")
      set.seed(1); null <- replicate(5000, { yy <- y; for (gg in unique(R$microbiome)) { ix <- which(R$microbiome == gg); yy[ix] <- y[sample(ix)] }; cor(x, yy, method = "spearman") })
      tibble(m1 = pr[1], m2 = pr[2], n = nrow(w), rho_all = rho_all, rho_within = rho_w, p_within = mean(abs(null) >= abs(rho_w) - 1e-12))
    })
  }) %>% ungroup()
  write_csv(coup, file.path(P$tbl, "integrated_animal_coupling.csv"))
  log_msg("Animal-level integration complete")
}

# ════════════════════════════════════════════════════════════════════════════
# Stage 4 · circuit enrichment across modalities (supertype units)
# ════════════════════════════════════════════════════════════════════════════
stage_circuit <- function() {
  # Every modality is judged the same, calibrated way: per supertype unit, the fraction of its
  # tests with p < 0.05 ("hit rate"). Main contrast: network vs other units, ratio of mean hit
  # rates, tissue-stratified permutation p (10,000) and bootstrap 95% CI. Secondary: a binary
  # "more hits than chance" flag (one-sided binomial at 0.05) in a logistic model adjusted for
  # unit size and tissue. (A flag like "any of 20–240 tests < 0.05" saturates by chance and is
  # not used.) MAST is nuclei-level: its rate is the fraction of tested genes with padj < 0.1.
  sup <- man %>% filter(level == "supertype") %>% mutate(circuit = in_circuit(name))
  pg <- readRDS(file.path(P$out, "pseudobulk_unit_genes.rds"))
  regu <- read_csv(file.path(P$tbl, "integrated_regulon_animal_unit.csv"), show_col_types = FALSE)
  rate_tbl <- function(d, mod, inference) d %>% mutate(modality = mod, inference = inference)
  pa <- catalog %>% filter(level == "supertype", robust, !contamination_flag, !is.na(p_state)) %>% group_by(unit_id) %>%
    summarise(n_tests = n(), n_hit = sum(p_state < 0.05), .groups = "drop") %>% rate_tbl("popAlign state", "animals")
  pb <- pg %>% group_by(unit_id) %>% summarise(n_tests = n(), n_hit = sum(p < 0.05), .groups = "drop") %>% rate_tbl("Pseudobulk expression", "animals")
  sc <- regu %>% filter(!is.na(p_perm)) %>% group_by(unit_id) %>% summarise(n_tests = n(), n_hit = sum(p_perm < 0.05), .groups = "drop") %>% rate_tbl("pySCENIC regulons", "animals")
  ds <- map_dfr(seq_len(nrow(sup)), function(i) { d <- load_dspin_unit(P, sup$tissue[i], sup$label[i]); if (is.null(d)) return(NULL)
    tibble(unit_id = sup$unit_id[i], n_tests = sum(!is.na(d$tests$p_welch)), n_hit = sum(d$tests$p_welch < 0.05, na.rm = TRUE),
           dspin_any_star = any(d$tests$star)) }) %>% rate_tbl("D-SPIN programs", "animals")
  ma <- map_dfr(seq_len(nrow(sup)), function(i) { m <- load_mast_unit(P, "supertype", sup$tissue[i], sup$label[i]); if (is.null(m)) return(NULL)
    tibble(unit_id = sup$unit_id[i], n_tests = nrow(m), n_hit = sum(m$padj < 0.1, na.rm = TRUE)) }) %>% rate_tbl("MAST DEGs (nuclei)", "nuclei")
  rates <- bind_rows(pa, pb, sc, ds %>% dplyr::select(-dspin_any_star), ma) %>%
    left_join(sup %>% dplyr::select(unit_id, tissue, name, circuit, n_cells_exported), by = "unit_id") %>%
    mutate(rate = n_hit / n_tests,
           alpha = if_else(modality == "MAST DEGs (nuclei)", NA_real_, 0.05),
           excess_flag = if_else(is.na(alpha), NA, n_hit > qbinom(0.95, n_tests, 0.05)),
           any_flag = n_hit > 0)
  write_csv(rates, file.path(P$tbl, "integrated_unit_rates.csv"))

  set.seed(42)
  res <- rates %>% group_by(modality, inference) %>% group_modify(function(d, k) {
    obs <- mean(d$rate[d$circuit]) / mean(d$rate[!d$circuit])
    null <- replicate(10000, { cc <- d$circuit; for (t in unique(d$tissue)) { ix <- which(d$tissue == t); cc[ix] <- sample(cc[ix]) }
      mean(d$rate[cc]) / mean(d$rate[!cc]) })
    boot <- replicate(2000, { a <- sample(d$rate[d$circuit], replace = TRUE); b <- sample(d$rate[!d$circuit], replace = TRUE); mean(a) / mean(b) })
    # size-adjusted test on the hit COUNTS (quasi-binomial: hits out of tests per unit), so units with
    # more nuclei (more power, esp. for MAST) do not drive the contrast. exp(coef) ≈ adjusted rate ratio
    # for small rates.
    m <- tryCatch(glm(cbind(n_hit, n_tests - n_hit) ~ circuit + log(n_cells_exported) + tissue,
                      family = quasibinomial, data = d), error = function(e) NULL)
    flag <- if (k$modality == "popAlign state") d$any_flag else d$excess_flag
    tibble(n_units = nrow(d), n_circuit = sum(d$circuit), n_other = sum(!d$circuit),
           rate_circuit = mean(d$rate[d$circuit]), rate_other = mean(d$rate[!d$circuit]),
           rate_ratio = obs, ci_lo = quantile(boot, 0.025, na.rm = TRUE), ci_hi = quantile(boot, 0.975, na.rm = TRUE),
           p_perm = (sum(abs(log(null)) >= abs(log(obs)) - 1e-12, na.rm = TRUE) + 1) / (sum(is.finite(null)) + 1),
           flag_circuit = sum(flag[d$circuit], na.rm = TRUE), flag_other = sum(flag[!d$circuit], na.rm = TRUE),
           ratio_adj = if (!is.null(m)) exp(coef(m)[["circuitTRUE"]]) else NA_real_,
           ratio_adj_lo = if (!is.null(m)) exp(coef(m)[["circuitTRUE"]] - 1.96 * coef(summary(m))["circuitTRUE", 2]) else NA_real_,
           ratio_adj_hi = if (!is.null(m)) exp(coef(m)[["circuitTRUE"]] + 1.96 * coef(summary(m))["circuitTRUE", 2]) else NA_real_,
           p_adj = if (!is.null(m)) coef(summary(m))["circuitTRUE", 4] else NA_real_)
  }) %>% ungroup()
  write_csv(res, file.path(P$tbl, "integrated_circuit_enrichment.csv"))

  # unit flags for the evidence table (calibrated)
  U <- sup %>% dplyr::select(unit_id, tissue, name, circuit, n_cells_exported) %>%
    left_join(pa %>% transmute(unit_id, pa_n_comp = n_tests, pa_state_hit = n_hit > 0), by = "unit_id") %>%
    left_join(pb %>% transmute(unit_id, n_genes = n_tests, pb_n_p05 = n_hit, pb_excess_hit = n_hit > qbinom(0.95, n_tests, 0.05)), by = "unit_id") %>%
    left_join(ds %>% transmute(unit_id, dspin_n_prog = n_tests, dspin_n_p05 = n_hit, dspin_any_star, dspin_hit = n_hit > qbinom(0.95, n_tests, 0.05)), by = "unit_id") %>%
    left_join(sc %>% transmute(unit_id, scenic_n_reg = n_tests, scenic_n_p05 = n_hit, scenic_animal_hit = n_hit > qbinom(0.95, n_tests, 0.05)), by = "unit_id") %>%
    left_join(ma %>% transmute(unit_id, mast_n_tested = n_tests, mast_n_deg = n_hit), by = "unit_id") %>%
    mutate(mast_burden_hit = (mast_n_deg / mast_n_tested) >= quantile(mast_n_deg / mast_n_tested, 0.75, na.rm = TRUE))
  write_csv(U, file.path(P$tbl, "integrated_unit_flags.csv"))
  log_msg("Circuit enrichment:"); print(res %>% dplyr::select(modality, rate_circuit, rate_other, rate_ratio, p_perm, ratio_adj, p_adj))
}

# ════════════════════════════════════════════════════════════════════════════
# Stage 5 · integrated evidence table
# ════════════════════════════════════════════════════════════════════════════
stage_evidence <- function() {
  U <- read_csv(file.path(P$tbl, "integrated_unit_flags.csv"), show_col_types = FALSE)
  regu <- read_csv(file.path(P$tbl, "integrated_regulon_animal_unit.csv"), show_col_types = FALSE)
  gb <- read_csv(file.path(P$tbl, "integrated_gene_bridge.csv"), show_col_types = FALSE)
  het <- read_csv(file.path(P$tbl, "integrated_heterogeneity_summary.csv"), show_col_types = FALSE)
  pbu <- read_csv(file.path(P$tbl, "integrated_pseudobulk_unit_summary.csv"), show_col_types = FALSE)
  top_reg <- regu %>% group_by(unit_id) %>% slice_min(p_perm, n = 3, with_ties = FALSE) %>%
    summarise(top_regulons = paste(glue("{regulon} ({if_else(effect > 0, '↑', '↓')}{signif(p_perm, 1)})"), collapse = ", "), .groups = "drop")
  gb_u <- gb %>% filter(level == "supertype") %>% group_by(unit_id) %>% slice_min(p_state, n = 1, with_ties = FALSE) %>%
    transmute(unit_id, lead_comp = comp, lead_p_state = p_state, shifted_program_vs_mast_z = mast_dirz_shift,
              program_regulon = scenic_top_shift, program_regulon_q = scenic_top_q_shift, program_dspin = dspin_best_shift,
              program_dspin_resp_p = dspin_best_resp_p_shift)
  ev <- U %>% left_join(top_reg, by = "unit_id") %>% left_join(gb_u, by = "unit_id") %>%
    left_join(het %>% dplyr::select(unit_id, n_het_fdr10, excess_het_p01), by = "unit_id") %>%
    left_join(pbu %>% dplyr::select(unit_id, rho_mast_pb), by = "unit_id") %>%
    mutate(n_animal_modalities = rowSums(across(c(pa_state_hit, pb_excess_hit, dspin_hit, scenic_animal_hit), ~ coalesce(.x, FALSE))),
           n_animal_tested = rowSums(across(c(pa_state_hit, pb_excess_hit, dspin_hit, scenic_animal_hit), ~ !is.na(.x)))) %>%
    arrange(desc(n_animal_modalities), desc(circuit), lead_p_state)
  write_csv(ev, file.path(P$tbl, "integrated_evidence_by_unit.csv"))
  leads <- read_csv(file.path(P$pa_tbl, "popAlign_followup_leads.csv"), show_col_types = FALSE)
  lead_ev <- leads %>% dplyr::select(tissue, lineage, level, unit_id, comp, effect, p, convergent) %>%
    left_join(gb %>% dplyr::select(unit_id, comp, shifted_feature, mast_dirz_shift, mast_dirz_control, mast_overlap_shift,
                                    dspin_best_shift, dspin_best_resp_p_shift, dspin_dir_concord_shift,
                                    scenic_top_shift, scenic_top_q_shift, scenic_top_diff_padj_shift, scenic_n_enriched_shift), by = c("unit_id", "comp")) %>%
    left_join(ev %>% dplyr::select(unit_id, n_animal_modalities, n_animal_tested, top_regulons, mast_n_deg, pb_n_p05, n_het_fdr10), by = "unit_id")
  write_csv(lead_ev, file.path(P$tbl, "integrated_evidence_leads.csv"))
  log_msg("Evidence tables written: ", nrow(ev), " units, ", nrow(lead_ev), " leads")
}

# ════════════════════════════════════════════════════════════════════════════
if (stage %in% c("gene", "all")) stage_gene()
if (stage %in% c("cell", "all")) stage_cell()
if (stage %in% c("animal", "all")) stage_animal()
if (stage %in% c("circuit", "all")) stage_circuit()
if (stage %in% c("evidence", "all")) stage_evidence()
if (stage %in% c("figures", "all")) source("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/notebooks/R_scripts/plot_integrated_analysis_sn.R")
log_msg("Stage '", stage, "' finished")
