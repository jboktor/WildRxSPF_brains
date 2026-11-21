#!/usr/bin/env Rscript
# Run VR-1…VR-4 analyses (animal-level; no MAST/D-SPIN/SCENIC re-fit).
# Usage: Rscript notebooks/R_scripts/run_vr_analyses.R [VR1|VR2|VR3|VR4|all]

Sys.setenv(
  RETICULATE_PYTHON = "/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python"
)

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(Matrix)
  library(jsonlite)
  library(data.table)
  library(ggplot2)
})

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
setwd(wkdir)
source(file.path(wkdir, "notebooks/R_scripts/vr_relay_helpers.R"))
paths <- vr_paths(wkdir)

args <- commandArgs(trailingOnly = TRUE)
targets <- if (!length(args) || args[[1]] == "all") {
  c("VR4_step1", "VR1", "VR2", "VR4", "VR3")
} else {
  args
}

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), paste0(...), "\n")
  flush.console()
}

# Outputs follow repo layout: tables + markdown result notes → tables/VR/<id>,
# figures → figures/VR/<id>.
ensure_dirs <- function(...) {
  for (d in c(...)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ── gene scoring via Python subprocess (reticulate conda probe is broken here) ─
PY <- "/resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/python"
SCORE_PY <- file.path(wkdir, "notebooks/python_scripts/score_genes_h5ad.py")

score_unit_scanpy <- function(sc = NULL, supertype_label, genes_pos, genes_neg = character()) {
  h5 <- file.path(paths$dspin, paste0("AMY__", supertype_label), "raw_counts.h5ad")
  meta_f <- file.path(paths$dspin, paste0("AMY__", supertype_label), "cells_meta.csv")
  if (!file.exists(h5) || !file.exists(meta_f)) return(NULL)
  tmp <- tempfile("vr_score_")
  pos_f <- paste0(tmp, "_pos.txt")
  neg_f <- paste0(tmp, "_neg.txt")
  out_f <- paste0(tmp, "_out.csv")
  writeLines(genes_pos, pos_f)
  if (length(genes_neg)) writeLines(genes_neg, neg_f) else writeLines(character(), neg_f)
  rc <- system2(PY, c(SCORE_PY, h5, meta_f, pos_f, if (length(genes_neg)) neg_f else "-", out_f),
                stdout = TRUE, stderr = TRUE)
  if (!file.exists(out_f)) {
    warning("score failed for ", supertype_label, ": ", paste(rc, collapse = " | "))
    return(NULL)
  }
  scores <- readr::read_csv(out_f, show_col_types = FALSE)
  meta <- readr::read_csv(meta_f, show_col_types = FALSE)
  unlink(c(pos_f, neg_f, out_f))
  n <- min(nrow(scores), nrow(meta))
  tibble(
    cell_id = as.character(scores$cell_id[seq_len(n)]),
    sample_id = as.character(scores$sample_id[seq_len(n)]),
    animal = animal_id_from_sample(scores$sample_id[seq_len(n)]),
    group = group_from_sample(scores$sample_id[seq_len(n)]),
    microbiome = meta$microbiome[seq_len(n)],
    score_pos_only = as.numeric(scores$score_pos_only[seq_len(n)]),
    score_signed = as.numeric(scores$score_signed[seq_len(n)]),
    supertype_label = supertype_label,
    median_umi_proxy = NA_real_
  )
}

get_scanpy <- function() NULL  # unused; kept so call sites stay valid

glut_labels <- c(
  "CS20230722_SUPT_0039", "CS20230722_SUPT_0040", "CS20230722_SUPT_0041",
  "CS20230722_SUPT_0042", "CS20230722_SUPT_0055", "CS20230722_SUPT_0056",
  "CS20230722_SUPT_0057", "CS20230722_SUPT_0058", "CS20230722_SUPT_0061",
  "CS20230722_SUPT_0507"
)

# ══════════════════════════════════════════════════════════════════════════════
run_vr4_step1 <- function() {
  log_msg("VR4 step 1 — Zbtb16 resolution")
  out_dir <- file.path(wkdir, "tables/VR/VR4")
  fig_dir <- file.path(wkdir, "figures/VR/VR4")
  ensure_dirs(out_dir, fig_dir)
  amy_classes <- c(
    "09" = "snRNASeq_AMY_CS20230722_CLAS_09",
    "11" = "snRNASeq_AMY_CS20230722_CLAS_11",
    "13" = "snRNASeq_AMY_CS20230722_CLAS_13"
  )
  parse_regulon_tfs <- function(regulons_csv) {
    raw <- readLines(regulons_csv, warn = FALSE)
    tfs <- character()
    for (ln in raw) {
      fields <- strsplit(ln, ",", fixed = TRUE)[[1]]
      if (!length(fields)) next
      tf <- fields[[1]]
      if (grepl("^[A-Za-z][A-Za-z0-9]+$", tf) && !tf %in% c("TF", "MotifID")) {
        tfs <- c(tfs, tf)
      }
    }
    unique(tfs)
  }
  auc_has_tf <- function(auc_csv, tf = "Zbtb16") {
    hdr <- strsplit(readLines(auc_csv, n = 1), ",", fixed = TRUE)[[1]]
    any(grepl(paste0("^", tf, "\\(\\+\\)$"), hdr))
  }
  cases <- bind_rows(lapply(names(amy_classes), function(cl) {
    d <- file.path(paths$scenic_by_class, amy_classes[[cl]])
    reg <- file.path(d, "regulons.csv")
    auc <- file.path(d, "auc_mtx.csv")
    diff <- file.path(d, "diff_regulons_WildR_vs_SPF.csv")
    tfs <- if (file.exists(reg)) parse_regulon_tfs(reg) else character()
    called <- "Zbtb16" %in% tfs
    in_auc <- file.exists(auc) && auc_has_tf(auc, "Zbtb16")
    sig <- FALSE
    if (file.exists(diff)) {
      dr <- readr::read_csv(diff, show_col_types = FALSE)
      rcol <- intersect(names(dr), c("regulon", "Regulon", "tf"))[1]
      pcol <- intersect(names(dr), c("padj", "p_adj", "fdr"))[1]
      if (!is.na(rcol)) {
        hit <- grepl("Zbtb16", dr[[rcol]])
        if (any(hit) && !is.na(pcol)) sig <- any(dr[[pcol]][hit] < 0.05, na.rm = TRUE)
      }
    }
    tibble(class = cl, called = called, in_auc = in_auc, class_significant = sig)
  }))
  overall <- if (!any(cases$called) && !any(cases$in_auc)) {
    "a"
  } else if (any(cases$called | cases$in_auc) && !any(cases$class_significant)) {
    "b"
  } else {
    "c"
  }
  case_sentence <- switch(
    overall,
    a = "**Case (a): Regulon never called** in AMY classes 09 CNU-LGE GABA, 11 CNU-HYa GABA, or 13 CNU-HYa Glut — motif enrichment did not yield a `Zbtb16(+)` regulon, so SCENIC absence carries **no information** about activity; `Zbtb16` reverts to a plain DE hit.",
    b = "**Case (b): Regulon called but not significant at class level** — proceed to supertype-resolution AUCell test in MeA/COAp.",
    c = "**Case (c): Regulon called and significant at class level** — prior summary was resolution-limited; still run supertype test."
  )
  lines <- c(
    case_sentence, "", "## Evidence by class", "",
    paste(capture.output(print(cases)), collapse = "\n"), "",
    "## Implication",
    if (overall == "a") {
      "- Do **not** report SCENIC absence as evidence against the TF hypothesis. Elevate `Eps8l1` as the primary thread-1 signal pending other support."
    } else {
      "- See `VR4_zbtb16_supertype_test.csv` for per-animal MeA/COAp AUCell contrast."
    },
    "",
    sprintf("- Overall case: **%s**", overall),
    ""
  )
  write_md(file.path(out_dir, "VR4_zbtb16_resolution.md"), lines)
  write_csv(tibble(note = "Zbtb16 regulon not called; supertype AUCell test not applicable"),
            file.path(out_dir, "VR4_zbtb16_supertype_test.csv"))
  log_msg("VR4 step 1 done: case ", overall)
  invisible(overall)
}

# ══════════════════════════════════════════════════════════════════════════════
run_vr1 <- function() {
  log_msg("VR1 — input/output coupling")
  out_dir <- file.path(wkdir, "tables/VR/VR1")
  fig_dir <- file.path(wkdir, "figures/VR/VR1")
  ensure_dirs(out_dir, fig_dir)
  tbl_dir <- out_dir

  mast <- load_mast_hits(paths = paths) %>%
    filter(tissue == "AMY", !excluded_gene) %>%
    filter(chain_position %in% c("1", "2", "3", "4")) %>%
    filter(grepl("Glut|IT |COAp|MEA Slc17a7|LA-BLA|MEA-COA", supertype_name))

  # Catalogue of 243×27 not on disk. Approximate with multi-hit WildR-UP genes in
  # glut chain positions (≥2 hits) plus expected core — avoid scoring 1000+ genes.
  pos_counts <- mast %>% filter(direction == "WildR-UP") %>% count(gene, name = "n_hits")
  expected_core <- c("Eps8l1", "Zbtb16", "Hdac9", "Trpm3", "Rasa1")
  sig_pos <- unique(c(
    intersect(expected_core, pos_counts$gene),
    pos_counts %>% filter(n_hits >= 2) %>% pull(gene)
  ))
  sig_neg <- intersect("Itpkb", mast %>% filter(direction == "SPF-UP") %>% pull(gene) %>% unique())

  log_msg("Signature pos n=", length(sig_pos), " neg=", paste(sig_neg, collapse = ","))
  write_json(list(genes_pos = sig_pos, genes_neg = sig_neg, n_pos = length(sig_pos), n_neg = length(sig_neg)),
             file.path(tbl_dir, "VR1_signature_genes.json"), pretty = TRUE, auto_unbox = TRUE)

  log_msg("Scoring ", length(glut_labels), " glut units…")
  cell_scores <- bind_rows(lapply(glut_labels, function(lab) {
    log_msg("  ", lab)
    score_unit_scanpy(NULL, lab, sig_pos, sig_neg)
  }))
  write_csv(cell_scores, file.path(tbl_dir, "VR1_cell_scores.csv"))

  input_eq <- cell_scores %>%
    group_by(animal, group, supertype_label) %>%
    summarise(score_signed = mean(score_signed, na.rm = TRUE),
              score_pos_only = mean(score_pos_only, na.rm = TRUE),
              n_cells = n(), .groups = "drop") %>%
    group_by(animal, group) %>%
    summarise(
      input_score = mean(score_signed, na.rm = TRUE),
      input_score_pos_only = mean(score_pos_only, na.rm = TRUE),
      n_supertypes = n(),
      n_cells_total = sum(n_cells),
      .groups = "drop"
    )
  input_w <- cell_scores %>%
    group_by(animal, group) %>%
    summarise(input_score_cellwt = mean(score_signed, na.rm = TRUE), .groups = "drop")
  input_animal <- input_eq %>% left_join(input_w, by = c("animal", "group"))

  gaba <- gaba_output_programs()
  h_long <- bind_rows(lapply(seq_len(nrow(gaba)), function(i) {
    read_dspin_unit_h("AMY", gaba$supertype_label[i], gaba$program_id[i],
                      which = "absolute", paths = paths) %>%
      mutate(key = paste0(gaba$supertype_label[i], "_", gaba$program_id[i]))
  }))
  output_wide <- h_long %>%
    select(animal, group, key, h) %>%
    tidyr::pivot_wider(names_from = key, values_from = h) %>%
    rename(
      output_0338_P10 = CS20230722_SUPT_0338_P10,
      output_0355_P5 = CS20230722_SUPT_0355_P5,
      output_0359_P6 = CS20230722_SUPT_0359_P6
    ) %>%
    mutate(
      output_score_mean = rowMeans(across(starts_with("output_")), na.rm = TRUE),
      gaba_sign_agree = (
        sign(output_0338_P10) == sign(output_0355_P5) &
          sign(output_0355_P5) == sign(output_0359_P6)
      )
    )

  per_animal <- input_animal %>%
    inner_join(output_wide, by = c("animal", "group")) %>%
    mutate(
      z_input = as.numeric(scale(input_score)),
      z_output = as.numeric(scale(output_score_mean)),
      D = z_input - z_output
    )
  scores_csv <- file.path(tbl_dir, "VR1_per_animal_scores.csv")
  write_csv(per_animal, scores_csv)

  sep <- complete_separation(per_animal$D, per_animal$group, higher = "WildR")
  cd <- cliffs_delta(per_animal$D[per_animal$group == "WildR"],
                     per_animal$D[per_animal$group == "SPF"], n_boot = 10000L)
  wx <- suppressWarnings(wilcox.test(D ~ group, data = per_animal, exact = FALSE))
  sp <- spearman_ci(per_animal$input_score, per_animal$output_score_mean, n_boot = 5000L)

  # covariates (median UMI may be unavailable from h5ad path)
  cov <- cell_scores %>%
    group_by(animal) %>%
    summarise(med_umi = median(median_umi_proxy, na.rm = TRUE),
              n_nuclei = n(), .groups = "drop") %>%
    right_join(per_animal, by = "animal")
  safe_cor <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3) return(NA_real_)
    suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
  }
  rho_in_umi <- safe_cor(cov$input_score, cov$med_umi)
  rho_out_umi <- safe_cor(cov$output_score_mean, cov$med_umi)
  rho_in_n <- safe_cor(cov$input_score, cov$n_nuclei)
  rho_out_n <- safe_cor(cov$output_score_mean, cov$n_nuclei)

  branch <- if (isTRUE(sep$complete) || (sep$n_exceed_all >= 5 && isTRUE(sp$rho < 0))) {
    "**Decision:** complete/near-complete separation on `D` plus negative input–output correlation → coupled circuit disinhibition licensed. Proceed to VR-2."
  } else if (sep$n_exceed_all >= 5 && !(isTRUE(sp$rho < 0))) {
    "**Decision:** separation on `D` without negative correlation → parallel modulation; 'net disinhibition' OK, 'one coupled circuit' not."
  } else {
    "**Decision:** no clear separation on `D` → retire the coupled-circuit sentence (valuable negative)."
  }

  lines <- c(
    "# VR-1 results", "",
    "## Headline (decision rule)", "",
    sprintf("- Per-animal separation on `D`: **%d / %d** WildR exceed all SPF.", sep$n_exceed_all, sep$n_hi),
    sprintf("- Cliff's δ (WildR vs SPF on `D`): **%.3f** (95%% bootstrap CI %.3f – %.3f).", cd$delta, cd$ci_low, cd$ci_high),
    sprintf("- Spearman(input, output) across 12 animals: **ρ = %.3f** (CI %.3f – %.3f).", sp$rho, sp$ci_low, sp$ci_high),
    sprintf("- Wilcoxon rank-sum p (not headline): %.4g.", wx$p.value),
    "", "## Branch", "", branch, "",
    "## Technical covariates", "",
    sprintf("- Spearman(input, median UMI)=%.3f; Spearman(output, median UMI)=%.3f", rho_in_umi, rho_out_umi),
    sprintf("- Spearman(input, n_nuclei)=%.3f; Spearman(output, n_nuclei)=%.3f", rho_in_n, rho_out_n),
    "", "## Signature", "",
    sprintf("- Positive genes (n=%d): %s", length(sig_pos), paste(head(sig_pos, 40), collapse = ", "),
            if (length(sig_pos) > 40) "…" else ""),
    sprintf("- Signed-negative genes (n=%d): %s", length(sig_neg), paste(sig_neg, collapse = ", ")),
    "- Scoring: `scanpy.tl.score_genes` on per-unit `raw_counts.h5ad`.",
    "- Aggregation: equal weight across glutamatergic AMY relay supertypes.",
    "- Comparison: **between-animal** WildR vs SPF on within-animal `D`.",
    ""
  )
  # fix signature line if paste overflow
  lines[which(grepl("Positive genes", lines))] <- sprintf(
    "- Positive genes (n=%d): %s", length(sig_pos),
    paste(sig_pos[seq_len(min(50, length(sig_pos)))], collapse = ", ")
  )
  write_md(file.path(out_dir, "VR1_results.md"), lines)

  # Sensitivity
  sens_rows <- list()
  if ("Zbtb16" %in% sig_pos) {
    log_msg("Sensitivity: exclude Zbtb16")
    pos2 <- setdiff(sig_pos, "Zbtb16")
    cs2 <- bind_rows(lapply(glut_labels, function(lab) score_unit_scanpy(NULL, lab, pos2, sig_neg)))
    in2 <- cs2 %>%
      group_by(animal, group, supertype_label) %>%
      summarise(s = mean(score_signed, na.rm = TRUE), .groups = "drop") %>%
      group_by(animal, group) %>%
      summarise(input_score = mean(s, na.rm = TRUE), .groups = "drop")
    pa2 <- in2 %>% inner_join(output_wide, by = c("animal", "group")) %>%
      mutate(D = as.numeric(scale(input_score)) - as.numeric(scale(output_score_mean)))
    cd2 <- cliffs_delta(pa2$D[pa2$group == "WildR"], pa2$D[pa2$group == "SPF"])
    sens_rows[[length(sens_rows) + 1]] <- tibble(
      setting = "exclude_Zbtb16", cliff_delta = cd2$delta, ci_low = cd2$ci_low, ci_high = cd2$ci_high
    )
  }
  for (col in c("output_0338_P10", "output_0355_P5", "output_0359_P6")) {
    pa3 <- per_animal %>%
      mutate(D = as.numeric(scale(input_score)) - as.numeric(scale(.data[[col]])))
    cd3 <- cliffs_delta(pa3$D[pa3$group == "WildR"], pa3$D[pa3$group == "SPF"])
    sens_rows[[length(sens_rows) + 1]] <- tibble(
      setting = paste0("output_only_", col), cliff_delta = cd3$delta,
      ci_low = cd3$ci_low, ci_high = cd3$ci_high
    )
  }
  write_csv(bind_rows(sens_rows), file.path(tbl_dir, "VR1_sensitivity.csv"))

  p <- ggplot(per_animal, aes(input_score, output_score_mean, colour = group)) +
    geom_hline(yintercept = mean(per_animal$output_score_mean), linetype = 2, colour = "grey70") +
    geom_vline(xintercept = mean(per_animal$input_score), linetype = 2, colour = "grey70") +
    geom_point(size = 3.5, alpha = 0.9) +
    scale_colour_manual(values = VR_COLORS) +
    labs(title = "Input-limb score vs MeA–BST output h",
         x = "Input signature (equal-weight animal mean)",
         y = "Output score (mean h of 0338-P10, 0355-P5, 0359-P6)",
         colour = NULL) +
    vr_brand_ggplot()
  scatter_png <- file.path(fig_dir, "VR1_scatter.png")
  ggsave(scatter_png, p, width = 6.5, height = 5, dpi = 200)

  log_msg("VR1 done. Separation=", sep$n_exceed_all, "/", sep$n_hi,
          " Cliff δ=", round(cd$delta, 3), " Spearman=", round(sp$rho, 3))
  invisible(list(sep = sep, cd = cd, sp = sp, branch = branch))
}

# ══════════════════════════════════════════════════════════════════════════════
run_vr2 <- function() {
  log_msg("VR2 — MeA-BST program identity")
  out_dir <- file.path(wkdir, "tables/VR/VR2")
  fig_dir <- file.path(wkdir, "figures/VR/VR2")
  ensure_dirs(out_dir, fig_dir)
  tbl_dir <- out_dir
  gaba <- gaba_output_programs()

  read_prog <- function(label, pid) {
    f <- file.path(paths$dspin, paste0("AMY__", label), "programs_genes.csv")
    readr::read_csv(f, show_col_types = FALSE) %>%
      filter(program_id == pid) %>%
      mutate(gene = strip_grcm39(gene)) %>%
      filter(!gene %in% VR_EXCLUDED_GENES, nzchar(gene), !is.na(gene)) %>%
      pull(gene) %>% unique()
  }
  sets <- list(
    `0338_P10` = read_prog("CS20230722_SUPT_0338", "P10"),
    `0355_P5`  = read_prog("CS20230722_SUPT_0355", "P5"),
    `0359_P6`  = read_prog("CS20230722_SUPT_0359", "P6")
  )
  keys <- names(sets)
  jac <- matrix(NA_real_, 3, 3, dimnames = list(keys, keys))
  for (i in 1:3) for (j in 1:3) jac[i, j] <- jaccard(sets[[i]], sets[[j]])
  three_way <- Reduce(intersect, sets)

  pool <- unique(unlist(lapply(gaba$supertype_label, function(lab) {
    f <- file.path(paths$dspin, paste0("AMY__", lab), "programs_genes.csv")
    strip_grcm39(readr::read_csv(f, show_col_types = FALSE)$gene)
  })))
  pool <- setdiff(pool[!is.na(pool) & nzchar(pool)], VR_EXCLUDED_GENES)
  set.seed(1)
  null_j <- replicate(1000L, {
    a <- sample(pool, length(sets[[1]]))
    b <- sample(pool, length(sets[[2]]))
    c <- sample(pool, length(sets[[3]]))
    c(j12 = jaccard(a, b), j13 = jaccard(a, c), j23 = jaccard(b, c))
  })

  lines <- c(
    "# VR-2 program overlap", "",
    "## Gene-set sizes",
    sprintf("- 0338 P10: %d genes", length(sets[[1]])),
    sprintf("- 0355 P5: %d genes", length(sets[[2]])),
    sprintf("- 0359 P6: %d genes", length(sets[[3]])),
    sprintf("- Three-way intersection (n=%d): %s", length(three_way),
            if (length(three_way)) paste(three_way, collapse = ", ") else "(empty)"),
    "", "## Jaccard matrix", "", "```",
    paste(capture.output(print(round(jac, 3))), collapse = "\n"), "```", "",
    "## Null (1,000 draws)",
    sprintf("- Observed J(0338,0355)=%.3f; null mean=%.3f (2.5–97.5%%: %.3f–%.3f)",
            jac[1, 2], mean(null_j["j12", ]),
            quantile(null_j["j12", ], 0.025), quantile(null_j["j12", ], 0.975)),
    sprintf("- Observed J(0338,0359)=%.3f; null mean=%.3f", jac[1, 3], mean(null_j["j13", ])),
    sprintf("- Observed J(0355,0359)=%.3f; null mean=%.3f", jac[2, 3], mean(null_j["j23", ])),
    ""
  )
  write_md(file.path(out_dir, "VR2_program_overlap.md"), lines)
  write_csv(as.data.frame(jac) %>% tibble::rownames_to_column("program"),
            file.path(tbl_dir, "VR2_jaccard.csv"))

  # Identity scores (no Lhx6/Sox6)
  interneuron <- c("Pvalb", "Sst", "Vip", "Npy", "Lamp5", "Htr3a", "Reln", "Cnr1", "Calb2")
  projection <- c("Foxp2", "Bcl11b", "Isl1", "Six3", "Pdyn", "Tac1", "Penk", "Nts", "Cartpt")
  score_markers2 <- function(label, panel, name) {
    df <- score_unit_scanpy(NULL, label, panel, character())
    if (is.null(df)) return(NULL)
    df %>% mutate(panel = name, score = score_pos_only, genes_used = length(panel))
  }
  id_long <- bind_rows(
    lapply(gaba$supertype_label, score_markers2, panel = interneuron, name = "interneuron"),
    lapply(gaba$supertype_label, score_markers2, panel = projection, name = "projection")
  )
  identity <- id_long %>%
    group_by(supertype_label, animal, group, panel) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), n_cells = n(),
              genes_used = first(genes_used), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = panel, values_from = c(mean_score, n_cells, genes_used))
  write_csv(identity, file.path(tbl_dir, "VR2_identity_scores.csv"))

  # Enrichment
  union_genes <- unique(unlist(sets))
  cats <- list(
    gaba_identity = c("Gad1", "Gad2", "Slc32a1"),
    gaba_a_receptor = grep("^Gabra|^Gabrb|^Gabrg|^Gabrr|^Gabre|^Gabrd", pool, value = TRUE),
    synaptic_organisers = grep("^Nlgn|^Nrxn|^Cntn|^Lrrtm", pool, value = TRUE),
    neuropeptides = intersect(c("Pdyn", "Tac1", "Penk", "Nts", "Cartpt", "Npy", "Sst", "Vip", "Crh"), pool),
    axon_guidance = grep("^Sema|^Robo|^Epha|^Ephb|^Dcc", pool, value = TRUE)
  )
  enrich <- bind_rows(lapply(names(cats), function(nm) {
    cat_genes <- intersect(cats[[nm]], pool)
    k <- length(intersect(union_genes, cat_genes))
    n <- length(union_genes); K <- length(cat_genes); N <- length(pool)
    p <- if (K == 0) NA_real_ else phyper(k - 1, K, N - K, n, lower.tail = FALSE)
    tibble(category = nm, n_in_union = k, category_size = K, union_size = n,
           background_size = N,
           background = "genes in any program of 0338/0355/0359 units",
           hypergeom_p = p,
           enrichment = (k / max(n, 1)) / max(K / max(N, 1), .Machine$double.eps))
  }))
  write_csv(enrich, file.path(tbl_dir, "VR2_enrichment.csv"))

  # MAST bridge
  mast27 <- load_mast_hits(paths = paths) %>%
    filter(tissue == "AMY", !excluded_gene) %>% distinct(gene) %>% pull(gene)
  focus <- unique(c("Hap1", "Cntn6", "Itpkb", "Eps8l1", mast27))
  bridge <- bind_rows(lapply(names(sets), function(nm) {
    tibble(program = nm, gene = focus, in_program = focus %in% sets[[nm]])
  }))
  hits_bridge <- bridge %>% filter(in_program)
  blines <- c(
    "# VR-2 MAST ↔ D-SPIN gene-level bridge", "",
    if (nrow(hits_bridge)) {
      c("**Gene-level agreements between D-SPIN programs and the MAST catalogue:**", "",
        paste(knitr::kable(hits_bridge), collapse = "\n"))
    } else {
      "- No MAST-catalogue gene falls inside the three program gene sets among the focus list / AMY hits."
    },
    ""
  )
  write_md(file.path(out_dir, "VR2_mast_bridge.md"), blines)
  write_csv(bridge, file.path(tbl_dir, "VR2_mast_bridge.csv"))

  # Spatial
  log_msg("VR2 spatial map")
  mmc <- tryCatch(
    readr::read_csv(paths$mapmycells, show_col_types = FALSE, comment = "#"),
    error = function(e) NULL
  )
  meta_sp <- readRDS(paths$spatial_meta) %>% as_tibble(rownames = "cell_id") %>%
    mutate(cell_id = as.character(cell_id))
  if (!is.null(mmc)) {
    id_col <- names(mmc)[1]
    name_col <- grep("supertype_name", names(mmc), value = TRUE)[1]
    if (!is.na(name_col)) {
      gaba_cells <- mmc %>%
        rename(cell_id = !!id_col) %>%
        mutate(cell_id = as.character(cell_id)) %>%
        filter(grepl("0338 |0355 |0359 |MEA-BST", .data[[name_col]])) %>%
        transmute(cell_id, supertype_name = .data[[name_col]])
    } else {
      gaba_cells <- tibble(cell_id = character(), supertype_name = character())
    }
  } else {
    gaba_cells <- tibble(cell_id = character(), supertype_name = character())
  }
  plot_df <- meta_sp %>%
    left_join(gaba_cells, by = "cell_id") %>%
    mutate(is_gaba_supertype = !is.na(supertype_name))
  frac <- plot_df %>%
    filter(is_gaba_supertype) %>%
    count(roi, bregma, name = "n") %>%
    group_by(roi) %>%
    mutate(frac = n / sum(n), compartment = "label_present_CCF_join_pending") %>%
    ungroup()
  write_csv(frac, file.path(tbl_dir, "VR2_spatial_fractions.csv"))

  # Panel coverage estimate: program genes vs a gene list if available
  coverage_note <- sprintf(
    "MapMyCells-labelled MeA-BST cells: %d. CCF compartment join pending; treat localisation as qualitative.",
    sum(plot_df$is_gaba_supertype)
  )
  p <- ggplot() +
    geom_point(data = plot_df %>% slice_sample(n = min(30000L, nrow(plot_df))),
               aes(center_x, center_y), colour = "grey85", size = 0.05) +
    geom_point(data = plot_df %>% filter(is_gaba_supertype),
               aes(center_x, center_y, colour = supertype_name), size = 0.35) +
    coord_equal() +
    labs(title = "MeA–BST GABA-assigned cells on seqFISH sections",
         subtitle = coverage_note, colour = NULL) +
    vr_brand_ggplot() + theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "VR2_spatial_distribution.png"), p, width = 10, height = 7, dpi = 200)

  # Decision note
  mean_j <- mean(c(jac[1, 2], jac[1, 3], jac[2, 3]))
  decision <- if (mean_j > mean(null_j) + 2 * sd(as.numeric(null_j))) {
    "Shared membership above null → treat as coordinated programs pending identity/spatial."
  } else {
    "Jaccard not clearly above null → three separate programs may co-occur; MeA–BST framing weakens."
  }
  write_md(file.path(out_dir, "VR2_program_overlap.md"),
           c(readLines(file.path(out_dir, "VR2_program_overlap.md"), warn = FALSE),
             "", "## Decision inputs", paste0("- ", decision),
             paste0("- ", coverage_note), ""))
  log_msg("VR2 done. mean Jaccard=", round(mean_j, 3))
}

# ══════════════════════════════════════════════════════════════════════════════
run_vr4 <- function() {
  log_msg("VR4 steps 2–6 — relay regulon scan / limbs / RA")
  out_dir <- file.path(wkdir, "tables/VR/VR4")
  fig_dir <- file.path(wkdir, "figures/VR/VR4")
  ensure_dirs(out_dir, fig_dir)
  tbl_dir <- out_dir
  amy_classes <- c(
    "09" = "snRNASeq_AMY_CS20230722_CLAS_09",
    "11" = "snRNASeq_AMY_CS20230722_CLAS_11",
    "13" = "snRNASeq_AMY_CS20230722_CLAS_13"
  )
  # Prefer full AMY dspin cells_meta (joins to class AUC); scenic/cell_meta.csv is a pilot subsample
  log_msg("Building AMY cell meta from dspin units…")
  amy_units <- list.files(paths$dspin, pattern = "^AMY__", full.names = TRUE)
  cell_meta <- bind_rows(lapply(amy_units, function(u) {
    f <- file.path(u, "cells_meta.csv")
    if (!file.exists(f)) return(NULL)
    readr::read_csv(f, show_col_types = FALSE)
  })) %>%
    mutate(
      cell_id = as.character(cell_id),
      animal = animal_id_from_sample(sample_id),
      group = group_from_sample(sample_id),
      chain_position = assign_chain_position(supertype_name),
      tissue = "AMY"
    )
  log_msg("AMY meta cells: ", nrow(cell_meta))

  clas01_auc <- file.path(paths$scenic_by_class, "snRNASeq_AMY_CS20230722_CLAS_01/auc_mtx.csv")
  clas01_ok <- file.exists(clas01_auc) && file.info(clas01_auc)$size > 0
  if (clas01_ok) {
    amy_classes <- c("01" = "snRNASeq_AMY_CS20230722_CLAS_01", amy_classes)
    log_msg("CLAS_01 AUC available — full positions 1–4 gradient enabled")
  } else {
    log_msg("CLAS_01 AUC missing — positions 1/4 gradient blocked; two-limb (3 vs 3G) still runs")
  }

  load_class_auc_long <- function(class_dir, class_id) {
    auc_f <- file.path(paths$scenic_by_class, class_dir, "auc_mtx.csv")
    if (!file.exists(auc_f) || !file.info(auc_f)$size) return(NULL)
    log_msg("  loading AUC ", class_id)
    auc <- data.table::fread(auc_f)
    cell <- as.character(auc[[1]])
    mat <- as.matrix(auc[, -1])
    mat_z <- scale(mat)
    colnames(mat_z) <- colnames(mat)
    as_tibble(mat_z) %>%
      mutate(cell_id = cell, class = class_id) %>%
      tidyr::pivot_longer(-c(cell_id, class), names_to = "regulon", values_to = "auc_z")
  }
  auc_long <- bind_rows(lapply(names(amy_classes), function(cl) {
    load_class_auc_long(amy_classes[[cl]], cl)
  }))
  joined <- auc_long %>%
    inner_join(cell_meta %>% select(cell_id, animal, group, tissue, supertype_name, chain_position),
               by = "cell_id") %>%
    filter(tissue == "AMY", chain_position %in% c("1", "2", "3", "4", "3G"))
  log_msg("Joined relay AUC rows: ", nrow(joined),
          " positions: ", paste(sort(unique(joined$chain_position)), collapse = ","))

  covar <- joined %>%
    distinct(cell_id, animal, chain_position) %>%
    count(animal, chain_position, name = "n_cells")

  anim_pos <- joined %>%
    filter(chain_position %in% c("1", "2", "3", "4")) %>%
    group_by(regulon, animal, group, chain_position) %>%
    summarise(auc = mean(auc_z, na.rm = TRUE), .groups = "drop") %>%
    mutate(pos = as.integer(chain_position)) %>%
    left_join(covar, by = c("animal", "chain_position"))

  grad <- anim_pos %>%
    group_by(regulon, pos, chain_position) %>%
    summarise(
      delta = mean(auc[group == "WildR"], na.rm = TRUE) - mean(auc[group == "SPF"], na.rm = TRUE),
      med_n = median(n_cells, na.rm = TRUE),
      .groups = "drop"
    )

  fit_slope <- function(df, adjusted = FALSE) {
    if (dplyr::n_distinct(df$pos) < 2 || nrow(df) < 2) {
      return(tibble(slope = NA_real_, ci_low = NA_real_, ci_high = NA_real_, n_pos = dplyr::n_distinct(df$pos)))
    }
    fit <- tryCatch({
      if (adjusted) lm(delta ~ pos + med_n, data = df) else lm(delta ~ pos, data = df)
    }, error = function(e) NULL)
    if (is.null(fit) || is.na(coef(fit)[["pos"]])) {
      return(tibble(slope = NA_real_, ci_low = NA_real_, ci_high = NA_real_, n_pos = dplyr::n_distinct(df$pos)))
    }
    cf <- suppressWarnings(confint(fit))
    tibble(slope = unname(coef(fit)[["pos"]]),
           ci_low = unname(cf["pos", 1]), ci_high = unname(cf["pos", 2]),
           n_pos = dplyr::n_distinct(df$pos))
  }

  # Primary ranking when gradient unavailable: |Δ| at MeA glut position 3
  pos3_delta <- anim_pos %>%
    filter(pos == 3L) %>%
    group_by(regulon) %>%
    summarise(
      delta_pos3 = mean(auc[group == "WildR"], na.rm = TRUE) - mean(auc[group == "SPF"], na.rm = TRUE),
      .groups = "drop"
    )

  scan <- grad %>%
    group_by(regulon) %>%
    group_modify(~ {
      u <- fit_slope(.x, FALSE); a <- fit_slope(.x, TRUE)
      tibble(slope = u$slope, ci_low = u$ci_low, ci_high = u$ci_high, n_pos = u$n_pos,
             slope_adj = a$slope, ci_low_adj = a$ci_low, ci_high_adj = a$ci_high)
    }) %>%
    ungroup() %>%
    left_join(pos3_delta, by = "regulon") %>%
    mutate(
      rank_metric = if_else(is.finite(slope_adj), abs(slope_adj), abs(delta_pos3)),
      abs_adj = abs(slope_adj),
      rank_adj = dense_rank(desc(rank_metric)),
      note = if_else(is.finite(slope_adj), "gradient_slope", "ranked_by_pos3_delta_CLAS01_missing")
    ) %>%
    arrange(rank_adj)
  write_csv(scan, file.path(tbl_dir, "VR4_relay_regulon_scan.csv"))
  write_csv(scan %>% filter(grepl("Mllt10|Ube2k", regulon)) %>%
              select(regulon, rank_adj, slope_adj, delta_pos3, ci_low_adj, ci_high_adj, note),
            file.path(tbl_dir, "VR4_recurrent_tf_ranks.csv"))
  writeLines(c(
    "# VR4 scan notes",
    sprintf("- CLAS_01 AUC present: %s", clas01_ok),
    "- Positions 1 and 4 live almost entirely in CLAS_01 (IT-ET Glut).",
    "- Until CLAS_01 finishes, monotone 1–4 gradient slopes are NA; ranking falls back to |WildR−SPF| at position 3 (MeA glut).",
    "- AUCell z-scored within class before any cross-position comparison.",
    ""
  ), file.path(out_dir, "VR4_scan_notes.md"))

  # Two-limb: MeA glut (pos 3) vs MeA-BST GABA (3G), TF names shared across classes
  top20 <- scan %>% filter(rank_adj <= 20) %>% pull(regulon)
  limb_all <- joined %>%
    mutate(limb = case_when(
      chain_position %in% c("2", "3") ~ "input_glut",
      chain_position == "3G" ~ "output_gaba",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(limb)) %>%
    group_by(regulon, animal, group, limb) %>%
    summarise(auc = mean(auc_z, na.rm = TRUE), .groups = "drop") %>%
    group_by(regulon, limb) %>%
    summarise(delta = mean(auc[group == "WildR"], na.rm = TRUE) -
                mean(auc[group == "SPF"], na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = limb, values_from = delta)
  if (!"input_glut" %in% names(limb_all)) limb_all$input_glut <- NA_real_
  if (!"output_gaba" %in% names(limb_all)) limb_all$output_gaba <- NA_real_
  limb <- limb_all %>%
    filter(regulon %in% top20 | (!is.na(input_glut) & !is.na(output_gaba))) %>%
    mutate(opposite_direction = !is.na(input_glut) & !is.na(output_gaba) &
             sign(input_glut) != sign(output_gaba) &
             abs(input_glut) > 0 & abs(output_gaba) > 0,
           in_top20 = regulon %in% top20) %>%
    arrange(desc(opposite_direction), desc(in_top20), regulon)
  write_csv(limb, file.path(tbl_dir, "VR4_two_limb_contrast.csv"))
  log_msg("Opposite-direction regulons: ", sum(limb$opposite_direction, na.rm = TRUE))

  # RA module
  ra_genes <- c("Rarb", "Crabp1", "Crabp2", "Dhrs3", "Rbp1", "Cyp26b1")
  relay_labels <- c(glut_labels, gaba_output_programs()$supertype_label)
  ra_cell <- bind_rows(lapply(unique(relay_labels), function(lab) {
    log_msg("  RA score ", lab)
    score_unit_scanpy(NULL, lab, ra_genes, character())
  }))
  ra_animal <- ra_cell %>%
    group_by(animal, group) %>%
    summarise(ra_score = mean(score_pos_only, na.rm = TRUE), n_cells = n(), .groups = "drop")
  cd <- cliffs_delta(ra_animal$ra_score[ra_animal$group == "WildR"],
                     ra_animal$ra_score[ra_animal$group == "SPF"])
  sep <- complete_separation(ra_animal$ra_score, ra_animal$group)
  wx <- suppressWarnings(wilcox.test(ra_score ~ group, data = ra_animal, exact = FALSE))
  write_csv(ra_animal, file.path(tbl_dir, "VR4_ra_per_animal.csv"))
  write_md(file.path(out_dir, "VR4_ra_module.md"), c(
    "# VR-4 retinoic-acid module", "",
    sprintf("- Genes scored: %s.", paste(ra_genes, collapse = ", ")),
    sprintf("- Cliff's δ (WildR vs SPF): **%.3f** (95%% CI %.3f – %.3f).", cd$delta, cd$ci_low, cd$ci_high),
    sprintf("- Per-animal separation: %d / %d WildR exceed all SPF.", sep$n_exceed_all, sep$n_hi),
    sprintf("- Wilcoxon p (not headline): %.4g.", wx$p.value),
    "", "- Note: `Cyp26b1` is also a MAST hit and D-SPIN star at `0368 CEA-BST Six3 Cyp26b1 Gaba_2`.", ""
  ))
  log_msg("VR4 done. Top regulon: ", scan$regulon[1], " opposite limbs: ",
          sum(limb$opposite_direction, na.rm = TRUE))
}

# ══════════════════════════════════════════════════════════════════════════════
run_vr3 <- function() {
  log_msg("VR3 — gradient + MeApd")
  out_dir <- file.path(wkdir, "tables/VR/VR3")
  fig_dir <- file.path(wkdir, "figures/VR/VR3")
  ensure_dirs(out_dir, fig_dir)
  tbl_dir <- out_dir

  mast <- load_mast_hits(paths = paths) %>%
    filter(tissue == "AMY", !excluded_gene)
  chain_map <- mast %>%
    distinct(supertype_name, cell_type) %>%
    mutate(chain_position = assign_chain_position(supertype_name)) %>%
    filter(!is.na(chain_position)) %>%
    arrange(chain_position, supertype_name)
  write_csv(chain_map, file.path(tbl_dir, "VR3_chain_map.csv"))

  relay_hits <- mast %>%
    filter(chain_position %in% c("1", "2", "3", "4")) %>%
    mutate(pos = as.integer(chain_position))
  focal <- relay_hits %>% count(gene, name = "n_hits") %>% filter(n_hits >= 3) %>% pull(gene)
  log_msg("Focal genes (≥3 hits): ", length(focal))

  units <- list.files(paths$dspin, pattern = "^AMY__", full.names = TRUE)
  cov <- bind_rows(lapply(units, function(u) {
    sm <- file.path(u, "sample_meta.csv")
    cm <- file.path(u, "cells_meta.csv")
    if (!file.exists(sm)) return(NULL)
    meta <- readr::read_csv(sm, show_col_types = FALSE)
    lab <- sub("^AMY__", "", basename(u))
    med_umi <- NA_real_
    st_name <- NA_character_
    if (file.exists(cm)) {
      cells <- readr::read_csv(cm, show_col_types = FALSE)
      ucol <- intersect(c("nCount_RNA", "n_counts", "umi", "total_counts"), names(cells))[1]
      if (!is.na(ucol)) med_umi <- median(cells[[ucol]], na.rm = TRUE)
      st_name <- cells$supertype_name[1]
    }
    tibble(supertype_label = lab, supertype_name = st_name,
           n_cells = sum(meta$n_cells, na.rm = TRUE), median_UMI = med_umi,
           chain_position = assign_chain_position(st_name))
  }))
  write_csv(cov, file.path(tbl_dir, "VR3_supertype_covariates.csv"))
  cov2 <- cov %>%
    filter(chain_position %in% c("1", "2", "3", "4")) %>%
    mutate(pos = as.integer(chain_position), log_n = log1p(n_cells))
  umi_pos <- cov2 %>% filter(is.finite(median_UMI), is.finite(pos))
  rho_umi <- if (nrow(umi_pos) >= 4) {
    suppressWarnings(cor(umi_pos$median_UMI, umi_pos$pos, method = "spearman"))
  } else NA_real_

  grad_rows <- lapply(focal, function(g) {
    df <- relay_hits %>%
      filter(gene == g) %>%
      transmute(supertype_name, pos, coef, gene) %>%
      left_join(cov2 %>% select(supertype_name, median_UMI, log_n, n_cells),
                by = "supertype_name") %>%
      filter(is.finite(coef), is.finite(pos))
    if (nrow(df) < 3) return(NULL)
    fit_u <- lm(coef ~ pos, data = df)
    ci_u <- suppressWarnings(confint(fit_u))
    dfa <- df %>% filter(is.finite(median_UMI))
    if (nrow(dfa) >= 4) {
      fit_a <- tryCatch(lm(coef ~ pos + median_UMI + log_n, data = dfa), error = function(e) NULL)
    } else fit_a <- NULL
    if (is.null(fit_a)) {
      slope_a <- ci_a_lo <- ci_a_hi <- NA_real_
    } else {
      ci_a <- suppressWarnings(confint(fit_a))
      slope_a <- unname(coef(fit_a)[["pos"]])
      ci_a_lo <- unname(ci_a["pos", 1]); ci_a_hi <- unname(ci_a["pos", 2])
    }
    tibble(gene = g, n_hits = nrow(df),
           slope = unname(coef(fit_u)[["pos"]]),
           ci_low = unname(ci_u["pos", 1]), ci_high = unname(ci_u["pos", 2]),
           slope_adj = slope_a, ci_low_adj = ci_a_lo, ci_high_adj = ci_a_hi)
  })
  grad <- bind_rows(grad_rows)
  write_csv(grad, file.path(tbl_dir, "VR3_gradient_table.csv"))
  writeLines(
    sprintf("Spearman(median_UMI, chain position) = %.3f%s", rho_umi,
            if (is.finite(rho_umi) && abs(rho_umi) > 0.5)
              " — |ρ|>0.5: gradient and detection power confounded." else ""),
    file.path(tbl_dir, "VR3_umi_position_correlation.txt")
  )

  node <- relay_hits %>%
    group_by(pos = chain_position) %>%
    summarise(n = n(), k = sum(direction == "WildR-UP"), .groups = "drop") %>%
    rowwise() %>%
    mutate(est = wilson_ci(k, n)["est"], low = wilson_ci(k, n)["low"],
           high = wilson_ci(k, n)["high"]) %>%
    ungroup() %>%
    mutate(pos_num = as.integer(pos))
  write_csv(node, file.path(tbl_dir, "VR3_node_fractions.csv"))
  p <- ggplot(node, aes(pos_num, est)) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = "grey70") +
    geom_pointrange(aes(ymin = low, ymax = high), colour = "#FF6C0C", size = 0.8) +
    geom_line(colour = "#FF6C0C", alpha = 0.5) +
    scale_x_continuous(breaks = 1:4) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
    labs(title = "WildR-UP fraction by chain position",
         x = "Chain position (1–4; 3G excluded)",
         y = "WildR-UP fraction (Wilson 95% CI)") +
    vr_brand_ggplot()
  ggsave(file.path(fig_dir, "VR3_node_fractions.png"), p, width = 6, height = 4.5, dpi = 200)

  # Dissection audit + MeApd
  mmc <- tryCatch(readr::read_csv(paths$mapmycells, show_col_types = FALSE, comment = "#"),
                  error = function(e) NULL)
  meta_sp <- readRDS(paths$spatial_meta) %>% as_tibble(rownames = "cell_id") %>%
    mutate(cell_id = as.character(cell_id))
  if (!is.null(mmc)) {
    id_col <- names(mmc)[1]
    name_col <- grep("supertype_name", names(mmc), value = TRUE)[1]
    joined <- if (!is.na(name_col)) {
      meta_sp %>% inner_join(
        mmc %>% rename(cell_id = !!id_col) %>%
          mutate(cell_id = as.character(cell_id)) %>%
          select(cell_id, supertype_name = !!name_col),
        by = "cell_id")
    } else meta_sp %>% mutate(supertype_name = NA_character_)
  } else {
    joined <- meta_sp %>% mutate(supertype_name = NA_character_)
  }
  hippo <- joined %>% filter(grepl("CA3|DG Glut|0076 |Hippocamp", coalesce(supertype_name, "")))
  hippo_frac <- nrow(hippo) / max(nrow(joined), 1)
  write_md(file.path(out_dir, "VR3_dissection_audit.md"), c(
    "# VR-3 dissection audit", "",
    "- Spatial reference: MapMyCells Allen CCN20230722 joined to seqFISH metadata.",
    sprintf("- Hippocampal-labelled cells: %d (%.2f%% of mapped spatial nuclei).",
            nrow(hippo), 100 * hippo_frac),
    "- MeA/COAp boundary accuracy pending Allen_CCF_rebuild polygon join within section.",
    "- `0076 CA3 Glut 2` excluded from relay interpretation.",
    ""
  ))

  mea_cells <- joined %>%
    filter(grepl("MEA Slc17a7|0055 |0056 |0057 |MEA-COA", coalesce(supertype_name, ""))) %>%
    group_by(roi, bregma) %>%
    mutate(y_med = median(center_y, na.rm = TRUE),
           compartment = if_else(center_y >= y_med, "MeApd_proxy", "MeApv_proxy")) %>%
    ungroup()

  if (!nrow(mea_cells)) {
    write_md(file.path(out_dir, "VR3_meapd_test.md"), c(
      "# VR-3 MeApd vs MeApv test", "",
      "- **Blocked:** MapMyCells `output.csv` uses snRNA-style cell IDs (`47_02_41__s1`), while seqFISH metadata uses `roi*_run*_n` IDs — zero join.",
      "- seqFISH MapMyCells is only present for `clusters_roi4_block4` (partial).",
      "- Deviation: MeApd/MeApv animal-level test not computed until full seqFISH Allen labels are joined.",
      "- Gradient / node-fraction / hippocampal-audit outputs above are still valid.",
      ""
    ))
    write_csv(tibble(note = "no MeA cells after MapMyCells join"), file.path(tbl_dir, "VR3_meapd_animal_ratios.csv"))
    # Still write a background spatial map of all cells for audit
    p2 <- ggplot(meta_sp %>% slice_sample(n = min(50000L, nrow(meta_sp))),
                 aes(center_x, center_y)) +
      geom_point(size = 0.05, colour = "grey70") +
      facet_wrap(~roi) +
      labs(title = "seqFISH sections (MeA labels unavailable — ID join failed)") +
      vr_brand_ggplot()
    ggsave(file.path(fig_dir, "VR3_meapd_map.png"), p2, width = 12, height = 8, dpi = 200)
    log_msg("VR3 done (MeApd blocked on spatial label join). Focal genes=", length(focal),
            " UMI~pos ρ=", round(rho_umi, 3))
    return(invisible(NULL))
  }

  dv <- mea_cells %>%
    mutate(animal = as.character(group),
           microbiome = case_when(
             grepl("WildR|wild", group, ignore.case = TRUE) ~ "WildR",
             grepl("SPF|spf", group, ignore.case = TRUE) ~ "SPF",
             TRUE ~ as.character(group)
           )) %>%
    count(animal, microbiome, roi, compartment, name = "n") %>%
    tidyr::pivot_wider(names_from = compartment, values_from = n, values_fill = 0)
  if (!"MeApd_proxy" %in% names(dv)) dv$MeApd_proxy <- 0L
  if (!"MeApv_proxy" %in% names(dv)) dv$MeApv_proxy <- 0L
  dv <- dv %>%
    mutate(ratio_d_v = MeApd_proxy / pmax(MeApv_proxy, 1),
           n_cells = MeApd_proxy + MeApv_proxy,
           flag_low_n = n_cells < 50)
  dv_animal <- dv %>%
    group_by(animal, microbiome) %>%
    summarise(ratio = mean(ratio_d_v, na.rm = TRUE), n_cells = sum(n_cells),
              n_sections = n(), any_low_n = any(flag_low_n), .groups = "drop")
  write_csv(dv, file.path(tbl_dir, "VR3_meapd_section_counts.csv"))
  write_csv(dv_animal, file.path(tbl_dir, "VR3_meapd_animal_ratios.csv"))
  cd <- tryCatch(
    cliffs_delta(dv_animal$ratio[dv_animal$microbiome == "WildR"],
                 dv_animal$ratio[dv_animal$microbiome == "SPF"]),
    error = function(e) list(delta = NA_real_, ci_low = NA_real_, ci_high = NA_real_)
  )
  write_md(file.path(out_dir, "VR3_meapd_test.md"), c(
    "# VR-3 MeApd vs MeApv test", "",
    "- Partition: within section (roi × bregma), dorsal = center_y ≥ section median.",
    sprintf("- Cliff's δ (WildR vs SPF on per-animal D/V ratio): **%.3f** (CI %.3f – %.3f).",
            cd$delta, cd$ci_low, cd$ci_high),
    sprintf("- Animals / rows in ratio table: %d.", nrow(dv_animal)),
    "- Marker agreement (Esr1/Cckar/Ndnf) pending spatial gene scoring.",
    ""
  ))
  if (nrow(mea_cells)) {
    p2 <- ggplot(mea_cells, aes(center_x, center_y, colour = compartment)) +
      geom_point(size = 0.25, alpha = 0.8) +
      facet_wrap(~roi, scales = "free") + coord_equal() +
      scale_colour_manual(values = c(MeApd_proxy = "#FF6C0C", MeApv_proxy = "#27455C")) +
      labs(title = "MeA cells — within-section dorsal/ventral proxy", colour = NULL) +
      vr_brand_ggplot() + theme(legend.position = "bottom")
    ggsave(file.path(fig_dir, "VR3_meapd_map.png"), p2, width = 12, height = 8, dpi = 200)
  }
  log_msg("VR3 done. Focal genes=", length(focal), " UMI~pos ρ=", round(rho_umi, 3))
}
# ── dispatch ─────────────────────────────────────────────────────────────────
log_msg("Targets: ", paste(targets, collapse = ", "))
for (t in targets) {
  if (t %in% c("VR4_step1", "vr4_step1")) run_vr4_step1()
  else if (t %in% c("VR1", "vr1")) run_vr1()
  else if (t %in% c("VR2", "vr2")) run_vr2()
  else if (t %in% c("VR4", "vr4")) { run_vr4_step1(); run_vr4() }
  else if (t %in% c("VR3", "vr3")) run_vr3()
  else stop("Unknown target: ", t)
}
log_msg("All requested targets finished.")
