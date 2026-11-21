#!/usr/bin/env Rscript
# Upgrade of the integrated MAST × pySCENIC × D-SPIN × popAlign analysis with the multi-level
# pySCENIC GRNs (class 36 / subclass 137 / supertype 364; animal-level tests 2026-10-06,
# permutation tests 2026-10-07). Report: notebooks/0X_Integrated-Analysis_SN.qmd.
#
# The base run (run_2026-10-06) linked every popAlign unit to its parent *class* GRN. Here each
# unit is linked to the GRN fitted on exactly its own nuclei (class unit → class GRN, subclass
# unit → subclass GRN, supertype unit → supertype GRN; 244/244 matched, barcodes identical).
# Nothing upstream is re-run: reads saved AUCell matrices, regulons, popAlign outputs and the
# base-run tables. Outputs go to a NEW run folder; run_2026-10-06 is left untouched.
#
# Stages
#   units      per-unit native-GRN regulon tests (from tables/SCENIC) + split permutation p
#   comps      per animal × popAlign subpopulation regulon AUC (native GRN): pooled t per regulon,
#              exact animal-relabelling null for the count of regulons at p < 0.05
#   bridge     gene bridge redone with native-GRN regulon targets
#   animal     pySCENIC leave-one-animal-out scores + cross-method coupling with native GRNs
#   circuit    network enrichment for native-GRN regulons (same model as the base run)
#   effects    high-level: animal-level effect-size distributions of every method vs the 6 vs 6
#              null, method-of-moments estimate of the true effect spread, power
#   concord    high-level: do methods agree on WHICH supertypes carry signal? (size-adjusted)
#   evidence   popAlign leads re-scored with native-GRN pySCENIC
#
# Run (≈ 10 min, 1 core, < 8 GB):
#   /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript \
#     notebooks/R_scripts/run_integrated_scenic_multilevel_sn.R [--stage all] [--out-tag 2026-10-07]

suppressPackageStartupMessages({ library(tidyverse); library(glue); library(data.table) })
wk <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
source(file.path(wk, "notebooks/R_scripts/integrated_sn_helpers.R"))
source(file.path(wk, "notebooks/R_scripts/scenic_animal_level_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[[i + 1]] }
stage <- get_arg("--stage", "all")
P <- int_paths("2026-10-06")                                  # base run (read only)
out_tag <- get_arg("--out-tag", "2026-10-07")
O <- list(tbl = file.path(wk, "tables/SN_Integrated", glue("run_{out_tag}")),
          out = file.path(wk, "data/interim/sn_integrated", glue("run_{out_tag}")))
walk(O, dir.create, recursive = TRUE, showWarnings = FALSE)
scenic_date <- "2026-10-06"; perm_date <- "2026-10-07"
run_stage <- function(s) stage %in% c(s, "all")
wtbl <- function(x, f) { write_csv(x, file.path(O$tbl, f)); log_msg("  wrote ", f, " (", nrow(x), " rows)"); invisible(x) }
rtbl <- function(f, dir = O$tbl) read_csv(file.path(dir, f), show_col_types = FALSE)

# ── unit ↔ native GRN map ────────────────────────────────────────────────────
man <- read_csv(P$pa_manifest, show_col_types = FALSE) %>%
  mutate(code = sub("^CS20230722_", "", label),
         stem = if_else(level == "class", glue("snRNASeq_{tissue}_{label}"), glue("snRNASeq_{tissue}_{code}")),
         grn_dir = file.path(wk, "data/interim/scenic", glue("by_{level}"), stem),
         circuit = in_circuit(name))
catalog <- read_csv(file.path(P$pa_tbl, "popAlign_subpopulation_catalog.csv"), show_col_types = FALSE)
suffix <- c(class = "", subclass = "_subclass", supertype = "_supertype")
scen_tests <- imap_dfr(suffix, ~ read_csv(file.path(wk, "tables/SCENIC", glue("SCENIC_animal_level_tests{.x}_{scenic_date}.csv")),
                                          show_col_types = FALSE, col_select = c(stem, regulon, tf, delta_animal, hedges_g, p_animal,
                                                                                 fdr_split, n_animals_WildR, n_animals_SPF, n_cells)) %>%
                         mutate(grn_level = .y))
perm_split <- read_csv(file.path(wk, "tables/SCENIC", glue("SCENIC_split_permutation_tests_{perm_date}.csv")), show_col_types = FALSE)

# pooled-t hit count (p < alpha) for every labelling column of L (animals × labellings, 1 = WildR)
nominal_counts <- function(m, L, alpha = 0.05) {
  nx <- colSums(L); ny <- nrow(L) - nx
  sx <- m %*% L; sy <- rowSums(m) - sx; qx <- (m^2) %*% L; qy <- rowSums(m^2) - qx
  mx <- sweep(sx, 2, nx, "/"); my <- sweep(sy, 2, ny, "/")
  ssx <- qx - sweep(sx^2, 2, nx, "/"); ssy <- qy - sweep(sy^2, 2, ny, "/")
  df <- nx + ny - 2; se <- sqrt(sweep(sweep(ssx + ssy, 2, df, "/"), 2, 1 / nx + 1 / ny, "*"))
  hit <- sweep(abs(mx - my) / se, 2, qt(1 - alpha / 2, df), ">"); hit[!is.finite(se) | se <= 1e-15] <- FALSE
  colSums(hit, na.rm = TRUE)
}
labellings <- function(n, k) { cb <- combn(n, k); L <- matrix(0L, n, ncol(cb)); for (j in seq_len(ncol(cb))) L[cb[, j], j] <- 1L; L }
hedges <- function(X, g) {   # animals × features → Hedges g (same J as scenic_animal_level_helpers::row_tests)
  nx <- sum(g); ny <- sum(!g); mx <- colMeans(X[g, , drop = FALSE]); my <- colMeans(X[!g, , drop = FALSE])
  sp <- sqrt(((nx - 1) * apply(X[g, , drop = FALSE], 2, var) + (ny - 1) * apply(X[!g, , drop = FALSE], 2, var)) / (nx + ny - 2))
  j <- 1 - 3 / (4 * (nx + ny) - 9); ifelse(sp > 0, j * (mx - my) / sp, NA_real_)
}
read_regulons <- function(f) {   # activating regulons, union of targets over motifs
  r <- fread(f, skip = 3, header = FALSE, sep = ",", quote = "\"", fill = TRUE,
             col.names = c("TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue", "OrthologousIdentity",
                           "Annotation", "Context", "TargetGenes", "RankAtMax"))
  r <- r[grepl("activating", Context)]
  r[, targets := str_match_all(TargetGenes, "\\('([^']+)'")]
  r %>% as_tibble() %>% group_by(TF) %>%
    summarise(targets = list(unique(unlist(map(targets, ~ .x[, 2])))), .groups = "drop") %>% mutate(regulon = glue("{TF}(+)"))
}

# ════════════════════════════════════════════════════════════════════════════
# units · native-GRN regulon tests per popAlign unit
# ════════════════════════════════════════════════════════════════════════════
if (run_stage("units")) {
  log_msg("units: ", nrow(man), " popAlign units")
  reg <- scen_tests %>% inner_join(man %>% dplyr::select(unit_id, stem, level, tissue, name, circuit), by = "stem")
  wtbl(reg %>% dplyr::select(unit_id, level, tissue, name, circuit, stem, regulon, tf, delta_animal, hedges_g, p_animal, fdr_split,
                             n_animals_WildR, n_animals_SPF), "integrated_scenic_native_regulon.csv")
  us <- man %>% dplyr::select(unit_id, level, tissue, name, circuit, stem, n_cells_exported) %>%
    left_join(reg %>% filter(!is.na(p_animal)) %>% group_by(unit_id) %>%
                summarise(n_regulons = n(), n_p05 = sum(p_animal < 0.05), n_fdr10 = sum(fdr_split < 0.1, na.rm = TRUE),
                          grn_nuclei = first(n_cells), .groups = "drop"), by = "unit_id") %>%
    left_join(perm_split %>% dplyr::select(stem, p_perm_split = p_perm, perm_mean_split = perm_mean, fdr_perm_level), by = "stem") %>%
    mutate(excess = n_p05 / (0.05 * n_regulons))
  wtbl(us, "integrated_scenic_native_unit.csv")
}

# ════════════════════════════════════════════════════════════════════════════
# comps · regulon AUC per animal × popAlign subpopulation, native GRN
# ════════════════════════════════════════════════════════════════════════════
if (run_stage("comps")) {
  log_msg("comps: per animal × subpopulation regulon AUC from native GRNs")
  comp_reg <- list(); comp_sum <- list(); animal_auc <- list()
  for (i in seq_len(nrow(man))) {
    u <- man[i, ]
    f <- file.path(u$grn_dir, "auc_mtx.csv"); if (!file.exists(f)) next
    a <- fread(f); setnames(a, 1, "cell"); setkey(a, cell)
    cells <- fread(file.path(u$output_dir, "cells.csv.gz"), select = c("cell_barcode", "sample_id", "microbiome", "comp"))
    a <- a[J(cells$cell_barcode), nomatch = 0L]; reg_cols <- setdiff(names(a), "cell")
    a <- merge(a, cells, by.x = "cell", by.y = "cell_barcode")
    robust_comps <- catalog %>% filter(unit_id == u$unit_id, robust, !contamination_flag) %>% pull(comp)
    # unit-level per-animal means over the popAlign nuclei (for the leave-one-animal-out scores)
    am <- a[, lapply(.SD, mean), by = .(sample_id, microbiome), .SDcols = reg_cols]
    am <- merge(am, a[, .(n_nuclei = .N), by = sample_id], by = "sample_id")
    animal_auc[[u$unit_id]] <- as_tibble(am) %>% mutate(unit_id = u$unit_id)
    for (cc in sort(unique(a$comp))) {
      sub <- a[a$comp == cc]
      ok_animals <- sub[, .(n_nuclei = .N), by = sample_id][n_nuclei >= 5]$sample_id
      ac <- sub[sub$sample_id %in% ok_animals, lapply(.SD, mean), by = .(sample_id, microbiome), .SDcols = reg_cols]
      g <- ac$microbiome == "WildR"
      if (sum(g) < 4 || sum(!g) < 4) next
      X <- as.matrix(ac[, ..reg_cols]); X <- X[, apply(X, 2, var) > 0, drop = FALSE]
      if (!ncol(X)) next
      rt <- row_tests(t(X), if_else(g, "WildR", "SPF"))
      L <- labellings(nrow(X), sum(g)); null <- nominal_counts(t(X), L)
      obs <- sum(rt$p_animal < 0.05, na.rm = TRUE)
      comp_reg[[glue("{u$unit_id}::{cc}")]] <- rt %>% transmute(unit_id = u$unit_id, comp = cc, regulon, hedges_g, p_animal)
      comp_sum[[glue("{u$unit_id}::{cc}")]] <- tibble(unit_id = u$unit_id, level = u$level, tissue = u$tissue, name = u$name,
        circuit = u$circuit, comp = cc, robust = cc %in% robust_comps, n_animals = nrow(X), n_WildR = sum(g),
        n_regulons = ncol(X), n_p05 = obs, perm_mean = mean(null), n_perm = length(null), p_perm_count = mean(null >= obs))
    }
    if (i %% 25 == 0) log_msg("  ", i, "/", nrow(man), " units")
  }
  cs <- bind_rows(comp_sum) %>% left_join(catalog %>% dplyr::select(unit_id, comp, p_state, q_state), by = c("unit_id", "comp"))
  wtbl(cs, "integrated_scenic_native_component.csv")
  wtbl(bind_rows(comp_reg), "integrated_scenic_native_component_regulon.csv")
  saveRDS(animal_auc, file.path(O$out, "animal_regulon_auc_native.rds"))
}

# ════════════════════════════════════════════════════════════════════════════
# bridge · do native-GRN regulons target the program popAlign sees moving?
# ════════════════════════════════════════════════════════════════════════════
if (run_stage("bridge")) {
  gb <- read_csv(file.path(P$tbl, "integrated_gene_bridge.csv"), show_col_types = FALSE)
  creg <- rtbl("integrated_scenic_native_component_regulon.csv")
  ureg <- rtbl("integrated_scenic_native_regulon.csv")
  units <- man %>% filter(unit_id %in% gb$unit_id)
  log_msg("bridge: ", nrow(units), " units, ", nrow(gb), " components")
  res <- map_dfr(seq_len(nrow(units)), function(i) {
    u <- units[i, ]; f <- file.path(u$grn_dir, "regulons.csv"); if (!file.exists(f)) return(NULL)
    regs <- read_regulons(f)
    W <- read_csv(file.path(u$output_dir, "onmf_W.csv.gz"), show_col_types = FALSE)
    mast <- load_mast_unit(P, u$level, u$tissue, u$label)
    U <- W$gene; if (!is.null(mast)) U <- intersect(U, mast$gene)
    top_genes <- function(fe) W %>% arrange(desc(.data[[fe]])) %>% filter(gene %in% U) %>% slice_head(n = 50) %>% pull(gene)
    rows <- gb %>% filter(unit_id == u$unit_id)
    map_dfr(seq_len(nrow(rows)), function(r) {
      row <- rows[r, ]; out <- tibble(unit_id = u$unit_id, comp = row$comp)
      cr <- creg %>% filter(unit_id == u$unit_id, comp == row$comp)
      for (arm in c("shift", "control")) {
        g <- top_genes(if (arm == "shift") row$shifted_feature else row$control_feature)
        en <- regs %>% mutate(hp = map(targets, ~ hyper_p(g, .x, U)), k = map_dbl(hp, "k"), p = map_dbl(hp, "p")) %>%
          mutate(q = p.adjust(p, "BH")) %>% arrange(p) %>%
          left_join(cr %>% dplyr::select(regulon, p_comp = p_animal), by = "regulon") %>%
          left_join(ureg %>% filter(unit_id == u$unit_id) %>% dplyr::select(regulon, p_unit = p_animal), by = "regulon")
        enr <- en %>% filter(q < 0.05)
        out[[glue("native_n_regulons_{arm}")]] <- nrow(en)
        out[[glue("native_n_enriched_{arm}")]] <- nrow(enr)
        out[[glue("native_top_{arm}")]] <- as.character(en$TF[1]); out[[glue("native_top_q_{arm}")]] <- en$q[1]
        out[[glue("native_top_p_unit_{arm}")]] <- en$p_unit[1]
        # of the regulons that target the program, how many shift per animal in its subpopulation?
        out[[glue("native_enriched_shifted_{arm}")]] <- sum(enr$p_comp < 0.05, na.rm = TRUE)
        out[[glue("native_enriched_tested_{arm}")]] <- sum(!is.na(enr$p_comp))
      }
      out
    })
  })
  wtbl(gb %>% dplyr::select(unit_id, level, tissue, name, comp, p_state, circuit, shifted_feature, control_feature, d_shift,
                            scenic_n_enriched_shift, scenic_n_enriched_control) %>% left_join(res, by = c("unit_id", "comp")),
       "integrated_gene_bridge_scenic_native.csv")
}

# ════════════════════════════════════════════════════════════════════════════
# animal · leave-one-animal-out pySCENIC scores (native GRNs) + coupling
# ════════════════════════════════════════════════════════════════════════════
loao_score <- function(X, is_wildr) {   # identical to run_integrated_analysis_sn.R
  X <- as.matrix(X); X <- X[, apply(X, 2, function(v) sd(v, na.rm = TRUE) > 0), drop = FALSE]
  if (!ncol(X)) return(rep(NA_real_, nrow(X)))
  Z <- scale(X)
  sapply(seq_len(nrow(Z)), function(a) {
    w <- colMeans(Z[-a, ][is_wildr[-a], , drop = FALSE]) - colMeans(Z[-a, ][!is_wildr[-a], , drop = FALSE])
    sum(Z[a, ] * sign(w) * abs(w)) / sum(abs(w))
  })
}
if (run_stage("animal")) {
  combos <- perm_combos()
  auc <- readRDS(file.path(O$out, "animal_regulon_auc_native.rds"))
  sup <- man %>% filter(level == "supertype")
  base_scores <- read_csv(file.path(P$tbl, "integrated_animal_scores.csv"), show_col_types = FALSE)
  new <- map_dfr(c("AMY", "HYP"), function(tis) map_dfr(c("network", "other"), function(sc) {
    us <- sup %>% filter(tissue == tis, circuit == (sc == "network"))
    per <- map_dfr(us$unit_id, function(uid) {
      a <- auc[[uid]]; if (is.null(a) || nrow(a) < 12) return(NULL)
      a <- arrange(a, microbiome, sample_id); rc <- setdiff(names(a), c("sample_id", "microbiome", "n_nuclei", "unit_id"))
      tibble(unit_id = uid, sample_id = a$sample_id, microbiome = a$microbiome, s = loao_score(a[, rc], a$microbiome == "WildR"))
    })
    per %>% group_by(unit_id) %>% mutate(s = as.numeric(scale(s))) %>% ungroup() %>%
      group_by(sample_id, microbiome) %>% summarise(score = mean(s, na.rm = TRUE), n_units = n_distinct(unit_id), .groups = "drop") %>%
      mutate(modality = "pySCENIC", scope = sc, tissue = tis)
  }))
  scores <- bind_rows(base_scores %>% filter(modality != "pySCENIC"), new)
  wtbl(scores, "integrated_animal_scores.csv")
  auc_of <- function(score, w) { r <- rank(score); (sum(r[w]) - sum(w) * (sum(w) + 1) / 2) / (sum(w) * sum(!w)) }
  pred <- scores %>% filter(!str_detect(modality, "^VR")) %>% group_by(tissue, scope, modality) %>%
    summarise(n = n(), accuracy = mean(sign(score) == if_else(microbiome == "WildR", 1, -1)),
              auc = auc_of(score, microbiome == "WildR"),
              p_one_sided = { w <- microbiome == "WildR"; obs <- auc_of(score, w)
                null <- apply(combos, 2, function(ix) auc_of(score, seq_along(score) %in% ix)); mean(null >= obs - 1e-12) },
              .groups = "drop")
  wtbl(pred, "integrated_animal_loao_prediction.csv")
  coup <- scores %>% group_by(tissue, scope) %>% group_modify(function(d, k) {
    w <- d %>% dplyr::select(sample_id, microbiome, modality, score) %>% pivot_wider(names_from = modality, values_from = score) %>% drop_na()
    mods <- setdiff(names(w), c("sample_id", "microbiome")); if (length(mods) < 2) return(tibble())
    R <- w %>% group_by(microbiome) %>% mutate(across(all_of(mods), ~ .x - mean(.x))) %>% ungroup()
    map_dfr(combn(mods, 2, simplify = FALSE), function(pr) {
      x <- R[[pr[1]]]; y <- R[[pr[2]]]; rho_w <- cor(x, y, method = "spearman")
      set.seed(1); null <- replicate(5000, { yy <- y; for (gg in unique(R$microbiome)) { ix <- which(R$microbiome == gg); yy[ix] <- y[sample(ix)] }; cor(x, yy, method = "spearman") })
      tibble(m1 = pr[1], m2 = pr[2], n = nrow(w), rho_all = cor(w[[pr[1]]], w[[pr[2]]], method = "spearman"),
             rho_within = rho_w, p_within = mean(abs(null) >= abs(rho_w) - 1e-12))
    })
  }) %>% ungroup()
  wtbl(coup, "integrated_animal_coupling.csv")
}

# ════════════════════════════════════════════════════════════════════════════
# circuit · network enrichment, native-GRN regulons (same model as the base run)
# ════════════════════════════════════════════════════════════════════════════
circuit_test <- function(d) {
  set.seed(42)
  obs <- mean(d$rate[d$circuit]) / mean(d$rate[!d$circuit])
  null <- replicate(10000, { cc <- d$circuit; for (t in unique(d$tissue)) { ix <- which(d$tissue == t); cc[ix] <- sample(cc[ix]) }
    mean(d$rate[cc]) / mean(d$rate[!cc]) })
  boot <- replicate(2000, mean(sample(d$rate[d$circuit], replace = TRUE)) / mean(sample(d$rate[!d$circuit], replace = TRUE)))
  m <- glm(cbind(n_hit, n_tests - n_hit) ~ circuit + log(n_cells_exported) + tissue, family = quasibinomial, data = d)
  cf <- coef(summary(m))["circuitTRUE", ]
  tibble(n_units = nrow(d), n_circuit = sum(d$circuit), n_other = sum(!d$circuit),
         rate_circuit = mean(d$rate[d$circuit]), rate_other = mean(d$rate[!d$circuit]), rate_ratio = obs,
         ci_lo = quantile(boot, 0.025, na.rm = TRUE), ci_hi = quantile(boot, 0.975, na.rm = TRUE),
         p_perm = (sum(abs(log(null)) >= abs(log(obs)) - 1e-12, na.rm = TRUE) + 1) / (sum(is.finite(null)) + 1),
         ratio_adj = exp(cf[[1]]), ratio_adj_lo = exp(cf[[1]] - 1.96 * cf[[2]]), ratio_adj_hi = exp(cf[[1]] + 1.96 * cf[[2]]), p_adj = cf[[4]])
}
if (run_stage("circuit")) {
  us <- rtbl("integrated_scenic_native_unit.csv") %>% filter(!is.na(n_regulons))
  base <- read_csv(file.path(P$tbl, "integrated_circuit_enrichment.csv"), show_col_types = FALSE)
  nat <- us %>% filter(level == "supertype") %>% transmute(unit_id, tissue, circuit, n_cells_exported, n_tests = n_regulons, n_hit = n_p05, rate = n_hit / n_tests)
  rows <- bind_rows(
    circuit_test(nat) %>% mutate(modality = "pySCENIC regulons (supertype GRN)", inference = "animals"),
    # all three levels: every popAlign unit with its own GRN, level as an extra covariate via stratified permutation
    { d <- us %>% transmute(unit_id, tissue = paste(tissue, level), circuit, n_cells_exported, n_tests = n_regulons, n_hit = n_p05, rate = n_hit / n_tests)
      circuit_test(d) %>% mutate(modality = "pySCENIC regulons (own GRN, all levels)", inference = "animals") })
  out <- bind_rows(base %>% mutate(modality = if_else(modality == "pySCENIC regulons", "pySCENIC regulons (class GRN)", modality)), rows)
  wtbl(out, "integrated_circuit_enrichment.csv")
  rates <- read_csv(file.path(P$tbl, "integrated_unit_rates.csv"), show_col_types = FALSE)
  wtbl(bind_rows(rates %>% mutate(modality = if_else(modality == "pySCENIC regulons", "pySCENIC regulons (class GRN)", modality)),
                 nat %>% left_join(man %>% dplyr::select(unit_id, name), by = "unit_id") %>%
                   mutate(modality = "pySCENIC regulons (supertype GRN)", inference = "animals", alpha = 0.05,
                          excess_flag = n_hit > qbinom(0.95, n_tests, 0.05), any_flag = n_hit > 0)),
       "integrated_unit_rates.csv")
}

# ════════════════════════════════════════════════════════════════════════════
# effects · animal-level effect sizes of every method vs the 6 vs 6 null
# ════════════════════════════════════════════════════════════════════════════
# For every animal-level test, Hedges g (WildR − SPF, pooled SD, small-sample J), restricted to
# 6 vs 6 so all methods share one design. Null = the same statistic under every relabelling of
# the tissue's 12 animals (924), applied jointly to all units and features of the tissue, so the
# correlation between genes / programs / regulons and between units is kept and no normality is
# assumed. Global statistic per method × tissue: mean g² over all tests. The excess of the observed
# mean g² over its permutation mean estimates the variance of the TRUE standardised effects
# (method of moments); its square root is the RMS true effect. Tissues have different animal
# sets, so the pooled null pairs AMY and HYP relabellings at random (10,000 pairs).
null_g2 <- function(nx, ny) { df <- nx + ny - 2; j <- 1 - 3 / (4 * (nx + ny) - 9); j^2 * df / (df - 2) * (1 / nx + 1 / ny) }   # theory, for reference
if (run_stage("effects")) {
  sup <- man %>% filter(level == "supertype")
  lc <- readRDS(file.path(P$out, "pseudobulk_unit_logcpm.rds"))
  mast_deg <- map_dfr(names(lc), function(uid) { p <- str_split(uid, "__")[[1]]; m <- load_mast_unit(P, "supertype", p[2], p[3])
    if (is.null(m)) NULL else tibble(unit_id = uid, gene = m$gene[m$padj < 0.1]) })
  # animals × features matrices (rownames = sample_id), one per method × unit
  mats <- list(); add <- function(method, uid, X) if (!is.null(X) && ncol(X) > 0 && nrow(X) == 12 &&
                                                        sum(str_detect(rownames(X), "^WildR")) == 6) mats[[length(mats) + 1]] <<- list(method = method, unit_id = uid, X = X)
  log_msg("effects: pseudobulk genes")
  for (uid in names(lc)) {
    X <- t(lc[[uid]]); dg <- mast_deg$gene[mast_deg$unit_id == uid]
    add("Pseudobulk: MAST DEGs", uid, X[, colnames(X) %in% dg, drop = FALSE])
    add("Pseudobulk: other genes", uid, X[, !colnames(X) %in% dg, drop = FALSE])
  }
  log_msg("effects: D-SPIN programs")
  for (i in seq_len(nrow(sup))) {
    d <- load_dspin_unit(P, sup$tissue[i], sup$label[i]); if (is.null(d)) next
    w <- d$responses %>% dplyr::select(sample_id, program_id, relative_h) %>% pivot_wider(names_from = program_id, values_from = relative_h)
    X <- as.matrix(w[, -1]); rownames(X) <- w$sample_id; add("D-SPIN programs", sup$unit_id[i], X)
  }
  log_msg("effects: popAlign subpopulation programs")
  for (i in seq_len(nrow(sup))) {
    rc <- catalog %>% filter(unit_id == sup$unit_id[i], robust, !contamination_flag) %>% pull(comp); if (!length(rc)) next
    ce <- read_csv(file.path(sup$output_dir[i], "animal_component_centroids.csv"), show_col_types = FALSE) %>% filter(comp %in% rc)
    f <- grep("^F[0-9]+$", names(ce), value = TRUE)
    w <- ce %>% dplyr::select(sample_id, comp, all_of(f)) %>% pivot_wider(names_from = comp, values_from = all_of(f), names_glue = "{comp}:{.value}")
    X <- as.matrix(w[, -1]); rownames(X) <- w$sample_id; X <- X[, colSums(is.na(X)) == 0, drop = FALSE]
    add("popAlign programs", sup$unit_id[i], X)
  }
  log_msg("effects: pySCENIC regulons (supertype GRNs)")
  pa_sc <- fread(file.path(wk, "data/interim/scenic/animal_level", glue("regulon_auc_per_animal_supertype_{scenic_date}.csv")))[n_cells >= SCENIC_MIN_CELLS]
  pa_sc[, unit_id := paste0("supertype__", tissue, "__CS20230722_", sub("^snRNASeq_(AMY|HYP)_", "", stem))]
  pa_sc[, sample_id := paste0(animal, "_", tissue)]
  for (uid in intersect(unique(pa_sc$unit_id), sup$unit_id)) {
    w <- dcast(pa_sc[unit_id == uid], sample_id ~ regulon, value.var = "mean_auc")
    X <- as.matrix(w[, -1]); rownames(X) <- w$sample_id; X <- X[, colSums(is.na(X)) == 0, drop = FALSE]
    add("pySCENIC regulons", uid, X)
  }
  log_msg("effects: ", length(mats), " method × unit matrices; permuting")
  tis_of <- function(uid) str_split(uid, "__")[[1]][2]
  # g for every feature (rows) × labelling (cols), same formula as hedges()
  g_lab <- function(X, L) {
    m <- t(X); nx <- colSums(L); ny <- nrow(L) - nx
    sx <- m %*% L; sy <- rowSums(m) - sx; qx <- (m^2) %*% L; qy <- rowSums(m^2) - qx
    mx <- sweep(sx, 2, nx, "/"); my <- sweep(sy, 2, ny, "/")
    ssx <- qx - sweep(sx^2, 2, nx, "/"); ssy <- qy - sweep(sy^2, 2, ny, "/")
    sp <- sqrt(pmax(sweep(ssx + ssy, 2, nx + ny - 2, "/"), 0)); j <- 1 - 3 / (4 * 12 - 9)
    g <- j * (mx - my) / sp; g[!is.finite(g)] <- NA; g
  }
  Ls <- list(); obs_rows <- list(); null_rows <- list()
  for (k in seq_along(mats)) {
    mm <- mats[[k]]; tis <- tis_of(mm$unit_id); X <- mm$X[order(rownames(mm$X)), , drop = FALSE]
    X <- X[, apply(X, 2, var) > 0, drop = FALSE]; if (!ncol(X)) next
    key <- paste(rownames(X), collapse = ",")
    if (is.null(Ls[[key]])) { L <- labellings(12, 6); rownames(L) <- rownames(X)
      truth <- as.integer(str_detect(rownames(X), "^WildR")); obs_col <- which(colSums(abs(L - truth)) == 0)
      # canonical labelling order per tissue: identify each labelling by its WildR set
      Ls[[key]] <- list(L = L, obs = obs_col, id = apply(L, 2, function(v) paste(sort(rownames(X)[v == 1]), collapse = ","))) }
    G <- g_lab(X, Ls[[key]]$L)
    obs_rows[[k]] <- tibble(method = mm$method, unit_id = mm$unit_id, tissue = tis, feature = colnames(X), g = G[, Ls[[key]]$obs])
    null_rows[[k]] <- tibble(method = mm$method, tissue = tis, lab = Ls[[key]]$id, s2 = colSums(G^2, na.rm = TRUE), n = colSums(!is.na(G)),
                             observed = seq_len(ncol(G)) == Ls[[key]]$obs)
  }
  eff <- bind_rows(obs_rows) %>% filter(is.finite(g))
  saveRDS(eff, file.path(O$out, "animal_effect_sizes_6v6.rds"))   # per-method observed g (no "all genes" duplicate)
  # per method × tissue: mean g² under every relabelling (summed over units sharing that labelling)
  # "Pseudobulk: all genes" = MAST DEGs + other genes. Only this pseudobulk row is a valid global
  # test: the MAST-DEG subset was selected with the true labels on the same animals' nuclei, so its
  # excess over relabellings is expected by construction (double dipping), not evidence of effect.
  nr <- bind_rows(null_rows)
  nr <- bind_rows(nr, nr %>% filter(str_detect(method, "^Pseudobulk")) %>% mutate(method = "Pseudobulk: all genes"))
  eff <- bind_rows(eff, eff %>% filter(str_detect(method, "^Pseudobulk")) %>% mutate(method = "Pseudobulk: all genes"))
  nl <- nr %>% group_by(method, tissue, lab) %>%
    summarise(mean_g2 = sum(s2) / sum(n), n = sum(n), observed = any(observed), .groups = "drop")
  write_csv(nl, file.path(O$out, "effect_size_permutation_null.csv"))
  set.seed(7)
  pooled <- nl %>% group_by(method) %>% group_modify(function(d, k) {
    a <- d %>% filter(tissue == "AMY"); h <- d %>% filter(tissue == "HYP")
    nA <- max(a$n); nH <- max(h$n)
    comb <- function(x, y) (x * nA + y * nH) / (nA + nH)
    obs <- comb(a$mean_g2[a$observed], h$mean_g2[h$observed])
    null <- comb(sample(a$mean_g2, 10000, TRUE), sample(h$mean_g2, 10000, TRUE))
    tibble(tissue = "Both", obs_g2 = obs, perm_g2 = mean(null), perm_q05 = quantile(null, 0.05), perm_q95 = quantile(null, 0.95),
           p_perm = (sum(null >= obs) + 1) / (length(null) + 1))
  }) %>% ungroup()
  per_tis <- nl %>% group_by(method, tissue) %>%
    summarise(obs_g2 = mean_g2[observed], perm_g2 = mean(mean_g2), perm_q05 = quantile(mean_g2, 0.05), perm_q95 = quantile(mean_g2, 0.95),
              p_perm = mean(mean_g2 >= mean_g2[observed]), .groups = "drop")
  n0 <- null_g2(6, 6); g_crit <- qt(0.975, 10) * sqrt(1 / 3) * (1 - 3 / (4 * 12 - 9))
  summ <- bind_rows(per_tis, pooled) %>%
    left_join(bind_rows(eff, eff %>% mutate(tissue = "Both")) %>% group_by(method, tissue) %>%
                summarise(n_tests = n(), n_units = n_distinct(unit_id), median_abs_g = median(abs(g)), share_p05 = mean(abs(g) > g_crit), .groups = "drop"),
              by = c("method", "tissue")) %>%
    mutate(theory_g2 = n0, ratio_perm = obs_g2 / perm_g2, var_excess = obs_g2 - perm_g2, rms_true = sqrt(pmax(var_excess, 0)),
           rms_upper95 = sqrt(pmax(obs_g2 - perm_q05, 0)),   # 95% upper bound: observed minus the null's 5th percentile
           n_per_group_80 = map_dbl(rms_true, ~ if (.x > 0.05) ceiling(power.t.test(delta = .x, power = 0.8, sig.level = 0.05)$n) else NA_real_))
  # cluster bootstrap over units for the RMS true effect (pooled tissues; permutation mean held fixed)
  set.seed(8)
  bs <- eff %>% group_by(method) %>% group_modify(function(d, k) {
    pu <- d %>% group_by(unit_id) %>% summarise(s2 = sum(g^2), n = n(), .groups = "drop")
    pm <- summ$perm_g2[summ$method == k$method & summ$tissue == "Both"]
    b <- replicate(1000, { ix <- sample(nrow(pu), replace = TRUE); sum(pu$s2[ix]) / sum(pu$n[ix]) - pm })
    tibble(tissue = "Both", rms_lo = sqrt(max(quantile(b, 0.025), 0)), rms_hi = sqrt(max(quantile(b, 0.975), 0)))
  }) %>% ungroup()
  summ <- summ %>% left_join(bs, by = c("method", "tissue")) %>% arrange(tissue, method)
  wtbl(summ, "integrated_effect_size_summary.csv")
  # |g| densities per method and the permutation-null |g| density (all relabellings except the observed one, pooled)
  grid <- seq(0, 4, length.out = 201)
  dens_of <- function(x) { dd <- density(abs(x), from = 0, to = 4, n = length(grid), bw = 0.08); dd$y }
  dens <- eff %>% group_by(method) %>% summarise(abs_g = list(grid), density = list(dens_of(g)), .groups = "drop") %>% unnest(c(abs_g, density))
  s <- sqrt(1 / 3) * (1 - 3 / (4 * 12 - 9))
  dens <- bind_rows(dens, tibble(method = "Null: t-distribution, 6 vs 6", abs_g = grid, density = 2 * dt(grid / s, 10) / s))
  wtbl(dens, "integrated_effect_size_density.csv")
  pw <- expand_grid(g = c(0.5, 0.8, 1.2, 1.8), alpha = c(0.05, 1e-4), n = c(3:30, seq(35, 100, 5))) %>%
    mutate(power = pmap_dbl(list(n, g, alpha), function(n, g, a) power.t.test(n = n, delta = g, sig.level = a)$power))
  wtbl(pw, "integrated_power_curves.csv")
}

# ════════════════════════════════════════════════════════════════════════════
# concord · do methods agree on WHICH supertypes carry signal?
# ════════════════════════════════════════════════════════════════════════════
# Per supertype unit, each method's hit rate (share of its tests at p < 0.05; MAST padj < 0.1).
# Spearman across units, raw and partial (ranks residualised on log nuclei and tissue, since
# unit size drives power for every method). Permutation p: shuffle one method's residual ranks
# within tissue (5,000).
if (run_stage("concord")) {
  rates <- rtbl("integrated_unit_rates.csv") %>% filter(modality != "pySCENIC regulons (class GRN)")
  W <- rates %>% dplyr::select(unit_id, tissue, n_cells_exported, modality, rate) %>% pivot_wider(names_from = modality, values_from = rate)
  mods <- setdiff(names(W), c("unit_id", "tissue", "n_cells_exported"))
  resid_rank <- function(x, d) { ok <- !is.na(x); r <- rep(NA_real_, length(x))
    r[ok] <- resid(lm(rank(x[ok]) ~ log(d$n_cells_exported[ok]) + d$tissue[ok])); r }
  set.seed(11)
  conc <- map_dfr(combn(mods, 2, simplify = FALSE), function(pr) {
    d <- W %>% filter(!is.na(.data[[pr[1]]]), !is.na(.data[[pr[2]]]))
    x <- resid_rank(d[[pr[1]]], d); y <- resid_rank(d[[pr[2]]], d); rp <- cor(x, y)
    null <- replicate(5000, { yy <- y; for (t in unique(d$tissue)) { ix <- which(d$tissue == t); yy[ix] <- y[sample(ix)] }; cor(x, yy) })
    tibble(m1 = pr[1], m2 = pr[2], n_units = nrow(d), rho_raw = cor(d[[pr[1]]], d[[pr[2]]], method = "spearman"),
           rho_size = cor(rank(d[[pr[1]]]), log(d$n_cells_exported)), rho_partial = rp, p_partial = mean(abs(null) >= abs(rp) - 1e-12))
  })
  wtbl(conc, "integrated_unit_concordance.csv")
  # size dependence of each method's hit rate (power scales with nuclei for nuclei-level methods)
  wtbl(rates %>% group_by(modality, inference) %>%
         summarise(n_units = n(), rho_size = cor(rate, log(n_cells_exported), method = "spearman", use = "complete.obs"),
                   p_size = cor.test(rate, log(n_cells_exported), method = "spearman", exact = FALSE)$p.value, .groups = "drop"),
       "integrated_unit_size_dependence.csv")
}

# ════════════════════════════════════════════════════════════════════════════
# evidence · popAlign leads re-scored with native-GRN pySCENIC
# ════════════════════════════════════════════════════════════════════════════
if (run_stage("evidence")) {
  leads <- read_csv(file.path(P$pa_tbl, "popAlign_followup_leads.csv"), show_col_types = FALSE)
  cs <- rtbl("integrated_scenic_native_component.csv"); us <- rtbl("integrated_scenic_native_unit.csv")
  br <- rtbl("integrated_gene_bridge_scenic_native.csv")
  creg <- rtbl("integrated_scenic_native_component_regulon.csv")
  top <- creg %>% group_by(unit_id, comp) %>% slice_min(p_animal, n = 3, with_ties = FALSE) %>%
    summarise(top_regulons = paste(glue("{regulon} (g = {sprintf('%+.1f', hedges_g)})"), collapse = ", "), .groups = "drop")
  ev <- leads %>% dplyr::select(tissue, lineage, level, name, unit_id, comp, effect, p, depth_flag) %>%
    left_join(cs %>% dplyr::select(unit_id, comp, comp_n_regulons = n_regulons, comp_n_p05 = n_p05, comp_perm_mean = perm_mean,
                                   comp_p_perm = p_perm_count), by = c("unit_id", "comp")) %>%
    left_join(us %>% dplyr::select(unit_id, grn_nuclei, unit_n_regulons = n_regulons, unit_n_p05 = n_p05, unit_p_perm = p_perm_split), by = "unit_id") %>%
    left_join(br %>% dplyr::select(unit_id, comp, native_n_enriched_shift, native_top_shift, native_enriched_shifted_shift,
                                   native_enriched_tested_shift), by = c("unit_id", "comp")) %>%
    left_join(top, by = c("unit_id", "comp"))
  wtbl(ev, "integrated_evidence_leads_scenic_native.csv")
}
log_msg("Stage '", stage, "' finished → ", O$tbl)
