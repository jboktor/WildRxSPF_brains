#!/usr/bin/env Rscript
# Animal-level summary of the D-SPIN snRNA-seq fits (WildR vs SPF) for 0X_DSPIN_Analysis.qmd.
# Light: reads the saved per-unit fit outputs in data/interim/dspin/<tissue>__<supertype_label>/
# (no refitting, no DAVID calls). ~1-2 min on a login node.
#
# What it computes
#   * Program-level test: for each unit × program, WildR vs SPF difference in relative h
#     (relative_h = h_animal − mean SPF h in that unit, from DSPIN.response_relative_to_control).
#     The animal is the unit of replication. Two-sided exact permutation test of |mean_WildR −
#     mean_SPF| over all choose(12, 6) = 924 relabelings of the animals (min p = 2/924 = 0.0022),
#     plus Hedges' g and "separation" (# WildR animals beyond the SPF range in the direction of
#     the shift). BH q across all unit × program tests.
#   * Unit-level omnibus: sum over the unit's programs of squared Welch t, permuted with the same
#     924 animal relabelings (keeps program-program correlation). BH q across units.
#   * Calibration: observed vs permutation-expected counts of small p (is there any excess signal?).
#   * Noise floor: mean |relative h| of SPF animals vs WildR animals per unit (the old ranking
#     metric used WildR |h| alone, which also grows when a unit has few nuclei per animal).
#     SPF animals are measured against a mean they contributed to, so their raw |h| is biased
#     low; abs_h_SPF_loo rescales SPF |h| by n/(n−1) (= deviation from the leave-one-out SPF
#     mean). Under no microbiome effect E[abs_h_ratio_loo] = sqrt((n+1)(n−1))/n ≈ 0.99 (n = 6)
#     and E[abs_h_ratio] = sqrt((n+1)/(n−1)) ≈ 1.18.
#   * Program annotation: top genes by oNMF weight; DAVID best curated term restricted to
#     GO BP/MF/CC + KEGG (UniProt sequence features / keywords dropped), Benjamini-adjusted.
#   * MAST concordance: member genes that are MAST hits (supertype level, padj < 0.1, the
#     2026-08-17 consolidated table) and whether their mean coef agrees in sign with Δh.
#
# Pre-specified lead rule (mirrors popAlign; nothing survives FDR at 6 vs 6 animals):
#   p_perm ≤ 0.0065 (≤ 6/924) and separation ≥ 4 WildR animals.
#
# Usage
#   Rscript notebooks/R_scripts/summarize_dspin_responses_sn.R --tag 2026-10-06
#
# Outputs (date-stamped, never overwrite)
#   data/interim/dspin/summary_<tag>/dspin_program_tests_<tag>.csv
#   data/interim/dspin/summary_<tag>/dspin_unit_tests_<tag>.csv
#   tables/DSPIN/summary_<tag>/DSPIN_*.csv   (display tables read by the notebook)

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1]]
}
tag <- get_arg("--tag", format(Sys.Date(), "%Y-%m-%d"))

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
dspin_dir <- file.path(wkdir, "data/interim/dspin")
out_dir <- file.path(dspin_dir, glue("summary_{tag}"))
tbl_dir <- file.path(wkdir, "tables/DSPIN", glue("summary_{tag}"))
mast_csv <- file.path(wkdir, "data/interim/mast/snRNASeq/mast_signif_hits_2026-08-17_p_adjust_0.1.csv")
for (d in c(out_dir, tbl_dir)) {
  if (dir.exists(d) && length(list.files(d))) stop("Output folder exists and is not empty: ", d, " (use a new --tag)")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)
curated_cats <- c("GOTERM_BP_DIRECT", "GOTERM_MF_DIRECT", "GOTERM_CC_DIRECT", "KEGG_PATHWAY")
p_lead <- 0.0065
sep_lead <- 4L

manifest <- read_csv(file.path(dspin_dir, "full_manifest.csv"), show_col_types = FALSE) %>%
  mutate(unit_id = glue("{tissue}__{supertype_label}") %>% as.character(),
         unit_dir = file.path(dspin_dir, unit_id))
fit_summary <- read_csv(file.path(dspin_dir, "fit_summary.csv"), show_col_types = FALSE)
# join on unit_id: fit_summary.csv has a blank supertype_label for the 5 pilot units
units <- manifest %>% semi_join(fit_summary %>% filter(status == "ok"), by = "unit_id")
log("Units with ok fits: ", nrow(units), " / ", nrow(manifest))

# ── helpers ─────────────────────────────────────────────────────────────────
read_weights <- function(f) {
  w <- read_csv(f, show_col_types = FALSE, col_types = cols(.default = "c"))
  ks <- names(w)[!str_starts(names(w), "Weights_")]
  map_dfr(ks, function(k) {
    tibble(program_id = paste0("P", k), gene = w[[k]], weight = as.numeric(w[[paste0("Weights_", k)]]))
  }) %>% filter(!is.na(gene), nzchar(gene)) %>% arrange(program_id, desc(weight))
}

read_david_list <- function(f) {
  x <- str_split(read_lines(f), "\t")
  x <- x[lengths(x) >= 13]
  x <- x[map_chr(x, 1) != "Category"]
  if (!length(x)) return(tibble())
  tibble(category = map_chr(x, 1), term = map_chr(x, 2), count = as.integer(map_chr(x, 3)),
         pvalue = as.numeric(map_chr(x, 5)), fold = as.numeric(map_chr(x, 10)),
         benjamini = as.numeric(map_chr(x, 12))) %>%
    distinct(category, term, .keep_all = TRUE)
}

clean_term <- function(x) x %>% str_remove("^GO:\\d+~") %>% str_remove("^mmu\\d+:")

# exact permutation over animal labels (shared across programs → keeps correlation)
perm_tests <- function(H, is_w) {
  # H: animals × programs matrix of relative h; is_w: logical (WildR)
  n <- nrow(H); nw <- sum(is_w)
  combos <- combn(n, nw)
  welch_t <- function(idx) {
    a <- H[idx, , drop = FALSE]; b <- H[-idx, , drop = FALSE]
    (colMeans(a) - colMeans(b)) / sqrt(apply(a, 2, var) / nrow(a) + apply(b, 2, var) / nrow(b))
  }
  dmean <- function(idx) colMeans(H[idx, , drop = FALSE]) - colMeans(H[-idx, , drop = FALSE])
  obs_d <- dmean(which(is_w)); obs_t <- welch_t(which(is_w))
  null_d <- apply(combos, 2, dmean)                   # programs × perms
  null_s <- apply(combos, 2, function(i) sum(welch_t(i)^2, na.rm = TRUE))
  list(
    p_prog = rowMeans(abs(null_d) >= abs(obs_d) - 1e-12),
    t_prog = obs_t,
    p_unit = mean(null_s >= sum(obs_t^2, na.rm = TRUE) - 1e-9),
    stat_unit = sum(obs_t^2, na.rm = TRUE),
    n_perm = ncol(combos)
  )
}

hedges_g <- function(a, b) {
  na <- length(a); nb <- length(b)
  sp <- sqrt(((na - 1) * var(a) + (nb - 1) * var(b)) / (na + nb - 2))
  (mean(a) - mean(b)) / sp * (1 - 3 / (4 * (na + nb) - 9))
}

# ── MAST (supertype level, padj < 0.1) ───────────────────────────────────────
mast <- read_csv(mast_csv, show_col_types = FALSE) %>%
  filter(cell_type_level == "supertype") %>%
  transmute(tissue, supertype_label = cell_type, gene = str_remove(primerid, "-GRCm39$"), coef, p_adjust)
log("MAST supertype hits: ", nrow(mast))

# ── per-unit loop ────────────────────────────────────────────────────────────
prog_rows <- list(); unit_rows <- list(); animal_rows <- list(); david_rows <- list()
for (i in seq_len(nrow(units))) {
  u <- units[i, ]
  rel <- read_csv(file.path(u$unit_dir, "relative_responses.csv"), show_col_types = FALSE)
  sm <- read_csv(file.path(u$unit_dir, "sample_meta.csv"), show_col_types = FALSE)
  rel <- rel %>% left_join(sm %>% select(sample_id, microbiome, n_cells), by = "sample_id")
  animal_rows[[i]] <- rel %>% mutate(unit_id = u$unit_id)

  H <- rel %>% select(sample_id, program_id, relative_h) %>%
    pivot_wider(names_from = program_id, values_from = relative_h)
  is_w <- sm$microbiome[match(H$sample_id, sm$sample_id)] == "WildR"
  Hm <- as.matrix(H[, -1]); rownames(Hm) <- H$sample_id
  pt <- perm_tests(Hm, is_w)

  per_prog <- map_dfr(colnames(Hm), function(p) {
    w <- Hm[is_w, p]; s <- Hm[!is_w, p]; d <- mean(w) - mean(s)
    tibble(program_id = p, mean_h_WildR = mean(w), mean_h_SPF = mean(s), delta_h = d,
           sd_h_WildR = sd(w), sd_h_SPF = sd(s), hedges_g = hedges_g(w, s),
           separation = if (d >= 0) sum(w > max(s)) else sum(w < min(s)),
           n_WildR_same_sign = sum(sign(w) == sign(d)),
           n_animals_WildR = length(w), n_animals_SPF = length(s))
  }) %>% mutate(welch_t = pt$t_prog[program_id], p_perm = pt$p_prog[program_id])

  # top genes by oNMF weight
  wts <- read_weights(file.path(u$unit_dir, "dspin_run/gene_programs_total_20_onmf_20_prior_0_weights.csv"))
  top_genes <- wts %>% group_by(program_id) %>%
    summarise(top_genes = paste(head(gene, 6), collapse = ", "), n_genes = n(), .groups = "drop")

  # DAVID: best curated (GO/KEGG) term per program
  dres <- file.path(u$unit_dir, "dspin_run/gene_programs_clean_for_david_david_results")
  dv <- map_dfr(list.files(dres, pattern = "^list_\\d+\\.csv$", full.names = TRUE), function(f) {
    read_david_list(f) %>% mutate(program_id = paste0("P", str_extract(basename(f), "\\d+")))
  })
  if (nrow(dv)) david_rows[[i]] <- dv %>% mutate(unit_id = u$unit_id, tissue = u$tissue, supertype_label = u$supertype_label)
  best <- if (nrow(dv)) {
    dv %>% filter(category %in% curated_cats, count >= 3) %>%
      arrange(program_id, benjamini, pvalue) %>% group_by(program_id) %>% slice_head(n = 1) %>% ungroup() %>%
      transmute(program_id, go_term = clean_term(term), go_category = str_remove(category, "_DIRECT|_PATHWAY"),
                go_benjamini = benjamini, go_pvalue = pvalue, go_count = count)
  } else tibble(program_id = character())

  # MAST concordance on all member genes
  mast_u <- mast %>% filter(tissue == u$tissue, supertype_label == u$supertype_label)
  mc <- wts %>% inner_join(mast_u, by = "gene") %>% group_by(program_id) %>%
    summarise(n_mast_hits = n(), mast_mean_coef = mean(coef),
              mast_genes = paste(head(gene[order(p_adjust)], 5), collapse = ", "), .groups = "drop")

  prog_rows[[i]] <- per_prog %>%
    left_join(top_genes, by = "program_id") %>% left_join(best, by = "program_id") %>% left_join(mc, by = "program_id") %>%
    mutate(unit_id = u$unit_id, tissue = u$tissue, supertype_label = u$supertype_label,
           supertype_name = u$supertype_name, n_mast_hits = coalesce(n_mast_hits, 0L),
           mast_sign_agrees = if_else(n_mast_hits > 0, sign(mast_mean_coef) == sign(delta_h), NA))

  abs_h <- rel %>% group_by(microbiome) %>% summarise(m = mean(abs(relative_h)), .groups = "drop") %>% deframe()
  n_spf <- sum(!is_w)
  unit_rows[[i]] <- tibble(
    unit_id = u$unit_id, tissue = u$tissue, supertype_label = u$supertype_label, supertype_name = u$supertype_name,
    n_cells = u$n_cells, n_SPF = u$n_SPF, n_WildR = u$n_WildR,
    min_cells_per_animal = min(sm$n_cells), median_cells_per_animal = median(sm$n_cells),
    n_mast_genes = u$n_mast_genes, augur_auc = u$augur_auc,
    mean_abs_h_WildR = abs_h[["WildR"]], mean_abs_h_SPF = abs_h[["SPF"]],
    abs_h_ratio = abs_h[["WildR"]] / abs_h[["SPF"]],
    mean_abs_h_SPF_loo = abs_h[["SPF"]] * n_spf / (n_spf - 1),
    abs_h_ratio_loo = abs_h[["WildR"]] / (abs_h[["SPF"]] * n_spf / (n_spf - 1)),
    rms_hedges_g = sqrt(mean(per_prog$hedges_g^2)), omnibus_stat = pt$stat_unit, p_unit = pt$p_unit,
    n_perm = pt$n_perm, min_p_program = min(per_prog$p_perm)
  )
  if (i %% 20 == 0) log("  processed ", i, " / ", nrow(units))
}

prog <- bind_rows(prog_rows) %>%
  mutate(q_perm = p.adjust(p_perm, "BH"),
         lead = p_perm <= p_lead & separation >= sep_lead,
         direction = if_else(delta_h > 0, "up in WildR", "down in WildR"),
         unit_short = paste(tissue, str_remove(supertype_name, "^\\d+ "), sep = " · "),
         program_label = glue("{program_id} · {word(top_genes, 1, 3, sep = fixed(', '))}") %>% as.character()) %>%
  relocate(unit_id, tissue, supertype_label, supertype_name, unit_short, program_id, program_label)
unit_t <- bind_rows(unit_rows) %>%
  mutate(q_unit = p.adjust(p_unit, "BH")) %>%
  left_join(prog %>% group_by(unit_id) %>% summarise(n_leads = sum(lead), .groups = "drop"), by = "unit_id") %>%
  arrange(p_unit, desc(rms_hedges_g))
animals <- bind_rows(animal_rows)
david_all <- bind_rows(david_rows)

# ── calibration: observed vs expected small-p counts (exact null of the discrete test) ──
# Under H0 each program's p is uniform over the 924-permutation support; for a symmetric
# balanced design P(p ≤ k/924) = k/924 (ties aside). Expected counts below use that.
cal <- tibble(threshold = c(2, 6, 12, 46, 92) / 924) %>%
  mutate(level = "program", n_tests = nrow(prog),
         observed = map_int(threshold, ~ sum(prog$p_perm <= .x + 1e-12)),
         expected = threshold * n_tests) %>%
  bind_rows(tibble(threshold = c(0.01, 0.05, 0.1)) %>%
              mutate(level = "unit", n_tests = nrow(unit_t),
                     observed = map_int(threshold, ~ sum(unit_t$p_unit <= .x + 1e-12)),
                     expected = threshold * n_tests)) %>%
  mutate(ratio = observed / expected)

# global sign test: are WildR |h| larger than SPF |h| across units?
noise_test <- wilcox.test(unit_t$mean_abs_h_WildR, unit_t$mean_abs_h_SPF, paired = TRUE)
noise_test_loo <- wilcox.test(unit_t$mean_abs_h_WildR, unit_t$mean_abs_h_SPF_loo, paired = TRUE)
noise_cor <- cor.test(log10(unit_t$median_cells_per_animal), unit_t$mean_abs_h_WildR, method = "spearman", exact = FALSE)
global <- tibble(
  n_units = nrow(unit_t), n_programs = nrow(prog),
  n_q_prog_lt_0.1 = sum(prog$q_perm < 0.1), n_q_unit_lt_0.1 = sum(unit_t$q_unit < 0.1),
  n_unit_p_lt_0.05 = sum(unit_t$p_unit < 0.05), n_leads = sum(prog$lead),
  n_lead_units = n_distinct(prog$unit_id[prog$lead]),
  median_abs_h_ratio = median(unit_t$abs_h_ratio),
  frac_units_WildR_gt_SPF = mean(unit_t$mean_abs_h_WildR > unit_t$mean_abs_h_SPF),
  paired_wilcox_p_absh = noise_test$p.value,
  median_abs_h_ratio_loo = median(unit_t$abs_h_ratio_loo),
  frac_units_WildR_gt_SPF_loo = mean(unit_t$abs_h_ratio_loo > 1),
  paired_wilcox_p_absh_loo = noise_test_loo$p.value,
  spearman_rho_absh_vs_cells = unname(noise_cor$estimate), spearman_p_absh_vs_cells = noise_cor$p.value,
  n_programs_with_go = sum(!is.na(prog$go_term) & prog$go_benjamini < 0.05),
  n_leads_mast_agree = sum(prog$lead & prog$mast_sign_agrees %in% TRUE),
  n_leads_with_mast = sum(prog$lead & prog$n_mast_hits > 0)
)

# ── write ────────────────────────────────────────────────────────────────────
write_csv(prog, file.path(out_dir, glue("dspin_program_tests_{tag}.csv")))
write_csv(unit_t, file.path(out_dir, glue("dspin_unit_tests_{tag}.csv")))
write_csv(animals, file.path(out_dir, glue("dspin_animal_relative_h_{tag}.csv")))
write_csv(david_all, file.path(out_dir, glue("dspin_david_all_terms_{tag}.csv")))

write_csv(prog, file.path(tbl_dir, "DSPIN_program_tests.csv"))
write_csv(unit_t, file.path(tbl_dir, "DSPIN_unit_tests.csv"))
write_csv(prog %>% filter(lead) %>% arrange(p_perm, desc(abs(hedges_g))), file.path(tbl_dir, "DSPIN_program_leads.csv"))
write_csv(cal, file.path(tbl_dir, "DSPIN_calibration.csv"))
write_csv(global, file.path(tbl_dir, "DSPIN_global_summary.csv"))

log("Programs: ", nrow(prog), " | leads: ", global$n_leads, " in ", global$n_lead_units, " units",
    " | program q<0.1: ", global$n_q_prog_lt_0.1, " | unit q<0.1: ", global$n_q_unit_lt_0.1)
print(cal)
print(glimpse(global))
log("Done → ", out_dir, " and ", tbl_dir)
