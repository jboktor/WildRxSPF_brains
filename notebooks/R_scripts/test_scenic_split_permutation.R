#!/usr/bin/env Rscript
# Animal-label permutation tests for an excess of small animal-level p-values
# in pySCENIC regulon tests, per split and per level × tissue.
#
# Why: counting regulons at p < 0.05 against a binomial(n, 0.05) null assumes
# regulons are independent. They are not (overlapping targets, shared cell-state
# axes), so a split-level binomial test is anti-conservative. Permuting the
# microbiome label over animals keeps the regulon correlation and is the exact
# null for the 6 vs 6 design.
#
# Method (reads saved per-animal means only; no GRN / AUCell recompute):
#   statistic = number of regulons with pooled-variance t p < 0.05 (same test
#   as scenic_animal_level_helpers.R::row_tests), animals with >= 10 nuclei.
#   - per split: enumerate every relabelling of the split's eligible animals
#     with the observed WildR count (<= 924); p_perm = share of relabellings
#     with a count >= observed (observed labelling included).
#   - per level × tissue: one relabelling of the tissue's 12 animals is applied
#     to every split jointly (animals are shared across splits), so the summed
#     count respects the dependence between splits. A split whose relabelled
#     groups fall below 3 animals contributes its null expectation 0.05 × n.
#
# Inputs:  data/interim/scenic/animal_level/regulon_auc_per_animal{,_subclass,_supertype}_<date>.csv
# Outputs (date-stamped, today):
#   tables/SCENIC/SCENIC_split_permutation_tests_<date>.csv
#   tables/SCENIC/SCENIC_level_permutation_tests_<date>.csv
#   data/interim/scenic/animal_level/level_permutation_null_<date>.csv
#
# Run:
#   /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript \
#     notebooks/R_scripts/test_scenic_split_permutation.R [animal_level_date]

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(glue)
})

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
source(file.path(wkdir, "notebooks/R_scripts/scenic_animal_level_helpers.R"))
al_dir <- file.path(wkdir, "data/interim/scenic/animal_level")
tbl_dir <- file.path(wkdir, "tables/SCENIC")
args <- commandArgs(trailingOnly = TRUE)
al_date <- if (length(args)) args[1] else
  max(str_extract(list.files(al_dir, "^regulon_auc_per_animal_\\d{4}-\\d{2}-\\d{2}\\.csv$"), "\\d{4}-\\d{2}-\\d{2}"))
run_date <- Sys.Date()
log_msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", paste0(...), "\n")
log_msg("Per-animal means dated ", al_date)

# Count of rows with pooled-t p < 0.05 for every labelling column of L (animals × perms, 1 = WildR)
nominal_counts <- function(m, L, alpha = 0.05) {
  nx <- colSums(L); ny <- nrow(L) - nx
  sx <- m %*% L; sy <- rowSums(m) - sx
  qx <- (m^2) %*% L; qy <- rowSums(m^2) - qx
  mx <- sweep(sx, 2, nx, "/"); my <- sweep(sy, 2, ny, "/")
  ssx <- qx - sweep(sx^2, 2, nx, "/"); ssy <- qy - sweep(sy^2, 2, ny, "/")
  df <- nx + ny - 2
  sp2 <- sweep(ssx + ssy, 2, df, "/")
  se <- sqrt(sweep(sp2, 2, 1 / nx + 1 / ny, "*"))
  tcrit <- qt(1 - alpha / 2, df)
  hit <- sweep(abs(mx - my) / se, 2, tcrit, ">")
  hit[!is.finite(se) | se <= 1e-15] <- FALSE
  colSums(hit, na.rm = TRUE)
}

levels_all <- c(class = "", subclass = "_subclass", supertype = "_supertype")
split_rows <- list(); null_rows <- list()

for (lv in names(levels_all)) {
  f <- file.path(al_dir, glue("regulon_auc_per_animal{levels_all[[lv]]}_{al_date}.csv"))
  if (!file.exists(f)) { log_msg("skip ", lv, ": no ", basename(f)); next }
  pa <- fread(f)[n_cells >= SCENIC_MIN_CELLS]
  log_msg(lv, ": ", uniqueN(pa$stem), " splits")
  for (tis in sort(unique(pa$tissue))) {
    pt <- pa[tissue == tis]
    animals <- sort(unique(pt$animal))
    truth <- setNames(unique(pt[, .(animal, microbiome)])$microbiome, unique(pt[, .(animal, microbiome)])$animal)[animals]
    n_w <- sum(truth == "WildR")
    combos <- combn(length(animals), n_w)
    Ljoint <- matrix(0L, length(animals), ncol(combos), dimnames = list(animals, NULL))
    for (k in seq_len(ncol(combos))) Ljoint[combos[, k], k] <- 1L
    obs_col <- which(colSums(abs(Ljoint - (truth == "WildR"))) == 0)
    joint_sum <- numeric(ncol(Ljoint))

    for (st in sort(unique(pt$stem))) {
      d <- pt[stem == st]
      w <- dcast(d, regulon ~ animal, value.var = "mean_auc")
      m <- as.matrix(w[, -1]); rownames(m) <- w$regulon
      m <- m[apply(m, 1, var) > 0, , drop = FALSE]
      grp <- truth[colnames(m)]
      nW <- sum(grp == "WildR"); nS <- sum(grp == "SPF")
      if (nW < SCENIC_MIN_ANIMALS || nS < SCENIC_MIN_ANIMALS) next
      obs <- nominal_counts(m, matrix(as.integer(grp == "WildR")))
      # exact split-level null: all relabellings of this split's animals
      cs <- combn(ncol(m), nW)
      Ls <- matrix(0L, ncol(m), ncol(cs))
      for (k in seq_len(ncol(cs))) Ls[cs[, k], k] <- 1L
      null_s <- nominal_counts(m, Ls)
      # joint relabelling for the level × tissue sum
      Lj <- Ljoint[colnames(m), , drop = FALSE]
      ok <- colSums(Lj) >= SCENIC_MIN_ANIMALS & colSums(1L - Lj) >= SCENIC_MIN_ANIMALS
      cnt <- rep(0.05 * nrow(m), ncol(Lj))
      cnt[ok] <- nominal_counts(m, Lj[, ok, drop = FALSE])
      joint_sum <- joint_sum + cnt
      meta <- d[1]
      split_rows[[length(split_rows) + 1]] <- tibble(
        level = lv, stem = st, tissue = tis, class_name = meta$class_name,
        subclass_name = meta$subclass_name, supertype_name = meta$supertype_name,
        n_animals_WildR = nW, n_animals_SPF = nS, n_tested = nrow(m),
        n_nominal_obs = obs, n_expected = 0.05 * nrow(m),
        perm_mean = mean(null_s), perm_q95 = unname(quantile(null_s, 0.95)),
        n_perm = length(null_s), p_perm = mean(null_s >= obs),
        p_binom = binom.test(obs, nrow(m), 0.05, alternative = "greater")$p.value)
    }
    null_rows[[length(null_rows) + 1]] <- tibble(level = lv, tissue = tis, perm = seq_along(joint_sum),
                                                 observed_labelling = seq_along(joint_sum) == obs_col,
                                                 n_nominal = joint_sum)
    log_msg("  ", tis, ": observed ", round(joint_sum[obs_col]), " vs permutation median ",
            round(median(joint_sum)), " (", ncol(Ljoint), " relabellings)")
  }
}

split_res <- bind_rows(split_rows) %>%
  group_by(level) %>%
  mutate(fdr_perm_level = p.adjust(p_perm, "BH"), fdr_binom_level = p.adjust(p_binom, "BH")) %>%
  ungroup() %>%
  arrange(factor(level, names(levels_all)), p_perm)
null_res <- bind_rows(null_rows)
level_res <- null_res %>%
  group_by(level, tissue) %>%
  summarise(n_nominal_obs = n_nominal[observed_labelling],
            perm_median = median(n_nominal), perm_q05 = quantile(n_nominal, 0.05),
            perm_q95 = quantile(n_nominal, 0.95), n_perm = n(),
            p_perm_excess = mean(n_nominal >= n_nominal[observed_labelling]),
            .groups = "drop")

write_csv(split_res, file.path(tbl_dir, glue("SCENIC_split_permutation_tests_{run_date}.csv")))
write_csv(level_res, file.path(tbl_dir, glue("SCENIC_level_permutation_tests_{run_date}.csv")))
write_csv(null_res, file.path(al_dir, glue("level_permutation_null_{run_date}.csv")))
print(level_res)
print(split_res %>% group_by(level) %>%
        summarise(splits = n(), p_perm_lt05 = sum(p_perm < 0.05), fdr_perm_lt10 = sum(fdr_perm_level < 0.10),
                  fdr_binom_lt10 = sum(fdr_binom_level < 0.10)))
log_msg("Done")
