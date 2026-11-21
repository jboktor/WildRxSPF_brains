# Animal-level (pseudobulk) WildR vs SPF test for pySCENIC AUCell matrices.
# Shared by run_scenic_one_split.R (per-split, at run time) and
# test_scenic_regulons_animal_level.R (post hoc across splits), so both use
# the identical test.
#
# Unit of replication = animal (6 WildR vs 6 SPF per tissue). Per animal ×
# regulon: mean AUCell over nuclei; animals need >= SCENIC_MIN_CELLS nuclei;
# a split is tested only with >= SCENIC_MIN_ANIMALS animals per group.
# Test: pooled-variance t on animal means (== lm(mean ~ group)), plus Hedges g,
# Cliff's delta, exact Wilcoxon and complete-separation count.

SCENIC_MIN_CELLS <- 10L
SCENIC_MIN_ANIMALS <- 3L

tf_from_regulon <- function(x) sub("(\\.pos|\\.neg|\\(\\+\\)|\\(-\\))$", "", x)
animal_from_sample <- function(s) gsub("_(AMY|HYP)$", "", as.character(s))

cliffs_delta <- function(x, y) {
  d <- outer(x, y, "-")
  (sum(d > 0) - sum(d < 0)) / length(d)
}

# Row-wise two-group statistics on a regulon × animal matrix of means.
row_tests <- function(m, grp) {
  x <- m[, grp == "WildR", drop = FALSE]
  y <- m[, grp == "SPF", drop = FALSE]
  nx <- ncol(x); ny <- ncol(y)
  mx <- rowMeans(x); my <- rowMeans(y)
  vx <- apply(x, 1, var); vy <- apply(y, 1, var)
  df <- nx + ny - 2
  sp <- sqrt(((nx - 1) * vx + (ny - 1) * vy) / df)
  se <- sp * sqrt(1 / nx + 1 / ny)
  tstat <- (mx - my) / se
  tq <- qt(0.975, df)
  j <- 1 - 3 / (4 * (nx + ny) - 9)
  idx <- seq_len(nrow(m))
  tibble::tibble(
    regulon = rownames(m),
    n_animals_WildR = nx, n_animals_SPF = ny,
    animal_mean_WildR = mx, animal_mean_SPF = my,
    delta_animal = mx - my,
    ci_low = (mx - my) - tq * se, ci_high = (mx - my) + tq * se,
    p_animal = ifelse(se > 0, 2 * pt(-abs(tstat), df), NA_real_),
    hedges_g = ifelse(sp > 0, j * (mx - my) / sp, NA_real_),
    cliffs_delta = vapply(idx, function(i) cliffs_delta(x[i, ], y[i, ]), numeric(1)),
    p_wilcox_animal = vapply(idx, function(i) suppressWarnings(tryCatch(
      wilcox.test(x[i, ], y[i, ], exact = TRUE)$p.value, error = function(e) NA_real_)), numeric(1)),
    animals_separated = vapply(idx, function(i)
      max(sum(x[i, ] > max(y[i, ])), sum(x[i, ] < min(y[i, ]))), numeric(1))
  )
}

#' Per-animal mean AUCell.
#' @param auc data.frame/data.table: cell_id + one column per regulon
#' @param meta data.frame with cell_id, sample, microbiome
#' @return list(per_animal = long table, wide = animal × regulon data.table, regulons)
animal_means <- function(auc, meta) {
  dt <- data.table::as.data.table(auc)
  data.table::setnames(dt, 1, "cell_id")
  dt[, cell_id := as.character(cell_id)]
  md <- data.table::as.data.table(meta)[, .(cell_id = as.character(cell_id), sample, microbiome)]
  dt <- merge(md, dt, by = "cell_id")
  dt[, animal := animal_from_sample(sample)]
  reg_cols <- setdiff(names(dt), c("cell_id", "sample", "microbiome", "animal"))
  wide <- dt[, c(list(n_cells = .N), lapply(.SD, mean)), by = .(animal, microbiome), .SDcols = reg_cols]
  long <- data.table::melt(wide, id.vars = c("animal", "microbiome", "n_cells"),
                           variable.name = "regulon", value.name = "mean_auc")
  list(per_animal = long, wide = wide, regulons = reg_cols, n_cells_matched = nrow(dt))
}

#' Animal-level test for one split.
#' @return list(res = tibble or NULL if underpowered, n_WildR, n_SPF, reason)
animal_level_test <- function(am, min_cells = SCENIC_MIN_CELLS, min_animals = SCENIC_MIN_ANIMALS) {
  keep <- am$wide[n_cells >= min_cells]
  n_w <- sum(keep$microbiome == "WildR"); n_s <- sum(keep$microbiome == "SPF")
  out <- list(res = NULL, n_WildR = n_w, n_SPF = n_s,
              reason = sprintf("animals with >= %d nuclei: WildR %d, SPF %d", min_cells, n_w, n_s))
  if (n_w < min_animals || n_s < min_animals) return(out)
  m <- t(as.matrix(keep[, am$regulons, with = FALSE]))
  colnames(m) <- keep$animal
  out$res <- row_tests(m, keep$microbiome) %>%
    dplyr::mutate(fdr_split = p.adjust(p_animal, "BH"),
                  tf = tf_from_regulon(regulon),
                  direction = dplyr::if_else(delta_animal > 0, "WildR higher", "SPF higher"))
  out
}
