#!/usr/bin/env Rscript
# Rank DSPIN (supertype_label × tissue) units from candidate_cellcounts + MAST + Augur.
# Writes:
#   pilot_candidates_draft.csv  — top 15
#   pilot_manifest.csv          — top ≤5 (interactive)
#   full_manifest.csv           — ALL hard_pass units (run_tier=full)
#   ranked_all.csv              — full ranking table
#
# Path key: tissue__{supertype_label} (never slugify supertype_name — "/" hazard).
#
# Usage:
#   Rscript notebooks/R_scripts/rank_dspin_pilots.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
})

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
data_path <- file.path(wkdir, "data/interim")
out_dir <- file.path(data_path, "dspin")
cand_path <- file.path(out_dir, "candidate_cellcounts.csv")
mast_path <- file.path(data_path, "mast/snRNASeq/supertype_mast_signif_hits_2026-01-04.rds")
augur_path <- file.path(
  data_path, "augur/snRNASeq/var_thresh_0.6/supertype_aucs_df_2026-01-10.rds"
)

stopifnot(file.exists(cand_path))
stopifnot(file.exists(mast_path))
stopifnot(file.exists(augur_path))

cand <- read_csv(cand_path, show_col_types = FALSE)
stopifnot("supertype_label" %in% names(cand))

mast <- readRDS(mast_path)
augur <- readRDS(augur_path)

# MAST / Augur keyed by display name; join via cand's name↔label map
mast_n <- mast %>%
  as_tibble() %>%
  filter(p_adjust < 0.1) %>%
  count(tissue, supertype_name, name = "n_mast_genes")

augur_n <- augur %>%
  as_tibble() %>%
  transmute(
    tissue = as.character(tissue),
    supertype_name = as.character(supertype_name),
    augur_auc = as.numeric(auc)
  ) %>%
  distinct(tissue, supertype_name, .keep_all = TRUE)

ranked <- cand %>%
  filter(hard_pass) %>%
  left_join(mast_n, by = c("tissue", "supertype_name")) %>%
  left_join(augur_n, by = c("tissue", "supertype_name")) %>%
  mutate(
    n_mast_genes = replace_na(n_mast_genes, 0L),
    augur_auc = replace_na(augur_auc, 0.5),
    z_auc = if (n() < 2) 0 else as.numeric(scale(augur_auc)),
    z_mast = if (n() < 2) 0 else as.numeric(scale(log1p(n_mast_genes))),
    score = 0.5 * replace_na(z_auc, 0) + 0.5 * replace_na(z_mast, 0)
  ) %>%
  arrange(desc(score), desc(min_cells_per_sample_id)) %>%
  mutate(
    rank = row_number(),
    run_tier = case_when(
      rank <= 5 ~ "interactive",
      rank <= 15 ~ "optional_sbatch",
      TRUE ~ "pool"
    ),
    unit_id = paste(tissue, supertype_label, sep = "__"),
    unit_path = file.path(out_dir, unit_id),
    # aliases kept for older readers
    slice_path = unit_path,
    h5ad = file.path(unit_path, "raw_counts.h5ad"),
    filtered = file.path(unit_path, "filtered.h5ad"),
    n_genes = NA_integer_,
    fit_status = "pending",
    prep_status = "pending"
  )

write_csv(ranked, file.path(out_dir, "ranked_all.csv"))

draft <- ranked %>% filter(rank <= 15)
write_csv(draft, file.path(out_dir, "pilot_candidates_draft.csv"))

pilot <- ranked %>% filter(run_tier == "interactive")
write_csv(pilot, file.path(out_dir, "pilot_manifest.csv"))

full <- ranked %>% mutate(run_tier = "full")
write_csv(full, file.path(out_dir, "full_manifest.csv"))

cat(
  "Wrote ranked_all (", nrow(ranked), "), draft (", nrow(draft),
  "), pilot (", nrow(pilot), "), full_manifest (", nrow(full), ") → ",
  out_dir, "\n",
  sep = ""
)
print(full %>% count(tissue))
print(summary(full$n_cells))
print(head(full %>% select(tissue, supertype_label, unit_id, supertype_name), 3))
