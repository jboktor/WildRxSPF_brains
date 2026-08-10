#!/usr/bin/env Rscript
# Aggregate completed tissue×class pySCENIC results (skip unfinished splits).
# Writes tables + figures under data/processed/scenic_by_class and figures/SCENIC/snRNASeq_by_class/_aggregate.

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(jsonlite)
  library(ggrepel)
  library(scales)
})

theme_set(theme_bw(base_size = 12))
wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
by_class <- file.path(wkdir, "data/interim/scenic/by_class")
manifest_path <- file.path(wkdir, "data/input/scenic/snRNASeq/manifest_class.csv")
out_proc <- file.path(wkdir, "data/processed/scenic_by_class")
out_fig <- file.path(wkdir, "figures/SCENIC/snRNASeq_by_class/_aggregate")
dir.create(out_proc, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", paste0(...), "\n")
}

tf_from_regulon <- function(x) {
  x %>%
    str_replace("\\.pos$", "") %>%
    str_replace("\\.neg$", "") %>%
    str_replace("\\(\\+\\)$", "") %>%
    str_replace("\\(-\\)$", "")
}

manifest <- read_csv(manifest_path, show_col_types = FALSE)
summary_files <- list.files(by_class, pattern = "^run_summary\\.json$", recursive = TRUE, full.names = TRUE)
diff_files <- list.files(by_class, pattern = "^diff_regulons_WildR_vs_SPF\\.csv$", recursive = TRUE, full.names = TRUE)

log_msg("Found ", length(summary_files), " summaries, ", length(diff_files), " diff tables")

summaries <- map_dfr(summary_files, function(p) {
  d <- fromJSON(p)
  tibble(
    stem = d$stem,
    tissue = d$tissue,
    class_label = d$class_label,
    class_name = d$class_name,
    n_cells = d$n_cells,
    n_regulons = d$n_regulons,
    n_sig_padj05 = d$n_sig_padj05,
    summary_path = p
  )
})

# Join expected splits; mark missing
status <- manifest %>%
  transmute(
    stem = str_replace(loom, "\\.loom$", ""),
    tissue, class_label, class_name,
    n_cells_manifest = n_cells,
    n_SPF, n_WildR
  ) %>%
  left_join(summaries, by = c("stem", "tissue", "class_label", "class_name")) %>%
  mutate(
    status = if_else(!is.na(n_regulons), "completed", "pending"),
    n_cells = coalesce(n_cells, n_cells_manifest)
  )

write_csv(status, file.path(out_proc, "run_status.csv"))
log_msg(
  "Completed ", sum(status$status == "completed"), " / ", nrow(status),
  "; pending: ", paste(status$stem[status$status == "pending"], collapse = ", ")
)

diff_all <- map_dfr(diff_files, read_csv, show_col_types = FALSE) %>%
  mutate(
    tf = tf_from_regulon(regulon),
    direction = case_when(
      delta > 0 ~ "WildR_up",
      delta < 0 ~ "SPF_up",
      TRUE ~ "zero"
    ),
    hit = !is.na(padj) & padj < 0.05
  )

write_csv(diff_all, file.path(out_proc, "diff_regulons_all_splits.csv"))

hits <- diff_all %>% filter(hit)
write_csv(hits, file.path(out_proc, glue("diff_regulons_padj05_{Sys.Date()}.csv")))

# Per-split hit counts
hit_counts <- diff_all %>%
  group_by(tissue, class_label, class_name) %>%
  summarise(
    n_regulons_tested = n(),
    n_sig = sum(hit),
    n_WildR_up = sum(hit & direction == "WildR_up"),
    n_SPF_up = sum(hit & direction == "SPF_up"),
    min_padj = suppressWarnings(min(padj, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(
    status %>% select(tissue, class_label, n_cells, status),
    by = c("tissue", "class_label")
  ) %>%
  arrange(tissue, desc(n_sig))

write_csv(hit_counts, file.path(out_proc, "hit_counts_by_split.csv"))

# TF recurrence across completed splits (hit lists only — not raw AUC)
tf_recurrence <- hits %>%
  distinct(tissue, class_label, class_name, tf, direction) %>%
  group_by(tf) %>%
  summarise(
    n_splits = n(),
    n_tissues = n_distinct(tissue),
    n_WildR_up = sum(direction == "WildR_up"),
    n_SPF_up = sum(direction == "SPF_up"),
    classes = paste(sort(unique(paste(tissue, class_name, sep = ":"))), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_splits), tf)

write_csv(tf_recurrence, file.path(out_proc, "tf_hit_recurrence.csv"))

# Top hits overall (by padj, then |delta|)
top_hits <- hits %>%
  arrange(padj, desc(abs(delta))) %>%
  select(tissue, class_name, class_label, regulon, tf, n_SPF, n_WildR,
         mean_SPF, mean_WildR, delta, p, padj, direction) %>%
  mutate(across(c(mean_SPF, mean_WildR, delta), ~ round(.x, 5)),
         across(c(p, padj), ~ signif(.x, 3)))

write_csv(top_hits, file.path(out_proc, "top_hits_padj05.csv"))

# ---- Figures ----

# 1) n_sig bars by class × tissue
p_nsig <- hit_counts %>%
  mutate(
    label = glue("{class_name}\n(n={comma(n_cells)})"),
    label = fct_reorder(label, n_sig)
  ) %>%
  ggplot(aes(x = label, y = n_sig, fill = tissue)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = n_sig), hjust = -0.15, size = 3) +
  coord_flip() +
  facet_wrap(~tissue, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c(AMY = "#E07A3D", HYP = "#3D6B8E"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Significant WildR vs SPF regulons per class (padj < 0.05)",
    subtitle = glue(
      "Completed class splits only (pending: {paste(status$stem[status$status=='pending'], collapse=', ')})"
    ),
    x = NULL, y = "# regulons (BH padj < 0.05)"
  )
ggsave(file.path(out_fig, "n_sig_by_class.png"), p_nsig, width = 10, height = 12, dpi = 160)

# 2) Direction stacked bars
p_dir <- hit_counts %>%
  filter(n_sig > 0) %>%
  pivot_longer(c(n_WildR_up, n_SPF_up), names_to = "dir", values_to = "n") %>%
  mutate(
    dir = recode(dir, n_WildR_up = "WildR higher", n_SPF_up = "SPF higher"),
    class_name = fct_reorder(class_name, n_sig)
  ) %>%
  ggplot(aes(x = class_name, y = n, fill = dir)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~tissue, scales = "free_y") +
  scale_fill_manual(values = c("WildR higher" = "#C44E52", "SPF higher" = "#4C72B0")) +
  labs(
    title = "Direction of significant regulon shifts",
    x = NULL, y = "# significant regulons", fill = NULL
  )
ggsave(file.path(out_fig, "sig_direction_by_class.png"), p_dir, width = 11, height = 8, dpi = 160)

# 3) Combined volcano (all completed splits) — density of tests
set.seed(1)
volc_df <- diff_all %>%
  mutate(
    plot_p = pmax(p, 1e-300),
    lab = if_else(
      hit & abs(delta) >= quantile(abs(delta[hit]), 0.95, na.rm = TRUE),
      paste(tf, class_name, sep = "\n"),
      NA_character_
    )
  )

p_volc <- ggplot(volc_df, aes(x = delta, y = -log10(plot_p))) +
  geom_point(aes(color = hit), alpha = 0.35, size = 0.9) +
  geom_text_repel(aes(label = lab), size = 2.6, max.overlaps = 25, lineheight = 0.85) +
  facet_wrap(~tissue) +
  scale_color_manual(values = c("FALSE" = "grey75", "TRUE" = "#D62728"), guide = "none") +
  labs(
    title = "WildR − SPF regulon AUC (all completed class splits)",
    x = "Δ mean AUC (WildR − SPF)", y = "−log10(p)"
  )
ggsave(file.path(out_fig, "volcano_all_splits.png"), p_volc, width = 11, height = 5.5, dpi = 160)

# 4) TF recurrence (TFs significant in ≥2 splits)
tf_plot <- tf_recurrence %>% filter(n_splits >= 2) %>% slice_head(n = 40)
p_tf <- tf_plot %>%
  mutate(tf = fct_reorder(tf, n_splits)) %>%
  ggplot(aes(x = tf, y = n_splits, fill = n_tissues)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "#9ECAE1", high = "#08306B", breaks = c(1, 2)) +
  labs(
    title = "TFs with significant regulon shifts in ≥2 class splits",
    subtitle = "Hit-list recurrence (regulons inferred separately per split — compare TF identity, not AUC)",
    x = NULL, y = "# class splits with padj < 0.05",
    fill = "# tissues"
  )
ggsave(file.path(out_fig, "tf_recurrence.png"), p_tf, width = 9, height = 8, dpi = 160)

# 5) Presence heatmap: top recurrent TFs × class
top_tfs <- tf_recurrence %>% filter(n_splits >= 2) %>% slice_head(n = 30) %>% pull(tf)
if (length(top_tfs)) {
  heat_df <- hits %>%
    filter(tf %in% top_tfs) %>%
    group_by(tf, tissue, class_name) %>%
    summarise(
      signed = sign(mean(delta)),
      .groups = "drop"
    ) %>%
    mutate(
      split = paste(tissue, class_name, sep = " | "),
      fill = case_when(
        signed > 0 ~ "WildR up",
        signed < 0 ~ "SPF up",
        TRUE ~ "zero"
      )
    )

  # order splits by tissue then n_sig
  split_order <- hit_counts %>%
    arrange(tissue, desc(n_sig)) %>%
    mutate(split = paste(tissue, class_name, sep = " | ")) %>%
    pull(split) %>%
    unique()

  p_heat <- heat_df %>%
    mutate(
      tf = factor(tf, levels = rev(top_tfs)),
      split = factor(split, levels = split_order)
    ) %>%
    ggplot(aes(x = split, y = tf, fill = fill)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_manual(values = c("WildR up" = "#C44E52", "SPF up" = "#4C72B0", "zero" = "grey90")) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
      panel.grid = element_blank()
    ) +
    labs(
      title = "Recurrent significant TFs across class splits",
      subtitle = "Tile = padj < 0.05 in that split; color = mean ΔAUC direction",
      x = NULL, y = NULL, fill = NULL
    )
  ggsave(file.path(out_fig, "tf_hit_heatmap.png"), p_heat, width = 14, height = 8, dpi = 160)
}

# 6) Completeness / run status strip
p_status <- status %>%
  mutate(
    label = glue("{tissue} | {class_name}"),
    label = fct_reorder(label, n_cells)
  ) %>%
  ggplot(aes(x = label, y = n_cells, fill = status)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c(completed = "#2CA02C", pending = "#FF7F0E")) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Class-split SCENIC run status",
    subtitle = "Aggregation excludes pending splits",
    x = NULL, y = "# cells in loom", fill = NULL
  )
ggsave(file.path(out_fig, "run_status.png"), p_status, width = 10, height = 10, dpi = 160)

# Conclusions text
n_done <- sum(status$status == "completed")
n_pend <- sum(status$status == "pending")
n_hit_tests <- nrow(hits)
n_splits_with_hits <- sum(hit_counts$n_sig > 0, na.rm = TRUE)
top_tf_line <- tf_recurrence %>%
  slice_head(n = 12) %>%
  mutate(s = glue("{tf} ({n_splits} splits)")) %>%
  pull(s) %>%
  paste(collapse = ", ")

concl <- c(
  glue("Class-split pySCENIC aggregation ({Sys.Date()})"),
  glue("Completed: {n_done}/{nrow(status)} tissue×class GRNs (pending: {paste(status$stem[status$status=='pending'], collapse=', ')})."),
  glue("Significant regulon tests (padj<0.05): {n_hit_tests} across {n_splits_with_hits} splits."),
  glue("Most recurrent significant TFs: {top_tf_line}."),
  "Note: regulons are inferred per split — cross-class comparisons use TF hit lists, not raw AUCell values.",
  "Input matrix: Seurat RNA counts (raw UMIs)."
)
writeLines(concl, file.path(out_proc, "conclusions.txt"))
write_json(
  list(
    date = as.character(Sys.Date()),
    n_completed = n_done,
    n_pending = n_pend,
    pending_stems = status$stem[status$status == "pending"],
    n_sig_tests = n_hit_tests,
    n_splits_with_hits = n_splits_with_hits,
    top_recurrent_tfs = head(tf_recurrence$tf, 20)
  ),
  file.path(out_proc, "aggregate_summary.json"),
  auto_unbox = TRUE, pretty = TRUE
)

log_msg("Wrote tables → ", out_proc)
log_msg("Wrote figures → ", out_fig)
print(hit_counts, n = 50)
print(head(tf_recurrence, 20))
