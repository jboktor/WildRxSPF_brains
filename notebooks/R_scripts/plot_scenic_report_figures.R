#!/usr/bin/env Rscript
# Report figures + display tables for 0X_SCENIC_pySCENIC.qmd.
# Reads saved outputs only (no GRN / AUCell recompute):
#   - data/processed/scenic_by_class/*            (cell-level aggregation)
#   - tables/SCENIC/SCENIC_animal_level_*_<date>  (test_scenic_regulons_animal_level.R)
#   - data/interim/scenic/animal_level/*          (per-animal regulon means)
#   - data/processed/scenic/regulon_rss_*.csv     (pilot RSS)
#   - tables/SCENIC/SCENIC_{split,level}_permutation_tests_<date>.csv (test_scenic_split_permutation.R)
# Writes figures/SCENIC/report/SCENIC_*.png and tables/SCENIC/SCENIC_*_display_<date>.csv.
#
# Style follows mapmycells_qc_functions.R: 12 in wide, 200 dpi, theme_minimal
# base 13, WildR orange / SPF blue, diverging blue-gray-orange for effects.
#
# Run:
#   /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript \
#     notebooks/R_scripts/plot_scenic_report_figures.R [animal_level_date]

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(scales)
  library(patchwork)
  library(ggrepel)
})

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
args <- commandArgs(trailingOnly = TRUE)
tbl_dir <- file.path(wkdir, "tables/SCENIC")
fig_dir <- file.path(wkdir, "figures/SCENIC/report")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

latest <- function(dir, pattern) {
  hits <- sort(list.files(dir, pattern = pattern, full.names = TRUE))
  if (!length(hits)) stop("No file matching ", pattern, " in ", dir)
  tail(hits, 1)
}
al_date <- if (length(args)) args[1] else
  str_extract(basename(latest(tbl_dir, "^SCENIC_animal_level_tests_\\d{4}-\\d{2}-\\d{2}\\.csv$")), "\\d{4}-\\d{2}-\\d{2}")

log_msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", paste0(...), "\n")
log_msg("Animal-level results dated ", al_date)

# ---- Style ----
group_cols <- c(WildR = "#E8590C", SPF = "#2E6FA8")
tissue_cols <- c(AMY = "#27455C", HYP = "#FF9240")
div_cols <- c(low = "#2E6FA8", mid = "#F0EFEC", high = "#E8590C")
fig_w <- 12
fig_dpi <- 200

theme_scenic <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#ECECEC", linewidth = 0.3),
      axis.text = element_text(colour = "#555555"),
      axis.title = element_text(colour = "#555555"),
      plot.title = element_text(face = "bold", colour = "#0A0A0A"),
      plot.subtitle = element_text(colour = "#555555"),
      strip.text = element_text(face = "bold", colour = "#0A0A0A", hjust = 0),
      legend.position = "top",
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.tag = element_text(face = "bold", size = base_size + 2),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}
save_fig <- function(p, file, h) {
  ggsave(file.path(fig_dir, file), p, width = fig_w, height = h, dpi = fig_dpi, bg = "white")
  log_msg("  wrote ", file)
}
short_class <- function(x) str_remove(x, "^\\d+\\s+")

# ---- Inputs ----
tests <- read_csv(file.path(tbl_dir, glue("SCENIC_animal_level_tests_{al_date}.csv")), show_col_types = FALSE)
split_sum <- read_csv(file.path(tbl_dir, glue("SCENIC_animal_level_split_summary_{al_date}.csv")), show_col_types = FALSE)
consensus <- read_csv(file.path(tbl_dir, glue("SCENIC_tf_consensus_animal_level_{al_date}.csv")), show_col_types = FALSE)
per_animal <- read_csv(file.path(wkdir, glue("data/interim/scenic/animal_level/regulon_auc_per_animal_{al_date}.csv")),
                       show_col_types = FALSE)
cell_all <- read_csv(file.path(wkdir, "data/processed/scenic_by_class/diff_regulons_all_splits.csv"), show_col_types = FALSE)
rss <- read_csv(latest(file.path(wkdir, "data/processed/scenic"), "^regulon_rss_.*\\.csv$"), show_col_types = FALSE)

split_sum <- split_sum %>%
  mutate(class_short = short_class(class_name),
         n_tested = coalesce(n_tested, 0L),
         n_cell_level_hits = coalesce(n_cell_level_hits, 0L))

# ---- Fig 1: coverage per split ----
# order classes by nuclei *within* each tissue; label carries a tissue key so
# facet_grid(free_y) keeps per-tissue ordering, stripped for display
cov <- split_sum %>%
  arrange(tissue, n_cells) %>%
  mutate(label = factor(paste0(tissue, "__", class_short), levels = unique(paste0(tissue, "__", class_short))))
strip_key <- function(x) str_remove(x, "^[A-Z]+__")
cov_long_animals <- cov %>%
  select(tissue, label, WildR = n_animals_WildR_ok, SPF = n_animals_SPF_ok, tested) %>%
  pivot_longer(c(WildR, SPF), names_to = "group", values_to = "n_animals")

pA <- ggplot(cov, aes(n_cells, label)) +
  geom_segment(aes(x = 100, xend = n_cells, yend = label), colour = "#C9C9C9", linewidth = 0.6) +
  geom_point(aes(colour = tissue), size = 2.8) +
  scale_x_log10(labels = label_comma(), breaks = c(100, 1000, 1e4, 5e4)) +
  scale_colour_manual(values = tissue_cols, guide = "none") +
  scale_y_discrete(labels = strip_key) +
  facet_grid(tissue ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Nuclei in class GRN (log scale)", y = NULL, tag = "A") +
  theme_scenic() + theme(strip.text.y = element_text(angle = 0))

pB <- ggplot(cov_long_animals, aes(n_animals, label, colour = group)) +
  geom_vline(xintercept = 3, linetype = "dashed", colour = "#A5A5A5") +
  geom_point(size = 2.6, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = group_cols, name = NULL) +
  scale_x_continuous(breaks = 0:6, limits = c(0, 6.3)) +
  facet_grid(tissue ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Animals with ≥10 nuclei", y = NULL, tag = "B") +
  theme_scenic() + theme(axis.text.y = element_blank(), strip.text.y = element_text(angle = 0))

pC <- ggplot(cov, aes(n_regulons, label)) +
  geom_col(aes(fill = tissue), width = 0.7) +
  scale_fill_manual(values = tissue_cols, guide = "none") +
  facet_grid(tissue ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Regulons called", y = NULL, tag = "C") +
  theme_scenic() + theme(axis.text.y = element_blank(), strip.text.y = element_text(angle = 0))

p1 <- (pA | pB | pC) + plot_layout(widths = c(1.3, 1, 0.8), guides = "collect") &
  theme(legend.position = "top")
save_fig(p1, "SCENIC_split_coverage.png", 11)

# ---- Fig 2: cell-level hits are driven by nucleus count, not animals ----
hit_df <- split_sum %>% filter(tested)
pA <- ggplot(hit_df, aes(n_cells, n_cell_level_hits, colour = tissue)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(aes(label = class_short), size = 3.4, colour = "#555555",
                  data = hit_df %>% filter(n_cell_level_hits >= 15), seed = 1,
                  min.segment.length = 0, box.padding = 0.4) +
  scale_x_log10(labels = label_comma()) +
  scale_colour_manual(values = tissue_cols, name = NULL) +
  labs(x = "Nuclei in class GRN (log scale)",
       y = "Regulons with cell-level padj < 0.05",
       title = "Cell-level hits track nuclei", tag = "A") +
  theme_scenic()

ch <- tests %>% filter(cell_level_hit)
pB <- ggplot(ch, aes(p_animal)) +
  geom_histogram(breaks = seq(0, 1, 0.05), fill = "#767676", colour = "white") +
  geom_hline(yintercept = nrow(ch) / 20, linetype = "dashed", colour = "#E8590C") +
  annotate("text", x = 1, y = nrow(ch) / 20, label = "expected if null", hjust = 1, vjust = -0.6,
           colour = "#E8590C", size = 3.6) +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) +
  labs(x = "Animal-level p (6 v 6)", y = "Cell-level hits",
       title = glue("{nrow(ch)} cell-level hits, per animal"), tag = "B") +
  theme_scenic()

# Example: strongest cell-level hit
ex <- cell_all %>% arrange(p) %>% slice(1)
ex_stem <- paste0("snRNASeq_", ex$tissue, "_", ex$class_label)
ex_reg <- str_replace(ex$regulon, "\\.pos$", "(+)")
ex_df <- per_animal %>% filter(stem == ex_stem, regulon == ex_reg)
ex_test <- tests %>% filter(stem == ex_stem, regulon == ex_reg)
pC <- ggplot(ex_df, aes(microbiome, mean_auc, colour = microbiome)) +
  geom_point(aes(size = n_cells), position = position_jitter(width = 0.12, seed = 1), alpha = 0.9) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.45, colour = "#0A0A0A", linewidth = 0.4) +
  scale_colour_manual(values = group_cols, guide = "none") +
  scale_size_area(max_size = 5, name = "Nuclei", labels = label_comma()) +
  labs(x = NULL, y = "Mean AUCell per animal",
       title = glue("{ex_reg} · {ex$tissue} {short_class(ex$class_name)}"),
       subtitle = glue("cell-level p = {signif(ex$p, 2)}  ·  animal-level p = {signif(ex_test$p_animal, 2)}"),
       tag = "C") +
  theme_scenic() + theme(legend.position = "right")

p2 <- (pA | pB | pC) + plot_layout(widths = c(1.15, 1, 0.9))
save_fig(p2, "SCENIC_cell_vs_animal_level.png", 5.6)

# ---- Fig 3: animal-level volcano + QQ ----
vdf <- tests %>%
  mutate(class_short = short_class(class_name),
         flag = case_when(fdr_split < 0.10 ~ "FDR < 0.10 (within class)",
                          p_animal < 0.05 ~ "nominal p < 0.05",
                          TRUE ~ "n.s."),
         flag = factor(flag, c("FDR < 0.10 (within class)", "nominal p < 0.05", "n.s.")))
lab_df <- vdf %>% group_by(tissue) %>% slice_min(p_animal, n = 6) %>% ungroup() %>%
  mutate(lab = glue("{regulon}\n{class_short}"))
pA <- ggplot(vdf, aes(hedges_g, -log10(p_animal))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#A5A5A5") +
  geom_point(aes(colour = flag, size = flag), alpha = 0.75) +
  geom_text_repel(data = lab_df, aes(label = lab), size = 3.1, lineheight = 0.85,
                  colour = "#0A0A0A", min.segment.length = 0, seed = 2, box.padding = 0.5,
                  max.overlaps = Inf) +
  scale_colour_manual(values = c("FDR < 0.10 (within class)" = "#E8590C",
                                 "nominal p < 0.05" = "#27455C", "n.s." = "#D4D4D4"), name = NULL) +
  scale_size_manual(values = c(2.6, 1.6, 0.9), guide = "none") +
  facet_wrap(~tissue) +
  labs(x = "Hedges g on animal means (WildR − SPF)", y = "−log10 animal-level p", tag = "A") +
  theme_scenic()

qq <- vdf %>% group_by(tissue) %>% arrange(p_animal) %>%
  mutate(expected = -log10(ppoints(n())), observed = -log10(p_animal)) %>% ungroup()
pB <- ggplot(qq, aes(expected, observed, colour = tissue)) +
  geom_abline(linetype = "dashed", colour = "#A5A5A5") +
  geom_point(size = 1.1, alpha = 0.8) +
  scale_colour_manual(values = tissue_cols, name = NULL) +
  coord_equal() +
  labs(x = "Expected −log10 p (uniform)", y = "Observed −log10 p", tag = "B") +
  theme_scenic()
p3 <- (pA | pB) + plot_layout(widths = c(2.3, 1))
save_fig(p3, "SCENIC_animal_level_volcano.png", 6.2)

# ---- Fig 4: superplots for the strongest animal-level regulons ----
top_reg <- tests %>%
  filter(n_animals_WildR >= 5, n_animals_SPF >= 5) %>%
  slice_min(p_animal, n = 8) %>%
  mutate(panel = glue("{regulon} · {tissue} {short_class(class_name)}\np = {signif(p_animal, 2)}, q = {signif(fdr_split, 2)}, g = {round(hedges_g, 1)}"),
         panel = fct_inorder(panel))
sp_df <- per_animal %>%
  filter(n_cells >= 10) %>%
  inner_join(top_reg %>% select(stem, regulon, panel), by = c("stem", "regulon"))
p4 <- ggplot(sp_df, aes(microbiome, mean_auc, colour = microbiome)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "#0A0A0A", linewidth = 0.35) +
  geom_point(aes(size = n_cells), position = position_jitter(width = 0.1, seed = 1), alpha = 0.9) +
  facet_wrap(~panel, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = group_cols, guide = "none") +
  scale_size_area(max_size = 4.5, name = "Nuclei per animal", labels = label_comma()) +
  labs(x = NULL, y = "Mean AUCell per animal") +
  theme_scenic(12) + theme(strip.text = element_text(size = 10.5))
save_fig(p4, "SCENIC_top_animal_level_regulons.png", 7.4)

# ---- Fig 5: recurrent TFs — cell-level hits vs animal-level effect ----
cell_rec <- tests %>% filter(cell_level_hit) %>% count(tf, name = "n_cell_hits") %>%
  arrange(desc(n_cell_hits))
top_tfs <- head(cell_rec$tf, 24)
split_levels <- split_sum %>% filter(tested) %>% arrange(tissue, desc(n_cells)) %>%
  pull(class_short) %>% unique()
heat <- tests %>%
  filter(tf %in% top_tfs, str_detect(regulon, "\\(\\+\\)$")) %>%
  mutate(class_short = factor(short_class(class_name), levels = split_levels),
         tf = factor(tf, levels = rev(top_tfs)),
         g = pmax(pmin(hedges_g, 2.5), -2.5))
p5 <- ggplot(heat, aes(class_short, tf)) +
  geom_tile(aes(fill = g), colour = "white", linewidth = 0.6) +
  geom_point(data = heat %>% filter(cell_level_hit), shape = 1, size = 2.3, colour = "#0A0A0A", stroke = 0.6) +
  geom_point(data = heat %>% filter(p_animal < 0.05), shape = 8, size = 2, colour = "#0A0A0A") +
  scale_fill_gradient2(low = div_cols["low"], mid = div_cols["mid"], high = div_cols["high"],
                       limits = c(-2.5, 2.5), name = "Animal-level\nHedges g\n(WildR − SPF)") +
  facet_grid(. ~ tissue, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL,
       caption = "Tile = regulon called in that class GRN (blank = not called). ○ = cell-level padj < 0.05.  ✱ = animal-level p < 0.05.") +
  theme_scenic(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10.5),
        axis.text.y = element_text(face = "italic", colour = "#0A0A0A"),
        panel.grid = element_blank(), legend.position = "right",
        plot.caption = element_text(colour = "#555555", hjust = 0))
save_fig(p5, "SCENIC_recurrent_tf_effects.png", 9.5)

# (Class-only TF consensus figure superseded by SCENIC_tf_consensus_levels.png below.)

# ---- Fig 7: pilot RSS identity sanity check ----
rss_df <- rss %>%
  mutate(regulon = str_replace(regulon, "\\.pos$", "(+)")) %>%
  group_by(subclass_name) %>% mutate(rank = rank(-rss)) %>% ungroup()
reg_order <- rss_df %>% group_by(regulon) %>% slice_max(rss, n = 1, with_ties = FALSE) %>%
  ungroup() %>% arrange(subclass_name, desc(rss)) %>% pull(regulon)
p7 <- rss_df %>%
  mutate(regulon = factor(regulon, levels = rev(unique(reg_order))),
         subclass_name = str_wrap(subclass_name, 18)) %>%
  ggplot(aes(subclass_name, regulon)) +
  geom_point(aes(size = rss, colour = rss)) +
  scale_colour_gradient(low = "#F0EFEC", high = "#E8590C", name = "RSS",
                        breaks = c(0.1, 0.2, 0.3, 0.4, 0.5), limits = c(0, 0.6)) +
  scale_size_area(max_size = 7, name = "RSS", breaks = c(0.1, 0.2, 0.3, 0.4, 0.5), limits = c(0, 0.6)) +
  guides(colour = "legend") +
  labs(x = NULL, y = NULL) +
  theme_scenic(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.text.y = element_text(face = "italic", size = 10),
        legend.position = "right")
save_fig(p7, "SCENIC_pilot_rss.png", 10)

# ---- Multi-level figures: class vs subclass vs supertype ----
# Needs test_scenic_regulons_animal_level.R --level {subclass,supertype} and,
# for the permutation panels, test_scenic_split_permutation.R.
lv_file <- function(stub, lv) file.path(tbl_dir, glue("{stub}{if (lv == 'class') '' else paste0('_', lv)}_{al_date}.csv"))
levels_have <- c("class", "subclass", "supertype")[file.exists(sapply(c("class", "subclass", "supertype"),
                                                                   function(lv) lv_file("SCENIC_animal_level_tests", lv)))]
if (length(levels_have) > 1) {
  lv_cols <- c(class = "#27455C", subclass = "#FF9240", supertype = "#E8590C")
  read_lv <- function(stub) map_dfr(levels_have, ~ read_csv(lv_file(stub, .x), show_col_types = FALSE) %>%
                                      mutate(level = .x)) %>% mutate(level = factor(level, levels_have))
  all_t <- read_lv("SCENIC_animal_level_tests")
  all_c <- read_lv("SCENIC_tf_consensus_animal_level")
  all_s <- read_lv("SCENIC_animal_level_split_summary") %>% filter(tested)
  perm_split_f <- list.files(tbl_dir, "^SCENIC_split_permutation_tests_\\d{4}-\\d{2}-\\d{2}\\.csv$", full.names = TRUE)
  perm_level_f <- list.files(tbl_dir, "^SCENIC_level_permutation_tests_\\d{4}-\\d{2}-\\d{2}\\.csv$", full.names = TRUE)
  perm_null_f <- list.files(file.path(wkdir, "data/interim/scenic/animal_level"),
                            "^level_permutation_null_\\d{4}-\\d{2}-\\d{2}\\.csv$", full.names = TRUE)
  have_perm <- length(perm_split_f) && length(perm_level_f) && length(perm_null_f)
  split_short <- function(level, class_name, subclass_name, supertype_name)
    str_remove(case_when(level == "class" ~ class_name, level == "subclass" ~ subclass_name,
                         TRUE ~ supertype_name), "^\\d+\\s+")

  # -- Fig: design and coverage at each resolution --
  cov_lv <- all_s %>%
    mutate(cover = case_when(n_animals_WildR_ok == 6 & n_animals_SPF_ok == 6 ~ "6 vs 6",
                             pmin(n_animals_WildR_ok, n_animals_SPF_ok) >= 5 ~ "≥ 5 per group",
                             TRUE ~ "3–4 per group"),
           cover = factor(cover, c("3–4 per group", "≥ 5 per group", "6 vs 6")))
  lv_lab <- cov_lv %>% count(level) %>% mutate(lab = glue("{level}\n({n} GRNs)"))
  lv_labeller <- setNames(lv_lab$lab, lv_lab$level)
  pA <- ggplot(cov_lv, aes(level, n_cells, colour = tissue)) +
    geom_boxplot(aes(fill = tissue), alpha = 0.12, outlier.shape = NA, width = 0.7,
                 position = position_dodge(width = 0.8)) +
    geom_point(size = 1, alpha = 0.55, position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8, seed = 1)) +
    scale_y_log10(labels = label_comma(), breaks = c(100, 300, 1000, 3000, 1e4, 3e4)) +
    scale_x_discrete(labels = lv_labeller) +
    scale_colour_manual(values = tissue_cols, name = NULL) +
    scale_fill_manual(values = tissue_cols, guide = "none") +
    labs(x = NULL, y = "Nuclei per GRN (log scale)", tag = "A") +
    theme_scenic()
  pB <- ggplot(cov_lv, aes(tissue, fill = cover)) +
    geom_bar(position = "fill", width = 0.75) +
    facet_wrap(~level, labeller = as_labeller(lv_labeller)) +
    scale_fill_manual(values = c("3–4 per group" = "#D4D4D4", "≥ 5 per group" = "#8FB3D3", "6 vs 6" = "#2E6FA8"),
                      name = NULL) +
    scale_y_continuous(labels = label_percent()) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(x = NULL, y = "Share of GRNs by animals with ≥ 10 nuclei", tag = "B") +
    theme_scenic() + theme(strip.text = element_text(hjust = 0.5, size = 11))
  pC <- ggplot(cov_lv, aes(level, n_regulons, colour = tissue)) +
    geom_boxplot(aes(fill = tissue), alpha = 0.12, outlier.shape = NA, width = 0.7,
                 position = position_dodge(width = 0.8)) +
    geom_point(size = 1, alpha = 0.55, position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8, seed = 1)) +
    scale_x_discrete(labels = lv_labeller) +
    scale_colour_manual(values = tissue_cols, name = NULL) +
    scale_fill_manual(values = tissue_cols, guide = "none") +
    labs(x = NULL, y = "Regulons called per GRN", tag = "C") +
    theme_scenic()
  p_cov <- (pA | pB | pC) + plot_layout(widths = c(1.15, 1, 1), guides = "collect") & theme(legend.position = "top")
  save_fig(p_cov, "SCENIC_levels_coverage.png", 5.6)

  # -- Fig: is there more signal than chance at any level? --
  qq <- all_t %>% filter(!is.na(p_animal)) %>% group_by(level) %>% arrange(p_animal) %>%
    mutate(expected = -log10(ppoints(n())), observed = -log10(p_animal)) %>% ungroup()
  pA <- ggplot(qq, aes(expected, observed, colour = level)) +
    geom_abline(linetype = "dashed", colour = "#A5A5A5") +
    geom_point(size = 0.9, alpha = 0.8) +
    scale_colour_manual(values = lv_cols, name = NULL) +
    coord_equal() +
    labs(x = "Expected −log10 p", y = "Observed −log10 p", tag = "A") +
    theme_scenic()
  cnt <- all_t %>% group_by(level) %>%
    summarise(tests = sum(!is.na(p_animal)), Observed = sum(p_animal < 0.05, na.rm = TRUE),
              `Expected (5%)` = round(0.05 * tests), .groups = "drop") %>%
    pivot_longer(c(Observed, `Expected (5%)`), names_to = "what", values_to = "n")
  pB <- ggplot(cnt, aes(level, n, fill = what)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(aes(label = comma(n)), position = position_dodge(width = 0.75), vjust = -0.4, size = 3.4,
              colour = "#555555") +
    scale_fill_manual(values = c(Observed = "#E8590C", `Expected (5%)` = "#D4D4D4"), name = NULL) +
    scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "Regulon tests with animal p < 0.05", tag = "B") +
    theme_scenic()
  if (have_perm) {
    pnull <- read_csv(max(perm_null_f), show_col_types = FALSE) %>% mutate(level = factor(level, levels_have))
    plev <- read_csv(max(perm_level_f), show_col_types = FALSE) %>%
      mutate(level = factor(level, levels_have), lab = glue("observed {n_nominal_obs}\np_perm = {round(p_perm_excess, 2)}"))
    pC <- ggplot(pnull, aes(n_nominal)) +
      geom_histogram(bins = 40, fill = "#BDBDBD", colour = NA) +
      geom_vline(data = plev, aes(xintercept = n_nominal_obs), colour = "#E8590C", linewidth = 0.9) +
      geom_text(data = plev, aes(x = Inf, y = Inf, label = lab), hjust = 1.05, vjust = 1.3, size = 3.1,
                colour = "#0A0A0A", lineheight = 0.9) +
      facet_wrap(level ~ tissue, scales = "free", ncol = 2, labeller = label_wrap_gen(multi_line = FALSE)) +
      scale_x_continuous(labels = label_comma(), n.breaks = 4) +
      labs(x = "Regulon tests at p < 0.05 under 924 animal relabellings", y = "Relabellings", tag = "C") +
      theme_scenic(11.5) + theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
    psplit <- read_csv(max(perm_split_f), show_col_types = FALSE) %>% mutate(level = factor(level, levels_have))
    pD <- ggplot(psplit, aes(p_perm, fill = level)) +
      geom_histogram(breaks = seq(0, 1, 0.1), colour = "white") +
      geom_hline(data = psplit %>% count(level) %>% mutate(e = n / 10), aes(yintercept = e),
                 linetype = "dashed", colour = "#555555") +
      facet_wrap(~level, ncol = 1, scales = "free_y") +
      scale_fill_manual(values = lv_cols, guide = "none") +
      scale_x_continuous(breaks = seq(0, 1, 0.25)) +
      labs(x = "Per-GRN permutation p (excess of p < 0.05)", y = "GRNs", tag = "D") +
      theme_scenic(11.5)
    p_null <- (pA | pB) / (pC | pD) + plot_layout(heights = c(1, 1.35))
    save_fig(p_null, "SCENIC_levels_null.png", 11.5)
  } else {
    save_fig((pA | pB), "SCENIC_levels_null.png", 5.5)
  }

  # -- Fig: superplots for the within-split FDR hits (all levels) --
  hits <- all_t %>% filter(fdr_split < 0.10) %>%
    mutate(split = split_short(level, class_name, subclass_name, supertype_name))
  top_hits <- hits %>% filter(n_animals_WildR >= 5, n_animals_SPF >= 5) %>%
    arrange(fdr_split, p_animal) %>% slice_head(n = 12) %>%
    mutate(panel = glue("{regulon} · {tissue} {level}\n{str_trunc(split, 30)}\nq = {signif(fdr_split, 2)}, g = {round(hedges_g, 1)}"),
           panel = fct_inorder(panel))
  pa_hits <- map_dfr(levels_have, function(lv) {
    f <- file.path(wkdir, glue("data/interim/scenic/animal_level/regulon_auc_per_animal{if (lv == 'class') '' else paste0('_', lv)}_{al_date}.csv"))
    st <- top_hits %>% filter(level == lv)
    if (!nrow(st) || !file.exists(f)) return(NULL)
    data.table::fread(f)[stem %in% st$stem] %>% as_tibble() %>%
      inner_join(st %>% select(stem, regulon, panel), by = c("stem", "regulon"))
  }) %>% filter(n_cells >= 10) %>% mutate(panel = factor(panel, levels(top_hits$panel)))
  p_sp <- ggplot(pa_hits, aes(microbiome, mean_auc, colour = microbiome)) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "#0A0A0A", linewidth = 0.35) +
    geom_point(aes(size = n_cells), position = position_jitter(width = 0.1, seed = 1), alpha = 0.9) +
    facet_wrap(~panel, scales = "free_y", ncol = 4) +
    scale_colour_manual(values = group_cols, guide = "none") +
    scale_y_continuous(labels = label_number(drop0trailing = TRUE),
                       n.breaks = 4) +
    scale_size_area(max_size = 4.5, name = "Nuclei per animal", labels = label_comma()) +
    labs(x = NULL, y = "Mean AUCell per animal") +
    theme_scenic(12) + theme(strip.text = element_text(size = 9.5, lineheight = 0.95))
  save_fig(p_sp, "SCENIC_fdr_hits_superplots.png", 11)

  # -- Fig: cross-split TF consensus at each level + candidate TFs --
  cand <- c("Mllt10", "Ube2k", "Zbtb20", "Klf12", "Bptf", "Rfx3", "Bclaf1", "Irf2")
  cons_lab <- all_c %>% group_by(level, tissue) %>% slice_min(p_animal, n = 4) %>% ungroup()
  pA <- ggplot(all_c, aes(hedges_g, -log10(p_animal))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#A5A5A5") +
    geom_point(aes(size = n_splits), colour = "#A5A5A5", alpha = 0.55) +
    geom_point(data = all_c %>% filter(tf %in% cand), aes(size = n_splits), colour = "#E8590C") +
    geom_text_repel(data = bind_rows(cons_lab, all_c %>% filter(tf %in% c("Mllt10", "Zbtb20", "Klf12"), tissue == "HYP")) %>%
                      distinct(level, tissue, tf, .keep_all = TRUE),
                    aes(label = tf), size = 3, fontface = "italic", seed = 3, min.segment.length = 0,
                    max.overlaps = Inf, box.padding = 0.35) +
    facet_grid(tissue ~ level) +
    scale_size_continuous(range = c(0.6, 3.5), name = "Splits pooled") +
    labs(x = "Consensus Hedges g (WildR − SPF)", y = "−log10 animal-level p", tag = "A") +
    theme_scenic(12) + theme(strip.text.y = element_text(angle = 0))
  cd <- all_c %>% filter(tf %in% cand) %>% mutate(tf = factor(tf, rev(cand)))
  pB <- ggplot(cd, aes(hedges_g, tf, colour = level)) +
    geom_vline(xintercept = 0, colour = "#A5A5A5") +
    geom_point(aes(size = n_splits), position = position_dodge(width = 0.6), alpha = 0.9) +
    scale_colour_manual(values = lv_cols, name = NULL) +
    scale_size_continuous(range = c(1.5, 4.5), name = "Splits pooled") +
    facet_wrap(~tissue) +
    labs(x = "Consensus Hedges g (WildR − SPF)", y = NULL, tag = "B") +
    theme_scenic() + theme(axis.text.y = element_text(face = "italic", colour = "#0A0A0A"))
  p_cons <- pA / pB + plot_layout(heights = c(1.5, 1))
  save_fig(p_cons, "SCENIC_tf_consensus_levels.png", 12)

  # -- Display tables --
  write_csv(
    hits %>%
      transmute(Level = level, Tissue = tissue, Split = coalesce(supertype_name, subclass_name, class_name),
                Regulon = regulon, Direction = direction,
                `Animals (W/S)` = glue("{n_animals_WildR}/{n_animals_SPF}"),
                `Mean AUC SPF` = signif(animal_mean_SPF, 2), `Mean AUC WildR` = signif(animal_mean_WildR, 2),
                `Hedges g` = round(hedges_g, 1), `Animals separated` = animals_separated,
                `p (animal)` = signif(p_animal, 2), `FDR (split)` = signif(fdr_split, 2), Nuclei = n_cells) %>%
      arrange(Level, `FDR (split)`),
    file.path(tbl_dir, glue("SCENIC_fdr_hits_all_levels_display_{al_date}.csv")))

  lvl_sum <- all_t %>% group_by(level, tissue) %>%
    summarise(tests = sum(!is.na(p_animal)), nominal = sum(p_animal < 0.05, na.rm = TRUE),
              fdr_split = sum(fdr_split < 0.10, na.rm = TRUE), fdr_global = sum(fdr_global < 0.10, na.rm = TRUE),
              # complete separation at 6 vs 6: chance = 2 / choose(12, 6) per test
              n66 = sum(n_animals_WildR == 6 & n_animals_SPF == 6 & !is.na(p_animal)),
              sep6 = sum(animals_separated == 6 & n_animals_WildR == 6 & n_animals_SPF == 6, na.rm = TRUE),
              .groups = "drop") %>%
    left_join(all_s %>% group_by(level, tissue) %>%
                summarise(grns = n(), med_nuclei = median(n_cells), med_reg = median(n_regulons),
                          full = sum(n_animals_WildR_ok == 6 & n_animals_SPF_ok == 6), .groups = "drop"),
              by = c("level", "tissue"))
  if (have_perm) lvl_sum <- lvl_sum %>%
    left_join(read_csv(max(perm_level_f), show_col_types = FALSE) %>% mutate(level = factor(level, levels_have)) %>%
                select(level, tissue, perm_median, p_perm_excess), by = c("level", "tissue"))
  write_csv(
    lvl_sum %>% arrange(level, tissue) %>%
      transmute(Level = level, Tissue = tissue, GRNs = grns, `Median nuclei` = round(med_nuclei),
                `GRNs at 6 vs 6` = full, `Median regulons` = round(med_reg), Tests = tests,
                `p < 0.05` = nominal, `Expected (5%)` = round(0.05 * tests),
                `Permutation median` = if (have_perm) round(perm_median) else NA,
                `Permutation p (excess)` = if (have_perm) round(p_perm_excess, 2) else NA,
                `Complete separation (6 vs 6)` = sep6, `Expected separations` = round(n66 * 2 / choose(12, 6)),
                `FDR < 0.10 (split)` = fdr_split, `FDR < 0.10 (global)` = fdr_global),
    file.path(tbl_dir, glue("SCENIC_level_summary_display_{al_date}.csv")))

  # Prespecified check: populations with popAlign state/abundance leads
  # (0X_popAlign_SN.qmd, MeA–BST–medial-hypothalamus network)
  if (have_perm) {
    circuit_pat <- "MEA-BST Lhx6 Nr2e1|DMH Hmx2|MEA-BST Lhx6 Sp9|TMd-DMH Foxd2|PVHd-DMH Lhx6|PVH-SO-PVa Otp|MEA-BST Otp Zic2"
    write_csv(
      read_csv(max(perm_split_f), show_col_types = FALSE) %>%
        filter(level != "class") %>%
        mutate(Population = if_else(level == "subclass", subclass_name, supertype_name)) %>%
        filter(str_detect(Population, circuit_pat)) %>%
        arrange(p_perm) %>%
        transmute(Level = level, Tissue = tissue, Population,
                  `Animals (W/S)` = glue("{n_animals_WildR}/{n_animals_SPF}"), `Regulons tested` = n_tested,
                  `p < 0.05 observed` = n_nominal_obs, `Permutation mean` = round(perm_mean, 1),
                  `Permutation p` = round(p_perm, 3)),
      file.path(tbl_dir, glue("SCENIC_popalign_circuit_permutation_display_{al_date}.csv")))
  }
}

# ---- Display tables ----
write_csv(
  split_sum %>%
    transmute(Tissue = tissue, Class = class_short, Nuclei = n_cells, Regulons = n_regulons,
              `Animals (W/S)` = glue("{n_animals_WildR_ok}/{n_animals_SPF_ok}"),
              `Cell-level hits` = n_cell_level_hits,
              `…also animal p<0.05` = n_cell_hits_nominal_animal,
              `Animal p<0.05` = n_nominal_p05, `Expected by chance` = n_expected_p05,
              `Animal FDR<0.10` = n_fdr10, Tested = tested),
  file.path(tbl_dir, glue("SCENIC_split_summary_display_{al_date}.csv")))

log_msg("Done → ", fig_dir)
