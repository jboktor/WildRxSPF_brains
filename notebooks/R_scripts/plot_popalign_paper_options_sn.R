#!/usr/bin/env Rscript
# Build paper-figure options from the completed snRNA-seq PopAlign run.
#
# This is a display/figure step only. It reads saved synthesis and per-unit
# CSV/CSV.GZ files; it does not refit PopAlign or submit jobs.
#
# Usage:
#   /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript \
#     notebooks/R_scripts/plot_popalign_paper_options_sn.R \
#     --run-tag 2026-10-04 --synth-tag 2026-10-05 --output-date 2026-11-06

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(patchwork)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1]]
}

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
run_tag <- get_arg("--run-tag", "2026-10-04")
synth_tag <- get_arg("--synth-tag", "2026-10-05")
output_date <- get_arg("--output-date", "2026-11-06")

pa_base <- file.path(wkdir, "data/interim/popalign/snRNASeq")
manifest_f <- file.path(pa_base, glue("manifest_{run_tag}.csv"))
tbl_in <- file.path(
  wkdir, "tables/popAlign/snRNASeq",
  glue("run_{run_tag}_synth-{synth_tag}")
)
fig_out <- file.path(
  wkdir, "figures/popAlign/snRNASeq",
  glue("paper_figure_options_{output_date}")
)
tbl_out <- file.path(
  wkdir, "tables/popAlign/snRNASeq",
  glue("paper_figure_options_{output_date}")
)
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_out, recursive = TRUE, showWarnings = FALSE)

required <- c(
  manifest_f,
  file.path(tbl_in, "popAlign_subpopulation_catalog.csv"),
  file.path(tbl_in, "popAlign_state_shift_tests.csv"),
  file.path(tbl_in, "popAlign_followup_leads.csv"),
  file.path(tbl_in, "popAlign_circuit_enrichment.csv")
)
stopifnot(all(file.exists(required)))

catalog <- read_csv(
  file.path(tbl_in, "popAlign_subpopulation_catalog.csv"),
  show_col_types = FALSE
)
state <- read_csv(
  file.path(tbl_in, "popAlign_state_shift_tests.csv"),
  show_col_types = FALSE
)
leads <- read_csv(
  file.path(tbl_in, "popAlign_followup_leads.csv"),
  show_col_types = FALSE
)
enrichment <- read_csv(
  file.path(tbl_in, "popAlign_circuit_enrichment.csv"),
  show_col_types = FALSE
)
manifest <- read_csv(manifest_f, show_col_types = FALSE)
gene_map <- read_csv(
  file.path(
    wkdir, "data/interim/seurat/snRNASeq",
    "seurat_gene_metadata_2026-10-05.csv"
  ),
  show_col_types = FALSE
) %>%
  distinct(export_symbol, .keep_all = TRUE)

gene_id <- setNames(gene_map$gene_id, gene_map$export_symbol)
gene_name <- setNames(gene_map$gene_name, gene_map$export_symbol)

# The same hypothesis-defined network used in build_popalign_report_assets_sn.R.
circuit_pattern <- "MEA|BST|PVH|DMH|AHN|TMd|MPN|MPO|VMH|SBPV"

COL <- c(
  SPF = "#27455C",
  WildR = "#E66101",
  circuit = "#E66101",
  other = "#A5A5A5",
  AMY = "#E66101",
  HYP = "#27455C",
  class = "#27455C",
  subclass = "#F39C56",
  supertype = "#C84C00",
  synaptic = "#E66101",
  adhesion = "#27455C",
  other_gene = "#A5A5A5"
)

theme_paper <- function(base_size = 10.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(colour = "#262626"),
      axis.text = element_text(colour = "#3A3A3A"),
      strip.background = element_rect(fill = "#F4F4F4", colour = NA),
      strip.text = element_text(face = "bold", colour = "#262626"),
      legend.position = "top",
      legend.justification = "left",
      plot.title = element_text(face = "bold", size = rel(1.2)),
      plot.subtitle = element_text(colour = "#555555", size = rel(0.92)),
      plot.caption = element_text(colour = "#666666", hjust = 0, size = rel(0.8)),
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )
}

save_option <- function(plot, stem, width, height) {
  ggsave(
    file.path(fig_out, paste0(stem, ".png")),
    plot, width = width, height = height, dpi = 320,
    bg = "white", limitsize = FALSE
  )
  ggsave(
    file.path(fig_out, paste0(stem, ".pdf")),
    plot, width = width, height = height,
    device = cairo_pdf, bg = "white", limitsize = FALSE
  )
}

short_name <- function(x) {
  x %>%
    str_remove("^\\d+\\s+") %>%
    str_remove("_\\d+$") %>%
    str_replace_all(" Gaba", " GABA")
}

short_lead <- function(tissue, name, level) {
  as.character(glue("{tissue} · {short_name(name)}\n{level}"))
}

mean_ci <- function(x) {
  n <- sum(!is.na(x))
  m <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(n)
  tibble(
    y = m,
    ymin = m - qt(0.975, df = n - 1) * se,
    ymax = m + qt(0.975, df = n - 1) * se
  )
}

unit_row <- function(unit_id) {
  out <- manifest %>% filter(.data$unit_id == !!unit_id)
  if (nrow(out) != 1) stop("Expected one manifest row for ", unit_id)
  out
}

read_unit <- function(unit_id, filename) {
  row <- unit_row(unit_id)
  read_csv(file.path(row$output_dir[[1]], filename), show_col_types = FALSE)
}

# ── Option 1: scope + all tests + circuit enrichment ─────────────────────────

scope_counts <- catalog %>%
  group_by(level, tissue) %>%
  summarise(
    units = n_distinct(unit_id),
    components = n(),
    robust_clean_components = sum(robust & !contamination_flag, na.rm = TRUE),
    .groups = "drop"
  )

p_scope <- bind_rows(
  scope_counts %>% transmute(level, tissue, measure = "Units modelled", n = units),
  scope_counts %>% transmute(level, tissue, measure = "Components tested", n = components)
) %>%
  mutate(
    level = factor(level, levels = c("class", "subclass", "supertype")),
    measure = factor(measure, levels = c("Units modelled", "Components tested"))
  ) %>%
  ggplot(aes(n, level, fill = tissue)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.72),
    hjust = -0.15, size = 3
  ) +
  facet_wrap(~ measure, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = COL[c("AMY", "HYP")]) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    x = NULL, y = NULL, fill = NULL,
    title = "Scope: two tissues × three Allen resolutions",
    subtitle = glue(
      "{n_distinct(catalog$unit_id)} units; {nrow(catalog)} components; ",
      "6 WildR vs 6 SPF animals per tissue"
    )
  ) +
  theme_paper() +
  theme(panel.grid.major.y = element_blank())

all_tests <- catalog %>%
  mutate(
    circuit = str_detect(name, circuit_pattern),
    clean_robust = robust & !contamination_flag,
    level = factor(level, levels = c("class", "subclass", "supertype"))
  ) %>%
  dplyr::select(
    unit_id, level, tissue, name, comp, circuit, clean_robust,
    p_state, q_state, p_abund, q_abund
  ) %>%
  pivot_longer(
    c(p_state, p_abund),
    names_to = "test",
    values_to = "p"
  ) %>%
  mutate(
    q = if_else(test == "p_state", q_state, q_abund),
    test = recode(
      test,
      p_state = "State shift",
      p_abund = "Abundance shift"
    )
  ) %>%
  filter(!is.na(p)) %>%
  group_by(test, level) %>%
  arrange(tissue, name, comp, .by_group = TRUE) %>%
  mutate(test_index = row_number()) %>%
  ungroup()

label_tests <- all_tests %>%
  filter(
    clean_robust,
    (test == "State shift" & circuit & p <= 0.0065) |
      (test == "Abundance shift" & str_detect(name, "MEA-BST Otp Zic2") & p < 0.05)
  ) %>%
  mutate(label = str_wrap(short_name(name), 18))

p_tests <- ggplot(
  all_tests,
  aes(test_index, -log10(p))
) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 2, linewidth = 0.4, colour = "#777777"
  ) +
  geom_point(
    aes(
      colour = if_else(circuit, "Circuit", "Other units"),
      shape = tissue,
      alpha = clean_robust
    ),
    size = 1.35, stroke = 0
  ) +
  geom_text_repel(
    data = label_tests,
    aes(label = label),
    size = 2.4,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.25,
    point.padding = 0.2,
    seed = 7,
    show.legend = FALSE
  ) +
  facet_grid(
    rows = vars(test), cols = vars(level),
    scales = "free_x", space = "free_x"
  ) +
  scale_colour_manual(
    values = c(Circuit = COL[["circuit"]], `Other units` = COL[["other"]])
  ) +
  scale_shape_manual(values = c(AMY = 16, HYP = 17)) +
  scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.16), guide = "none") +
  scale_x_continuous(breaks = NULL) +
  labs(
    x = "Every tested component, ordered by tissue and cell type",
    y = expression(-log[10]("exact-permutation p")),
    colour = NULL, shape = NULL,
    title = "All component-level tests",
    subtitle = "Orange: hypothesis-defined MeA/BST/medial-hypothalamus circuit; dashed line: nominal p = 0.05",
    caption = "Faded points are components that fail robustness or contamination criteria. No individual state or abundance test survives BH FDR."
  ) +
  theme_paper(9.5) +
  theme(
    panel.spacing.x = unit(0.7, "lines"),
    panel.grid.major.x = element_blank()
  )

rate_ci <- function(hits, total) {
  ci <- binom.test(hits, total)$conf.int
  tibble(rate = hits / total, low = ci[[1]], high = ci[[2]])
}

enrich_plot <- enrichment %>%
  filter(
    outcome == "hit",
    grouping == "circuit",
    scope %in% c("all levels", "subclass only", "supertype only")
  ) %>%
  mutate(
    scope = recode(
      scope,
      `all levels` = "All levels",
      `subclass only` = "Subclass",
      `supertype only` = "Supertype"
    )
  )

enrich_rates <- map_dfr(seq_len(nrow(enrich_plot)), function(i) {
  r <- enrich_plot[i, ]
  bind_rows(
    rate_ci(r$in_group_hits, r$in_group_units) %>%
      mutate(group = "Circuit", scope = r$scope, odds_ratio = r$odds_ratio, p = r$p),
    rate_ci(r$out_group_hits, r$out_group_units) %>%
      mutate(group = "Other units", scope = r$scope, odds_ratio = r$odds_ratio, p = r$p)
  )
}) %>%
  mutate(
    scope = factor(scope, levels = c("All levels", "Subclass", "Supertype")),
    group = factor(group, levels = c("Circuit", "Other units"))
  )

enrich_labels <- enrich_plot %>%
  transmute(
    scope = factor(scope, levels = c("All levels", "Subclass", "Supertype")),
    label = glue("OR {sprintf('%.1f', odds_ratio)}; p = {format(signif(p, 2), scientific = FALSE)}")
  )

p_enrichment <- ggplot(enrich_rates, aes(rate, scope, colour = group)) +
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0.16, linewidth = 0.55) +
  geom_point(size = 2.6) +
  geom_text(
    data = enrich_labels,
    aes(x = 0.74, y = scope, label = label),
    inherit.aes = FALSE,
    hjust = 1, size = 3
  ) +
  scale_colour_manual(
    values = c(Circuit = COL[["circuit"]], `Other units` = COL[["other"]])
  ) +
  scale_x_continuous(
    labels = scales::percent,
    limits = c(0, 0.75),
    breaks = c(0, 0.2, 0.4, 0.6)
  ) +
  labs(
    x = "Units with ≥1 nominal state shift",
    y = NULL, colour = NULL,
    title = "State shifts concentrate in the circuit",
    subtitle = "Exact binomial 95% CIs; Fisher odds ratio and p shown"
  ) +
  theme_paper() +
  theme(panel.grid.major.y = element_blank())

option1 <- p_scope + p_tests + p_enrichment +
  plot_layout(
    design = "
    AB
    CB
    ",
    widths = c(1.05, 2.45),
    heights = c(1, 1)
  ) +
  plot_annotation(
    title = "PopAlign tests broad cell-state space and localizes the strongest pattern to one circuit",
    subtitle = "All statistics use animals as replicates (6 WildR vs 6 SPF); individual component signals are exploratory",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11, colour = "#555555")
    )
  )

save_option(
  option1,
  glue("popAlign_option1_scope_all_tests_{output_date}"),
  width = 15.8, height = 9.4
)

write_csv(
  all_tests,
  file.path(tbl_out, glue("popAlign_scope_all_tests_{output_date}.csv"))
)

# ── Option 2: animal-level MeA–BST leads + signatures ─────────────────────────

state_focus_ids <- c(
  "supertype__AMY__CS20230722_SUPT_0355",
  "subclass__AMY__CS20230722_SUBC_075",
  "subclass__HYP__CS20230722_SUBC_074"
)
abundance_focus_id <- "subclass__AMY__CS20230722_SUBC_121"

state_focus <- leads %>%
  filter(unit_id %in% state_focus_ids, effect == "state shift") %>%
  left_join(
    state %>% dplyr::select(unit_id, comp, top_shift_feature, top_shift_value),
    by = c("unit_id", "comp")
  ) %>%
  mutate(
    display = case_when(
      unit_id == "supertype__AMY__CS20230722_SUPT_0355" ~
        "AMY · Lhx6 Nr2e1 GABA\nsupertype 0355",
      unit_id == "subclass__AMY__CS20230722_SUBC_075" ~
        "AMY · Lhx6 Nr2e1 GABA\nsubclass 075",
      unit_id == "subclass__HYP__CS20230722_SUBC_074" ~
        "HYP · Lhx6 Sp9 GABA\nsubclass 074",
      TRUE ~ short_lead(tissue, name, level)
    )
  )

state_animals <- map_dfr(seq_len(nrow(state_focus)), function(i) {
  r <- state_focus[i, ]
  centroids <- read_unit(r$unit_id, "animal_component_centroids.csv")
  value <- centroids %>%
    filter(comp == r$comp) %>%
    pull(all_of(r$top_shift_feature))
  centroids %>%
    filter(comp == r$comp) %>%
    transmute(
      unit_id = r$unit_id,
      tissue = r$tissue,
      level = r$level,
      cell_type = r$name,
      display = r$display,
      sample_id,
      microbiome,
      comp,
      feature = r$top_shift_feature,
      raw_value = value,
      z_value = as.numeric(scale(value)),
      p_perm_multivariate = r$p,
      q_perm_multivariate = r$q
    )
})

state_annotations <- state_focus %>%
  transmute(
    display = factor(display, levels = state_focus$display),
    label = glue("multivariate p = {signif(p, 2)}\nq = {signif(q, 2)}")
  )

p_state_animals <- state_animals %>%
  mutate(
    microbiome = factor(microbiome, levels = c("SPF", "WildR")),
    display = factor(display, levels = state_focus$display)
  ) %>%
  ggplot(aes(microbiome, z_value, colour = microbiome)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "#D0D0D0") +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 2.2, alpha = 0.82
  ) +
  stat_summary(fun.data = mean_ci, geom = "errorbar", width = 0.12, linewidth = 0.65) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3.4) +
  geom_text(
    data = state_annotations,
    aes(x = 1.5, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.15, size = 2.8
  ) +
  facet_wrap(~ display, nrow = 1) +
  scale_colour_manual(values = COL[c("SPF", "WildR")]) +
  labs(
    x = NULL,
    y = "Top shifted oNMF feature (z-score)",
    colour = NULL,
    title = "State shifts recur across related MeA–BST inhibitory populations",
    subtitle = "Each point is one animal; diamonds and bars are means and t-based 95% CIs"
  ) +
  theme_paper() +
  theme(legend.position = "none")

abundance_focus <- leads %>%
  filter(unit_id == abundance_focus_id, effect == "abundance shift") %>%
  slice_head(n = 1)
stopifnot(nrow(abundance_focus) == 1)

abundance_animals <- read_unit(
  abundance_focus_id,
  "animal_component_fractions.csv"
) %>%
  filter(comp == abundance_focus$comp) %>%
  mutate(
    unit_id = abundance_focus_id,
    cell_type = abundance_focus$name,
    p_perm = abundance_focus$p,
    q_perm = abundance_focus$q,
    p_depth_adjusted = abundance_focus$p_abund_glm_depthadj
  )

p_abundance_animals <- abundance_animals %>%
  mutate(microbiome = factor(microbiome, levels = c("SPF", "WildR"))) %>%
  ggplot(aes(microbiome, soft_frac, colour = microbiome)) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 2.4, alpha = 0.85
  ) +
  stat_summary(fun.data = mean_ci, geom = "errorbar", width = 0.12, linewidth = 0.65) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3.6) +
  scale_colour_manual(values = COL[c("SPF", "WildR")]) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  labs(
    x = NULL,
    y = "Soft component fraction",
    colour = NULL,
    title = "A glutamatergic MeA–BST population expands in WildR",
    subtitle = glue(
      "AMY · MEA–BST Otp Zic2 Glut · component {abundance_focus$comp}\n",
      "exact p = {signif(abundance_focus$p, 2)}; q = {signif(abundance_focus$q, 2)}; ",
      "depth-adjusted GLM p = {signif(abundance_focus$p_abund_glm_depthadj, 2)}"
    )
  ) +
  theme_paper() +
  theme(legend.position = "none")

program_genes_focus <- map_dfr(seq_len(nrow(state_focus)), function(i) {
  r <- state_focus[i, ]
  read_unit(r$unit_id, "onmf_top_genes.csv") %>%
    filter(feature == r$top_shift_feature) %>%
    arrange(rank) %>%
    slice_head(n = 10) %>%
    mutate(
      unit_id = r$unit_id,
      display = r$display,
      gene_name = coalesce(unname(gene_name[gene]), gene)
    )
})

synaptic_genes <- c(
  "Gria1", "Gria2", "Grm1", "Grm7", "Grm8", "Grik1", "Grik2", "Grid1",
  "Grid2", "Erbb4", "Nrg1", "Htr2c", "Gabrg3", "Slc8a1", "Cadps2",
  "Cacna1e", "Kcnc2"
)
adhesion_genes <- c(
  "Ptprd", "Ptprt", "Lsamp", "Tenm1", "Tenm2", "Ctnna2", "Cadm1", "Cdh18",
  "Pcdh15", "Cntnap5a", "Cntnap5b", "Nrxn1", "Nlgn1", "Il1rapl1",
  "Il1rapl2", "Prickle2", "Brinp3"
)

program_genes_focus <- program_genes_focus %>%
  mutate(
    theme = case_when(
      gene_name %in% synaptic_genes ~ "Synaptic receptor / signalling",
      gene_name %in% adhesion_genes ~ "Adhesion / axon organization",
      TRUE ~ "Other top-ranked gene"
    ),
    theme = factor(
      theme,
      levels = c(
        "Synaptic receptor / signalling",
        "Adhesion / axon organization",
        "Other top-ranked gene"
      )
    ),
    display = factor(display, levels = state_focus$display)
  )

p_program_genes <- ggplot(
  program_genes_focus,
  aes(rank, fct_reorder(gene_name, -rank), colour = theme)
) +
  geom_segment(aes(x = 10.5, xend = rank, yend = gene_name), colour = "#DDDDDD") +
  geom_point(size = 2.5) +
  facet_wrap(~ display, nrow = 1, scales = "free_y") +
  scale_x_reverse(breaks = c(1, 5, 10)) +
  scale_colour_manual(
    values = c(
      `Synaptic receptor / signalling` = COL[["synaptic"]],
      `Adhesion / axon organization` = COL[["adhesion"]],
      `Other top-ranked gene` = COL[["other_gene"]]
    )
  ) +
  labs(
    x = "Rank in shifted program (1 = highest weight)",
    y = NULL,
    colour = NULL,
    title = "The shifted programs prioritize synaptic and adhesion genes",
    subtitle = "Top 10 genes from the oNMF feature with the largest group shift; colours are descriptive gene classes"
  ) +
  theme_paper(9.5) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_line(colour = "#EEEEEE", linewidth = 0.3)
  )

option2_bottom <- p_abundance_animals + p_program_genes +
  plot_layout(widths = c(0.82, 2.18))
option2 <- p_state_animals / option2_bottom +
  plot_layout(heights = c(0.95, 1.25)) +
  plot_annotation(
    title = "MeA–BST populations offer an animal-level cellular and molecular story",
    subtitle = "Convergent inhibitory state shifts plus an exploratory glutamatergic abundance shift; all panels show saved raw outputs",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11, colour = "#555555")
    )
  )

save_option(
  option2,
  glue("popAlign_option2_mea_bst_story_{output_date}"),
  width = 15.8, height = 9.5
)

write_csv(
  bind_rows(
    state_animals %>% mutate(readout = "state"),
    abundance_animals %>%
      transmute(
        unit_id, tissue = "AMY", level = "subclass",
        cell_type, display = "AMY · Otp Zic2 Glut\nsubclass 121",
        sample_id, microbiome, comp,
        feature = "soft_frac", raw_value = soft_frac,
        z_value = as.numeric(scale(soft_frac)),
        p_perm_multivariate = p_perm,
        q_perm_multivariate = q_perm,
        readout = "abundance"
      )
  ),
  file.path(tbl_out, glue("popAlign_mea_bst_animal_values_{output_date}.csv"))
)
write_csv(
  program_genes_focus,
  file.path(tbl_out, glue("popAlign_mea_bst_shifted_program_genes_{output_date}.csv"))
)

# ── Option 3: molecular themes across strong circuit state leads ─────────────

circuit_state <- leads %>%
  filter(
    effect == "state shift",
    p <= 0.0065,
    str_detect(name, circuit_pattern)
  ) %>%
  left_join(
    state %>% dplyr::select(
      unit_id, comp, top_shift_feature, top_shift_value
    ),
    by = c("unit_id", "comp")
  ) %>%
  mutate(
    display = short_lead(tissue, name, level),
    display = str_wrap(display, width = 25),
    direction = if_else(
      top_shift_value > 0,
      "WildR enriched (program higher)",
      "WildR depleted (program lower)"
    )
  )

all_program_genes <- map_dfr(seq_len(nrow(circuit_state)), function(i) {
  r <- circuit_state[i, ]
  read_unit(r$unit_id, "onmf_top_genes.csv") %>%
    filter(feature == r$top_shift_feature) %>%
    arrange(rank) %>%
    slice_head(n = 50) %>%
    mutate(
      unit_id = r$unit_id,
      comp = r$comp,
      display = r$display,
      tissue = r$tissue,
      level = r$level,
      cell_type = r$name,
      feature = r$top_shift_feature,
      top_shift_value = r$top_shift_value,
      direction = r$direction,
      gene_name = coalesce(unname(gene_name[gene]), gene),
      gene_id = unname(gene_id[gene])
    )
})

go_cache_f <- file.path(
  tbl_out,
  glue("popAlign_circuit_shifted_program_GO_{output_date}.csv")
)
go_all <- if (file.exists(go_cache_f) && !("--recompute-go" %in% args)) {
  message("Reusing shifted-program GO cache: ", go_cache_f)
  read_csv(go_cache_f, show_col_types = FALSE)
} else {
  go_result <- map_dfr(seq_len(nrow(circuit_state)), function(i) {
    r <- circuit_state[i, ]
    row <- unit_row(r$unit_id)
    genes <- all_program_genes %>%
      filter(unit_id == r$unit_id, comp == r$comp) %>%
      pull(gene_id) %>%
      na.omit() %>%
      unique()
    universe <- read_csv(
      file.path(row$output_dir[[1]], "onmf_W.csv.gz"),
      show_col_types = FALSE,
      col_select = "gene"
    ) %>%
      pull(gene) %>%
      { unname(gene_id[.]) } %>%
      na.omit() %>%
      unique()

    enriched <- enrichGO(
      gene = genes,
      universe = universe,
      OrgDb = org.Mm.eg.db,
      keyType = "ENSEMBL",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      minGSSize = 5,
      maxGSSize = 500,
      readable = FALSE
    )
    if (is.null(enriched) || !nrow(as.data.frame(enriched))) return(tibble())

    as_tibble(as.data.frame(enriched)) %>%
      mutate(
        unit_id = r$unit_id,
        comp = r$comp,
        display = r$display,
        tissue = r$tissue,
        level = r$level,
        cell_type = r$name,
        feature = r$top_shift_feature
      )
  })
  write_csv(go_result, go_cache_f)
  go_result
}

write_csv(
  all_program_genes,
  file.path(tbl_out, glue("popAlign_circuit_shifted_program_genes_{output_date}.csv"))
)

# Representative mechanistic terms are shown instead of a long list of highly
# redundant GO parent/child terms. Selection is biological and descriptive;
# significance still requires within-program BH < 0.10.
representative_terms <- c(
  "synapse assembly",
  "regulation of postsynaptic membrane potential",
  "neurotransmitter secretion",
  "ERBB4 signaling pathway",
  "cell-cell adhesion via plasma-membrane adhesion molecules",
  "dendrite development",
  "learning or memory",
  "glutamate receptor signaling pathway"
)
go_plot <- go_all %>%
  filter(Description %in% representative_terms, p.adjust < 0.10)

if (nrow(go_plot)) {
  go_plot <- go_plot %>%
    separate_wider_delim(GeneRatio, delim = "/", names = c("gene_n", "gene_d")) %>%
    mutate(
      across(c(gene_n, gene_d), as.numeric),
      gene_ratio = gene_n / gene_d,
      term = str_wrap(Description, 38)
    )
} else {
  go_plot <- tibble()
}

display_order <- circuit_state %>%
  arrange(tissue, p, level, name) %>%
  pull(display)

p_go <- if (nrow(go_plot)) {
  go_plot %>%
    mutate(
      display = factor(display, levels = display_order),
      term = fct_reorder(term, -log10(p.adjust), .fun = max)
    ) %>%
    ggplot(aes(display, term)) +
    geom_point(aes(size = gene_ratio, colour = -log10(p.adjust)), alpha = 0.9) +
    scale_colour_gradient(low = "#F3C6AA", high = "#B33C00") +
    scale_size_continuous(range = c(2, 6), labels = scales::percent) +
    labs(
      x = NULL, y = NULL,
      size = "Gene ratio",
      colour = expression(-log[10]("BH FDR")),
      title = "Shifted-program GO terms",
      subtitle = "Representative mechanistic terms from top-50 program genes; within-program BH < 0.10"
    ) +
    theme_paper(9.2) +
    theme(
      axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
      panel.grid.major = element_line(colour = "#EEEEEE", linewidth = 0.3)
    )
} else {
  ggplot() +
    annotate("text", x = 0, y = 0, label = "No shifted-program GO terms at BH < 0.10") +
    theme_void() +
    labs(title = "Shifted-program GO terms")
}

gene_recurrence <- all_program_genes %>%
  filter(rank <= 20) %>%
  count(gene_name, name = "n_programs") %>%
  filter(
    n_programs >= 2 |
      gene_name %in% c(synaptic_genes, adhesion_genes)
  ) %>%
  arrange(desc(n_programs), gene_name) %>%
  slice_head(n = 24)

gene_matrix <- all_program_genes %>%
  semi_join(gene_recurrence, by = "gene_name") %>%
  filter(rank <= 20) %>%
  left_join(gene_recurrence, by = "gene_name") %>%
  mutate(
    theme = case_when(
      gene_name %in% synaptic_genes ~ "Synaptic receptor / signalling",
      gene_name %in% adhesion_genes ~ "Adhesion / axon organization",
      TRUE ~ "Recurring top-program gene"
    ),
    display = factor(display, levels = display_order),
    gene_name = fct_reorder(gene_name, n_programs, .desc = FALSE)
  )

p_genes <- ggplot(gene_matrix, aes(display, gene_name)) +
  geom_point(aes(size = 21 - rank, colour = theme), alpha = 0.9) +
  scale_size_continuous(
    range = c(1.8, 5.5),
    breaks = c(1, 6, 11, 16, 20),
    labels = c("20", "15", "10", "5", "1")
  ) +
  scale_colour_manual(
    values = c(
      `Synaptic receptor / signalling` = COL[["synaptic"]],
      `Adhesion / axon organization` = COL[["adhesion"]],
      `Recurring top-program gene` = COL[["other_gene"]]
    )
  ) +
  labs(
    x = NULL, y = NULL,
    size = "Gene rank",
    colour = NULL,
    title = "Recurring genes and pre-specified synaptic/adhesion genes",
    subtitle = "Top 20 genes per shifted program; larger points are more highly weighted"
  ) +
  theme_paper(9.2) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
    panel.grid.major = element_line(colour = "#EEEEEE", linewidth = 0.3)
  )

p_leads <- circuit_state %>%
  mutate(
    display = factor(display, levels = rev(display_order)),
    share = share_pct / 100
  ) %>%
  ggplot(aes(-log10(p), display)) +
  geom_vline(
    xintercept = -log10(0.0065),
    linetype = 2, linewidth = 0.4, colour = "#777777"
  ) +
  geom_segment(aes(x = 0, xend = -log10(p), yend = display), colour = "#D5D5D5") +
  geom_point(aes(size = share, colour = tissue, shape = level), alpha = 0.9) +
  scale_colour_manual(values = COL[c("AMY", "HYP")]) +
  scale_shape_manual(values = c(subclass = 16, supertype = 17)) +
  scale_size_continuous(labels = scales::percent, range = c(2.5, 6)) +
  labs(
    x = expression(-log[10]("exact-permutation p")),
    y = NULL,
    colour = NULL,
    shape = "Allen level",
    size = "Component share",
    title = "Strong circuit state-shift leads",
    subtitle = "All shown components are robust and contamination-clean"
  ) +
  theme_paper(9.2) +
  theme(panel.grid.major.y = element_blank())

option3_right <- p_go / p_genes + plot_layout(heights = c(1, 1))
option3 <- p_leads + option3_right +
  plot_layout(widths = c(0.95, 2.05)) +
  plot_annotation(
    title = "Synaptic receptor and adhesion programs provide a molecular frame for the circuit leads",
    subtitle = "Exploratory program annotation across the strongest MeA/BST/medial-hypothalamus state shifts",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11, colour = "#555555")
    )
  )

save_option(
  option3,
  glue("popAlign_option3_circuit_molecular_themes_{output_date}"),
  width = 16.5, height = 11.2
)

# ── Options 4–5: direction-separated GO program annotations ─────────────────

# GO is over-representation among the top-loading genes of a shifted oNMF
# feature. Direction comes from that feature's WildR-minus-SPF centroid shift;
# it is not a gene-level differential-expression or directional GSEA statistic.
go_directed <- go_all %>%
  dplyr::select(-any_of(c("top_shift_value", "direction"))) %>%
  left_join(
    circuit_state %>%
      dplyr::select(
        unit_id, comp, top_shift_value, direction,
        state_p = p, state_q = q
      ),
    by = c("unit_id", "comp")
  ) %>%
  mutate(
    direction = factor(
      direction,
      levels = c(
        "WildR enriched (program higher)",
        "WildR depleted (program lower)"
      )
    )
  )

write_csv(
  go_directed,
  file.path(
    tbl_out,
    glue("popAlign_circuit_shifted_program_GO_directed_{output_date}.csv")
  )
)

go_direction_plot <- go_directed %>%
  filter(Description %in% representative_terms, p.adjust < 0.10) %>%
  separate_wider_delim(
    GeneRatio, delim = "/", names = c("gene_n", "gene_d")
  ) %>%
  mutate(
    across(c(gene_n, gene_d), as.numeric),
    gene_ratio = gene_n / gene_d,
    term = factor(
      Description,
      levels = rev(representative_terms)
    ),
    display = factor(display, levels = display_order),
    evidence = -log10(p.adjust)
  )

p_go_direction <- ggplot(
  go_direction_plot,
  aes(display, term)
) +
  geom_point(
    aes(size = gene_ratio, fill = evidence),
    shape = 21, colour = "white", stroke = 0.35
  ) +
  facet_grid(
    cols = vars(direction),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_gradient(low = "#F2D8C8", high = "#A63603") +
  scale_size_continuous(range = c(2.5, 7), labels = scales::percent) +
  labs(
    x = NULL,
    y = NULL,
    size = "Gene ratio",
    fill = expression(-log[10]("GO BH FDR")),
    title = "GO annotations separate into WildR-high and WildR-low cell-state programs",
    subtitle = paste(
      "Direction is the sign of the program's WildR − SPF centroid shift;",
      "GO tests use its top-50 loading genes"
    ),
    caption = paste(
      "Program direction is not directional GSEA and does not imply that every",
      "gene in the GO term changes in the same direction."
    )
  ) +
  theme_paper(10) +
  theme(
    axis.text.x = element_text(angle = 52, hjust = 1, vjust = 1),
    panel.grid.major = element_line(colour = "#EEEEEE", linewidth = 0.3),
    strip.text = element_text(size = 10.5)
  )

save_option(
  p_go_direction,
  glue("popAlign_option4_directed_GO_by_lead_{output_date}"),
  width = 16.5, height = 7.8
)

direction_summary <- go_directed %>%
  filter(Description %in% representative_terms, p.adjust < 0.10) %>%
  distinct(unit_id, comp, Description, direction) %>%
  count(Description, direction, name = "n_programs") %>%
  complete(
    Description = representative_terms,
    direction,
    fill = list(n_programs = 0L)
  ) %>%
  mutate(
    signed_n = if_else(
      direction == "WildR depleted (program lower)",
      -n_programs,
      n_programs
    ),
    total = ave(n_programs, Description, FUN = sum),
    Description = fct_reorder(Description, total)
  )

p_go_direction_summary <- ggplot(
  direction_summary,
  aes(signed_n, Description, fill = direction)
) +
  geom_vline(xintercept = 0, linewidth = 0.45, colour = "#777777") +
  geom_col(width = 0.68) +
  geom_text(
    aes(
      label = if_else(n_programs > 0, as.character(n_programs), ""),
      hjust = if_else(signed_n < 0, 1.15, -0.15)
    ),
    size = 3.4, colour = "#262626"
  ) +
  scale_fill_manual(
    values = c(
      `WildR enriched (program higher)` = COL[["WildR"]],
      `WildR depleted (program lower)` = COL[["SPF"]]
    )
  ) +
  scale_x_continuous(
    breaks = seq(-8, 8, by = 2),
    labels = abs,
    expand = expansion(mult = c(0.10, 0.10))
  ) +
  labs(
    x = "Number of independently fitted shifted programs",
    y = NULL,
    fill = NULL,
    title = "The same broad pathways occur in both directions across cell types",
    subtitle = paste(
      "Orange: program higher in WildR; blue: program lower in WildR.",
      "Counts include GO BH FDR < 0.10."
    ),
    caption = paste(
      "Counts summarize programs, not animals, and related subclass/supertype",
      "models are not statistically independent."
    )
  ) +
  theme_paper(10.5) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank()
  )

save_option(
  p_go_direction_summary,
  glue("popAlign_option5_directed_GO_summary_{output_date}"),
  width = 10.8, height = 6.3
)

message("Wrote paper figure options to: ", fig_out)
message("Wrote plotted data tables to: ", tbl_out)
