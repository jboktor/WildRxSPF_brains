#!/usr/bin/env Rscript
# Publication figure options for the completed single-nucleus D-SPIN analysis.
#
# Reads saved model summaries and raw per-animal responses; never refits D-SPIN.
# The output directory must be absent or empty so existing figures/tables cannot
# be overwritten accidentally.
#
# Usage:
#   /resnick/groups/MazmanianLab/jboktor/software/miniforge3/envs/spatialomics/bin/Rscript \
#     notebooks/R_scripts/plot_dspin_paper_options_sn.R \
#     --tag 2026-10-06 --output-tag 2026-11-06

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(patchwork)
  library(ggrepel)
  library(scico)
  library(ggbeeswarm)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1]]
}

tag <- get_arg("--tag", "2026-10-06")
output_tag <- get_arg("--output-tag", format(Sys.Date(), "%Y-%m-%d"))

wkdir <- "/resnick/groups/mthomson/jboktor/WILDRxSPF_brains"
dspin_dir <- file.path(wkdir, "data/interim/dspin")
sum_dir <- file.path(dspin_dir, glue("summary_{tag}"))
source_tbl_dir <- file.path(wkdir, "tables/DSPIN", glue("summary_{tag}"))
fig_dir <- file.path(wkdir, "figures/DSPIN/snRNASeq", glue("paper_options_{output_tag}"))
tbl_dir <- file.path(wkdir, "tables/DSPIN", glue("paper_options_{output_tag}"))
popalign_f <- file.path(
  wkdir, "tables/popAlign/snRNASeq",
  "run_2026-10-04_synth-2026-10-05",
  "popAlign_state_shift_tests.csv"
)

inputs <- c(
  program_tests = file.path(source_tbl_dir, "DSPIN_program_tests.csv"),
  unit_tests = file.path(source_tbl_dir, "DSPIN_unit_tests.csv"),
  calibration = file.path(source_tbl_dir, "DSPIN_calibration.csv"),
  global_summary = file.path(source_tbl_dir, "DSPIN_global_summary.csv"),
  animal_responses = file.path(sum_dir, glue("dspin_animal_relative_h_{tag}.csv")),
  popalign = popalign_f
)
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs)) {
  stop("Missing required inputs:\n", paste(missing_inputs, collapse = "\n"))
}

nonempty <- function(path) dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE))
if (nonempty(fig_dir) || nonempty(tbl_dir)) {
  stop(
    "Refusing to overwrite a non-empty output directory. Use a new --output-tag.\n",
    fig_dir, "\n", tbl_dir
  )
}

prog <- read_csv(inputs[["program_tests"]], show_col_types = FALSE)
unit_t <- read_csv(inputs[["unit_tests"]], show_col_types = FALSE)
cal <- read_csv(inputs[["calibration"]], show_col_types = FALSE)
global <- read_csv(inputs[["global_summary"]], show_col_types = FALSE)
animals <- read_csv(inputs[["animal_responses"]], show_col_types = FALSE)
popalign <- read_csv(inputs[["popalign"]], show_col_types = FALSE)

stopifnot(
  nrow(prog) == 2800L,
  nrow(unit_t) == 140L,
  all(c("unit_id", "program_id", "relative_h", "sample_id", "microbiome") %in% names(animals))
)

# ── House style ─────────────────────────────────────────────────────────────
mb_cols <- c(SPF = "#27455C", WildR = "#FF6C0C")
tis_cols <- c(AMY = "#3B7EA1", HYP = "#A36AAB")
dir_cols <- c("Down in WildR" = "#27455C", "Up in WildR" = "#FF6C0C")
role_cols <- c(
  "Biological anchor" = "#FF6C0C",
  "Nominal benchmark" = "#27455C",
  "Excluded" = "#A5A5A5"
)

theme_dspin_paper <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#EDEDED", linewidth = 0.3),
      panel.border = element_rect(colour = "#D0D0D0"),
      strip.background = element_rect(fill = "#FFF3EB", colour = NA),
      strip.text = element_text(face = "bold", size = rel(0.85), colour = "#222222"),
      plot.title = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle = element_text(colour = "#555555", size = rel(0.82)),
      plot.title.position = "plot",
      legend.position = "top",
      legend.justification = "left",
      axis.title = element_text(size = rel(0.9)),
      plot.tag = element_text(face = "bold", size = rel(1.15))
    )
}
theme_set(theme_dspin_paper())

short_unit <- function(x, n = 42) str_trunc(str_remove(x, "^\\d+ "), n)
fmt_p <- function(x) ifelse(
  x < 0.001,
  formatC(x, format = "e", digits = 1),
  formatC(x, format = "f", digits = 3)
)
circuit_pattern <- "MEA|BST|COA|CEA|PAA"

is_circuit_name <- function(x) str_detect(x, circuit_pattern)

# ── Candidate audit table ───────────────────────────────────────────────────
candidate_spec <- tribble(
  ~unit_id, ~program_id, ~role, ~story_short, ~story_label, ~rationale, ~caution,
  "AMY__CS20230722_SUPT_0056", "P7", "Biological anchor",
  "MEA Glut P7",
  "MEA Slc17a7 Glut_2 · P7\nAr / Esr1 / Htr7",
  "Steroid- and serotonin-receptor genes in a MeA glutamatergic unit provide a social/sex-behavior-relevant hypothesis.",
  "Social behavior was not measured; D-SPIN q is null and the popAlign unit is not nominal.",
  "AMY__CS20230722_SUPT_0338", "P10", "Biological anchor",
  "MEA–BST Sox6 P10",
  "MEA–BST Sox6 Gaba_2 · P10\nCntnap5a / Pam / Grin3a",
  "Synaptic/excitability program in a MeA–BST GABA unit, with Cntnap5a, Grin3a and Thsd7b in the MAST bridge.",
  "Nominal only; popAlign support uses the same animals and is not independent replication.",
  "AMY__CS20230722_SUPT_0355", "P5", "Biological anchor",
  "MEA–BST Lhx6 P5",
  "MEA–BST Lhx6 Nr2e1 Gaba_2 · P5\nPcdh7 / Fgf13 / Cdh13",
  "Cell-adhesion/excitability program with Thrb, Prkcb, Cdh13 and Pcdh7 in the MAST bridge; its unit also has a nominal popAlign state shift.",
  "Only one WildR animal lies beyond the SPF range for P5; the popAlign component is not the D-SPIN program.",
  "AMY__CS20230722_SUPT_0511", "P19", "Nominal benchmark",
  "COAa–PAA–MEA P19",
  "COAa–PAA–MEA Barhl2 Glut_3 · P19\nstrongest nominal D-SPIN lead",
  "Smallest attainable exact p and complete animal separation make this a useful benchmark for the scale of the nominal signals.",
  "Top genes and DAVID annotation do not provide a compelling social-behavior mechanism.",
  "AMY__CS20230722_SUPT_0384", "P7", "Excluded",
  "CEA–BST Pdyn P7",
  "CEA–BST Ebf1 Pdyn Gaba_1 · P7",
  "The parent cell-type label is relevant to extended-amygdala biology.",
  "Pdgfra, Itpr2 and Vcan indicate an OPC/glial program inside a neuronal supertype; exclude from the biological story.",
  "AMY__CS20230722_SUPT_0359", "P6", "Excluded",
  "MEA–BST Nfib P6",
  "MEA–BST Lhx6 Nfib Gaba_3 · P6",
  "This was one of the original VR-2 output programs.",
  "Slc1a2, Slc1a3 and Slc4a4 indicate astrocytic contamination; VR-2 found no shared program across the three proposed GABA outputs."
)

popalign_best <- popalign %>%
  filter(level == "supertype", tissue == "AMY") %>%
  group_by(tissue, label) %>%
  slice_min(p_perm, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    tissue,
    supertype_label = label,
    popalign_component = comp,
    popalign_shift = shift,
    popalign_p_perm = p_perm,
    popalign_q_perm = q_perm,
    popalign_top_feature = top_shift_feature
  )

candidate_evidence <- candidate_spec %>%
  left_join(prog, by = c("unit_id", "program_id")) %>%
  left_join(popalign_best, by = c("tissue", "supertype_label")) %>%
  mutate(
    evidence_tier = case_when(
      role == "Biological anchor" ~ "Exploratory nominal hypothesis",
      role == "Nominal benchmark" ~ "Strongest nominal benchmark; weak annotation",
      TRUE ~ "Excluded from interpretation"
    ),
    dspin_fdr_result = "No program survives BH FDR; q_perm approximately 0.994",
    popalign_context = case_when(
      popalign_p_perm <= 0.05 ~ glue(
        "Nominal unit-level state shift: component {popalign_component}, ",
        "p={fmt_p(popalign_p_perm)}, q={fmt_p(popalign_q_perm)}"
      ),
      TRUE ~ glue(
        "No nominal unit-level state shift: minimum p={fmt_p(popalign_p_perm)}, ",
        "q={fmt_p(popalign_q_perm)}"
      )
    ),
    social_behavior_scope = paste(
      "MeA/BST relevance is anatomical and gene-based;",
      "social behavior was not measured in this experiment."
    )
  ) %>%
  select(
    role, evidence_tier, unit_id, tissue, supertype_label, supertype_name,
    program_id, story_short, story_label, delta_h, hedges_g, separation,
    p_perm, q_perm, n_animals_SPF, n_animals_WildR, top_genes, go_term,
    go_category, go_benjamini, n_mast_hits, mast_genes, mast_sign_agrees,
    popalign_component, popalign_shift, popalign_p_perm, popalign_q_perm,
    popalign_top_feature, popalign_context, rationale, caution,
    dspin_fdr_result, social_behavior_scope
  )

if (nrow(candidate_evidence) != nrow(candidate_spec) || anyNA(candidate_evidence$hedges_g)) {
  stop("Candidate definitions did not join cleanly to DSPIN_program_tests.csv")
}

# ── Option A: complete calibrated landscape ─────────────────────────────────
unit_land <- unit_t %>%
  mutate(
    is_mea_bst = tissue == "AMY" & is_circuit_name(supertype_name),
    unit_short = short_unit(supertype_name, 34),
    label_unit = if_else(
      row_number() %in% order(p_unit)[1:4] |
        (is_mea_bst & p_unit <= 0.10),
      unit_short,
      NA_character_
    )
  )

unit_order <- unit_t %>%
  arrange(tissue, supertype_name) %>%
  transmute(unit_id, unit_order = row_number())

program_land <- prog %>%
  left_join(unit_order, by = "unit_id") %>%
  arrange(tissue, unit_order, as.integer(str_remove(program_id, "^P"))) %>%
  group_by(tissue) %>%
  mutate(test_index = row_number()) %>%
  ungroup() %>%
  mutate(
    is_mea_bst = tissue == "AMY" & is_circuit_name(supertype_name),
    effect_direction = if_else(hedges_g >= 0, "Up in WildR", "Down in WildR")
  )

p_scope <- unit_land %>%
  ggplot(aes(n_cells, median_cells_per_animal, colour = tissue)) +
  geom_point(size = 2.2, alpha = 0.78) +
  geom_point(
    data = ~ filter(.x, is_mea_bst),
    shape = 21, fill = NA, colour = "#0A0A0A", stroke = 0.8, size = 3.2
  ) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = tis_cols) +
  annotate(
    "label", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2,
    label = glue(
      "{global$n_units} units × 20 programs\n",
      "{scales::comma(global$n_programs)} exact tests\n",
      "12 animals (6 SPF / 6 WildR)"
    ),
    size = 3.5, label.size = 0, fill = "white"
  ) +
  labs(
    x = "Nuclei in unit (log scale)",
    y = "Median nuclei per animal (log scale)",
    colour = NULL,
    title = "Breadth and sampling of the fitted cohort",
    subtitle = "Outlined: AMY MeA / COA / CEA / BST-labelled units"
  )

p_programs <- ggplot(program_land, aes(test_index, -log10(p_perm))) +
  geom_point(colour = "#B8B8B8", size = 0.75, alpha = 0.55) +
  geom_point(
    data = ~ filter(.x, is_mea_bst),
    aes(fill = effect_direction),
    shape = 21, colour = "#0A0A0A", stroke = 0.25, size = 1.65, alpha = 0.9
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#767676") +
  geom_hline(yintercept = -log10(6 / 924), linetype = "dotted", colour = "#FF6C0C") +
  facet_wrap(~tissue, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = dir_cols) +
  scale_x_continuous(breaks = NULL) +
  labs(
    x = "All unit × program tests, ordered by cell type",
    y = expression(-log[10](" exact permutation p")),
    fill = NULL,
    title = "All 2,800 program tests",
    subtitle = "MeA/BST-family points are outlined; dashed p=0.05, dotted p=6/924"
  )

p_units <- ggplot(unit_land, aes(rms_hedges_g, -log10(p_unit), colour = tissue)) +
  geom_point(size = 2.2, alpha = 0.7) +
  geom_point(
    data = ~ filter(.x, is_mea_bst),
    shape = 21, fill = NA, colour = "#0A0A0A", stroke = 0.8, size = 3.2
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#767676") +
  geom_text_repel(
    data = ~ filter(.x, !is.na(label_unit)),
    aes(label = label_unit),
    size = 2.7, min.segment.length = 0, max.overlaps = 12,
    box.padding = 0.25, seed = 11, show.legend = FALSE
  ) +
  scale_colour_manual(values = tis_cols) +
  labs(
    x = "RMS Hedges' g across 20 programs",
    y = expression(-log[10](" unit omnibus p")),
    colour = NULL,
    title = "All 140 unit-level omnibus tests",
    subtitle = "No unit survives FDR; labels show the strongest and nominal MeA/BST units"
  )

expected_per_bin <- nrow(prog) * 0.05
p_calibration <- ggplot(prog, aes(p_perm)) +
  geom_histogram(
    breaks = seq(0, 1, 0.05),
    fill = "#9DB3C4", colour = "white", linewidth = 0.25
  ) +
  geom_hline(
    yintercept = expected_per_bin,
    linetype = "dashed", colour = "#FF6C0C", linewidth = 0.8
  ) +
  annotate(
    "label", x = 0.52, y = expected_per_bin * 1.18,
    label = "Expected count per bin under the null",
    colour = "#FF6C0C", size = 3.2, label.size = 0, fill = "white"
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), expand = expansion(mult = 0.01)) +
  labs(
    x = "Exact permutation p",
    y = "Programs",
    title = "The full test set is calibrated to the null",
    subtitle = glue(
      "{global$n_q_prog_lt_0.1} programs and {global$n_q_unit_lt_0.1} units at FDR < 0.1; ",
      "{global$n_leads} nominal leads, a count consistent with chance"
    )
  )

option_a <- (
  p_scope + p_programs + plot_layout(widths = c(0.82, 1.55))
) / (
  p_units + p_calibration + plot_layout(widths = c(1.15, 1))
) +
  plot_annotation(
    title = "D-SPIN tested coordinated program responses across the single-nucleus atlas",
    subtitle = paste(
      "The analysis is broad and animal-replicated, but the complete result set follows the permutation null.",
      "MeA/BST-family units are highlighted as a hypothesis-generating subset."
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 17),
      plot.subtitle = element_text(size = 11, colour = "#555555")
    )
  )

# ── Option B: MeA–BST hypothesis landscape ──────────────────────────────────
mea_programs <- prog %>%
  filter(tissue == "AMY", is_circuit_name(supertype_name)) %>%
  group_by(unit_id) %>%
  mutate(unit_min_p = min(p_perm)) %>%
  ungroup() %>%
  mutate(
    unit_plot = fct_reorder(short_unit(supertype_name, 38), unit_min_p, .fun = min, .desc = TRUE),
    effect_direction = if_else(hedges_g >= 0, "Up in WildR", "Down in WildR"),
    candidate_key = paste(unit_id, program_id, sep = "::")
  )

candidate_keys <- candidate_spec %>%
  filter(role != "Excluded") %>%
  transmute(candidate_key = paste(unit_id, program_id, sep = "::"), role, story_short)

mea_programs <- mea_programs %>%
  left_join(candidate_keys, by = "candidate_key") %>%
  mutate(role = coalesce(role, "Context"))

p_mea_context <- ggplot(mea_programs, aes(hedges_g, unit_plot)) +
  geom_vline(xintercept = 0, colour = "#767676", linewidth = 0.4) +
  geom_point(
    aes(size = -log10(p_perm)),
    colour = "#B8B8B8", alpha = 0.55, position = position_jitter(height = 0.13, seed = 7)
  ) +
  geom_point(
    data = ~ filter(.x, role != "Context"),
    aes(fill = effect_direction, size = -log10(p_perm)),
    shape = 21, colour = "#0A0A0A", stroke = 0.65
  ) +
  geom_text_repel(
    data = ~ filter(.x, role != "Context"),
    aes(label = story_short),
    size = 3, min.segment.length = 0, seed = 13,
    box.padding = 0.3, nudge_y = 0.2, show.legend = FALSE
  ) +
  scale_fill_manual(values = dir_cols) +
  scale_size_continuous(range = c(1, 4.7)) +
  labs(
    x = "Hedges' g (WildR − SPF relative h)",
    y = NULL,
    fill = NULL,
    size = expression(-log[10](p)),
    title = "All programs in AMY MeA / COA / CEA / BST-labelled units",
    subtitle = "Highlighted programs are hypotheses or a nominal benchmark; all other programs remain visible"
  ) +
  theme(axis.text.y = element_text(size = 8.2))

focus_spec <- candidate_spec %>% filter(role != "Excluded")
focus_animals <- animals %>%
  inner_join(
    focus_spec %>% select(unit_id, program_id, role, story_short, story_label),
    by = c("unit_id", "program_id")
  ) %>%
  left_join(
    candidate_evidence %>% select(unit_id, program_id, hedges_g, p_perm, q_perm, separation),
    by = c("unit_id", "program_id")
  ) %>%
  mutate(
    microbiome = factor(microbiome, c("SPF", "WildR")),
    facet_label = glue(
      "{story_label}\n",
      "g={round(hedges_g, 2)} · p={fmt_p(p_perm)} · q={fmt_p(q_perm)} · ",
      "{separation}/6 beyond SPF"
    ) %>% as.character()
  )

p_animals <- ggplot(focus_animals, aes(microbiome, relative_h, colour = microbiome)) +
  geom_hline(yintercept = 0, colour = "#B8B8B8", linewidth = 0.35) +
  geom_quasirandom(width = 0.14, size = 2.7, alpha = 0.9) +
  stat_summary(
    fun = mean, geom = "crossbar", width = 0.5,
    colour = "#0A0A0A", linewidth = 0.4
  ) +
  facet_wrap(~facet_label, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = mb_cols) +
  labs(
    x = NULL,
    y = "SPF-centred program response (relative h)",
    colour = NULL,
    title = "Raw animal-level responses",
    subtitle = "Each point is one animal; crossbars are group means"
  ) +
  theme(strip.text = element_text(face = "plain", size = 8.2, hjust = 0))

weights_file <- "gene_programs_total_20_onmf_20_prior_0_weights.csv"
read_program_weights <- function(unit_id, program_id, story_short) {
  k <- str_remove(program_id, "^P")
  f <- file.path(dspin_dir, unit_id, "dspin_run", weights_file)
  if (!file.exists(f)) stop("Missing oNMF weight table: ", f)
  x <- read_csv(f, show_col_types = FALSE)
  tibble(
    story_short,
    program_id,
    gene = x[[k]],
    weight = as.numeric(x[[paste0("Weights_", k)]])
  ) %>%
    filter(!is.na(gene), nzchar(gene), !is.na(weight))
}

weights_long <- pmap_dfr(
  focus_spec %>%
    filter(role == "Biological anchor") %>%
    select(unit_id, program_id, story_short),
  read_program_weights
)

story_genes <- c(
  "Ar", "Esr1", "Htr7", "Cntnap5a", "Pam", "Thsd7b", "Grin3a", "Cacna2d1",
  "Pcdh7", "Fgf13", "Cdh13", "Prkcb", "Thrb", "Cdh12"
)

weights_plot <- weights_long %>%
  group_by(story_short) %>%
  slice_max(weight, n = 12, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(story_short, weight) %>%
  mutate(
    story_gene = factor(
      paste(story_short, gene, sep = "___"),
      levels = unique(paste(story_short, gene, sep = "___"))
    ),
    highlighted_gene = gene %in% story_genes
  )

p_weights <- ggplot(weights_plot, aes(weight, story_gene, fill = highlighted_gene)) +
  geom_col(width = 0.72) +
  facet_wrap(~story_short, scales = "free_y", ncol = 1) +
  scale_y_discrete(labels = function(x) str_remove(x, "^.*___")) +
  scale_fill_manual(values = c(`TRUE` = "#FF6C0C", `FALSE` = "#C8D3DB")) +
  labs(
    x = "oNMF weight",
    y = NULL,
    fill = NULL,
    title = "Highest-weight program genes",
    subtitle = "Orange genes anchor the proposed biological labels"
  ) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 8.2),
    strip.text = element_text(hjust = 0)
  )

option_b <- p_mea_context / (
  p_animals + p_weights + plot_layout(widths = c(1.55, 0.9))
) +
  plot_layout(heights = c(0.92, 1.45)) +
  plot_annotation(
    title = "MeA–BST programs nominate independent social-behavior-related hypotheses",
    subtitle = paste(
      "These are anatomical and gene-based hypotheses, not FDR-significant findings or a coupled circuit.",
      "Social behavior was not measured."
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 17),
      plot.subtitle = element_text(size = 11, colour = "#555555")
    )
  )

# ── Write once, after every plot and table has been assembled ────────────────
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

candidate_f <- file.path(tbl_dir, "DSPIN_mea_bst_candidate_evidence.csv")
registry_f <- file.path(tbl_dir, "DSPIN_paper_figure_source_registry.csv")
write_csv(candidate_evidence, candidate_f)
write_csv(
  enframe(inputs, name = "input", value = "path") %>%
    mutate(
      rows = map_int(path, ~ nrow(read_csv(.x, show_col_types = FALSE))),
      modified = as.character(file.info(path)$mtime)
    ),
  registry_f
)

save_plot_pair <- function(plot, stem, width, height) {
  png_f <- file.path(fig_dir, paste0(stem, ".png"))
  pdf_f <- file.path(fig_dir, paste0(stem, ".pdf"))
  if (file.exists(png_f) || file.exists(pdf_f)) stop("Refusing to overwrite: ", stem)
  ggsave(png_f, plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(
    pdf_f, plot, width = width, height = height,
    device = grDevices::cairo_pdf, bg = "white"
  )
  message("wrote ", png_f)
  message("wrote ", pdf_f)
}

save_plot_pair(option_a, "DSPIN_paper_option_A_calibrated_landscape", 16, 11)
save_plot_pair(option_b, "DSPIN_paper_option_B_mea_bst_hypotheses", 16, 13)

message("Candidate evidence: ", candidate_f)
message("Source registry: ", registry_f)
message("Done: ", fig_dir)
