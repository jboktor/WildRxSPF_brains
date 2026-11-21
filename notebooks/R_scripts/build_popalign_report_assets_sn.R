#!/usr/bin/env Rscript
# Report assets for 0X_popAlign_SN.qmd, built from a finished popAlign synthesis.
# Light (reads saved tables only; seconds). Writes into the synthesis' own
# tables/figures folders:
#   popAlign_followup_leads.csv        curated leads (robust, uncontaminated components)
#   popAlign_meta_program_summary.csv  one row per meta-program (size, levels, top genes)
#   popAlign_meta_programs_top.png     spaced heatmap: top meta-programs × units
#
# Lead rule (pre-specified, no FDR survivor exists at 6 vs 6 animals):
#   robust & !contamination & (state p <= 0.0065  [≤ 6/924]
#                               | abundance p < 0.05, non-redundant, ≥5/6 animals
#                                 separated, and depth-adjusted GLM p < 0.05)
#   Ranking: replicated lineages first, then lineages whose hits are all depth-clean, then best p.
#   "convergent" = the same Allen lineage (name without numeric prefix / _N suffix)
#   in the same tissue yields leads in ≥ 2 distinct units (e.g. subclass + its supertype).
#
# Usage:
#   Rscript notebooks/R_scripts/build_popalign_report_assets_sn.R --run-tag 2026-10-04 --synth-tag 2026-10-05

suppressPackageStartupMessages({ library(tidyverse); library(glue) })
source("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/notebooks/R_scripts/popalign_sn_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[[i + 1]] }
P <- pa_paths(get_arg("--run-tag", "2026-10-04"))
sfx <- glue("_synth-{get_arg('--synth-tag', '2026-10-05')}")
tbl <- paste0(P$tbl, sfx); fig <- paste0(P$fig, sfx)
stopifnot(dir.exists(tbl), dir.exists(fig))

cat_ <- read_csv(file.path(tbl, "popAlign_subpopulation_catalog.csv"), show_col_types = FALSE)
man <- read_csv(P$manifest, show_col_types = FALSE)
gmap <- load_gene_map(P)
lk_nm <- setNames(gmap$gene_name, gmap$export_symbol)
lineage_key <- function(x) x %>% str_remove("^\\d+ ") %>% str_remove("_\\d+$")
first_n <- function(x, n) map_chr(str_split(coalesce(x, ""), ",\\s*|\\s*\\|\\s*"), ~ paste(head(.x[nzchar(.x)], n), collapse = ", "))

# ── follow-up leads ──────────────────────────────────────────────────────────
ok <- cat_ %>% filter(robust, !contamination_flag)
state_leads <- ok %>% filter(p_state <= 0.0065) %>% mutate(effect = "state shift")
abund_leads <- ok %>%
  filter(p_abund < 0.05, !coalesce(redundant, FALSE), separation >= 5, p_abund_glm_depthadj < 0.05) %>%
  mutate(effect = "abundance shift")
leads <- bind_rows(state_leads, abund_leads)

# translate the oNMF program that moves most (state leads) into its top genes
prog_genes <- leads %>% filter(effect == "state shift") %>% distinct(unit_id, top_shift_feature) %>%
  left_join(man %>% dplyr::select(unit_id, output_dir), by = "unit_id") %>%
  mutate(genes = map2_chr(output_dir, top_shift_feature, function(d, f) {
    tg <- read_csv(file.path(d, "onmf_top_genes.csv"), show_col_types = FALSE) %>% filter(feature == f) %>% arrange(rank)
    paste(head(coalesce(unname(lk_nm[tg$gene]), tg$gene), 6), collapse = ", ")
  })) %>% dplyr::select(unit_id, top_shift_feature, shifted_program_genes = genes)

leads <- leads %>%
  left_join(prog_genes, by = c("unit_id", "top_shift_feature")) %>%
  mutate(lineage = lineage_key(name)) %>%
  group_by(tissue, lineage) %>%
  mutate(n_hits_lineage = n_distinct(unit_id), levels_hit = paste(sort(unique(level)), collapse = " + ")) %>%
  ungroup() %>%
  mutate(convergent = n_hits_lineage >= 2,          # ≥ 2 distinct units, not 2 comps of one unit
         separation = if_else(effect == "state shift", NA_real_, separation),
         p = if_else(effect == "state shift", p_state, p_abund),
         q = if_else(effect == "state shift", q_state, q_abund),
         direction = case_when(effect == "abundance shift" & log2fc_frac > 0 ~ "↑ in WildR",
                               effect == "abundance shift" ~ "↓ in WildR",
                               TRUE ~ "state moved"),
         share_pct = round(100 * atlas_weight),
         markers = first_n(top_markers, 6),
         go = first_n(top_go, 1)) %>%
  # rank lineages: replicated across units > all hits depth-clean > best p; rows within lineage by p
  group_by(tissue, lineage) %>%
  mutate(lineage_clean = !any(depth_flag), lineage_best_p = min(p)) %>%
  ungroup() %>%
  arrange(desc(convergent), desc(lineage_clean), lineage_best_p, tissue, lineage, depth_flag, p, desc(stability)) %>%
  dplyr::select(tissue, lineage, level, name, comp, effect, direction, p, q, separation, depth_flag, lineage_clean,
                mean_frac_SPF, mean_frac_WildR, p_abund_glm_depthadj, share_pct, stability, purity,
                shifted_program_genes, markers, go, convergent, levels_hit, n_hits_lineage, unit_id)
write_csv(leads, file.path(tbl, "popAlign_followup_leads.csv"))
message("Leads: ", nrow(leads), " (", sum(leads$convergent), " in convergent lineages)")

# ── meta-programs: summary + spaced heatmap ─────────────────────────────────
mm <- read_csv(file.path(tbl, "popAlign_meta_program_membership.csv"), show_col_types = FALSE) %>%
  filter(!is.na(meta_program)) %>%
  left_join(cat_ %>% dplyr::select(unit_id, comp, level, tissue, name, contamination_flag), by = c("unit_id", "comp"))
mg <- read_csv(file.path(tbl, "popAlign_meta_program_genes.csv"), show_col_types = FALSE)
mp_sum <- mm %>%
  group_by(meta_program) %>%
  summarise(n_components = n(), n_units = n_distinct(unit_id),
            n_class = sum(level == "class"), n_subclass = sum(level == "subclass"), n_supertype = sum(level == "supertype"),
            tissues = paste(sort(unique(tissue)), collapse = "+"),
            frac_contaminated = mean(contamination_flag), .groups = "drop") %>%
  left_join(mg %>% group_by(meta_program) %>% arrange(desc(frac_members), .by_group = TRUE) %>%
              summarise(top_genes = paste(head(gene_name, 8), collapse = ", ")), by = "meta_program") %>%
  arrange(desc(n_units))
write_csv(mp_sum, file.path(tbl, "popAlign_meta_program_summary.csv"))

top_mp <- mp_sum %>% filter(n_units >= 4) %>% slice_head(n = 20)
lab_mp <- top_mp %>% mutate(lab = glue("MP{meta_program}  ·  {first_n(top_genes, 4)}")) %>%
  dplyr::select(meta_program, lab)
hm <- mm %>% semi_join(top_mp, by = "meta_program") %>%
  distinct(meta_program, unit_id, level, tissue, name) %>%
  left_join(lab_mp, by = "meta_program") %>%
  mutate(lab = factor(lab, levels = rev(lab_mp$lab)),
         level = factor(level, levels = pa_levels),
         unit = glue("{str_trunc(name, 30)} ({substr(level, 1, 3)})"))
# order units by tissue, then by the first meta-program they belong to (block-diagonal look)
ord <- hm %>% group_by(unit, tissue) %>% summarise(first_mp = max(as.integer(lab)), .groups = "drop") %>%
  arrange(tissue, desc(first_mp)) %>% pull(unit)
hm <- hm %>% mutate(unit = factor(unit, levels = unique(ord)))
n_units_shown <- n_distinct(hm$unit)
# one panel per tissue, stacked, panel height proportional to its number of programs
suppressPackageStartupMessages(library(patchwork))
hm_panel <- function(d, tis, show_legend) {
  d <- d %>% filter(tissue == tis) %>% mutate(lab = droplevels(lab), unit = droplevels(unit))
  ggplot(d, aes(unit, lab, fill = level)) +
    geom_tile(colour = "white", linewidth = 1.1, width = 0.92, height = 0.86) +
    scale_fill_manual(values = c(class = "#27455C", subclass = "#FF9240", supertype = "#FF6C0C"),
                      name = "Allen level", drop = FALSE) +
    labs(x = NULL, y = NULL, title = tis) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8.5),
          axis.text.y = element_text(size = 11.5, colour = "#0A0A0A"),
          panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 14, colour = "#FF6C0C"),
          legend.position = if (show_legend) "top" else "none")
}
tissues <- sort(unique(hm$tissue))
n_rows <- map_int(tissues, ~ n_distinct(hm$lab[hm$tissue == .x]))
p <- wrap_plots(map2(tissues, seq_along(tissues), ~ hm_panel(hm, .x, .y == 1)), ncol = 1, heights = n_rows) +
  plot_annotation(
    title = glue("Top {nrow(top_mp)} recurrent meta-programs (each present in ≥ 4 units)"),
    subtitle = "Columns = units containing a member component · rows = meta-program ID · most-shared marker genes",
    theme = theme(plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 12)))
ggsave(file.path(fig, "popAlign_meta_programs_top.png"), p,
       width = max(14, 0.24 * max(table(distinct(hm, unit, tissue)$tissue)) + 7),
       height = 0.42 * sum(n_rows) + 2 * 3.0 + 1.2, dpi = 200, limitsize = FALSE, bg = "white")
message("Meta-program heatmap: ", nrow(top_mp), " programs × ", n_units_shown, " units")

# ── circuit enrichment of state shifts (unit level) ─────────────────────────
# Hypothesis from the VR relay analyses: the medial amygdala → BST → medial
# hypothalamus network. Defined once by Allen region names, not tuned:
circuit_pattern <- "MEA|BST|PVH|DMH|AHN|TMd|MPN|MPO|VMH|SBPV"
# A unit counts once ("any robust, uncontaminated component with state p < 0.05"),
# so components of one unit are not treated as independent. Tested per level as
# well, because subclass and supertype units share nuclei.
rc <- cat_ %>% filter(robust, !contamination_flag, !is.na(p_state))
units <- rc %>% group_by(unit_id, level, tissue, name) %>%
  summarise(hit = any(p_state < 0.05), strong = any(p_state <= 0.0065), .groups = "drop") %>%
  mutate(circuit = str_detect(name, circuit_pattern),
         transmitter = case_when(str_detect(name, "Gaba|GABA") ~ "GABA", str_detect(name, "Glut") ~ "Glut", TRUE ~ "non-neuronal"))
fisher_row <- function(d, what, by, scope) {
  a <- d[[by]]; h <- d[[what]]
  tab <- matrix(c(sum(a & h), sum(a & !h), sum(!a & h), sum(!a & !h)), 2, byrow = TRUE)
  ft <- fisher.test(tab)
  tibble(scope = scope, outcome = what, grouping = by,
         in_group_hits = tab[1, 1], in_group_units = sum(a), out_group_hits = tab[2, 1], out_group_units = sum(!a),
         rate_in = tab[1, 1] / sum(a), rate_out = tab[2, 1] / sum(!a),
         odds_ratio = unname(ft$estimate), p = ft$p.value)
}
units_hy <- units %>% mutate(hyp_gaba = tissue == "HYP" & transmitter == "GABA", gaba = transmitter == "GABA")
enr <- bind_rows(
  fisher_row(units_hy, "hit", "circuit", "all levels"),
  fisher_row(units_hy, "strong", "circuit", "all levels"),
  fisher_row(filter(units_hy, level == "subclass"), "hit", "circuit", "subclass only"),
  fisher_row(filter(units_hy, level == "supertype"), "hit", "circuit", "supertype only"),
  fisher_row(filter(units_hy, tissue == "AMY"), "hit", "circuit", "AMY only"),
  fisher_row(filter(units_hy, tissue == "HYP"), "hit", "circuit", "HYP only"),
  fisher_row(units_hy, "hit", "gaba", "all levels (GABA vs rest)"),
  fisher_row(units_hy, "hit", "hyp_gaba", "all levels (HYP GABA vs rest)"))
# confound check: "any component shifted" is easier to hit with more components, so
# adjust for component count, unit size, level and tissue (unit level), and repeat at
# component level.
u_cov <- rc %>% group_by(unit_id, level, tissue, name) %>%
  summarise(hit = any(p_state < 0.05), n_comp = n(), n_cells = sum(n_cells_hard), .groups = "drop") %>%
  mutate(circuit = str_detect(name, circuit_pattern))
m_unit <- glm(hit ~ circuit + log(n_comp) + log(n_cells) + level + tissue, family = binomial, data = u_cov)
m_comp <- glm(I(p_state < 0.05) ~ circuit + log(n_cells_hard) + level + tissue, family = binomial,
              data = rc %>% mutate(circuit = str_detect(name, circuit_pattern)))
adj <- tibble(scope = c("unit level, adjusted (n_comp, n_cells, level, tissue)", "component level, adjusted (n_cells, level, tissue)"),
              outcome = "hit", grouping = "circuit",
              odds_ratio = exp(c(coef(m_unit)[["circuitTRUE"]], coef(m_comp)[["circuitTRUE"]])),
              p = c(coef(summary(m_unit))["circuitTRUE", 4], coef(summary(m_comp))["circuitTRUE", 4]))
enr <- bind_rows(enr, adj)
write_csv(enr, file.path(tbl, "popAlign_circuit_enrichment.csv"))
write_csv(units %>% arrange(desc(circuit), desc(hit)), file.path(tbl, "popAlign_circuit_units.csv"))

# bars = the five raw comparisons only (adjusted-model rows have no counts; shown in the subtitle)
scopes_plot <- c("all levels", "subclass only", "supertype only", "AMY only", "HYP only")
adj_u <- enr %>% filter(scope == "unit level, adjusted (n_comp, n_cells, level, tissue)")
pe <- enr %>% filter(outcome == "hit", grouping == "circuit", scope %in% scopes_plot) %>%
  mutate(scope = factor(scope, levels = rev(c("all levels", "subclass only", "supertype only", "AMY only", "HYP only")))) %>%
  pivot_longer(c(rate_in, rate_out), names_to = "grp", values_to = "rate") %>%
  mutate(grp = if_else(grp == "rate_in", "MeA / BST / medial hypothalamus", "all other units"),
         n_lab = if_else(grp == "MeA / BST / medial hypothalamus", glue("{in_group_hits}/{in_group_units}"), glue("{out_group_hits}/{out_group_units}")),
         p_lab = glue("OR {sprintf('%.1f', odds_ratio)}, p = {signif(p, 2)}"))
p_enr <- ggplot(pe, aes(rate, scope, fill = grp)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = n_lab), position = position_dodge(width = 0.75), hjust = -0.15, size = 3.6) +
  geom_text(data = distinct(pe, scope, p_lab), aes(x = 0.62, y = scope, label = p_lab), inherit.aes = FALSE,
            hjust = 1, size = 3.8, fontface = "bold", colour = "#0A0A0A") +
  scale_x_continuous(labels = scales::percent, limits = c(0, 0.62), expand = c(0, 0)) +
  scale_fill_manual(values = c(`MeA / BST / medial hypothalamus` = "#FF6C0C", `all other units` = "#A5A5A5"), name = NULL) +
  labs(x = "Units with a nominal microbiome state shift (p < 0.05)", y = NULL,
       title = "State shifts concentrate in the MeA → BST → medial hypothalamus network",
       subtitle = glue("Each unit counted once · Fisher's exact test per row · robust, uncontaminated components only\n",
                       "Adjusted for subpopulation count, nuclei, level and tissue: OR {sprintf('%.1f', adj_u$odds_ratio)}, p = {signif(adj_u$p, 2)}")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top", panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold"),
        plot.title.position = "plot")
ggsave(file.path(fig, "popAlign_circuit_enrichment.png"), p_enr, width = 11, height = 5.2, dpi = 200, bg = "white")

# ── do DEG-rich units change composition? ────────────────────────────────────
om <- read_csv(file.path(tbl, "popAlign_abundance_omnibus.csv"), show_col_types = FALSE)
deg <- map_dfr(names(P$mast), function(lv) {
  if (!file.exists(P$mast[[lv]])) return(NULL)
  readRDS(P$mast[[lv]]) %>% as_tibble() %>% filter(p_adjust < 0.1) %>%
    count(tissue, label = as.character(cell_type), name = "n_deg") %>% mutate(level = lv)
})
dc <- om %>% filter(!is.na(p_omnibus)) %>% left_join(deg, by = c("tissue", "label", "level")) %>%
  mutate(n_deg = replace_na(n_deg, 0L))
write_csv(dc, file.path(tbl, "popAlign_deg_vs_composition.csv"))
sp <- cor.test(log1p(dc$n_deg), -log10(dc$p_omnibus), method = "spearman", exact = FALSE)
p_dc <- ggplot(dc, aes(n_deg + 1, -log10(p_omnibus), colour = level)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
  geom_point(alpha = 0.75, size = 2.3) + scale_x_log10() +
  scale_colour_manual(values = c(class = "#27455C", subclass = "#FF9240", supertype = "#FF6C0C"), name = "Allen level") +
  labs(x = "MAST DEGs in the unit (padj < 0.1, + 1, log scale)", y = "Composition shift, −log10 omnibus p",
       title = "Units with many DEGs do not change their subpopulation mix",
       subtitle = glue("Spearman ρ = {sprintf('%.2f', sp$estimate)} (p = {signif(sp$p.value, 2)}); dashed line p = 0.05")) +
  theme_minimal(base_size = 13) + theme(legend.position = "top", plot.title = element_text(face = "bold"), plot.title.position = "plot")
ggsave(file.path(fig, "popAlign_deg_vs_composition.png"), p_dc, width = 9, height = 5.5, dpi = 200, bg = "white")
message("Circuit enrichment + DEG-vs-composition written")
