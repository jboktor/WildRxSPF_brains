#!/usr/bin/env Rscript
# Synthesis of per-unit popAlign runs (snRNA-seq; class, subclass, supertype).
# Reads data/interim/popalign/snRNASeq/results_<run-tag>/<level>/<tissue>__<label>/,
# writes tables → tables/popAlign/snRNASeq/run_<tag>/, figures → figures/popAlign/snRNASeq/run_<tag>/,
# heavier intermediates → data/interim/popalign/snRNASeq/synthesis_<tag>/.
#
# Usage:
#   Rscript notebooks/R_scripts/run_popalign_synthesis_sn.R [--run-tag 2026-10-04] [--synth-tag 2026-10-05] [--skip-go]
#   sbatch notebooks/shell_scripts/run_popalign_synthesis_sn.sh [--run-tag ...] [--skip-go]
# Test overrides: --manifest <csv> --out <dir>

suppressPackageStartupMessages({
  library(tidyverse); library(glue); library(patchwork); library(ggrepel)
})
source("/resnick/groups/mthomson/jboktor/WILDRxSPF_brains/notebooks/R_scripts/popalign_sn_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[[i + 1]] }
skip_go <- "--skip-go" %in% args
no_go <- "--no-go" %in% args   # testing only: skip GO entirely (empty GO table)
P <- pa_paths(get_arg("--run-tag", "2026-10-04"))
# --synth-tag: version the synthesis outputs (tables / figures / synthesis cache) so a
# re-synthesis of the same popAlign run never overwrites an earlier one
if (!is.null(get_arg("--synth-tag"))) {
  sfx <- glue("_synth-{get_arg('--synth-tag')}")
  P$synth <- paste0(P$synth, sfx); P$fig <- paste0(P$fig, sfx); P$tbl <- paste0(P$tbl, sfx)
}
if (!is.null(get_arg("--manifest"))) P$manifest <- get_arg("--manifest")
if (!is.null(get_arg("--out"))) { o <- get_arg("--out"); P$synth <- file.path(o, "synthesis"); P$fig <- file.path(o, "figures"); P$tbl <- file.path(o, "tables") }
walk(c(P$synth, P$fig, P$tbl), dir.create, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "-", paste0(...), "\n")
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                           strip.background = element_rect(fill = "#FFF3EB", colour = NA)))
save_fig <- function(p, name, w, h) ggsave(file.path(P$fig, glue("popAlign_{name}.png")), p, width = w, height = h, dpi = 200)
lvl <- function(x) factor(x, levels = pa_levels)

man <- read_csv(P$manifest, show_col_types = FALSE)
st <- read_status(man)
write_csv(st %>% dplyr::select(array_index, unit_id, level, tissue, name, n_cells_exported, status, beta, atlas_K,
                               atlas_pool, boot_K_agreement, n_features, n_features_at_bound,
                               atlas_K_at_upper_bound, plot_deltas, runtime_min),
          file.path(P$tbl, "popAlign_unit_status.csv"))
done <- st %>% filter(status == "ok")
log_msg("Units done: ", nrow(done), "/", nrow(man))
stopifnot(nrow(done) > 0)

# ── load ─────────────────────────────────────────────────────────────────────
comps <- read_unit_csv(done, "atlas_components.csv")
fracs <- read_unit_csv(done, "animal_component_fractions.csv")
cents <- read_unit_csv(done, "animal_component_centroids.csv")
markers <- read_unit_csv(done, "atlas_markers.csv")
al <- read_unit_csv(done, "alignments.csv")
universe_by_unit <- done %>%
  mutate(g = map(input_dir, ~ read_tsv(file.path(.x, "genes.tsv"), col_names = FALSE, show_col_types = FALSE)$X1)) %>%
  { setNames(.$g, .$unit_id) }
# gene identity from the Seurat feature metadata: export symbols → Ensembl IDs
log_msg("Gene map: ", P$gene_meta)
gmap <- load_gene_map(P)
lk_id <- setNames(gmap$gene_id, gmap$export_symbol)
lk_nm <- setNames(gmap$gene_name, gmap$export_symbol)
universe_ens <- map(universe_by_unit, ~ unname(na.omit(lk_id[.x])))
log_msg("Unmapped universe genes (median per unit): ",
        median(map2_int(universe_by_unit, universe_ens, ~ length(.x) - length(.y))))

# ── statistics ───────────────────────────────────────────────────────────────
log_msg("Per-animal depth"); adepth <- animal_depth(done)
log_msg("Abundance tests"); ab <- test_abundance(fracs, adepth)
abund <- ab$comp; abund_unit <- ab$unit
log_msg("State-shift tests"); state <- test_state_shift(cents)
log_msg("Alignment summaries"); align_sum <- summarise_alignments(al)
log_msg("Allen composition"); allen <- allen_composition(done)
log_msg("Depth QC"); depth <- depth_qc(done)

sets <- component_marker_sets(markers) %>% annotate_sets(gmap)
log_msg("Marker sets: ", nrow(sets), " | contamination flags")
contam <- contamination_flags(sets, done)
log_msg("D-SPIN / MAST overlap")
dsp <- dspin_overlap(sets, done, P)
mst <- mast_overlap(sets, done, universe_by_unit, P)
go_cache <- file.path(P$synth, "go_components.rds")
go <- if (no_go) go_empty else if (skip_go && file.exists(go_cache)) readRDS(go_cache) else {
  log_msg("GO enrichment (Ensembl IDs)"); go_enrich_components(sets, universe_ens) }
if (!no_go) saveRDS(go, go_cache)
log_msg("Meta-programs"); mp <- meta_programs(sets)
saveRDS(mp, file.path(P$synth, "meta_programs.rds"))
log_msg("Cross-level maps")
xl <- bind_rows(cross_level_map(done, "class", "subclass"), cross_level_map(done, "subclass", "supertype"))
write_csv(xl, file.path(P$tbl, "popAlign_cross_level_map.csv"))

# ── catalog: one row per discovered subpopulation ────────────────────────────
top_go <- go %>% filter(!is.na(p.adjust)) %>% group_by(unit_id, comp) %>%
  slice_min(p.adjust, n = 3, with_ties = FALSE) %>%
  summarise(top_go = paste(Description, collapse = " | "), go_min_padj = min(p.adjust), .groups = "drop")
top_mk <- sets %>% mutate(top_markers = map_chr(names, ~ paste(head(.x, 10), collapse = ", "))) %>%
  dplyr::select(unit_id, comp, n_markers, top_markers)
unitK <- done %>% dplyr::select(unit_id, atlas_K, boot_K_agreement, beta)

catalog <- comps %>%
  dplyr::select(unit_id, level, tissue, label, name, comp, atlas_weight, stability,
                frac_boots_jaccard_ge_075, n_cells_hard) %>%
  left_join(unitK, by = "unit_id") %>%
  left_join(allen %>% dplyr::select(unit_id, comp, allen_level, dominant_allen, purity, entropy_bits, n_allen_10pct),
            by = c("unit_id", "comp")) %>%
  left_join(depth %>% dplyr::select(unit_id, comp, median_log10_umi, delta_log10_umi, unit_depth_eta2),
            by = c("unit_id", "comp")) %>%
  left_join(top_mk, by = c("unit_id", "comp")) %>%
  left_join(top_go, by = c("unit_id", "comp")) %>%
  left_join(abund %>% dplyr::select(unit_id, comp, mean_frac_SPF, mean_frac_WildR, log2fc_frac, g_abund = g,
                                    separation, redundant, p_abund = p_perm, q_abund = q_perm,
                                    p_abund_glm = p_glm, q_abund_glm = q_glm, p_abund_glm_depthadj = p_glm_depth,
                                    p_omnibus, q_omnibus),
            by = c("unit_id", "comp")) %>%
  left_join(state %>% dplyr::select(unit_id, comp, shift, top_shift_feature, p_state = p_perm, q_state = q_perm),
            by = c("unit_id", "comp")) %>%
  left_join(align_sum %>% dplyr::select(unit_id, comp, present_SPF, present_WildR, log10_jd_SPF, log10_jd_WildR, p_perm_jd),
            by = c("unit_id", "comp")) %>%
  left_join(dsp, by = c("unit_id", "comp")) %>%
  left_join(mst %>% dplyr::select(unit_id, comp, n_markers_mast, p_mast_enrich, mast_genes), by = c("unit_id", "comp")) %>%
  left_join(mp$membership %>% dplyr::select(unit_id, comp, meta_program), by = c("unit_id", "comp")) %>%
  left_join(contam, by = c("unit_id", "comp")) %>%
  mutate(
    contamination_flag = coalesce(contamination_flag, FALSE),
    substructure = atlas_K > 1,
    # |Δ median log10 UMI| >= 0.1 (~26% depth difference) → possibly technical
    depth_flag = coalesce(abs(delta_log10_umi) >= 0.1, FALSE),
    # robustness is about the component being real; depth is judged on the effect below
    robust = substructure & stability >= 0.75 & atlas_weight >= 0.02 & coalesce(n_markers, 0L) >= 10,
    kind = case_when(!substructure ~ "no substructure (K = 1)",
                     contamination_flag ~ glue("possible contamination ({contamination_type})"),
                     purity >= 0.7 ~ "Allen-aligned subtype",
                     purity < 0.5 ~ "cross-type state",
                     TRUE ~ "mixed"),
    abundance_hit = coalesce(q_abund < 0.1, FALSE) | coalesce(separation == 6 & p_abund < 0.05, FALSE),
    state_hit = coalesce(q_state < 0.1, FALSE),
    # a depth-associated component's abundance hit must survive the depth covariate
    depth_survives = !depth_flag | coalesce(p_abund_glm_depthadj < 0.05, FALSE),
    microbiome_effect = case_when(
      (abundance_hit | state_hit) & !substructure & state_hit ~ "state (whole population)",
      (abundance_hit | state_hit) & !robust ~ "hit, component not robust",
      (abundance_hit | state_hit) & contamination_flag ~ "excluded (possible contamination)",
      abundance_hit & !depth_survives ~ "excluded (depth-confounded)",
      abundance_hit & state_hit ~ "abundance + state",
      abundance_hit ~ "abundance", state_hit ~ "state",
      coalesce(p_abund < 0.05, FALSE) | coalesce(p_state < 0.05, FALSE) ~ "nominal",
      TRUE ~ "none"),
    interpretability = (robust & !contamination_flag) * 1 + !is.na(top_go) * 1 + (coalesce(n_markers, 0L) >= 25) * 0.5 +
      (!is.na(meta_program)) * 0.5 + (coalesce(dspin_jaccard, 0) >= 0.1) * 0.5
  ) %>%
  arrange(lvl(level), dplyr::desc(interpretability), q_abund)

write_csv(catalog, file.path(P$tbl, "popAlign_subpopulation_catalog.csv"))
write_csv(abund, file.path(P$tbl, "popAlign_abundance_tests.csv"))
write_csv(abund_unit, file.path(P$tbl, "popAlign_abundance_omnibus.csv"))
write_csv(state, file.path(P$tbl, "popAlign_state_shift_tests.csv"))
write_csv(depth, file.path(P$tbl, "popAlign_depth_qc.csv"))
write_csv(go, file.path(P$tbl, "popAlign_component_GO.csv"))
mp$genes <- mp$genes %>% mutate(gene_name = coalesce(unname(lk_nm[gene]), gene), gene_id = unname(lk_id[gene])) %>%
  dplyr::rename(export_symbol = gene) %>% relocate(gene_name, gene_id, .after = n_members)
write_csv(mp$genes, file.path(P$tbl, "popAlign_meta_program_genes.csv"))
write_csv(contam, file.path(P$tbl, "popAlign_contamination_flags.csv"))
write_csv(mp$membership %>% dplyr::select(unit_id, comp, meta_program, n_markers), file.path(P$tbl, "popAlign_meta_program_membership.csv"))
log_msg("Catalog rows: ", nrow(catalog))

# ── figures ──────────────────────────────────────────────────────────────────
tissue_cols <- c(AMY = "#FF6C0C", HYP = "#27455C")
effect_cols <- c(`abundance + state` = "#B33C00", abundance = "#FF6C0C", state = "#FF9240",
                 `state (whole population)` = "#FFC08A", `hit, component not robust` = "#767676",
                 `excluded (depth-confounded)` = "#27455C", `excluded (possible contamination)` = "#555555",
                 nominal = "#A5A5A5", none = "#E5E5E5")

# 1. model overview: K per unit (K = 1 = no substructure) and bootstrap K agreement
p1 <- done %>% ggplot(aes(n_cells_exported, atlas_K, colour = tissue)) +
  geom_jitter(width = 0, height = 0.15, alpha = 0.8) + scale_x_log10() +
  facet_wrap(~ lvl(level), scales = "free") + scale_colour_manual(values = tissue_cols) +
  labs(x = "Nuclei in unit (log)", y = "Atlas components (K, held-out BIC)", title = "Subpopulations per unit")
p2 <- done %>% ggplot(aes(lvl(level), boot_K_agreement, colour = tissue)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.7) + scale_colour_manual(values = tissue_cols) +
  labs(x = NULL, y = "Bootstraps re-selecting the same K", title = "Is K itself stable?") +
  theme(legend.position = "none")
save_fig(p1 + p2 + plot_layout(widths = c(3, 1)), "unit_overview", 13, 4.5)

# 2. abundance volcano (non-redundant component tests)
ab_plot <- abund %>% filter(!redundant) %>% mutate(level = lvl(level))
p <- ab_plot %>% ggplot(aes(delta_logit, -log10(p_perm), colour = q_perm < 0.1)) +
  geom_jitter(alpha = 0.7, width = 0, height = 0.03) + facet_wrap(~ level) +
  geom_hline(yintercept = -log10(2 / 924), linetype = 3, colour = "grey50") +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#FF6C0C"), name = "q < 0.1", na.value = "grey70") +
  geom_text_repel(data = ab_plot %>% group_by(level) %>% slice_min(p_perm, n = 6, with_ties = FALSE),
                  aes(label = glue("{tissue} {str_trunc(name, 20)} c{comp}")),
                  size = 2.5, colour = "black", max.overlaps = 30) +
  labs(x = "Δ logit fraction (WildR − SPF), animal level", y = "−log10 exact permutation p",
       title = "Subpopulation abundance shifts", caption = "Dotted line: smallest attainable exact p (2/924)")
save_fig(p, "abundance_volcano", 13, 5)

# 3. per-animal fractions for top abundance hits
top_hits <- ab_plot %>% arrange(p_perm, dplyr::desc(abs(g))) %>% slice_head(n = 8)
if (nrow(top_hits)) {
  p <- fracs %>% semi_join(top_hits, by = c("unit_id", "comp")) %>%
    mutate(panel = glue("{level} · {tissue} · {str_trunc(name, 22)} · c{comp}")) %>%
    ggplot(aes(microbiome, soft_frac, colour = microbiome)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) + geom_jitter(width = 0.1, size = 2) +
    facet_wrap(~ panel, scales = "free_y", ncol = 4) + scale_colour_manual(values = pa_colors) +
    labs(x = NULL, y = "Fraction of unit nuclei", title = "Top abundance shifts (one point = one animal)") +
    theme(legend.position = "none")
  save_fig(p, "top_abundance_animals", 13, 6.5)
}

# 4. state shift
p <- state %>% mutate(level = lvl(level)) %>%
  ggplot(aes(shift, -log10(p_perm), colour = q_perm < 0.1)) +
  geom_point(alpha = 0.7) + facet_wrap(~ level) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#FF6C0C"), name = "q < 0.1") +
  labs(x = "Standardised centroid shift (oNMF space)", y = "−log10 exact permutation p",
       title = "Within-subpopulation expression-state shifts")
save_fig(p, "state_shift", 12, 4.5)

# 5. interpretability landscape (units with substructure)
p <- catalog %>% filter(substructure, !is.na(purity)) %>% mutate(level = lvl(level)) %>%
  ggplot(aes(purity, stability, size = atlas_weight, colour = microbiome_effect, shape = depth_flag)) +
  geom_hline(yintercept = 0.75, linetype = 3, colour = "grey50") +
  geom_point(alpha = 0.75) + facet_wrap(~ level) +
  scale_colour_manual(values = effect_cols, drop = FALSE) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), name = "Depth-flagged") +
  labs(x = "Allen purity (dominant next-finer label)", y = "clusterboot stability (mean Jaccard)",
       size = "Atlas weight", colour = "Microbiome effect", title = "Known types or new states?")
save_fig(p, "purity_stability", 13, 5)

# 6. depth QC
p <- depth %>% left_join(man %>% dplyr::select(unit_id, level), by = "unit_id") %>%
  mutate(level = lvl(level)) %>%
  ggplot(aes(delta_log10_umi)) + geom_histogram(bins = 40, fill = "#27455C") +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = 2, colour = "#FF6C0C") + facet_wrap(~ level, scales = "free_y") +
  labs(x = "Component median log10 UMI − rest of unit", y = "Components",
       title = "Sequencing-depth association of components (dashed: flag threshold)")
save_fig(p, "depth_qc", 12, 4)

# 7. meta-programs across units
if (nrow(mp$genes)) {
  mp_lab <- mp$genes %>% group_by(meta_program) %>% slice_max(frac_members, n = 6, with_ties = FALSE) %>%
    summarise(lbl = paste(gene_name, collapse = ", "))
  p <- mp$membership %>% filter(!is.na(meta_program)) %>%
    left_join(man %>% dplyr::select(unit_id, level, tissue, name), by = "unit_id") %>%
    distinct(meta_program, unit_id, level, tissue, name) %>%
    left_join(mp_lab, by = "meta_program") %>%
    mutate(mp = glue("MP{meta_program}: {lbl}"), level = lvl(level)) %>%
    ggplot(aes(str_trunc(glue("{tissue} {name}"), 32), mp, fill = level)) + geom_tile(colour = "white") +
    scale_fill_manual(values = c(class = "#27455C", subclass = "#FF9240", supertype = "#FF6C0C")) +
    labs(x = NULL, y = NULL, title = "Recurrent meta-programs across cell types") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5))
  save_fig(p, "meta_programs", 16, 7)
}

# 8. cross-level consistency
if (nrow(xl)) {
  p <- xl %>% group_by(coarse_level, fine_level, coarse_unit, coarse_comp) %>%
    summarise(top_map = max(p), n_fine = n_distinct(fine_unit[p >= 0.1]), .groups = "drop") %>%
    mutate(pair = glue("{coarse_level} → {fine_level}"), n_fine = factor(pmin(n_fine, 4), levels = 0:4)) %>%
    ggplot(aes(top_map, fill = n_fine)) + geom_histogram(bins = 20) + facet_wrap(~ pair) +
    scale_fill_manual(values = c(`0` = "#E5E5E5", `1` = "#27455C", `2` = "#FF9240", `3` = "#FF6C0C", `4` = "#B33C00"),
                      name = "Finer units ≥10%", drop = FALSE) +
    labs(x = "Largest share of a coarse component in one finer-level component",
         y = "Coarse-level components", title = "Do coarse subpopulations resolve into finer-level states?")
  save_fig(p, "cross_level", 11, 4.5)
}
log_msg("Synthesis complete")
